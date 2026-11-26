#!/usr/bin/env python
"""GROW hybrid worker: runs ONE KBase-app step locally with the app's own tool versions and command lines, then saves
the result into the same narrative by direct Workspace save (user decision 2026-11-23, both BRIEFs).

The driver (assembly_binning_2026 / multisample_binning_2026 pipeline.py) writes <workdir>/spec.json and starts
    setsid ~/Documents/py_venv/bin/python -u hw.py <workdir>
The worker never touches a ledger. It writes <workdir>/result.json at the end:
    {"status": "done"|"empty"|"failed", "ref": UPA|None, "error": str, "intrinsic": bool, "versions": {...}, ...}
and <workdir>/worker.log. Large intermediates live under <workdir>/tmp and are removed on success.

Stages (spec["stage"]) and the KBase app each one replaces (commands: see ~/Documents/metag_env/VERSIONS.md):
  assemble  MEGAHIT/run_megahit e2e6d557      -> KBaseGenomeAnnotations.Assembly (AssemblyUtil 59855bcb shape)
  metabat   metabat/run_metabat bd317234       -> KBaseMetagenomes.BinnedContigs (MetagenomeUtils 49489d11 shape)
  maxbin    kb_maxbin/run_maxbin2 3574f0b2     -> BinnedContigs
  concoct   kb_concoct/run_kb_concoct 936c9bfe -> BinnedContigs
  extract   MetagenomeUtils extract 54c6bc35 + kb_Msuite CheckM 7495e02e: bin Assemblies (CheckM completeness /
            contamination as object metadata) + KBaseSets.AssemblySet + a KBaseReport.Report whose file_links carry
            CheckM_summary_table.tsv.zip and full_output.zip exactly as kb_Msuite's report does (so the driver's
            collect_checkm_select.py reads it unchanged). The driver's checkm stage then adopts that report.
"""
import csv
import glob
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import traceback
import uuid
import zipfile
from collections import Counter
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kbwrite as K  # noqa: E402
from kbrpc import rpc  # noqa: E402  (read-only)

ENV = os.path.expanduser("~/Documents/metag_env")
ENVS = os.path.join(ENV, "mamba_root", "envs")
MB = os.path.join(ENV, "src", "kbase-metabat", "metabat_binary")
MAXBIN = os.path.join(ENV, "maxbin")
MEGAHIT = os.path.join(ENV, "megahit", "MEGAHIT-1.2.9-Linux-x86_64-static", "bin", "megahit")
CHECKM_DATA = os.path.join(ENV, "checkm_data")
HOST = os.uname().nodename

VERSIONS = {
    "assemble": {"MEGAHIT": "1.2.9 (official static Linux build, as kb_megahit e2e6d557)"},
    "metabat": {"MetaBAT": "2.12.1", "bbmap": "37.90", "samtools": "1.8", "jgi_summarize_bam_contig_depths": "2.12.1",
                "source": "kbase-metabat bd317234 metabat_binary/ (same binaries)"},
    "maxbin": {"MaxBin": "2.2.4", "bowtie2": "2.2.3", "HMMER": "3.1b1", "FragGeneScan": "1.30",
               "source": "MaxBin-2.2.4 source + its auxiliary versions (as kb_maxbin 3574f0b2)"},
    "concoct": {"CONCOCT": "1.1.0", "bowtie2": "2.3.0", "samtools": "1.3.1",
                "source": "as kb_concoct 936c9bfe (CONCOCT master 2021-03 = 1.1.0; Debian stretch bowtie2/samtools)"},
    "extract": {"CheckM": "1.1.2", "HMMER": "3.1b2", "Prodigal": "2.6.3", "pplacer": "1.1.alpha19",
                "checkm_data": "2015_01_16", "source": "as kb_Msuite 7495e02e (lineage_wf --reduced_tree)"},
}
VERSIONS["checkm"] = VERSIONS["extract"]
SERVICE = {  # the KBase app each local step stands in for (recorded in provenance)
    "assemble": ("MEGAHIT", "run_megahit", "e2e6d557da8a92d8b9a6cb4893a2a7293a3068bf"),
    "metabat": ("metabat", "run_metabat", "bd317234abfc301ff4d0fd63b1abe7615052b78f"),
    "maxbin": ("kb_maxbin", "run_max_bin", "3574f0b2cc343f4f530bd52b86c8cf8b6b1427f3"),
    "concoct": ("kb_concoct", "run_kb_concoct", "936c9bfed7e07dd1ddcb76ac4fa775b73d60afbb"),
    "extract": ("MetagenomeUtils", "extract_binned_contigs_as_assembly", "54c6bc35fa899823966a3b0ab972e2cdfd2a36a2"),
    "checkm": ("kb_Msuite", "run_checkM_lineage_wf", "7495e02e3380354255c859823aa52028e6e0a433"),
}


class Intrinsic(Exception):
    """A failure of the tool on this input (not of the machine): the driver may move to the next variant."""


class Superseded(Exception):
    """Another version of this step's output arrived (from KBase) while it ran: not pushed (no-duplicate rule)."""

    def __init__(self, name, ref, saved):
        super().__init__(f"{name} already has version {ref} saved {saved}, during this run: local result not pushed")
        self.name, self.ref = name, ref


WORKER_START = datetime.now(timezone.utc)


def pre_push(spec, names):
    """No-duplicate rule (user decision 2026-11-25): before a push, no other version of these output names may have
    been saved since this worker started. If one was, the first (validated by the driver) is kept: raise Superseded."""
    for n in [x for x in names if x]:
        infos = rpc("ws", "Workspace.get_object_info3", [{"objects": [{"wsid": int(spec["ws_id"]), "name": n}],
                                                          "ignoreErrors": 1}])[0]["infos"]
        i = infos[0] if infos else None
        if i and datetime.strptime(i[3], "%Y-%m-%dT%H:%M:%S%z") >= WORKER_START:
            raise Superseded(n, f"{i[6]}/{i[0]}/{i[4]}", i[3])


class Empty(Exception):
    """The tool ran and produced nothing (no bins): the stage ends 'empty'."""


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


LOG = None


def log(msg):
    line = f"[{now()}] {msg}"
    print(line, flush=True)


def sh(cmd, cwd=None, env_bin=None, check=True, log_to=None):
    """Run a bash command (logged verbatim). env_bin: extra PATH dirs first."""
    env = dict(os.environ)
    if env_bin:
        env["PATH"] = ":".join(env_bin) + ":" + env["PATH"]
    log(f"RUN: {cmd}")
    t = time.time()
    out = open(log_to, "ab") if log_to else None
    p = subprocess.run(["/bin/bash", "-o", "pipefail", "-c", cmd], cwd=cwd, env=env,
                       stdout=out or None, stderr=subprocess.STDOUT if out else None)
    if out:
        out.close()
    log(f"exit {p.returncode} after {time.time() - t:.0f} s")
    if check and p.returncode != 0:
        raise RuntimeError(f"command failed ({p.returncode}): {cmd[:300]}")
    return p.returncode


# ------------------------------------------------------------------------------------------------ inputs
def get_obj(ref, included=None):
    spec = {"ref": ref}
    if included:
        spec["included"] = included
    return rpc("ws", "Workspace.get_objects2", [{"objects": [spec]}])[0]["data"][0]


def fetch_reads(ref, cache):
    """A PairedEndLibrary's stored FASTQ (interleaved .gz in these narratives), md5-verified, cached by ref."""
    os.makedirs(cache, exist_ok=True)
    o = get_obj(ref, ["/lib1", "/lib2", "/interleaved", "/read_count", "/total_bases"])
    d = o["data"]
    if d.get("lib2") or not d.get("interleaved"):
        raise RuntimeError(f"{ref}: not a single interleaved file; not handled")
    f = d["lib1"]["file"]
    path = os.path.join(cache, ref.replace("/", "_") + "." + f["file_name"].split(".", 1)[1])
    if not (os.path.exists(path) and os.path.exists(path + ".md5ok")):
        log(f"downloading reads {ref} ({d['lib1'].get('size', 0) / 1e9:.1f} GB) -> {path}")
        K.download_blob(f["id"], path, expect_md5=f.get("remote_md5"))
        open(path + ".md5ok", "w").write(f.get("remote_md5") or "")
    return path, d


def fetch_assembly_fasta(ref, path):
    """The Assembly's blob FASTA as stored (what AssemblyUtil.get_assembly_as_fasta hands the apps), md5-verified."""
    o = get_obj(ref, ["/fasta_handle_info", "/fasta_handle_ref", "/md5", "/num_contigs", "/dna_size"])
    h = o["data"]["fasta_handle_info"]["handle"]
    if not os.path.exists(path):
        log(f"downloading assembly {ref} -> {path}")
        K.download_blob(h["id"], path, expect_md5=h.get("remote_md5"))
    return path, o["data"]


def read_fasta(path):
    """[(id, description, SEQUENCE)] in file order, Bio.SeqIO semantics (id = first token)."""
    from Bio import SeqIO
    return [(r.id, r.description[len(r.id):].strip(), str(r.seq)) for r in SeqIO.parse(path, "fasta")]


# ------------------------------------------------------------------------------------------------ object builders
def assembly_data(fasta_path, assembly_id):
    """AssemblyUtil 59855bcb FastaToAssembly._parse_fasta + _build_assembly_object (fasta_handle_info added later)."""
    from Bio import SeqIO
    total, base_counts, md5s, contigs = 0, {"A": 0, "G": 0, "C": 0, "T": 0}, [], {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        ci = {"contig_id": rec.id, "name": rec.id, "description": rec.description[len(rec.id):].strip(),
              "length": len(rec.seq)}
        total += ci["length"]
        cnt = dict(Counter(seq))
        for ch, n in cnt.items():
            base_counts[ch] = base_counts.get(ch, 0) + n
            if ch not in "ACGTNRYSWKMBDHV.-":
                raise ValueError(f"non nucleic acid character {ch}")
        if "N" in cnt:
            ci["Ncount"] = cnt["N"]
        ci["md5"] = hashlib.md5(seq.encode()).hexdigest()
        md5s.append(ci["md5"])
        ci["gc_content"] = round(float(cnt.get("G", 0) + cnt.get("C", 0)) / float(ci["length"]), 5)
        if rec.id in contigs:
            raise ValueError(f"duplicate contig id {rec.id}")
        contigs[rec.id] = ci
    if not contigs:
        raise Empty("no contigs")
    return {"md5": hashlib.md5(",".join(sorted(md5s)).encode()).hexdigest(), "base_counts": base_counts,
            "dna_size": total, "gc_content": round(float(base_counts["G"] + base_counts["C"]) / float(total), 5)
            if total else None, "contigs": contigs, "num_contigs": len(contigs), "assembly_id": assembly_id,
            "type": "Unknown"}


def save_assembly(spec, fasta_path, name, prov, meta=None, upload_name=None):
    data = assembly_data(fasta_path, name)
    fh = K.blob_upload(fasta_path, upload_name or os.path.basename(fasta_path))
    data["fasta_handle_ref"] = fh["handle"]["hid"]
    data["fasta_handle_info"] = fh
    info = K.save_object(spec["ws_id"], "KBaseGenomeAnnotations.Assembly", name, data, prov, meta)
    return info, data


def provenance(spec, stage, method_params, input_refs, extra=""):
    svc, meth, commit = SERVICE[stage]
    where = HOST + {"TF": " (CheckM on TF)", "poplar": " (MEGAHIT on poplar)"}.get(spec.get("remote"), "")
    desc = (f"Run LOCALLY on host {where} (GROW hybrid local/KBase mode, user decision 2026-11-23) in place of KBase "
            f"app {svc}.{meth} ({commit[:8]}), same tool versions and parameters: "
            + "; ".join(f"{k} {v}" for k, v in VERSIONS.get(stage if stage != "checkm" else "extract", {}).items())
            + (f". {extra}" if extra else ""))
    return [{"time": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+0000"),
             "service": svc, "service_ver": commit, "method": meth, "method_params": [method_params],
             "input_ws_objects": list(input_refs), "script": f"grow/hybrid/hw.py {stage} (local, host={HOST})",
             "script_ver": spec.get("script_ver", ""), "description": desc[:1000],
             "script_command_line": f"hw.py {spec['workdir']}"}]


def binned_contigs_from_dir(bin_dir, assembly_ref, contig_info):
    """MetagenomeUtils 49489d11 MetagenomeFileUtils.file_to_binned_contigs, minus the save."""
    from Bio import SeqIO
    files = os.listdir(bin_dir)
    bin_ids = [f for f in files if re.match(r'.*\.\d{3}\.fasta', f)]
    summaries = [os.path.join(bin_dir, f) for f in files if f.endswith(".summary")]
    lines_all = []
    for s in summaries:
        lines_all.append(open(s).readlines())
    bins = []
    for bid in bin_ids:
        gc = sum_len = cov = None
        for lines in lines_all:                           # _generate_contig_bin_summary: last summary file wins
            for line in lines:
                ll = line.split("\t")
                if ll[0] == bid:
                    if len(ll) == 5:
                        gc, sum_len, cov = round(float(ll[4]) / 100, 5), int(ll[3]), round(float(ll[2].partition('%')[0]) / 100, 5)
                    elif len(ll) == 4:
                        gc, sum_len, cov = round(float(ll[3]) / 100, 5), int(ll[2]), round(float(ll[1].partition('%')[0]) / 100, 5)
        if gc is None:
            raise RuntimeError(f"bin {bid} missing from the summary file (MetagenomeUtils would fail)")
        contigs = {}
        for rec in SeqIO.parse(os.path.join(bin_dir, bid), "fasta"):
            c = contig_info.get(rec.id)
            if c is None:
                raise RuntimeError(f"contig {rec.id} of {bid} is not in the assembly")
            contigs[rec.id] = {"gc": c["gc_content"], "len": c["length"]}
        bins.append({"bid": bid, "contigs": contigs, "n_contigs": len(contigs), "gc": gc, "sum_contig_len": sum_len,
                     "cov": cov})
    total = 0
    for lines in lines_all:
        for line in lines[1:]:
            ll = line.split("\t")
            if len(ll) == 5:
                total += int(ll[3])
            elif len(ll) == 4:
                total += int(ll[2])
    return {"assembly_ref": assembly_ref, "bins": bins, "total_contig_len": total}


def contig_info_of(fasta):
    info = {}
    from Bio import SeqIO
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper()
        info[rec.id] = {"gc_content": round(float(seq.count("G") + seq.count("C")) / float(len(seq)), 5),
                        "length": len(seq)}
    return info


def check_against_assembly(contig_info, asm_ref):
    """The local FASTA must be the Assembly object: same contig ids and lengths (spot-check of gc)."""
    o = get_obj(asm_ref, ["/num_contigs", "/dna_size"])["data"]
    if o["num_contigs"] != len(contig_info) or o["dna_size"] != sum(c["length"] for c in contig_info.values()):
        raise RuntimeError(f"local FASTA does not match {asm_ref}: {len(contig_info)} contigs vs {o['num_contigs']}")


# ------------------------------------------------------------------------------------------------ stages
class Overflow(Exception):
    """FM's memory budget would be exceeded: stopped before swap / oomd; the driver requeues the step to poplar,
    which resumes from the kept MEGAHIT checkpoints (user decision 2026-11-25)."""


def megahit_cmd(spec, reads, out, tmp, mg=None):
    p = spec["params"]
    cmd = f"{mg or MEGAHIT} --12 {reads}"
    if p.get("megahit_parameter_preset"):
        cmd += f" --presets {p['megahit_parameter_preset']}"
    if p.get("min_count"):
        cmd += f" --min-count {p['min_count']}"
    if p.get("k_list"):
        cmd += f" --k-list {','.join(str(k) for k in p['k_list'])}"
    for k in ("k_min", "k_max", "k_step"):
        if p.get(k):
            cmd += f" --{k.replace('_', '-')} {p[k]}"
    return cmd + (f" --min-contig-len {p.get('min_contig_length', 1000)} --num-cpu-threads {spec['threads']}"
                  f" -m {int(spec['mem_gb'] * 1e9)} -o {out} --tmp-dir {tmp}")


def session_rss_gb(sid):
    tot = 0
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as fh:
                f = fh.read().rsplit(")", 1)[1].split()
            if int(f[3]) == sid:                                   # field 6 (session id)
                tot += int(f[21]) * os.sysconf("SC_PAGE_SIZE")     # field 24 (rss pages)
        except (OSError, IndexError, ValueError):
            pass
    return tot / 1e9


def run_watched(cmd, log_to, rss_limit_gb, floor_gb):
    """Run cmd in its own session; stop it (SIGTERM to the session, then SIGKILL) when its total RSS passes
    rss_limit_gb or MemAvailable drops below floor_gb. Returns (rc, reason or None)."""
    log(f"RUN (memory watchdog: RSS <= {rss_limit_gb} GB, MemAvailable >= {floor_gb} GB): {cmd}")
    t = time.time()
    out = open(log_to, "ab")
    pr = subprocess.Popen(["/bin/bash", "-o", "pipefail", "-c", cmd], stdout=out, stderr=subprocess.STDOUT,
                          start_new_session=True)
    reason, peak = None, 0.0
    while pr.poll() is None:
        time.sleep(10)
        rss = session_rss_gb(pr.pid)
        peak = max(peak, rss)
        with open("/proc/meminfo") as fh:
            avail = next(int(x.split()[1]) for x in fh if x.startswith("MemAvailable:")) / 1048576
        if rss > rss_limit_gb or avail < floor_gb:
            reason = (f"memory watchdog: RSS {rss:.0f} GB (limit {rss_limit_gb}), MemAvailable {avail:.0f} GB "
                      f"(floor {floor_gb}) - stopped before swap/oomd")
            log(reason)
            os.killpg(pr.pid, signal.SIGTERM)
            for _ in range(30):
                if pr.poll() is not None:
                    break
                time.sleep(2)
            if pr.poll() is None:
                os.killpg(pr.pid, signal.SIGKILL)
            pr.wait()
            break
    out.close()
    log(f"exit {pr.returncode} after {time.time() - t:.0f} s (peak RSS {peak:.0f} GB)")
    return pr.returncode, reason, peak


# ------------------------------------------------------------------------------------------------ poplar (overflow)
POPLAR_ROOT = "/scratch/afreiburger/grow_hybrid"
POPLAR_MEGAHIT = f"{POPLAR_ROOT}/env/megahit/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit"
POPLAR_RUNNER = f"{POPLAR_ROOT}/bin/poplar_runner.py"
POPLAR_MIN_FREE_GB = 200
POPLAR_STAGE = "/scratch1/afreiburger/grow_hybrid_stage"   # NFS: reads, control files, outputs (2026-11-25, 3rd decision)


def poplar(cmd, timeout=180, check=True, stdin=None):
    p = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "poplar", cmd], capture_output=True,
                       text=True, timeout=timeout, input=stdin)
    if check and p.returncode != 0:
        raise RuntimeError(f"ssh poplar failed ({p.returncode}): {p.stderr[-300:]}")
    return p.stdout


def megahit_resume_plan(old_tmp):
    """Check that a stopped MEGAHIT 1.2.9 run can be continued elsewhere: options.json + checkpoints.txt present, the
    checkpoint numbering maps onto the step sequence, and the paths it records are the old run's. Returns
    (options, out_dir, temp_dir, next_step, next_k) or raises."""
    out = os.path.join(old_tmp, "megahit_out")
    opt = json.load(open(os.path.join(out, "options.json")))
    done = [int(l.split()[0]) for l in open(os.path.join(out, "checkpoints.txt")) if l.strip().endswith("done")]
    if done != list(range(len(done))) or len(done) < 4:
        raise RuntimeError(f"checkpoints not a complete prefix past the first graph: {done}")
    if os.path.realpath(opt["out_dir"]) != os.path.realpath(out) or not os.path.isdir(opt["temp_dir"]) or \
            not opt["temp_dir"].startswith(old_tmp) or opt.get("no_local") or not opt.get("kmin_1pass"):
        raise RuntimeError("options.json paths/flags not as expected for a resume")
    # 0 library file, 1 build library, 2 first graph (1-pass), 3 assemble(k_min); then per k: local, iterate,
    # build_graph, assemble
    n = len(done)
    i = (n - 4) // 4 + 1
    step = ("local", "iterate", "build_graph", "assemble")[(n - 4) % 4]
    if i >= len(opt["k_list"]):
        raise RuntimeError("all k done; merge step only - rerun the merge instead")
    return opt, opt["out_dir"], opt["temp_dir"], step, opt["k_list"][i]


def run_megahit_poplar(spec, tmp, reads_ref):
    """MEGAHIT on poplar (overflow host). Either resumes a stopped run from its checkpoints (spec['resume_from'], the
    old attempt's tmp on FM: copied over, recorded paths rewritten, nothing else changed) or runs from the KBase
    library, which poplar pulls itself (token on stdin, never in a file). FM polls; a dropped MFA ControlMaster only
    delays polling (up to 72 h). final.contigs.fa comes back md5-checked; FM saves the Assembly."""
    rid = os.path.basename(spec["workdir"].rstrip("/")) + "_" + spec["unit"]
    # 3rd decision 2026-11-25: run dir (job, logs, reads, finished outputs) on /scratch1 (NFS, room); MEGAHIT's work
    # dir (output + tmp) on local /scratch when the gate found room there, else on /scratch1
    rd = f"{POPLAR_STAGE}/runs/{rid}"
    on_local = spec.get("poplar_work_fs", "/scratch") == "/scratch"
    work = f"{POPLAR_ROOT}/runs/{rid}" if on_local else f"{POPLAR_STAGE}/work/{rid}"
    free = float(poplar(f"df -BG --output=avail {work.rsplit('/runs/', 1)[0].rsplit('/work/', 1)[0]} | tail -1")
                 .strip().rstrip("G"))
    floor = POPLAR_MIN_FREE_GB if on_local else 1000
    poplar(f"mkdir -p {rd}/reads {rd}/out {work}")
    final = f"{rd}/out/final.contigs.fa"
    job = {"downloads": [], "outputs": [final], "min_free_gb": 150 if on_local else 500, "free_path": work,
           "moves": [[f"{work}/megahit_out/final.contigs.fa", final], [f"{work}/megahit_out/log", f"{rd}/out/megahit_log"],
                     [f"{work}/megahit_out/options.json", f"{rd}/out/options.json"]],
           "rm_after": [work]}
    log(f"poplar run dir {rd} (/scratch1); MEGAHIT work dir {work}")
    resume = spec.get("resume_from")
    note = ""
    if resume:
        opt, out_dir, temp_dir, step, next_k = megahit_resume_plan(resume)
        need = sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(resume) for f in fs) / 1e9
        if free - need < floor + 100:
            raise RuntimeError(f"poplar {work}: {free:.0f} GB free; resume needs {need:.0f} GB + growth")
        new_out, new_tmp = f"{work}/megahit_out", f"{work}/{os.path.basename(temp_dir)}"
        reads_new = [f"{rd}/reads/{os.path.basename(x)}" for x in opt["pe12"]]
        log(f"resume on poplar: {len(open(os.path.join(out_dir, 'checkpoints.txt')).readlines())} checkpoints done, "
            f"next step {step} toward k{next_k}; copying {need:.0f} GB")
        sh(f"rsync -a --partial {out_dir}/ poplar:{new_out}/ && rsync -a --partial {temp_dir}/ poplar:{new_tmp}/"
           f" && ssh -o BatchMode=yes poplar mkdir -p {rd}/reads && rsync -a --partial {' '.join(opt['pe12'])} poplar:{rd}/reads/")
        # verify the copy: sizes + md5 of every file under 1 GB, size of the big ones
        for src, dst in ((out_dir, new_out), (temp_dir, new_tmp)):
            l = sorted((os.path.relpath(os.path.join(r, f), src), os.path.getsize(os.path.join(r, f)))
                       for r, _, fs in os.walk(src) for f in fs)
            rem = sorted((x.split("\t")[0], int(x.split("\t")[1])) for x in
                         poplar(f"cd {dst} && find . -type f -printf '%P\\t%s\\n'").splitlines() if x)
            if l != rem:
                raise RuntimeError(f"copy of {src} to poplar differs (files/sizes)")
        lsum = sh_out(f"md5sum {' '.join(opt['pe12'])} {temp_dir}/reads.lib.bin | cut -d' ' -f1")
        rsum = poplar(f"md5sum {' '.join(reads_new)} {new_tmp}/reads.lib.bin | cut -d' ' -f1")
        if lsum.split() != rsum.split():
            raise RuntimeError("md5 of reads / reads.lib.bin differ after the copy")
        # rewrite only the recorded paths (options.json; reads.lib text) - everything else is byte-identical
        o2 = dict(opt)
        o2.update(out_dir=new_out, temp_dir=new_tmp, pe12=reads_new)
        open(os.path.join(tmp, "options.json"), "w").write(json.dumps(o2))
        lib = open(os.path.join(temp_dir, "reads.lib")).read().replace(temp_dir, new_tmp)
        for a, b in zip(opt["pe12"], reads_new):
            lib = lib.replace(a, b)
        open(os.path.join(tmp, "reads.lib"), "w").write(lib)
        info = open(os.path.join(temp_dir, "reads.lib.lib_info")).read()
        for a, b in zip(opt["pe12"], reads_new):
            info = info.replace(a, b)
        open(os.path.join(tmp, "reads.lib.lib_info"), "w").write(info)
        sh(f"scp -q -o BatchMode=yes {tmp}/options.json poplar:{new_out}/options.json && "
           f"scp -q -o BatchMode=yes {tmp}/reads.lib {tmp}/reads.lib.lib_info poplar:{new_tmp}/")
        if step in ("local", "iterate"):                    # a killed iterate leaves partial edges for k_next that
            poplar(f"rm -rf {new_tmp}/k{next_k}")           # build_graph would read: drop them (recreated by iterate)
        poplar(f"grep -rl /home/freiburger {new_out}/options.json {new_tmp}/reads.lib {new_tmp}/reads.lib.lib_info; true")
        job["cmd"] = f"nice -n 10 {POPLAR_MEGAHIT} --continue -o {new_out} > megahit.log 2>&1"
        note = (f"Resumed on poplar from FM's checkpoints ({step} toward k{next_k}; options.json/reads.lib paths "
                f"rewritten, data byte-identical; MEGAHIT --continue reuses the recorded options: threads "
                f"{opt['num_cpu_threads']}, -m {opt['memory'] / 1e9:.0f} GB)")
    else:
        o = get_obj(reads_ref, ["/lib1", "/interleaved"])["data"]
        f = o["lib1"]["file"]
        need = o["lib1"].get("size", 0) / 1e9 * 7                 # work dir (reads are on /scratch1)
        if free - need < floor:
            raise RuntimeError(f"poplar {work}: {free:.0f} GB free < {floor} + {need:.0f} GB")
        rp = f"{rd}/reads/{reads_ref.replace('/', '_')}.{f['file_name'].split('.', 1)[1]}"
        poplar(f"mkdir -p {rd}/reads")
        job["downloads"] = [{"node": f["id"], "md5": f.get("remote_md5"), "path": rp}]
        job["cmd"] = "nice -n 10 " + megahit_cmd(spec, rp, f"{work}/megahit_out", work, POPLAR_MEGAHIT) + " > megahit.log 2>&1"
        note = f"Run on poplar from the KBase library (pulled there, md5-checked); threads {spec['threads']}, -m {spec['mem_gb']} GB"
    open(os.path.join(tmp, "job.json"), "w").write(json.dumps(job, indent=1))
    sh(f"scp -q -o BatchMode=yes {tmp}/job.json poplar:{rd}/job.json")
    poplar(f"python3 {POPLAR_RUNNER} {rd}", stdin=(K.token() + "\n") if job["downloads"] else "\n")
    log(f"MEGAHIT started on poplar in {rd}")
    misses = 0
    while True:
        time.sleep(300)
        try:
            st = poplar(f"cat {rd}/STATE 2>/dev/null; true", timeout=90).strip()
            misses = 0
        except Exception as e:  # noqa: BLE001
            misses += 1
            if misses in (1, 12) or misses % 72 == 0:
                log(f"poplar poll failed ({misses}; the MFA ControlMaster may need a new login): {e}")
            if misses > 864:                                    # 72 h unreachable
                raise RuntimeError("poplar unreachable for 72 h; the run there may still finish - check by hand")
            continue
        if st == "done":
            break
        if st == "failed":
            det = poplar(f"cat {rd}/state.json; tail -c 1500 {rd}/megahit.log; true")
            raise Intrinsic(f"MEGAHIT failed on poplar: {det[-1500:]}")
    js = json.loads(poplar(f"cat {rd}/state.json"))
    out = os.path.join(tmp, "megahit_out")
    os.makedirs(out, exist_ok=True)
    fa = os.path.join(out, "final.contigs.fa")
    sh(f"scp -q -o BatchMode=yes poplar:{final} {fa} && "
       f"scp -q -o BatchMode=yes poplar:{rd}/out/megahit_log {spec['workdir']}/megahit_poplar.log")
    if md5(fa) != js["md5"][final]:
        raise RuntimeError("final.contigs.fa md5 differs after the copy from poplar")
    log(f"final.contigs.fa copied back from poplar (md5 ok); MEGAHIT took {js.get('seconds')} s there")
    poplar(f"rm -rf {rd} {work}", timeout=900, check=False)   # reads, outputs, any work-dir remains
    log(f"poplar run dir {rd} removed (contigs on FM, md5-verified)")   # (user decisions 2026-11-25, 2nd/3rd)
    return fa, note, None


def sh_out(cmd):
    return subprocess.run(["/bin/bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True, check=True).stdout


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def stage_assemble(spec, tmp):
    p = spec["params"]
    prd = None
    if spec.get("remote") == "poplar":
        fa, where_note, prd = run_megahit_poplar(spec, tmp, p["read_library_ref"])
        host_note = "MEGAHIT on poplar.cels.anl.gov (overflow host), saved from FM. " + where_note
    else:
        reads, rd = fetch_reads(p["read_library_ref"], spec["reads_cache"])
        out = os.path.join(tmp, "megahit_out")
        shutil.rmtree(out, ignore_errors=True)
        rc, over, peak = run_watched(megahit_cmd(spec, reads, out, tmp), os.path.join(spec["workdir"], "megahit.log"),
                                     spec.get("rss_limit_gb", spec["mem_gb"]), spec.get("mem_floor_gb", 20))
        if over:
            raise Overflow(over)
        fa = os.path.join(out, "final.contigs.fa")
        if rc != 0 or not os.path.exists(fa):
            tail = open(os.path.join(spec["workdir"], "megahit.log"), errors="replace").read()[-1500:]
            if re.search(r"No space left", tail):
                raise RuntimeError(f"megahit: disk full (rc {rc}): {tail[-400:]}")
            # memory failures are intrinsic, as the driver treats them for KBase: the ladder's min-count 2 needs less
            raise Intrinsic(f"megahit failed (rc {rc}): {tail[-400:]}")
        host_note = (f"Reads: the stored interleaved FASTQ of {p['read_library_ref']} given as --12 (KBase passes the "
                     f"same pairs as -1/-2); threads {spec['threads']} and -m {spec['mem_gb']} GB instead of 63 / 0.9 "
                     f"(resource settings only); peak RSS {peak:.0f} GB")
    if os.path.getsize(fa) == 0:
        raise Intrinsic("megahit produced an empty final.contigs.fa")
    prov = provenance(spec, "assemble", {**p, "workspace_name": spec["ws_name"]}, [p["read_library_ref"]], host_note)
    pre_push(spec, [spec["out_name"]])
    info, data = save_assembly(spec, fa, spec["out_name"], prov, upload_name="final.contigs.fa")
    if prd:                                               # saved and validated by the save: free poplar's /scratch
        poplar(f"rm -rf {prd}", timeout=600, check=False)
        log(f"poplar run dir {prd} removed")
    if spec.get("resume_from") and os.path.isdir(spec["resume_from"]):   # the FM checkpoints it resumed from
        shutil.rmtree(spec["resume_from"], ignore_errors=True)
        log(f"FM checkpoints {spec['resume_from']} removed (assembly saved)")
    return {"ref": K.upa(info), "size": data["dna_size"], "n_contigs": data["num_contigs"]}


def stage_metabat(spec, tmp):
    p = spec["params"]
    asm = os.path.join(tmp, "assembly.fa")
    fetch_assembly_fasta(p["assembly_ref"], asm)
    ci = contig_info_of(asm)
    check_against_assembly(ci, p["assembly_ref"])
    res = os.path.join(tmp, "Metabat_results")
    os.makedirs(res, exist_ok=True)
    env = [MB, os.path.join(MB, "bbmap"), os.path.join(ENVS, "env_map", "bin")]
    dicts = []
    for i, ref in enumerate(p["reads_list"]):
        reads, _ = fetch_reads(ref, spec["reads_cache"])
        bam = os.path.join(res, f"lib{i}.sam.sorted.bam")
        # bbmap.sh ... out=<sam> then samtools view | sort  ==  the same records streamed (no SAM on disk)
        sh(f"bash {MB}/bbmap/bbmap.sh -Xmx30g fast threads={spec['threads']} ref={asm} in={reads} out=stdout.sam "
           f"mappedonly nodisk overwrite 2>{res}/lib{i}.bbmap.log | samtools view -F 0x04 -uS - | "
           f"samtools sort -@ 4 -m 2G -T {res}/lib{i}.sorttmp - -o {bam}", cwd=res, env_bin=env)
        depth = os.path.join(res, f"lib{i}.depth.txt")
        sh(f"jgi_summarize_bam_contig_depths --outputDepth {depth} --minContigLength {p['min_contig_length']} "
           f"--minContigDepth 1 {bam}", cwd=res, env_bin=env)
        os.remove(bam)
        d = {}                                            # MetabatUtil.createDictFromDepthfile
        with open(depth) as f:
            header = f.readline().rstrip().split("\t")
            for line in f:
                vals = line.rstrip().split("\t")
                d[vals[0].split()[0] if " " in vals[0] else vals[0]] = vals[1:]
            d["header"] = header
        dicts.append(json.loads(json.dumps(d)))           # KBase passes them through JSON files
    merged = os.path.join(tmp, "merged.depth.txt")        # CombineDepthTables.createFinalDepthfile
    header_list, main = ["contigName", "contigLen", "totalAvgDepth"], {}
    for dd in dicts:
        for c in dd:
            if c != "header":
                main[c] = [dd[c][0], 0]
    for dd in dicts:
        header_list += [dd["header"][3], dd["header"][4]]
        for c in main:
            if c in dd:
                try:
                    old = main[c][1]
                except IndexError:
                    old = 0
                    main[c].append(0)
                main[c][1] = float(old) + float(dd[c][2])
                main[c].append(dd[c][2])
                main[c].append(dd[c][3])
            else:
                main[c].append(float(0))
                main[c].append(float(0))
    with open(merged, "w") as f:
        f.write("\t".join(header_list) + "\n")
        for c in ci:                                      # assembly order (SeqIO ids)
            if c in main:
                f.write("{}\t{}\n".format(c, "\t".join(str(x) for x in main[c])))
    if sum(1 for _ in open(merged)) < 2:
        raise Intrinsic("merged depth file is empty")
    sh(f"/bin/bash {MB}/metabat_wrapper.sh -a {asm} -m {p['min_contig_length']} -o {res} -d {merged} -t {spec['threads']}",
       cwd=tmp, env_bin=env, log_to=os.path.join(spec["workdir"], "metabat.log"))
    bins_dir = os.path.join(res, "Bins")
    return save_binned(spec, "metabat", bins_dir, p, ci)


def stage_maxbin(spec, tmp):
    p = spec["params"]
    asm = os.path.join(tmp, "assembly.fa")
    fetch_assembly_fasta(p["assembly_ref"], asm)
    ci = contig_info_of(asm)
    check_against_assembly(ci, p["assembly_ref"])
    if spec.get("remote") == "TF":
        return maxbin_tf(spec, tmp, p, asm, ci)
    res = os.path.join(tmp, "result")
    os.makedirs(res, exist_ok=True)
    bt = os.path.join(MAXBIN, "bowtie2-2.2.3")
    # run_MaxBin.pl -reads_list: formatBowtie2 (bowtie2-build <contig> <out>.idx), then per reads file
    # 'bowtie2 -p T -x <out>.idx -U <reads> -S <sam>' and getsam(<sam>) -> abundance. Same commands here, with the SAM
    # streamed through a FIFO into MaxBin's own getsam (it reads the SAM once, sequentially); then the abundance files
    # go to run_MaxBin.pl -abund_list, which is what -reads_list does internally after mapping.
    sh(f"{bt}/bowtie2-build {asm} {res}/Bin.idx > {res}/bowtie2build.out 2> {res}/bowtie2build.err", cwd=res)
    abunds = []
    for i, ref in enumerate(p["reads_list"]):
        reads, _ = fetch_reads(ref, spec["reads_cache"])
        fifo = os.path.join(res, f"Bin.sam{i}")
        if os.path.exists(fifo):
            os.remove(fifo)
        os.mkfifo(fifo)
        ab = os.path.join(res, f"reads.abund{i + 1}")
        getsam = (f"perl -e 'require \"{MAXBIN}/MaxBin-2.2.4/_getabund.pl\"; getsam(\"{fifo}\", \"{ab}\");'")
        sh(f"{bt}/bowtie2 -p {spec['threads']} -x {res}/Bin.idx -U {reads} -S {fifo} 1>{fifo}.bowtie2.out "
           f"2>{fifo}.bowtie2.err & bpid=$! ; {getsam} ; gs=$? ; wait $bpid ; bs=$? ; "
           f"echo bowtie2=$bs getsam=$gs ; [ $bs -eq 0 ] && [ $gs -eq 0 ]", cwd=res)
        if os.path.getsize(ab) == 0:
            raise RuntimeError(f"empty abundance file {ab}")
        os.remove(fifo)
        abunds.append(ab)
    with open(os.path.join(res, "abund_list.txt"), "w") as f:
        f.write("\n".join(abunds) + "\n")
    for f in glob.glob(os.path.join(res, "Bin.idx*")):
        os.remove(f)
    out = os.path.join(tmp, "maxbin_out")
    os.makedirs(out, exist_ok=True)
    mlog = os.path.join(spec["workdir"], "maxbin.log")
    opts = (f"-prob_threshold {p['prob_threshold']} -markerset {p['markerset']} -min_contig_length "
            f"{p['min_contig_length']}")
    rc = sh(f"perl {MAXBIN}/MaxBin-2.2.4/run_MaxBin.pl -contig {asm} -out Bin -abund_list {res}/abund_list.txt "
            f"-thread {p.get('thread', 32) if spec['threads'] >= p.get('thread', 32) else spec['threads']} {opts}",
            cwd=out, log_to=mlog, check=False)
    return finish_maxbin(spec, p, out, mlog, rc, ci)


def finish_maxbin(spec, p, out, mlog, rc, ci):
    txt = open(mlog, errors="replace").read()
    if "cannot be binned" in txt:       # kb_maxbin fails with this; the driver maps it to 'empty' (deterministic)
        raise Empty("MaxBin2: too few marker genes in this assembly to bin it (dataset cannot be binned)")
    if rc != 0:
        raise RuntimeError(f"run_MaxBin.pl failed ({rc}): {txt[-400:]}")
    return save_binned(spec, "maxbin", out, p, ci)


def concoct_tf(spec, tmp, p, clean, res):
    """CONCOCT's mapping, coverage table and clustering on TF (2026-11-25), the same commands as the local path with
    TF's env_bt2 (bowtie2 2.3.0 / samtools 1.3.1) and env_concoct (CONCOCT 1.1.0), both created from FM's exact
    package lists. merge_cutup_clustering / extract_fasta_bins then run on FM as before."""
    reads = [fetch_reads(ref, spec["reads_cache"])[0] for ref in p["reads_list"]]
    t = spec.get("remote_threads", 10)
    seed = int.from_bytes(os.urandom(4), "little") % 1000000000      # kb_concoct draws a random bowtie2-build seed
    E = f"{TF_ENV}/mamba_root/envs"
    steps = ["O=$PATH", "mkdir -p concoct_output_dir", f"export PATH={E}/env_bt2/bin:$O",
             f"nice -n 10 bowtie2-build -f assembly_clean.fa --threads {t} --seed {seed} idx.bt2 > bt2build.log 2>&1"]
    for i, r in enumerate(reads):
        n = f"lib{i}_" + os.path.basename(r)
        fw, rv = f"lib{i}_forward.fastq.gz", f"lib{i}_reverse.fastq.gz"
        steps += [f"zcat reads/{n} | paste - - - - - - - - | tee >(cut -f 1-4 | tr '\\t' '\\n' | gzip -1 > {fw}) | "
                  f"cut -f 5-8 | tr '\\t' '\\n' | gzip -1 > {rv}",
                  f"nice -n 10 bowtie2 -x idx.bt2 -1 {fw} -2 {rv} --threads {t} 2> concoct_output_dir/lib{i}.bt2.log | "
                  f"samtools view -F 0x04 -uS - | samtools sort -@ 4 -m 2G -T concoct_output_dir/lib{i}.sorttmp - "
                  f"-o concoct_output_dir/lib{i}_sorted.bam",
                  f"samtools index concoct_output_dir/lib{i}_sorted.bam", f"rm -f {fw} {rv} reads/{n}"]
    cmd = (f"concoct --composition_file concoct_output_dir/split_contigs.fa -l {p['min_contig_length']} -b "
           f"concoct_output_dir --coverage_file concoct_output_dir/coverage_table.tsv -t {t} "
           f"-k {p['kmer_size']} -c {p['max_clusters_for_vgmm']} -i {p['max_iterations_for_vgmm']} "
           f"--total_percentage_pca {p['total_percentage_pca']}")
    if int(p.get("no_cov_normalization") or 0):
        cmd += " --no_cov_normalization"
    if int(p.get("no_total_coverage") or 0):
        cmd += " --no_total_coverage"
    steps += ["rm -f idx.bt2*", f"export PATH={E}/env_concoct/bin:$O",
              f"cut_up_fasta.py assembly_clean.fa -c {p['contig_split_size']} -o {p['contig_split_overlap']} "
              f"--merge_last -b temp.bed > concoct_output_dir/split_contigs.fa",
              f"export PATH={E}/env_concoct/bin:{E}/env_bt2/bin:$O",
              "concoct_coverage_table.py temp.bed concoct_output_dir/*_sorted.bam > concoct_output_dir/coverage_table.tsv",
              "rm -f concoct_output_dir/*_sorted.bam*", f"export PATH={E}/env_concoct/bin:$O",
              f"nice -n 10 {cmd} > concoct.log 2>&1"]
    rc = run_tf_job(spec, tmp, "CONCOCT", [clean], " && ".join(steps), ["concoct_output_dir", "concoct.log",
                                                                          "bt2build.log"],
                    big=[(r, f"lib{i}_" + os.path.basename(r)) for i, r in enumerate(reads)])
    lg = os.path.join(tmp, "concoct.log")
    if os.path.exists(lg):
        shutil.copy(lg, os.path.join(spec["workdir"], "concoct.log"))
    if rc != 0:
        txt = open(lg, errors="replace").read()[-400:] if os.path.exists(lg) else "(mapping / coverage step)"
        raise RuntimeError(f"CONCOCT on TF failed ({rc}): {txt}")


def stage_concoct(spec, tmp):
    p = spec["params"]
    asm = os.path.join(tmp, "assembly.fa")
    fetch_assembly_fasta(p["assembly_ref"], asm)
    ci = contig_info_of(asm)
    check_against_assembly(ci, p["assembly_ref"])
    # retrieve_and_clean_assembly: headers renamed to unique space-free ids and reverted in the bins afterwards
    # (FastaHeaderRename) -> equivalent to writing a copy with id-only headers, the ids themselves are unchanged.
    clean = os.path.join(tmp, "assembly_clean.fa")
    with open(asm) as fi, open(clean, "w") as fo:
        for line in fi:
            fo.write(line.split()[0] + "\n" if line.startswith(">") else line)
    res = os.path.join(tmp, "concoct_output_dir")
    os.makedirs(res, exist_ok=True)
    cc = [os.path.join(ENVS, "env_concoct", "bin")]
    if spec.get("remote") == "TF":                        # mapping, coverage table and clustering on TF
        concoct_tf(spec, tmp, p, clean, res)
    else:
        bt = [os.path.join(ENVS, "env_bt2", "bin")]
        cc = [os.path.join(ENVS, "env_concoct", "bin")]
        seed = int.from_bytes(os.urandom(4), "little") % 1000000000      # kb_concoct draws a random bowtie2-build seed
        sh(f"bowtie2-build -f {clean} --threads {spec['threads']} --seed {seed} {tmp}/idx.bt2 > {tmp}/bt2build.log 2>&1",
           cwd=tmp, env_bin=bt)
        for i, ref in enumerate(p["reads_list"]):
            reads, _ = fetch_reads(ref, spec["reads_cache"])
            fw, rv = os.path.join(tmp, f"lib{i}_forward.fastq.gz"), os.path.join(tmp, f"lib{i}_reverse.fastq.gz")
            # deinterleave_fastq.sh (kb_concoct): line-group split 1-4 / 5-8
            sh(f"zcat {reads} | paste - - - - - - - - | tee >(cut -f 1-4 | tr '\\t' '\\n' | gzip -1 > {fw}) | "
               f"cut -f 5-8 | tr '\\t' '\\n' | gzip -1 > {rv}", cwd=tmp)
            sh(f"bowtie2 -x {tmp}/idx.bt2 -1 {fw} -2 {rv} --threads {spec['threads']} 2> {res}/lib{i}.bt2.log | "
               f"samtools view -F 0x04 -uS - | samtools sort -@ 4 -m 2G -T {res}/lib{i}.sorttmp - -o {res}/lib{i}_sorted.bam",
               cwd=tmp, env_bin=bt)
            sh(f"samtools index {res}/lib{i}_sorted.bam", cwd=tmp, env_bin=bt)
            os.remove(fw)
            os.remove(rv)
        for f in glob.glob(os.path.join(tmp, "idx.bt2*")):
            os.remove(f)
        sh(f"cut_up_fasta.py {clean} -c {p['contig_split_size']} -o {p['contig_split_overlap']} --merge_last -b temp.bed "
           f"> concoct_output_dir/split_contigs.fa", cwd=tmp, env_bin=cc)
        sh("concoct_coverage_table.py temp.bed concoct_output_dir/*_sorted.bam > concoct_output_dir/coverage_table.tsv",
           cwd=tmp, env_bin=cc + bt)
        for f in glob.glob(os.path.join(res, "*_sorted.bam*")):
            os.remove(f)
        cmd = (f"concoct --composition_file concoct_output_dir/split_contigs.fa -l {p['min_contig_length']} -b "
               f"concoct_output_dir --coverage_file concoct_output_dir/coverage_table.tsv -t {spec['threads']} "
               f"-k {p['kmer_size']} -c {p['max_clusters_for_vgmm']} -i {p['max_iterations_for_vgmm']} "
               f"--total_percentage_pca {p['total_percentage_pca']}")
        if int(p.get("no_cov_normalization") or 0):
            cmd += " --no_cov_normalization"
        if int(p.get("no_total_coverage") or 0):
            cmd += " --no_total_coverage"
        sh(cmd, cwd=tmp, env_bin=cc, log_to=os.path.join(spec["workdir"], "concoct.log"))
    sh(f"merge_cutup_clustering.py concoct_output_dir/clustering_gt{p['min_contig_length']}.csv > "
       f"concoct_output_dir/clustering_merged.csv", cwd=tmp, env_bin=cc)
    fb = os.path.join(res, "final_bins")
    os.makedirs(fb, exist_ok=True)                        # generate_concoct_extract_fasta_bins_command: mkdir first
    sh(f"extract_fasta_bins.py {clean} concoct_output_dir/clustering_merged.csv --output_path {fb}", cwd=tmp,
       env_bin=cc)
    for f in os.listdir(fb):                              # rename_and_standardize_bin_names
        if f.endswith(".fa"):
            os.rename(os.path.join(fb, f), os.path.join(fb, "bin." + f.split(".fa")[0].zfill(3) + ".fasta"))
    with open(os.path.join(fb, "binned_contig.summary"), "w") as f:   # make_binned_contig_summary_file_...
        f.write("Bin name\tCompleteness\tGenome size\tGC content\n")
        for fn in sorted(os.listdir(fb)):
            if fn.endswith(".fasta"):
                seqs = [s for _, _, s in read_fasta(os.path.join(fb, fn))]
                acgt = sum(sum(s.upper().count(b) for b in "ACGT") for s in seqs)
                gc = sum(s.upper().count("G") + s.upper().count("C") for s in seqs)
                # bbtools stats.sh format=3: contig_bp = all bases; gc_avg = GC/ACGT; the app writes gc_avg*100
                f.write(f"{fn}\t0\t{sum(len(s) for s in seqs)}\t{(gc / acgt if acgt else 0) * 100}\n")
    return save_binned(spec, "concoct", fb, p, ci)


def maxbin_tf(spec, tmp, p, asm, ci):
    """MaxBin2 entirely on TF (2026-11-25: FM is often loaded past its cores by other users' jobs): the same commands
    as the local path (bowtie2 2.2.3 -U per library into MaxBin's getsam through a FIFO, then run_MaxBin.pl
    -abund_list), with TF's md5-identical MaxBin2 tree. FM only fetches the reads (KBase, md5-checked), ships them and
    saves the result."""
    reads = [fetch_reads(ref, spec["reads_cache"])[0] for ref in p["reads_list"]]
    t = spec.get("remote_threads", 10)
    M, B = f"{TF_ENV}/maxbin", f"{TF_ENV}/maxbin/bowtie2-2.2.3"
    steps = [f"mkdir -p result maxbin_out && cd result",
             f"{B}/bowtie2-build ../assembly.fa Bin.idx > bowtie2build.out 2> bowtie2build.err"]
    for i, r in enumerate(reads):
        n = f"lib{i}_" + os.path.basename(r)
        steps.append(f"mkfifo Bin.sam{i} && {{ nice -n 10 {B}/bowtie2 -p {t} -x Bin.idx -U ../reads/{n} -S Bin.sam{i} "
                     f"1>Bin.sam{i}.bowtie2.out 2>Bin.sam{i}.bowtie2.err & bpid=$! ; perl -e 'require \"$ENV{{HOME}}/metag_env/maxbin/"
                     f"MaxBin-2.2.4/_getabund.pl\"; getsam(\"Bin.sam{i}\", \"reads.abund{i + 1}\");' ; gs=$? ; wait $bpid ; bs=$? ; "
                     f"echo bowtie2=$bs getsam=$gs ; [ $bs -eq 0 ] && [ $gs -eq 0 ] ; }} && test -s reads.abund{i + 1} "
                     f"&& rm -f Bin.sam{i}")
    steps += [f"printf '%s\\n' {' '.join(f'$PWD/reads.abund{i + 1}' for i in range(len(reads)))} > abund_list.txt",
              "rm -f Bin.idx*", "cd ../maxbin_out",
              f"nice -n 10 perl {M}/MaxBin-2.2.4/run_MaxBin.pl -contig ../assembly.fa -out Bin -abund_list "
              f"../result/abund_list.txt -thread {min(p.get('thread', 32), t)} -prob_threshold {p['prob_threshold']} "
              f"-markerset {p['markerset']} -min_contig_length {p['min_contig_length']} > ../maxbin.log 2>&1"]
    rc = run_tf_job(spec, tmp, "MaxBin2", [asm], " && ".join(steps), ["maxbin_out", "maxbin.log", "result"],
                    big=[(r, f"lib{i}_" + os.path.basename(r)) for i, r in enumerate(reads)])
    mlog = os.path.join(spec["workdir"], "maxbin.log")
    if os.path.exists(os.path.join(tmp, "maxbin.log")):
        shutil.copy(os.path.join(tmp, "maxbin.log"), mlog)
    else:                                                 # the mapping failed before MaxBin ran
        errs = "".join(open(f, errors="replace").read()[-400:] for f in glob.glob(os.path.join(tmp, "result", "*.err")))
        open(mlog, "w").write(f"mapping on TF failed (rc {rc}):\n{errs}")
    return finish_maxbin(spec, p, os.path.join(tmp, "maxbin_out"), mlog, rc, ci)


def save_binned(spec, stage, bin_dir, p, ci):
    data = binned_contigs_from_dir(bin_dir, p["assembly_ref"], ci)
    if not data["bins"]:
        raise Empty(f"{stage}: no bins")
    mp = {**p, "workspace_name": spec["ws_name"]}
    prov = provenance(spec, stage, mp, [p["assembly_ref"]] + list(p["reads_list"]),
                      "Mapping output streamed straight into samtools instead of an intermediate SAM file (same records)")
    pre_push(spec, [spec["out_name"]])
    info = K.save_object(spec["ws_id"], "KBaseMetagenomes.BinnedContigs", spec["out_name"], data, prov)
    return {"ref": K.upa(info), "n_bins": len(data["bins"])}


def stage_extract(spec, tmp):
    """extract_binned_contigs_as_assembly + run_checkM_lineage_wf, locally, in one step."""
    p = spec["params"]
    bc_ref, suffix, set_name = p["binned_contig_obj_ref"], p["assembly_suffix"], p["assembly_set_name"]
    bc = get_obj(bc_ref)["data"]
    asm = os.path.join(tmp, "assembly.fa")
    fetch_assembly_fasta(bc["assembly_ref"], asm)
    seqs = {i: s.upper() for i, _, s in read_fasta(asm)}
    wanted = p["extracted_assemblies"].split(",")
    bins = {b["bid"]: b for b in bc["bins"]}
    bdir = os.path.join(tmp, "bins")
    os.makedirs(bdir, exist_ok=True)
    files = {}
    for bid in wanted:                                    # binned_contigs_to_file: '>id\nSEQ\n' per contig
        name = bid + suffix
        fp = os.path.join(bdir, name + ".fna")
        with open(fp, "w") as f:
            for cid in bins[bid]["contigs"]:
                f.write(f">{cid}\n{seqs[cid]}\n")
        files[bid] = fp
    ext, qa, out = run_checkm(spec, tmp, bdir, {b + suffix for b in wanted})
    # bin Assemblies (AssemblyUtil save, as extract does: name = bid + suffix, upload file name = bid)
    refs, created = [], []
    prov = provenance(spec, "extract", {"workspace_name": spec["ws_name"], "binned_contig_obj_ref": bc_ref,
                                        "extracted_assemblies": p["extracted_assemblies"], "assembly_suffix": suffix,
                                        "assembly_set_name": set_name}, [bc_ref],
                      "CheckM 1.1.2 lineage_wf --reduced_tree (as kb_Msuite 7495e02e) run locally on these bins; "
                      "completeness/contamination are in each bin Assembly's metadata")
    objs = []                                             # AssemblyUtil: all blobs first, then ONE save_objects
    for bid in wanted:
        e = ext[bid + suffix]
        meta = {"CheckM completeness": f"{float(e['Completeness']):.2f}",
                "CheckM contamination": f"{float(e['Contamination']):.2f}",
                "CheckM marker lineage": str(e.get("marker lineage")),
                "CheckM version": "1.1.2 reduced_tree (local)"}
        data = assembly_data(files[bid], bid + suffix)
        fh = K.blob_upload(files[bid], bid)
        data["fasta_handle_ref"] = fh["handle"]["hid"]
        data["fasta_handle_info"] = fh
        objs.append({"type": "KBaseGenomeAnnotations.Assembly", "name": bid + suffix, "data": data,
                     "provenance": prov, "meta": meta})
    pre_push(spec, [set_name] + [o["name"] for o in objs])
    for info in K.save_objects(spec["ws_id"], objs):
        refs.append(K.upa(info))
        created.append({"ref": K.upa(info), "description": info[1]})
    if len(refs) > 1:
        sinfo = K.save_object(spec["ws_id"], "KBaseSets.AssemblySet", set_name,
                              {"description": f"binned assemblies from {p.get('binned_contig_name', bc_ref)}",
                               "items": [{"ref": r, "label": ""} for r in refs]},
                              provenance(spec, "extract", {"workspace": spec["ws_name"], "output_object_name": set_name},
                                         refs, "AssemblySet saved directly (SetAPI.save_assembly_set_v1 shape)"))
        out_ref = K.upa(sinfo)
    else:
        out_ref = refs[0]
    rref, n_pass = save_checkm_report(spec, tmp, out, qa, ext, out_ref, set_name, created)
    return {"ref": out_ref, "n_bins": len(wanted), "checkm_report_ref": rref, "n_pass_50_10": n_pass,
            "bin_refs": refs}


def run_checkm(spec, tmp, bdir, expect):
    """CheckM 1.1.2, as kb_Msuite: 'checkm lineage_wf -t 4 --reduced_tree <bins> <out>' (bin id = Assembly name).
    --tab_table -f only adds a copy of the qa table; threads only change speed."""
    import ast
    out = os.path.join(tmp, "checkm_out")
    shutil.rmtree(out, ignore_errors=True)
    if spec.get("remote") == "TF":
        run_checkm_tf(spec, tmp, bdir, out)
    else:
        ck = [os.path.join(ENVS, "env_checkm", "bin")]
        sh(f"checkm lineage_wf -t {spec['threads']} --pplacer_threads {min(spec['threads'], 8)} --reduced_tree -x fna "
           f"--tab_table -f {tmp}/qa.tsv {bdir} {out}", cwd=tmp, env_bin=ck,
           log_to=os.path.join(spec["workdir"], "checkm.log"))
    ext = {}
    for line in open(os.path.join(out, "storage", "bin_stats_ext.tsv")):
        if line.strip():
            k, js = line.split("\t", 1)
            ext[k] = ast.literal_eval(js)
    qa = list(csv.DictReader(open(os.path.join(tmp, "qa.tsv")), delimiter="\t"))
    if set(ext) != set(expect) or len(qa) != len(expect):
        raise RuntimeError(f"CheckM reported {len(ext)}/{len(qa)} of {len(expect)} bins")
    return ext, qa, out


# ------------------------------------------------------------------------------------------------ TF (remote compute)
TF_ROOT = "tf_runs/grow_hybrid"                       # on TF, under ~andrew; no KBase token ever goes there
TF_ENV = "$HOME/metag_env"


def tf(cmd, timeout=120, check=True):
    p = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "TF", cmd], capture_output=True,
                       text=True, timeout=timeout)
    if check and p.returncode != 0:
        raise RuntimeError(f"ssh TF failed ({p.returncode}): {p.stderr[-300:]}")
    return p.stdout


def run_checkm_tf(spec, tmp, bdir, out):
    """CheckM on TF (user plan 2026-11-23: small input, CPU-heavy): bins (md5 list) -> TF; a detached run there
    writes state files; FM polls over ssh (a dropped connection only delays polling); outputs come back as a tar
    whose md5 TF records; host=TF goes into result.json and the provenance."""
    rid = os.path.basename(spec["workdir"].rstrip("/")) + "_" + spec["unit"]
    rdir = f"{TF_ROOT}/{rid}"
    sh(f"cd {bdir} && md5sum *.fna > ../bins.md5 && tar cf ../bins.tar *.fna", cwd=tmp)
    tf(f"mkdir -p {rdir}")
    sh(f"scp -q -o BatchMode=yes {tmp}/bins.tar {tmp}/bins.md5 TF:{rdir}/", cwd=tmp)
    thr = spec.get("remote_threads", 8)
    script = (f"cd ~/{rdir} && mkdir -p bins && tar xf bins.tar -C bins && (cd bins && md5sum -c --quiet ../bins.md5) "
              f"&& export PATH={TF_ENV}/mamba_root/envs/env_checkm/bin:$PATH && "
              f"nice -n 10 checkm lineage_wf -t {thr} --pplacer_threads {min(thr, 8)} --reduced_tree -x fna "
              f"--tab_table -f qa.tsv bins checkm_out > checkm.log 2>&1 && "
              f"tar cf result.tar qa.tsv checkm_out && md5sum result.tar > result.md5 && echo done > STATE "
              f"|| echo failed > STATE")
    tf(f"cd ~/{rdir} && rm -f STATE && echo '{script}' > run.sh && (setsid nohup bash run.sh > run.out 2>&1 < /dev/null &)")
    log(f"CheckM started on TF in ~/{rdir}")
    misses = 0
    while True:
        time.sleep(60)
        try:
            st = tf(f"cat ~/{rdir}/STATE 2>/dev/null; true", timeout=60).strip()
            misses = 0
        except Exception as e:  # noqa: BLE001
            misses += 1
            log(f"TF poll failed ({misses}): {e}")
            if misses > 180:                                  # ~3 h unreachable
                raise RuntimeError("TF unreachable for 3 h")
            continue
        if st == "done":
            break
        if st == "failed":
            raise RuntimeError("CheckM failed on TF: " + tf(f"tail -c 1500 ~/{rdir}/checkm.log; true"))
    sh(f"scp -q -o BatchMode=yes TF:{rdir}/result.tar TF:{rdir}/result.md5 {tmp}/", cwd=tmp)
    sh(f"md5sum -c --quiet result.md5 && tar xf result.tar", cwd=tmp)
    shutil.copy(os.path.join(tmp, "checkm_out", "checkm.log") if os.path.exists(os.path.join(tmp, "checkm_out", "checkm.log"))
                else os.path.join(tmp, "qa.tsv"), os.path.join(spec["workdir"], "checkm_tf.log"))
    tf(f"rm -rf ~/{rdir}")
    log("CheckM results copied back from TF (md5 verified); TF run dir removed")


def shrink_slot(spec, threads=1, mem_gb=2):
    """After the FM part of a TF step (mapping / coverage) the FM worker only waits: give the gate its room back."""
    try:
        import gate as G
        with G.Gate() as g:
            w = g.d["workers"].get(spec["workdir"])
            if w:
                w.update(threads=threads, mem_gb=mem_gb)
        log(f"FM gate slot reduced to {threads} thread / {mem_gb} GB while TF computes")
    except Exception as e:  # noqa: BLE001
        log(f"gate slot not reduced ({type(e).__name__}: {e})")


def run_tf_job(spec, tmp, label, files, script, results, big=()):
    """A detached job on TF (user plan 2026-11-23, binners added 2026-11-25): `files` (paths under tmp) go over as a
    tar with an md5 list, checked there; `script` runs in the TF run dir (nice 10, the pinned envs in TF_ENV, the same
    builds as FM: MaxBin2 tree md5-identical, env_concoct package set identical); `results` (paths relative to the run
    dir) come back as a tar whose md5 TF records. No KBase token ever goes to TF."""
    rid = os.path.basename(spec["workdir"].rstrip("/")) + "_" + spec["unit"]
    rdir = f"{TF_ROOT}/{rid}"
    rel = [os.path.relpath(f, tmp) for f in files]
    sh(f"md5sum {' '.join(rel)} > in.md5 && tar cf in.tar {' '.join(rel)} in.md5", cwd=tmp)
    tf(f"rm -rf ~/{rdir} && mkdir -p ~/{rdir}")
    sh(f"scp -q -o BatchMode=yes {tmp}/in.tar TF:{rdir}/", cwd=tmp)
    if big:                                               # reads: one copy each into reads/, md5-checked there
        tf(f"mkdir -p ~/{rdir}/reads")
        with open(os.path.join(tmp, "reads.md5"), "w") as f:
            for src, name in big:
                f.write(f"{md5(src)}  reads/{name}\n")
        for src, name in big:
            sh(f"scp -q -o BatchMode=yes {src} TF:{rdir}/reads/{name}", cwd=tmp)
        sh(f"scp -q -o BatchMode=yes {tmp}/reads.md5 TF:{rdir}/reads.md5", cwd=tmp)
    chk = "md5sum -c --quiet in.md5" + (" && md5sum -c --quiet reads.md5" if big else "")
    body = (f"cd ~/{rdir} && tar xf in.tar && {chk} && ( {script} ) ; echo $? > rc ; "
            f"tar cf result.tar rc {' '.join(results)} 2>/dev/null ; md5sum result.tar > result.md5 && echo done > STATE "
            f"|| echo failed > STATE")
    open(os.path.join(tmp, "tf_run.sh"), "w").write(body + "\n")
    sh(f"scp -q -o BatchMode=yes {tmp}/tf_run.sh TF:{rdir}/run.sh", cwd=tmp)
    tf(f"cd ~/{rdir} && rm -f STATE && (setsid nohup bash run.sh > run.out 2>&1 < /dev/null &)")
    log(f"{label} started on TF in ~/{rdir}")
    shrink_slot(spec)
    misses = 0
    while True:
        time.sleep(60)
        try:
            st = tf(f"cat ~/{rdir}/STATE 2>/dev/null; true", timeout=60).strip()
            misses = 0
        except Exception as e:  # noqa: BLE001
            misses += 1
            log(f"TF poll failed ({misses}): {e}")
            if misses > 180:
                raise RuntimeError("TF unreachable for 3 h")
            continue
        if st == "done":
            break
        if st == "failed":
            raise RuntimeError(f"{label} on TF: transfer/check failed: " + tf(f"tail -c 800 ~/{rdir}/run.out; true"))
    sh(f"scp -q -o BatchMode=yes TF:{rdir}/result.tar TF:{rdir}/result.md5 {tmp}/ && md5sum -c --quiet result.md5 && "
       f"tar xf result.tar && rm -f result.tar", cwd=tmp)
    rc = int(open(os.path.join(tmp, "rc")).read().strip() or 1)
    tf(f"rm -rf ~/{rdir}")
    log(f"{label} results copied back from TF (md5 verified, rc {rc}); TF run dir removed")
    return rc


def save_checkm_report(spec, tmp, out, qa, ext, out_ref, set_name, created):
    # CheckM report: CheckM_summary_table.tsv.zip + full_output.zip as file_links (kb_Msuite's names)
    summ = io.StringIO()
    cols = ["Bin Name", "Marker Lineage", "# Genomes", "# Markers", "# Marker Sets", "0", "1", "2", "3", "4", "5+",
            "Completeness", "Contamination"]
    w = csv.writer(summ, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for r in sorted(qa, key=lambda r: r["Bin Id"]):
        w.writerow([r["Bin Id"], r["Marker lineage"], r["# genomes"], r["# markers"], r["# marker sets"], r["0"], r["1"],
                    r["2"], r["3"], r["4"], r["5+"], round(float(r["Completeness"]), 2),
                    round(float(r["Contamination"]), 2)])
    z1 = os.path.join(tmp, "CheckM_summary_table.tsv.zip")
    with zipfile.ZipFile(z1, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("CheckM_summary_table.tsv", summ.getvalue())
    z2 = os.path.join(tmp, "full_output.zip")
    with zipfile.ZipFile(z2, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as z:
        for root, _, fs in os.walk(out):
            for fn in fs:
                full = os.path.join(root, fn)
                z.write(full, os.path.relpath(full, out))
    links = []
    for path, desc in ((z1, "TSV Summary Table from CheckM"), (z2, "Full output of CheckM")):
        fh = K.blob_upload(path, os.path.basename(path))
        links.append({"URL": f"{K.HANDLE_URL}/node/{fh['shock_id']}", "description": desc, "handle": fh["handle"]["hid"],
                      "label": "", "name": os.path.basename(path)})
    n_pass = sum(1 for e in ext.values() if float(e["Completeness"]) >= 50 and float(e["Contamination"]) < 10)
    rep = {"text_message": (f"CheckM 1.1.2 lineage_wf --reduced_tree run locally on host {HOST} (GROW hybrid mode) on "
                            f"{len(ext)} bins of {set_name}: {n_pass} at >= 50 % completeness and < 10 % "
                            f"contamination. Per-bin values are in CheckM_summary_table.tsv.zip"
                            + (" and in each bin Assembly's metadata." if created else ".")),
           "warnings": [], "objects_created": created, "file_links": links, "html_links": []}
    rname = f"{spec['unit']}.CheckM_local.report"
    rinfo = K.save_object(spec["ws_id"], "KBaseReport.Report", rname, rep,
                          provenance(spec, "checkm", {"workspace_name": spec["ws_name"], "input_ref": out_ref,
                                                      "reduced_tree": "1", "threads": str(spec["threads"])}, [out_ref]))
    return K.upa(rinfo), n_pass


def stage_checkm(spec, tmp):
    """CheckM only, for a unit whose extract already ran in KBase: the bin Assemblies of the set (or the one bare
    Assembly) are fetched from their blobs and CheckM'd; a report object is saved (the existing Assemblies' metadata
    cannot be changed after the fact, so here the values are in the report only)."""
    p = spec["params"]
    target = p["input_ref"]
    head = K.get_info(target)
    if head[2].startswith("KBaseGenomeAnnotations.Assembly"):
        items = [target]
    else:
        items = [it["ref"] for it in get_obj(target)["data"]["items"]]
    infos = rpc("ws", "Workspace.get_object_info3", [{"objects": [{"ref": r} for r in items]}])[0]["infos"]
    bdir = os.path.join(tmp, "bins")
    os.makedirs(bdir, exist_ok=True)
    for r, i in zip(items, infos):
        fetch_assembly_fasta(r, os.path.join(bdir, i[1] + ".fna"))
    ext, qa, out = run_checkm(spec, tmp, bdir, {i[1] for i in infos})
    created = []
    rref, n_pass = save_checkm_report(spec, tmp, out, qa, ext, target, head[1], created)
    return {"ref": rref, "n_bins": len(items), "n_pass_50_10": n_pass}


STAGES = {"assemble": stage_assemble, "metabat": stage_metabat, "maxbin": stage_maxbin, "concoct": stage_concoct,
          "extract": stage_extract, "checkm": stage_checkm}


def main():
    wd = os.path.abspath(sys.argv[1])
    spec = json.load(open(os.path.join(wd, "spec.json")))
    spec["workdir"] = wd
    tmp = spec.get("tmp") or os.path.join(wd, "tmp")
    os.makedirs(tmp, exist_ok=True)
    # TMPDIR only for small files and multiprocessing sockets (AF_UNIX paths must stay < 108 chars, which the deep
    # work dirs exceed: CheckM's mp.Manager failed with EOFError); every large intermediate has an explicit path in tmp
    short = f"/tmp/grow_hw_{os.getpid()}"
    os.makedirs(short, exist_ok=True)
    os.environ["TMPDIR"] = short
    res = {"stage": spec["stage"], "unit": spec["unit"], "host": HOST + {"TF": "+TF (CheckM)", "poplar": "+poplar (MEGAHIT)"}.get(spec.get("remote"), ""), "started_utc": now(),
           "versions": VERSIONS.get(spec["stage"]), "pid": os.getpid()}
    json.dump({**res, "status": "running"}, open(os.path.join(wd, "status.json"), "w"))
    try:
        out = STAGES[spec["stage"]](spec, tmp)
        res.update(status="done", **out)
    except Empty as e:
        res.update(status="empty", error=str(e))
    except Superseded as e:
        res.update(status="superseded", superseded_by=e.ref, error=str(e))
    except Overflow as e:
        res.update(status="failed", overflow=True, intrinsic=False, error=str(e))
    except Intrinsic as e:
        res.update(status="failed", intrinsic=True, error=str(e))
    except BaseException as e:  # noqa: BLE001
        res.update(status="failed", intrinsic=False, error=f"{type(e).__name__}: {e}",
                   traceback=traceback.format_exc()[-3000:])
    res["finished_utc"] = now()
    with open(os.path.join(wd, "result.json.tmp"), "w") as f:
        json.dump(res, f, indent=1)
    os.replace(os.path.join(wd, "result.json.tmp"), os.path.join(wd, "result.json"))
    if res["status"] in ("done", "empty") and not spec.get("keep_tmp"):
        shutil.rmtree(tmp, ignore_errors=True)
    shutil.rmtree(short, ignore_errors=True)
    log(f"result: {res['status']} {res.get('ref') or res.get('error', '')}")


if __name__ == "__main__":
    main()
