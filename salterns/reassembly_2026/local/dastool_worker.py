#!/usr/bin/env python
"""Local DAS Tool 1.1.7 refinement of one assembly unit, launched detached by ../pipeline.py (stage 'dastool').

Steps: (1) download the Assembly FASTA from the Blobstore and verify it against the Workspace object
(fetch_assembly.py); (2) write one contig2bin TSV per binner from its KBaseMetagenomes.BinnedContigs
(binnedcontigs_to_contig2bin.py); (3) run DAS Tool 1.1.7 from ~/Documents/dastool_env (diamond, score threshold 0.5,
duplicate penalty 0.6, megabin penalty 0.5, the same thresholds as KBase's DAS Tool app); (4) build and save one
KBaseMetagenomes.BinnedContigs '<unit>.DASTool.BinnedContigs' in ws 275645 (dastool_to_binnedcontigs.py builds and
schema-validates it), which the KBase extract and CheckM stages then consume unchanged.
Writes <workdir>/status.json: {"state": running|done|empty|error, ...}; the driver never edits it.

Usage: python dastool_worker.py --unit U --assembly-ref W/O/V --binner metabat2=REF [--binner ...] --workdir DIR
"""
import argparse
import json
import os
import subprocess
import sys
import time
import traceback
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
KBUTILLIB = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "KBUtilLib", "src"))
if os.path.isdir(KBUTILLIB):
    sys.path.insert(0, KBUTILLIB)
DASTOOL_BIN = os.path.expanduser("~/Documents/dastool_env/env/bin")
PY = sys.executable
THREADS = 8


def status(workdir, **kw):
    kw["utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    tmp = os.path.join(workdir, "status.json.part")
    with open(tmp, "w") as fh:
        json.dump(kw, fh, indent=1)
    os.replace(tmp, os.path.join(workdir, "status.json"))


def run(cmd, log, env=None):
    with open(log, "a") as fh:
        fh.write(f"\n$ {' '.join(cmd)}\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode


def tail(path, n=1200):
    """The end of a log, one line, so the driver can tell a KBase outage (transient) from a real failure."""
    try:
        with open(path, "rb") as fh:
            fh.seek(0, 2)
            fh.seek(max(0, fh.tell() - n))
            return fh.read().decode(errors="replace").strip().replace("\n", " | ")
    except OSError:
        return ""


SAVE_WAITS = (60, 120, 300, 600, 1200)                   # ~38 min of retries for a Workspace outage


def saved_since(ws, ws_id, obj_name, t0):
    """The info of obj_name if a version of it was saved after t0 (a save that succeeded behind a timeout)."""
    try:
        i = ws.get_object_info3({"objects": [{"wsid": ws_id, "name": obj_name}], "ignoreErrors": 1})["infos"][0]
    except Exception:
        return None
    if not i:
        return None
    return i if datetime.strptime(i[3], "%Y-%m-%dT%H:%M:%S%z").timestamp() >= t0 - 60 else None


def save_with_retry(params, log):
    from kbutillib.kb_genome_utils import KBGenomeUtils
    ws = KBGenomeUtils(kb_version="prod").ws_client()
    obj_name, t0 = params["objects"][0]["name"], time.time()
    for k, wait in enumerate(SAVE_WAITS + (None,)):
        try:
            return ws.save_objects(params)[0]
        except Exception as e:
            with open(log, "a") as fh:
                fh.write(f"save_objects try {k + 1} failed: {type(e).__name__}: {str(e)[:500]}\n")
            if wait is None:
                raise
            time.sleep(wait)
            info = saved_since(ws, params["id"], obj_name, t0)
            if info:                                     # do not save a second version of the same object
                return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit", required=True)
    ap.add_argument("--assembly-ref", required=True)
    ap.add_argument("--binner", action="append", required=True, help="label=BinnedContigs ref")
    ap.add_argument("--workdir", required=True)
    a = ap.parse_args()
    wd = os.path.abspath(a.workdir)
    os.makedirs(wd, exist_ok=True)
    log = os.path.join(wd, "worker.log")
    try:
        status(wd, state="running", pid=os.getpid(), step="fetch")
        from fetch_assembly import fetch
        fna = os.path.join(wd, f"{a.unit}.assembly.fna")
        rep = fetch(a.assembly_ref, fna, contig_check=True)
        json.dump(rep, open(os.path.join(wd, "fetch.json"), "w"), indent=1)

        status(wd, state="running", pid=os.getpid(), step="contig2bin")
        c2b_dir = os.path.join(wd, "c2b")
        if run([PY, os.path.join(HERE, "binnedcontigs_to_contig2bin.py"), fna, a.assembly_ref, *a.binner, "--out", c2b_dir],
               log) != 0:
            raise RuntimeError(f"contig2bin conversion failed; worker.log ends: {tail(log)}")
        labels, tsvs = [], []
        for spec in a.binner:                                   # skip a binner that produced no binned contig
            lab = spec.split("=", 1)[0]
            tsv = os.path.join(c2b_dir, f"{lab}.contigs2bin.tsv")
            if os.path.getsize(tsv) > 0:
                labels.append(lab)
                tsvs.append(tsv)
        if not tsvs:
            status(wd, state="empty", reason="no binner produced a binned contig")
            return

        status(wd, state="running", pid=os.getpid(), step="das_tool", inputs=labels)
        out_prefix = os.path.join(wd, "out", a.unit)
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        env = dict(os.environ, PATH=DASTOOL_BIN + os.pathsep + os.environ.get("PATH", ""))
        rc = run(["nice", "-n", "10", os.path.join(DASTOOL_BIN, "DAS_Tool"), "-i", ",".join(tsvs), "-l", ",".join(labels),
                  "-c", fna, "-o", out_prefix, "--search_engine", "diamond", "--score_threshold", "0.5",
                  "--duplicate_penalty", "0.6", "--megabin_penalty", "0.5", "--write_bin_evals", "--write_bins",
                  "-t", str(THREADS)], log, env=env)
        dlog = out_prefix + "_DASTool.log"
        dtext = open(dlog).read() if os.path.exists(dlog) else ""
        if "No bins with bin-score" in dtext:
            status(wd, state="empty", reason="DAS Tool: no bins with bin-score above the threshold")
            return
        if rc != 0:
            raise RuntimeError(f"DAS_Tool exit {rc} (see worker.log and {dlog})")

        status(wd, state="running", pid=os.getpid(), step="build_binnedcontigs")
        save_json = os.path.join(wd, f"{a.unit}.save_objects.json")
        cmd = [PY, os.path.join(HERE, "dastool_to_binnedcontigs.py"), "--fasta", fna, "--assembly-upa", a.assembly_ref,
               "--dastool-prefix", out_prefix, "--sample", a.unit, "--out", save_json]
        for spec in a.binner:
            if spec.split("=", 1)[0] in labels:
                cmd += ["--binner", spec]
        if run(cmd, log) != 0:
            raise RuntimeError(f"building the BinnedContigs object failed; worker.log ends: {tail(log)}")
        built = json.load(open(save_json))

        status(wd, state="running", pid=os.getpid(), step="save")
        info = save_with_retry(built["save_objects_params"][0], log)
        ref = f"{info[6]}/{info[0]}/{info[4]}"
        for f in (fna, out_prefix + "_proteins.faa"):            # large intermediates; bins FASTA are kept
            if os.path.exists(f):
                os.remove(f)
        status(wd, state="done", output_ref=ref, output_name=info[1], n_bins=len(built["predicted_assemblies"]),
               predicted=save_json, bins_dir=out_prefix + "_DASTool_bins", inputs=labels,
               assembly_bp=rep["bp_fasta"], assembly_contigs=rep["n_fasta"])
    except Exception as e:
        try:
            with open(log, "a") as fh:
                fh.write(traceback.format_exc())
            status(wd, state="error", error=f"{type(e).__name__}: {str(e)[:1500]}")
        except Exception:        # e.g. a full disk: the driver sees a dead worker without a final status and retries
            pass


if __name__ == "__main__":
    main()
