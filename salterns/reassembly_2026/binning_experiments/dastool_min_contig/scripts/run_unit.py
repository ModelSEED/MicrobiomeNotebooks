#!/usr/bin/env python
"""Local DAS Tool 1.1.7 re-runs of one re-assembly unit with different minimum contig lengths (no KBase writes).

Variants (all with the exact DAS Tool settings of ../../local/dastool_worker.py: diamond, --score_threshold 0.5,
--duplicate_penalty 0.6, --megabin_penalty 0.5, --write_bin_evals, --write_bins):
  V0  the three binners as run (must reproduce work/dastool/<unit>/attempt*/out/<unit>_DASTool_contig2bin.tsv)
  V1  MaxBin2 bins with contigs < 2500 bp removed (post-hoc approximation of MaxBin2 at min_contig_length 2500)
  V2  all three binners with contigs < 3000 bp removed
  V3  all three binners with contigs < 5000 bp removed
  V4  MetaBAT2 + CONCOCT only (no MaxBin2)

Single-copy-gene (SCG) calls are reused: DAS Tool 1.1.7 with --proteins P --resume only checks that P exists and is
non-empty and then reads P.bacteria.scg / P.archaea.scg (the prodigal + diamond step is skipped). The stored SCG
tables of the production run (work/dastool/<unit>/attempt*/out/<unit>_proteins.faa.{bacteria,archaea}.scg) are copied
next to a one-line stub P, so every variant (V0 included) is a selection-only run on identical SCG calls.
With --fresh-v0, one extra run 'V0fresh' repeats prodigal + diamond from scratch (-t 24) to show that the stored SCG
tables themselves reproduce.

Steps: fetch the Assembly FASTA (../../local/fetch_assembly.py, md5-verified, read-only on KBase) into
assemblies/<unit>/; build the variant contig2bin inputs from the stored c2b/*.contigs2bin.tsv and <unit>.seqlength;
run the variants concurrently under nice -n 10; skani the same-sample original MAGs and all 'lost' originals against
the whole assembly (is the genome in the assembly at all?); delete the assembly FASTA (md5s kept in runs/<unit>/fetch.json).

  V5  ensemble: the bins as run PLUS length-filtered copies as extra candidate sets (MaxBin2 >= 2500, all >= 3000,
      all >= 5000), so DAS Tool picks a length-filtered bin only where it scores higher (bin ids suffixed .ge<min>)
  V6  ensemble: the bins as run PLUS all three binners >= 5000 as extra candidates

Usage: ~/Documents/py_venv/bin/python -I run_unit.py <unit> [--fresh-v0] [--keep-assembly] [--only=V5,V6]
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)                                         # binning_experiments/dastool_min_contig
RA = os.path.dirname(os.path.dirname(EXP))                          # reassembly_2026
S = os.path.dirname(RA)                                             # salterns
PY = os.path.expanduser("~/Documents/py_venv/bin/python")
DASTOOL_BIN = os.path.expanduser("~/Documents/dastool_env/env/bin")
SKANI = os.path.expanduser("~/Documents/gtdbtk_venv/bin/skani")
BINNERS = ["metabat2", "maxbin2", "concoct"]
ALL = lambda m: [(b if not m else f"{b}_ge{m}", b, m) for b in BINNERS]
VARIANTS = {   # name: [(DAS Tool label, source binner, minimum contig length kept in bp)]
    "V0": ALL(0),
    "V1": [("metabat2", "metabat2", 0), ("maxbin2", "maxbin2", 2500), ("concoct", "concoct", 0)],
    "V2": [(b, b, 3000) for b in BINNERS],
    "V3": [(b, b, 5000) for b in BINNERS],
    "V4": [("metabat2", "metabat2", 0), ("concoct", "concoct", 0)],
    # ensembles: the length-filtered bins are offered as EXTRA candidates next to the bins as run
    "V5": ALL(0) + [("maxbin2_ge2500", "maxbin2", 2500)] + ALL(3000) + ALL(5000),
    "V6": ALL(0) + ALL(5000),
}
DAS_ARGS = ["--search_engine", "diamond", "--score_threshold", "0.5", "--duplicate_penalty", "0.6",
            "--megabin_penalty", "0.5", "--write_bin_evals", "--write_bins"]


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def attempt_dir(unit):
    for a in sorted(glob.glob(os.path.join(RA, "work/dastool", unit, "attempt*")), reverse=True):
        st = json.load(open(os.path.join(a, "status.json")))
        if st.get("state") == "done":
            return a
    raise SystemExit(f"no finished DAS Tool attempt for {unit}")


def read_lengths(seqlength):
    L = {}
    with open(seqlength) as fh:
        for line in fh:
            c, n = line.rstrip("\n").split("\t")
            L[c.split()[0]] = int(n)
    return L


def build_inputs(unit, att, rundir, L, variants):
    """Variant contig2bin TSVs; returns {variant: (labels, tsvs)} and per-input filter counts. Bin ids of an extra
    (length-filtered) candidate set get the suffix '.ge<min>' so they stay distinct from the bins as run."""
    stats, spec = {}, {}
    for v in variants:
        d = os.path.join(rundir, v, "inputs")
        os.makedirs(d, exist_ok=True)
        labels, tsvs = [], []
        for lab, b, m in VARIANTS[v]:
            src = os.path.join(att, "c2b", f"{b}.contigs2bin.tsv")
            dst = os.path.join(d, f"{lab}.contigs2bin.tsv")
            suffix = f".ge{m}" if lab != b else ""
            n_in = n_out = bp_in = bp_out = 0
            with open(src) as fi, open(dst, "w") as fo:
                for line in fi:
                    c, bn = line.rstrip("\n").split("\t")[:2]
                    c = c.split()[0]
                    n_in += 1
                    bp_in += L[c]
                    if L[c] >= m:
                        fo.write(f"{c}\t{bn}{suffix}\n")
                        n_out += 1
                        bp_out += L[c]
            stats[f"{v}:{lab}"] = {"source": b, "min_len": m, "contigs_in": n_in, "contigs_kept": n_out, "bp_in": bp_in,
                                   "bp_kept": bp_out}
            if n_out:
                labels.append(lab)
                tsvs.append(dst)
        spec[v] = (labels, tsvs)
    return spec, stats


def das_cmd(labels, tsvs, fna, prefix, proteins=None, threads=2):
    cmd = ["nice", "-n", "10", os.path.join(DASTOOL_BIN, "DAS_Tool"), "-i", ",".join(tsvs), "-l", ",".join(labels),
           "-c", fna, "-o", prefix, *DAS_ARGS, "-t", str(threads)]
    if proteins:
        cmd += ["--proteins", proteins, "--resume"]
    return cmd


def read_c2b(path):
    s = set()
    with open(path) as fh:
        for line in fh:
            c, b = line.rstrip("\n").split("\t")[:2]
            s.add((c, b))
    return s


def compare(stored, new):
    a, b = read_c2b(stored), read_c2b(new)
    ab, bb = {}, {}
    for c, x in a:
        ab.setdefault(x, set()).add(c)
    for c, x in b:
        bb.setdefault(x, set()).add(c)
    same_bins = sum(1 for k in ab if bb.get(k) == ab[k])
    return {"stored_pairs": len(a), "new_pairs": len(b), "pairs_only_stored": len(a - b), "pairs_only_new": len(b - a),
            "stored_bins": len(ab), "new_bins": len(bb), "bins_identical": same_bins, "identical": a == b}


def main():
    unit = sys.argv[1]
    fresh = "--fresh-v0" in sys.argv
    keep = "--keep-assembly" in sys.argv
    only = next((a.split("=", 1)[1].split(",") for a in sys.argv if a.startswith("--only=")), list(VARIANTS))
    att = attempt_dir(unit)
    rundir = os.path.join(EXP, "runs", unit)
    os.makedirs(rundir, exist_ok=True)
    stored = json.load(open(os.path.join(att, "fetch.json")))
    t0 = time.time()

    # 1. assembly FASTA, in its own directory (downloaded data), md5-verified against the Workspace object
    adir = os.path.join(EXP, "assemblies", unit)
    os.makedirs(adir, exist_ok=True)
    fna = os.path.join(adir, f"{unit}.assembly.fna")
    if not os.path.exists(fna):
        log(f"{unit}: fetching {stored['ref']}")
        out = subprocess.run([PY, "-I", os.path.join(RA, "local/fetch_assembly.py"), stored["ref"], fna],
                             check=True, capture_output=True, text=True).stdout
        rep = json.loads(out)
        assert rep["download_md5"] == stored["download_md5"], ("download md5 differs from the production run", rep)
        assert rep["assembly_md5_recomputed"] == stored["assembly_md5_recomputed"], rep
        json.dump(rep, open(os.path.join(rundir, "fetch.json"), "w"), indent=1)
    t_fetch = time.time() - t0

    # 2. inputs; stored SCG tables next to a stub protein file
    L = read_lengths(os.path.join(att, "out", f"{unit}.seqlength"))
    spec, fstats = build_inputs(unit, att, rundir, L, only)
    sdir = os.path.join(rundir, "scg")
    os.makedirs(sdir, exist_ok=True)
    prot = os.path.join(sdir, f"{unit}_proteins.faa")
    with open(prot, "w") as fh:
        fh.write(">stub_1\nM\n")    # never read with --resume; DAS Tool only checks that it exists and is non-empty
    for k in ("bacteria", "archaea"):
        shutil.copy2(os.path.join(att, "out", f"{unit}_proteins.faa.{k}.scg"), f"{prot}.{k}.scg")

    # 3. DAS Tool, variants concurrently (selection-only: R is single-threaded)
    procs, timing = {}, {}
    runs = dict(spec)
    if fresh:
        runs["V0fresh"] = spec["V0"]
    for v, (labels, tsvs) in runs.items():
        od = os.path.join(rundir, v, "out")
        if os.path.isdir(od):
            shutil.rmtree(od)
        os.makedirs(od)
        prefix = os.path.join(od, unit)
        cmd = das_cmd(labels, tsvs, fna, prefix, None if v == "V0fresh" else prot, 24 if v == "V0fresh" else 2)
        env = dict(os.environ, PATH=DASTOOL_BIN + os.pathsep + os.environ.get("PATH", ""))
        fh = open(os.path.join(rundir, v, "das_tool.stdout.log"), "w")
        fh.write("$ " + " ".join(cmd) + "\n")
        fh.flush()
        procs[v] = (subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env), time.time(), fh)
        log(f"{unit}: started {v} ({','.join(labels)})")
    rcs = {}
    for v, (p, ts, fh) in procs.items():
        rcs[v] = p.wait()
        fh.close()
        timing[v] = round(time.time() - ts, 1)
        log(f"{unit}: {v} exit {rcs[v]} after {timing[v]} s")

    # 4. checks: --resume really skipped the SCG search; V0 reproduces the stored contig2bin
    checks = {}
    for v in runs:
        dlog = os.path.join(rundir, v, "out", f"{unit}_DASTool.log")
        txt = open(dlog).read() if os.path.exists(dlog) else ""
        checks[v] = {"exit": rcs[v], "seconds": timing[v],
                     "scg_skipped": txt.count("Skipping SCG prediction") if v != "V0fresh" else None,
                     "no_bins": "No bins with bin-score" in txt}
    sc = os.path.join(att, "out", f"{unit}_DASTool_contig2bin.tsv")
    if "V0" in runs:
        checks["V0"]["vs_stored"] = compare(sc, os.path.join(rundir, "V0", "out", f"{unit}_DASTool_contig2bin.tsv"))
    if fresh:
        checks["V0fresh"]["vs_stored"] = compare(sc, os.path.join(rundir, "V0fresh", "out", f"{unit}_DASTool_contig2bin.tsv"))
        scg = {}
        for k in ("bacteria", "archaea"):
            a = open(os.path.join(att, "out", f"{unit}_proteins.faa.{k}.scg")).read().splitlines()
            b = open(os.path.join(rundir, "V0fresh", "out", f"{unit}_proteins.faa.{k}.scg")).read().splitlines()
            scg[k] = {"stored": len(a), "fresh": len(b), "identical_as_sets": set(a) == set(b),
                      "only_stored": len(set(a) - set(b)), "only_fresh": len(set(b) - set(a))}
        checks["V0fresh"]["scg_vs_stored"] = scg
        pf = os.path.join(rundir, "V0fresh", "out", f"{unit}_proteins.faa")
        if os.path.exists(pf):
            checks["V0fresh"]["proteins_faa_bytes"] = os.path.getsize(pf)
            os.remove(pf)                                         # large; the .scg tables are kept

    # 5. is each same-sample original / each 'lost' original in the assembly at all? (skani vs the whole assembly)
    ql = os.path.join(rundir, "skani_originals_query.txt")
    if os.path.exists(ql) and "V0" in runs:
        sk = os.path.join(rundir, "skani_originals_vs_assembly.tsv")
        t1 = time.time()
        r = subprocess.run(["nice", "-n", "10", SKANI, "dist", "--ql", ql, "-r", fna, "-t", "16", "--min-af", "1",
                            "-o", sk], capture_output=True, text=True)
        checks["skani_vs_assembly"] = {"exit": r.returncode, "seconds": round(time.time() - t1, 1),
                                       "stderr_tail": r.stderr[-400:]}

    cf = os.path.join(rundir, "run_checks.json")
    if os.path.exists(cf):                       # a partial run (--only) adds to the checks of the earlier runs
        prev = json.load(open(cf))
        fstats = {**prev.get("filter_stats", {}), **fstats}
        checks = {**prev, **checks}
    checks["fetch_seconds"] = round(t_fetch, 1)
    checks["filter_stats"] = fstats
    checks["production_attempt"] = os.path.relpath(att, RA)
    checks["assembly"] = {k: stored[k] for k in ("ref", "name", "download_md5", "assembly_md5_recomputed",
                                                  "n_fasta", "bp_fasta")}
    json.dump(checks, open(cf, "w"), indent=1)
    if not keep:
        os.remove(fna)
        log(f"{unit}: deleted {fna}")
    log(f"{unit}: done in {round(time.time() - t0)} s")


if __name__ == "__main__":
    main()
