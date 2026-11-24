#!/usr/bin/env python
"""Where did the 'lost' original MAGs go? Compare them with every raw binner bin (MetaBAT2, MaxBin2, CONCOCT, before
DAS Tool) of the experiment unit that should contain them.

For each unit with in-scope lost originals (results/lost_originals.tsv): re-fetch the assembly (md5-verified, read-only),
write each raw binner bin as FASTA (temporary), skani the lost originals against them, and join DAS Tool's own SCG
evaluation of every input bin (production out/<unit>_allBins.eval) plus each raw bin's share of bp in contigs < 2500 bp.
Output: results/lost_trace.tsv (one row per lost original x matching raw bin, ANI >= 95). Temporary FASTA are deleted.
Usage: ~/Documents/py_venv/bin/python -I trace_lost.py
"""
import glob
import json
import os
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
RA = os.path.dirname(os.path.dirname(EXP))
RES = os.path.join(EXP, "results")
PY = os.path.expanduser("~/Documents/py_venv/bin/python")
SKANI = os.path.expanduser("~/Documents/gtdbtk_venv/bin/skani")
sys.path.insert(0, HERE)
from analyze import production_attempt  # noqa: E402


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def main():
    orig = pd.read_csv(os.path.join(RES, "originals.tsv"), sep="\t", keep_default_na=False)
    lost = pd.read_csv(os.path.join(RES, "lost_originals.tsv"), sep="\t", keep_default_na=False)
    lost = lost[lost.experiment_unit != ""].merge(orig[["genome", "fasta"]], on="genome")
    rows = []
    for u in sorted({x for s in lost.experiment_unit for x in s.split(";")}):
        lu = lost[[u in s.split(";") for s in lost.experiment_unit]]
        a = production_attempt(u)
        stored = json.load(open(os.path.join(a, "fetch.json")))
        adir = os.path.join(EXP, "assemblies", u)
        os.makedirs(adir, exist_ok=True)
        fna = os.path.join(adir, f"{u}.assembly.fna")
        rep = json.loads(subprocess.run([PY, "-I", os.path.join(RA, "local/fetch_assembly.py"), stored["ref"], fna],
                                        check=True, capture_output=True, text=True).stdout)
        assert rep["download_md5"] == stored["download_md5"]
        seqs = read_fasta(fna)
        bdir = os.path.join(EXP, "runs", u, "raw_bins")
        os.makedirs(bdir, exist_ok=True)
        lt = {}
        for b in ("metabat2", "maxbin2", "concoct"):
            t = pd.read_csv(os.path.join(a, "c2b", f"{b}.contigs2bin.tsv"), sep="\t", header=None, names=["c", "bin"], dtype=str)
            t["c"] = t.c.str.split().str[0]
            for bn, g in t.groupby("bin"):
                with open(os.path.join(bdir, f"{bn}.fa"), "w") as fh:
                    for c in g.c:
                        fh.write(f">{c}\n{seqs[c]}\n")
                L = pd.Series([len(seqs[c]) for c in g.c])
                lt[bn] = (L[L < 2500].sum() / L.sum(), L.sum(), len(L))
        del seqs
        rl = os.path.join(bdir, "refs.txt")
        pd.Series(sorted(glob.glob(os.path.join(bdir, "*.fa")))).to_csv(rl, index=False, header=False)
        ql = os.path.join(bdir, "queries.txt")
        lu.fasta.to_csv(ql, index=False, header=False)
        out = os.path.join(bdir, "skani.tsv")
        subprocess.run(["nice", "-n", "10", SKANI, "dist", "--ql", ql, "--rl", rl, "-t", "16", "--min-af", "10", "-o", out],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        s = pd.read_csv(out, sep="\t")
        ev = pd.read_csv(os.path.join(a, "out", f"{u}_allBins.eval"), sep="\t").set_index("bin")
        f2g = dict(zip(lu.fasta, lu.genome))
        for r in s[s.ANI >= 95].itertuples():
            bn = os.path.basename(r.Ref_file)[:-3]
            e = ev.loc[bn] if bn in ev.index else None
            rows.append({"unit": u, "genome": f2g[r.Query_file], "raw_bin": bn, "binner": bn.split(".")[0], "ANI": r.ANI,
                         "af_of_original": r.Align_fraction_query, "af_of_raw_bin": r.Align_fraction_ref,
                         "raw_bin_mbp": round(lt[bn][1] / 1e6, 2), "raw_bin_contigs": lt[bn][2],
                         "raw_bin_frac_bp_lt2500": round(lt[bn][0], 3),
                         "das_completeness": e.SCG_completeness if e is not None else None,
                         "das_redundancy": e.SCG_redundancy if e is not None else None,
                         "das_score": round(e.bin_score, 3) if e is not None else None})
        shutil.rmtree(bdir)
        os.remove(fna)
        print(u, "done", len(lu), "lost originals", flush=True)
    t = pd.DataFrame(rows).merge(lost[["genome", "genus", "sample", "mimag", "lost30", "in_assembly_af"]], on="genome")
    t.sort_values(["genome", "af_of_original"], ascending=[True, False]).to_csv(os.path.join(RES, "lost_trace.tsv"),
                                                                              sep="\t", index=False)


if __name__ == "__main__":
    main()
