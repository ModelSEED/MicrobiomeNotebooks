#!/usr/bin/env python3
"""Fold the brute-force re-scan of the census's reference searches into the census outputs (2026-11-17).

The census's skani searches vs novel_all_representatives, novel_hq_v2 and gtdb_r220_consolidated were single
12,859-query calls, so skani 0.3.1 used its marker-hash index, which loses hits (run_census_rescan.sh explains and
measures it). run_census_rescan.sh repeated them with --no-marker-index and no -n cap. This script:
  1. builds a compute dir = the original one (triangle, CheckM2, versions) with the three re-scanned search tables;
  2. re-runs census.py into a scratch output dir;
  3. ASSERTS that only reference-search fields changed: in mag_census.tsv the <db>_species_hit/_species_ani/_best_ani/
     _best_af columns, in clusters.tsv the *_in_novel_* columns, in headline.json the <db>_* and r220_vs_r232 keys;
     every other column/key (clusters, strains, quality, R232, adjacency) must be identical;
  4. with --apply, copies those columns/keys into census/outputs (keeping the iterativeID column and everything else
     there) and into ~/gtdbtk_kbase/census_final, and prints the before/after numbers for RESULTS.md.
Usage: ~/Documents/py_venv/bin/python -I apply_rescan.py <dir with re-scanned skani_vs_*.tsv> [--apply]
"""
import json, os, shutil, subprocess, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESCAN = sys.argv[1]
APPLY = "--apply" in sys.argv
ORIG_CD = os.path.expanduser("~/gtdbtk_kbase/census_compute")
CD = os.path.expanduser("~/gtdbtk_kbase/census_compute_rescan")
SCRATCH = os.path.expanduser("~/gtdbtk_kbase/census_rescan_out")
FINAL = os.path.expanduser("~/gtdbtk_kbase/census_final")
REPO = os.path.join(HERE, "outputs")
DBS = ("novel_all_representatives", "novel_hq_v2", "gtdb_r220_consolidated")

os.makedirs(CD, exist_ok=True)
for f in ("skani_triangle.tsv", "checkm2", "skani_binary.md5", "skani_version.txt"):
    dst = os.path.join(CD, f)
    if not os.path.lexists(dst):
        os.symlink(os.path.join(ORIG_CD, f), dst)
for n in DBS:
    src = os.path.join(RESCAN, f"skani_vs_{n}.tsv")
    assert os.path.getsize(src) > 0, src
    shutil.copyfile(src, os.path.join(CD, f"skani_vs_{n}.tsv"))
subprocess.run([sys.executable, "-I", os.path.join(HERE, "census.py"), CD, SCRATCH], check=True,
               stdout=open(SCRATCH + ".log", "w"), stderr=subprocess.STDOUT)


def ref_col(c):
    return any(c.startswith(n + "_") for n in DBS) or "_in_novel_" in c


changes = {}
for fname, key in (("mag_census.tsv", "mag"), ("clusters.tsv", "species_cluster")):
    old = pd.read_csv(os.path.join(FINAL, fname), sep="\t", dtype=str, keep_default_na=False).set_index(key)
    new = pd.read_csv(os.path.join(SCRATCH, fname), sep="\t", dtype=str, keep_default_na=False).set_index(key)
    assert set(old.index) == set(new.index), f"{fname}: row sets differ"
    new = new.loc[old.index]
    extra = set(old.columns) - set(new.columns)          # e.g. iterativeID, added downstream
    assert set(new.columns) <= set(old.columns), f"{fname}: new columns {set(new.columns) - set(old.columns)}"
    for c in new.columns:
        diff = int((old[c] != new[c]).sum())
        if diff and not ref_col(c):
            raise SystemExit(f"{fname}: non-reference column {c} changed in {diff} rows; not applying")
        if diff:
            changes[f"{fname}:{c}"] = diff
    print(f"{fname}: columns kept as is (added downstream): {sorted(extra)}")
    if APPLY:
        for target in (os.path.join(REPO, fname), os.path.join(FINAL, fname)):
            t = pd.read_csv(target, sep="\t", dtype=str, keep_default_na=False).set_index(key)
            assert set(t.index) == set(new.index)
            for c in new.columns:
                if ref_col(c):
                    t[c] = new.loc[t.index, c]
            t.reset_index().to_csv(target, sep="\t", index=False)

ho = json.load(open(os.path.join(FINAL, "headline.json")))
hn = json.load(open(os.path.join(SCRATCH, "headline.json")))
for k in set(ho) | set(hn):
    if ho.get(k) != hn.get(k):
        if not (any(k.startswith(n) for n in DBS) or k == "r220_vs_r232"):
            raise SystemExit(f"headline.json: non-reference key {k} changed: {ho.get(k)} -> {hn.get(k)}")
        changes[f"headline:{k}"] = {"before": ho.get(k), "after": hn.get(k)}
print(json.dumps(changes, indent=1))
if APPLY:
    for target in (os.path.join(REPO, "headline.json"), os.path.join(FINAL, "headline.json")):
        h = json.load(open(target))
        for k in hn:
            if any(k.startswith(n) for n in DBS) or k == "r220_vs_r232":
                h[k] = hn[k]
        json.dump(h, open(target, "w"), indent=1)
    print("applied to", REPO, "and", FINAL)
else:
    print("dry run; re-run with --apply")
