#!/usr/bin/env python3
"""Pull the KBDLSKANI sweep's per-batch hit tables from poplar and merge them (local; run any time, idempotent).

Writes, under kbdl_build/skani_sweep/:
  kbdl_skani_hits.tsv.gz   mag, iterativeID, archive, db, ref_genome, ani (0-1 fraction; distance = 1 - ani). Every MAG
                           swept has >= 1 row per DB; a row with empty ref_genome/ani means "searched, no hit".
                           NO alignment fraction (KBDLSKANI drops it): ani >= 0.95 is NOT a species call.
  failed_mags.tsv          MAGs whose search failed twice (sweep ledger table failed_mags)
  summary.json             coverage, per-DB no-hit counts, and a PAIR-LEVEL check of every (MAG, reference) novel hit
                           against the census's brute-force re-scan (run_census_rescan.sh, --no-marker-index, no -n):
                           both must hold the same pairs with the same ANI, which is what shows the sweep's jobs
                           were not screened by skani's lossy marker index
Asserts: no MAG x DB pair comes from two batch files; every MAG is a census MAG.
Usage: ~/Documents/py_venv/bin/python -I collect_skani_sweep.py
"""
import glob, gzip, json, os, sqlite3, subprocess, tempfile
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "skani_sweep")
CEN = os.path.join(HERE, "..", "census", "outputs")
RESCAN_NOVEL = os.path.expanduser("~/gtdbtk_kbase/census_rescan_noindex/skani_vs_novel_all_representatives.tsv")
REMOTE = "poplar:/scratch/afreiburger/wetland_kbdl_build"
os.makedirs(OUT, exist_ok=True)
stage = tempfile.mkdtemp(prefix="skani_sweep_")
subprocess.run(["rsync", "-a", "--include=batch_*.tsv", "--exclude=*", f"{REMOTE}/skani_hits/", stage + "/"], check=True)
subprocess.run(["rsync", "-a", f"{REMOTE}/skani_ledger.sqlite", stage + "/"], check=True)

files = sorted(glob.glob(os.path.join(stage, "batch_*.tsv")))
parts = []
for f in files:
    d = pd.read_csv(f, sep="\t", dtype={"ref_genome": str})
    d["batch"] = os.path.basename(f)
    parts.append(d)
h = pd.concat(parts, ignore_index=True)
src = h.groupby(["mag", "db"]).batch.nunique()
assert (src == 1).all(), f"MAG x DB pairs from > 1 batch: {src[src > 1].index[:5].tolist()}"
ids = pd.read_csv(os.path.join(CEN, "mag_iterativeIDs.tsv"), sep="\t").set_index("mag").iterativeID
cen = pd.read_csv(os.path.join(CEN, "mag_census.tsv"), sep="\t").set_index("mag")
assert set(h.mag) <= set(cen.index), "non-census MAG in the sweep"
h.insert(1, "iterativeID", h.mag.map(ids))
h = h.drop(columns="batch").sort_values(["mag", "db", "ani"], ascending=[True, True, False])
with gzip.open(os.path.join(OUT, "kbdl_skani_hits.tsv.gz"), "wt", compresslevel=9) as fh:
    h.to_csv(fh, sep="\t", index=False, float_format="%.6f")

led = sqlite3.connect(os.path.join(stage, "skani_ledger.sqlite"))
failed = pd.read_sql("SELECT * FROM failed_mags", led) if led.execute(
    "SELECT 1 FROM sqlite_master WHERE name='failed_mags'").fetchone() else pd.DataFrame(columns=["mag"])
failed.to_csv(os.path.join(OUT, "failed_mags.tsv"), sep="\t", index=False)
batches = dict(led.execute("SELECT state, count(*) FROM batches GROUP BY state").fetchall())

# best hit per MAG x DB; pair-level agreement with the brute-force census re-scan vs novel_all_representatives
best = h.dropna(subset=["ani"]).groupby(["mag", "db"]).ani.max().unstack()
swept = sorted(set(h.mag))
r = pd.read_csv(RESCAN_NOVEL, sep="\t", usecols=["Ref_file", "Query_file", "ANI"])
r["mag"] = r.Query_file.map(lambda x: os.path.basename(x)[:-4])
r["ref"] = r.Ref_file.map(lambda x: os.path.splitext(os.path.basename(x))[0])
r = r[r.mag.isin(swept)]
rp = dict(zip(zip(r.mag, r.ref), r.ANI / 100))
n = h[(h.db == "novel_all_representatives") & h.ref_genome.notna()]
kp = dict(zip(zip(n.mag, n.ref_genome), n.ani))
shared = set(kp) & set(rp)
summary = {
    "batch_files": len(files), "batches_by_state": batches, "mags_swept": len(swept), "of_census_mags": len(cen),
    "failed_mags": len(failed), "rows": len(h),
    "mags_with_hit": {d: int(best[d].notna().sum()) if d in best else 0 for d in sorted(set(h.db))},
    "mags_without_hit": {d: int(len(swept) - (best[d].notna().sum() if d in best else 0)) for d in sorted(set(h.db))},
    "novel_pairs_vs_bruteforce_census": {
        "kbdl_pairs": len(kp), "census_pairs": len(rp), "only_kbdl": len(set(kp) - set(rp)),
        "only_census": len(set(rp) - set(kp)),
        "max_abs_ani_diff": float(max((abs(kp[k] - rp[k]) for k in shared), default=0.0)),
    },
}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
