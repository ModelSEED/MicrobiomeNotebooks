#!/usr/bin/env python3
"""Pull the KBDLSKANI sweep's per-batch hit tables from poplar and merge them (local; run any time, idempotent).

Writes, under kbdl_build/skani_sweep/:
  kbdl_skani_hits.tsv.gz   mag, iterativeID, archive, db, ref_genome, ani (0-1 fraction; distance = 1 - ani). Every MAG
                           swept has >= 1 row per DB; a row with empty ref_genome/ani means "searched, no hit".
                           NO alignment fraction (KBDLSKANI drops it): ani >= 0.95 is NOT a species call.
  failed_mags.tsv          MAGs whose search failed twice (sweep ledger table failed_mags)
  summary.json             coverage, per-DB no-hit counts, and agreement with the census's local skani search vs
                           novel_all_representatives (census used -n 5 per query over one 12,859-query run)
Asserts: no MAG x DB pair comes from two batch files; every MAG is a census MAG.
Usage: ~/Documents/py_venv/bin/python -I collect_skani_sweep.py
"""
import glob, gzip, json, os, sqlite3, subprocess, tempfile
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "skani_sweep")
CEN = os.path.join(HERE, "..", "census", "outputs")
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

# best hit per MAG x DB; agreement with the census's local novel_all_representatives search
best = h.dropna(subset=["ani"]).groupby(["mag", "db"]).ani.max().unstack()
swept = sorted(set(h.mag))
nov = best.get("novel_all_representatives", pd.Series(dtype=float)).reindex(swept)
cnov = cen.loc[swept, "novel_all_representatives_best_ani"] / 100
both = nov.notna() & cnov.notna()
summary = {
    "batch_files": len(files), "batches_by_state": batches, "mags_swept": len(swept), "of_census_mags": len(cen),
    "failed_mags": len(failed), "rows": len(h),
    "mags_with_hit": {d: int(best[d].notna().sum()) if d in best else 0 for d in sorted(set(h.db))},
    "mags_without_hit": {d: int(len(swept) - (best[d].notna().sum() if d in best else 0)) for d in sorted(set(h.db))},
    "novel_vs_census": {
        "both_have_hit": int(both.sum()),
        "best_ani_equal_within_1e-4": int(((nov - cnov).abs() < 1e-4)[both].sum()),
        "kbdl_higher": int(((nov - cnov) >= 1e-4)[both].sum()), "census_higher": int(((cnov - nov) >= 1e-4)[both].sum()),
        "kbdl_hit_census_none": int((nov.notna() & cnov.isna()).sum()),
        "census_hit_kbdl_none": int((nov.isna() & cnov.notna()).sum()),
    },
}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
