#!/usr/bin/env python3
"""Census update for run 4 (13,164 MAGs: the 12,859 of 2026-11-13 + 305 from the 11 re-assembled libraries).

1. Assemble the run-4 compute dir ~/gtdbtk_kbase/census_compute_run4 from
   - poplar wetland_census/run4 (run_census_compute_run4.sh): CheckM2 --allmodels for the 305, the skani triangle over
     all 13,164, and the brute-force searches of the 305 vs the three reference DBs;
   - the existing census inputs for the 12,859: CheckM2 (default + allmodels tables) and the brute-force re-scan.
   The 305's default-mode CheckM2 rows are derived from their --allmodels output by CheckM2's own choice
   (Completeness = Completeness_Specific if the Neural Network (Specific) model was used, else _General); for the 12,859
   this rule reproduces the default run exactly (asserted here).
   ASSERTED: the new triangle restricted to the 12,859 is identical to the census triangle (same pairs, ANI and AF).
2. Run census.py (inputs switched to run-4 paths) -> ~/gtdbtk_kbase/census_run4_out, then select_reps.py, then
   iterative_ids.py seeded with the current mag_iterativeIDs.tsv (existing IDs frozen; the 305 get the next numbers).
3. Compare with the current census for the 12,859: species-cluster partition (membership) changes vs label-only
   changes, strain clusters, iterativeIDs (must be identical), representative worklist; write run4_vs_run3.json.
Nothing in census/outputs is changed here; apply with --apply (copies the run-4 outputs over census/outputs and
~/gtdbtk_kbase/census_final, after the checks).
Usage: ~/Documents/py_venv/bin/python -I census_run4.py [--apply]
"""
import json, os, shutil, subprocess, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
G = os.path.expanduser("~/gtdbtk_kbase")
RAW = os.path.join(G, "census_compute_run4_raw")
CD4 = os.path.join(G, "census_compute_run4")
OUT4 = os.path.join(G, "census_run4_out")
OLD_CD, RESCAN, FINAL = os.path.join(G, "census_compute"), os.path.join(G, "census_rescan_noindex"), os.path.join(G, "census_final")
REPO = os.path.join(HERE, "outputs")
DBS = ("novel_all_representatives", "novel_hq_v2", "gtdb_r220_consolidated")
APPLY = "--apply" in sys.argv
PY = [sys.executable, "-I"]
P = "poplar:/scratch1/afreiburger/wetland_census/run4"
for d in (RAW, CD4, OUT4, os.path.join(CD4, "checkm2")):
    os.makedirs(d, exist_ok=True)
subprocess.run(["rsync", "-a", f"{P}/checkm2_new305/quality_report.tsv", f"{P}/skani_triangle.tsv",
                *[f"{P}/skani_vs_{n}.new305.tsv" for n in DBS], f"{P}/census_run4.events", RAW + "/"], check=True)
new305 = set(pd.read_csv(os.path.join(G, "run4_new_mags.tsv"), sep="\t").name)

# --- CheckM2 ---
def chosen(df):
    return df.Completeness_Specific.where(df.Completeness_Model_Used.str.startswith("Neural"), df.Completeness_General)

old_def = pd.read_csv(os.path.join(OLD_CD, "checkm2", "quality_report.tsv"), sep="\t")
old_all = pd.read_csv(os.path.join(OLD_CD, "checkm2", "quality_report_allmodels_merged.tsv"), sep="\t")
chk = old_all.set_index("Name").loc[old_def.Name]
assert (chosen(chk.reset_index()).round(2).values == old_def.Completeness.round(2).values).all(), "chosen-model rule"
new_all = pd.read_csv(os.path.join(RAW, "quality_report.tsv"), sep="\t")
assert set(new_all.Name) == new305 and len(new_all) == 305, "CheckM2 must cover exactly the 305"
new_def = new_all.assign(Completeness=chosen(new_all))[old_def.columns]
pd.concat([old_def, new_def]).to_csv(os.path.join(CD4, "checkm2", "quality_report.tsv"), sep="\t", index=False)
pd.concat([old_all, new_all[old_all.columns]]).to_csv(os.path.join(CD4, "checkm2", "quality_report_allmodels_merged.tsv"),
                                                      sep="\t", index=False)

# --- skani triangle: the old part must be unchanged ---
def edges(path, keep=None):
    t = pd.read_csv(path, sep="\t", usecols=["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"])
    a = t.Ref_file.map(lambda x: os.path.basename(x)[:-4]); b = t.Query_file.map(lambda x: os.path.basename(x)[:-4])
    lo, hi = a.where(a < b, b), b.where(a < b, a)
    afl = t.Align_fraction_ref.where(a < b, t.Align_fraction_query); afh = t.Align_fraction_query.where(a < b, t.Align_fraction_ref)
    e = pd.DataFrame({"a": lo, "b": hi, "ANI": t.ANI, "af_a": afl, "af_b": afh})
    if keep is not None:
        e = e[e.a.isin(keep) & e.b.isin(keep)]
    return e.sort_values(["a", "b"]).reset_index(drop=True)

old_mags = set(old_def.Name)
e_old = edges(os.path.join(OLD_CD, "skani_triangle.tsv"))
e_new = edges(os.path.join(RAW, "skani_triangle.tsv"), keep=old_mags)
same = len(e_old) == len(e_new) and (e_old.values == e_new.values).all()
print(f"triangle: old {len(e_old)} edges; new run restricted to the 12,859: {len(e_new)}; identical: {same}")
assert same, "the new triangle does not reproduce the census triangle for the 12,859"
shutil.copyfile(os.path.join(RAW, "skani_triangle.tsv"), os.path.join(CD4, "skani_triangle.tsv"))
for f in ("skani_binary.md5", "skani_version.txt"):
    shutil.copyfile(os.path.join(OLD_CD, f), os.path.join(CD4, f))

# --- reference searches: brute-force re-scan (12,859) + the 305 ---
for n in DBS:
    a = pd.read_csv(os.path.join(RESCAN, f"skani_vs_{n}.tsv"), sep="\t")
    b = pd.read_csv(os.path.join(RAW, f"skani_vs_{n}.new305.tsv"), sep="\t")
    qa = set(a.Query_file.map(os.path.basename)); qb = set(b.Query_file.map(os.path.basename))
    assert not qa & qb, n
    pd.concat([a, b]).to_csv(os.path.join(CD4, f"skani_vs_{n}.tsv"), sep="\t", index=False)

# --- census, reps, iterativeIDs ---
subprocess.run(PY + [os.path.join(HERE, "census.py"), CD4, OUT4], check=True, stdout=open(OUT4 + ".census.log", "w"),
               stderr=subprocess.STDOUT)
subprocess.run(PY + [os.path.join(HERE, "select_reps.py"), OUT4], check=True, stdout=open(OUT4 + ".reps.log", "w"),
               stderr=subprocess.STDOUT)
for f in ("mag_iterativeIDs.tsv", "dram2_hq_multisample.tsv"):
    shutil.copyfile(os.path.join(REPO, f), os.path.join(OUT4, f))
subprocess.run(PY + [os.path.join(HERE, "iterative_ids.py"), OUT4], check=True, stdout=open(OUT4 + ".ids.log", "w"),
               stderr=subprocess.STDOUT)

# --- compare with the current census for the 12,859 ---
o = pd.read_csv(os.path.join(REPO, "mag_census.tsv"), sep="\t").set_index("mag")
n4 = pd.read_csv(os.path.join(OUT4, "mag_census.tsv"), sep="\t").set_index("mag")
assert set(n4.index) == old_mags | new305 and len(n4) == 13164
oo, nn = o.loc[sorted(old_mags)], n4.loc[sorted(old_mags)]
label_changed = int((oo.species_cluster != nn.species_cluster).sum())
# partition change: two old MAGs together before and apart now, or vice versa
po = oo.groupby("species_cluster").groups; pn = nn.groupby("species_cluster").groups
def partition(gr):
    return {frozenset(v) for v in gr.values()}
before, after = partition(po), partition(pn)
after_old_only = {frozenset(s) for s in after}
moved = set()
for s in before ^ after_old_only:
    moved |= set(s)
idchg = int((oo.iterativeID != nn.iterativeID).sum())
reps_o = set(pd.read_csv(os.path.join(REPO, "rep_worklist.tsv"), sep="\t").mag)
reps_n = set(pd.read_csv(os.path.join(OUT4, "rep_worklist.tsv"), sep="\t").mag)
h_o = json.load(open(os.path.join(REPO, "headline.json")))
h_n = json.load(open(os.path.join(OUT4, "headline.json")))
new = n4.loc[sorted(new305)]
summary = {
    "mags": len(n4), "new_mags": len(new305),
    "old_mags_cluster_label_changed": label_changed,
    "old_mags_in_clusters_whose_old_membership_changed": len(moved),
    "old_mags_strain_label_changed": int((oo.strain_cluster != nn.strain_cluster).sum()),
    "old_mags_iterativeID_changed": idchg,
    "new_mags_joining_existing_clusters": int(new.species_cluster.isin(set(oo.species_cluster)).sum()),
    "new_mags_in_new_clusters": int((~new.species_cluster.isin(set(oo.species_cluster))).sum()),
    "new_mags_by_source": new.cluster_source.value_counts().to_dict(),
    "clusters_before_after": [int(o.species_cluster.nunique()), int(n4.species_cluster.nunique())],
    "rep_worklist_before_after": [len(reps_o), len(reps_n)], "reps_added": len(reps_n - reps_o),
    "reps_dropped": sorted(reps_o - reps_n),
    "headline_changed_keys": sorted(k for k in set(h_o) | set(h_n) if h_o.get(k) != h_n.get(k)),
}
json.dump(summary, open(os.path.join(OUT4, "run4_vs_run3.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
assert idchg == 0, "existing iterativeIDs changed"
if APPLY:
    for f in sorted(os.listdir(OUT4)):
        if f.endswith((".tsv", ".json")):
            for dst in (REPO, FINAL):
                shutil.copyfile(os.path.join(OUT4, f), os.path.join(dst, f))
    print("applied to", REPO, "and", FINAL)
else:
    print("dry run: outputs in", OUT4, "; re-run with --apply to copy them over census/outputs")
