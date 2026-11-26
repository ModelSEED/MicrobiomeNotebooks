#!/usr/bin/env python
"""READ-ONLY w.r.t. ~/Documents/EmilyKin. Builds the authoritative R220 -> R232 maps for the 291 EmilyKin MAGs.

Outputs (rename_plan/data/):
  mag_id_map_r220_r232.tsv        all 291: bin, old (committed) ID, new ID, case old/new, flags (ASV collision, chains)
  mag_id_map_r220_r232.json       {bin: new_id} for all 291 (drop-in replacement for scripts/_new_ids.json 'mag_new')
  mag_id_changes_r220_r232.csv    the changed MAGs only (bin, old, new, cases)
  mag_lineage_r220_r232.tsv       bin, R220 / R232 GTDB lineage, per-rank R220 / R232 names, changed ranks,
                                  agreement of the repo's four R220 taxonomy copies
  taxon_name_changes_r220_r232.tsv rank, R220 name -> R232 name over the 291 MAGs, n MAGs, whether the R220 name
                                  survives elsewhere in the R232 MAG taxonomy, and whether it is also an ASV
                                  (MiDAS/SILVA) taxon name in taxonomy.json (=> plain text replacement unsafe)
  build_maps.log
Verifies: simulation baseline reproduces the committed map 291/291 (re-run here, not just read from the log),
the committed map equals mag_abundance_summary.{json,tsv} / intersection CSV, and the new map is collision free.
"""
import collections, csv, json, os, re, sys, importlib.util
sys.dont_write_bytecode = True

EK = os.path.expanduser("~/Documents/EmilyKin")
PREP = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
OUT = os.path.join(PREP, "rename_plan", "data")
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
log_lines = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s); log_lines.append(s)

def lin2d(lin):
    parts = lin.split(";")
    return {R: (parts[i][3:] if i < len(parts) else "") for i, R in enumerate(RANKS)}

# --- re-run the simulation (read-only on the repo) ------------------------------------------------------------
spec = importlib.util.spec_from_file_location("sim", os.path.join(PREP, "simulate_iterativeids_r232.py"))
sim = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim)
taxonomy = json.load(open(os.path.join(EK, "taxonomy.json")))
ids = json.load(open(os.path.join(EK, "iterativeIDs.json")))           # ASV iterativeID -> hash
hash2id = {h: i for i, h in ids.items()}   # hash -> ASV iterativeID
asv_ids = set(ids.keys())
relab = json.load(open(os.path.join(EK, "relative_abundance.json")))
mean_ab = {h: (sum(d.values()) / len(d) if d else 0.0) for h, d in relab.items()}
a2m = json.load(open(os.path.join(EK, "asv_to_mag_mapping.json")))["asv_to_mag"]
mags = json.load(open(os.path.join(EK, "mag_abundance_summary.json")))
mapped = collections.defaultdict(list)
for h, rec in a2m.items():
    if rec.get("confidence_tier") in {"species", "genus", "family"} and rec.get("best_mag"):
        mapped[rec["best_mag"]].append(h)
order = [r[0] for r in csv.reader(open(os.path.join(EK, "mag_abundance_by_day_intersection_with_ids.csv")))][1:]
order = [m for m in order if m]
committed_json = json.load(open(os.path.join(EK, "mag_iterativeID_old_to_new.json")))
committed_csv = {r["MAG"]: r for r in csv.DictReader(open(os.path.join(EK, "mag_iterativeID_old_to_new.csv")))}
inter = {r["MAG"]: r["iterativeID"] for r in csv.DictReader(open(os.path.join(EK, "mag_abundance_by_day_intersection_with_ids.csv")))}
mas_tsv = {r["MAG"]: r for r in csv.DictReader(open(os.path.join(EK, "mag_abundance_summary.tsv")), delimiter="\t")}
mqs = {r["MAG"]: r for r in csv.DictReader(open(os.path.join(EK, "mag_quality_summary.tsv")), delimiter="\t")}
g220 = {r["user_genome"]: r for r in csv.DictReader(open(os.path.join(EK, "gtdbtk.bac120.summary.tsv")), delimiter="\t")}
g232 = {r["user_genome"]: r for r in csv.DictReader(open(os.path.join(PREP, "results", "gtdbtk.r232.summary.tsv")), delimiter="\t")}

base, base_case = sim.name_ids(order, lambda m: mags.get(m, {}).get("taxonomy", {}), taxonomy, mapped, mean_ab, hash2id)
new, new_case = sim.name_ids(order, lambda m: sim.lineage_to_dict(g232[m]["classification"]), taxonomy, mapped, mean_ab, hash2id)
log(f"MAGs in order: {len(order)}; R220 summary {len(g220)}; R232 summary {len(g232)}")
log(f"baseline (R220) reproduces committed mag_iterativeID_old_to_new.json: {sum(base[m]==committed_json.get(m) for m in order)}/{len(order)}")
log(f"committed json == csv new_mag_iterativeID: {sum(committed_json[m]==committed_csv[m]['new_mag_iterativeID'] for m in order)}/{len(order)}")
log(f"committed == intersection CSV iterativeID: {sum(committed_json[m]==inter[m] for m in order)}/{len(order)}")
log(f"committed == mag_abundance_summary.tsv mag_iterativeID: {sum(committed_json[m]==mas_tsv[m]['mag_iterativeID'] for m in order)}/{len(order)}")
log(f"committed == mag_abundance_summary.json mag_iterativeID: {sum(committed_json[m]==mags[m]['mag_iterativeID'] for m in order)}/{len(order)}")
log(f"baseline case == committed csv assignment_case: {sum(base_case[m]==committed_csv[m]['assignment_case'] for m in order)}/{len(order)}")
simtsv = {r["MAG"]: r for r in csv.DictReader(open(os.path.join(PREP, "results", "iterativeid_simulation_r232.tsv")), delimiter="\t")}
log(f"re-run == results/iterativeid_simulation_r232.tsv (r232 id): {sum(new[m]==simtsv[m]['simulated_r232_id'] for m in order)}/{len(order)}")
changed = [m for m in order if new[m] != base[m]]
log(f"changed MAG iterativeIDs: {len(changed)}; root changed: {sum(base[m].rsplit('.',1)[0]!=new[m].rsplit('.',1)[0] for m in changed)}")
dups = {v: c for v, c in collections.Counter(new.values()).items() if c > 1}
log(f"duplicate new IDs: {dups}")

old_owner = {base[m]: m for m in order}           # old id -> bin
new_owner = {new[m]: m for m in order}
rows = []
for m in order:
    o, n = base[m], new[m]
    chain_in = old_owner.get(n) if (n in old_owner and old_owner[n] != m) else ""   # new id was the old id of another MAG
    chain_out = new_owner.get(o) if (o in new_owner and new_owner[o] != m) else ""   # old id becomes the new id of another MAG
    rows.append(dict(MAG=m, old_id=o, new_id=n, case_old=base_case[m], case_new=new_case[m], changed=n != o,
                     root_changed=o.rsplit(".", 1)[0] != n.rsplit(".", 1)[0],
                     old_is_asv_id=o in asv_ids, new_is_asv_id=n in asv_ids,
                     new_id_was_old_id_of=chain_in, old_id_becomes_new_id_of=chain_out,
                     old_display=o.replace("Ca_", "Ca. ").replace("_m", "") if o.startswith("Ca_") else "",
                     r220_genus=lin2d(g220[m]["classification"])["Genus"], r232_genus=lin2d(g232[m]["classification"])["Genus"]))
with open(os.path.join(OUT, "mag_id_map_r220_r232.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
json.dump({m: new[m] for m in order}, open(os.path.join(OUT, "mag_id_map_r220_r232.json"), "w"), indent=0)
with open(os.path.join(OUT, "mag_id_changes_r220_r232.csv"), "w") as fh:
    w = csv.writer(fh, lineterminator="\n"); w.writerow(["MAG", "old_mag_iterativeID", "new_mag_iterativeID", "case_r220", "case_r232",
                                                         "old_is_asv_id", "new_is_asv_id", "chain_note"])
    for r in rows:
        if r["changed"]:
            note = []
            if r["new_id_was_old_id_of"]: note.append(f"new id was old id of {r['new_id_was_old_id_of']}")
            if r["old_id_becomes_new_id_of"]: note.append(f"old id becomes new id of {r['old_id_becomes_new_id_of']}")
            w.writerow([r["MAG"], r["old_id"], r["new_id"], r["case_old"], r["case_new"], r["old_is_asv_id"], r["new_is_asv_id"], "; ".join(note)])
log("old IDs that are also current ASV iterativeIDs (ASV keeps it):", sum(r["old_is_asv_id"] for r in rows if r["changed"]))
log("new IDs that are also current ASV iterativeIDs:", [r["new_id"] for r in rows if r["changed"] and r["new_is_asv_id"]])
log("chains:", [(r["MAG"], r["old_id"], r["new_id"], r["new_id_was_old_id_of"], r["old_id_becomes_new_id_of"]) for r in rows
               if r["changed"] and (r["new_id_was_old_id_of"] or r["old_id_becomes_new_id_of"])])
# unchanged MAG IDs equal to a changed MAG's new id? (would be a duplicate; checked above)
# substring hazards: old id that is a prefix of another current id (MAG or ASV) at an ID boundary
allcur = asv_ids | set(base.values())
haz = collections.defaultdict(list)
for r in rows:
    if not r["changed"]: continue
    o = r["old_id"]
    for x in allcur:
        if x != o and x.startswith(o):
            haz[o].append(x)
log("old IDs that are string prefixes of other current IDs (substring hazards):",
    {k: sorted(v)[:8] + (["..."] if len(v) > 8 else []) for k, v in haz.items()})
json.dump({k: sorted(v) for k, v in haz.items()}, open(os.path.join(OUT, "substring_hazards.json"), "w"), indent=1)

# --- lineage map --------------------------------------------------------------------------------------------
lrows = []
cons = collections.Counter()
for m in order:
    a, b = lin2d(g220[m]["classification"]), lin2d(g232[m]["classification"])
    tj = mags[m]["taxonomy"]; tt = mas_tsv[m]; tq = mqs.get(m, {})
    ok_json = all((tj.get(R, "") or "") == a[R] for R in RANKS[:6])
    ok_tsv = all((tt.get(R, "") or "") == a[R] for R in RANKS)
    ok_q = all((tq.get(R, "") or "") == a[R] for R in RANKS)
    cons["json_eq_summary"] += ok_json; cons["tsv_eq_summary"] += ok_tsv; cons["quality_eq_summary"] += ok_q
    ch = [R for R in RANKS if a[R] != b[R]]
    d = dict(MAG=m, r220_classification=g220[m]["classification"], r232_classification=g232[m]["classification"],
             changed_ranks=",".join(ch), first_changed_rank=ch[0] if ch else "",
             r232_method=g232[m]["classification_method"], r232_closest_ref=g232[m]["closest_genome_reference"],
             r232_ani=g232[m]["closest_genome_ani"], r232_af=g232[m]["closest_genome_af"], r232_red=g232[m]["red_value"],
             r232_msa_percent=g232[m]["msa_percent"], r232_translation_table=g232[m]["translation_table"],
             r232_warnings=g232[m]["warnings"], repo_json_matches_r220=ok_json, repo_tsv_matches_r220=ok_tsv,
             repo_quality_matches_r220=ok_q)
    for R in RANKS:
        d[f"r220_{R}"] = a[R]; d[f"r232_{R}"] = b[R]
    lrows.append(d)
log("R220 taxonomy copies agree with repo gtdbtk.bac120.summary.tsv:", dict(cons))
with open(os.path.join(OUT, "mag_lineage_r220_r232.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(lrows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(lrows)
log("MAGs with any lineage change:", sum(1 for r in lrows if r["changed_ranks"]),
    "| first changed rank:", dict(collections.Counter(r["first_changed_rank"] or "identical" for r in lrows)))

# --- taxon-name changes ----------------------------------------------------------------------------------------
asv_names = collections.defaultdict(set)   # name -> ASV ranks it occurs at
for h, t in taxonomy.items():
    for k, v in t.items():
        if v: asv_names[str(v)].add(k)
r232_names = {R: collections.Counter(r[f"r232_{R}"] for r in lrows) for R in RANKS}
r220_names = {R: collections.Counter(r[f"r220_{R}"] for r in lrows) for R in RANKS}
tc = collections.Counter()
for r in lrows:
    for R in RANKS:
        if r[f"r220_{R}"] != r[f"r232_{R}"]:
            tc[(R, r[f"r220_{R}"], r[f"r232_{R}"])] += 1
with open(os.path.join(OUT, "taxon_name_changes_r220_r232.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["rank", "r220_name", "r232_name", "n_mags", "r220_name_still_used_in_r232_mag_taxonomy",
                "r220_name_n_mags_r220", "r220_name_is_asv_taxon_name", "asv_ranks", "r232_name_is_asv_taxon_name",
                "r220_name_in_any_mag_id_old"])
    oldroots = collections.Counter(v.rsplit(".", 1)[0].replace("Ca_", "") for v in base.values())
    for (R, a, b), n in sorted(tc.items(), key=lambda x: (RANKS.index(x[0][0]), -x[1])):
        still = any(r232_names[RR].get(a, 0) for RR in RANKS) if a else False
        w.writerow([R, a, b, n, still, r220_names[R].get(a, 0), bool(a and a in asv_names),
                    ",".join(sorted(asv_names.get(a, []))), bool(b and b in asv_names), oldroots.get(a, 0)])
log("distinct (rank, r220->r232) name changes:", len(tc), "| species-rank:", sum(1 for k in tc if k[0] == "Species"))
open(os.path.join(OUT, "build_maps.log"), "w").write("\n".join(log_lines) + "\n")
