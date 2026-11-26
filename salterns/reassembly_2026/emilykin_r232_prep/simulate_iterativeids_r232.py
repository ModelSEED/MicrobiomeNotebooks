#!/usr/bin/env python
"""READ-ONLY simulation: which EmilyKin MAG iterativeIDs would change if the GTDB taxonomy moved from R220 to R232?

Replicates the current MAG naming rule of ~/Documents/EmilyKin/scripts/recompute_mag_ids.py (commit 408cb39 map:
agreeing mapped non-weak ASV first -> inherit its id; else GTDB genus -> '{Genus}.{N}_m'; else deepest non-numeric
GTDB label -> '{label}.{N}_m'), first with the R220 taxonomy in mag_abundance_summary.json (must reproduce the
committed map 291/291), then with a substituted R232 lineage per MAG (from tables/emilykin_r232_projection.tsv or,
after the real run, from the R232 GTDB-Tk summary via --r232-summary).
Never writes into the EmilyKin repo (the real generators there are one-shot and unsafe to re-run).
Usage: python -I simulate_iterativeids_r232.py [--emilykin ~/Documents/EmilyKin] [--r232-summary gtdbtk.bac120.summary.tsv]
"""
import argparse, collections, csv, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
HASDIGIT = re.compile(r"\d")
GTDB = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]
G2A = {"Domain": "Kingdom", "Phylum": "Phylum", "Class": "Class", "Order": "Order", "Family": "Family", "Genus": "Genus"}


def numeric(s): return bool(HASDIGIT.search(str(s)))
def clean(v):
    s = str(v).strip(); return "" if s.lower() in ("", "nan", "none") else s
def norm(v):
    v = clean(v); return max(v.split("_"), key=len).lower() if v else ""


def name_ids(order, tax_of, taxonomy, mapped, mean_ab, hash2id):
    mag_new, mag_case, mc = {}, {}, collections.defaultdict(int)
    for m in order:
        t = tax_of(m)
        rr, rv = "", ""
        for L in GTDB:
            v = clean(t.get(L, ""))
            if v and not numeric(v):
                rr, rv = L, v
        agree = []
        if rr and rr != "Domain" and rv:
            for h in mapped.get(m, []):
                av = clean(taxonomy.get(h, {}).get(G2A[rr], ""))
                if av and norm(av) == norm(rv):
                    agree.append(h)
        if agree:
            best = max(agree, key=lambda h: (mean_ab.get(h, 0.0), hash2id.get(h, "")))
            mag_new[m], mag_case[m] = hash2id[best], "asv_inherit"; continue
        g = clean(t.get("Genus", ""))
        if g and not numeric(g):
            mc[g] += 1; mag_new[m], mag_case[m] = f"{g}.{mc[g]}_m", "gtdb_genus"
        else:
            lbl = rv or "Bacteria"; mc[lbl] += 1
            mag_new[m], mag_case[m] = f"{lbl}.{mc[lbl]}_m", "gtdb_lowest"
    return mag_new, mag_case


def lineage_to_dict(lin):
    parts = lin.split(";")
    return {L: (parts[i][3:] if i < len(parts) else "") for i, L in enumerate(GTDB)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emilykin", default=os.path.expanduser("~/Documents/EmilyKin"))
    ap.add_argument("--projection", default=os.path.join(HERE, "tables", "emilykin_r232_projection.tsv"))
    ap.add_argument("--r232-summary", default=None, help="real R232 gtdbtk.bac120.summary.tsv (overrides --projection)")
    ap.add_argument("--out", default=os.path.join(HERE, "tables", "iterativeid_simulation_r232.tsv"))
    a = ap.parse_args()
    R = a.emilykin
    taxonomy = json.load(open(os.path.join(R, "taxonomy.json")))
    ids = json.load(open(os.path.join(R, "iterativeIDs.json")))
    hash2id = {h: i for i, h in ids.items()}
    relab = json.load(open(os.path.join(R, "relative_abundance.json")))
    mean_ab = {h: (sum(d.values()) / len(d) if d else 0.0) for h, d in relab.items()}
    a2m = json.load(open(os.path.join(R, "asv_to_mag_mapping.json")))["asv_to_mag"]
    mags = json.load(open(os.path.join(R, "mag_abundance_summary.json")))
    mapped = collections.defaultdict(list)
    for h, rec in a2m.items():
        if rec.get("confidence_tier") in {"species", "genus", "family"} and rec.get("best_mag"):
            mapped[rec["best_mag"]].append(h)
    order = [r[0] for r in csv.reader(open(os.path.join(R, "mag_abundance_by_day_intersection_with_ids.csv")))][1:]
    order = [m for m in order if m]
    committed = {r["MAG"]: r["new_mag_iterativeID"] for r in csv.DictReader(open(os.path.join(R, "mag_iterativeID_old_to_new.csv")))}

    base, base_case = name_ids(order, lambda m: mags.get(m, {}).get("taxonomy", {}), taxonomy, mapped, mean_ab, hash2id)
    repro = sum(1 for m in order if base[m] == committed.get(m))
    src = a.r232_summary or a.projection
    col = "classification" if a.r232_summary else "projected_r232_classification"
    keycol = "user_genome" if a.r232_summary else "mag"
    lin232 = {r[keycol]: r[col] for r in csv.DictReader(open(src), delimiter="\t")}
    new, new_case = name_ids(order, lambda m: lineage_to_dict(lin232[m]) if m in lin232 else
                             mags.get(m, {}).get("taxonomy", {}), taxonomy, mapped, mean_ab, hash2id)
    ch = [m for m in order if new[m] != base[m]]
    with open(a.out, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["MAG", "committed_id", "simulated_r220_id", "case_r220", "simulated_r232_id", "case_r232", "changed",
                    "r220_genus", "r232_genus"])
        for m in order:
            w.writerow([m, committed.get(m, ""), base[m], base_case[m], new[m], new_case[m], new[m] != base[m],
                        clean(mags.get(m, {}).get("taxonomy", {}).get("Genus", "")),
                        lineage_to_dict(lin232[m])["Genus"] if m in lin232 else ""])
    print(f"baseline reproduces committed map: {repro}/{len(order)}")
    print(f"cases R220: {dict(collections.Counter(base_case.values()))}")
    print(f"cases R232: {dict(collections.Counter(new_case.values()))}")
    print(f"MAG iterativeIDs that change with {os.path.basename(src)}: {len(ch)}")
    by = collections.Counter((base_case[m], new_case[m]) for m in ch)
    print("  by case (r220 -> r232):", dict(by))
    print("  root changed (not just renumbered):", sum(1 for m in ch if base[m].rsplit('.', 1)[0] != new[m].rsplit('.', 1)[0]))
    print("  duplicates in simulated R232 map:", {v: c for v, c in collections.Counter(new.values()).items() if c > 1})
    for m in ch:
        print(f"   {m:16s} {base[m]:26s} -> {new[m]:26s} [{base_case[m]} -> {new_case[m]}]")


if __name__ == "__main__":
    sys.exit(main())
