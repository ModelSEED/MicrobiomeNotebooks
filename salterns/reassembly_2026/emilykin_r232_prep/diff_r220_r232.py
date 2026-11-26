#!/usr/bin/env python
"""Per-MAG diff of the EmilyKin GTDB R220 classification against the new GTDB-Tk 2.8.0 / R232 run.

Inputs:  --r220     the R220 summary in use (default ~/Documents/EmilyKin/gtdbtk.bac120.summary.tsv; NOTE it is a hybrid:
                    for the 134 best-of-both reverted MAGs only the 'classification' column is theirs)
         --r232     one or more R232 summary files (gtdbtk.bac120.summary.tsv and, if present, gtdbtk.ar53.summary.tsv)
         --projection  tables/emilykin_r232_projection.tsv (the desk prediction; scored against the real run)
Outputs (--out, default results/):
         gtdbtk.r232.summary.tsv     bac120 + ar53 rows in one file (input for simulate_iterativeids_r232.py)
         r220_vs_r232_per_mag.tsv    one row per MAG: both lineages, method/ANI/AF/ref on both sides, change class
         r220_vs_r232_summary.tsv    counts per change class, species gained/lost/renamed, prediction accuracy
Change classes (first rank, domain->species, whose name differs): identical; <rank>_suffix_only (e.g. g__Luteimonas_D ->
g__Luteimonas); <rank>_placeholder_to_latin (g__JAAKGP01 -> g__Pseudogemmobacter); <rank>_renamed; <rank>_gained /
<rank>_lost (an empty rank on one side, e.g. species_gained = R220 's__' -> named R232 species).
Usage: python -I diff_r220_r232.py --r232 <out>/gtdbtk.bac120.summary.tsv [<out>/gtdbtk.ar53.summary.tsv]
"""
import argparse, collections, csv, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def strip_suffix(x):
    return " ".join(re.sub(r"_[A-Z]+$", "", t) for t in x.split(" "))


def placeholder(name):
    core = name[3:]
    if name.startswith("s__"):
        return bool(re.search(r" sp\d+", core)) or bool(re.search(r"\d", core.split(" ")[0]))
    return bool(re.search(r"\d", core))


def lineage_change(a, b):
    a = (a + [f"{r[0]}__" for r in RANKS])[:7]
    b = (b + [f"{r[0]}__" for r in RANKS])[:7]
    for i, (x, y) in enumerate(zip(a, b)):
        if x == y:
            continue
        if x.endswith("__") and not y.endswith("__"):
            return RANKS[i], "gained"
        if y.endswith("__") and not x.endswith("__"):
            return RANKS[i], "lost"
        if strip_suffix(x) == strip_suffix(y):
            return RANKS[i], "suffix_only"
        if placeholder(x) and not placeholder(y):
            return RANKS[i], "placeholder_to_latin"
        return RANKS[i], "renamed"
    return "none", "identical"


def read_summary(paths):
    rows = {}
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[r["user_genome"]] = r
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r220", default=os.path.expanduser("~/Documents/EmilyKin/gtdbtk.bac120.summary.tsv"))
    ap.add_argument("--r232", nargs="+", required=True)
    ap.add_argument("--projection", default=os.path.join(HERE, "tables", "emilykin_r232_projection.tsv"))
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    old, new = read_summary([a.r220]), read_summary(a.r232)
    proj = {}
    if a.projection and os.path.exists(a.projection):
        proj = {r["mag"]: r for r in csv.DictReader(open(a.projection), delimiter="\t")}
    missing = sorted(set(old) - set(new))
    extra = sorted(set(new) - set(old))

    # one combined R232 summary (bac120 columns; ar53 has the same header)
    cols = None
    with open(os.path.join(a.out, "gtdbtk.r232.summary.tsv"), "w") as fh:
        for p in a.r232:
            with open(p) as fi:
                rd = csv.DictReader(fi, delimiter="\t")
                if cols is None:
                    cols = rd.fieldnames
                    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
                    w.writeheader()
                for r in rd:
                    w.writerow(r)

    cnt = collections.Counter()
    out_rows = []
    for m in sorted(old):
        o = old[m]
        n = new.get(m)
        lo = o["classification"].split(";")
        ln = n["classification"].split(";") if n else []
        if not n:
            cls = "missing_in_r232"
        elif n["classification"].startswith("Unclassified"):
            cls = "unclassified_in_r232"
        else:
            rank, kind = lineage_change(lo, ln)
            cls = "identical" if kind == "identical" else f"{rank}_{kind}"
        cnt[f"change_class={cls}"] += 1
        sp_o = len(lo) == 7 and not lo[6].endswith("__")
        sp_n = len(ln) == 7 and not ln[6].endswith("__")
        sp = ("kept_same" if sp_o and sp_n and lo[6] == ln[6] else "kept_renamed" if sp_o and sp_n else
              "lost" if sp_o else "gained" if sp_n else "none_either")
        cnt[f"species={sp}"] += 1
        g_same = len(ln) >= 6 and lo[:6] == ln[:6]
        cnt[f"same_lineage_to_genus={g_same}"] += 1
        if len(ln) >= 6:
            cnt[f"same_phylum_and_genus_names={lo[1] == ln[1] and lo[5] == ln[5]}"] += 1
        row = {"mag": m, "r220_classification": o["classification"], "r232_classification": n["classification"] if n else "",
               "change_class": cls, "species_status": sp,
               "r220_method": o.get("classification_method", ""), "r232_method": n.get("classification_method", "") if n else "",
               "r220_ref": o.get("closest_genome_reference", ""), "r232_ref": n.get("closest_genome_reference", "") if n else "",
               "r232_ani": n.get("closest_genome_ani", "") if n else "", "r232_af": n.get("closest_genome_af", "") if n else "",
               "r232_red": n.get("red_value", "") if n else "", "r232_warnings": n.get("warnings", "") if n else ""}
        if m in proj and n:
            p = proj[m]["projected_r232_classification"].split(";")
            row["predicted_r232"] = proj[m]["projected_r232_classification"]
            row["prediction_basis"] = proj[m]["projection_basis"]
            row["prediction_exact"] = p == ln
            row["prediction_species_ok"] = (p[6] == ln[6]) if len(ln) == 7 else ""
            row["prediction_to_genus_ok"] = p[:6] == ln[:6]
            cnt[f"prediction[{proj[m]['projection_basis']}] species_ok={row['prediction_species_ok']}"] += 1
            cnt[f"prediction[{proj[m]['projection_basis']}] to_genus_ok={row['prediction_to_genus_ok']}"] += 1
        out_rows.append(row)
    fields = list(dict.fromkeys(k for r in out_rows for k in r))
    with open(os.path.join(a.out, "r220_vs_r232_per_mag.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    with open(os.path.join(a.out, "r220_vs_r232_summary.tsv"), "w") as fh:
        fh.write("metric\tn_mags\n")
        fh.write(f"mags_r220\t{len(old)}\nmags_r232\t{len(new)}\nmissing_in_r232\t{len(missing)}\nextra_in_r232\t{len(extra)}\n")
        for k, v in sorted(cnt.items()):
            fh.write(f"{k}\t{v}\n")
    print(open(os.path.join(a.out, "r220_vs_r232_summary.tsv")).read())
    if missing or extra:
        print("WARNING: genome sets differ:", missing[:5], extra[:5])


if __name__ == "__main__":
    sys.exit(main())
