#!/usr/bin/env python
"""Predict, MAG by MAG, what GTDB R232 will call the 291 EmilyKin MAGs (classified so far against GTDB R220).

Two independent lines of evidence, both cheap:
  1. Name mapping: each R220 taxon of a MAG -> the R232 name that the same reference genomes carry in R232
     (dominant share over all genomes of that taxon present in both releases; gtdb_downloads/ taxonomy files).
     For MAGs assigned to an R220 species, the R220 closest reference genome itself is looked up in R232.
  2. skani preview: skani_preview/skani_search_r232.tsv = `skani search -d <R232 GTDB-Tk skani db>` on the 291 local
     MAG FASTAs (the exact search GTDB-Tk 2.8.0 runs first in classify_wf), with GTDB-Tk 2.8.0's species rule:
     AF = max(AF_ref, AF_query) >= 0.5; take the closest AF-passing reference (ANI, then AF); assign its species if
     ANI >= that representative's radius (release232/radii/gtdb_radii.tsv). This is the species call GTDB-Tk will make;
     ranks above species still come from pplacer and are only approximated here (by the name mapping).

Output: tables/emilykin_r232_projection.tsv (one row per MAG) and tables/emilykin_r232_projection_summary.tsv.
Usage: python -I project_emilykin_r232.py [--summary ~/Documents/EmilyKin/gtdbtk.bac120.summary.tsv]
"""
import argparse, collections, csv, gzip, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
R232 = os.path.expanduser("~/Documents/gtdbtk_data/release232")


def key(acc):
    m = re.search(r"GC[AF]_(\d+)", acc)
    return m.group(1) if m else acc


def read_tax(path_dir, r):
    tax = {}
    for dom in ("ar53", "bac120"):
        with gzip.open(os.path.join(path_dir, f"r{r}", f"{dom}_taxonomy_r{r}.tsv.gz"), "rt") as fh:
            for line in fh:
                gid, lin = line.rstrip("\n").split("\t")
                tax[key(gid)] = (gid, lin.split(";"))
    return tax


def strip_suffix(x):
    return " ".join(re.sub(r"_[A-Z]+$", "", t) for t in x.split(" "))


def placeholder(name):
    """GTDB placeholder (alphanumeric, e.g. g__JAEXTR01, s__UBA5532 sp002426885) vs Latin name"""
    core = name[3:]
    if name.startswith("s__"):
        return bool(re.search(r" sp\d+", core)) or bool(re.search(r"\d", core.split(" ")[0]))
    return bool(re.search(r"\d", core))


def lineage_change(a, b):
    """first rank whose name differs, and how"""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.expanduser("~/Documents/EmilyKin/gtdbtk.bac120.summary.tsv"))
    ap.add_argument("--skani", default=os.path.join(HERE, "skani_preview", "skani_search_r232.tsv"))
    ap.add_argument("--downloads", default=os.path.join(HERE, "gtdb_downloads"))
    ap.add_argument("--out", default=os.path.join(HERE, "tables"))
    a = ap.parse_args()

    t220, t232 = read_tax(a.downloads, "220"), read_tax(a.downloads, "232")
    radius = {}
    with open(os.path.join(R232, "radii", "gtdb_radii.tsv")) as fh:
        for line in fh:
            sp, gid, rad = line.rstrip("\n").split("\t")[:3]
            radius[key(gid)] = float(rad)

    # R220 taxon -> dominant R232 lineage prefix (to the same rank) over the genomes in both releases
    by_taxon = collections.defaultdict(collections.Counter)
    for k, (_, lin) in t220.items():
        if k in t232:
            l2 = t232[k][1]
            for i in range(1, 6):
                by_taxon[(i, lin[i])][tuple(l2[: i + 1])] += 1

    def project(lin220):
        """deepest named R220 rank (<= genus) -> dominant R232 prefix; returns lineage, share, n, n_alternatives"""
        for i in range(5, 0, -1):
            if not lin220[i].endswith("__") and (i, lin220[i]) in by_taxon:
                c = by_taxon[(i, lin220[i])]
                pre, n = c.most_common(1)[0]
                tot = sum(c.values())
                rest = list(pre) + [RANKS[j][0] + "__" for j in range(i + 1, 7)]
                return rest, n / tot, tot, len(c)
        return None, 0, 0, 0

    # skani hits per MAG
    hits = collections.defaultdict(list)
    with open(a.skani) as fh:
        next(fh)
        for line in fh:
            t = line.rstrip("\n").split("\t")
            mag = os.path.basename(t[1]).rsplit(".fa", 1)[0]
            ref = key(os.path.basename(t[0]))
            hits[mag].append((float(t[2]), max(float(t[3]), float(t[4])) / 100, ref))

    rows, cnt = [], collections.Counter()
    with open(a.summary) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mag, lin220 = r["user_genome"], r["classification"].split(";")
            sp220 = not lin220[6].endswith("__")
            ref220 = r["closest_genome_reference"] if r["closest_genome_reference"] not in ("", "N/A") else ""
            # (1) R220 reference genome looked up in R232
            ref_in_232 = t232.get(key(ref220), ("", [""] * 7)) if ref220 else ("", [""] * 7)
            # (2) skani preview with the GTDB-Tk 2.8.0 rule
            h = sorted(hits.get(mag, []), key=lambda x: (-x[0], -x[1]))
            afp = [x for x in h if x[1] >= 0.5]
            call = None
            if afp and afp[0][2] in radius and afp[0][0] >= radius[afp[0][2]]:
                call = afp[0]
            best_any = h[0] if h else None
            best_af = afp[0] if afp else None
            proj, share, ntax, nalt = project(lin220)
            if call:
                pred = t232[call[2]][1]
                basis = "skani_species_call"
            elif proj:
                pred = proj
                basis = "name_map_from_r220_taxon"
            else:
                pred = lin220
                basis = "unmapped"
            rank, kind = lineage_change(lin220, pred)
            # genus-level change ignoring the species column
            grank, gkind = lineage_change(lin220[:6], pred[:6])
            cls = ("identical" if kind == "identical" else
                   f"species_{kind}" if rank == "species" else f"{rank}_{kind}")
            gcls = "same_to_genus" if gkind == "identical" else f"{grank}_{gkind}"
            sp_status = ("kept_same_name" if sp220 and call and pred[6] == lin220[6] else
                         "kept_renamed" if sp220 and call else
                         "lost" if sp220 else
                         "gained" if call else "none_either")
            cnt[f"class={cls}"] += 1
            cnt[f"to_genus={gcls}"] += 1
            cnt[f"species={sp_status}"] += 1
            cnt[f"basis={basis}"] += 1
            if not call and proj and share < 0.9:
                cnt["uncertain_split_taxon(share<0.9)"] += 1
            if not call and best_af and t232.get(best_af[2]) and t232[best_af[2]][1][5] != pred[5]:
                cnt["closest_AFpass_ref_in_other_genus_than_projected"] += 1
            rows.append({
                "mag": mag, "r220_method": r["classification_method"], "r220_classification": ";".join(lin220),
                "r220_species_ref": ref220, "r220_ref_lineage_in_r232": ";".join(ref_in_232[1]) if ref220 else "",
                "r232_skani_species_call": t232[call[2]][1][6] if call else "",
                "r232_call_ref": t232[call[2]][0] if call else "", "r232_call_ani": call[0] if call else "",
                "r232_call_af": round(call[1], 4) if call else "",
                "r232_closest_af_pass_ref": t232[best_af[2]][0] if best_af and best_af[2] in t232 else "",
                "r232_closest_af_pass_ani": best_af[0] if best_af else "",
                "r232_closest_af_pass_af": round(best_af[1], 4) if best_af else "",
                "r232_closest_af_pass_radius": radius.get(best_af[2], "") if best_af else "",
                "r232_closest_af_pass_lineage": ";".join(t232[best_af[2]][1]) if best_af and best_af[2] in t232 else "",
                "r232_best_any_ani": best_any[0] if best_any else "", "r232_best_any_af": round(best_any[1], 4) if best_any else "",
                "n_skani_hits": len(h),
                "projected_r232_classification": ";".join(pred), "projection_basis": basis,
                "r220_taxon_share_to_dominant": round(share, 3) if proj else "", "r220_taxon_n_genomes": ntax,
                "r220_taxon_n_r232_targets": nalt,
                "change_class": cls, "change_to_genus": gcls, "species_status": sp_status,
            })
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "emilykin_r232_projection.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(a.out, "emilykin_r232_projection_summary.tsv"), "w") as fh:
        fh.write("metric\tn_mags\n")
        for m, c in sorted(cnt.items()):
            fh.write(f"{m}\t{c}\n")
    print(open(os.path.join(a.out, "emilykin_r232_projection_summary.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
