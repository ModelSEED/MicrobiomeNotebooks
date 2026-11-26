#!/usr/bin/env python
"""Join the measured prFBA skani search (310 MAGs vs GTDB R226, mag_skani_hits.py, -n 20
--min-af 0) to the R232-based forecast; map the measured best hit and the selector-relevant
hits to BV-BRC genome_ids and the prFBA caches. Writes measured_r226_per_mag.tsv,
measured_r226_classes.tsv, forecast_vs_measured.tsv and measured_r226_summary.json.
Usage: compare_r226_measured.py FORECAST_DIR
"""
import collections
import csv
import gzip
import json
import os
import statistics
import sys

OUT = sys.argv[1]
SK = os.path.join(OUT, "skani_r226", "out")
PRFBA = "/home/freiburger/Documents/prFBA"
MAP = f"{PRFBA}/reference_maps/gtdb_r226_acc_to_bvbrc.tsv.gz"
KBC = f"{PRFBA}/bvbrc_cache/kbase_genome_cache"
GUILD_COLS = ["mepn_consumer", "betaine_reducer", "tma_producer", "methanogen", "archaea"]
FLOOR = 0.83
BAND = 0.005   # ANI_KNOBS.band (mag_skani_hits.py L107): 0.5 ANI points = a tie

fc = {r["mag"]: r for r in csv.DictReader(open(os.path.join(OUT, "forecast_per_mag.tsv")), delimiter="\t")}
hits = json.load(open(os.path.join(SK, "mag_top20_skani_hits.gtdb_r226.json")))
mp = {}
with gzip.open(MAP, "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        mp[r["gtdb_accession"]] = r
pgf = {}
with open(os.path.join(OUT, "pgfam_cache_index.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        pgf[r["bvbrc_genome_id"]] = int(r["n_pgfams"])
kbc = {f[:-5] for f in os.listdir(KBC) if f.endswith(".json")}


def tier(bi):
    if bi is None:
        return "abstain"
    return "species" if bi >= 0.95 else "genus" if bi >= 0.85 else "family" if bi >= FLOOR else "abstain"


def gid_info(acc):
    r = mp.get(acc)
    gid = r["bvbrc_genome_id"] if r else ""
    n = pgf.get(gid) if gid else None
    return dict(in_map=int(r is not None), gid=gid, match_type=r["match_type"] if r else "",
                is_mag_ref=(r["ncbi_genome_category"] != "none") if r else "",
                pgfam="" if not gid else ("absent" if n is None else "empty" if n == 0 else "nonempty"),
                n_pgfams="" if n is None else n, kbase=int(gid in kbc) if gid else "")


rows = []
band_accs = set()
top_accs = set()
for key, rec in sorted(hits.items()):
    mag = key.removesuffix(".contigs__")
    f = fc[mag]
    bi = rec["best_identity"]
    t = tier(bi)
    gi = gid_info(rec["best_ref"]) if rec["best_ref"] else gid_info("")
    # AF-passing hits within the tie band of the best (the anchor handful the selector starts from)
    band = [h for h in rec["top20"] if h["af_pass"] and bi is not None and h["identity"] >= bi - BAND]
    for h in band:
        band_accs.add(h["genome_id"])
    for h in rec["top20"]:
        top_accs.add(h["genome_id"])
    best_any = max(rec["top20"], key=lambda h: h["ani"]) if rec["top20"] else None
    rows.append(dict(
        mag=mag, domain=f["domain"], **{c: f[c] for c in GUILD_COLS}, completeness=f["completeness"],
        mag_len=rec["mag_len"], n_contigs=rec["n_contigs"], n_hits=rec["n_hits"],
        r226_ani_class=rec["ani_class"], r226_best_ref=rec["best_ref"] or "",
        r226_best_ani=rec["best_ani"] if rec["best_ani"] is not None else "",
        r226_best_af_query=rec["best_af_query"] if rec["best_af_query"] is not None else "",
        r226_best_af_ref=rec["best_af_ref"] if rec["best_af_ref"] is not None else "",
        r226_af_pass=int(bool(rec["af_pass"])), r226_selector_tier=t,
        r226_max_af_query_any_hit=max((h["af_query"] for h in rec["top20"]), default=""),
        r226_max_ani_any_hit=best_any["ani"] if best_any else "",
        n_band_hits=len(band),
        best_in_map=gi["in_map"], best_bvbrc_genome_id=gi["gid"], best_match_type=gi["match_type"],
        best_ref_is_mag=gi["is_mag_ref"], best_pgfam_cache=gi["pgfam"], best_n_pgfams=gi["n_pgfams"],
        best_kbase_cache=gi["kbase"],
        r232_task_class=f["task_ani_class"], r232_fine_class=f["fine_ani_class"], r232_best_ref=f["best_r232_ref"],
        r232_best_ani=f["best_r232_ani"], r232_best_af=f["best_r232_af"], r232_r226_status=f["r226_status"],
        r232_prfba_tier=f["prfba_tier_r232"], forecast_tier_r226=f["forecast_tier_r226"],
        same_best_ref_as_r232=int(bool(rec["best_ref"]) and rec["best_ref"] == f["best_r232_ref"]),
    ))

with open(os.path.join(OUT, "measured_r226_per_mag.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

groups = [("all", lambda r: True), ("Bacteria", lambda r: r["domain"] == "Bacteria"),
          ("Archaea", lambda r: r["domain"] == "Archaea")] + \
         [(c, (lambda c: (lambda r: r[c] == "True"))(c)) for c in GUILD_COLS]
CL = [">=95", "85-95", "<85", "low_AF_unresolved", "no_hit"]
TI = ["species", "genus", "family", "abstain"]


def quart(xs):
    xs = sorted(xs)
    if not xs:
        return ["", "", "", "", ""]
    if len(xs) == 1:
        return [xs[0]] * 5
    q = statistics.quantiles(xs, n=4, method="inclusive")
    return [xs[0], round(q[0], 2), round(q[1], 2), round(q[2], 2), xs[-1]]


with open(os.path.join(OUT, "measured_r226_classes.tsv"), "w") as fh:
    fh.write("group\tn\t" + "\t".join(CL) + "\t|\t" + "\t".join(f"tier_{t}" for t in TI) +
             "\t|\taugmentable(tier!=abstain)\tbest_has_genome_id\tbest_pgfam_cached\tbest_kbase_cached\n")
    for g, fn in groups:
        sub = [r for r in rows if fn(r)]
        cc = collections.Counter(r["r226_ani_class"] for r in sub)
        tc = collections.Counter(r["r226_selector_tier"] for r in sub)
        aug = [r for r in sub if r["r226_selector_tier"] != "abstain"]
        fh.write(f"{g}\t{len(sub)}\t" + "\t".join(str(cc.get(c, 0)) for c in CL) + "\t|\t" +
                 "\t".join(str(tc.get(t, 0)) for t in TI) + "\t|\t" +
                 f"{len(aug)}\t{sum(1 for r in aug if r['best_bvbrc_genome_id'])}\t"
                 f"{sum(1 for r in aug if r['best_pgfam_cache'] == 'nonempty')}\t"
                 f"{sum(1 for r in aug if r['best_kbase_cache'] == 1)}\n")
    fh.write("\nclass\tgroup\tn\tcompl_min\tcompl_Q1\tcompl_median\tcompl_Q3\tcompl_max\t"
             "af_query_min\taf_query_Q1\taf_query_median\taf_query_Q3\taf_query_max\n")
    for c in CL:
        for g, fn in groups[:3]:
            sub = [r for r in rows if r["r226_ani_class"] == c and fn(r)]
            if not sub and g != "all":
                continue
            cq = quart([float(r["completeness"]) for r in sub])
            aq = quart([float(r["r226_best_af_query"]) for r in sub if r["r226_best_af_query"] != ""])
            fh.write(f"{c}\t{g}\t{len(sub)}\t" + "\t".join(map(str, cq)) + "\t" + "\t".join(map(str, aq)) + "\n")

xt_class = collections.Counter((r["r232_fine_class"], r["r226_ani_class"]) for r in rows)
xt_tier = collections.Counter((r["forecast_tier_r226"], r["r226_selector_tier"]) for r in rows)
xt_status = collections.Counter((r["r232_r226_status"] or "no_r232_hit", r["r226_ani_class"]) for r in rows)
with open(os.path.join(OUT, "forecast_vs_measured.tsv"), "w") as fh:
    fh.write("table\tforecast_or_r232\tmeasured_r226\tn\n")
    for (a, b), n in sorted(xt_class.items()):
        fh.write(f"r232_fine_class_x_r226_class\t{a}\t{b}\t{n}\n")
    for (a, b), n in sorted(xt_tier.items()):
        fh.write(f"forecast_tier_x_measured_tier\t{a}\t{b}\t{n}\n")
    for (a, b), n in sorted(xt_status.items()):
        fh.write(f"r232_best_ref_r226_status_x_r226_class\t{a}\t{b}\t{n}\n")

aug = [r for r in rows if r["r226_selector_tier"] != "abstain"]
band_gids = {mp[a]["bvbrc_genome_id"] for a in band_accs if a in mp and mp[a]["bvbrc_genome_id"]}
top_gids = {mp[a]["bvbrc_genome_id"] for a in top_accs if a in mp and mp[a]["bvbrc_genome_id"]}


def cache_split(gids):
    c = collections.Counter("absent" if g not in pgf else "empty" if pgf[g] == 0 else "nonempty" for g in gids)
    return dict(c, kbase_cache=sum(1 for g in gids if g in kbc))


lowaf = [r for r in rows if r["r226_ani_class"] == "low_AF_unresolved"]
summ = dict(
    n=len(rows),
    r226_classes=dict(collections.Counter(r["r226_ani_class"] for r in rows)),
    r226_selector_tiers=dict(collections.Counter(r["r226_selector_tier"] for r in rows)),
    augmentable=len(aug),
    augmentable_by_domain={d: sum(1 for r in aug if r["domain"] == d) for d in ("Bacteria", "Archaea")},
    augmentable_by_guild={c: f"{sum(1 for r in aug if r[c] == 'True')}/{sum(1 for r in rows if r[c] == 'True')}"
                          for c in GUILD_COLS},
    aug_best_in_map=sum(r["best_in_map"] for r in aug),
    aug_best_has_genome_id=sum(1 for r in aug if r["best_bvbrc_genome_id"]),
    aug_best_pgfam=dict(collections.Counter(r["best_pgfam_cache"] for r in aug)),
    aug_best_kbase_cache=sum(1 for r in aug if r["best_kbase_cache"] == 1),
    aug_best_ref_is_mag=sum(1 for r in aug if r["best_ref_is_mag"] is True),
    unique_best_refs_aug=len({r["r226_best_ref"] for r in aug}),
    unique_best_gids_aug=len({r["best_bvbrc_genome_id"] for r in aug if r["best_bvbrc_genome_id"]}),
    unique_best_gids_aug_cache=cache_split({r["best_bvbrc_genome_id"] for r in aug if r["best_bvbrc_genome_id"]}),
    band_hit_accessions=len(band_accs),
    band_unmapped=sum(1 for a in band_accs if not (a in mp and mp[a]["bvbrc_genome_id"])),
    band_gids=len(band_gids), band_gids_cache=cache_split(band_gids),
    top20_accessions=len(top_accs),
    top20_unmapped=sum(1 for a in top_accs if not (a in mp and mp[a]["bvbrc_genome_id"])),
    top20_gids=len(top_gids), top20_gids_cache=cache_split(top_gids),
    low_af_n=len(lowaf),
    low_af_r232_class=dict(collections.Counter(r["r232_fine_class"] for r in lowaf)),
    low_af_best_af_query_quartiles=quart([float(r["r226_best_af_query"]) for r in lowaf]),
    low_af_best_ani_quartiles=quart([float(r["r226_best_ani"]) for r in lowaf]),
    low_af_max_af_query_any_hit_quartiles=quart([float(r["r226_max_af_query_any_hit"]) for r in lowaf]),
    same_best_ref_as_r232=sum(r["same_best_ref_as_r232"] for r in rows),
    forecast_tier_vs_measured={f"{a}->{b}": n for (a, b), n in sorted(xt_tier.items())},
)
json.dump(summ, open(os.path.join(OUT, "measured_r226_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
