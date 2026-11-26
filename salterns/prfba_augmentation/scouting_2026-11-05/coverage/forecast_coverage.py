#!/usr/bin/env python
"""Coverage forecast for the prFBA MAG arm on the 310 dereplicated saltern MAGs.

Read-only inputs (nothing outside the scratch forecast dir is written):
  * GTDB-Tk 2.8.0 / R232 summaries + ANI-screen unassigned skani hits (salterns/MAG_QC)
  * guilds_r232.csv (guild flags + CheckM completeness)
  * prFBA reference_maps/gtdb_r226_acc_to_bvbrc.tsv.gz (R226 rep -> BV-BRC genome_id)
  * GTDB R226 metadata (MD5-verified copy in an earlier prFBA session's scratch) and
    GTDB R232 metadata (prFBA region_validation copy) for R232 -> R226 cluster membership
  * pgfam_cache_index.tsv (stream-scan of bvbrc_cache/genome_gene_families.json)
  * bvbrc_cache/kbase_genome_cache/ filenames

prFBA rules mirrored (prFBA/mag_skani_hits.py): AF-aware best hit = highest ANI with
AF >= 0.15 (best_hit, L269-274; ANI_KNOBS.cov_min, L101-102); classes from classify()
(L256-266); the selector abstains below ANI_KNOBS.family_floor = 0.83 (L105-106).
AF caveat: every GTDB-Tk AF here is max(AF_ref, AF_query) (gtdbtk/external/skani.py
L431/462/487), and GTDB-Tk ran skani with its default --min-af 15, so the 0.15 gate is a
proxy that every reported hit passes by construction.
"""
import collections
import csv
import gzip
import json
import os
import re
import statistics
import sys

SAL = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns"
GT = f"{SAL}/MAG_QC/gtdbtk_local_r232/classify"
GUILDS = f"{SAL}/correlational/r232/guilds_r232.csv"
PRFBA = "/home/freiburger/Documents/prFBA"
MAP = f"{PRFBA}/reference_maps/gtdb_r226_acc_to_bvbrc.tsv.gz"
R226 = "/tmp/claude-1000/-home-freiburger-Documents-prFBA/0c6f5bc6-8224-5a34-9b20-7d5643487901/scratchpad/gtdb_r226_meta"
R232 = f"{PRFBA}/region_validation/picrust2_comparison/data/mag_arm/gtdb"
KBC = f"{PRFBA}/bvbrc_cache/kbase_genome_cache"
OUT = sys.argv[1] if len(sys.argv) > 1 else "."
PGF_INDEX = os.path.join(OUT, "pgfam_cache_index.tsv")

GUILD_COLS = ["mepn_consumer", "betaine_reducer", "tma_producer", "methanogen", "archaea"]
AF_MIN = 0.15
FLOOR = 83.0


def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def strip_pref(acc):
    return re.sub(r"^(GB_|RS_)", "", acc or "")


def core(acc):
    m = re.match(r"GC[AF]_(\d+)\.\d+$", acc or "")
    return m.group(1) if m else None


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------ GTDB-Tk R232
mags = {}
for dom, fname in (("Bacteria", "gtdbtk.bac120.summary.tsv"), ("Archaea", "gtdbtk.ar53.summary.tsv")):
    for r in csv.DictReader(open(f"{GT}/{fname}"), delimiter="\t"):
        name = r["user_genome"].removesuffix(".contigs__")
        hits = {}

        def add(acc, ani, af, src):
            acc = strip_pref(acc)
            ani, af = fnum(ani), fnum(af)
            if not acc or acc == "N/A" or ani is None:
                return
            cur = hits.get(acc)
            if cur is None or (ani, af or 0) > (cur["ani"], cur["af"] or 0):
                hits[acc] = dict(ani=ani, af=af, src=src)

        add(r["closest_genome_reference"], r["closest_genome_ani"], r["closest_genome_af"], "closest_genome")
        add(r["closest_placement_reference"], r["closest_placement_ani"], r["closest_placement_af"], "closest_placement")
        orr = r["other_related_references(genome_id,species_name,radius,ANI,AF)"]
        if orr and orr != "N/A":
            for ent in orr.split(";"):
                parts = [p.strip() for p in ent.split(",")]
                if len(parts) >= 5:
                    add(parts[0], parts[-2], parts[-1], "other_related")
        mags[name] = dict(mag=name, user_genome=r["user_genome"], gtdbtk_file_domain=dom,
                          classification=r["classification"], method=r["classification_method"],
                          closest_genome_reference=strip_pref(r["closest_genome_reference"]).replace("N/A", ""),
                          closest_genome_reference_radius=r["closest_genome_reference_radius"].replace("N/A", ""),
                          closest_genome_ani=r["closest_genome_ani"].replace("N/A", ""),
                          closest_genome_af=r["closest_genome_af"].replace("N/A", ""),
                          closest_placement_reference=strip_pref(r["closest_placement_reference"]).replace("N/A", ""),
                          closest_placement_ani=r["closest_placement_ani"].replace("N/A", ""),
                          closest_placement_af=r["closest_placement_af"].replace("N/A", ""),
                          hits=hits, n_raw_unassigned=0)

unassigned_users = set()
with gzip.open(f"{GT}/ani_screen/gtdbtk.ani_screen_unassigned_hits.tsv.gz", "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        name = r["user_genome"].removesuffix(".contigs__")
        unassigned_users.add(name)
        m = mags[name]
        m["n_raw_unassigned"] += 1
        acc = strip_pref(r["reference_genome"])
        ani, af = float(r["ani"]), float(r["af"])
        cur = m["hits"].get(acc)
        if cur is None or (ani, af) > (cur["ani"], cur["af"] or 0):
            m["hits"][acc] = dict(ani=ani, af=af, src="unassigned_raw")
        elif cur["src"] != "unassigned_raw":
            cur["src"] += "+unassigned_raw"

method_counts = collections.Counter(m["method"] for m in mags.values())
topo_ani = {k for k, m in mags.items() if m["method"] == "taxonomic classification defined by topology and ANI"}

# ------------------------------------------------------------------ guilds
guild = {r["mag"]: r for r in csv.DictReader(open(GUILDS))}
assert set(guild) == set(mags), "guild/GTDB-Tk name mismatch"

# ------------------------------------------------------------------ prFBA best-hit rule
def best(hits):
    if not hits:
        return None, "no_hit"
    good = {a: h for a, h in hits.items() if (h["af"] or 0) >= AF_MIN}
    pool = good or hits
    acc = max(pool, key=lambda a: (pool[a]["ani"], pool[a]["af"] or 0))
    if not good:
        return acc, "low_AF_unresolved"
    ani = pool[acc]["ani"]
    return acc, (">=95" if ani >= 95 else "85-95" if ani >= 85 else "<85")


def task_class(cls):
    return {"no_hit": "none", "low_AF_unresolved": "none(low_AF)"}.get(cls, cls)


def fine_class(cls, ani):
    if cls == "<85":
        return "83-85" if ani >= FLOOR else "<83"
    return cls


def tier(cls, ani):
    if cls in ("no_hit", "low_AF_unresolved"):
        return "abstain"
    if ani >= 95:
        return "species"
    if ani >= 85:
        return "genus"
    if ani >= FLOOR:
        return "family"
    return "abstain"


def af_class(af):
    if af is None:
        return "n/a"
    return ">=0.5" if af >= 0.5 else "0.15-0.5" if af >= 0.15 else "<0.15"


# ------------------------------------------------------------------ R226 map + metadata
mp, mp_core = {}, collections.defaultdict(list)
with gzip.open(MAP, "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        mp[r["gtdb_accession"]] = r
        mp_core[r["acc_core"]].append(r["gtdb_accession"])

all_accs = set()
for m in mags.values():
    all_accs |= set(m["hits"])
want_cores = {core(a) for a in all_accs if core(a)}


def scan_meta(paths, want_exact, want_cores_):
    """accession (stripped) -> dict(rep, is_rep, taxonomy) for rows matching exact or core."""
    out_exact, out_core = {}, collections.defaultdict(list)
    for p in paths:
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ia, irep, iisrep, itax = (hdr.index(c) for c in
                                      ("accession", "gtdb_genome_representative", "gtdb_representative", "gtdb_taxonomy"))
            for line in fh:
                f = line.rstrip("\n").split("\t")
                acc = strip_pref(f[ia])
                c = core(acc)
                if acc in want_exact or c in want_cores_:
                    rec = dict(acc=acc, rep=strip_pref(f[irep]), is_rep=f[iisrep] == "t", tax=f[itax])
                    if acc in want_exact:
                        out_exact[acc] = rec
                    out_core[c].append(rec)
    return out_exact, out_core


r226_exact, r226_core = scan_meta([f"{R226}/bac120_metadata_r226.tsv.gz", f"{R226}/ar53_metadata_r226.tsv.gz"],
                                  all_accs, want_cores)
# R226 reps reached via membership; need their R232 status
r226_reps_needed = {rec["rep"] for rec in r226_exact.values()} | {rec["rep"] for v in r226_core.values() for rec in v}
r232_exact, r232_core = scan_meta([f"{R232}/bac120_metadata_r232.tsv.gz", f"{R232}/ar53_metadata_r232.tsv.gz"],
                                  all_accs | r226_reps_needed, {core(a) for a in all_accs | r226_reps_needed if core(a)})

# ------------------------------------------------------------------ caches
pgf = {}
with open(PGF_INDEX) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        pgf[r["bvbrc_genome_id"]] = int(r["n_pgfams"])
kbc = {f[:-5] for f in os.listdir(KBC) if f.endswith(".json")}


def resolve_r226(acc):
    """R232 reference accession -> (status, r226_target_accession, detail)."""
    if acc in mp:
        return "A_r226_rep_exact", acc, ""
    c = core(acc)
    if c and mp_core.get(c):
        return "A_r226_rep_core", mp_core[c][0], "accession differs in prefix/version"
    rec = r226_exact.get(acc)
    if rec is None and c:
        recs = r226_core.get(c) or []
        rec = recs[0] if recs else None
    if rec is not None:
        rep = rec["rep"]
        r232rep = (r232_exact.get(rep) or {}).get("rep", "")
        same = (r232rep == acc) or (core(r232rep) == c)
        return ("B_r226_member_same_r232_cluster" if same else "B_r226_member_cluster_changed"), rep, \
               f"R226 rep {rep}; that rep's R232 rep = {r232rep or 'not in R232'}"
    return "C_absent_from_r226", "", "genome not in GTDB R226 metadata"


def bv(acc):
    r = mp.get(acc) if acc else None
    if not r:
        return dict(in_map=0, bvbrc_genome_id="", match_type="", map_cds="")
    return dict(in_map=1, bvbrc_genome_id=r["bvbrc_genome_id"], match_type=r["match_type"], map_cds=r["cds"])


def cache_cols(gid):
    if not gid:
        return dict(pgfam_cache="", n_pgfams="", kbase_cache="")
    n = pgf.get(gid)
    return dict(pgfam_cache="absent" if n is None else ("empty" if n == 0 else "nonempty"),
                n_pgfams="" if n is None else n, kbase_cache=int(gid in kbc))


# ------------------------------------------------------------------ per-MAG rows
rows = []
for name in sorted(mags):
    m, g = mags[name], guild[name]
    acc, cls = best(m["hits"])
    h = m["hits"].get(acc) if acc else None
    ani = h["ani"] if h else None
    af = h["af"] if h else None
    status, target, detail = resolve_r226(acc) if acc else ("", "", "")
    # measured R226-restricted best: best AF-pass R232 hit that is itself an R226 rep (map)
    r226hits = {a: x for a, x in m["hits"].items() if a in mp or (core(a) and mp_core.get(core(a)))}
    racc, rcls = best(r226hits)
    rh = r226hits.get(racc) if racc else None
    b1, c1 = bv(target), None
    c1 = cache_cols(b1["bvbrc_genome_id"])
    racc_map = racc if (racc in mp) else (mp_core[core(racc)][0] if racc else "")
    b2 = bv(racc_map)
    c2 = cache_cols(b2["bvbrc_genome_id"])
    # forecast tier against the R226 skani DB
    if status.startswith("A_"):
        f_tier, f_basis = tier(cls, ani), "best R232 hit is an R226 rep (measured ANI)"
    elif status == "B_r226_member_same_r232_cluster":
        f_tier = tier(cls, ani)
        f_basis = "best R232 hit's genome sits in R226 cluster whose rep shares its R232 cluster (ANI to the R226 rep not measured)"
        if rh is not None:
            f_basis += f"; measured R226-rep hit {racc} {rh['ani']}"
    elif rh is not None:
        f_tier = tier(rcls, rh["ani"])
        f_basis = f"measured ANI to an R226 rep among the R232 hits ({racc})"
        if status.startswith("B_"):
            f_basis += "; the R226 cluster of the best hit also changed in R232"
    elif status.startswith("B_"):
        f_tier, f_basis = "unknown", "best hit's R226 cluster changed in R232 and no R226 rep among measured hits"
    elif acc:
        f_tier, f_basis = "unknown", "best hit absent from R226 and no R226 rep among measured hits"
    else:
        f_tier, f_basis = "abstain", "no skani hit to any R232 rep (R226 is a subset)"
    row = dict(
        mag=name, domain=g["domain"], gtdbtk_domain_file=m["gtdbtk_file_domain"],
        **{c: g[c] for c in GUILD_COLS},
        completeness=g["completeness"], contamination=g["contamination"],
        gtdbtk_method=m["method"], gtdbtk_classification=m["classification"],
        closest_genome_reference=m["closest_genome_reference"], closest_genome_ani=m["closest_genome_ani"],
        closest_genome_af=m["closest_genome_af"], closest_genome_radius=m["closest_genome_reference_radius"],
        closest_placement_reference=m["closest_placement_reference"], closest_placement_ani=m["closest_placement_ani"],
        closest_placement_af=m["closest_placement_af"],
        n_r232_hits=len(m["hits"]), n_unassigned_raw_hits=m["n_raw_unassigned"],
        best_r232_ref=acc or "", best_r232_ani=ani if ani is not None else "", best_r232_af=af if af is not None else "",
        best_r232_src=h["src"] if h else "",
        prfba_ani_class=cls, task_ani_class=task_class(cls), fine_ani_class=fine_class(cls, ani) if ani is not None else cls,
        af_class=af_class(af), prfba_tier_r232=tier(cls, ani) if ani is not None else "abstain",
        r226_status=status, r226_target=target, r226_detail=detail,
        target_in_map=b1["in_map"], target_bvbrc_genome_id=b1["bvbrc_genome_id"], target_match_type=b1["match_type"],
        target_pgfam_cache=c1["pgfam_cache"], target_n_pgfams=c1["n_pgfams"], target_kbase_cache=c1["kbase_cache"],
        r226rep_best_measured_ref=racc_map or "", r226rep_best_measured_ani=rh["ani"] if rh else "",
        r226rep_best_measured_af=rh["af"] if rh else "",
        r226rep_bvbrc_genome_id=b2["bvbrc_genome_id"], r226rep_pgfam_cache=c2["pgfam_cache"],
        r226rep_kbase_cache=c2["kbase_cache"],
        forecast_tier_r226=f_tier, forecast_basis=f_basis,
    )
    rows.append(row)

with open(os.path.join(OUT, "forecast_per_mag.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# ------------------------------------------------------------------ class tables
CLS_ORDER = [">=95", "85-95", "<85", "none"]
FINE_ORDER = [">=95", "85-95", "83-85", "<83", "no_hit", "low_AF_unresolved"]
groups = [("all", lambda r: True), ("Bacteria", lambda r: r["domain"] == "Bacteria"),
          ("Archaea", lambda r: r["domain"] == "Archaea")] + \
         [(c, (lambda c: (lambda r: r[c] == "True"))(c)) for c in GUILD_COLS]

with open(os.path.join(OUT, "forecast_ani_classes.tsv"), "w") as fh:
    fh.write("group\tn\t" + "\t".join(CLS_ORDER) + "\t|\t" + "\t".join(f"fine:{c}" for c in FINE_ORDER) +
             "\t|\tAF>=0.5\tAF0.15-0.5\t|\ttier_species\ttier_genus\ttier_family\ttier_abstain\n")
    for gname, fn in groups:
        sub = [r for r in rows if fn(r)]
        cc = collections.Counter(r["task_ani_class"] for r in sub)
        fc = collections.Counter(r["fine_ani_class"] for r in sub)
        ac = collections.Counter(r["af_class"] for r in sub)
        tc = collections.Counter(r["prfba_tier_r232"] for r in sub)
        fh.write(f"{gname}\t{len(sub)}\t" + "\t".join(str(cc.get(c, 0)) for c in CLS_ORDER) + "\t|\t" +
                 "\t".join(str(fc.get(c, 0)) for c in FINE_ORDER) + "\t|\t" +
                 f"{ac.get('>=0.5', 0)}\t{ac.get('0.15-0.5', 0)}\t|\t" +
                 "\t".join(str(tc.get(t, 0)) for t in ("species", "genus", "family", "abstain")) + "\n")


def quart(xs):
    xs = sorted(xs)
    if not xs:
        return ["", "", "", "", ""]
    if len(xs) == 1:
        return [xs[0]] * 5
    q = statistics.quantiles(xs, n=4, method="inclusive")
    return [xs[0], round(q[0], 2), round(q[1], 2), round(q[2], 2), xs[-1]]


with open(os.path.join(OUT, "forecast_completeness_by_class.tsv"), "w") as fh:
    fh.write("class_scheme\tclass\tgroup\tn\tcompl_min\tcompl_Q1\tcompl_median\tcompl_Q3\tcompl_max\t"
             "af_min\taf_Q1\taf_median\taf_Q3\taf_max\n")
    for scheme, key, order in (("task", "task_ani_class", CLS_ORDER), ("fine", "fine_ani_class", FINE_ORDER)):
        for c in order:
            for gname, fn in groups[:3]:
                sub = [r for r in rows if r[key] == c and fn(r)]
                if not sub and gname != "all":
                    continue
                cq = quart([float(r["completeness"]) for r in sub])
                aq = quart([float(r["best_r232_af"]) for r in sub if r["best_r232_af"] != ""])
                fh.write(f"{scheme}\t{c}\t{gname}\t{len(sub)}\t" + "\t".join(map(str, cq)) + "\t" +
                         "\t".join(map(str, aq)) + "\n")

# ------------------------------------------------------------------ unique references
refrows = []
per_acc = collections.defaultdict(list)
for r in rows:
    if r["best_r232_ref"]:
        per_acc[r["best_r232_ref"]].append(r)
for acc, rs in sorted(per_acc.items()):
    status, target, detail = resolve_r226(acc)
    b = bv(target)
    c = cache_cols(b["bvbrc_genome_id"])
    refrows.append(dict(r232_ref=acc, n_mags=len(rs), max_ani=max(float(x["best_r232_ani"]) for x in rs),
                        classes=",".join(sorted({x["fine_ani_class"] for x in rs})),
                        direct_in_map=int(acc in mp), r226_status=status, r226_target=target, detail=detail,
                        bvbrc_genome_id=b["bvbrc_genome_id"], match_type=b["match_type"], map_cds=b["map_cds"],
                        pgfam_cache=c["pgfam_cache"], n_pgfams=c["n_pgfams"], kbase_cache=c["kbase_cache"]))
with open(os.path.join(OUT, "forecast_best_reference_mapping.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(refrows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(refrows)

# all R232 hit accessions (mapping burden for a multi-reference selector)
allrows = []
for acc in sorted(all_accs):
    status, target, detail = resolve_r226(acc)
    b = bv(target)
    c = cache_cols(b["bvbrc_genome_id"])
    n_m = sum(1 for m in mags.values() if acc in m["hits"])
    allrows.append(dict(r232_ref=acc, n_mags_hit=n_m, r226_status=status, r226_target=target,
                        bvbrc_genome_id=b["bvbrc_genome_id"], pgfam_cache=c["pgfam_cache"],
                        kbase_cache=c["kbase_cache"]))
with open(os.path.join(OUT, "forecast_all_hit_reference_mapping.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(allrows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(allrows)

# ------------------------------------------------------------------ summary JSON
def cnt(it):
    return dict(collections.Counter(it))


with_hit = [r for r in rows if r["best_r232_ref"]]
aug_r232 = [r for r in rows if r["prfba_tier_r232"] != "abstain"]
summ = dict(
    n_mags=len(rows),
    gtdbtk_methods=dict(method_counts),
    unassigned_hit_mags=len(unassigned_users),
    unassigned_equals_topology_and_ani=unassigned_users == topo_ani,
    ani_screen_radius_values=cnt(m["closest_genome_reference_radius"] for m in mags.values() if m["closest_genome_reference"]),
    task_classes=cnt(r["task_ani_class"] for r in rows),
    fine_classes=cnt(r["fine_ani_class"] for r in rows),
    af_classes=cnt(r["af_class"] for r in rows),
    prfba_tier_r232=cnt(r["prfba_tier_r232"] for r in rows),
    best_src=cnt(r["best_r232_src"] for r in rows),
    funnel_r232=dict(
        n=len(rows), has_r232_skani_hit=len(with_hit),
        af_pass_and_ani_ge_83=len(aug_r232),
        best_ref_direct_in_r226_map=sum(1 for r in aug_r232 if r["best_r232_ref"] in mp),
        best_ref_resolved_to_r226_target=sum(1 for r in aug_r232 if r["r226_target"]),
        target_has_bvbrc_genome_id=sum(1 for r in aug_r232 if r["target_bvbrc_genome_id"]),
        target_pgfam_nonempty=sum(1 for r in aug_r232 if r["target_pgfam_cache"] == "nonempty"),
        target_kbase_cache=sum(1 for r in aug_r232 if r["target_kbase_cache"] == 1),
    ),
    r226_status_all=cnt(r["r226_status"] for r in rows if r["best_r232_ref"]),
    r226_status_by_class={c: cnt(r["r226_status"] for r in rows if r["fine_ani_class"] == c) for c in FINE_ORDER},
    forecast_tier_r226=cnt(r["forecast_tier_r226"] for r in rows),
    forecast_tier_r226_by_domain={d: cnt(r["forecast_tier_r226"] for r in rows if r["domain"] == d) for d in ("Bacteria", "Archaea")},
    forecast_tier_r226_by_guild={c: cnt(r["forecast_tier_r226"] for r in rows if r[c] == "True") for c in GUILD_COLS},
    unique_best_refs=len(refrows),
    unique_best_refs_status=cnt(x["r226_status"] for x in refrows),
    unique_best_refs_direct_in_map=sum(x["direct_in_map"] for x in refrows),
    unique_best_refs_with_genome_id=sum(1 for x in refrows if x["bvbrc_genome_id"]),
    unique_best_refs_pgfam=cnt(x["pgfam_cache"] for x in refrows if x["bvbrc_genome_id"]),
    unique_best_refs_kbase_cache=sum(1 for x in refrows if x["kbase_cache"] == 1),
    all_hit_accessions=len(allrows),
    all_hit_status=cnt(x["r226_status"] for x in allrows),
    all_hit_unique_genome_ids=len({x["bvbrc_genome_id"] for x in allrows if x["bvbrc_genome_id"]}),
    all_hit_pgfam=cnt(x["pgfam_cache"] for x in {x["bvbrc_genome_id"]: x for x in allrows if x["bvbrc_genome_id"]}.values()),
    all_hit_kbase_cache=sum(1 for x in {x["bvbrc_genome_id"]: x for x in allrows if x["bvbrc_genome_id"]}.values() if x["kbase_cache"] == 1),
    map_rows=len(mp), map_rows_without_genome_id=sum(1 for r in mp.values() if not r["bvbrc_genome_id"]),
    pgfam_index_keys=len(pgf), kbase_cache_files=len(kbc),
    r226_meta_exact_hits=len(r226_exact), r232_meta_exact_hits=len(r232_exact),
)
json.dump(summ, open(os.path.join(OUT, "forecast_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
