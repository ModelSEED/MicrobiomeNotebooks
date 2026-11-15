#!/usr/bin/env python3
"""Weisenhorn wetland MAG census (handoff Phase 3, as revised after verification on 2026-11-13).

Population: the 12,859 MAGs of KBase narrative 272312 (DAS Tool SCG >= 50% complete / < 10% contaminated, extracted by
MetagenomeUtils; FASTAs = ~/gtdbtk_kbase/genomes, manifest all_manifest_run3_merged.tsv).

Species rule (GTDB's, as applied by GTDB-Tk): ANI >= the reference's radius (95.0 here) AND max(AF_ref, AF_query) >= 50%.
  1. GTDB side: local GTDB-Tk 2.8.0 / GTDB R232 (classify summaries). A MAG with an R232 species belongs to that cluster.
  2. MAGs with no R232 species are clustered de novo, greedy and representative-centred: MAGs in descending quality
     score join the best-matching existing centroid at ANI >= 95 / max AF >= 50, otherwise they found a new cluster.
     No MAG is moved into an R232 species through another MAG ("member join"): GTDB-Tk already tested these MAGs against
     the species representative and rejected them (review 2026-11-13: 109 of 120 such joins contradicted GTDB-Tk).
     Instead each de novo MAG is LABELLED with the R232 species it touches (skani >= 95 / >= 50 to a member) and with
     GTDB-Tk's own ANI/AF to that species' representative: below_radius, low_af (ANI >= radius, AF < 0.5: usually a
     fragment of a known species), or no_gtdbtk_value. qc_known_novel_edges.tsv lists every such boundary edge.
  3. Novel-side label (not a key): each cluster's centroid vs KBDL novel_all_representatives (SPIRE/MGnify MAGs) and
     novel_hq_v2 (isolates, the handoff's choice) with the same rule.
  4. The handoff's GTDB r220 route is reported for comparison only (skani vs gtdb_r220_consolidated, same rule).
  5. Strain clusters at 99% ANI / max AF >= 50 within each species cluster (for representative selection).
Quality score = completeness - 5 x contamination, from CheckM2 if present, else CheckM v1, else DAS Tool SCG.
skani 0.3.1 outputs (poplar run_census_compute.sh); AF in percent. Writes mag_census.tsv, clusters.tsv, headline.json.
Cluster ids for de novo clusters are a hash of the cluster's lexicographically smallest member (stable across reruns).
Sites = sampling coordinates merged by single linkage within SITE_KM (MAG-bearing runs only).
Usage: python -I census.py <compute dir with skani_*.tsv and checkm2/> <out dir>
"""
import csv, hashlib, json, math, os, re, sys
from collections import Counter, defaultdict

import pandas as pd

CD, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
H = "/home/freiburger/gtdbtk_kbase/"
SUMM = [H + "merged_run3/gtdbtk_run3/classify/gtdbtk.bac120.summary.tsv",
        H + "merged_run3/gtdbtk_run3/classify/gtdbtk.ar53.summary.tsv"]
MANIFEST = H + "all_manifest_run3_merged.tsv"
CHECKM1 = H + "narrative_run/checkm_all_MAGs.tsv"
POOL = H + "narrative_run/mag_annotation_pool.json"
META = "/home/freiburger/Documents/MicrobiomeNotebooks/wetland/Sample_metadata_wetland_only.csv"
ANI_SP, AF_MIN, ANI_STRAIN = 95.0, 50.0, 99.0
SITE_KM = 2.0
RANKS = "dpcofgs"


def stem(p):
    return re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", os.path.basename(str(p)))


def lineage_ranks(c):
    d = dict(x.split("__", 1) for x in str(c).split(";") if "__" in x)
    return {r: d.get(r, "") for r in RANKS}


# ---------- MAGs, taxonomy, quality, metadata ----------
man = pd.read_csv(MANIFEST, sep="\t", usecols=["name", "ref", "size_bp", "contigs"]).set_index("name")
mags = sorted(man.index)
summ = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in SUMM]).set_index("user_genome")
assert set(summ.index) == set(mags), "GTDB-Tk summaries do not cover exactly the manifest MAGs"
rank = {m: lineage_ranks(summ.at[m, "classification"]) for m in mags}

q = pd.DataFrame(index=mags)
ck2p = os.path.join(CD, "checkm2", "quality_report.tsv")
if os.path.exists(ck2p):
    ck2 = pd.read_csv(ck2p, sep="\t").set_index("Name")
    q = q.join(ck2[["Completeness", "Contamination", "Completeness_Model_Used", "Contig_N50", "Total_Contigs",
                    "Additional_Notes"]].add_prefix("checkm2_"))
ck1 = pd.read_csv(CHECKM1, sep="\t").set_index("Bin Id")
q = q.join(ck1[["Completeness", "Contamination", "Marker Lineage"]].add_prefix("checkm1_"))
pool = json.load(open(POOL))
das = {}
for bins in pool["sets"].values():
    for b in bins:
        das[f"{b['bid']}__{b['acc']}"] = (b["completeness"], b["contamination"], b["n50"])
q["das_completeness"] = [das.get(m, (None,) * 3)[0] for m in mags]
q["das_contamination"] = [das.get(m, (None,) * 3)[1] for m in mags]
q["das_n50"] = [das.get(m, (None,) * 3)[2] for m in mags]
for c in ("checkm2_Completeness", "checkm2_Contamination"):
    if c not in q:
        q[c] = float("nan")


def quality(m):
    r = q.loc[m]
    for src, cc, xc in (("checkm2", "checkm2_Completeness", "checkm2_Contamination"),
                        ("checkm1", "checkm1_Completeness", "checkm1_Contamination"),
                        ("das_scg", "das_completeness", "das_contamination")):
        if pd.notna(r.get(cc)) and pd.notna(r.get(xc)):
            return src, float(r[cc]), float(r[xc])
    return "none", float("nan"), float("nan")


qual = {m: quality(m) for m in mags}
score = {m: qual[m][1] - 5 * qual[m][2] for m in mags}
meta = pd.read_csv(META, dtype=str).drop_duplicates("SRA run").set_index("SRA run")
acc = {m: m.split("__")[-1] for m in mags}


def _coord(a):
    try:
        r = meta.loc[a]
        return (round(float(r["Latitude"]), 4), round(float(r["Longitude"]), 4))
    except Exception:
        return None


def _hav(p, q_):
    p1, p2, dl = math.radians(p[0]), math.radians(q_[0]), math.radians(q_[1] - p[1])
    h = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 6371 * math.asin(min(1, math.sqrt(h)))


coords = sorted({c for c in (_coord(a) for a in set(acc.values())) if c})
parent = {c: c for c in coords}
def _find(c):
    while parent[c] != c:
        parent[c] = parent[parent[c]]
        c = parent[c]
    return c
for i, a in enumerate(coords):                 # single linkage within SITE_KM
    for b in coords[i + 1:]:
        if abs(a[0] - b[0]) < 0.1 and _hav(a, b) <= SITE_KM:
            parent[_find(a)] = _find(b)
groups = defaultdict(list)
for c in coords:
    groups[_find(c)].append(c)
site_label = {c: "%.3f,%.3f" % min(g) for g in groups.values() for c in g}
site_of = {m: site_label.get(_coord(acc[m]), "") for m in mags}


def site(m):
    return site_of[m]


# ---------- skani tables ----------
def read_skani(path):
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path, sep="\t", usecols=["Ref_file", "Query_file", "ANI", "Align_fraction_ref",
                                             "Align_fraction_query"])
    d["ref"], d["qry"] = d.Ref_file.map(stem), d.Query_file.map(stem)
    d["maxaf"] = d[["Align_fraction_ref", "Align_fraction_query"]].max(axis=1)
    d["minaf"] = d[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
    return d


tri = read_skani(os.path.join(CD, "skani_triangle.tsv"))
assert tri is not None, "skani_triangle.tsv missing"
edges = tri[(tri.ANI >= ANI_SP) & (tri.maxaf >= AF_MIN)]
nbr = defaultdict(dict)                       # MAG -> {MAG: ANI} for species-level edges
for a, b, ani in zip(edges.ref, edges.qry, edges.ANI):
    if a != b:
        nbr[a][b] = max(ani, nbr[a].get(b, 0))
        nbr[b][a] = max(ani, nbr[b].get(a, 0))
s_edges = tri[(tri.ANI >= ANI_STRAIN) & (tri.maxaf >= AF_MIN)]
snbr = defaultdict(dict)
for a, b, ani in zip(s_edges.ref, s_edges.qry, s_edges.ANI):
    if a != b:
        snbr[a][b] = max(ani, snbr[a].get(b, 0))
        snbr[b][a] = max(ani, snbr[b].get(a, 0))

# ---------- species clusters ----------
cluster, source = {}, {}
for m in mags:
    if rank[m]["s"]:
        cluster[m], source[m] = "R232:" + rank[m]["s"], "R232"
# R232 species -> (representative accession, radius), from the species' own members
sp_rep = {}
for m in mags:
    if source.get(m) == "R232":
        sp_rep.setdefault(cluster[m], (str(summ.at[m, "closest_genome_reference"]),
                                       float(summ.at[m, "closest_genome_reference_radius"])))
ORR = "other_related_references(genome_id,species_name,radius,ANI,AF)"
WARN_RE = re.compile(r"Closest reference (\S+) has ANI=([\d.]+)% but alignment fraction AF=([\d.]+)")


def gtdbtk_to_ref(m, accn):
    """GTDB-Tk's own (ANI, AF) from MAG m to reference accn, or None."""
    s_ = summ.loc[m]
    for rc, ac, fc in (("closest_genome_reference", "closest_genome_ani", "closest_genome_af"),
                       ("closest_placement_reference", "closest_placement_ani", "closest_placement_af")):
        if str(s_.get(rc)) == accn and pd.notna(s_.get(ac)):
            return float(s_[ac]), float(s_[fc])
    for t in str(s_.get(ORR, "")).split(";"):
        f = [x.strip() for x in t.split(",")]
        if len(f) == 5 and f[0] == accn:
            return float(f[3]), float(f[4])
    for w in WARN_RE.finditer(str(s_.get("warnings", ""))):
        if w.group(1) == accn:
            return float(w.group(2)), float(w.group(3))
    return None


adjacent = {}                                  # de novo MAG -> label of the R232 species it touches
for m in mags:
    if m in cluster:
        continue
    best = max(((ani, o) for o, ani in nbr[m].items() if source.get(o) == "R232"), default=None)
    if not best:
        continue
    sp = cluster[best[1]]
    accn, radius = sp_rep[sp]
    g = gtdbtk_to_ref(m, accn)
    cls = ("no_gtdbtk_value" if g is None else "below_radius" if g[0] < radius else
           "low_af" if g[1] < 0.5 else "passes_rule")       # passes_rule should not occur (GTDB-Tk would assign)
    adjacent[m] = {"species": sp, "member": best[1], "edge_ani": best[0], "rep": accn,
                   "rep_ani": g[0] if g else None, "rep_af": g[1] if g else None, "class": cls}
novel = sorted((m for m in mags if m not in cluster), key=lambda m: (-score[m], m))
centroids = []                                 # greedy, representative-centred
cent_of = {}
for m in novel:
    cands = [(nbr[m][c], c) for c in centroids if c in nbr[m]]
    if cands:
        cent_of[m] = max(cands)[1]
    else:
        centroids.append(m)
        cent_of[m] = m
cl_members = defaultdict(list)
for m in novel:
    cl_members[cent_of[m]].append(m)
cid = {c: "WN_" + hashlib.sha1(min(cl_members[c]).encode()).hexdigest()[:8] for c in centroids}
assert len(set(cid.values())) == len(cid), "cluster id hash collision"
centroid_of_cluster = {}
for m in novel:
    cluster[m], source[m] = "novel:" + cid[cent_of[m]], "de_novo"
for c in centroids:
    centroid_of_cluster[cluster[c]] = c

# ---------- strain clusters (99%) within each species cluster ----------
members = defaultdict(list)
for m in mags:
    members[cluster[m]].append(m)
strain = {}
for cl, ms in members.items():
    cents = []
    for m in sorted(ms, key=lambda m: (-score[m], m)):
        c = max(((snbr[m][x], x) for x in cents if x in snbr[m]), default=None)
        if c:
            strain[m] = strain[c[1]]
        else:
            cents.append(m)
            strain[m] = f"{cl}|S{len(cents)}"


# ---------- reference-DB labels (novel side; r220 comparison) ----------
def best_hits(name):
    d = read_skani(os.path.join(CD, f"skani_vs_{name}.tsv"))
    if d is None:
        return None
    ok = d[(d.ANI >= ANI_SP) & (d.maxaf >= AF_MIN)].sort_values("ANI", ascending=False).drop_duplicates("qry")
    anyhit = d.sort_values("ANI", ascending=False).drop_duplicates("qry")
    return ({r.qry: (r.ref, r.ANI, r.maxaf) for r in ok.itertuples()},
            {r.qry: (r.ref, r.ANI, r.maxaf) for r in anyhit.itertuples()})


dbs = {n: best_hits(n) for n in ("novel_all_representatives", "novel_hq_v2", "gtdb_r220_consolidated")}

# ---------- per-MAG table ----------
first_unnamed = {}
for m in mags:
    r = rank[m]
    first_unnamed[m] = next((k for k in RANKS if not r[k]), "none")   # 'none' = has species
rows = []
for m in mags:
    s = summ.loc[m]
    src, comp, cont = qual[m]
    row = {"mag": m, "kbase_ref": man.at[m, "ref"], "accession": acc[m],
           "bioproject": meta.at[acc[m], "Bioproject"] if acc[m] in meta.index else "",
           "habitat": meta.at[acc[m], "Habitat"] if acc[m] in meta.index else "", "site": site(m),
           "size_bp": man.at[m, "size_bp"], "contigs": man.at[m, "contigs"],
           "r232_classification": s["classification"], "r232_method": s["classification_method"],
           "r232_closest_reference": s["closest_genome_reference"], "r232_closest_ani": s["closest_genome_ani"],
           "r232_closest_af": s["closest_genome_af"], "r232_red": s["red_value"], "r232_warnings": s["warnings"],
           "first_unnamed_rank": first_unnamed[m], "domain": rank[m]["d"], "phylum": rank[m]["p"],
           "species_cluster": cluster[m], "cluster_source": source[m], "strain_cluster": strain[m],
           "quality_source": src, "completeness": comp, "contamination": cont, "score": round(score[m], 2),
           "checkm2_completeness": q.at[m, "checkm2_Completeness"], "checkm2_contamination": q.at[m, "checkm2_Contamination"],
           "checkm2_model": q.at[m, "checkm2_Completeness_Model_Used"] if "checkm2_Completeness_Model_Used" in q else "",
           "checkm2_n50": q.at[m, "checkm2_Contig_N50"] if "checkm2_Contig_N50" in q else "",
           "checkm1_completeness": q.at[m, "checkm1_Completeness"], "checkm1_contamination": q.at[m, "checkm1_Contamination"],
           "das_completeness": q.at[m, "das_completeness"], "das_contamination": q.at[m, "das_contamination"],
           "das_n50": q.at[m, "das_n50"]}
    a = adjacent.get(m)
    row.update({"r232_adjacent_species": a["species"] if a else "", "r232_adjacent_member": a["member"] if a else "",
                "r232_adjacent_edge_ani": a["edge_ani"] if a else "", "r232_adjacent_rep_ani": a["rep_ani"] if a else "",
                "r232_adjacent_rep_af": a["rep_af"] if a else "", "r232_adjacent_class": a["class"] if a else ""})
    for n, v in dbs.items():
        if v is None:
            continue
        hit = v[0].get(m)
        anyh = v[1].get(m)
        row[f"{n}_species_hit"] = hit[0] if hit else ""
        row[f"{n}_species_ani"] = hit[1] if hit else ""
        row[f"{n}_best_ani"] = anyh[1] if anyh else ""
        row[f"{n}_best_af"] = anyh[2] if anyh else ""
    rows.append(row)
df = pd.DataFrame(rows)
hab = dict(zip(df.mag, df.habitat))
df.to_csv(os.path.join(OUT, "mag_census.tsv"), sep="\t", index=False)
# boundary QC: species-level (>=95/50) and strain-level (>=99/50) edges between an R232 cluster and a de novo cluster
qc = edges[(edges.ref != edges.qry)].copy()
qc["cl_a"], qc["cl_b"] = qc.ref.map(cluster), qc.qry.map(cluster)
qc = qc[(qc.cl_a.str.startswith("R232:")) != (qc.cl_b.str.startswith("R232:"))]
qc["strain_level"] = qc.ANI >= ANI_STRAIN
qc[["ref", "qry", "ANI", "Align_fraction_ref", "Align_fraction_query", "cl_a", "cl_b", "strain_level"]].to_csv(
    os.path.join(OUT, "qc_known_novel_edges.tsv"), sep="\t", index=False)
# within-cluster and cross-cluster strain-level edges (select_reps.py uses them to avoid near-identical representatives)
se = s_edges[s_edges.ref != s_edges.qry]
se[["ref", "qry", "ANI", "maxaf"]].rename(columns={"ref": "mag_a", "qry": "mag_b"}).to_csv(
    os.path.join(OUT, "strain_edges.tsv"), sep="\t", index=False)

# ---------- clusters ----------
crow = []
for cl, ms in members.items():
    rep = min(ms, key=lambda m: (-score[m], m))            # same order as the greedy pass: = centroid for de novo
    lab = centroid_of_cluster.get(cl, rep)
    adj = Counter(adjacent[m]["class"] for m in ms if m in adjacent)
    r = {"species_cluster": cl, "source": "R232" if cl.startswith("R232:") else "de_novo",
         "n_mags": len(ms), "n_strains": len({strain[m] for m in ms}), "n_samples": len({acc[m] for m in ms}),
         "n_sites": len({site(m) for m in ms} - {""}), "habitats": ";".join(sorted({hab[m] for m in ms} - {""})),
         "best_mag": rep, "best_score": round(score[rep], 2), "lineage": summ.at[rep, "classification"],
         "first_unnamed_rank": first_unnamed[rep], "n_mags_gt70_lt10": sum(1 for m in ms if qual[m][1] > 70 and qual[m][2] < 10),
         "centroid": lab if cl.startswith("novel:") else "",
         "r232_adjacent_species": ";".join(sorted({adjacent[m]["species"] for m in ms if m in adjacent})),
         "n_adjacent_below_radius": adj.get("below_radius", 0), "n_adjacent_low_af": adj.get("low_af", 0),
         "n_adjacent_no_value": adj.get("no_gtdbtk_value", 0),
         "all_members_low_af_adjacent": bool(ms) and adj.get("low_af", 0) == len(ms)}
    for n, v in dbs.items():
        if v is not None and n.startswith("novel"):
            r[f"centroid_in_{n}"] = bool(v[0].get(lab))
            r[f"any_member_in_{n}"] = any(m in v[0] for m in ms)
    crow.append(r)
cdf = pd.DataFrame(crow).sort_values(["source", "n_mags"], ascending=[True, False])
cdf.to_csv(os.path.join(OUT, "clusters.tsv"), sep="\t", index=False)

# ---------- headline ----------
known = cdf.source.eq("R232")
qclu = cdf.n_mags_gt70_lt10 > 0
head = {
    "population_mags": len(mags), "assemblies_with_mags": df.accession.nunique(),
    "quality_source_counts": dict(Counter(qual[m][0] for m in mags)),
    "mags_with_r232_species": int((df.cluster_source == "R232").sum()),
    "mags_de_novo_r232_adjacent": dict(Counter(df.loc[df.r232_adjacent_class != "", "r232_adjacent_class"])),
    "clusters_de_novo_all_members_low_af_adjacent": int(cdf.all_members_low_af_adjacent.sum()),
    "clusters_de_novo_excluding_low_af_fragment_clusters": int((~known).sum() - cdf.all_members_low_af_adjacent.sum()),
    "qc_known_novel_edges": {"species_level": int(len(qc)), "strain_level": int(qc.strain_level.sum())},
    "mags_de_novo": int((df.cluster_source == "de_novo").sum()),
    "species_clusters_total": len(cdf), "clusters_r232": int(known.sum()), "clusters_de_novo": int((~known).sum()),
    "clusters_with_a_mag_gt70_lt10": int(qclu.sum()),
    "clusters_r232_with_a_mag_gt70_lt10": int((known & qclu).sum()),
    "strain_clusters_total": df.strain_cluster.nunique(),
    "mag_share_known_species_pct": round(100 * (df.cluster_source != "de_novo").mean(), 1),
    "cluster_share_known_pct": round(100 * known.mean(), 1),
    "cluster_share_known_quality_clusters_pct": round(100 * (known & qclu).sum() / max(qclu.sum(), 1), 1),
    "novelty_depth_of_de_novo_mags": dict(Counter(df.loc[df.cluster_source == "de_novo", "first_unnamed_rank"])),
    "novelty_depth_of_de_novo_clusters": dict(Counter(cdf.loc[~known, "first_unnamed_rank"])),
    "domains": dict(Counter(df.domain)),
    "by_habitat_pct_known_species": {h: round(100 * (g.cluster_source != "de_novo").mean(), 1)
                                     for h, g in df.groupby("habitat")},
    "sites": len(set(site_of.values()) - {""}),
    "thresholds": {"ani_species": ANI_SP, "max_af_min_pct": AF_MIN, "ani_strain": ANI_STRAIN, "site_km": SITE_KM, "skani": open(
        os.path.join(CD, "skani_version.txt")).read().strip() if os.path.exists(os.path.join(CD, "skani_version.txt")) else "?"},
}
for n, v in dbs.items():
    if v is None:
        continue
    hits = set(v[0])
    head[f"{n}_mags_with_species_level_hit"] = len(hits)
    if n.startswith("novel"):
        dn = cdf[~known]
        head[f"{n}_de_novo_clusters_seen"] = int(dn[f"centroid_in_{n}"].sum())
    else:
        r232sp = set(df.loc[df.cluster_source == "R232", "mag"])
        head["r220_vs_r232"] = {"r220_species_hits": len(hits), "r232_species": len(r232sp),
                                "both": len(hits & r232sp), "r232_only": len(r232sp - hits), "r220_only": len(hits - r232sp)}
json.dump(head, open(os.path.join(OUT, "headline.json"), "w"), indent=1, default=str)
print(json.dumps(head, indent=1, default=str))
