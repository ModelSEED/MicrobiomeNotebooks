#!/usr/bin/env python3
"""Weisenhorn wetland MAGs: representative worklist for the functional census (handoff Phase 5, revised 2026-11-13).

Input: mag_census.tsv from census.py. Rules (fixed here so the worklist is reproducible):
  eligible    completeness > 70 and contamination < 10 (quality_source: CheckM2 > CheckM v1 > DAS Tool SCG), and no
              GTDB-Tk "markers with multiple hits" warning (a contamination/chimerism flag)
  score       completeness - 5 x contamination + 0.5 x log10(N50); ties: fewer contigs, then name
  per species cluster: walk the eligible MAGs by score; take the best MAG of each 99%-ANI strain cluster not yet
              represented, preferring (as a tie-break, not a quota) a sampling site not yet represented; cap MAX_REPS.
              A candidate with a >= 99% ANI / max AF >= 50 edge to ANY representative already chosen (in any cluster;
              strain_edges.tsv) is skipped, so near-identical genomes are never annotated twice. R232 clusters are
              processed first, so a de novo MAG next to a known species yields to that species' representative.
  relaxation  only for a cluster with NO eligible MAG: one rep at >= 50 / < 10 without the multiple-hit flag,
              relaxed=True. Clusters with nothing above that floor are listed as census-only (no annotation).
Outputs: rep_worklist.tsv (one row per chosen MAG), cluster_reps.tsv (per cluster: n chosen, relaxed, reason),
rep_summary.json. Usage: python -I select_reps.py <census out dir> [MAX_REPS=10]
"""
import json, math, os, sys
from collections import Counter, defaultdict

import pandas as pd

D = sys.argv[1]
MAX_REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 10
df = pd.read_csv(os.path.join(D, "mag_census.tsv"), sep="\t", dtype={"site": str})
df["multihit"] = df.r232_warnings.fillna("").str.contains("multiple hits", case=False)
n50 = pd.to_numeric(df.checkm2_n50, errors="coerce").fillna(pd.to_numeric(df.das_n50, errors="coerce")).fillna(1)
df["rep_score"] = df.completeness - 5 * df.contamination + 0.5 * n50.clip(lower=1).map(math.log10)
df["eligible"] = (df.completeness > 70) & (df.contamination < 10) & ~df.multihit
df["floor_ok"] = (df.completeness >= 50) & (df.contamination < 10) & ~df.multihit
df = df.sort_values(["rep_score", "contigs", "mag"], ascending=[False, True, True])
adj = defaultdict(set)
se = pd.read_csv(os.path.join(D, "strain_edges.tsv"), sep="\t")
for a, b in zip(se.mag_a, se.mag_b):
    adj[a].add(b)
    adj[b].add(a)
chosen_all = set()
order = sorted(df.species_cluster.unique(), key=lambda c: (not c.startswith("R232:"), c))
groups = dict(tuple(df.groupby("species_cluster", sort=False)))

picks, crow = [], []
for cl in order:
    g = groups[cl]
    free = ~g.mag.map(lambda m: bool(adj[m] & chosen_all))
    el = g[g.eligible & free]
    relaxed = False
    n_str, n_el_str = int(g.strain_cluster.nunique()), int(g.loc[g.eligible, "strain_cluster"].nunique())
    base = {"species_cluster": cl, "n_mags": len(g), "n_eligible": int(g.eligible.sum()), "n_strains": n_str,
            "n_eligible_strains": n_el_str}
    if el.empty:
        fl = g[g.floor_ok & free]
        if fl.empty:
            why = ("census-only: every eligible MAG is >= 99% ANI to a representative already chosen"
                   if g.eligible.any() or g.floor_ok.any() else "census-only: no MAG >= 50/<10 without the multiple-hit flag")
            crow.append({**base, "n_reps": 0, "relaxed": False, "reason": why})
            continue
        chosen, relaxed = [fl.iloc[0]], True
        chosen_all.add(fl.iloc[0].mag)
    else:
        chosen, strains, sites = [], set(), set()
        pool = list(el.itertuples(index=False))
        while pool and len(chosen) < MAX_REPS:
            fresh = [r for r in pool if r.strain_cluster not in strains and not (adj[r.mag] & chosen_all)]
            if not fresh:
                break
            top = fresh[0].rep_score
            # among strain-new candidates within 1 score point of the best, prefer a new site
            near = [r for r in fresh if top - r.rep_score <= 1.0]
            pick = next((r for r in near if r.site and r.site not in sites), fresh[0])
            chosen.append(pick)
            chosen_all.add(pick.mag)
            strains.add(pick.strain_cluster)
            if pick.site:
                sites.add(pick.site)
            pool = [r for r in pool if r is not pick]
        chosen = [pd.Series(r._asdict()) for r in chosen]
    for i, r in enumerate(chosen, 1):
        d = r.to_dict()
        d.update(rep_rank=i, relaxed=relaxed)
        picks.append(d)
    crow.append({**base, "n_reps": len(chosen), "relaxed": relaxed,
                 "reason": "relaxed to the >= 50/<10 floor: no eligible MAG" if relaxed else
                 (f"capped at {MAX_REPS} of {n_el_str} eligible strains" if n_el_str > MAX_REPS else
                  "all eligible strains represented" + (f" ({n_str - n_el_str} strains have no eligible MAG)"
                                                        if n_str > n_el_str else ""))})

wl = pd.DataFrame(picks)
cols = ["mag", "kbase_ref", "species_cluster", "cluster_source", "strain_cluster", "rep_rank", "relaxed",
        "quality_source", "completeness", "contamination", "rep_score", "checkm2_completeness", "checkm2_contamination",
        "checkm1_completeness", "checkm1_contamination", "das_completeness", "das_contamination", "size_bp", "contigs",
        "r232_classification", "first_unnamed_rank", "r232_adjacent_species", "r232_adjacent_class", "accession",
        "bioproject", "habitat", "site"]
wl[[c for c in cols if c in wl]].to_csv(os.path.join(D, "rep_worklist.tsv"), sep="\t", index=False)
cr = pd.DataFrame(crow)
cr.to_csv(os.path.join(D, "cluster_reps.tsv"), sep="\t", index=False)
s = {"max_reps_per_cluster": MAX_REPS, "clusters": len(cr), "clusters_with_reps": int((cr.n_reps > 0).sum()),
     "clusters_relaxed": int(cr.relaxed.sum()), "clusters_census_only": int((cr.n_reps == 0).sum()),
     "worklist_mags": len(wl), "worklist_relaxed_mags": int(wl.relaxed.sum()),
     "worklist_by_source": dict(Counter(wl.cluster_source)), "worklist_by_domain": dict(Counter(wl.r232_classification.str.split(";").str[0])),
     "worklist_quality_source": dict(Counter(wl.quality_source)),
     "worklist_mimag_hq_by_completeness_contamination": int(((wl.completeness > 90) & (wl.contamination < 5)).sum()),
     "excluded_multihit_mags": int(df.multihit.sum()),
     "reps_per_cluster_distribution": {int(k): v for k, v in sorted(Counter(cr.n_reps).items())},
     "worklist_r232_adjacent": dict(Counter(wl.r232_adjacent_class.fillna("").replace("", "none")))}
json.dump(s, open(os.path.join(D, "rep_summary.json"), "w"), indent=1, default=str)
print(json.dumps(s, indent=1, default=str))
