#!/usr/bin/env python
"""Three-way comparison of the salterns MAG sets, with shared species clusters:

  original    the original study's 310 dereplicated MAGs (2020). FASTA ../MAG_contigs/.
              QC: CheckM2 1.1.0 + GUNC 1.1.1 run on poplar (vs2020/checkm2*, gunc*; copies in manuscript_numbers_tmp/),
              and KBase CheckM 1.1.2 (2026 re-QC) from comparison/mags_old_vs_new.tsv.
  first_pass  the 1,324 MAGs of the 2026 KBase pipeline (metaSPAdes --meta per sample + MEGAHIT site co-assemblies,
              KBase CheckM >= 50 / < 10). QC: KBase CheckM for all; CheckM2 + GUNC for the 776 per-sample MAGs
              (rollout_final/results/old2026_bins_checkm2_gunc.tsv; co-assembly MAGs were not re-scored).
  rollout     the 935 kept MAGs of the Arm C rollout (merge + Tadpole + default-mode SPAdes; CheckM2 >= 50 / < 10 and
              GUNC pass). rollout_final/results/mags_kept.tsv.
  coassembly  (added 2026-11-19) all 718 DAS Tool bins of the 4 first-pass MEGAHIT site co-assemblies re-scored with the
              rollout's QC (annotation_runs/rollout/scripts/coasm_qc.sh: CheckM2 + GUNC + CheckM1, rolib.py keep,
              GTDB-Tk R232); the kept ones (coasm_final/*_sediment/keep/bins.tsv, coasm_final/gtdbtk/).
The FINAL MAG set = rollout + coassembly, dereplicated jointly (final_dereplicated.tsv). The original / first pass /
rollout comparisons keep their 3-set meaning (cluster_sets3; the co-assembly genomes duplicate first-pass ones).
Taxonomy: GTDB-Tk 2.8.0 / GTDB R232 for all sets.

Species clusters: skani 0.4.0 `triangle -E --min-af 30` over all 2,569 genomes together (run on poplar with the
gtdbtk280 env; comparison_rollout/three_way_skani_edges.tsv), edges with ANI >= 95 % and max aligned fraction >= 50 %,
greedy centroid clustering, best genome first by quality = completeness - 5 x contamination (CheckM2 where the MAG has
it, else KBase CheckM). The dereplicated rollout set = the best rollout MAG of each cluster that has one.
Recovery of an original MAG (set-wide) = a direct edge (ANI >= 95, max AF >= 50) to a MAG of the other set, any
sample; manuscript_numbers.py reports the same-sample version.

Run: ~/Documents/py_venv/bin/python compare_rollout.py  -> comparison_rollout/
"""
import csv, glob, os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.dirname(HERE)
OUT = os.path.join(HERE, "comparison_rollout")
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
SETS3 = ["original", "first_pass", "rollout"]
SETS = SETS3 + ["coassembly"]
FINAL = ["rollout", "coassembly"]          # the final MAG set
LABEL = {"original": "original (310, dereplicated)", "first_pass": "2026 first pass (1,324)",
         "rollout": "Arm C rollout (935 kept)", "coassembly": "site co-assemblies, rollout QC (520 kept)"}
TAXA = [("Archaea", "domain", "Archaea"), ("Halobacteriota (haloarchaea + methanogens)", "phylum", "Halobacteriota"),
        ("Methanogens (Methanosarcinales)", "order", "Methanosarcinales"), ("Methanolobus", "genus", "Methanolobus"),
        ("Methanohalophilus", "genus", "Methanohalophilus"), ("Rhodobacteraceae", "family", "Rhodobacteraceae"),
        ("Roseovarius", "genus", "Roseovarius"), ("Marivita", "genus", "Marivita"), ("Albimonas", "genus", "Albimonas"),
        ("Halomonadaceae", "family", "Halomonadaceae"), ("Halomonas", "genus", "Halomonas"),
        ("Halanaerobiaceae", "family", "Halanaerobiaceae"), ("Halanaerobium", "genus", "Halanaerobium"),
        ("Desulfobacterota (sulfate reducers)", "phylum", "Desulfobacterota")]
SITE = lambda s: ("R2_restored" if s.startswith("R2_restored") else s.split("_")[0])


def split_tax(c):
    d = dict(x.split("__", 1) for x in str(c).split(";") if "__" in x)
    return [d.get(r[0], "") for r in RANKS]


def gunc_key(fname):  # GUNC 1.1.1 names a genome by its file name cut at the first '.fa'
    return fname.split(".fa")[0]


def load():
    m = pd.read_csv(os.path.join(HERE, "comparison/mags_old_vs_new.tsv"), sep="\t")
    idx = pd.read_csv(os.path.join(S, "MAG_contigs/index.tsv"), sep="\t").set_index("assembly_name")
    tmp = os.path.join(HERE, "manuscript_numbers_tmp")
    C2, GU = {}, {}
    for f in ("checkm2__quality_report.tsv", "checkm2_rest__quality_report.tsv"):
        C2.update({r["Name"]: r for r in csv.DictReader(open(os.path.join(tmp, f)), delimiter="\t")})
    for f in ("gunc__GUNC.progenomes_2.1.maxCSS_level.tsv", "gunc_rest__GUNC.progenomes_2.1.maxCSS_level.tsv"):
        GU.update({r["genome"]: r for r in csv.DictReader(open(os.path.join(tmp, f)), delimiter="\t")})
    rows = []
    for r in m[m.genome.str.startswith("old|")].itertuples():
        name = r.genome[4:]; f = idx.loc[name, "file"]; stem = f[:-4]          # poplar copy: <stem>.fa
        c2, gu = C2.get(stem), GU.get(gunc_key(stem + ".fa"), {})
        rows.append(dict(genome=f"original|{name}", set="original", sample=r.sample, classification=r.classification,
                         kbase_completeness=r.completeness, kbase_contamination=r.contamination,
                         checkm2_completeness=float(c2["Completeness"]) if c2 else None,
                         checkm2_contamination=float(c2["Contamination"]) if c2 else None,
                         gunc_pass=gu.get("pass.GUNC", ""), n_contigs=r.n_contigs, size_bp=r.size_bp, n50=r.n50))
    o2 = pd.read_csv(os.path.join(HERE, "rollout_final/results/old2026_bins_checkm2_gunc.tsv"), sep="\t")
    o2 = o2.set_index(o2["sample"] + "|" + o2["bin"])
    for r in m[m.genome.str.startswith("new|")].itertuples():
        name = r.genome[4:]; key = f"{r.sample}|{name.split('__')[0].removesuffix('.fasta').removesuffix('.fa')}"
        c = o2.loc[key] if key in o2.index else None
        rows.append(dict(genome=f"first_pass|{name}", set="first_pass", sample=r.sample,
                         classification=r.classification, kbase_completeness=r.completeness,
                         kbase_contamination=r.contamination,
                         checkm2_completeness=c.checkm2_completeness if c is not None else None,
                         checkm2_contamination=c.checkm2_contamination if c is not None else None,
                         gunc_pass=str(c.gunc_pass) if c is not None else "",
                         n_contigs=r.n_contigs, size_bp=r.size_bp, n50=r.n50))
    k = pd.read_csv(os.path.join(HERE, "rollout_final/results/mags_kept.tsv"), sep="\t")
    for r in k.itertuples():
        rows.append(dict(genome=f"rollout|{r.mag_id}", set="rollout", sample=r.sample,
                         classification=r.gtdb_classification, kbase_completeness=None, kbase_contamination=None,
                         checkm2_completeness=r.checkm2_completeness, checkm2_contamination=r.checkm2_contamination,
                         gunc_pass=str(r.gunc_pass), n_contigs=r.n_contigs, size_bp=r.size_bp, n50=r.n50))
    # the 4 first-pass MEGAHIT site co-assemblies, all 718 DAS Tool bins re-scored with the rollout's QC (coasm_qc.sh)
    cb = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(HERE, "coasm_final/*_sediment/keep/bins.tsv")))])
    cb = cb[cb.keep == "yes"]
    gt = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(HERE, "coasm_final/gtdbtk/gtdbtk.*.summary.tsv"))])
    gt = gt.set_index("user_genome").classification
    for r in cb.itertuples():
        rows.append(dict(genome=f"coassembly|{r.mag_id}", set="coassembly", sample=r.sample, classification=gt[r.mag_id],
                         kbase_completeness=None, kbase_contamination=None,
                         checkm2_completeness=r.checkm2_completeness, checkm2_contamination=r.checkm2_contamination,
                         gunc_pass=str(r.gunc_pass), n_contigs=r.n_contigs, size_bp=r.size_bp, n50=r.n50))
    d = pd.DataFrame(rows)
    d[RANKS] = pd.DataFrame(d.classification.map(split_tax).tolist(), index=d.index)
    d["site"] = d["sample"].map(lambda s: s.split("_sediment")[0] + " co-assembly" if str(s).endswith("_sediment")
                                else SITE(str(s)))
    has2 = d.checkm2_completeness.notna()
    d["qc_tool"] = ["CheckM2" if h else "KBase CheckM" for h in has2]
    d["completeness"] = d.checkm2_completeness.where(has2, d.kbase_completeness)
    d["contamination"] = d.checkm2_contamination.where(has2, d.kbase_contamination)
    d["quality"] = d.completeness - 5 * d.contamination
    d["mimag"] = ["high (>=90/<5)" if c >= 90 and x < 5 else "medium (>=50/<10)" if c >= 50 and x < 10 else "below MQ"
                  for c, x in zip(d.completeness, d.contamination)]
    # the rollout's KEEP rule (CheckM2 >= 50 / < 10 and GUNC pass), where CheckM2 + GUNC exist
    d["keep_rule"] = [("pass" if (c >= 50 and x < 10 and g == "True") else "fail") if q == "CheckM2" else "n/a"
                      for c, x, g, q in zip(d.checkm2_completeness, d.checkm2_contamination, d.gunc_pass, d.qc_tool)]
    assert d.genome.is_unique and d.groupby("set").size().to_dict() == {"first_pass": 1324, "original": 310,
                                                                        "rollout": 935, "coassembly": len(cb)}
    return d


def edges(d):
    inp = pd.read_csv(os.path.join(OUT, "skani_inputs.tsv"), sep="\t")
    byf = dict(zip(inp.poplar_fasta, inp.genome))
    e = pd.read_csv(os.path.join(OUT, "four_set_skani_edges.tsv"), sep="\t")
    e["a"], e["b"] = e.Ref_file.map(byf), e.Query_file.map(byf)
    assert e.a.notna().all() and e.b.notna().all()
    e["af"] = e[["Align_fraction_ref", "Align_fraction_query"]].max(axis=1)
    return e


def cluster(d, e):
    s = e[(e.ANI >= 95) & (e.af >= 50)]
    nb = {}
    for a, b in zip(s.a, s.b):
        nb.setdefault(a, set()).add(b); nb.setdefault(b, set()).add(a)
    rep = {}
    for g in d.sort_values(["quality", "completeness", "genome"], ascending=[False, False, True]).genome:
        if g in rep:
            continue
        rep[g] = g
        for h in nb.get(g, ()):
            rep.setdefault(h, g)
    d["species_cluster"] = d.genome.map(rep)
    return d, nb


def main():
    d = load()
    e = edges(d)
    d, nb = cluster(d, e)
    sets4 = d.groupby("species_cluster").set.agg(lambda s: "+".join(x for x in SETS if x in set(s)))
    d["cluster_sets"] = d.species_cluster.map(sets4)
    # the 3-set view (original / first pass / rollout) for the comparisons that predate the co-assembly re-scoring
    sets_of = d.groupby("species_cluster").set.agg(lambda s: "+".join(x for x in SETS3 if x in set(s)))
    d["cluster_sets3"] = d.species_cluster.map(sets_of)
    sets_of = sets_of[sets_of != ""]
    fin = d[d.set.isin(FINAL)]
    best_final = fin.sort_values(["quality", "completeness", "genome"], ascending=[False, False, True]) \
        .drop_duplicates("species_cluster")
    d["final_species_rep"] = d.genome.isin(set(best_final.genome))
    roll = d[d.set == "rollout"]
    best_roll = roll.sort_values(["quality", "completeness", "genome"], ascending=[False, False, True]) \
        .drop_duplicates("species_cluster")
    d["rollout_species_rep"] = d.genome.isin(set(best_roll.genome))
    # set-wide recovery of each original MAG by each other set (direct edge >= 95 / >= 50)
    for other in ("first_pass", "rollout", "coassembly"):
        d[f"recovered_by_{other}"] = [(g.startswith("original|") and any(h.startswith(other + "|") for h in nb.get(g, ())))
                                      for g in d.genome]
    d["recovered_by_final"] = d.recovered_by_rollout | d.recovered_by_coassembly
    os.makedirs(OUT, exist_ok=True)
    d.to_csv(os.path.join(OUT, "mags_three_way.tsv"), sep="\t", index=False)
    best_roll.to_csv(os.path.join(OUT, "rollout_dereplicated.tsv"), sep="\t", index=False)
    best_final.to_csv(os.path.join(OUT, "final_dereplicated.tsv"), sep="\t", index=False)

    # the two "lost" sets of originals (manuscript_numbers.py): 37 species-lost, the project's 30
    cm = pd.read_csv(os.path.join(HERE, "comparison/mags_old_vs_new.tsv"), sep="\t")
    newc = set(cm.species_cluster[cm.genome.str.startswith("new|")])
    lost37 = {g[4:] for g, c, q in zip(cm.genome, cm.species_cluster, cm.mimag)
              if g.startswith("old|") and c not in newc and q != "below MQ"}
    gl = pd.read_csv(os.path.join(HERE, "annotation_runs/tadpole/scripts/genome_check_mags.tsv"), sep="\t")
    lost30 = set(gl.mag[gl.population == "lost30"])
    assert len(lost37) == 37 and len(lost30) == 30 and lost30 <= lost37
    orig = d[d.set == "original"].assign(name=lambda x: x.genome.str[len("original|"):])

    L = []; w = L.append
    w("# Salterns MAG sets: original (2020) vs 2026 first pass vs Arm C rollout\n")
    w("Generated by `compare_rollout.py`. One species clustering over all three sets: skani 0.4.0, ANI >= 95 % and "
      "aligned fraction >= 50 %, greedy centroids by completeness - 5 x contamination. QC = CheckM2 1.1.0 where the MAG "
      "has it (all originals and rollout MAGs; the 776 per-sample first-pass MAGs), else KBase CheckM 1.1.2 (the 548 "
      "first-pass co-assembly MAGs). Taxonomy GTDB-Tk 2.8.0 / R232 everywhere.\n")
    w("Not like-for-like by design: the original set is dereplicated across samples, the first pass is KBase-CheckM-"
      "selected (>= 50 / < 10) and includes 4 site co-assemblies, and the rollout is CheckM2 + GUNC-filtered per-sample "
      "MAGs. Species clusters are the fair unit across sets.\n")

    def q(x):
        return {"MAGs": len(x), "species clusters": x.species_cluster.nunique(),
                "HQ (>=90/<5)": int(x.mimag.str.startswith("high").sum()),
                "MQ": int(x.mimag.str.startswith("medium").sum()), "below MQ": int((x.mimag == "below MQ").sum()),
                "pass the rollout KEEP rule (CheckM2 + GUNC)": f"{(x.keep_rule == 'pass').sum()} of "
                                                               f"{(x.keep_rule != 'n/a').sum()} scored",
                "median completeness %": round(x.completeness.median(), 1),
                "median contamination %": round(x.contamination.median(), 2),
                "median contigs": int(x.n_contigs.median()), "median N50 (kb)": round(x.n50.median() / 1e3, 1),
                "median size (Mb)": round(x.size_bp.median() / 1e6, 2),
                "Bacteria / Archaea": f"{(x.domain == 'Bacteria').sum()} / {(x.domain == 'Archaea').sum()}",
                "with a GTDB species": int((x.species != "").sum())}
    tq = pd.DataFrame({LABEL[s]: q(d[d.set == s]) for s in SETS})
    tq["rollout, dereplicated (best MAG per species)"] = pd.Series(q(best_roll))
    tq["FINAL: rollout + co-assembly, dereplicated"] = pd.Series(q(best_final))
    w("## 1. Numbers and quality\n"); w(tq.to_markdown()); w("")

    combo = sets_of.value_counts()
    w("## 2. Species clusters shared between sets\n")
    w(f"{len(sets_of)} species clusters contain an original, first-pass or rollout MAG (this section is that 3-set "
      f"view; {len(sets4)} clusters in all, with the co-assembly set).\n")
    w(pd.DataFrame({"species clusters": combo}).to_markdown()); w("")
    for s in SETS3:
        n = sets_of.str.contains(s).sum()
        w(f"- {LABEL[s]}: {n} species clusters.")
    w(f"- {LABEL['coassembly']}: {sets4.str.contains('coassembly').sum()} species clusters.")
    w("")
    fp = d[d.set == "first_pass"].assign(coasm=lambda x: x["sample"].str.endswith("_sediment"))
    g = fp[fp.cluster_sets3 == "first_pass"].groupby("species_cluster").agg(
        per_sample=("coasm", lambda s: (~s).any()), keep=("keep_rule", lambda s: (s == "pass").any()),
        best_c=("completeness", "max"))
    ps = fp[(~fp.coasm) & (fp.keep_rule == "pass")]
    psc = ps.groupby("species_cluster").cluster_sets3.first()
    w(f"**Why the first pass has more species ({sets_of.str.contains('first_pass').sum()}) than the rollout "
      f"({sets_of.str.contains('rollout').sum()}).** Of the {len(g)} species clusters found only by the first pass, "
      f"{(~g.per_sample).sum()} contain only MAGs of the 4 MEGAHIT site co-assemblies, which the rollout did not redo "
      f"(it replaced the 39 per-sample assemblies only); {(~g.per_sample & (g.best_c >= 90)).sum()} of those reach "
      f"KBase-CheckM completeness >= 90. Only {g.per_sample.sum()} involve per-sample MAGs, and {(g.per_sample & g.keep).sum()} "
      f"of those have a MAG passing the CheckM2 + GUNC KEEP rule. Like for like, the {len(ps)} per-sample first-pass "
      f"MAGs that pass the KEEP rule fall in {len(psc)} species clusters, {psc.str.contains('rollout').sum()} of which "
      f"have a rollout MAG ({ps.cluster_sets3.str.contains('rollout').sum()} of the {len(ps)} MAGs).\n")
    # ---- the final MAG set
    fs = best_final.assign(src=lambda x: x.species_cluster.map(fin.groupby("species_cluster").set.agg(
        lambda s: "+".join(y for y in FINAL if y in set(s)))))
    fpc = set(d.species_cluster[d.set == "first_pass"])
    w("## 2b. The final MAG set: rollout + site co-assemblies (rollout QC), dereplicated jointly\n")
    w(f"{len(fin)} MAGs ({(fin.set == 'rollout').sum()} rollout, {(fin.set == 'coassembly').sum()} co-assembly) form "
      f"{len(best_final)} species clusters; the representative is the best MAG of the cluster "
      f"(`final_dereplicated.tsv`). Species found by the rollout only: {(fs.src == 'rollout').sum()}; by the co-assemblies "
      f"only: {(fs.src == 'coassembly').sum()}; by both: {(fs.src == 'rollout+coassembly').sum()}. Representatives from "
      f"the co-assemblies: {(fs.set == 'coassembly').sum()}. HQ representatives: "
      f"{fs.mimag.str.startswith('high').sum()}; archaeal species: {(fs.domain == 'Archaea').sum()}. Of the "
      f"{len(fpc)} first-pass species, {len(fpc & set(fs.species_cluster))} are in the final set.\n")

    w("## 3. Recovery of the original MAGs (set-wide: a MAG of the other set at ANI >= 95 % and AF >= 50 %, any sample)\n")
    rec = pd.DataFrame({
        "all 310": [int(orig.recovered_by_first_pass.sum()), int(orig.recovered_by_rollout.sum())],
        "MQ+ by CheckM2 (of " + str(int((orig.mimag != 'below MQ').sum())) + ")":
            [int(orig[orig.mimag != "below MQ"].recovered_by_first_pass.sum()),
             int(orig[orig.mimag != "below MQ"].recovered_by_rollout.sum())],
        "species-lost 37": [int(orig[orig.name.isin(lost37)].recovered_by_first_pass.sum()),
                            int(orig[orig.name.isin(lost37)].recovered_by_rollout.sum())],
        "project's lost 30": [int(orig[orig.name.isin(lost30)].recovered_by_first_pass.sum()),
                              int(orig[orig.name.isin(lost30)].recovered_by_rollout.sum())]},
        index=["recovered by the 2026 first pass", "recovered by the Arm C rollout"])
    rec.loc["recovered by the FINAL set (rollout + co-assembly)"] = [
        int(orig.recovered_by_final.sum()), int(orig[orig.mimag != "below MQ"].recovered_by_final.sum()),
        int(orig[orig.name.isin(lost37)].recovered_by_final.sum()), int(orig[orig.name.isin(lost30)].recovered_by_final.sum())]
    w(rec.to_markdown()); w("")
    w("The first-pass row for the two lost sets is 0 by construction where their definition matches "
      "(species-lost: no first-pass MAG in the original's 2-way species cluster); it can be > 0 here because this "
      "clustering is joint over all sets and recovery is a direct edge. Same-sample recovery (the manuscript's numbers) is in "
      "`manuscript_numbers_<date>.json`.\n")
    l30 = orig[orig.name.isin(lost30)]
    w("Project's lost 30 recovered by the rollout: " + ", ".join(
        f"{r.genus or r.family} ({r.sample})" for r in l30[l30.recovered_by_rollout].itertuples()) + ".\n")
    w("Still not recovered: " + ", ".join(
        f"{r.genus or r.family or r.order} ({r.sample})" for r in l30[~l30.recovered_by_rollout].itertuples()) + ".\n")

    w("## 4. Best MAG per shared species cluster\n")
    def best(s):
        return d[d.set == s].sort_values(["quality", "completeness", "genome"], ascending=[False, False, True]) \
            .drop_duplicates("species_cluster").set_index("species_cluster")
    br = best("rollout")
    for s in ("original", "first_pass"):
        b = best(s); sh = b.index.intersection(br.index)
        dq = br.loc[sh, "quality"] - b.loc[sh, "quality"]
        w(f"- Rollout vs {LABEL[s]}, {len(sh)} shared clusters: the best rollout MAG has the higher quality score in "
          f"{(dq > 0).sum()}, lower in {(dq < 0).sum()}; median completeness change "
          f"{(br.loc[sh, 'completeness'] - b.loc[sh, 'completeness']).median():+.1f} points; median N50 ratio "
          f"{(br.loc[sh, 'n50'] / b.loc[sh, 'n50']).median():.2f}x.")
    w("")

    rows = []
    for lab, rank, val in TAXA:
        row = {"taxon": lab}
        for s in SETS:
            x = d[(d.set == s) & (d[rank] == val)]
            row[f"{s} MAGs (HQ)"] = f"{len(x)} ({x.mimag.str.startswith('high').sum()})"
            row[f"{s} species"] = x.species_cluster.nunique()
        x = d[d[rank] == val]
        row["species only in rollout"] = int((x.groupby("species_cluster").set.agg(set) == {"rollout"}).sum())
        row["FINAL species"] = int(best_final[best_final[rank] == val].species_cluster.nunique())
        rows.append(row)
    w("## 5. Taxa of interest (GTDB R232)\n"); w(pd.DataFrame(rows).set_index("taxon").to_markdown()); w("")

    ph = pd.DataFrame({f"{s} species": d[d.set == s].groupby("phylum").species_cluster.nunique() for s in SETS}) \
        .fillna(0).astype(int)
    w("## 6. Phyla (species clusters per set)\n")
    w(ph.sort_values("rollout species", ascending=False).head(25).to_markdown()); w("")
    si = pd.DataFrame({f"{s} MAGs": d[d.set == s].site.value_counts() for s in SETS}).fillna(0).astype(int)
    w("## 7. By site\n"); w(si.to_markdown()); w("")
    w("Files: `mags_three_way.tsv` (every MAG: set, sample, both QC values where available, KEEP-rule status, R232 "
      "lineage, contigs, N50, species cluster, which sets share the cluster, rollout species representative, set-wide "
      "recovery of originals); `rollout_dereplicated.tsv` (best rollout MAG per species cluster); "
      "`three_way_skani_edges.tsv`, `skani_inputs.tsv` (poplar paths).\n")
    open(os.path.join(OUT, "README.md"), "w").write("\n".join(L))
    # machine-readable summary for the manuscript builder (clauded_insights/make_r232_review_docx.py)
    import json
    bo, bf = best("original"), best("first_pass")
    sh_o, sh_f = bo.index.intersection(br.index), bf.index.intersection(br.index)
    arch = {s: int(d[(d.set == s) & (d.domain == "Archaea")].species_cluster.nunique()) for s in SETS}
    summ = {
        "n_mags": {s: int((d.set == s).sum()) for s in SETS},
        "n_species": {s: int(sets_of.str.contains(s).sum()) for s in SETS3},
        "n_species_total": int(len(sets_of)),
        "rollout_bacteria": int(((d.set == "rollout") & (d.domain == "Bacteria")).sum()),
        "rollout_archaea": int(((d.set == "rollout") & (d.domain == "Archaea")).sum()),
        "archaeal_species": arch,
        "shared_species_rollout_original": int(len(sh_o)), "shared_species_rollout_first_pass": int(len(sh_f)),
        "rollout_better_than_original": int(((br.loc[sh_o, "quality"] - bo.loc[sh_o, "quality"]) > 0).sum()),
        "rollout_n50_ratio_vs_original": round(float((br.loc[sh_o, "n50"] / bo.loc[sh_o, "n50"]).median()), 2),
        "rollout_better_than_first_pass": int(((br.loc[sh_f, "quality"] - bf.loc[sh_f, "quality"]) > 0).sum()),
        "rollout_n50_ratio_vs_first_pass": round(float((br.loc[sh_f, "n50"] / bf.loc[sh_f, "n50"]).median()), 2),
        "originals_recovered_setwide": {"first_pass": int(orig.recovered_by_first_pass.sum()),
                                        "rollout": int(orig.recovered_by_rollout.sum())},
        "lost30_recovered_setwide_rollout": int(l30.recovered_by_rollout.sum()),
        "lost37_recovered_setwide_rollout": int(orig[orig.name.isin(lost37)].recovered_by_rollout.sum()),
        "first_pass_only_species": int(len(g)), "first_pass_only_species_coassembly_only": int((~g.per_sample).sum()),
        "first_pass_per_sample_keep_mags": int(len(ps)),
        "first_pass_per_sample_keep_mags_with_rollout_species": int(ps.cluster_sets3.str.contains("rollout").sum()),
        "coassembly": {"bins": int(len(pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(HERE, "coasm_final/*_sediment/keep/bins.tsv"))]))),
                       "kept": int((d.set == "coassembly").sum()),
                       "kept_hq": int(((d.set == "coassembly") & d.mimag.str.startswith("high")).sum())},
        "final": {"mags": int(len(fin)), "species": int(len(best_final)),
                  "species_rollout_only": int((fs.src == "rollout").sum()),
                  "species_coassembly_only": int((fs.src == "coassembly").sum()),
                  "species_both": int((fs.src == "rollout+coassembly").sum()),
                  "hq_representatives": int(fs.mimag.str.startswith("high").sum()),
                  "archaeal_species": int((fs.domain == "Archaea").sum()),
                  "originals_recovered_setwide": int(orig.recovered_by_final.sum()),
                  "lost30_recovered_setwide": int(orig[orig.name.isin(lost30)].recovered_by_final.sum()),
                  "lost37_recovered_setwide": int(orig[orig.name.isin(lost37)].recovered_by_final.sum()),
                  "first_pass_species_in_final": int(len(fpc & set(fs.species_cluster))),
                  "first_pass_species": int(len(fpc))},
    }
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print("\n".join(L))


if __name__ == "__main__":
    main()
