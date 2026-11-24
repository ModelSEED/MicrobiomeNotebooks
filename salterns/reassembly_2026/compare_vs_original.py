#!/usr/bin/env python
"""Diff report: the original study's 310 dereplicated MAGs vs the 2026 re-assembly's 1,324 MAGs.

Both sides use GTDB-Tk 2.8.0 / GTDB R232 and KBase CheckM (lineage_wf):
  original : ../MAG_QC/comparison_vs_2020.csv (2026 re-QC of the 310: Completeness_2026, Contamination_2026, gtdb_r232);
             FASTA in ../MAG_contigs/ (index.tsv)
  new      : work/checkm/<unit>/selected_<unit>.tsv (pass_50_10 == 1) and the GTDB-Tk summaries in gtdbtk/
The original 310 are dereplicated; the new MAGs are not (39 per-sample assemblies + 4 per-site sediment
co-assemblies recover the same organisms repeatedly). Species-level clusters make the two comparable:
skani 0.4.0 triangle over all genomes, then greedy centroid clustering (highest completeness - 5 x contamination
first) at ANI >= 95 % with the larger aligned fraction >= 50 % (GTDB's species thresholds).
Contig counts and N50 are computed from the FASTA for both sides.

Run: ~/Documents/py_venv/bin/python compare_vs_original.py  -> comparison/
"""
import glob, os, subprocess
from collections import Counter

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.dirname(HERE)
OUT = os.path.join(HERE, "comparison")
SKANI = os.path.expanduser("~/Documents/gtdbtk_venv/bin/skani")
SITE = lambda s: ("R2_restored" if s.startswith("R2_restored") else s.split("_")[0])
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
TAXA = [  # (label, rank, value) for the paper's taxa of interest
    ("Archaea", "domain", "Archaea"), ("Halobacteriota (haloarchaea + methanogens)", "phylum", "Halobacteriota"),
    ("Methanogens (Methanosarcinales)", "order", "Methanosarcinales"), ("Methanolobus", "genus", "Methanolobus"),
    ("Methanohalophilus", "genus", "Methanohalophilus"),
    ("Rhodobacteraceae", "family", "Rhodobacteraceae"), ("Roseovarius", "genus", "Roseovarius"),
    ("Marivita", "genus", "Marivita"), ("Albimonas", "genus", "Albimonas"),
    ("Halomonadaceae", "family", "Halomonadaceae"), ("Halomonas", "genus", "Halomonas"),
    ("Halanaerobiaceae", "family", "Halanaerobiaceae"), ("Halanaerobium", "genus", "Halanaerobium"),
    ("Desulfobacterota (sulfate reducers)", "phylum", "Desulfobacterota"),
]


def split_tax(c):
    d = dict(x.split("__", 1) for x in str(c).split(";") if "__" in x)
    return [d.get(r[0], "") for r in RANKS]


def fasta_stats(path):
    lens, n = [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lens.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lens.append(n)
    lens.sort(reverse=True)
    tot, acc = sum(lens), 0
    for L in lens:
        acc += L
        if acc >= tot / 2:
            return len(lens), tot, L
    return 0, 0, 0


def load():
    o = pd.read_csv(os.path.join(S, "MAG_QC/comparison_vs_2020.csv"))
    idx = pd.read_csv(os.path.join(S, "MAG_contigs/index.tsv"), sep="\t").set_index("assembly_name")
    old = pd.DataFrame({"genome": "old|" + o.assembly_name, "set": "original (2020, dereplicated)",
                        "sample": o.mag_key.str.split("|").str[0],
                        "completeness": o.Completeness_2026, "contamination": o.Contamination_2026,
                        "classification": o.gtdb_r232,
                        "fasta": [os.path.join(S, "MAG_contigs", idx.loc[n, "file"]) for n in o.assembly_name]})
    rows = []
    summ = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(HERE, "gtdbtk/gtdbtk.*.summary.tsv"))])
    summ = summ.set_index("user_genome").classification
    for f in sorted(glob.glob(os.path.join(HERE, "work/checkm/*/selected_*.tsv"))):
        unit = os.path.basename(os.path.dirname(f))
        t = pd.read_csv(f, sep="\t")
        for r in t[t.pass_50_10 == 1].itertuples():
            rows.append({"genome": "new|" + r.assembly_name, "set": "re-assembly (2026)", "sample": unit,
                         "completeness": r.completeness, "contamination": r.contamination,
                         "classification": summ[r.assembly_name],
                         "fasta": os.path.join(HERE, f"work/checkm/{unit}/fasta/{r.assembly_name}.fna")})
    d = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    d[RANKS] = pd.DataFrame(d.classification.map(split_tax).tolist(), index=d.index)
    d["site"] = d["sample"].map(lambda s: s.split("_sediment")[0] + " co-assembly" if s.endswith("_sediment") else SITE(s))
    d["quality"] = d.completeness - 5 * d.contamination
    d["mimag"] = ["high (>=90/<5)" if c >= 90 and x < 5 else "medium (>=50/<10)" if c >= 50 and x < 10 else "below MQ"
                  for c, x in zip(d.completeness, d.contamination)]
    st = [fasta_stats(p) for p in d.fasta]
    d["n_contigs"], d["size_bp"], d["n50"] = zip(*st)
    assert (d.size_bp > 0).all(), "missing FASTA"
    return d


def cluster(d):
    os.makedirs(OUT, exist_ok=True)
    lst = os.path.join(OUT, "skani_inputs.txt")
    d.fasta.to_csv(lst, index=False, header=False)
    edges = os.path.join(OUT, "skani_edges.tsv")
    if not os.path.exists(edges):
        subprocess.run([SKANI, "triangle", "-l", lst, "-t", "48", "-E", "--min-af", "30", "-o", edges], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    e = pd.read_csv(edges, sep="\t")
    byf = dict(zip(d.fasta, d.genome))
    e = e[(e.ANI >= 95) & (e[["Align_fraction_ref", "Align_fraction_query"]].max(axis=1) >= 50)]
    nb = {}
    for a, b in zip(e.Ref_file.map(byf), e.Query_file.map(byf)):
        nb.setdefault(a, set()).add(b); nb.setdefault(b, set()).add(a)
    rep = {}
    for g in d.sort_values(["quality", "completeness"], ascending=False).genome:   # greedy centroids
        if g in rep:
            continue
        rep[g] = g
        for h in nb.get(g, ()):
            rep.setdefault(h, g)
    d["species_cluster"] = d.genome.map(rep)
    return d


def main():
    d = cluster(load())
    d.drop(columns=["fasta"]).to_csv(os.path.join(OUT, "mags_old_vs_new.tsv"), sep="\t", index=False)
    old, new = d[d.genome.str.startswith("old|")], d[d.genome.str.startswith("new|")]
    cl = d.groupby("species_cluster").genome.agg(lambda s: ("old" if s.str.startswith("old|").any() else "") +
                                                 ("new" if s.str.startswith("new|").any() else ""))
    L = []
    w = L.append
    w("# Salterns MAGs: original study vs 2026 re-assembly\n")
    w("Both sets: same 39 metagenomes, KBase CheckM (lineage_wf) and GTDB-Tk 2.8.0 / GTDB R232. The original 310 "
      "MAGs are dereplicated; the 1,324 new MAGs are not (39 per-sample assemblies + 4 per-site sediment co-assemblies), "
      "so species clusters (skani, ANI >= 95 %, aligned fraction >= 50 %) are the fair comparison.\n")
    def q(x):
        return {"MAGs": len(x), "species clusters": x.species_cluster.nunique(),
                "high quality (>=90/<5)": int((x.mimag.str.startswith("high")).sum()),
                "medium quality": int((x.mimag.str.startswith("medium")).sum()),
                "below MQ (<50 or >=10)": int((x.mimag == "below MQ").sum()),
                "median completeness %": round(x.completeness.median(), 1),
                "median contamination %": round(x.contamination.median(), 2),
                "median contigs": int(x.n_contigs.median()), "median N50 (kb)": round(x.n50.median() / 1e3, 1),
                "median size (Mb)": round(x.size_bp.median() / 1e6, 2), "total size (Gb)": round(x.size_bp.sum() / 1e9, 2),
                "Bacteria / Archaea": f"{(x.domain == 'Bacteria').sum()} / {(x.domain == 'Archaea').sum()}",
                "with a GTDB species": int((x.species != "").sum()), "no GTDB genus": int((x.genus == "").sum())}
    tq = pd.DataFrame({"original (310, dereplicated)": q(old), "re-assembly (all 1,324)": q(new)})
    # best new MAG per species cluster (a dereplicated view of the new set)
    best_new = new.sort_values("quality", ascending=False).drop_duplicates("species_cluster")
    tq["re-assembly, best MAG per species cluster"] = pd.Series(q(best_new))
    w("## 1. Numbers and quality\n"); w(tq.to_markdown()); w("")
    w(f"Species clusters across both sets: {len(cl)}; in both: {(cl == 'oldnew').sum()}, original only: "
      f"{(cl == 'old').sum()}, re-assembly only: {(cl == 'new').sum()}.\n")
    # shared clusters: is the best new MAG better than the best old one?
    sh = cl[cl == "oldnew"].index
    bo = old.sort_values("quality", ascending=False).drop_duplicates("species_cluster").set_index("species_cluster")
    bn = best_new.set_index("species_cluster")
    cmp_ = pd.DataFrame({"dq": bn.loc[sh, "quality"] - bo.loc[sh, "quality"],
                         "dc": bn.loc[sh, "completeness"] - bo.loc[sh, "completeness"],
                         "dn50": bn.loc[sh, "n50"] / bo.loc[sh, "n50"]})
    w(f"For the {len(sh)} species clusters found by both, the best re-assembly MAG vs the best original MAG: "
      f"better quality score (completeness - 5 x contamination) in {(cmp_.dq > 0).sum()}, worse in {(cmp_.dq < 0).sum()}; "
      f"median change in completeness {cmp_.dc.median():+.1f} points; median N50 ratio {cmp_.dn50.median():.2f}x.\n")
    lost = old[old.species_cluster.map(cl) == "old"]
    w(f"Original species clusters NOT recovered by the re-assembly ({len(lost)} MAGs): "
      + (", ".join(f"{r.genus or r.family or r.order} ({r.completeness:.0f}/{r.contamination:.1f})"
                   for r in lost.sort_values(['phylum', 'genus']).itertuples()) or "none") + "\n")
    # what the lost original species were, and whether a fragment of them was recovered (ANI >= 95, AF < 50)
    e = pd.read_csv(os.path.join(OUT, "skani_edges.tsv"), sep="\t")
    byf = dict(zip(d.fasta, d.genome))
    near = {}
    for a, b, ani in zip(e.Ref_file.map(byf), e.Query_file.map(byf), e.ANI):
        for x, y in ((a, b), (b, a)):
            if x in set(lost.genome) and str(y).startswith("new|"):
                near[x] = max(near.get(x, 0), ani)
    good = lost[lost.mimag != "below MQ"]
    frag = sum(near.get(g, 0) >= 95 for g in good.genome)
    w(f"Of these {len(lost)} original MAGs, {(lost.mimag == 'below MQ').sum()} were below MIMAG medium quality in the "
      f"2026 re-QC. Of the {len(good)} medium- or high-quality ones, {frag} have a re-assembly MAG at >= 95 % ANI that "
      f"covers < 50 % of the genome (fragment recovered), and {sum(near.get(g, 0) < 80 for g in good.genome)} have no "
      f"re-assembly MAG above 80 % ANI. By sample type: "
      + ", ".join(f"{k} {v}" for k, v in Counter("water" if s.endswith("H2O") else "R2_restored sediment"
                                                  if s.startswith("R2_restored") else "sediment"
                                                  for s in good["sample"]).most_common()) + ".\n")
    # fragmentation: binner mix and contig-length floor
    d["binner"] = d.genome.str.extract(r"(maxbin|metabat|concoct)", expand=False)
    old, new = d[d.genome.str.startswith("old|")], d[d.genome.str.startswith("new|")]
    bt = pd.DataFrame({"original MAGs": old.binner.value_counts(), "re-assembly MAGs": new.binner.value_counts(),
                       "re-assembly median N50 (kb)": (new.groupby("binner").n50.median() / 1e3).round(1),
                       "re-assembly median contigs": new.groupby("binner").n_contigs.median()})
    w("**Why the re-assembly MAGs are more fragmented.** Every original MAG contig is >= 2,500 bp. In the re-assembly, "
      "MetaBAT2 and CONCOCT also used 2,500 bp, but MaxBin2 used 1,000 bp (pipeline.py), and DAS Tool kept far more "
      "MaxBin2 bins than in 2020. MaxBin2 bins have the lowest N50. MetaBAT2/CONCOCT bins are also below the original "
      "N50 (19.1 kb), so the metaSPAdes assemblies themselves are more fragmented than the 2020 ones; the cause is "
      "not established (read processing differs: Trimmomatic here vs JGI's pipeline in 2020).\n")
    w(bt.to_markdown()); w("")
    # taxa of interest
    rows = []
    for lab, rank, val in TAXA:
        o, n = old[old[rank] == val], new[new[rank] == val]
        rows.append({"taxon": lab, "original MAGs": len(o), "original species clusters": o.species_cluster.nunique(),
                     "re-assembly MAGs": len(n), "re-assembly species clusters": n.species_cluster.nunique(),
                     "new species clusters not in original": len(set(n.species_cluster) - set(o.species_cluster)
                                                                 - set(old.species_cluster)),
                     "original HQ / re-assembly HQ MAGs": f"{(o.mimag.str.startswith('high')).sum()} / "
                                                          f"{(n.mimag.str.startswith('high')).sum()}"})
    w("## 2. Taxa of interest (GTDB R232 on both sides)\n"); w(pd.DataFrame(rows).set_index("taxon").to_markdown()); w("")
    for lab, rank, val in [("Methanogens", "order", "Methanosarcinales"), ("Albimonas", "genus", "Albimonas"),
                           ("Halomonas", "genus", "Halomonas")]:
        x = d[d[rank] == val].sort_values(["set", "species_cluster"])
        w(f"**{lab}:** " + "; ".join(f"{r.set.split(' ')[0]} {r.genome.split('|')[1][:60]} ({r.sample}, "
                                       f"{r.completeness:.0f}/{r.contamination:.1f}, {r.species or r.genus})"
                                       for r in x.itertuples()) + "\n")
    # phyla and sites
    ph = pd.DataFrame({"original MAGs": old.phylum.value_counts(), "re-assembly MAGs": new.phylum.value_counts(),
                       "original species": old.groupby("phylum").species_cluster.nunique(),
                       "re-assembly species": new.groupby("phylum").species_cluster.nunique()}).fillna(0).astype(int)
    w("## 3. Phyla\n"); w(ph.sort_values("re-assembly species", ascending=False).head(20).to_markdown()); w("")
    si = pd.DataFrame({"original MAGs": old.site.value_counts(), "re-assembly MAGs": new.site.value_counts()}).fillna(0).astype(int)
    w("## 4. By site (original = dereplicated, so counts are after cross-sample dereplication)\n"); w(si.to_markdown()); w("")
    w("Files: comparison/mags_old_vs_new.tsv (per MAG: set, sample, CheckM, R232 lineage, contigs, N50, species "
      "cluster); comparison/skani_edges.tsv.\n")
    open(os.path.join(OUT, "README.md"), "w").write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
