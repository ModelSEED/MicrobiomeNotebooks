#!/usr/bin/env python3
"""One DRAM2 gene table with BOTH coding and non-coding genes, plus a per-MAG summary (per Andrew, 2026-11-18:
"make sure that both coding and non-coding genes are annotated").

DRAM2 v2.0.0-beta17 annotates the two kinds in separate outputs:
  coding      results/raw-annotations.all.tsv.gz      Prodigal CDS x KOfam, dbCAN, MEROPS, VOG, CAMPER, CANT-HYD,
                                                        methyl, sulfur, FeGenie (one row per CDS)
  non-coding  results/noncoding/combined_trna_scan.tsv.gz tRNAscan-SE tRNAs (amino acid, anticodon, score)
              results/noncoding/combined_rrna_scan.tsv barrnap 5S/16S/23S rRNAs
(DRAM2 calls no other non-coding RNA: no tmRNA, Rfam ncRNA or CRISPR. The KBDL Bakta builds of the same MAGs carry
those; see kbdl_build.)
Outputs:
  results/dram2_genes_coding_noncoding.tsv.gz  mag, gene_type (CDS|tRNA|rRNA), gene_id, contig, start, end, strand,
      product (CDS: the KOfam description if any, else dbCAN/MEROPS/CAMPER/...; tRNA: "tRNA-<aa> (<anticodon>)";
      rRNA: barrnap product), kofam_id, kofam_passes_ko_threshold, any_db_hit, score  (large; not in git)
  results/dram2_per_mag_genes.tsv  per MAG: CDS, CDS with any database hit, with a KO, with a threshold-passing KO,
      tRNA, distinct tRNA amino acids, 5S/16S/23S rRNA
Usage: ~/Documents/py_venv/bin/python -I merge_coding_noncoding.py
"""
import os
import pandas as pd

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
HIT_COLS = ["kofam_id", "dbcan_id", "merops_id", "camper_id", "vogdb_id", "cant_hyd_id", "sulfur_id", "fegenie_id",
            "methyl_id"]
DESC = ["kofam_description", "dbcan_description", "merops_description", "camper_description", "cant_hyd_description",
        "sulfur_id", "fegenie_id", "methyl_id", "vogdb_id"]

cds = pd.read_csv(os.path.join(R, "raw-annotations.all.tsv.gz"), sep="\t", low_memory=False,
                  usecols=["query_id", "input_fasta", "mag", "start_position", "stop_position", "strandedness",
                           "kofam_passes_ko_threshold"] + sorted(set(HIT_COLS + DESC)))
f2m = cds.drop_duplicates("input_fasta").set_index("input_fasta").mag
assert f2m.is_unique and len(f2m) == 1184, len(f2m)
prod = cds[DESC[0]].copy()
for c in DESC[1:]:
    prod = prod.fillna(cds[c])
coding = pd.DataFrame({
    "mag": cds.mag, "gene_type": "CDS", "gene_id": cds.query_id, "contig": cds.query_id.str.rsplit("_", n=1).str[0],
    "start": cds.start_position, "end": cds.stop_position,
    "strand": cds.strandedness.map({1: "+", -1: "-", "1": "+", "-1": "-"}), "product": prod,
    "kofam_id": cds.kofam_id, "kofam_passes_ko_threshold": cds.kofam_passes_ko_threshold,
    "any_db_hit": cds[HIT_COLS].notna().any(axis=1), "score": pd.NA})

t = pd.read_csv(os.path.join(R, "noncoding", "combined_trna_scan.tsv.gz"), sep="\t", dtype=str)
t["query_id"] = t.query_id.str.strip()
assert set(t.input_fasta) <= set(f2m.index)
b, e = t.begin.astype(int), t.end.astype(int)
trna = pd.DataFrame({
    "mag": t.input_fasta.map(f2m), "gene_type": "tRNA",
    "gene_id": t.query_id + "_tRNA" + t["tRNA #"].str.strip(), "contig": t.query_id,
    "start": b.where(b <= e, e), "end": e.where(b <= e, b), "strand": (b <= e).map({True: "+", False: "-"}),
    "product": "tRNA-" + t.type.str.strip() + " (" + t.codon.str.strip() + ")",
    "kofam_id": pd.NA, "kofam_passes_ko_threshold": pd.NA, "any_db_hit": pd.NA, "score": t.score})

r = pd.read_csv(os.path.join(R, "noncoding", "combined_rrna_scan.tsv"), sep="\t", dtype=str)
assert set(r.input_fasta) <= set(f2m.index)
rprod = r.note.str.extract(r"product=([^;]+)")[0]
rrna = pd.DataFrame({
    "mag": r.input_fasta.map(f2m), "gene_type": "rRNA",
    "gene_id": r.query_id + "_rRNA_" + r.begin + "_" + r.end, "contig": r.query_id,
    "start": r.begin.astype(int), "end": r.end.astype(int), "strand": r.strand, "product": rprod,
    "kofam_id": pd.NA, "kofam_passes_ko_threshold": pd.NA, "any_db_hit": pd.NA, "score": r["e-value"]})

allg = pd.concat([coding, trna, rrna], ignore_index=True).sort_values(["mag", "contig", "start"], kind="stable")
allg.to_csv(os.path.join(R, "dram2_genes_coding_noncoding.tsv.gz"), sep="\t", index=False, compression="gzip")

g = coding.groupby("mag")
summ = pd.DataFrame({
    "CDS": g.size(), "CDS_any_db_hit": g.any_db_hit.sum(), "CDS_with_KO": g.kofam_id.count(),
    "CDS_KO_passing_threshold": g.kofam_passes_ko_threshold.apply(lambda s: (s.astype(str) == "True").sum()),
    "tRNA": trna.groupby("mag").size(),
    "tRNA_amino_acids": trna.assign(aa=t.type.str.strip()).query("aa not in ['Undet', 'Pseudo', 'Sup']")
                            .groupby("mag").aa.nunique(),
})
for kind in ("5S", "16S", "23S"):
    summ[f"rRNA_{kind}"] = rrna[rrna["product"].str.startswith(kind, na=False)].groupby("mag").size()
summ = summ.reindex(sorted(f2m)).fillna(0).astype(int)
summ.index.name = "mag"
summ.to_csv(os.path.join(R, "dram2_per_mag_genes.tsv"), sep="\t")
tot = summ.sum()
print(f"{len(summ)} MAGs; genes: CDS {tot.CDS:,} (any DB hit {tot.CDS_any_db_hit:,}, KO {tot.CDS_with_KO:,}, "
      f"KO passing threshold {tot.CDS_KO_passing_threshold:,}); tRNA {tot.tRNA:,}; rRNA 5S {tot.rRNA_5S:,}, "
      f"16S {tot.rRNA_16S:,}, 23S {tot.rRNA_23S:,}")
print("MAGs with >= 1:", {c: int((summ[c] > 0).sum()) for c in summ.columns})
print("median tRNA amino acids per MAG:", summ.tRNA_amino_acids.median())
