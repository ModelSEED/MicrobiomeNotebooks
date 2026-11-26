#!/usr/bin/env python
"""Table S4 (2026 revision): the 39 metagenomes used (24 sediment, 15 water), with public accessions and the read and
assembly numbers of this study. Sources: reassembly_2026/manifest.tsv (SRA/IMG/GOLD/BioProject/BioSample, spots,
bases), exported_img_data.tsv (IMG assemblies: original and SPAdes version 2), abundance_final/qc_read_totals.tsv
(BBDuk-trimmed reads), manuscript_numbers_2026-11-19.json (second-pass assembly and MAGs per sample), and the CH4 flux
of the sediment cores (data/Cliff_Sample_Metadata_BGC_NMR.csv).
Run: ~/Documents/py_venv/bin/python make_table_s4.py -> ../../Table S4 - metagenomes, accessions and re-assembly (2026).xlsx
"""
import json
from pathlib import Path

import pandas as pd

H = Path(__file__).resolve().parent
R = H.parent
S = R.parent
m = pd.read_csv(R / "manifest.tsv", sep="\t")
img = pd.read_csv(S / "exported_img_data.tsv", sep="\t")
img["sample"] = img["Genome Name / Sample Name"].str.extract(r"Salt Pond MetaG (\S+?)_MG")[0]
v2 = img[img["Genome Name / Sample Name"].str.contains("SPAdes")].set_index("sample")["taxon_oid"]
q = pd.read_csv(R / "abundance_final/qc_read_totals.tsv", sep="\t", header=None, names=["sample", "qc_reads"]).set_index("sample").qc_reads
N = json.load(open(R / "manuscript_numbers_2026-11-19.json"))
rows = pd.DataFrame(N["rows"]).set_index("sample")
meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
rest = {"R1": "unrestored saltern", "R2": "unrestored saltern", "SF2": "restored saltern", "R2A": "reference wetland"}
site_type = m.apply(lambda r: rest.get(r["site"], "restored pond (additional sample)") if "restored" not in r["sample"] else "restored pond R2 (additional sample)", axis=1)
out = pd.DataFrame({
    "Sample": m["sample"], "Type": m["type"], "Site": m["site"], "Site type": site_type, "Layer": m["layer"],
    "NCBI BioProject": m["bioproject"], "NCBI BioSample": m["biosample"], "SRA run": m["srr"], "GOLD": m["gold"],
    "IMG/M taxon OID (JGI assembly)": m["img"].map(lambda x: "" if pd.isna(x) else str(int(x))), "IMG/M taxon OID (JGI SPAdes v2)": m["sample"].map(v2).map(lambda x: "" if pd.isna(x) else str(int(x))),
    "Read pairs (SRA spots)": m["spots"], "Bases (SRA)": m["bases"], "Reads after trimming (this study)": m["sample"].map(q),
    "Second-pass MAGs kept (high-quality)": m["sample"].map(lambda s: f"{rows.loc[s, 'new_kept']} ({rows.loc[s, 'new_hq']})" if s in rows.index and pd.notna(rows.loc[s, "new_kept"]) else ""),
    "CH4 flux (umol m-2 d-1)": m["sample"].map(lambda s: meta.loc[s, "CH4_umol_m2_d"] if s in meta.index else None),
})
assert len(out) == 39 and out["SRA run"].notna().all()
assert (out.Type == "Sediment").sum() == 26 and (out.Type == "Water").sum() == 13
out.to_excel(S / "Table S4 - metagenomes, accessions and re-assembly (2026).xlsx", index=False)
print(out.head(3).T.to_string()); print(out.Type.value_counts().to_dict())
