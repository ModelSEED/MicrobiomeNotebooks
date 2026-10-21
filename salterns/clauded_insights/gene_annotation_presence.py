#!/usr/bin/env python
"""
Recover reaction presence from RAST GENE ANNOTATION (not model reactions), using the
local RASTified genome JSONs in jplfaria/Cliff MAGs/RASTified_MAG_json/.  This captures
reactions the base ModelSEED models never built (e.g. rxn00270, rxn05011) and corrects
under/over-counting.

For each reaction we match the RAST 'functions' strings (EC number + enzyme name); we
parse ONLY the functions field (raw text regex) so protein-sequence substrings can't
create false positives.  Presence is deduplicated to the abundance-table MAGs, then the
per-reaction read abundance (% of mapped community) is recomputed and compared to the
model-based version.

Outputs (suffix _geneannot):
  reaction_presence_model_vs_geneannot.csv
  read_abundance_by_reaction_totMAGs_geneannot.csv
  read_abundance_by_reaction_totMAGs_byType_geneannot.csv
"""
from __future__ import annotations
import glob, os, re
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
RAST = REPO / "jplfaria/Cliff MAGs/RASTified_MAG_json"

MEPN = ["rxn26456", "rxn26457", "rxn43120", "rxn26447", "rxn04651"]
BET  = ["rxn00270", "rxn05010", "rxn05011", "rxn17220", "rxn24612", "rxn24613", "rxn03127"]

# RAST function substrings (lower-case) that indicate each reaction's enzyme
PATTERNS = {
    "rxn26456": ["subunit phng", "phng "],
    "rxn26457": ["subunit phnh", "phnh "],
    "rxn26447": ["subunit phni", "phni "],
    "rxn43120": ["c-p lyase", "ec 4.7.1.1"],
    "rxn04651": ["phnn "],
    "rxn00270": ["glycine n-methyltransferase"],            # gly -> sarcosine (GSMT)
    "rxn05010": ["sarcosine n-methyltransferase"],          # sarcosine -> DMG (SDMT)
    "rxn05011": ["dimethylglycine n-methyltransferase"],    # DMG -> betaine
    "rxn17220": ["betaine reductase component b"],          # betaine reductase (EC 1.21.4.4)
    "rxn24612": ["trimethylamine methyltransferase"],       # mttB
    "rxn24613": ["trimethylamine methyltransferase"],
    "rxn03127": ["methyl-coenzyme m reductase", "methyl coenzyme m reductase"],
}

def norm(x):
    x = x.strip().strip('"').lower()
    for suf in [".contigs__.rast.json", ".contigs__.rast", ".contigs.fa_assembly.rast.json",
                ".contigs.fa_assembly.rast", ".contigs.fa", ".contigs", ".rast.json",
                ".json", "__.rast", "__", ".rast"]:
        x = x.replace(suf, "")
    return x.replace("_bins.", "_bins_").rstrip("_.")

# ---- gather each MAG's RAST function strings (dedup by normalised name) ------------
fn_re = re.compile(r'"functions"\s*:\s*\[(.*?)\]', re.S)
mag_functions: dict[str, str] = {}
for fp in glob.glob(str(RAST / "*.json")):
    mag = norm(os.path.basename(fp))
    text = Path(fp).read_text(errors="ignore")
    funcs = " ".join(fn_re.findall(text)).lower()
    mag_functions[mag] = mag_functions.get(mag, "") + " " + funcs
print(f"parsed RAST annotations for {len(mag_functions)} distinct MAGs")

# ---- abundance table + dedup universe ---------------------------------------------
rel = pd.read_csv(REPO / "data/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
rel.index = [norm(i) for i in rel.index]
SAMPLES = list(rel.columns)
idx = set(rel.index)

# ---- gene-annotation presence per reaction ----------------------------------------
def present(rxn):
    pats = PATTERNS[rxn]
    return {m for m, fx in mag_functions.items() if any(p in fx for p in pats)}

geneannot = {r: present(r) & idx for r in MEPN + BET}

# model-based n_MAGs (from the already-generated file) for comparison
model = pd.read_csv(HERE / "read_abundance_by_reaction_totMAGs.csv", index_col=0)["n_MAGs"]

cmp = pd.DataFrame({
    "pathway": ["MePn"]*len(MEPN) + ["Bet"]*len(BET),
    "n_MAGs_model": [int(model.get(r, 0)) for r in MEPN + BET],
    "n_MAGs_geneannot": [len(geneannot[r]) for r in MEPN + BET],
}, index=MEPN + BET)
cmp.to_csv(HERE / "reaction_presence_model_vs_geneannot.csv")
print("\n== reaction presence: model-based vs gene-annotation (# abundance MAGs) ==")
print(cmp.to_string())

# ---- recompute read abundance (%) with gene-annotation presence -------------------
def abund(mags):
    m = [x for x in mags if x in idx]
    return rel.loc[m, SAMPLES].sum()*100 if m else pd.Series(0.0, index=SAMPLES)

rows = {r: abund(geneannot[r]) for r in MEPN + BET}
tbl = pd.DataFrame(rows).T
mepn_u = set().union(*[geneannot[r] for r in MEPN]); bet_u = set().union(*[geneannot[r] for r in BET])
tbl.loc["MePn_pathway_union"] = abund(mepn_u); tbl.loc["Bet_pathway_union"] = abund(bet_u)
tbl.insert(0, "n_MAGs", [len(geneannot[r]) for r in MEPN+BET] + [len(mepn_u & idx), len(bet_u & idx)])
tbl.insert(0, "pathway", ["MePn"]*len(MEPN)+["Bet"]*len(BET)+["MePn","Bet"])
tbl["MEAN_across_samples"] = tbl[SAMPLES].mean(axis=1)
tbl.to_csv(HERE / "read_abundance_by_reaction_totMAGs_geneannot.csv")

# by sample type (D1/D2 within each site type)
def stype(c):
    depth = "D1" if c.endswith("_D1") else "D2" if c.endswith("_D2") else None
    if depth is None: return None
    if c.startswith("R2A"): t="Reference_R2A"
    elif c.startswith("SF2"): t="Restored_SF2"
    elif c.startswith(("R1_","R2_")): t="Unrestored_R1R2"
    else: return None
    return f"{t}_{depth}"
grp={}
for c in SAMPLES:
    g=stype(c)
    if g: grp.setdefault(g,[]).append(c)
bt=pd.DataFrame({g: tbl[cols].mean(axis=1) for g,cols in grp.items()})
bt.insert(0,"n_MAGs",tbl["n_MAGs"]); bt.insert(0,"pathway",tbl["pathway"])
order=["pathway","n_MAGs","Unrestored_R1R2_D1","Unrestored_R1R2_D2","Restored_SF2_D1",
       "Restored_SF2_D2","Reference_R2A_D1","Reference_R2A_D2"]
bt[order].to_csv(HERE / "read_abundance_by_reaction_totMAGs_byType_geneannot.csv")
print("\n== by-type (gene-annotation), signature reactions ==")
print(bt[order].loc[["rxn43120","rxn17220","rxn05011","rxn03127","MePn_pathway_union","Bet_pathway_union"]].round(2).to_string())
print("\nSaved *_geneannot.csv (presence comparison, totMAGs, byType)")
