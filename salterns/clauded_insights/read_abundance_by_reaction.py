#!/usr/bin/env python
"""
Read→MAG mapping → per-reaction read abundance for the MePn and betaine pathways.

read->MAG : data/Cliff_310MAG_relabund.tsv  (the reconciled CoverM read-recruitment
            result: 311 MAGs x 24 samples, relative abundance among mapped reads,
            columns sum to 1). This is the product of the bam_bins mapping in ./reads,
            keyed by MAG name (the raw BAMs use un-reconcilable KBH_/sha256 IDs).
MAG->rxn  : (a) methane_producing_mags.csv — CURATED presence of the 5 MePn reactions
            (incl. rxn43120, which base RAST models lack); (b) every local MAG model
            (models/, models/MePn_consumers/, betaine_models/final/) — union of
            modelreactions per MAG, which supplies the curated betaine reactions
            (rxn17220 lives in betaine_models/final; rxn03127 in the methanogen model).
sum       : for each reaction, sum the relative abundance of MAGs carrying it, per
            sample; pathway score = union of MAGs carrying ANY pathway reaction
            (avoids double-counting) plus the literal per-reaction sum.

Reaction IDs validated against ../../ModelSEEDDatabase/Objects/Reactions.json.gz.
Run: ~/Documents/py_venv/bin/python read_abundance_by_reaction.py
"""
from __future__ import annotations
import csv, glob, gzip, json, os, re
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[1]
MSDB = Path("/home/freiburger/Documents/ModelSEEDDatabase/Objects/Reactions.json.gz")

MEPN = ["rxn26456", "rxn26457", "rxn43120", "rxn26447", "rxn04651"]
BET  = ["rxn00270", "rxn05010", "rxn05011", "rxn17220", "rxn24612", "rxn24613", "rxn03127"]

# ---- validate reaction IDs against ModelSEED DB ------------------------------------
_db = json.load(gzip.open(MSDB))
_db = {r["id"]: r for r in _db} if isinstance(_db, list) else _db
print("== reaction ID validation (ModelSEEDDatabase) ==")
for r in MEPN + BET:
    d = _db.get(r)
    print(f"  {r:10s} {'OK ' if d else 'NOT IN DB'} {d.get('name','')[:52] if d else '(custom/curated or malformed)'}")

# ---- name normalisation (join model / csv names to the relabund index) ------------
def norm(x: str) -> str:
    x = x.strip().strip('"').lower()
    for suf in [".contigs__.rast.json", ".contigs__.rast", ".contigs.fa_assembly.rast.json",
                ".contigs.fa_assembly.rast", ".contigs.fa", ".contigs", ".rast.json",
                "_betaine.json", ".json", "__.rast", "__", ".rast"]:
        x = x.replace(suf, "")
    x = x.replace("_bins.", "_bins_")            # unify the bins.<binner> / bins_<binner> separator
    return x.rstrip("_.")

# ---- read->MAG abundance ----------------------------------------------------------
rel = pd.read_csv(REPO / "data/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
rel.index = [norm(i) for i in rel.index]
SAMPLES = list(rel.columns)
print(f"\nread->MAG abundance: {rel.shape[0]} MAGs x {len(SAMPLES)} samples (cols sum to {rel[SAMPLES[0]].sum():.2f})")

# ---- MAG -> reaction presence (from RAST GENE ANNOTATION) --------------------------
# Presence is taken from the RAST 'functions' strings in the local RASTified genome
# JSONs, matched by EC number + enzyme name.  This captures reactions the base ModelSEED
# models never built (e.g. rxn00270, rxn05011) and doesn't depend on model curation.
# Only the functions field is parsed, so protein-sequence substrings can't create
# false positives.  (rxn24612/rxn24613 both annotate as the mttB TMA-methyltransferase
# family and cannot be separated at the gene level.)
RAST = REPO / "jplfaria/Cliff MAGs/RASTified_MAG_json"
PATTERNS = {
    "rxn26456": ["subunit phng", "phng "],
    "rxn26457": ["subunit phnh", "phnh "],
    "rxn26447": ["subunit phni", "phni "],
    "rxn43120": ["c-p lyase", "ec 4.7.1.1"],
    "rxn04651": ["phnn "],
    "rxn00270": ["glycine n-methyltransferase"],
    "rxn05010": ["sarcosine n-methyltransferase"],
    "rxn05011": ["dimethylglycine n-methyltransferase"],
    "rxn17220": ["betaine reductase component b"],
    "rxn24612": ["trimethylamine methyltransferase"],
    "rxn24613": ["trimethylamine methyltransferase"],
    "rxn03127": ["methyl-coenzyme m reductase", "methyl coenzyme m reductase"],
}
_fn_re = re.compile(r'"functions"\s*:\s*\[(.*?)\]', re.S)
mag_functions: dict[str, str] = {}
for fp in glob.glob(str(RAST / "*.json")):
    mag = norm(os.path.basename(fp))
    funcs = " ".join(_fn_re.findall(Path(fp).read_text(errors="ignore"))).lower()
    mag_functions[mag] = mag_functions.get(mag, "") + " " + funcs
rxn_mags = {r: {m for m, fx in mag_functions.items() if any(p in fx for p in PATTERNS[r])}
            for r in MEPN + BET}
print(f"parsed RAST gene annotation for {len(mag_functions)} distinct MAGs")

# ---- join diagnostics -------------------------------------------------------------
idx = set(rel.index)
print("\n== per-reaction MAG support (in curated table/models  ->  matched to abundance table) ==")
for grp, rxns in [("MePn", MEPN), ("Bet", BET)]:
    print(f"  [{grp}]")
    for r in rxns:
        mm = rxn_mags[r]
        matched = mm & idx
        print(f"    {r:10s} MAGs={len(mm):3d}  matched_to_abundance={len(matched):3d}"
              + ("   <-- malformed/absent" if not mm else ""))

# ---- sum read abundance per reaction, per sample ----------------------------------
def abundance_of(mags):
    m = [x for x in mags if x in idx]
    return rel.loc[m, SAMPLES].sum() if m else pd.Series(0.0, index=SAMPLES)

per_rxn = {r: abundance_of(rxn_mags[r]) for r in MEPN + BET}
rxn_table = pd.DataFrame(per_rxn).T          # reactions x samples
rxn_table["MEAN_across_samples"] = rxn_table[SAMPLES].mean(axis=1)
rxn_table["n_MAGs"] = [len(rxn_mags[r] & idx) for r in rxn_table.index]

# pathway score = union of MAGs carrying ANY pathway reaction (no double counting)
mepn_union = set().union(*[rxn_mags[r] for r in MEPN])
bet_union  = set().union(*[rxn_mags[r] for r in BET])
path = pd.DataFrame({"MePn_pathway_union": abundance_of(mepn_union),
                     "Bet_pathway_union":  abundance_of(bet_union)}).T
path["MEAN_across_samples"] = path[SAMPLES].mean(axis=1)
path["n_MAGs"] = [len(mepn_union & idx), len(bet_union & idx)]
# literal "sum of the reactions" (double-counts multi-reaction MAGs)
path.loc["MePn_sum_of_reactions"] = rxn_table.loc[MEPN, SAMPLES].sum()
path.loc["Bet_sum_of_reactions"]  = rxn_table.loc[BET,  SAMPLES].sum()
path.loc["MePn_sum_of_reactions","MEAN_across_samples"] = path.loc["MePn_sum_of_reactions",SAMPLES].mean()
path.loc["Bet_sum_of_reactions","MEAN_across_samples"]  = path.loc["Bet_sum_of_reactions",SAMPLES].mean()

out = pd.concat([rxn_table, path])
# classify each reaction / pathway row as MePn or Bet
out.insert(0, "pathway", ["MePn" if (r in MEPN or str(r).startswith("MePn"))
                          else "Bet" if (r in BET or str(r).startswith("Bet")) else ""
                          for r in out.index])
out.to_csv(REPO / "clauded_insights/read_abundance_by_reaction.csv")

# analogous file: total % abundance (the summed relative abundance of the possessing
# MAGs, expressed as a percentage of the mapped community) -- values x100.
abund_cols = SAMPLES + ["MEAN_across_samples"]
pct = out.copy()
pct[abund_cols] = pct[abund_cols] * 100.0
pct.to_csv(REPO / "clauded_insights/read_abundance_by_reaction_totMAGs.csv")

pd.set_option("display.width", 200, "display.max_columns", 30)
print("\n== per-reaction read abundance (fraction of mapped community) ==")
print(rxn_table[["n_MAGs", "MEAN_across_samples"]].round(4).to_string())
print("\n== pathway read abundance ==")
print(path[["n_MAGs", "MEAN_across_samples"]].round(4).to_string())
print("\n== per-sample pathway union scores (first 8 samples) ==")
print(path.loc[["MePn_pathway_union","Bet_pathway_union"], SAMPLES[:8]].round(4).to_string())
print(f"\nSaved read_abundance_by_reaction.csv  (fractions, reactions+pathways x {len(SAMPLES)} samples)")
print("Saved read_abundance_by_reaction_totMAGs.csv  (same, as total %% abundance of possessing MAGs)")
print("\n== total %% abundance of MAGs possessing each reaction (mean across samples) ==")
print(pct.loc[MEPN + BET, ["n_MAGs", "MEAN_across_samples"]].round(2).to_string())
