"""
Load individual organism models from models/, add betaine reductase (rxn17220)
and glycine reductase (rxn07207) reactions with required metabolites,
and export to betaine_models/final/.
"""
import json
import os

MODELS_DIR = os.path.join(os.path.dirname(__file__), "models")
OUT_DIR = os.path.join(os.path.dirname(__file__), "betaine_models", "final")

# Betaine-consuming bins mapped to their sample-specific model files
BETAINE_BINS = {
    "concoct_out.9 (R1_A_D1)": "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs__.RAST.json",
    "metabat.47 (R1_A_D1)": "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs__.RAST.json",
    "metabat.51 (R1_A_D1)": "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs__.RAST.json",
    "metabat.27 (R1_A_D2)": "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs__.RAST.json",
    "metabat.45 (R1_A_D2)": "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs__.RAST.json",
    "metabat.48 (R1_A_D2)": "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs__.RAST.json",
    "metabat.50 (R1_A_D2)": "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs__.RAST.json",
    "concoct_out.59 (R1_B_D1)": "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs__.RAST.json",
    "metabat.58 (R1_B_D1)": "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs__.RAST.json",
    "concoct_out.73 (R1_B_D2)": "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs__.RAST.json",
    "maxbin.047 (R1_C_D1)": "Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs__.RAST.json",
    "concoct_out.85 (R1_C_D2)": "Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs__.RAST.json",
    "metabat.40 (R1_C_D2)": "Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs__.RAST.json",
    "concoct_out.98 (R2_A_D1)": "Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs__.RAST.json",
    "metabat.52 (R2_B_D2)": "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs__.RAST.json",
}

# Reactions to add
# rxn17220: Betaine reductase (EC 1.21.4.4)
#   BET + Red-Thioredoxin + Pi + 2H+ -> TMA + Acetylphosphate + Ox-Thioredoxin + H2O
RXN17220 = {
    "id": "rxn17220_c0",
    "name": "Betaine reductase",
    "metabolites": {
        "cpd00001_c0": 1,    # H2O
        "cpd00009_c0": -1,   # Pi
        "cpd00067_c0": -2,   # H+
        "cpd00196_c0": 1,    # Acetylphosphate
        "cpd00441_e0": 1,    # Trimethylamine (extracellular)
        "cpd00540_e0": -1,   # Glycine betaine (extracellular)
        "cpd27735_c0": 1,    # Oxidized Thioredoxin
        "cpd28060_c0": -1,   # Reduced Thioredoxin
    },
    "lower_bound": 0,
    "upper_bound": 1000,
    "gene_reaction_rule": "",
    "objective_coefficient": 1.0,
    "subsystem": "Betaine metabolism",
}

# rxn07207: Glycine reductase (EC 1.21.4.2)
#   Glycine + Red-Thioredoxin + Pi + 2H+ -> Acetylphosphate + NH3 + Ox-Thioredoxin
RXN07207 = {
    "id": "rxn07207_c0",
    "name": "Glycine reductase",
    "metabolites": {
        "cpd00033_c0": -1,   # Glycine
        "cpd28060_c0": -1,   # Reduced Thioredoxin
        "cpd00009_c0": -1,   # Pi
        "cpd00067_c0": -2,   # H+
        "cpd00196_c0": 1,    # Acetylphosphate
        "cpd00013_c0": 1,    # NH3
        "cpd27735_c0": 1,    # Oxidized Thioredoxin
    },
    "lower_bound": 0,
    "upper_bound": 1000,
    "gene_reaction_rule": "",
    "objective_coefficient": 1.0,
    "subsystem": "Betaine metabolism",
}

# Metabolite that may be absent from individual models
TMA_E0 = {
    "id": "cpd00441_e0",
    "name": "(CH3)3N [e0]",
    "compartment": "e0",
    "charge": 1,
    "formula": "C3H10N",
    "notes": {"modelseed_template_id": "cpd00441_e"},
}

# Exchange reaction for TMA
EX_TMA = {
    "id": "EX_cpd00441_e0",
    "name": "Exchange for (CH3)3N [e0]",
    "metabolites": {"cpd00441_e0": -1},
    "lower_bound": 0,
    "upper_bound": 1000,
    "gene_reaction_rule": "",
    "notes": {},
    "annotation": {},
}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    for label, fname in BETAINE_BINS.items():
        filepath = os.path.join(MODELS_DIR, fname)
        with open(filepath) as fh:
            model = json.load(fh)

        met_ids = {m["id"] for m in model["metabolites"]}
        rxn_ids = {r["id"] for r in model["reactions"]}
        added = []

        # add TMA extracellular metabolite if missing
        if "cpd00441_e0" not in met_ids:
            model["metabolites"].append(TMA_E0)
            added.append("cpd00441_e0")

        # add betaine reductase
        if "rxn17220_c0" not in rxn_ids:
            model["reactions"].append(RXN17220)
            added.append("rxn17220_c0")

        # add glycine reductase
        if "rxn07207_c0" not in rxn_ids:
            model["reactions"].append(RXN07207)
            added.append("rxn07207_c0")

        # add TMA exchange reaction
        if "EX_cpd00441_e0" not in rxn_ids:
            model["reactions"].append(EX_TMA)
            added.append("EX_cpd00441_e0")

        outpath = os.path.join(OUT_DIR, fname)
        with open(outpath, "w") as fh:
            json.dump(model, fh)

        size_kb = os.path.getsize(outpath) / 1024
        print(f"{label:>30}: {', '.join(added) if added else 'no changes'} ({size_kb:.0f} KB)")

    print(f"\nDone. {len(BETAINE_BINS)} models saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
