"""B11. Hand C1 balance from the biomass stoichiometry of the original Roseovarius.9 model (bio1 + Protein/mRNA/DNA
synthesis pseudo-reactions), independent of any FBA: glycine made by SHMT (Ser -> Gly + CH2-THF) releases one C1;
the biomass consumes C1 for purines (2 formyl each, de novo), Met, His (AICAR recycling), dTMP, pantoate, THF-C1
cofactors and SAM methylations. If the model must make Thr (and hence Ile) from glycine (no homoserine kinase),
the C1 surplus that something must dispose of is predicted below and compared with the FBA MePn flux per unit
biomass (0.4697, with 1 MePn -> 1 drained 5'-dA -> 1 purine rebuilt de novo = net -1 C1)."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import cobra
cobra.Configuration().solver = "glpk"
inv = {iid(m): m for m in model_ids()}
for t in ["Roseovarius.9", "Albimonas.1"]:
    m = cobra.io.load_json_model(str(MEPN_MODELS / f"{inv[t]}.json"))
    b = m.reactions.bio1
    c = {x.id: -v for x, v in b.metabolites.items()}
    prot, rna, dna = c["cpd11463_c0"], c["cpd11462_c0"], c["cpd11461_c0"]
    P = {x.id: -v for x, v in m.reactions.rxn05296_c0.metabolites.items()}
    R = {x.id: -v for x, v in m.reactions.rxn05295_c0.metabolites.items()}
    D = {x.id: -v for x, v in m.reactions.rxn05294_c0.metabolites.items()}
    cof = c.get("cpd00003_c0", 0)  # all small cofactors share one coefficient in this template biomass
    gly_prot = prot * P["cpd00033_c0"]; thr = prot * P["cpd00161_c0"]; ile = prot * P["cpd00322_c0"]
    met = prot * P["cpd00060_c0"]; his = prot * P["cpd00119_c0"]
    pur_rna = rna * (R["cpd00002_c0"] + R["cpd00038_c0"]); pur_dna = dna * (D["cpd00115_c0"] + D["cpd00241_c0"])
    pur_cof = cof * 5        # NAD, NADP, CoA, FAD, SAM adenine
    purines = pur_rna + pur_dna + pur_cof
    dtmp = dna * D["cpd00357_c0"]
    heme_gly = 8 * (c.get("cpd00028_c0", 0) + c.get("cpd00557_c0", 0))   # ALA synthase: 8 Gly per tetrapyrrole
    gsh_gly = c.get("cpd00042_c0", 0)
    gly_noThr = gly_prot + purines + heme_gly + gsh_gly
    gly_withThr = gly_noThr + thr + ile
    c1_pur = 2 * purines
    c1_met = met + cof                    # protein Met + SAM's Met
    c1_his = his                          # AICAR from His synthesis needs one 10-formyl-THF to return to IMP
    c1_misc = dtmp + 2 * cof + 2 * cof + (3 + 1 + 2) * cof   # dTMP; pantoate (CoA, ACP); 10-CHO-/5-CH3-THF; SAM methyls UQ8(3), MK8(1), siroheme(2)
    c1_dem = c1_pur + c1_met + c1_his + c1_misc
    print(f"===== {t} (per unit biomass, mmol/gDW)")
    print(f"  glycine for protein {gly_prot:.4f}, purines {purines:.4f} (RNA {pur_rna:.4f}, DNA {pur_dna:.4f}, cofactors {pur_cof:.4f}), tetrapyrroles {heme_gly:.4f}, GSH {gsh_gly:.4f}")
    print(f"  Thr {thr:.4f} + Ile {ile:.4f} = {thr+ile:.4f} (made from glycine when there is no homoserine kinase)")
    print(f"  C1 demand: purines {c1_pur:.4f} + Met {c1_met:.4f} + His {c1_his:.4f} + misc {c1_misc:.4f} = {c1_dem:.4f}")
    print(f"  SHMT glycine (=C1 released) with Thr from Gly : {gly_withThr:.4f} -> C1 surplus {gly_withThr - c1_dem:+.4f}")
    print(f"  SHMT glycine (=C1 released) with Thr from Asp : {gly_noThr:.4f} -> C1 surplus {gly_noThr - c1_dem:+.4f}")
    print(f"  total P demand of biomass (for scale): ~0.6605 mmol P/gDW (r/(1+r) x 0.660485 is the ratio-bound methane)")
