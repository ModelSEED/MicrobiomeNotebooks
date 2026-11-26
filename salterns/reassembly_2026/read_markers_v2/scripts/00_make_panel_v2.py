#!/usr/bin/env python
"""v2 panel (meta/ko_panel_v2.tsv) = v1 panel rows (unchanged, incl. kegg_cap) + one counted label per KO of
community_ko/ko_panel_methylated.tsv not counted in v1 + new decoy KOs (close paralogs of the new targets; KO numbers
checked against KEGG REST list/ko of 2026-11-25).
Role flips v1 -> v2: K00499 cmo, K00479 gbcA, K21832 gbcB, K12905 hepD become counted targets.
Label = panel gene symbol, made filesystem-safe and unique per KO; KEGG symbols confirmed 2026-11-25 for the five DMSP
lyases (dddP/dddW/dddY/dddX/dddD), K23996 (phnY), K28979 (tdm), K24071 (gsmt-sdmt -> label gsdmt), K18897 (sdmt).
kegg_cap for new targets: 1500 bacterial genes, keeping every archaeal gene (the v1 marker rule), per-KO seed;
for new decoys: 1000 random prokaryote genes, per-KO seed. euk_cap: eukaryote genes fetched for counted targets
(<label>_euk decoy_reported class), random, per-KO seed.
Usage: python -I 00_make_panel_v2.py <v1 ko_panel.tsv> <ko_panel_methylated.tsv> <out.tsv>"""
import sys, csv
v1p, newp, out = sys.argv[1:4]
v1 = list(csv.DictReader(open(v1p), delimiter="\t"))
newrows = list(csv.DictReader(open(newp), delimiter="\t"))
FLIP = {"K00499": "CMO", "K00479": "gbcA", "K21832": "gbcB", "K12905": "hepD"}   # panel symbols (cmo was a v1 decoy label)
# Output contract (coordinator, 2026-11-25): label = panel `gene` verbatim, unique per KO; a repeated symbol -> <gene>_<KO>
# (pta K13788). Exception: the panel's explicit placeholders ('ddd (K28067)' ... 'gene symbol to confirm from KEGG', and
# 'mpnH/phnY*') get the KEGG symbols confirmed on 2026-11-25. label_safe = filesystem-safe form for file names.
LABEL = {"K28067": "dddP", "K28068": "dddW", "K28069": "dddY", "K28072": "dddX", "K28073": "dddD", "K23996": "phnY",
         "K13788": "pta_K13788"}
DECOYS = [  # ko, label, group (paralog family), what it guards
    ("K02567", "napA", "MoBisPGD_decoy", "torA/torZ/dmsA/ddhA/fwdB"), ("K00370", "narG", "MoBisPGD_decoy", ""),
    ("K22015", "fdhF", "MoBisPGD_decoy", ""), ("K00123", "fdoG", "MoBisPGD_decoy", ""), ("K08351", "bisC", "MoBisPGD_decoy", ""),
    ("K08352", "phsA", "MoBisPGD_decoy", ""), ("K07309", "ynfE", "MoBisPGD_decoy", ""), ("K08348", "fdnG", "MoBisPGD_decoy", ""),
    ("K22516", "fdhA_F420", "MoBisPGD_decoy", ""), ("K05299", "fdhA_NADP", "MoBisPGD_decoy", ""),
    ("K07307", "dmsB", "FeS_subunit_decoy", "ddhB"), ("K00371", "narH", "FeS_subunit_decoy", ""), ("K00124", "fdoH", "FeS_subunit_decoy", ""),
    ("K00374", "narI", "cytb_decoy", "ddhC"),
    ("K00001", "ADH", "MDR_decoy", "frmA/fdhA/mdo"), ("K13953", "adhP", "MDR_decoy", ""), ("K17068", "fdm", "MDR_decoy", "mdo"),
    ("K00086", "dhaT", "FeADH_decoy", "mdh/gbsB"), ("K08325", "yqhD", "FeADH_decoy", ""), ("K04072", "adhE", "FeADH_decoy", ""),
    ("K00048", "fucO", "FeADH_decoy", ""), ("K13921", "pduQ", "FeADH_decoy", ""),
    ("K17755", "codA", "GMC_decoy", "AOX/betA"),
    ("K00128", "ALDH", "ALDH_decoy", "betB"), ("K00135", "gabD", "ALDH_decoy", ""), ("K00140", "mmsA", "ALDH_decoy", ""),
    ("K00146", "feaB", "ALDH_decoy", ""), ("K00138", "aldB", "ALDH_decoy", ""), ("K14519", "aldH", "ALDH_decoy", ""),
    ("K00285", "dadA", "DAAO_decoy", "soxB/dmgo/soxA_mono"), ("K03153", "thiO", "DAAO_decoy", ""), ("K00111", "glpA", "DAAO_decoy", ""),
    ("K15736", "lhgO", "DAAO_decoy", ""), ("K09471", "puuB", "DAAO_decoy", ""), ("K00306", "PIPOX", "DAAO_decoy", "soxA_mono"),
    ("K19191", "mabO", "DAAO_decoy", ""),
    ("K00605", "gcvT", "GcvT_decoy", "dmdA/tdm/dmgdh/sardh/soxA"), ("K06980", "ygfZ", "GcvT_decoy", ""), ("K20155", "mao", "GcvT_decoy", ""),
    ("K01915", "glnA", "GS_decoy", "gmaS"),
    ("K00265", "gltB", "GltB_decoy", "mgsABC/fwdC"), ("K00266", "gltD", "GltB_decoy", ""), ("K00284", "gltS", "GltB_decoy", ""),
    ("K00548", "metH", "B12_MT_decoy", "corrinoid proteins/bhmt"), ("K16955", "mtsB", "B12_MT_decoy", "corrinoid proteins"),
    ("K25221", "mtgA", "UROD_MT_decoy", "mtaA/mtbA/mtsA"), ("K26930", "mtqA", "UROD_MT_decoy", ""),
    ("K25915", "mtvA", "UROD_MT_decoy", ""), ("K25917", "mtoA", "UROD_MT_decoy", ""),
    ("K15023", "acsE", "MeTr_decoy", "mtrH"),
    ("K00485", "FMO", "FMO_decoy", "tmm"), ("K03379", "chnB", "FMO_decoy", ""),
    ("K00219", "fadH", "OYE_decoy", "tmd"), ("K10680", "nemA", "OYE_decoy", ""), ("K00354", "OYE", "OYE_decoy", ""),
    ("K01897", "fadD", "AMPbinding_decoy", "acs/dmdB"), ("K01908", "prpE", "AMPbinding_decoy", ""), ("K02182", "caiC", "AMPbinding_decoy", ""),
    ("K01907", "AACS", "AMPbinding_decoy", ""),
    ("K00249", "acd", "ACAD_decoy", "dmdC"), ("K06445", "fadE", "ACAD_decoy", ""), ("K00248", "bcd", "ACAD_decoy", ""),
    ("K00253", "ivd", "ACAD_decoy", ""),
    ("K01692", "paaF", "crotonase_decoy", "dmdD"), ("K01715", "crt", "crotonase_decoy", ""), ("K13766", "liuC", "crotonase_decoy", ""),
    ("K15866", "paaG", "crotonase_decoy", ""),
    ("K00929", "buk", "ASKHA_decoy", "ackA"), ("K00932", "tdcD", "ASKHA_decoy", ""),
    ("K00634", "ptb", "PTA_PTB_decoy", "pta/eutD"),
    ("K00198", "cooS", "CODH_decoy", "cdhA"), ("K14138", "acsB", "ACS_decoy", "cdhC"),
    ("K00058", "serA", "2HADH_decoy", "fdh"), ("K00090", "ghrB", "2HADH_decoy", ""), ("K12972", "ghrA", "2HADH_decoy", ""),
    ("K03778", "ldhA", "2HADH_decoy", ""), ("K03473", "pdxB", "2HADH_decoy", ""),
    ("K11175", "purN", "formyltransferase_decoy", "purU"), ("K00604", "fmt", "formyltransferase_decoy", ""),
    ("K00547", "mmuM", "HMT_decoy", "bhmt"),
    ("K00114", "exaA", "PQQ_ADH_decoy", "mxaF/xoxF"), ("K17760", "qheDH", "PQQ_ADH_decoy", ""),
    ("K01262", "pepP", "M24_decoy", "dddP"), ("K01271", "pepQ", "M24_decoy", ""), ("K01265", "map", "M24_decoy", ""),
    ("K01905", "acdA", "AcdA_decoy", "dddX"), ("K22224", "acdB", "AcdA_decoy", ""), ("K24012", "acdAB", "AcdA_decoy", ""),
    ("K18593", "hbcl_ADP", "AcdA_decoy", ""), ("K18594", "hpcs_ADP", "AcdA_decoy", ""),
    ("K08298", "caiB", "CoAT3_decoy", "dddD"), ("K07749", "frc", "CoAT3_decoy", ""), ("K15871", "baiF", "CoAT3_decoy", ""),
    ("K18702", "uctC", "CoAT3_decoy", ""), ("K13607", "fldA", "CoAT3_decoy", ""), ("K20882", "hadA", "CoAT3_decoy", ""),
    ("K18896", "gsmt", "GB_MT_decoy", "gsdmt/sdmt"), ("K13042", "bsmB", "GB_MT_decoy", ""),
    ("K21195", "phnY_AEP", "Pn_dioxygenase_decoy", "phnY/phnZ"), ("K21196", "phnZ_AHEP", "Pn_dioxygenase_decoy", ""),
    ("K24290", "tmpB", "Pn_dioxygenase_decoy", ""), ("K24291", "tmpA", "Pn_dioxygenase_decoy", ""),
    ("K13812", "fae_hps", "fae_decoy", "fae"),
    ("K00494", "luxA", "LLM_decoy", "dmoA/ssuD/sfnG/mer"), ("K20938", "ladA", "LLM_decoy", ""), ("K15510", "fgd1", "LLM_decoy", ""),
    ("K05337", "fer", "ferredoxin_decoy", "fwdH/fwdF/fwdG"), ("K03616", "rnfB", "ferredoxin_decoy", ""),
    ("K00176", "korD", "ferredoxin_decoy", ""), ("K14127", "mvhD", "ferredoxin_decoy", ""), ("K03855", "fixX", "ferredoxin_decoy", ""),
    ("K00196", "cooF", "ferredoxin_decoy", ""),
    ("K02040", "pstS", "SBP_decoy", "phnD"), ("K11081", "phnS", "SBP_decoy", ""),
]
rows = []
v1kos = set()
for r in v1:
    r = dict(r)
    v1kos.add(r["ko"])
    if r["ko"] in FLIP:
        r["role"] = "target"; r["gene"] = FLIP[r["ko"]]
        r["note"] = (r["note"] + "; " if r["note"] else "") + "v2: counted target (was a v1 decoy)"
    r["source"] = "v1"
    r["euk_cap"] = "300" if r["role"] == "target" else "0"
    rows.append(r)
for r in newrows:
    if r["ko"] in v1kos:
        continue
    lab = LABEL.get(r["ko"], r["gene"])
    rows.append(dict(ko=r["ko"], gene=lab, role="target", group=f'{r["compound"]}|{r["process"]}', kegg_cap="1500",
                     note=f'{r["status"]}; {r["specificity"]}; {r["note"]}'.strip("; "), source="v2_panel", euk_cap="300"))
for ko, lab, grp, guards in DECOYS:
    assert ko not in v1kos and ko not in {r["ko"] for r in newrows}, ko
    rows.append(dict(ko=ko, gene=lab, role="decoy", group=grp, kegg_cap="1000", note=f"v2 decoy vs {guards}" if guards else "v2 decoy",
                     source="v2_decoy", euk_cap="0"))
labs = {}
for r in rows:
    if r["gene"] in labs and r["gene"] not in ("uL4", "uL16"):
        raise SystemExit(f"duplicate label {r['gene']} {r['ko']} {labs[r['gene']]}")
    labs[r["gene"]] = r["ko"]
kos = [r["ko"] for r in rows]
assert len(kos) == len(set(kos))
import re
for r in rows:
    r["label_safe"] = re.sub(r"[^A-Za-z0-9_.-]+", "_", r["gene"]).strip("_")
safe = [r["label_safe"] for r in rows if r["gene"] not in ("uL4", "uL16")]
assert len(safe) == len(set(safe)), "label_safe collision"
cols = ["ko", "gene", "label_safe", "role", "group", "kegg_cap", "euk_cap", "source", "note"]
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
newk = {r["ko"] for r in newrows}
print("rows", len(rows), "targets", sum(r["role"] == "target" for r in rows), "decoys", sum(r["role"] == "decoy" for r in rows),
      "panel KOs covered as counted labels:", sum(1 for r in rows if r["ko"] in newk and r["role"] in ("target",)), "of", len(newk))
missing = [k for k in newk if not any(r["ko"] == k and r["role"] == "target" for r in rows)]
print("panel KOs not counted:", missing)
