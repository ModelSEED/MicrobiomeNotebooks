#!/usr/bin/env python
"""Build marker_ko_table.tsv and the two custom DRAM distillate sheets from the cached sources in sources/.

marker_ko_table.tsv: one row per marker KO this paper's guilds and pathway figures rely on (markers_input.tsv), with
  - the KEGG name/EC (KEGG REST `list/ko`, 2026-11-11) and the KOfam score threshold + score type in two releases
    (2026-10-26 current; 2026-05-28, the profiles used for the local control test),
  - whether DRAM 1.5.0 / DRAM2 v2.0.0-beta42 would *report* it in their distillates (the KOfam call itself only needs
    a threshold; DRAM skips KOs whose threshold is '-'),
  - what the KBase/ModelSEED translation tables turn the KO or its EC into (cb_annotation_ontology_api data),
  - which guilds_r232.csv evidence column used the corresponding RAST string.
custom_distillate_salterns.dram1.tsv / .dram2.tsv: the same markers as a DRAM custom distillate sheet
  (DRAM1: `DRAM.py distill --custom_distillate`; DRAM2: `--sum_custom`).
Run: ~/Documents/py_venv/bin/python build_marker_table.py
"""
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = HERE / "sources"
ONTO = Path.home() / "Documents/cb_annotation_ontology_api/data"

m = pd.read_csv(HERE / "markers_input.tsv", sep="\t", dtype=str)
kegg = pd.read_csv(SRC / "kegg/kegg_rest_list_ko_2026-11-11.tsv", sep="\t", header=None, names=["ko", "kegg_name"],
                   index_col=0)
kof_new = pd.read_csv(SRC / "kegg/kofam_ko_list_2026-10-26.gz", sep="\t", index_col=0, dtype=str)
kof_old = pd.read_csv(SRC / "kegg/kofam_ko_list_2026-05-28_profiles_used_locally.gz", sep="\t", index_col=0, dtype=str)

# DRAM1 1.5.0 distillate sheets
d1 = pd.read_csv(SRC / "dram1_v1.5.0/genome_summary_form.tsv", sep="\t", dtype=str)
d1h = pd.read_csv(SRC / "dram1_v1.5.0/function_heatmap_form.tsv", sep="\t", dtype=str)
d1m = pd.read_csv(SRC / "dram1_v1.5.0/module_step_form.tsv", sep="\t", dtype=str)

# DRAM2 beta42 sheets; which are loaded without --sum_custom
D2_TOPIC = {"assimilation_and_cofactor_metabolism": "topic (default)", "cellular_machinery": "topic (default)",
            "energy_acquisition_bioenergetics": "topic (default)",
            "environmental_interaction_and_adaptation": "topic (default)",
            "distill_ag": "eco (--sum_ecos ag)", "distill_aquatic": "eco (--sum_ecos aquatic)",
            "distill_bgc": "eco (--sum_ecos bgc)", "distill_eng_sys": "eco (--sum_ecos eng_sys)",
            "distill_gut": "eco (--sum_ecos gut)", "methylotrophy_distillate": "shipped, load with --sum_custom"}
d2 = {}
for f in sorted((SRC / "dram2_v2.0.0-beta42").glob("*.tsv")):
    txt = f.read_text()
    d2[f.stem] = set(re.findall(r"\bK\d{5}\b", txt))

# KBase / ModelSEED translation (KEGG release 95.0 KO->rxn; EC->rxn)
ko_seed = pd.read_csv(ONTO / "kegg_95_0_ko_seed.tsv", sep="\t", dtype=str).set_index("ko_id")
ec_tr = pd.read_csv(ONTO / "EC_translation.tsv", sep="\t", dtype=str)
ec_rxn = ec_tr.groupby("External ID")["ModelSEED ID"].apply(lambda s: ";".join(sorted(set(s))))

GUILD_COL = {"K06163": "ev_phnJ (guild: mepn_consumer)", "K06164": "ev_phn_accessory", "K06165": "ev_phn_accessory",
             "K06166": "ev_phn_accessory", "K05780": "ev_phn_accessory", "K06162": "ev_phn_accessory",
             "K05774": "ev_phn_accessory", "K21579": "ev_grdH (guild: betaine_reducer)",
             "K21578": "ev_grdH (guild: betaine_reducer)", "K22443": "ev_cntA (guild: tma_producer)",
             "K20038": "ev_cutC_cluster matches only BMC genes 'clustered with' CutC (guild: tma_producer)",
             "K00399": "ev_mcrA (guild: methanogen)", "K16178": "ev_mtbB_mtmB", "K16176": "ev_mtbB_mtmB",
             "K16179": "ev_mtbB_mtmB (RAST 'corrinoid protein' strings also match)",
             "K16177": "ev_mtbB_mtmB (RAST 'corrinoid protein' strings also match)"}

rows = []
for _, r in m.iterrows():
    k = r.ko
    name = kegg.kegg_name.get(k, "")
    ecs = re.findall(r"\[EC:([^\]]+)\]", name)
    ecs = ecs[0].split() if ecs else []
    d1_rows = d1[d1.gene_id == k]
    d1_desc = "; ".join(sorted({f"{a} > {b} > {c}".replace(" > nan", "") for a, b, c in
                                zip(d1_rows.sheet, d1_rows.header, d1_rows.subheader)}))
    d1_hm = "; ".join(sorted(d1h[d1h.function_ids.fillna("").str.contains(k)].function_name.str.strip()))
    d1_mod = ";".join(sorted(set(d1m[d1m.ko == k].module)))
    d2_hits = [s for s, kos in d2.items() if k in kos]
    d2_default = [s for s in d2_hits if D2_TOPIC.get(s, "").startswith("topic")]
    rows.append({
        "ko": k, "gene": r.gene, "group": r.group, "role": r.role, "kegg_name": name, "ec": " ".join(ecs),
        "kofam_threshold_2026-10": kof_new.threshold.get(k, "absent"),
        "kofam_score_type_2026-10": kof_new.score_type.get(k, "absent"),
        "kofam_threshold_2026-05": kof_old.threshold.get(k, "absent"),
        "kofam_score_type_2026-05": kof_old.score_type.get(k, "absent"),
        "dram_kofam_callable": "no: KOfam threshold '-' (DRAM skips)" if kof_new.threshold.get(k) == "-" else "yes",
        "dram1_metabolism_summary_row": d1_desc or "none",
        "dram1_product_heatmap_row": d1_hm or "none",
        "dram1_module_step_modules": d1_mod or "none",
        "dram2_sheets": ", ".join(f"{s} [{D2_TOPIC.get(s, '?')}]" for s in d2_hits) or "none",
        "dram2_in_default_summary": "yes" if d2_default else "no",
        "ko_to_modelseed_rxn_kegg95": (ko_seed.seed_ids.get(k) if k in ko_seed.index and
                                       isinstance(ko_seed.seed_ids.get(k), str) else "none"),
        "ec_to_modelseed_rxn": "; ".join(f"{e}: {ec_rxn.get(e, 'none')}" for e in ecs) or "no EC",
        "guilds_r232_evidence_column": GUILD_COL.get(k, ""),
    })
t = pd.DataFrame(rows)
t.to_csv(HERE / "marker_ko_table.tsv", sep="\t", index=False)

# custom distillate sheets (KOfam-callable markers only)
c = t[t.dram_kofam_callable == "yes"]
d1c = pd.DataFrame({"gene_id": c.ko, "gene_description": c.gene + "; " + c.role, "module": c.group,
                    "sheet": "Salterns methane", "header": c.group, "subheader": c.gene, "potential_amg": False})
d1c.to_csv(HERE / "custom_distillate_salterns.dram1.tsv", sep="\t", index=False)
d2c = pd.DataFrame({"gene_id": c.ko, "gene_description": c.gene + "; " + c.role, "module": c.group,
                    "topic_ecosystem": "salterns_methane", "header": c.group, "subheader": c.gene})
d2c.to_csv(HERE / "custom_distillate_salterns.dram2.tsv", sep="\t", index=False)

print(f"{len(t)} markers; KOfam threshold '-': {', '.join(t.loc[t.dram_kofam_callable != 'yes', 'gene'])}")
print("in DRAM1 metabolism summary:", ", ".join(t.loc[t.dram1_metabolism_summary_row != "none", "gene"]))
print("in DRAM1 product heatmap:", ", ".join(t.loc[t.dram1_product_heatmap_row != "none", "gene"]))
print("in DRAM2 default topic sheets:", ", ".join(t.loc[t.dram2_in_default_summary == "yes", "gene"]))
print("in DRAM2 methylotrophy sheet only via --sum_custom:",
      ", ".join(t.loc[t.dram2_sheets.str.contains("methylotrophy") & (t.dram2_in_default_summary == "no"), "gene"]))
print("in no DRAM2 sheet:", ", ".join(t.loc[t.dram2_sheets == "none", "gene"]))
print("no KO->ModelSEED rxn:", ", ".join(t.loc[t.ko_to_modelseed_rxn_kegg95 == "none", "gene"]))
print("threshold changed 2026-05 -> 2026-10 by >15%:", ", ".join(
    f"{g} ({a}->{b})" for g, a, b in zip(t.gene, t["kofam_threshold_2026-05"], t["kofam_threshold_2026-10"])
    if a not in ("-", "absent") and b not in ("-", "absent") and abs(float(b) - float(a)) / float(a) > 0.15))
