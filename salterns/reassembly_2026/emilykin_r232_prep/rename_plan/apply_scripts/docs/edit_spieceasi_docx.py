"""Docs stage, step D: tracked changes in ~/Documents/EmilyKin/SPIEC-EASI_network_analysis.docx (plan paras 140/357 =
python-docx paragraphs 53/72: 'PAO-like MAG (Azospira.3 ...)' -> ASV node wording; plus paragraphs 14/16, which state
the MAG-link rule that decision D2 replaced by the frozen R220-era list).  Node labels stay ASV IDs (D2)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from docx_track import Doc

U = sys.argv[1]; DATE = sys.argv[2]
P = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
d = Doc(U, DATE)

D2 = ("GTDB R232 update (2026-11-12), decision D2: renaming MAGs must not change the network. The MAG link is now "
      "the frozen list of the 140 ASV-MAG pairs whose IDs matched under the r220-based naming "
      "(network/mag_linked_asvs_r220_frozen.tsv), so the 53 MAG-linked nodes, edges, modules and tests are unchanged. "
      "Node labels remain ASV iterativeIDs; the MAG's current (R232) iterativeID is reported by bin in the output "
      "tables (e.g. spieceasi_mag_nodes.csv, column MAG_iterativeID). For 17 network nodes the ASV label and the MAG "
      "iterativeID now differ (e.g. ASV Azospira.3 -> MAG CAN_4_bin.208 = Proximibacter.2_m; ASV Ca_Competibacter.32 "
      "-> MAG coasm_bin.365 = Contendibacter.1_m).")

p = d.find_para("The shared iterativeID nomenclature makes the linkage explicit")
d.replace(p, "an ASV whose iterativeID equals a MAG iterativeID is a genome-resolved node.",
          "an ASV whose iterativeID equalled a MAG iterativeID under the GTDB r220-based naming is a genome-resolved "
          "node; these 140 ASV–MAG pairs are frozen, so the GTDB R232 renaming of 51 MAG iterativeIDs does not change "
          "the network.", label="para14 link rule", comment=D2)
d.replace(p, "(`mag_iterativeID_old_to_new.json`;", "(`network/mag_linked_asvs_r220_frozen.tsv`;",
          label="para14 link source")

p = d.find_para("Each ASV node is joined to its guild through its MAG")
d.replace(p, "the current naming map (`mag_iterativeID_old_to_new.json`), not through the guild table's own iterativeID "
             "column, which carries an older naming map.",
          "the frozen ASV–MAG list (`network/mag_linked_asvs_r220_frozen.tsv`), never through an iterativeID string.",
          label="para16 guild join")

p = d.find_para("all of them join a non-canonical gene-content PAO-like MAG (Azospira.3")
d.replace(p, "gene-content PAO-like MAG (Azospira.3,",
          "gene-content PAO-like MAG-linked ASV node (Azospira.3 [MAG CAN_4_bin.208, Proximibacter.2_m under GTDB R232],",
          label="para53 (plan 140) Azospira.3 node",
          comment="Wording only: these are ASV nodes (labels unchanged, D2). Under GTDB R232 the MAG behind node "
                  "Azospira.3 (CAN_4_bin.208) is Proximibacter.2_m and the MAG behind node Ca_Competibacter.32 "
                  "(coasm_bin.365) is Contendibacter.1_m; Pseudoxanthomonas.23, Zoogloea.2, Ca_Competibacter.26 and .51 "
                  "keep their IDs. coasm_bin.365 stays canonical GAO (family Competibacteraceae).")
d.replace(p, "(.26, .32, .51)", "(.26, .32 [MAG coasm_bin.365, Contendibacter.1_m under GTDB R232], .51)",
          label="para53 (plan 140) .32 node")

p = d.find_para("between gene-content PAO-like MAGs (Azospira.3")
d.replace(p, "gene-content PAO-like MAGs (Azospira.3,",
          "gene-content PAO-like MAG-linked ASV nodes (Azospira.3 [MAG CAN_4_bin.208, Proximibacter.2_m under GTDB "
          "R232],", label="para72 (plan 357) Azospira.3 node")

d.save()
for x in d.log:
    print(json.dumps(x, ensure_ascii=False)[:300])
json.dump(d.log, open(f"{P}/apply_scratch/docs/spieceasi_docx_changes.json", "w"), indent=1, ensure_ascii=False)
