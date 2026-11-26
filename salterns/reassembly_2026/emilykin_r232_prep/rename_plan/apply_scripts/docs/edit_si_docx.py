"""Docs stage, step C2: tracked changes (author 'Claude (GTDB R232 update)') in
~/Documents/EmilyKin/SI-paper3 (3)_correction.docx, per docs/PLAN.md 4.3 + CRITIC + decisions D3/D7."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from docx_track import Doc, W, fix_settings_track_revisions

U = sys.argv[1]; DATE = sys.argv[2]
R = os.path.expanduser("~/Documents/EmilyKin")
P = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
d = Doc(U, DATE)

# ---- para 209 (ppk1 text): genome-level GTDB R232 species (D7)
p = d.find_para("All six Accumulibacter MAGs carried a recoverable, full-length ppk1 sequence")
d.insert_after(p, "and coasm_bin.476 to clade ID (100/100).",
               " At genome level, GTDB R232 assigns coasm_bin.185, coasm_bin.250 and coasm_bin.347 to *Accumulibacter "
               "regalis*, *A. delftensis* and *A. phosphatis*, matching these clades, and coasm_bin.476 to "
               "*Accumulibacter* sp032229525 (98.68% ANI).", label="para209 GTDB R232 species (D7)",
               comment="Added for GTDB R232 (scientific statement). The R232 genome-level species now agree with the "
                       "ppk1 clades: coasm_bin.185 A. regalis (97.58% ANI), coasm_bin.250 A. delftensis (98.60%), "
                       "coasm_bin.347 A. phosphatis (98.14%). coasm_bin.476 is Accumulibacter sp032229525 (98.68% ANI, "
                       "AF 0.84; a representative added to GTDB in R226), so it is not a novel species, although its "
                       "clade ID ppk1 remains nearest only to sparsely characterized references. Not changed here: "
                       "CAN_1_bin.98 is called a 'candidate-novel Type I lineage' below; that refers to its ppk1 only - "
                       "its genome belongs to GTDB species Accumulibacter sp014879755 (99.48% ANI) in both r220 and "
                       "R232. Consider 'divergent, unclassified Type I ppk1 lineage'.")

# ---- para 212: CAN_5_bin.40 genus; Figure S4 reference labels
p = d.find_para("A quality audit of the 188 overwrites")
d.replace(p, "Giesbergeria", "Simplicispira", label="para212 CAN_5_bin.40 genus")
d.insert_after(p, "among their closest GTDB r220 reference genomes",
               ", labelled with their GTDB R232 names,", label="para212 Fig S4 reference names")

# ---- para 216: naming cases
p = d.find_para("On the best-of-both 291 catalogue this yields 155 ASV-inherited")
d.insert_after(p, "On the best-of-both 291 catalogue", " classified with GTDB R232", label="para216 R232")
d.replace(p, "155 ASV-inherited, 62 GTDB-genus, and 74 GTDB-lowest",
          "138 ASV-inherited, 70 GTDB-genus, and 83 GTDB-lowest", label="para216 138/70/83",
          comment="GTDB R232 update: the same naming rule, re-run on the R232 taxonomy with the ASV iterativeIDs held "
                  "fixed, gives 138/70/83 (r220: 155/62/74). 51 MAG iterativeIDs changed (50 at the root); no ASV "
                  "iterativeID changed. Old -> new pairs by bin: mag_iterativeID_r220_to_r232.csv.")

# ---- para 252: Figure S4 legend
p = d.find_para("Figure S4. Maximum-likelihood genome tree of 65 selected MAGs")
d.insert_after(p, "61 GTDB r220 species-representative reference genomes",
               " (labelled with their GTDB R232 names; GCA_017302355.1 is not in R232)", label="para252 R232 names")
d.replace(p, "rooted on the Patescibacteria clade", "rooted on the Patescibacteriota clade", label="para252 root")
d.replace(p, "for CAN_4_bin.208, whose largest-cluster Azospira_A representative",
          "for CAN_4_bin.208 (Proximibacter in GTDB R232), whose largest-cluster r220 Azospira_A representative",
          label="para252 CAN_4_bin.208")
d.replace(p, "the Azospira_A representative sharing", "the r220 Azospira_A representative sharing",
          label="para252 Azospira_A (2)")
d.replace(p, "(Azospira_A restricta, GCF_016858125.1)",
          "(Azospira_A restricta, GCF_016858125.1; Proximibacter restrictus in R232)", label="para252 reference")
d.replace(p, "best skani ANI to a GTDB r220 representative (Table S6)",
          "best skani ANI to a GTDB R232 representative (Table S6)", label="para252 skani release")

# ---- Figure S4 / S5 images
d.replace_image("image5.png", f"{R}/fig1_tree/fig1_tree.png", label="Figure S4 image",
                comment="Figure S4 re-rendered with GTDB R232 data (scripts/render_fig1_tree.py); same RAxML tree and "
                        "rooting (116 of 123 nodes >= 70% support). 8 MAG tips relabelled: Ca_Competibacter.32 -> "
                        "Contendibacter.1_m, Azospira.3 -> Proximibacter.2_m, Flavipsychrobacter.6_m -> Taibaiella.5, "
                        "Ignavibacteria.7 -> Tepidaquicellaceae.3_m, Ignavibacteria.5 -> Tepidaquicellaceae.1_m, "
                        "Ignavibacteria.21 -> Tepidaquicellales.1_m, Ca_Falkowbacteria.7 -> Patescibacterium.1_m, "
                        "Paceibacteria.1_m -> Minisyncoccia.1_m. Reference tips show R232 names (e.g. Accumulibacter "
                        "regalis, A. delftensis, Proximibacter restrictus, Contendibacter sp018262135, Lysobacter_I "
                        "sp019510125, Patescibacterium sp002344425, JAFIGV01 sp937871155, Niabella defluvii); "
                        "GCA_017302355.1 left GTDB and is marked 'not in GTDB R232'. Red circles follow the rebuilt "
                        "Table S6: 10 MAGs change ANI class (e.g. CAN_5_bin.70 >=95 -> 85-95 because its r220 reference "
                        "left GTDB; coasm_bin.53, coasm_bin.10, coasm_bin.579, coasm_bin.489, coasm_bin.321 and "
                        "coasm_bin.521 now >=95). Phylum colours: Patescibacteriota.")
d.replace_image("image4.png", f"{R}/clade_classify/out/nosz_top5_tree.png", label="Figure S5 image",
                comment="Figure S5 re-rendered (clade_classify/render_nosz_top5_tree.py): the only change is the "
                        "rank-5 MAG CAN_3_bin.221, now labelled Tepidaquicellales.1_m (GTDB R232 iterativeID; was "
                        "Ignavibacteria.21; genus OLB5 unchanged). Tree, ranking and abundances are unchanged.")

d.save()
nbad = fix_settings_track_revisions(U)
print("settings: removed invalid w:trackChanges x", nbad, "; w:trackRevisions present")
for x in d.log:
    print(json.dumps(x, ensure_ascii=False)[:300])
json.dump(d.log, open(f"{P}/apply_scratch/docs/si_docx_changes.json", "w"), indent=1, ensure_ascii=False)
