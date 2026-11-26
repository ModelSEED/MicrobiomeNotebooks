"""Docs stage, step C1: tracked changes (author 'Claude (GTDB R232 update)') in
~/Documents/EmilyKin/papert3_070626 (2)_correction.docx, per docs/PLAN.md 4.2 + CRITIC + decisions D3/D7.
Run: edit_main_docx.py <unpacked dir>  (the caller unzips the backup and rezips)."""
import csv, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from docx_track import Doc, W, fix_settings_track_revisions

U = sys.argv[1]
DATE = sys.argv[2]
R = os.path.expanduser("~/Documents/EmilyKin")
P = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
d = Doc(U, DATE)

# ---- para 45: GTDB-Tk version / release
p = d.find_para("Taxonomic assignments were generated using GTDB-Tk v2.7.2 against Genome Taxonomy Database release r220")
d.replace(p, "release r220", "release R232", label="para45 release")
d.replace(p, "v2.7.2", "v2.8.0 (classify_wf)", label="para45 version",
          comment="GTDB R232 update (2026-11-12): all 291 MAGs were re-classified locally with GTDB-Tk v2.8.0 "
                  "classify_wf against GTDB release R232, the reference data this GTDB-Tk version expects (the earlier "
                  "v2.7.2 run used r220 data and logged a version-mismatch warning). Consider adding a citation for the "
                  "R232 release; the plan suggests Parks et al. 2025, Nucleic Acids Res., doi:10.1093/nar/gkaf1040 "
                  "(please verify before citing). Per-MAG changes: mag_iterativeID_r220_to_r232.csv and Tables S1/S6.")

# ---- para 57: phyla
p = d.find_para("These 291 MAGs spanned 14 bacterial phyla")
d.replace(p, "14 bacterial phyla", "15 bacterial phyla", label="para57 phylum count",
          comment="GTDB R232 renames Patescibacteria to Patescibacteriota (the same 21 MAGs) and places one former "
                  "Myxococcota MAG (coasm_bin.258) in Myxococcota_A, so the catalogue has 15 GTDB phylum names. If you "
                  "prefer to merge GTDB suffixed names (Figure 1 does this for genera), keep 14 phyla and write "
                  "Myxococcota (3).")
d.replace(p, "Patescibacteria (21)", "Patescibacteriota (21)", label="para57 Patescibacteriota")
d.replace(p, "Myxococcota (3)", "Myxococcota (2), Myxococcota_A (1)", label="para57 Myxococcota")

# ---- para 59: Table 1 legend
p = d.find_para("GTDB genus, GTDB-Tk (r220) genus, or family")
d.replace(p, "GTDB-Tk (r220) genus", "GTDB-Tk (R232) genus", label="para59 Table 1 legend",
          comment="Table 1 updated by MAG bin to GTDB R232. iterativeID changes: CAN_3_bin.221 Ignavibacteria.21 -> "
                  "Tepidaquicellales.1_m; CAN_4_bin.208 Azospira.3 -> Proximibacter.2_m; coasm_bin.365 "
                  "Ca_Competibacter.32 -> Contendibacter.1_m; CAN_1_bin.2 Ignavibacteria.5 -> Tepidaquicellaceae.1_m; "
                  "CAN_1_bin.195 Ca_Falkowbacteria.7 -> Patescibacterium.1_m; CAN_1_bin.48 Paceibacteria.1_m -> "
                  "Minisyncoccia.1_m; coasm_bin.40 Flavipsychrobacter.6_m -> Taibaiella.5; coasm_bin.337 "
                  "Ignavibacteria.7 -> Tepidaquicellaceae.3_m. GTDB genus changes: SCMT01 -> Lysobacter_I, Azospira_A "
                  "-> Proximibacter, JARSSO01 -> Contendibacter, UBA919 -> Patescibacterium, Flavipsychrobacter -> "
                  "JAFIGV01, and new genus names where only a family was given before (f__B-1AR -> JAFDZO01, f__UBA2400 "
                  "-> JBDHTQ01, f__GWB1-50-10 -> CAYHLH01). The old IDs Ignavibacteria.21/.5/.7, Azospira.3, "
                  "Ca_Competibacter.32 and Ca_Falkowbacteria.7 are still ASV names elsewhere; only the MAG labels "
                  "changed. coasm_bin.365 keeps its guild (family Competibacteraceae).")

# ---- Table 1 cells, keyed by the bin in column 1
idmap = json.load(open(f"{R}/mag_iterativeID_old_to_new.json"))
gt = {r["user_genome"]: r["classification"].split(";") for r in
      csv.DictReader(open(f"{R}/gtdbtk.bac120.summary.tsv", newline=""), delimiter="\t")}
plan = [r for r in csv.DictReader(open(f"{P}/docs/table1_plan.tsv"), delimiter="\t")
        if r["file"] == "papert3_070626 (2)_correction.docx"]
tbl = [t for t in d.body.iter(W + "tbl")][0]
rows = tbl.findall(W + "tr")
hdr = [d.text(tc.find(W + "p")) for tc in rows[0].findall(W + "tc")]
assert hdr[:3] == ["MAG", "iterativeID", "GTDB genus"], hdr
byb = {}
for tr in rows[1:]:
    cells = tr.findall(W + "tc")
    b = "".join(d.text(x) for x in cells[0].iter(W + "p"))
    byb[b] = cells
assert len(byb) == 53
for r in plan:
    cells = byb[r["bin"]]
    col = 1 if "iterativeID" in r["location"] else 2
    ps = cells[col].findall(W + "p")
    assert len(ps) == 1
    assert d.text(ps[0]) == r["current"], (r["bin"], d.text(ps[0]), r["current"])
    ital = None
    if col == 2:  # GTDB genus: Latin names italic, alphanumeric placeholders roman (the table's convention)
        ital = not any(ch.isdigit() for ch in r["proposed"])
    d.replace(ps[0], r["current"], r["proposed"], label=f"Table 1 {r['bin']} col{col+1}", italic=ital)
# D1 gate on Table 1: every row's ID and genus (or f__family) equal R232
bad = []
for b, cells in byb.items():
    idv = d.text(cells[1].find(W + "p"))
    g = d.text(cells[2].find(W + "p"))
    lin = gt[b]
    want_g = lin[5][3:] if lin[5][3:] else "f__" + lin[4][3:]
    if idv != idmap[b] or g != want_g:
        bad.append((b, idv, idmap[b], g, want_g))
print("Table 1 D1 gate: rows", len(byb), "mismatches", bad)
assert not bad

# ---- Figure 1 image (main image1) -> regenerated per-cell composition
d.replace_image("image1.png", f"{R}/figures/percell_composition.png", label="Figure 1 image",
                comment="Figure 1 re-rendered with GTDB R232 genera (scripts/render_percell_composition.py; same "
                        "script, data and colours). Legend changes: SCMT01 -> Lysobacter (GTDB Lysobacter_I), Azospira "
                        "-> Proximibacter, Giesbergeria -> Simplicispira, UBA5532 (Patescibacteria) -> UBA5532 "
                        "(Patescibacteriota); each keeps the same MAGs and colour. JJ008 is now CAN_5_bin.112 only "
                        "(coasm_bin.492 is CADEDV01 in R232 and is pooled in Others). Unclassified falls from 1.89% to "
                        "0.57% (five-phase mean; 11 MAGs gained a GTDB genus) and becomes the last band. Competibacter "
                        "and Accumulibacter values are unchanged.")

# ---- para 656: D7 claim revision (Table S6 rebuild + coasm_bin.476)
p = d.find_para("Of the 291 MAGs, 48 were substantially divergent from currently available genomes")
d.replace(p, "48", "28", label="para656 48->28")
d.replace(p, "77", "66", label="para656 77->66")
d.replace(p, "three of these were closely related to existing references (≥97% ANI), whereas the near-complete "
             "coasm_bin.476 (99.98% completeness, 0.20% contamination) aligned poorly to any reference and likely "
             "represents a novel Accumulibacter lineage.",
          "all six belong to GTDB R232 species, and five are closely related to their species representative (≥97% "
          "ANI), including the near-complete coasm_bin.476 (99.98% completeness, 0.20% contamination; *Accumulibacter* "
          "sp032229525, 98.68% ANI).", label="para656 coasm_bin.476 novelty (D7)",
          comment="Scientific change from GTDB R232, not a rename. Table S6 was rebuilt with skani v0.4.0 against the "
                  "GTDB R232 species representatives with the same settings as before (--min-af 0 -n 3): <85% ANI 48 -> "
                  "28 (23 of them with no hit), unresolved (<15% aligned) 77 -> 66, >=95% 96 -> 121. All six "
                  "Accumulibacter MAGs now have a GTDB species. coasm_bin.476 is Accumulibacter sp032229525 (GTDB-Tk "
                  "98.68% ANI, AF 0.84 to GCF_032229525.1, a representative added to GTDB in R226; skani AF_query 72.9%), "
                  "so it is no longer a candidate novel species; its clade ID ppk1 is still nearest only to sparsely "
                  "characterized references (BPBW2884/BPBW2899). >=97% ANI: coasm_bin.185 A. regalis 97.58, coasm_bin.250 "
                  "A. delftensis 98.60, coasm_bin.347 A. phosphatis 98.14, CAN_1_bin.98 sp014879755 99.48, coasm_bin.476 "
                  "98.68; CAN_4_bin.64 (sp003538495) is 96.05. The R232 species of the three named-clade MAGs now match "
                  "their ppk1 clades (IA regalis, IC delftensis, IIA phosphatis).")

# ---- Figure 4 caption: label map comment (no generator; not re-rendered)
fm = list(csv.DictReader(open(f"{P}/apply_scratch/docs/fig4_label_map.tsv"), delimiter="\t"))
m232 = json.load(open(f"{R}/mag_iterativeID_r220_to_r232.json"))
ch = [r for r in fm if r["status"] == "changes"]
parts = []
for r in ch:
    tag = "R232" if m232[r["bin"]]["old_id"] != m232[r["bin"]]["new_id"] else "older map"
    parts.append(f"{r['label_in_figure']} -> {r['R232_ID']} ({r['bin']}, {tag})")
p = d.find_para("Figure 4. Denitrifier functional classes and their phase dynamics.")
cid = d.comment_on_para(p, "Figure 4: Emily's 2026-10-21 version (paper3/papert3_071626.pdf p. 24) labels the "
                           "x-axis by MAG bin, so it needs no change for GTDB R232. The older image embedded here "
                           "could not be re-rendered (no generator for this 5-class version in the repository); its "
                           "x-axis 'MAG iteration ID' labels come from the workbook's older ID map; "
                           f"with the current GTDB R232 IDs (keyed by bin) {len(ch)} of its {len(fm)} labels change: "
                           + "; ".join(parts) + ". Classes, genes and values are unchanged (the 5-class calls do not "
                           "depend on taxonomy). If this ID-labelled version is kept, re-render it from the updated 'MAG "
                           "iteration ID' column of Tables S1/S7 (MAG_target_gene_filled_UPDATED (1).xlsx).")
d.log.append(dict(label="Figure 4 label-map comment", old="", new=f"{len(ch)} label changes", comment=cid))

d.save()
nbad = fix_settings_track_revisions(U)
print("settings: removed invalid w:trackChanges x", nbad, "; w:trackRevisions present")
for x in d.log:
    print(json.dumps(x, ensure_ascii=False)[:300])
json.dump(d.log, open(f"{P}/apply_scratch/docs/main_docx_changes.json", "w"), indent=1, ensure_ascii=False)
