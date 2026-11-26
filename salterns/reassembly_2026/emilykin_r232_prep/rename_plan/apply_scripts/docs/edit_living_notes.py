"""Docs stage, step E: living notes and the novel-bacteria module text for GTDB R232 (D5: living files are updated
where the string means the MAG; logs, run records, hand-off notes, frozen files are not).  Exact, unique substring
edits (asserted count == 1) plus bin-keyed cell patches; line endings preserved per file."""
import csv, io, json, os, re, sys
R = os.path.expanduser("~/Documents/EmilyKin")
sys.path.insert(0, f"{R}/scripts")
from _binkeyed import patch_table  # noqa: E402

DRY = "--dry-run" in sys.argv
LOG = []


def edit(rel, pairs):
    path = f"{R}/{rel}"
    raw = open(path, "rb").read()
    crlf = b"\r\n" in raw
    t = raw.decode("utf-8")
    if crlf:
        t = t.replace("\r\n", "\n")
    for old, new in pairs:
        n = t.count(old)
        assert n == 1, (rel, old[:70], n)
        t = t.replace(old, new, 1)
    if crlf:
        t = t.replace("\n", "\r\n")
    if not DRY:
        open(path, "wb").write(t.encode("utf-8"))
    LOG.append((rel, len(pairs)))


R232_NOTE_FIG1 = """## 7. GTDB R232 relabel (2026-11-12)

The figure was re-rendered with GTDB R232 data (the R232 rename of the EmilyKin MAGs; decision D3 "relabel, do
not rebuild"). The RAxML tree, the alignment and the 61 reference genomes are unchanged (they remain the r220
selection described above, so the r220 names in sections 2-6 are the names at the time of the build).
- MAG tips: 8 labels change with the R232 iterativeIDs (Ca_Competibacter.32 -> Contendibacter.1_m, Azospira.3 ->
  Proximibacter.2_m, Flavipsychrobacter.6_m -> Taibaiella.5, Ignavibacteria.7 -> Tepidaquicellaceae.3_m,
  Ignavibacteria.5 -> Tepidaquicellaceae.1_m, Ignavibacteria.21 -> Tepidaquicellales.1_m, Ca_Falkowbacteria.7 ->
  Patescibacterium.1_m, Paceibacteria.1_m -> Minisyncoccia.1_m); MAG taxonomy from the R232
  `gtdbtk.bac120.summary.tsv`.
- Reference tips: R232 names from the GTDB-Tk release232 taxonomy (`render_fig1_tree.py --ref-taxonomy`); 9 labels
  change (e.g. *Azospira_A restricta* -> *Proximibacter restrictus*, *Accumulibacter* sp000585075 -> *A. regalis*,
  *A. phosphatis_G* -> *A. delftensis*, JARSSO01 sp018262135 -> *Contendibacter* sp018262135). GB_GCA_017302355.1
  (*Azonexus* sp017302355) is not in R232 and keeps its r220 name with "(not in GTDB R232)".
- Red circles: from the R232 Table S6 skani rebuild (`reprocess_291/skani_novelty_r232/`); 10 MAGs change ANI class.
- Rooting: unchanged 7-tip clade, now named Patescibacteriota (subtitle "rooted on Patescibacteriota clade (7 tips)").
- Re-running `build_fig1_alignment.py` or `score_ref_candidates.py` now reads R232 inputs and would not reproduce the
  r220 alignment above; use the 9e087c3 files for that.

"""

# ---- fig1_tree/README.md
edit("fig1_tree/README.md", [
    ("CAN_4_bin.208 (Azospira.3) is now 100 % complete",
     "CAN_4_bin.208 (Proximibacter.2_m; Azospira.3 before the GTDB R232 rename) is now 100 % complete"),
    ("- The tree is rooted at the midpoint of the Patescibacteria stem.",
     "- The tree is rooted at the midpoint of the Patescibacteria stem (Patescibacteriota in GTDB R232; the renderer\n"
     "  accepts both names and its subtitle names the clade by its tips' phylum)."),
    ("- Labels: MAG tips show the current iterativeID (bold); reference tips show the GTDB r220 species\n  (grey italic).",
     "- Labels: MAG tips show the current iterativeID (bold); reference tips show the GTDB species (grey italic):\n"
     "  the GTDB R232 name since 2026-11-12 (section 7), the r220 species before."),
    ("## Reproduce\n", R232_NOTE_FIG1 + "## Reproduce\n"),
])

# ---- clade_classify/NOSZ_TOP5_NOTES.md
edit("clade_classify/NOSZ_TOP5_NOTES.md", [
    ("PNG 5309 × 2416 px,\naspect 2.197, the same pixel size as the copy embedded in `SI-paper3 (3)_correction.docx`)",
     "PNG 5311 × 2416 px\nsince the GTDB R232 relabel of 2026-11-12 (5309 × 2416 before; longer rank-5 label), embedded in\n"
     "`SI-paper3 (3)_correction.docx` as a tracked image replacement)"),
    ("`CAN_3_bin.221` (Ignavibacteria.21 / OLB5, Bacteroidota, II.B/C low)",
     "`CAN_3_bin.221` (Tepidaquicellales.1_m, Ignavibacteria.21 before GTDB R232 / OLB5, Bacteroidota, II.B/C low)"),
    ("| 5 | CAN_3_bin.221 | Ignavibacteria.21       |", "| 5 | CAN_3_bin.221 | Tepidaquicellales.1_m   |"),
    ("Ignavibacteria.21 / OLB5 (22.2)", "Tepidaquicellales.1_m / OLB5 (22.2)"),
])

# ---- clade_classify/CLADE_I_II_DISCREPANCY.md
edit("clade_classify/CLADE_I_II_DISCREPANCY.md", [
    ("| **CAN_5_bin.40** | *Giesbergeria* |", "| **CAN_5_bin.40** | *Simplicispira* (r220 *Giesbergeria*) |"),
    ("(the 6th, *Giesbergeria*, is a borderline", "(the 6th, *Simplicispira* CAN_5_bin.40, is a borderline"),
    ("and the *Giesbergeria* MAG).", "and the *Simplicispira* (r220 *Giesbergeria*) MAG CAN_5_bin.40)."),
])

# ---- gene_ab_figure/METHODS_nosz_clades.md
edit("gene_ab_figure/METHODS_nosz_clades.md", [
    ("`CAN_5_bin.40` *Giesbergeria* low-confidence", "`CAN_5_bin.40` *Simplicispira* (r220 *Giesbergeria*) low-confidence"),
    ("| **II.G** | 2 | *Giesbergeria*, UBA5066 | *Desulfitobacterium* + env. | *Giesbergeria* high (60), UBA5066 medium (14) |",
     "| **II.G** | 2 | *Simplicispira* (r220 *Giesbergeria*), UBA5066 | *Desulfitobacterium* + env. | *Simplicispira* high (60), UBA5066 medium (14) |"),
])

# ---- clade_classify/out/summary.md (classify_nosz.py output; bin-keyed genus cell, script not re-run)
edit("clade_classify/out/summary.md", [
    ("| CAN_5_bin.40 | Giesbergeria  |", "| CAN_5_bin.40 | Simplicispira |"),
])

# ---- phase_rpkm_chord/README.md
edit("phase_rpkm_chord/README.md", [
    ("→ 2 Clade I (Ca_Competibacter.51, Giesbergeria.1) + 11 Clade II.",
     "→ 2 Clade I (Ca_Competibacter.51, Simplicispira.1 = CAN_5_bin.40 under GTDB R232; the chord figures,\n"
     "  not re-rendered, still label it Giesbergeria_m.1) + 11 Clade II."),
])

# ---- r13 annotation notes (the MAGs' GTDB phylum / genus; KBase object renames)
R13_NOTE = ("> **GTDB R232 update (2026-11-12).** The MAG taxonomy is now GTDB R232 (GTDB-Tk v2.8.0): phylum "
            "Patescibacteria is Patescibacteriota (the same 21 MAGs; the CPR flag and envelope calls are unchanged), "
            "and 51 MAG iterativeIDs changed (`mag_iterativeID_r220_to_r232.csv`). In KBase ws 262319 the 51 "
            "Assemblies and Genomes were renamed by object ID (versions and checksums unchanged) and the 291 "
            "Assemblies carry the R232 lineage (v2); the RAST feature IDs inside each Genome keep the old name, so join "
            "on bin. Member (ASV) names are unchanged. `taxonomy.tsv` and `panel/mag_panel.tsv` were rebuilt and "
            "`calls/summary.tsv` and `rast/genomes/_manifest.tsv` patched by bin; see the dated freeze-2a record in "
            "`panel/PANEL_FROZEN.md`. Notes below that quote r220 names describe runs made before the update.\n\n")
# ANNOTATIONS.md and plan_r13/track2/GENOME_OBJECT_SPEC.md are NOT edited: both are in the ASV md5 baseline
# (apply_log/asv_md5_before.tsv) and hold no MAG ID; the R232 note goes into crosswalk_NOTES.md instead.
txt = open(f"{R}/community_modeling/analysis/r13/annotations/crosswalk_NOTES.md", encoding="utf-8").read()
first_hdr_end = txt.index("\n") + 1
edit("community_modeling/analysis/r13/annotations/crosswalk_NOTES.md", [(txt[:first_hdr_end], txt[:first_hdr_end] + "\n" + R13_NOTE.rstrip("\n") + "\n")])
edit("community_modeling/analysis/r13/annotations/crosswalk_NOTES.md", [
    ("as expected for Patescibacteria.", "as expected for Patescibacteriota (CPR)."),
    ("(36 Chloroflexota + 21 Patescibacteria + coasm_bin.207)", "(36 Chloroflexota + 21 Patescibacteriota + coasm_bin.207)"),
    ("**21 Patescibacteria, `unclear_CPR`**", "**21 Patescibacteriota (r220: Patescibacteria), `unclear_CPR`**"),
    ("Genome name = `<setname>_<iterativeID>`.",
     "Genome name = `<setname>_<iterativeID>`. (GTDB R232 update, 2026-11-12: 51 Genome objects were renamed to their "
     "R232 iterativeID by object ID; the feature IDs inside them keep the old name. Join on bin.)"),
])
edit("community_modeling/analysis/r13/annotations/kofam_NOTES.md", [
    ("The 16 code-11 Patescibacteria are all ok.", "The 16 code-11 Patescibacteriota (r220: Patescibacteria) are all ok."),
    ("coasm_bin.479 (a 0.5 Mb UBA9973 Patescibacteria genome)",
     "coasm_bin.479 (a 0.5 Mb Patescibacteriota genome; order UBA9973 in r220, Phycocordibacterales in R232)"),
])
edit("community_modeling/analysis/r13/annotations/panel_NOTES.md", [
    ("Competibacter CAN_4_bin.116 and JARSSO01 have hemY only",
     "Competibacter CAN_4_bin.116 and Contendibacter coasm_bin.365 (r220 JARSSO01) have hemY only"),
])
# ---- novel-bacteria module: bin-keyed RED / method / novelty from the R232 summary + notes
gt = {r["user_genome"]: r for r in csv.DictReader(open(f"{R}/gtdbtk.bac120.summary.tsv", newline=""), delimiter="\t")}
idmap = json.load(open(f"{R}/mag_iterativeID_old_to_new.json"))
nb = list(csv.DictReader(open(f"{R}/gene_ab_figure/data/novel_bacteria_summary.tsv", newline=""), delimiter="\t"))
upd = {"RED": {}, "novelty": {}, "method": {}}
for r in nb:
    b = r["bin"]; g = gt[b]
    sp = g["classification"].split(";")[6][3:]
    upd["method"][b] = g["classification_method"]
    upd["RED"][b] = f"{float(g['red_value']):.3f}" if g["red_value"] not in ("", "N/A") else "N/A"
    upd["novelty"][b] = (f"not species-novel in GTDB R232 (s__{sp}; ANI {float(g['closest_genome_ani']):.2f}%)" if sp
                         else "species-novel (s__; no GTDB R232 species; topology and ANI)")
if not DRY:
    res = patch_table(f"{R}/gene_ab_figure/data/novel_bacteria_summary.tsv", "\t", "bin", upd)
    LOG.append(("gene_ab_figure/data/novel_bacteria_summary.tsv", res))

# .md tables: N1 rows keyed by bin (ID + RED cells); N2/N3 rows keyed by ID through N1's ID -> bin map
md = f"{R}/gene_ab_figure/data/novel_bacteria_summary.md"
t = open(md, encoding="utf-8").read()
id2bin = {}
out = []
for line in t.split("\n"):
    if line.startswith("| ") and not line.startswith("|---"):
        cells = line.split("|")
        c1, c2 = cells[1].strip(), (cells[2].strip() if len(cells) > 2 else "")
        if c2 in idmap:                       # Table N1: | iterativeID | bin | ... | RED |
            id2bin[c1] = c2
            cells[1] = f" {idmap[c2]} "
            cells[-2] = f" {upd['RED'][c2]} "
        elif c1 in id2bin:                    # Tables N2/N3: | iterativeID | ...
            cells[1] = f" {idmap[id2bin[c1]]} "
        line = "|".join(cells)
    out.append(line)
t2 = "\n".join(out)
NB_NOTE = ("\n> **GTDB R232 update (2026-11-12).** The seven MAGs were chosen on GTDB r220 (no species, placement fully "
           "defined by topology). Under GTDB R232 (GTDB-Tk v2.8.0) coasm_bin.260 (Ca_Sarcinithrix.1) is assigned to the "
           "species JAAEKA01 sp049294015 (ANI 98.33%), so six of the seven remain species-novel; the other six are now "
           "placed by \"topology and ANI\" with new RED values (Table N1). Conversely CAN_5_bin.70 (Azonexus.12) lost "
           "its r220 species, because its reference genome left GTDB. Figure N1 (`10_novel_denitrification_bubble.py`) "
           "now draws the R232 set (with CAN_5_bin.70, without coasm_bin.260), while Figure N2 and Tables N2-N3 keep "
           "the r220 set; which set the module uses is an author decision. IDs follow the current map (CAN_1_bin.210 "
           "is Dokdonella.8; Dokdonella.14 was its older ID).\n")
anchor = "i.e. each represents a candidate novel species.\n"
assert t2.count(anchor) == 1
t2 = t2.replace(anchor, anchor + NB_NOTE, 1)
if not DRY:
    open(md, "w", encoding="utf-8").write(t2)
LOG.append(("gene_ab_figure/data/novel_bacteria_summary.md", "N1 ID/RED + N2/N3 ID by bin, R232 note"))

SUPP_NOTE = ("\n> **GTDB R232 update (2026-11-12).** This module was written on the 276-MAG catalogue with GTDB r220 and is "
             "kept as written. Under GTDB R232 coasm_bin.260 (*Ca_Sarcinithrix.1*) belongs to the species JAAEKA01 "
             "sp049294015 (98.33% ANI), so six of the seven remain species-novel, and CAN_5_bin.70 (*Azonexus.12*) lost "
             "its r220 species. None of the seven MAG IDs changes in R232; *Dokdonella.14* here is CAN_1_bin.210, now "
             "*Dokdonella.8*. See `gene_ab_figure/data/novel_bacteria_summary.md` for the R232 values; whether the "
             "module keeps the r220 set of seven is an author decision.\n")
edit("novel_bacteria_supplementary.md", [
    ("project's novel organisms.*\n", "project's novel organisms.*\n" + SUPP_NOTE),
])
for x in LOG:
    print(x)
