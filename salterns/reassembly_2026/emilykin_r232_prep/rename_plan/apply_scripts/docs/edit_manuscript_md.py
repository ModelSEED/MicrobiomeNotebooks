"""Docs stage, step B: direct edits of MANUSCRIPT_methods_results_discussion.md for GTDB R232 (docs/PLAN.md 4.4).
Each edit is an exact, unique sentence-level substring (asserted to occur exactly once), never a global string
replace.  Numbers come from the regenerated repo tables (compute_manuscript_stats_291 re-run in scratch,
Table_S6_skani_novelty.tsv, gtdbtk.bac120.summary.tsv, the per-cell render log)."""
import os, sys
R = os.path.expanduser("~/Documents/EmilyKin")
F = f"{R}/MANUSCRIPT_methods_results_discussion.md"
raw = open(F, "rb").read()
assert b"\r\n" not in raw
t = raw.decode("utf-8")

E = []  # (label, old, new)

# --- line 45: GTDB-Tk version, command, method counts; delete the version-mismatch caveat
old45_start = "The 291 dereplicated MAGs were taxonomically classified with GTDB-Tk v2.7.2 (`classify_wf`;"
i = t.index(old45_start); j = t.index("\n", i)
old45 = t[i:j]
assert old45.endswith("but the resulting taxonomy corresponds to GTDB r220 and is reported as such throughout.")
new45 = ("The 291 dereplicated MAGs were taxonomically classified with GTDB-Tk v2.8.0 (`classify_wf`; "
         "`gtdbtk classify_wf --batchfile <291 MAG FASTA files> --out_dir <out> --cpus 24 --pplacer_cpus 8 "
         "--min_perc_aa 10 --scratch_dir <scratch>`) against GTDB release R232 (the GTDB-Tk `release232` reference "
         "data; ANI pre-screen with skani v0.4.0). Classification was performed locally rather than deferred to any "
         "external service (2026-11-11; summary `gtdbtk.bac120.summary.tsv`; the run outputs are archived in KBase "
         "report 262319/892/1). It supersedes the earlier GTDB-Tk v2.7.2 classification against r220; the MAG "
         "iterativeIDs that changed with it are listed in `mag_iterativeID_r220_to_r232.csv` (§2.7). No archaeal "
         "placements were produced (only a bac120 summary was written), and no genomes failed. Across the 291 MAGs, "
         "122 were classified by the ANI pre-screen (\"ani_screen\"), 97 by topology and ANI (\"defined by topology "
         "and ANI\"), 51 by topology alone (\"fully defined by topology\"), and 21 by relative evolutionary divergence "
         "(RED, \"taxonomic novelty determined using RED\").")
E.append(("L45 GTDB-Tk v2.8.0 / R232, method counts, caveat deleted", old45, new45))

# --- line 61: naming-case counts
E.append(("L61 case counts",
          "producing 155 ASV-inherited, 62 GTDB-genus, and 74 lowest-GTDB-label MAG identifiers (291 total, no blanks).",
          "producing, on the GTDB R232 taxonomy (§2.5), 138 ASV-inherited, 70 GTDB-genus, and 83 lowest-GTDB-label MAG "
          "identifiers (291 total, no blanks). Re-running the unchanged rule on R232 changed 51 MAG identifiers relative "
          "to the r220-based map (155/62/74), 50 of them at the root; no ASV identifier changed, and the old-to-new pairs "
          "are listed by bin in `mag_iterativeID_r220_to_r232.csv`."))

# --- line 87 and 153: CAN_5_bin.40 genus
E.append(("L87 CAN_5_bin.40 genus", "(CAN_5_bin.40, *Giesbergeria*), retained as a low-confidence Clade I call",
          "(CAN_5_bin.40, *Simplicispira*), retained as a low-confidence Clade I call"))
E.append(("L153 CAN_5_bin.40 genus (1)", "and CAN_5_bin.40 (*Giesbergeria*), the latter a low-confidence call",
          "and CAN_5_bin.40 (*Simplicispira*), the latter a low-confidence call"))
E.append(("L153 CAN_5_bin.40 genus (2)", "being *Giesbergeria* (CAN_5_bin.40, 9 genes)",
          "being *Simplicispira* (CAN_5_bin.40, 9 genes)"))

# --- line 103: skani on R232
E.append(("L103 skani R232",
          "Genome-wide average nucleotide identity of each MAG to its closest GTDB r220 representative was computed with "
          "skani v0.3.2 (`skani search` against the prebuilt GTDB r220 skani database, `--min-af 0 -n 3`).",
          "Genome-wide average nucleotide identity of each MAG to its closest GTDB R232 species representative was "
          "computed with skani v0.4.0 (`skani search` against the GTDB-Tk `release232` skani database of 199,923 species "
          "representatives, `--min-af 0 -n 3`; `reprocess_291/skani_novelty_r232/`), with best-hit taxonomy from the "
          "R232 `gtdb_taxonomy.tsv`."))
E.append(("L103 -n 3 caveat (coasm_bin.535)",
          "These calls populate the standalone `Table_S6_skani_novelty.tsv` (291 MAGs) and the `Table S6_ANI` sheet of "
          "the supplementary workbook.",
          "These calls populate the standalone `Table_S6_skani_novelty.tsv` (291 MAGs) and the `Table S6_ANI` sheet of "
          "the supplementary workbook. Because `-n 3` keeps only the three highest-ANI hits per MAG, a MAG whose top "
          "hits are short high-identity fragments is called `low_AF_unresolved` even when a full-length hit exists; this "
          "affects one MAG, coasm_bin.535, which GTDB-Tk assigns to *Novosphingobium* sp017302615 (98.4% ANI, AF 0.90)."))

# --- line 107: Fig S4 methods (relabel-only option; tree inputs stay r220)
E.append(("L107 CAN_4_bin.208 genus",
          "and the now 100%-complete *Azospira_A* MAG CAN_4_bin.208,",
          "and the now 100%-complete *Proximibacter* (r220 *Azospira_A*) MAG CAN_4_bin.208,"))
E.append(("L107 61 r220 references (R232 names)",
          "of 65 selected MAGs plus 61 GTDB r220 reference genomes over the concatenated bac120 marker alignment",
          "of 65 selected MAGs plus 61 GTDB r220 reference genomes (labelled with their GTDB R232 names; GCA_017302355.1, "
          "which is not in R232, keeps its r220 name) over the concatenated bac120 marker alignment"))
E.append(("L107 reference choice (r220 Azospira_A -> R232 Proximibacter restrictus)",
          "the largest-cluster *Azospira_A* representative (GCA_030917225.1) has a fragmentary bac120 row (65% of masked "
          "columns occupied), so the *Azospira_A* representative sharing the most identical aligned residues with the MAG "
          "was used instead: *Azospira_A restricta* (GCF_016858125.1;",
          "the largest-cluster r220 *Azospira_A* representative (GCA_030917225.1) has a fragmentary bac120 row (65% of "
          "masked columns occupied), so the r220 *Azospira_A* representative sharing the most identical aligned residues "
          "with the MAG was used instead: *Azospira_A restricta* (GCF_016858125.1; *Proximibacter restrictus* in GTDB R232;"))

# --- line 123: taxonomy results
old123_start = "GTDB-Tk (r220) classified all 291 MAGs as Bacteria"
i = t.index(old123_start); j = t.index("\n", i)
old123 = t[i:j]
assert old123.endswith("i.e. only about a third of the catalogue is confidently conspecific with a named GTDB genome.")
new123 = ("GTDB-Tk (R232) classified all 291 MAGs as Bacteria (no Archaea, no failed genomes), spanning 15 phyla (14 if "
          "*Myxococcota_A* is merged with *Myxococcota*). The catalogue was dominated by Pseudomonadota (131 MAGs), "
          "Bacteroidota (62), and Chloroflexota (36), which together accounted for 229 of 291 MAGs (79%); the remaining "
          "phyla were Patescibacteriota (21; Patescibacteria in r220), Verrucomicrobiota (10), Planctomycetota (10), "
          "Acidobacteriota (8), Spirochaetota (4), Myxococcota (2), Bdellovibrionota (2), and one MAG each of "
          "Myxococcota_A, Gemmatimonadota, Hydrogenedentota, Armatimonadota and Elusimicrobiota. Rank resolution was "
          "complete through order (291/291), with 290 MAGs resolved at family, 283 at genus (183 distinct genera), and "
          "122 at species (169 unnamed). The high genus-level but lower species-level resolution, together with the 21 "
          "MAGs flagged as taxonomic novelty by RED, is typical of GTDB-Tk placement of environmental MAGs. Genome-wide "
          "ANI to the closest GTDB R232 representative (skani; §2.13, Table S6) corroborated this novelty: of the 291 "
          "MAGs, 121 shared ≥95% ANI with a representative (same species), 76 fell in the 85–95% band (putatively "
          "species-distinct), 28 were <85% (stronger divergence; 23 of them with no skani hit), and 66 were unresolved "
          "because <15% of the MAG aligned to any representative — i.e. about two-fifths of the catalogue (42%) is "
          "confidently conspecific with a GTDB species representative.")
E.append(("L123 taxonomy results (15 phyla, ranks, RED, S6 buckets)", old123, new123))

# --- line 141
E.append(("L141 case counts",
          "The per-MAG iterativeID assignment produced 155 ASV-inherited, 62 GTDB-genus, and 74 lowest-GTDB-label identifiers.",
          "The per-MAG iterativeID assignment produced 138 ASV-inherited, 70 GTDB-genus, and 83 lowest-GTDB-label "
          "identifiers on GTDB R232 (155/62/74 on r220)."))

# --- line 147: coasm_bin.185 species
E.append(("L147 coasm_bin.185 species", "coasm_bin.185 (*Accumulibacter* sp000585075;",
          "coasm_bin.185 (*Accumulibacter regalis*;"))

# --- line 157: genome-level species of the ppk1-typed MAGs (D7: coasm_bin.476 is not species-novel)
E.append(("L157 GTDB R232 species agree with the ppk1 clades; coasm_bin.476 not species-novel",
          "and coasm_bin.476 in clade ID (100/100). The remaining two",
          "and coasm_bin.476 in clade ID (100/100). At genome level GTDB R232 assigns these four to *Accumulibacter "
          "regalis* (coasm_bin.185, 97.58% ANI), *A. delftensis* (coasm_bin.250, 98.60%), *A. phosphatis* "
          "(coasm_bin.347, 98.14%) and *Accumulibacter* sp032229525 (coasm_bin.476, 98.68% ANI, AF 0.84; a species "
          "representative added to GTDB in R226), so the three named-clade calls agree with the GTDB species and "
          "coasm_bin.476 is not novel at species level. The remaining two"))

# --- line 181: delete the Discussion caveat
E.append(("L181 delete the r220/v2.7.2 caveat",
          " Finally, the GTDB taxonomy is r220 although the GTDB-Tk binary (v2.7.2) natively expects r232; the runs "
          "completed without failure, but all GTDB labels should be cited as r220.",
          ""))

# --- line 191: Figure 1 legend genus list (from the regenerated figure's printed table)
E.append(("L191 Fig 1 band list",
          "(*Competibacter*, *Accumulibacter*, *Thauera*, UBA5532, *Pseudoxanthomonas*, *Azonexus*, JAEULV01, Unclassified, "
          "OLB5, *Leptovillus*, UBA1573, UBA5066, *Dokdonella*, JJ008, SCMT01, SpSt-398, *Desulfobacillus*, *Zoogloea*, "
          "*Azospira*, *Giesbergeria*)",
          "(*Competibacter*, *Accumulibacter*, *Thauera*, UBA5532, *Pseudoxanthomonas*, *Azonexus*, JAEULV01, OLB5, "
          "*Leptovillus*, UBA1573, UBA5066, *Dokdonella*, *Lysobacter* (GTDB *Lysobacter_I*), JJ008, SpSt-398, "
          "*Desulfobacillus*, *Zoogloea*, *Proximibacter*, *Simplicispira*, Unclassified; GTDB R232 genera)"))

# --- line 225: Table S5 legend
E.append(("L225 Table S5 coasm_bin.476 species",
          "coasm_bin.476 = ID (100/100);",
          "coasm_bin.476 = ID (100/100; genome in the GTDB R232 species *Accumulibacter* sp032229525, 98.68% ANI);"))

# --- line 229: Figure S4 legend
E.append(("L229 Fig S4 references (R232 names)",
          "of 65 selected MAGs plus 61 GTDB r220 references — each MAG's",
          "of 65 selected MAGs plus 61 GTDB r220 references (labelled with their GTDB R232 names; GCA_017302355.1, absent "
          "from R232, keeps its r220 name) — each MAG's"))
E.append(("L229 Fig S4 root", "median 100), rooted on the Patescibacteria clade.",
          "median 100), rooted on the Patescibacteriota clade."))
E.append(("L229 Fig S4 skani release", "no skani hit to any GTDB r220 representative, i.e. <~80% ANI)",
          "no skani hit to any GTDB R232 representative, i.e. <~80% ANI)"))
E.append(("L229 Fig S4 reference labels", "and references with their GTDB species.",
          "and references with their GTDB R232 species."))

# --- line 231: Table S6 legend
E.append(("L231 Table S6 release + buckets",
          "Per-MAG (291) closest GTDB r220 representative, skani ANI, query/reference alignment fractions, ANI class "
          "(≥95 / 85–95 / <85 / low_AF_unresolved) and a novelty interpretation. Distribution: 96 same-species (≥95%), "
          "70 putatively species-distinct (85–95%), 48 stronger divergence (<85%), 77 unresolved by low alignment fraction.",
          "Per-MAG (291) closest GTDB R232 species representative (skani v0.4.0), skani ANI, query/reference alignment "
          "fractions, ANI class (≥95 / 85–95 / <85 / low_AF_unresolved) and a novelty interpretation. Distribution: 121 "
          "same-species (≥95%), 76 putatively species-distinct (85–95%), 28 stronger divergence (<85%; 23 with no skani "
          "hit), 66 unresolved by low alignment fraction (r220: 96 / 70 / 48 / 77)."))

dry = "--dry-run" in sys.argv
log = []
for label, old, new in E:
    n = t.count(old)
    assert n == 1, (label, n)
    t = t.replace(old, new, 1)
    log.append((label, len(old), len(new)))
for l in log:
    print("%-80s %6d -> %6d chars" % l)
if not dry:
    open(F, "wb").write(t.encode("utf-8"))
    print("written", F)
