# R232 review of the saltern manuscript (2026-10-21, revised 2026-10-26, rebased 2026-11-04, revised 2026-11-08, 2026-11-16, rebased and revised 2026-11-19 and 2026-11-20, revised 2026-11-23 and 2026-11-24)

**File:** `Metabolic Manuscript - R232 review (tracked changes).docx`.
- **Base:** since 2026-11-20 it is built on the Google Doc as of 20 November 2026, through the verified local mirror `Metabolic Manuscript - Google Doc mirror (2026-11-20).docx` (2026-11-19: the 2026-11-17 mirror; 2026-11-04 to 2026-11-18: the 2026-11-03 mirror).
- **Tracked changes:** every edit of this review is a Word tracked change by "Claude (R232 review)". Susannah Tringe's open suggestions appear as her own tracked changes.
- **Undoing the review:** reject only the "Claude (R232 review)" changes (Word: Show Markup > Specific People, then Reject All Changes Shown). That restores the Google Doc text, tables and figures.

**Accompanying file:** `Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx`, one row per MAG (310).

**Rebuild** in two steps, from `clauded_insights/`:
1. **The mirror:** `make_gdoc_mirror_docx.py <merged.docx>`, where `<merged.docx>` is the 2 July export `Metabolic Manuscript.docx` after the docx skill's `merge_runs.py`.
   - The Google Doc changes are specified in `gdoc_sync/gdoc_sync_2026-11-03.json`.
   - The script checks every paragraph against `gdoc_sync/gdoc_views_2026-11-03.json`, the Docs API text with suggestions rejected and accepted.
   - `--date 2026-11-17` builds the newer mirror from `gdoc_sync/gdoc_{sync,views}_2026-11-17.json` (see "Google Doc sync 2026-11-17" below); without it the 11-03 mirror is rebuilt byte-identically.
2. **The review:** `make_r232_review_docx.py "../Metabolic Manuscript - Google Doc mirror (2026-11-03).docx"`. Its inputs are:
   - `correlational/r232/` and `data/`;
   - `r232_review_edits_2026-10-26.json`, with 41 tracked edits and comments, each with its rationale and source (edits revised on 2026-11-08 carry that date, which their tracked changes show);
   - `figure7_r232.png`.

Paragraph numbers (P…) in these notes are those of the July export, which the mirror keeps.

> **Rebased on the Google Doc (2026-11-04).** The Drive connector cannot export the Google Doc as .docx because it is over the size limit. The mirror therefore rebuilds it from the July export plus every change the Docs API reports:
> - the 8 suggestions accepted on 27 August;
> - Susannah's 37 suggestions of 3 November (22 edits), as her tracked changes;
> - the comment changes since July: the thread resolved on 27 August, and the 13 threads and 1 reply added since then.
> - the 4 threads Susannah added on 5 November (Figure 1 colours, Figure 2 scale and detection, and two on "betaine" vs "glycine betaine"), added on 2026-11-09. The Google Doc text itself has not changed since 3 November.
>
> With the suggestions rejected, and again with them accepted, all 491 paragraphs match the Google Doc. The 20 figures are identical. The review of the 70 open comment threads is in `Google Doc comments - review (2026-11-04).md`.

> **Author fixes (2026-10-21).** The author reports fixing the three draft errors (Table 3, P69, P35) in their own copy. That copy is neither the local .docx nor the Google Doc, which still has the original text on 3 November. The tracked edits at P35 and P69 now duplicate those fixes. The R232 Table 3 replacement is a separate update. A hand fix of the old Table 3 also needs new Unrestored values, because that column pooled the R2A cores (see the Table 3 row below).

## Google Doc sync 2026-11-17

The Google Doc changed on 2026-11-16 and 2026-11-17 (Drive modifiedTime 2026-11-17T23:20:24Z; read through the Docs and Drive APIs on 2026-11-19). The Docs API attributes every change to A. Freiburger. All changes are in the title block, the Abstract and the first five Introduction paragraphs. From P25 ("Methods") on, all 466 paragraphs are unchanged.
- **Suggestions accepted (27 of Susannah's):**
  - P8: Kyle Hartman added to the author list (her 11-03 suggestion).
  - P23–P24: all 24 of her 2026-05-07 suggestions, plus her 11-03 "in marine-influenced systems" (P23) and "enabling" (P24).
  - P24 now reads exactly as her suggested text. The 19 other 11-03 edits (Methods P32–P57) are still open.
- **Abstract (P18) rewritten and split** into two paragraphs, at "Here we explore…". This includes his own suggestions of 11-16 (2 accepted, 1 rejected).
  - The methanogenesis/MePn mechanism clause moves from the first part into "centered around …".
  - "increase carbon sequestration" → "maintain sequestered carbon".
  - "Key experimental findings include: site-specific metabolomes reflect adaptations to salinity, and…" → "In addition to experimentally revealing salinity-dependent metabolic changes, a key *in vitro* finding is that…".
  - The "primary methanogenic substrate" claim now follows the in vitro result ("…increases N₂O production, which supports the role of methylphosphonate being a primary methanogenic substrate in these environments").
  - "These findings were corroborated by metabolic modeling, whereby the trace consumption of methylphosphonic acid reproduced the experimentally measured methane emissions, providing further evidence for…" → "This was corroborated by metabolic modeling, where trace consumption of methylphosphonic acid was sufficient to reproduce the experimental methane emissions."
- **Introduction P20–P22 restructured** by direct edits:
  - P20 absorbs "(i.e., > 7% salt) (15)" and the compatible-solute list, as parentheticals on salt-in/salt-out. It contains a typo in the Google Doc: "use the the".
  - The compatible-solute paragraph (P22) now opens with its former last sentence. It then reads "For this study, there are two notable compatible solutes. Firstly, glycine betaine (GB)… Secondly, the multifunctional disaccharide trehalose…", and it moves ahead of P21.
  - P21 is cut to "The biogeochemistry metabolism of organisms in these hypersaline environments (16) is notably different in pathways such as … do not generate enough energy …". As written in the Google Doc, this sentence is ungrammatical.
- **P23 opening rewritten** after the acceptances: "Previous studies (31–33), including our own of hypersaline unrestored salterns (34), showing that methane emissions from marine-influenced environments are inversely related to salinity in a log-linear relationship. This log-linear decline is attributed to the sulfate in seawater which enables sulfate reducers…".
  - The first sentence has no main verb.
  - The earlier contrast, "we've previously observed elevated methane emissions in hypersaline unrestored salterns (34)", is gone. (34) is now cited for the log-linear decline.
- **Comments:** now 242 threads: 72 open (60 anchored, 12 orphaned) and 170 resolved. Drive's lastCommentedTime (11-05) does not reflect these changes.
  - Resolved on 11-17: the trehalose thread (Word ids 2, 3) and "Need to confirm that this reference supports my statement about sulfate" on "15" (id 6).
  - Orphaned: "Methylphosphonates are also organic molecules, not minerals" (id 8) is still open, but its anchor was deleted by the accepted suggestion.
  - Opened and resolved on 11-16: "This sounds tangential … can it be safely removed from the abstract", on the metabolome clause, which is now gone.
  - The mirror goes from 83 to 79 comments. The "methanogens" thread is re-anchored in the rewritten P23.
  - Checked against the export's markers: the 66 comments on plain text match exactly. The other 13 (on images, equations or ranges across paragraphs) could not be compared automatically; they are unchanged from the 11-03 mirror.
- **New files:**
  - `clauded_insights/gdoc_sync/gdoc_views_2026-11-17.json` (492 paragraphs);
  - `clauded_insights/gdoc_sync/gdoc_sync_2026-11-17.json`, the 11-03 spec plus `accept_all_revisions`, `accepted_suggestions` and `direct_edits` (rewrites carry the Docs API runs, including italics, sub/superscripts and Zotero citation links);
  - `Metabolic Manuscript - Google Doc mirror (2026-11-17).docx`, built with `make_gdoc_mirror_docx.py <merged.docx> --date 2026-11-17`.
- **Validation:** all 492 paragraphs match the Google Doc with suggestions rejected and with them accepted. The rewritten paragraphs' formatting matches the API. Outside P8 and P18–P24, the document XML, the 3 tables and the 20 figures are identical to the 11-03 mirror. The docx validator reports no new errors. The 11-03 mirror still rebuilds byte-identically.
- **Anchor check for `make_r232_review_docx.py`** (as found by the sync; resolved 2026-11-19, see the next section):
  - **Indices shift.** The Abstract split adds a paragraph, so every paragraph from P19 on moves down by one: P[29]→30, P[30]→31, P[35]→36 … P[175]→176, including the equations P[70]/[72]/[75]→71/73/76, P[128]→129, P[162]→163, and every `paragraph` index in `r232_review_edits_2026-10-26.json` except 18. The paragraphs themselves are unchanged.
  - **Unaffected:** `tables[2]`, `first_par` and P[18]'s position.
  - **Run unchanged on the new mirror,** the script picks its style templates from the wrong paragraphs and inserts the Methods heading in the wrong place. It then stops at the first `tracked_replace` (P[35], "anchor occurs 0x").
  - **Tracked-replace check:** of the 44 `old` strings (5 in `main()`, 39 in the edits JSON), 43 still match once the index is shifted. One is broken: the edits-JSON replace at P18, "These findings were corroborated by metabolic modeling, whereby…". That sentence no longer exists. Its rewritten form is in the new P19, and the "primary methanogenic substrate" claim it softened now sits on the in vitro sentence.
  - **P18 review comment:** it still anchors, but it discusses that modeling clause, which has moved to P19.
  - **Other review edits:** none touch P19–P24. Paragraph numbers in these notes stay July-export numbers, so add 1 from P19 on when reading the 11-17 mirror.

## 2026-11-24: per-gram framing of the models; curated methanogen; PrmA block

**Why the framing changed.** The model-vs-measured comparison multiplied CH4 per g of biomass formed by standing biomass and by 24 h d⁻¹. That assumes every member grows at 1 h⁻¹.
- **Affected:** the magnitudes (14 % / 55 % of measured, "one seventh to one half", betaine 2.6×).
- **Unaffected:** the rankings.

**The new framing** (`reassembly_2026/modelling_2026/simulations/`: `per_gram_yields.py`, `specific_rates.py`, `figure7_pergram.py`).
- **Model output:** CH4 per g of biomass formed.
- **Comparison:** the growth rate each guild's measured biomass would need to supply the measured flux, and for the methanogen the required specific methanogenic activity.
- **Biomass:** relative abundance × soil DNA × layer thickness, with DNA at 10 % of dry biomass. It is used as measured, per Andrew. The text states that at 3 % every required rate is 3.3× smaller.

**Yields** (mmol CH4 g⁻¹ biomass formed):

| Guild | Yield |
|---|---|
| MPn consumers, ratio 0.07 | 0.043 |
| MPn consumers, ratio 0.33 | 0.164 |
| MPn consumers, MPn the only P source | 0.66 |
| Betaine reducers (0.06 CH4 per TMA) | 2.9–11.5 (48–192 mmol TMA g⁻¹) |
| *Methanolobus*, curated, on methylamines | 126 (120–159 across ion couplings) |

- **Why MPn consumers yield so little:** one CH4 per MPn assimilated, so methane follows phosphorus demand. A methanogen's methane is its catabolic product, about 2,900× more per gram.

**Doubling time needed to supply the measured CH4:**

| Guild | Doubling time |
|---|---|
| MPn consumers, ratio 0.07 | 2–15 min |
| MPn consumers, ratio 0.33 | 8–56 min |
| MPn consumers, MPn the only P source | 0.5–3.8 h |
| Betaine reducers (0.06) | 1–35 h |
| *Methanolobus* | R1 4–38 h (2–21 mmol g⁻¹ h⁻¹); R2 0.9–1.5 h (60–99 mmol g⁻¹ h⁻¹) |

- **Methylamine pools:** they could supply the measured flux at 0.9–3.4 d turnover.
- **Betaine → TMA chain:** for complete conversion the betaine reducers would need to double every 1.5–54 d; the methanogen is the bottleneck.
- **Ranking:** no guild ranks the cores (ρ −0.43 to 0.03).

**Text** (`r232_review_edits_2026-11-24.json`, generated by `make_modelling_edits_2026-11-24.py` on top of the frozen 11-23 JSON; assertions guard the wording):
- **Rewritten:** Abstract (yields and required rates; the inference is in the P18 and P100 comments); Methods P35 (curation, PrmA), P36, P38 and P40 (the explicit-rate comparison and the DNA constant); Results P69–P76 (MPn, betaine and the methanogen chain, with the loophole that MPn degraded beyond assimilatory demand is not represented); Discussion P95; Conclusions P99; the Figure 7 caption.
- **New Figure 7:** `figure7_pergram.png`.
- **Tracked deletions:** equations P39, P70, P72 and P75, plus their lead-ins.
- **New P73 comment:** suggests citing Hoehler and Jørgensen 2013 (verified: Nat Rev Microbiol 11:83–94) for implausible sediment growth.

**Models:**
- ***Methanolobus* curated** (models/methanogen/CURATION.md): 19 gene-justified changes; methanogenic on mono-, di- and trimethylamine; no MtaB/MtaC, so no methanol.
- **The KOfam-arm methanol growth was an artefact.**
- **PrmA:** rxn00844 closed in both *Roseovarius* models. The gene is propane 2-monooxygenase.

**Still open:**
- the read-based gene survey (`reassembly_2026/read_markers/`, running): whether unbinned methanogens raise the methanogen biomass, especially in R2;
- the Google Doc re-sync (the Doc changed 2026-11-23);
- the suggestion batch.

**Build:** "All validations PASSED" (492 → 1,076 paragraphs).

## 2026-11-23: modelling text and Figure 7 from the 2026 guild models

**Models:** `reassembly_2026/modelling_2026/models/README.md`.
- One model per MPn-consumer and betaine-reducer species of the final set, from the curated templates and ModelSEEDpy.
- RAST alone, chosen by the pre-registered comparison of RAST, +KOfam, +eggNOG and +Bakta.
- Kept: 38 MePn and 18 betaine models. JASFMS01.1 is dropped.
- No MePn model needs MePn to grow.
- The *Methanolobus* model is not methanogenic as built and is not simulated.

**Simulations:** `reassembly_2026/modelling_2026/simulations/README.md` (original cell-27 and betaine protocols, run on TF).

**Edits:** `clauded_insights/r232_review_edits_2026-11-23.json`, generated by `make_modelling_edits_2026-11-23.py` from the simulation outputs, with the wording's preconditions asserted. It replaces all 21 modelling entries of the 11-19 file and adds:
- the P35 Methods rewrite (model building; the KBase MS2 and methanogen-medium sentences go);
- a P37 comment that Figure S10 comes from the 2020 models;
- a P73 comment recommending deletion of the *Albimonas* sentence (no *Albimonas* *phnJ* in the final set). It is not deleted, because it ends inside an existing comment range and a citation link.

**Builder (`make_r232_review_docx.py`):**
- Figure 7 is `simulations/figure7_final.png`.
- The P70, P72 and P75 equations are swapped to the new fits, read from `consistent_basis_final.json`:
  - MePn 0.07: −0.02x + 198;
  - MePn 0.33: −0.09x + 764;
  - betaine: 0.58x + 2347.
- The interim P35 comment and the "310 dereplicated MAGs" P35 replace are removed.
- The originals paragraph no longer mentions the 2020 curated modelling sets.
- Table 4 notes give the modelled counts.
- The overview and re-assembly comments are updated.
- `delete_contents` now also track-deletes the runs of other authors' open insertions (nested `w:del` in `w:ins`). Before, Susannah's "(2022)" suggestion in P35 survived the paragraph rewrite.

**Numbers that change** (consistent areal basis, six unrestored cores):

| | Old models | New models |
|---|---|---|
| MePn, ratio 0.07 | 19 % of measured | **14 %** (5–33 % among cores; ρ = 0.03) |
| MePn, ratio 0.33 | 72 % | **55 %** (18–126 %) |
| Betaine | 9 *grdH* models; below measured in one core; "follows DNA (R² = 0.97)" | 18 models; **2.6× measured, above it in all six cores** (1.6–5.5×); ρ = 0.60, exact p = 0.24; DNA R² = 0.15, so the DNA argument is withdrawn |

- **New MePn statements:**
  - Methane per unit biomass is set by the ratio and the template biomass P content: 0.043–0.077 mmol g⁻¹ at 0.07.
  - Rhodobacteraceae give 75 % of the simulated methane, *Roseovarius* 49 %.
- **Abstract:** "about one seventh to one half".

**Decisions for Andrew** (also in the P36 comment):
1. ***Roseovarius*.2:** its annotated methane monooxygenase oxidises its MePn methane. It is kept as built; closing it gives 15 % instead of 14 %.
2. **JASFMS01.1** is dropped. rxn01345 would rescue it but fails two ATP tests.
3. **The methanogen model** is acetogenic as built and was not simulated.
4. **Template curation queue:** betaine and glycine reductase directions, Fpo, N-cycle reversibility, and the missing betaine reductase complex.
5. **Figure S10** (2020 models) is not regenerated.
6. **The *Albimonas* sentence** and ref 78.

**Build:** `make_r232_review_docx.py "../Metabolic Manuscript - Google Doc mirror (2026-11-20).docx"`. Validation: "All validations PASSED" (492 → 1,075 paragraphs). The Google Doc is unchanged since 2026-11-20T13:56Z (checked 2026-11-23).

**Figure audit (2026-11-23, later).** Every MAG- or model-dependent figure is now from the final MAG set or the 2026 models.
- **Figures 7, 8, 9, S8 and S11:** each embedded image is byte-identical to the current PNG. Figure S11 was regenerated and came out byte-identical, because its guilds input had been rewritten after it with the same content.
- **Figure S10 is regenerated** (`reassembly_2026/modelling_2026/figureS10/`): the core gap-filling per ATP test medium of the 2026 drafts. Its caption is rewritten.
- **P37 claim corrected:** aerobic L-lactate is the cheapest medium for the MePn consumers but ranks 16th of 22 for the betaine reducers. A comment notes that the betaine models are aerobic although the organisms are anaerobic fermenters.
- **Google Doc suggestion test (2026-11-24):** one test suggestion was made through the Docs API with `writeMode: SUGGEST` from Andrew's commenter-only account. It changes Methods P37 "(see Figure **S9**)" to "**S10**": the 9 is a suggested deletion and "10" a suggested insertion, both bold, id `suggest.hzsyv3lhfhtb`. It was verified in a fresh `read_doc`, so suggestions work.
  - The Doc changed again on 2026-11-23 at 22:40Z (4:40 PM CST): Andrew rewrote Introduction P24 and added "in triplicate" to the Methods field-sites paragraph (P27). The mirror must be re-synced before the batch.
  - The batch is held until the modelling magnitudes are reframed with an explicit growth rate (`community_eval/EVALUATION.md` §5).
- **Community-modelling evaluation** (`reassembly_2026/modelling_2026/community_eval/EVALUATION.md`): porting the EmilyKin or codiffusion_bioreactor pipelines is not recommended.
  - No measured input ranks the six cores' methane.
  - Two DNA depths read 0.16 µg g⁻¹ and need checking.
  - No manuscript text was added from this evaluation.
- **Not MAG-dependent:**
  - Figures 1 and S4 (IMG KO counts);
  - Figures 2, 3, S2 and S3 (NMR);
  - Figures 4, S5 and S9 (microcosm gases);
  - Figures 5 and 6 (isolates; Figure 5's caption has R232 edits);
  - Figure S6 (16S and IMG read mapping);
  - Figure S7 (DNA vs CH4);
  - Figure S1 (photo).


## 2026-11-20: the abundance, guild and co-occurrence analyses use the final MAG set; Google Doc sync 2026-11-20

**Google Doc sync 2026-11-20.** The Doc changed on 2026-11-20 (Drive modifiedTime 13:56:07Z, 7:56 AM CST); it was read through the Docs and Drive APIs the same day.
- **Text:** A. Freiburger rewrote three Introduction paragraphs directly:
  - P21 gains the citation (38);
  - the last three sentences of P23 are rewritten (the "(38)" moves to P21);
  - the first two sentences of P24 are rewritten. The second is unfinished in the Doc: "…although other methylated substrates such as TMA were ."
- **Unchanged:** the other 489 paragraphs and the paragraph count (492). No review edit touches P19–P24.
- **Comments:** Andrew resolved three threads:
  - "methanogen" (Word id 7);
  - "Most abundant genes…" (id 9);
  - "what was the evidence for betaine…", with Clifton's and Susannah's replies (ids 10, 11, 84). Its anchor sentence was deleted.

  The open "restoration efforts" thread (ids 12, 13) is re-anchored on its unchanged text. The mirror has 74 comments, down from 79.
- **New files:**
  - `clauded_insights/gdoc_sync/gdoc_{sync,views}_2026-11-20.json`, built by `gdoc_sync/tools/make_spec_1120.py` from the 11-17 spec plus the API read;
  - `Metabolic Manuscript - Google Doc mirror (2026-11-20).docx`, built with `make_gdoc_mirror_docx.py <merged.docx> --date 2026-11-20`.
- **Validation:**
  - All 492 paragraphs match the Doc with suggestions rejected and with them accepted.
  - `check_comments.py`: 61 comments match the export anchors exactly. The other 13 sit on images, equations or cross-paragraph ranges, the same as on 11-17. All 57 open anchored threads are present.
  - The 11-17 mirror still rebuilds byte-identically.

**The final MAG set becomes the primary set** (`correlational/final585/`, README there). The analyses of `correlational/r232/` were rerun with the same statistics on the 585 species of the re-assembled set:
- **Abundance:** CoverM read counts (`reassembly_2026/abundance_final/`).
- **Guilds:** KOfam; a species is in a guild when any of its MAGs carries the marker.
- **Labels:** a species keeps the R232 label of its original MAG (260 species).

The review now takes these from the final set:
- **Tables:** 3, 4 (rebuilt: species and MAG counts, members, notes on originals and RAST calls), 5 (guild rows: MPn, betaine, methanogen, TMA candidates, archaea) and S6.
- **Figures:** 8, 9, S8 and S11 (`*_final585.png`).
- **Methods:**
  - the correlation and co-occurrence paragraphs (species instead of MAGs);
  - a new abundance and labelling paragraph in the re-assembly Methods (CoverM v0.7.0, minimap2 v2.30, ≥ 95 % identity, ≥ 75 % aligned);
  - CoverM and minimap2 added to the citation comment.
- **Results:**
  - the guild–methane paragraph, the methanogen paragraph and the co-occurrence paragraph;
  - the re-assembly Results moved to directly after the MAG-quality Results, ahead of the analyses that use the set;
  - the originals paragraph now cites Table S5 and says the analyses below use the final set.
- **Discussion:**
  - the *Methanolobus* sentence (species, four MAGs);
  - *Nitrosopumilus*: two species, four MAGs, reference wetland only, about 0.07 % of reads each. No MAG has a KOfam hit to *mpnS* (K18049, E ≤ 10^-5^; KOfam has no threshold for it);
  - the P95 and P99 guild sentences;
  - the P96 unbinned shares (89 / 85 / 60 / 51 % for SF2 / R2A / R2 / R1);
  - the data-availability folders.
- **Comments:** the Abstract comment, the overview comment and the Table 3, Figure 8 and S11 comments updated. A new interim comment on P35 says the models are being rebuilt.
- **Table S8:** each MAG gets its species label; representatives get mean abundance by site type and ρ/q with CH₄. The caption says so.

**Claims that changed (not just numbers):**
- **Guild specificity is gone.** No guild beats the random-species null on either basis:
  - areal: 27–48 % of random sets ≥ ρ;
  - DNA-normalised: MPn ρ = 0.78 (9 %) and betaine ρ = 0.76 (14 %).

  The earlier "MPn consumers stood out on the DNA-normalized basis (higher than 99 % of random 22-MAG sets)" is withdrawn in the P78, P95 and P99 text, the S8 caption and the Abstract comment.
- **New nominal result within the six unrestored cores:** the DNA-normalised MPn sum is nominally correlated (ρ = 0.83, p = 0.042, n = 6; 19 % of random sets do as well). The text reports it as nominal.
- **Archaea in the network:** eight of the ten archaeal taxa in the network fall in one module, which also holds eight bacterial taxa (16 in all). In the R232 network, all nine archaeal taxa formed one purely haloarchaeal module.
- **Permutation null:** shuffled residuals give denser networks (236–617 edges against 194), not networks "of similar size". Modularity is within the null, and only the same-phylum excess (21 % against at most 13 %) exceeds it.
- **Methanogen:** one species (four MAGs), still detected only in the six unrestored cores.

**Safeguards:** the builder asserts the conditions this wording depends on. These are the null fractions, the p values within the unrestored cores, *Methanolobus* detection, the no-edge guild taxa, the prefilter result and *Nitrosopumilus* detection. A refresh that changes them fails instead of writing stale text.

**Interim (models):** the metabolic models are being rebuilt for the final set's 57 MPn, betaine and methanogen species. Until then, the modelling Methods and Results (P35–P40, P69–P77), Figure 7, and the modelling sentences of the Abstract, Discussion (P95) and Conclusions (P99) still describe the original models.

**Build:** `make_r232_review_docx.py "../Metabolic Manuscript - Google Doc mirror (2026-11-20).docx"`. Validation: "All validations PASSED" (492 → 1,075 paragraphs).

**To do when the models are rebuilt:**
- Rewrite the modelling edits in `r232_review_edits_2026-11-19.json` (P18/19, P36–P40, P69–P77, P95, P99, P100, P129) for the new models.
- Remove the dangling cross-references to the old sets: P36 and P74 cite "the nine MAGs carrying *grdH*; Table 4", but Table 4 now lists the 18 betaine species.
- Replace Figure 7 and the fitted-line equations (`tracked_math_swap` at P70/P72/P75).
- Remove the interim P35 comment.
- Check the CoverM denominator for every sample: for R1_A_D1 (2026-11-20), the CoverM total of 57,655,962 reads equals `qc_read_totals.tsv`.

## 2026-11-19 (late): the re-assembly section describes the final MAG set

At Andrew's request, the re-assembly section now covers the **final MAG set**: the 935 second-pass MAGs plus the 520 MAGs of the four first-pass site co-assemblies that pass the same CheckM2 + GUNC criteria. All 718 co-assembly bins were re-scored with the rollout's QC (`reassembly_2026/annotation_runs/rollout/scripts/coasm_qc.sh`; `coasm_final/`). Every number in the text comes from these generated files:
- `reassembly_2026/comparison_rollout/summary.json`, written by `compare_rollout.py` (four sets, joint skani clustering);
- `final_guilds_summary.json`, written by `final_set_guilds.py`.

**KOfam calls for all 1,455 MAGs** (2026-10-26 release, DRAM rule, Sec rescue):
- 935 from `annotation_runs/kofam_rollout`;
- 482 co-assembly MAGs from the first-pass `kofam` track;
- 38 co-assembly MAGs from `annotation_runs/kofam_coasm38`. Its validation found 0 differing calls.

**Changes in `add_reassembly_2026()`:**
- **Methods, second pass:** one added sentence. The 718 co-assembly bins were assessed and classified in the same way, and the 520 that met the criteria form, with the second-pass MAGs, the final MAG set (Table S8).
- **Methods, annotation:** genes of the 1,324 first-pass, 935 second-pass and 38 co-assembly-only MAGs were annotated. DRAM is stated as run on the first-pass species representatives (558), which is where it was run.
- **Results, new paragraph:** 1,455 MAGs in 585 species clusters, 275 with an HQ representative and 57 archaeal.
  - 212 species come only from the co-assemblies and 148 only from the second pass.
  - The set contains 504 of the 557 first-pass species.
  - It recovers 268 of the 310 originals (any sample), including 16 of the lost 30.
- **Results, guild paragraph:** now the final set. If `final_guilds_summary.json` is missing, the script falls back to the first-pass text.

  | Guild | MAGs | Species |
  |---|---|---|
  | *phnJ* (MePn consumers) | 109 | 39: 19 Rhodobacteraceae, 6 haloarchaea |
  | betaine reductase | 53 | 18: 17 Bacillota, 1 Cloacimonadota |
  | *mcrA* | 4 *Methanolobus* | 1 |
  | TMA producer | 0 | 0 |

  - Sec read-through takes *grdH* genes from 10 to 58.
  - One-sided TMA genes: *cutC* without *cutD* in 13 MAGs, *cntA* without *cntB* in 70.
  - MttB family without *mcrA*: 355 MAGs, including 138 of the 139 Rhodobacteraceae.
  - RAST MePn consumers: 22 MAGs in 21 species; 15 of those species are recovered with *phnJ*, and 20 further species have no original MAG.
  - The *Albimonas*.1 and BM714.4 statements still hold: the 6 *Albimonas* MAGs have no *phnJ*, and the BM714.4 cluster has no betaine call.
- **Table S8:** a new caption after Table S7, for a separate file like Table S5: `Table S8 - re-assembled MAG set taxonomy, quality and guilds (GTDB R232).xlsx`, built by `reassembly_2026/make_table_s8.py`. It has one row per MAG (1,455) with source, sample or co-assembly, site, GTDB R232, CheckM2/CheckM/GUNC, quality, size/contigs/N50, species cluster and its representative, overlap with an original species, and guild calls.
- **Comment on the rollout paragraph:** now names these sources. It notes that the main analyses (Tables 3–5, Figures 7–8, metabolic models) still use the 310 original MAGs.

**Refresh order:**
1. `compare_rollout.py`
2. `final_set_guilds.py`
3. `make_table_s8.py`
4. `make_r232_review_docx.py "../Metabolic Manuscript - Google Doc mirror (2026-11-17).docx"`, then validate.

Validation: "All validations PASSED" (492 → 1,081 paragraphs). The previous build is in the session scratchpad as `review_before_finalset.docx`.

## 2026-11-19: rebased on the 2026-11-17 Google Doc; re-assembly section made final

**Rebase.** `make_r232_review_docx.py "../Metabolic Manuscript - Google Doc mirror (2026-11-17).docx"` now builds the review on the 2026-11-17 mirror.
- **Index shift:** `main()` wraps `doc.paragraphs` in `Paras`. It keeps July-export indices and adds 1 from P19 on when the input has the split Abstract (asserted: 282 paragraphs, and P19 starts "Here we explore"). On the 11-03 mirror it changes nothing; that build still passes validation.
- **Edits file:** with the 11-17 mirror the script reads `clauded_insights/r232_review_edits_2026-11-19.json`, a copy of the 10-26 file with two changes (the 10-26 file is untouched):
  - **Abstract modeling replace, retargeted** (`mirror_index` 19) to Andrew's rewritten sentence, which still claims full reproduction: "This was corroborated by metabolic modeling, where trace consumption of methylphosphonic acid was sufficient to reproduce the experimental methane emissions." It now reads "Metabolic modeling partly corroborated this: … could supply about one fifth to three quarters of the measured methane emissions, depending on the assumed share of methylphosphonate in phosphate uptake, although the simulated emissions did not track the differences among cores."
  - **New comment** on "a primary methanogenic substrate", the claim the Doc moved onto the in vitro sentence: the modeling supports methylphosphonate as a source, not as the primary substrate.
- **Abstract comment:** "The modeling clause is now qualified" follows the clause to P19. Its wording now matches the repaired-model edit.
- **Check:** apart from the Abstract and the five Introduction paragraphs Andrew edited (P20–P24), every body element of the new build equals the 11-03 build's element one position earlier, including all tracked changes. Validation: "All validations PASSED" (492 → 1,079 paragraphs).
- **Previous build:** kept in the session scratchpad as `review_before_20261119.docx`.

**Re-assembly section made final.** The rollout finished 2026-11-19 12:27 UTC with 39/39 samples. The section is built by `add_reassembly_2026()`; its tracked insertions are now dated 2026-11-19.
- **Methods:**
  - The second-pass paragraph now states the outcome of the two-sample comparison: R1_A_D1 met every criterion and SF2_C_D1 exceeded the contamination criterion, so the second pass was applied to all 39 samples with the CheckM2 + GUNC filter.
  - It also states the GTDB-Tk R232 classification and the joint three-set species clustering.
- **Results:**
  - **Rollout paragraph:** switches to final wording automatically once `n_complete == n_samples`. Kept MAGs (HQ; median N50): 935 (473; 18.0 kb), against 603 (145; 5.4 kb) for the first pass and 345 (174; 17.2 kb) for the originals. 334 of the 450 originals are recovered in the same sample, including 14 of the 37 species-lost originals.
  - **New paragraph** from `reassembly_2026/comparison_rollout/summary.json` (made by `compare_rollout.py`):
    - 373 species, 45 of them archaeal;
    - set-wide recovery of 263/310 originals, against 256 for the first pass, including 16 of the 30 lost;
    - best-MAG quality and N50 against the originals and against the first pass;
    - the first pass's extra species: 212 of its 230 species not shared with the rollout come only from co-assemblies.
- **Comment:** a "final numbers" comment replaces the PRELIMINARY one. It names the sources and lists what is pending:
  - KOfam/guild calls for the second-pass MAGs (running), so the guild paragraph still describes the first-pass MAGs;
  - the co-assembly MAGs, which are being re-scored with the rollout's QC at Andrew's request.
- **Table S7:** the caption now covers all 39 samples.
- **Refresh:**
  1. `bash /scratch1/afreiburger/salterns_2026/vs2020/skani_completed.sh` (poplar).
  2. `manuscript_numbers.py`.
  3. `compare_rollout.py` (needs `comparison_rollout/three_way_skani_edges.tsv` from poplar).
  4. `make_r232_review_docx.py "../Metabolic Manuscript - Google Doc mirror (2026-11-17).docx"`, then validate with `--original` set to the same mirror.

## 2026-11-16 addition: re-assembly, re-binning and annotation of the metagenomes (preliminary)

Tracked insertions dated 2026-11-16, built by `add_reassembly_2026()` in `make_r232_review_docx.py`:
- **Methods:** a new subsection, "Re-assembly, re-binning, and functional annotation of the metagenomes", placed after the MAG-quality methods. It has three paragraphs:
  - the first pass (KBase: Trimmomatic → metaSPAdes, site co-assembly with MEGAHIT, three binners with site multi-sample coverage, DAS Tool 1.1.7, CheckM, GTDB-Tk 2.8.0/R232);
  - the second pass (poplar: BBDuk → BBMerge with extension → Tadpole → standard-mode SPAdes 3.15.3, the same binners and DAS Tool, then CheckM2 + GUNC keep filter, with the two-sample pre-registered comparison);
  - annotation (Prodigal + KOfam 2026-10-26 with the DRAM rule, the selenocysteine read-through for Grd genes, DRAM 1.5 on the 558 best-per-species MAGs, and the revised guild rules).
  - A comment lists the citations to add.
- **Results:** a new subsection, "Re-assembly and re-binning of the saltern metagenomes", placed after Co-occurrence. Its four paragraphs cover:
  - the first-pass outcome (1,324 MAGs, 558 species vs 296, fragmentation, and the reconstruction of the original assemblies);
  - the pilot and binning test (R1_A_D1 and SF2_C_D1);
  - the second-pass partial rollout, with numbers taken from the latest `reassembly_2026/manuscript_numbers_<date>.json`. A "PRELIMINARY" comment gives the snapshot date and how to regenerate it;
  - the KOfam/DRAM guild calls for the first-pass MAGs.
- **Table S7:** per completed sample, the second pass vs the first pass vs the original MAGs (kept, HQ, N50), plus the originals recovered. All three sets are scored with the same CheckM2 1.1.0 + GUNC 1.1.1.
- **To refresh as the rollout progresses:**
  1. On poplar, run `bash /scratch1/afreiburger/salterns_2026/vs2020/skani_completed.sh`.
  2. Run `~/Documents/py_venv/bin/python reassembly_2026/manuscript_numbers.py`.
  3. Run `make_r232_review_docx.py <mirror>`, then validate.
- **Sources:**
  - `reassembly_2026/comparison/README.md`
  - `binning_experiments/assembly_vs_jgi/README.md`
  - `annotation_runs/{tadpole,binning_test,kofam,dram,rollout}/` (RESULTS.md / REPORT.md)

## 2026-11-08 revision: the three MePn-obligate models repaired

**What changed.** Three of the 20 simulated MePn-consumer models (*Roseovarius*.9, *Roseovarius*.11 and *Albimonas*.1) could not grow without methylphosphonate. That dependence is a modelling artefact (`clauded_insights/mepn_obligacy_check/README.md`, checked 2026-11-04/05 and confirmed independently on 2026-11-08). Their stored methane per unit of biomass was 11–13× that of every other model at ratio 0.07 (3–3.4× at 0.33), and they produced 72 % of the simulated methane at ratio 0.07.
- **Causes:**
  - *Roseovarius*.9 and *Albimonas*.1 have no outlet for surplus one-carbon units except MePn degradation, because homoserine kinase is unannotated. It is also missing from the other Rhodobacteraceae MAGs and from reference genomes that grow without MePn.
  - *Roseovarius*.11 lacks PRPP synthetase. Its MAG is 81 % complete and 20 % contaminated.
  - The models were gapfilled on a phosphate-free MePn medium, which hid these gaps.
  - *Albimonas*.1's C–P lyase genes are frameshifted.
- **Repair:** one reaction each: homoserine kinase (rxn01300) for *Roseovarius*.9 and *Albimonas*.1, PRPP synthetase (rxn00770) for *Roseovarius*.11.
  - Afterwards every model gives r/(1+r) × 0.660485 methane per unit biomass, and the ratio constraint binds in all 20 models.
  - All 29 metabolic repair variants give the same methane.

**Corrected numbers** (consistent areal basis; `modelling_consistent_basis_r232.json`, sets `MePn_*_20_repaired`, asserted against the cell-27 re-run in `mepn_obligacy_check/analysis/step5_consistent_basis.json`):

| | as built (2026-10-26 text) | repaired (now in the text) |
|---|---|---|
| Model / measured at 0.07 | 56 % (18–141 % by core) | **19 % (6–46 %)** |
| Model / measured at 0.33 | 100 % (30–245 %) | **72 % (22–174 %)** |
| Spearman ρ (exact p), R² at 0.07 | −0.20 (0.71), 0.08 | −0.09 (0.92), 0.05 |
| Fit P70 / P72 | −0.20x + 907 / −0.32x + 1563 | **−0.06x + 292 / −0.21x + 1106** |
| Share of the three models | 72 % | 18 % |
| Rhodobacteraceae share | 99.6 % (*Roseovarius* 81.5 %) | 98.9 % (*Roseovarius* 67.7 %) |

**Edits revised** (dated 2026-11-08):
- **Abstract P18:** "about half" becomes "part (about one fifth to three quarters, depending on the assumed share of MePn in phosphate uptake)", and "important source" becomes "a source".
- **Methods P36:**
  - new sentences describe the obligacy, its cause and the one-reaction repairs;
  - the ratio-constraint rationale now says the constraint binds in all 20 models;
  - a new comment lists the caveats: local template, no re-gapfilling, frameshifts not read-checked, and export repairs (at most 31 %).
- **Results P69–P73:**
  - the new totals, ρ and R²;
  - the equations P70/P72;
  - P73 now says the "non-binding" three were an artefact, and the earlier "did not bind … natural outcome" wording is withdrawn;
  - the Rhodobacteraceae share.
- **Discussion P95:** "about one fifth (three quarters at 0.33)". The last sentence now says the models translate MPn-consumer biomass and the assumed ratio into methane rather than testing their metabolism.
- **Conclusions P99:** 19 % / 72 %. The comment on P100 is updated to match.
- **Figure 7 and caption:** the repaired set is filled orange with its fit. Open orange circles show the as-built models, with no fit.

**Still open:**
- **No re-gapfilling:** the drafts were not re-gapfilled from scratch on a Pi-containing, MePn-free medium. That would be the cleanest test.
- **Template:** fliu's KBase template (ws 12218) was not available locally.
- **Ratio:** the magnitude now depends only on the consumer biomass and the assumed ratio r. Nothing in the data fixes r.

## 2026-10-26 revision

What this revision changed, and where the evidence is:

- **Modelling text and numbers (P18, P36–P40, P69–P77, P95, P99, Figure 7 caption; comment on P100).** The headline finding of this revision: **on a dimensionally consistent basis, the simulated methane does not track the measured emissions among the six unrestored cores.**
  - **Why the old fits looked good.** The stored simulations fixed each model's biomass to relative abundance ÷ DNA, averaged over depths. Two depths have almost no DNA (R2_B_D1 and R2_C_D2, 0.16 µg g⁻¹), so every simulated value carries a 1/DNA term. All earlier fits came from that term:
    - the draft's R² = 0.978 / 0.968 / 0.862, with both axes divided by DNA (0.966 matches no stored fit);
    - an intermediate "areal" correction from July (`figure7_corrected.py`, ρ = 0.60), where multiplying by the core's DNA did not cancel the term.

    The factor VW DNA × mean(1/DNA) alone gives ρ = 0.60 against measured methane.
  - **Consistent basis.** Each MAG's flux per unit of simulated abundance is rescaled to its biomass in the core: relative abundance × DNA × layer thickness, summed over D1 and D2.
  - **MePn models:**
    - ρ = −0.20, exact p = 0.71, n = 6, R² = 0.08 at ratio 0.07, with 56 % of measured emissions in total (18–141 % by core);
    - ρ = −0.09, R² = 0.06 at 0.33, with 100 % in total (30–245 %);
    - so the models can supply methane of the measured magnitude but not the pattern among cores.
  - **Betaine models** (TMA mostly at the uptake bound, so converted as before): ρ = −0.09, R² = 0.10.
  - **Other numbers:**
    - CH₄ : MePn is exactly 1 : 1, not 1.4 : 1.
    - Figure S7 (DNA vs CH₄) is R² = 0.22, ρ = 0.52, p = 0.08, so "significantly" became "weakly and not significantly".
    - The "unconstrained" R² = 0.970 run is not stored. It matches the three MAGs (*Roseovarius*.9, .11 and *Albimonas*.1) for which the ratio constraint does not bind; they supply 72 % of the simulated methane.
  - **Sources:** `clauded_insights/modelling_consistent_basis_r232.py` → `.json`. `modelling_numbers_r232.json` traces the draft's values but its "areal" basis is superseded.
- **Equations P70, P72, P75 (Office Math).** Each is replaced as a tracked deletion of the old equation plus a tracked insertion of the new fit:
  - P70: −0.20x + 907;
  - P72: −0.32x + 1563;
  - P75: −0.71x + 3192.

  The units drop the 1/(µg DNA g⁻¹ soil) term. Word shows these as tracked changes; LibreOffice ignores revisions on equations and shows both versions.
- **Betaine models.** These are restricted to the nine *grdH* carriers. A re-simulation reproduces the stored fluxes exactly, because no model fix was supported by the annotations. The fit is ρ = −0.09, p = 0.87, R² = 0.10. Other findings:
  - The betaine reductase reaction (rxn17220_c0) was added by hand to all 15 models with no gene association.
  - 33 of 54 fluxes sit at the 0.1 uptake bound, so modelled betaine methane tracks DNA density (R² = 0.97), not abundance.
  - Only *Halanaerobium*.2 is forced to zero TMA. The old review comment that "the *Halanaerobium* MAGs make no TMA" was wrong for *Halanaerobium*.1 and .3.
  - Evidence: `correlational/betaine_grdH_r232/`.
- **Figure 7 replaced** by `clauded_insights/figure7_r232.png` (`figure7_r232.py`). It is a single panel on the consistent areal basis with a 1:1 line, MePn (20 models) and betaine (9 *grdH* models). The old panel b is superseded by Figure 9, Figure S8 and Table 5.
- **Other corrections from the verification pass:**
  - Methods cited Figure S9 for gap-filling; it now cites Figure S10.
  - A Discussion sentence attributing the community-wide correlation to "higher metabolic activity" is deleted; the revised Results attribute it to the site gradient.
  - Table 4's last column is now "In the curated modeling set but lacking the marker", noting which of those MAGs were simulated.
  - The new co-occurrence wording is tightened (filters, edge counts, maxT).
  - The deleted and new Table 3 are separated by a paragraph so that Word does not merge them.
  - A comment on the Conclusions suggests hedging "primary source … precise mechanism".
- **Halomonas isolates (P47, P67, Figure 5 caption).** GTDB-Tk 2.8.0 / R232 keeps both isolates in *Halomonas*.
  - I1 is *Halomonas* sp. D, the C–P lyase carrier. It is deposited as GCA_032334135.1 and is its own species cluster, *Halomonas* sp032334135.
  - I4 is *H. saccharevitans*: 98.5 % ANI to the type strain.
  - R232 places 11 of the 23 Figure 5 references labelled *Halomonas* in other genera. Three of the six MPn-positive references remain *Halomonas*.
  - Evidence: `MAG_QC/gtdbtk_isolates_r232/` (`run_gtdbtk_isolates.sh`, README).
- **Co-occurrence paragraph and Figure 8 caption (R232 insertions) qualified** by `correlational/r232/prefilter_sensitivity/within_REPORT.md`:
  - The network's size and modularity fall within Graphical Lasso fits to permuted residuals; only the excess of edges between related taxa does not.
  - The key taxa have no edges because the lasso acts on unstandardized residuals (degree tracks residual variance, ρ = 0.90).
  - Relaxing the prevalence filter gives no edge to a guild taxon that lacked one; new edges mostly join rarely detected taxa.
  - At most two pairs reach the maxT threshold, depending on the number of permutations.
- **Smaller corrections:** "Saligentibacter" → *Salegentibacter* (P67).

## What the new analysis changed

- **Quality (CheckM):** the CheckM re-run reproduced the original values exactly. 275 of 310 MAGs meet MIMAG medium quality; 35 exceed 10 % contamination, including *Roseovarius*.11 and .12.
- **Taxonomy (GTDB-Tk 2.8.0, R232):**
  - 280 bacterial and 30 archaeal MAGs.
  - 290 MAGs resolved to genus (originally 201) and 176 to species (originally 19).
  - The genera central to the MPn argument are unchanged. One *Halomonas* MAG is now *Vreelandella*.
- **Gene-based guilds:**
  - MPn consumers (*phnJ*): 22, of the 27 in the modeling set.
  - Betaine reducers (*grdH*): 9, of the 15 in the modeling set.
  - Methanogens (*mcrA*): 1, *Methanolobus* sp020055655.
  - *cutC* is not annotated in any MAG.
  - The before/after membership is in `correlational/r232/guild_membership_2020_vs_r232.csv` / `.md`.
- **Guild abundance vs methane:**
  - On the areal basis, every guild correlates across the 12 cores, but none beats random MAG sets.
  - On the DNA-normalized basis, only MPn consumers stand out (ρ = 0.92, above 99 % of random sets).
  - Nothing is significant within the six unrestored cores.
- **Methanogen:** *Methanolobus*.1 occurs only in the unrestored cores, and within them it does not track methane (ρ = −0.54, p = 0.27).
- **Co-occurrence:** the archaea now enter the network (169 taxa, 13 archaeal). No pair passes FDR. The network has 6 Louvain modules (Q = 0.747), including one haloarchaeal module. Its size and modularity are within the permutation null, and the key taxa have no edges, for the variance reason above.

## Tracked text changes (P-numbers refer to the original paragraphs)

| Where | Change |
|---|---|
| Abstract (P18) | Modelling clause: "reproduced the rank order, though not the magnitude" |
| Methods, after Metagenomics (P30) | **New subsection** "MAG quality, taxonomy, and functional guilds": CheckM, GTDB-Tk R232, iterativeID naming, gene rules for the guilds |
| Methods, Metabolic Modeling (P35) | "dereplicated medium and high quality MAGs" → "310 dereplicated MAGs (275 meet MIMAG medium quality; Table S5)" |
| Methods (P36) | Phosphate ratio = a minimum MePn : Pi uptake ratio; 20 of 27 MePn consumers simulated (why); betaine models = the nine *grdH* carriers; 0.1 = maximum uptake (mmol g~DW~⁻¹ h⁻¹, unit to confirm) |
| Methods (P38, P40) | Regression and Spearman; conversion with each core's DNA density (areal basis) |
| Methods, isolates (P47) | GTDB-Tk v2.8.0 / R232 reclassification sentence; BioProject and assembly accessions |
| Methods, Analyses (P52) | Replaced: correlation design (12 cores and 6 unrestored, two bases, 4,000-set random comparison, BH-FDR) |
| Methods, Analyses (after P53) | **New paragraph:** co-occurrence methods (CLR, within-site residuals, permutation FDR, graphical lasso, Louvain) |
| Results, after Metagenomics (P57) | **New subsection** "MAG quality, taxonomy, and functional guilds" |
| Results, culturing (P67) | R232 placement of both isolates and of the six MPn-positive references; *Salegentibacter* spelling |
| Results, modeling (P69–P77) | Areal-basis statistics replace every old R² (see above); equations P70/P72/P75 swapped; P76 rewritten for the nine *grdH* betaine models; P77 stoichiometry 1 : 1 |
| Results (before P78) | **New heading** "Guild abundances and methane flux" |
| Results (P78) | Replaced: the old R² values (0.974/0.966/0.852/0.798) give way to the R232 guild correlations and random-set comparison |
| Results (after P78) | **New paragraphs:** methanogen vs methane; new heading "Co-occurrence" with the network results and their permutation/variance caveats |
| Discussion (P83) | Inserted sentence: the only *mcrA* MAG is *Methanolobus* sp020055655; it is restricted to the unrestored salterns and does not track methane |
| Discussion (P86) | Inserted sentence: the *Nitrosopumilus* MAGs occur only in the reference wetland, and RAST annotated no *mpnS* in any MAG |
| Discussion (P95) | First sentence replaced (MPn specificity only on the DNA-normalized basis); modelling sentences replaced (rank order, 1 : 1, betaine models not correlated; models translate abundance into methane) |
| Conclusions (P99) | "Secondly…" and "Thirdly…" sentences replaced to match |
| Data Availability (P102) | Added the repository location of the MAG QC and correlation analyses |

## Figures and tables

| Item | Status |
|---|---|
| Table 3 | **Replaced.** Recomputed with R232 families for all 310 MAGs (archaea included; † marks archaeal families) and correct unbinned shares (52 / 88 / 89 %). **Cause of the old table's errors:** it came from `correlational/abundances.ipynb` cell 1, which had three faults:<br>• `'R1' in col or 'R2' in col` also matches the six R2A cores, so "unrestored" pooled R1 + R2 + R2A.<br>• The 30 archaea had no taxonomy entry and were filled as `unbinned_reads`, so they were dropped with the unbinned reads; "68 %" is unbinned plus archaea over those 18 cores.<br>• The notebook outputs unrestored / reference / restored, which was copied under the headers Unrestored / Restored / Reference.<br>The notebook is now fixed (2026-10-21) and prints the corrected 2020-family values. They are also in `Table 3 - corrected values (2020 GTDB families).csv`. |
| Table 4 | **New.** Functional guilds: marker gene, members, and what was dropped from the modeling sets. |
| Table 5 | **New.** Guild–methane Spearman correlations on both bases, 12 cores and 6 unrestored, with the random-set fractions. |
| Figure 5 | Caption: R232 genus placement of the reference genomes labelled *Halomonas* (comment on the one inferred genome). |
| Figure 7 | **Replaced** (2026-10-26) by a single panel on the areal basis (`clauded_insights/figure7_r232.png`); caption rewritten. |
| Figure 8 | **Replaced.** Both old panels and the caption are deleted. Panel b (the Louvain view) was re-rendered with R232 labels and archaea and is now the single panel; panel a (the plain graphical-lasso layout of the same edges) was not re-rendered. |
| Figure 9 | **New.** Guild abundance vs methane, areal basis (`correlational/r232/guild_methane_raw_r232.png`). |
| Figure S8 | **Replaced** with the DNA-normalized version of Figure 9; the stray caption fragment (old P151) is deleted. |
| Figure S11 | **Replaced** (2026-10-30) by `correlational/r232/figureS11_r232.png` (`figureS11_r232.py`): 169 GTDB R232 taxa (13 archaeal, †), 186 graphical-lasso asterisks, gene-based highlights (*phnJ* orange, *grdH* blue, *mcrA* purple). Caption rewritten (the old one misdescribed the asterisks as ρ > 0.7); now cited in the Co-occurrence results. |
| Table S5 | **New caption.** Content is in the Excel file. |
| Table S6 | **New.** Per-core methane flux, *Methanolobus*.1 abundance, and total archaeal abundance. |

## Comments (18, by "Claude (R232 review)")

1. **Start of document:** summary of this review and revision.
2. **Abstract:** the modelling clause is now qualified; consider qualifying any guild-correlation claim too.
3. **Methods P47:** KBase narrative 139158 is not readable by the author's account and has no public read entry, so it does not appear to be public.
4. **Methods P47:** reference 68 is the GTDB-Tk v1 paper; add the GTDB-Tk v2 paper.
5. **Results P67:** "*Halomonas* sp. D" is isolate I1 (NCBI "Halomonas sp. I1"); use one name.
6. **Results P67:** the C–P lyase *Halomonas* is supported by the isolate, not by the MAGs (neither Halomonadaceae MAG carries *phn* genes).
7. **Results P69, P71:** what the replaced equations were and why the old R² (and the intermediate ρ = 0.60) were inflated by the 1/DNA term.
8. **Results P74:** the replaced betaine equation and the betaine-model caveats (hand-added reaction, O₂ dependence, non-unique TMA optima, contamination).
9. **Discussion P83:** reconcile *Ca.* Methanosalis SBSPR1A (ref. 58) with GTDB's *Methanolobus* sp020055655.
10. **Discussion P93:** the 2021 *Halomonas* species count predates the split into *Halomonas* and *Vreelandella*.
11. **Table 3:** the three confirmed causes of the old table's errors.
12. **Figure 5 caption:** the R232 counts assume the *H. campaniensis* tip is strain SK03.
13. **Figure 8 caption:** the figure is now a single panel; the old panel a (plain glasso layout) was not re-rendered.
14. **Figure S11 caption:** what was regenerated and why the old asterisk description was wrong.
15. **Table S3:** not in the document; consider adding the isolates' genome-based names and accessions.
16. **Conclusions:** consider hedging "primary source … precise mechanism" in light of Figure 7 and Table 5.
17. **Methods P36 (2026-11-08):** the repair of the three MePn-obligate models, its limits, and the as-built numbers.

(Topic 7 carries two comments, which makes 18 in total.) Resolved and removed in this revision: the isolate-reclassification request, the modelling-R² comment, the betaine-model comment, and the Figure 7 panel-b comment.

## Not changed

- **FBA simulations of the MPn models** were not re-run, except for the three repaired models (cell 27 re-run in `mepn_obligacy_check/`; the other 17 are unchanged by the repair). They use the depth-mean DNA-normalised abundances of commit 22f4629, not the current DNA-weighted merge; since each model's methane is proportional to its abundance, a re-run would only rescale them.
- **Figure 5:** tip labels are the original KBase names (the tree object in narrative 139158 could not be read).
- **Citation formatting:** citations in inserted text, e.g. (51), (53), (67), (68), are plain text, not reference-manager links. Relink them if you use one.
