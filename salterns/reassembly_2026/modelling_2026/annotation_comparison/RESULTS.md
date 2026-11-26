# Annotation source for the 2026 salterns models: results of the pre-registered comparison

Run 2026-11-20/21 (TF), analysed 2026-11-22. Rule and metrics: `../ANNOTATION_CHOICE_CRITERIA.md`, fixed before any compared model was built. The arithmetic used to apply the rule was written into `../pipeline/compare_arms.py` (docstring) after 7 test builds and before the full run was inspected. Nothing was written to KBase and nothing was committed.

## Decision: arm A (RAST)

No merged arm meets all five criteria, so the rule returns arm A.

- **Arm B (RAST + KOfam)** meets (a)–(d) but fails (e): the median blocked-reaction fraction rises by 2.69 percentage points (limit 2).
- **Arms C (+ eggNOG) and D (+ Bakta)** fail (d) and (e).
  - (d): one model, *RI-231* `R1_sediment__metabat2.bin.142_sub`, makes ATP from NH3 + O2 in C and D but not in A.
  - (e): the blocked fraction rises by 6.5 and 6.4 points.
- All three merged arms clear (a) and (b) by wide margins. C and D give 72 % more gene-associated reactions and 36 % fewer gap-filled reactions than A. The added evidence enlarges the models: about 600 extra reactions per model in C and D carry only non-RAST evidence. Much of that addition does not connect, so the blocked fraction rises. About 250–290 of those reactions per model (median) rest on a single weak source (EC-only, or eggNOG without a KO; metric 6).

## What was compared

**Genomes:** the 57 species representatives of `../../comparison_rollout/modelling_species.tsv`.
- *UBA8670* `R2_sediment__metabat2.bin.434_sub` is both a MePn consumer and a betaine reducer, so it is built under both protocols. That gives **58 models per arm**, all paired.
- Guilds: MePn consumers 39, betaine reducers 18, methanogen 1 (*Methanolobus*).

**Genes:** one gene set for every arm, the Prodigal 2.6.3 genes (`-p meta`) of the KOfam tracks.
- RAST (service; `../rast_service/`) and KOfam annotate these genes.
- eggNOG-mapper 2.1.13 ran on the same proteins, so its annotations join by protein id.
- Bakta 1.12.0 calls its own genes. They were attached to Prodigal genes at ≥ 80 % reciprocal coordinate overlap, same contig and strand.
  - 161,756 of 169,139 non-pseudo Bakta CDS (95.6 %) matched.
  - The other 7,383 entered arm D as extra features `bakta_<locus_tag>`.
- Model gene ids are the Prodigal protein ids without the `<genome>|` prefix, because cobra parses `|` in a GPR as OR.

**Evidence per arm** (`../pipeline/evidence.py` → `evidence/<mag>.json.gz`; per-genome counts in `evidence_stats.tsv`):

| Arm | Gene → reaction evidence added to RAST |
|---|---|
| A | none: RAST roles through the templates' complexes |
| B | KOfam best KO (DRAM rule), with the 53 accepted Sec read-through rescues (`annotations/kofam/sec_readthrough_57.tsv`) replacing the KO of their 106 piece genes |
| C | B + eggNOG `KEGG_ko`, `EC`, `BiGG_Reaction` |
| D | C + Bakta EC numbers, and Bakta products matched to template role names (`product` and `product (EC x)`, normalised as ModelSEEDpy normalises roles; triggering roles only) |

**Mapping terms to reactions** (`../pipeline/maps.py` → `../pipeline/maps/term_maps.json`):

| Source | Path to ModelSEED reactions |
|---|---|
| KO | KEGG REST `link/reaction/ko` (cached 2026-11-20 21:09 UTC as `../pipeline/kegg_cache/link_reaction_ko_2026-11-20.tsv`, md5 `864f3a82…`, 12,416 links), then the ModelSEEDDatabase `KEGG` reaction aliases |
| EC | ModelSEEDDatabase `Unique_ModelSEED_Reaction_ECs.txt`, complete 4-level ECs only |
| BiGG | ModelSEEDDatabase `BiGG` reaction aliases |

- The ModelSEEDDatabase checkout was `tautomer-adenine-n1_20261113` at `41345b04`, read only.
- Only reactions present in the template were used.
- Extra genes enter as a gene association (an extra GPR term) of the template reaction, before the spontaneous and universal reactions are chosen. The template is never changed by an arm.

**Pipeline:** identical in every arm (`../pipeline/build_model.py`; protocols in its docstring and in `../models/README.md`).

- **Templates:** the combined curated templates (`../templates/`, md5s in `../templates/PROVENANCE.md`), chosen by file name: Archaea → ArchaeaTemplateV6; Bacillota (and its splits), Actinomycetota, Chloroflexota or Patescibacteriota → GramPosModelTemplateV7; any other phylum → GramNegModelTemplateV7.
- **ATP correction:** MSATPCorrection (Core-V6, ModelSEEDpy's 54 ATP media). Its tests are passed to every MSGapfill.
- **Guild protocols:**
  - MePn: the gapfill.ipynb protocol, with a first gap-fill on the MePn-free lactate + Pi medium (scouting advice).
  - Betaine: aerobic lactate, then the cell-2 pathway additions.
  - Methanogen: methanol + methylamines, then methane export.
- **Software:**
  - ModelSEEDpy `sync-cshenry-main` at `311f720`; TF installed a `git archive` of that commit. Python sources are identical to the local editable install.
  - cobra 0.29.1, optlang 1.9.1.
  - Gurobi 13.0.1, 1 thread per build (`tf_env/gurobi.env`).
  - Environment `~/venvs/salterns_modelling` on TF (`tf_env/pip_freeze.txt`, `tf_env/COMMITS.tsv`).

**Code versions:** 210 of the 232 builds ran `build_model.py` v1 (md5 `7f8f7f65…`, kept as `tf_env/build_model_v1_main_run.py`).
- The other 22, all MePn consumers, failed in v1 at the second gap-fill. MSGapfill's prefilter raised "the model itself violates them", because the first gap-fill's integrated solution already broke an ATP test. In 15 of the 22 that test was Pyr.SO3: assimilatory sulfite reductase makes sulfite an electron sink for pyruvate fermentation.
- They were re-run with v2 (md5 `99652b0d…`), which drops the tests the model already violates for that stage only and records them (`atp_tests_violated_at_start`). v2 behaves identically to v1 whenever no test is violated, which was the case for all 210 v1 builds; otherwise their prefilter would have raised.
- The final ATP-test pass rate still counts the violated tests as failures.
- The 22 v1 error reports are in `errors_main_run_v1/`.

## Results

The four tables (metric medians, per guild, decision rule, per model) are at the end of this file. Every value is in `per_model_metrics.tsv` and `summary.json`; the raw build reports and models are in `runs/<arm>/`.

## Points needed to read the numbers

1. **Betaine guild support is 0 in arm A by construction.**
   - All three templates carry the betaine reductase roles (ftr33907, ftr01398, ftr11139), but no complex links them to a reaction.
   - The template's betaine reductase rxn11678 therefore cannot be triggered by RAST. KOfam K21579 (grdH) → R07228 → rxn11678 triggers it.
   - Criterion (c) is met by B–D for this guild, but only because arm A cannot score.
2. **The template writes betaine reductase and glycine reductase backwards.**
   - rxn11678 (betaine reductase) is written acetyl-P + TMA + oxidised thioredoxin → betaine, direction `>`. rxn11676 (glycine reductase) has the same problem.
   - Both reactions can only run in the synthesis direction, the reverse of Stickland reduction.
   - When KO evidence adds rxn11678, it forms a flux-neutral null cycle with the hand-added rxn17220_c0 of the betaine protocol. Betaine-reduction capacity is therefore measured as TMA export, not as rxn17220 flux.
   - This is a template curation item.
3. **The ATP tests differ by arm.** Each arm's own ATP correction selects them, so the test count is reported next to the pass rate (median 7 tests in A, 10 in B–D).
4. **"Free ATP" here means ATP from NH3 + O2 (or NH3 + nitrate), not from nothing.**
   - The salts medium of the free-ATP test contains NH3 as the N source.
   - Every failing model (A: 3; B: 3; C and D: 1) oxidises NH3 to nitrite with a nitrite reductase run in the oxidising direction: rxn13647 (ubiquinol) or rxn03978 (cytochrome c), plus rxn14173 to nitrate in C and D.
   - In A and B the reaction comes from the energy-limited second MePn gap-fill (lactate 1). In C and D it comes from annotation.
   - None of these organisms is a nitrifier. ModelSEEDpy's ATP-test media contain no NH3, so ATP correction cannot detect this.
   - Robustness check: the gap-fill solutions (any stage) that contain rxn13647 or rxn03978 number 6 in A, 7 in B, 4 in C and 3 in D. A restriction on them could not move the (e) medians by the 0.7–4.5 points needed.
5. **Paired view.** For arm B, the ratio of medians gives m2 = 0.773 (passes (b)), but the median per-model change is −17.8 %, which would fail. The ratio of medians decides, as fixed in advance. B fails (e) on either view, so this does not matter.
6. **The added reactions are mostly template-conditional, but in C and D about 106 per model are template `gapfilling`-type reactions.** These are reactions the template normally offers only to the gap-filler; EC numbers add them as gene-associated. That accounts for much of C and D's blocked-fraction rise.
7. **eggNOG's BiGG column adds almost nothing.** 1,565 genes carry a value in total, and the column holds BiGG gene ids rather than reactions, so it maps to nothing. KO mappings sometimes miss the template's own reaction ids:
   - K06163 (phnJ) maps to rxn40081, which is not in the templates. The template's C–P lyase is rxn26459, reached by the RAST role.
   - K21580 now maps to R08791 (corytuberine N-methyltransferase, rxn14200).
8. **Non-growing models after gap-filling:** A 2, B 1, C 4, D 3 (`summary.json`, `non_growing_per_arm`). In each case the gap-filler reported a solution, but the integrated model did not grow; see `../models/README.md` for the arm-A cases.
9. **The methanogen behaves differently by arm**; this is the main modelling caveat for the final set.
   - **Arm A:** *Methanolobus* grows 0.42 h⁻¹ on methylamines as a methyl → acetate acetogen. CH4 at maximum growth is 0, its FVA range is [0, 1] at 99 % growth, and with acetate export blocked growth is 0. The arm-A draft lacks a working route that reduces methanophenazine: the template's Fpo (rxn24614) is `<` only, and the ATP correction filtered the `>` direction of Rnf (rxa45615). The stage-2 methane gap-fill only made CH4 possible apart from growth (formate hydrogenlyase, rxn08518).
   - **Arms B–D:** the KO-mapped methanol:CoM methyltransferase (rxn10568) lets the model grow methanogenically on methanol: CH4 52–64 at maximum growth, and growth with acetate export blocked equals full growth. These models do not grow on any methylamine alone.

   | Arm | Growth (h⁻¹) | CH4 at max growth | Growth, acetate export blocked | Grows on methanol only | Grows on MMA / DMA / TMA only | Growth on the manuscript's DMA + H2 |
   |---|---:|---:|---:|---|---|---:|
   | A | 0.418 | 0 | 0 | no | yes / yes / yes | 1.32 |
   | B | 0.158 | 64.2 | 0.158 | yes | no / no / no | 0.81 |
   | C | 0.254 | 52.3 | 0.186 | yes | no / no / no | 0.87 |
   | D | 0.255 | 52.2 | 0.186 | yes | no / no / no | 0.87 |

10. **Methane can be oxidised instead of released.**
    - The template's C–P lyase rxn26459 releases CH4 in the cytosol. A methane monooxygenase (rxn00843/rxn00844) can then consume it.
    - That reaction was gap-filled in 1 arm-A MePn model (*Haloquadratum*) and gene-associated in 2 (A, B) or 5 (C, D).
    - In 1–2 growing MePn models per arm, CH4 export at maximum growth on the Pi-free MePn medium is 0.
    - Like the NH3 routes, this affects all arms alike and is handled in the final set, not in the comparison.

## Final models

The final set is arm A's builds after final-set repairs that apply to all 58 models alike, run after this decision. The repairs and the reasons for each are in `../models/README.md`; the scripts are `../pipeline/repair_model.py` and `../pipeline/finalize.py`. The comparison builds in `runs/` are unchanged.

## Results tables

### Metric medians (means in brackets) over the paired models

| Metric | Arm A | Arm B | Arm C | Arm D |
|---|---:|---:|---:|---:|
| Draft reactions with a gene association | 802 (821) | 1005 (1036) | 1380 (1434) | 1392 (1452) |
| 1. Gene-associated, non-gap-filled reactions (final model) | 802 (820) | 1004 (1035) | 1378 (1427) | 1370 (1443) |
| 2. Gap-filled reactions, guild-medium gap-fills | 97 (104) | 75 (87) | 62 (73) | 62 (72) |
| 2'. ATP-correction core gap-fill (not in rule) | 0.0 (0.6) | 0.0 (0.4) | 0.0 (0.4) | 0.0 (0.4) |
| 3. Guild-pathway gene support (fraction) | 1.000 (0.595) | 1.000 (0.905) | 1.000 (0.924) | 1.000 (0.924) |
| 4. ATP-test pass rate | 1.000 (0.959) | 1.000 (0.968) | 1.000 (0.962) | 1.000 (0.962) |
| 4'. ATP tests per model | 7.0 (8.3) | 10.0 (8.9) | 10.0 (10.2) | 10.0 (10.3) |
| 5. Blocked reaction fraction | 0.409 (0.405) | 0.435 (0.435) | 0.474 (0.471) | 0.472 (0.467) |
| 5'. Dead-end metabolite fraction | 0.197 (0.196) | 0.178 (0.182) | 0.181 (0.186) | 0.180 (0.185) |
| 6. Reactions with weak single-source non-RAST evidence only | 0 (0) | 0 (0) | 288 (288) | 256 (254) |
| Reactions with non-RAST evidence only (final model) | 0 (0) | 203 (215) | 590 (609) | 616 (626) |
| Final model reactions | 1025 (1035) | 1194 (1236) | 1592 (1626) | 1606 (1647) |
| Growth on guild medium (1/h) | 0.048 (0.346) | 0.048 (0.326) | 0.050 (0.423) | 0.050 (0.545) |
| 4. Models passing the free-ATP test | 55/58 | 55/58 | 57/58 | 57/58 |
| Models growing on the guild medium | 56/58 | 57/58 | 54/58 | 55/58 |
| MePn-obligate models | 0 | 0 | 0 | 0 |

### Per guild (medians; m3 mean in brackets)

| Guild | n | m3 A | m3 B | m3 C | m3 D | m1 A | m1 B | m1 C | m1 D | m2 A | m2 B | m2 C | m2 D | blocked A | blocked B | blocked C | blocked D |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| betaine_reducer | 18 | 0.00 (0.00) | 1.00 (0.94) | 1.00 (1.00) | 1.00 (1.00) | 686 | 856 | 1214 | 1266 | 154 | 132 | 115 | 112 | 0.382 | 0.401 | 0.439 | 0.438 |
| mepn_consumer | 39 | 1.00 (0.86) | 1.00 (0.89) | 1.00 (0.89) | 1.00 (0.89) | 984 | 1263 | 1652 | 1680 | 71 | 55 | 47 | 46 | 0.419 | 0.452 | 0.489 | 0.480 |
| methanogen | 1 | 0.82 (0.82) | 0.82 (0.82) | 0.91 (0.91) | 0.91 (0.91) | 525 | 666 | 914 | 915 | 104 | 73 | 57 | 57 | 0.396 | 0.412 | 0.462 | 0.458 |

### Decision rule, each merged arm against arm A

| Criterion | Arm B | Arm C | Arm D |
|---|---|---|---|
| (a) m1 ratio to A >= 1.10 | 1.252 yes | 1.720 yes | 1.709 yes |
| (b) m2 ratio to A <= 0.80 | 0.773 yes | 0.639 yes | 0.639 yes |
| (c) guild support not lower, every guild | betaine yes, mepn yes, methanogen yes | betaine yes, mepn yes, methanogen yes | betaine yes, mepn yes, methanogen yes |
| (d) no model gains free ATP; ATP-test pass rate not lower | 0 gain; 1.000 vs 1.000 yes | 1 gain; 1.000 vs 1.000 **no** | 1 gain; 1.000 vs 1.000 **no** |
| (e) blocked fraction change <= +2 pp | +2.69 pp **no** | +6.55 pp **no** | +6.36 pp **no** |
| **All met** | **no** | **no** | **no** |
| Paired view (median per-model change): m1 / m2 / blocked | +25.6% / -17.8% / +2.95 pp | +75.1% / -31.2% / +6.60 pp | +76.5% / -32.4% / +6.12 pp |

**Decision: arm A.**

### Per model

m1 = gene-associated non-gap-filled reactions; m2 = gap-filled; m3 = guild-pathway gene support; blk = blocked fraction; ATP = free-ATP test passed / ATP tests passed; grow = growth on guild medium.

| Guild | MAG | Genus | m1 A | m1 B | m1 C | m1 D | m2 A | m2 B | m2 C | m2 D | m3 A | m3 B | m3 C | m3 D | blk A | blk B | blk C | blk D | ATP A | ATP B | ATP C | ATP D | grow A | grow B | grow C | grow D |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| betaine | `R1_A_D2__maxbin2.Bin.123_sub` | BM714 | 686 | 847 | 1179 | 1203 | 141 | 123 | 104 | 103 | 0.00 | 1.00 | 1.00 | 1.00 | 0.333 | 0.364 | 0.458 | 0.443 | ok 2/2 | ok 2/2 | ok 2/2 | ok 2/2 | 1.377 | 1.331 | 1.996 | 1.996 |
| betaine | `R1_A_D2__metabat2.bin.020` | DATFJG01 | 847 | 1057 | 1411 | 1417 | 113 | 95 | 79 | 78 | 0.00 | 1.00 | 1.00 | 1.00 | 0.383 | 0.416 | 0.413 | 0.413 | ok 5/6 | ok 6/6 | ok 6/6 | ok 6/6 | 1.468 | 1.115 | 1.119 | 1.121 |
| betaine | `R1_B_D2__concoct.bin.112_sub` | Halanaerobium | 686 | 827 | 939 | 1255 | 157 | 147 | 144 | 118 | 0.00 | 1.00 | 1.00 | 1.00 | 0.405 | 0.429 | 0.367 | 0.405 | ok 7/7 | ok 10/10 | ok 10/10 | ok 10/10 | 0.570 | 0.576 | 1.087 | 2.663 |
| betaine | `R1_C_D1__concoct.bin.031` | UBA8670 | 739 | 915 | 1264 | 1283 | 151 | 132 | 112 | 111 | 0.00 | 1.00 | 1.00 | 1.00 | 0.343 | 0.397 | 0.441 | 0.440 | ok 8/9 | ok 10/10 | ok 10/10 | ok 10/10 | 1.081 | 1.076 | 0.598 | 1.108 |
| betaine | `R1_C_D1__maxbin2.Bin.003_sub` | M55B127 | 683 | 847 | 1159 | 1173 | 161 | 129 | 123 | 124 | 0.00 | 1.00 | 1.00 | 1.00 | 0.370 | 0.394 | 0.455 | 0.452 | ok 5/5 | ok 5/5 | ok 5/5 | ok 5/5 | 0.994 | 1.039 | 0.525 | 2.608 |
| betaine | `R1_sediment__concoct.bin.052_sub` | CAYRVO01 | 528 | 638 | 899 | 924 | 206 | 200 | 192 | 174 | 0.00 | 1.00 | 1.00 | 1.00 | 0.406 | 0.434 | 0.437 | 0.436 | ok 7/7 | ok 7/7 | ok 6/7 | ok 7/7 | 0.530 | 0.601 | **0** | 1.661 |
| betaine | `R1_sediment__metabat2.bin.186` | DAOWES01 | 566 | 706 | 1122 | 1129 | 137 | 115 | 92 | 91 | 0.00 | 1.00 | 1.00 | 1.00 | 0.323 | 0.378 | 0.454 | 0.455 | ok 5/5 | ok 5/5 | ok 5/5 | ok 5/5 | 1.062 | 1.092 | 0.071 | 0.071 |
| betaine | `R1_sediment__metabat2.bin.206` | Halanaerobium | 829 | 1028 | 1360 | 1381 | 116 | 107 | 89 | 83 | 0.00 | 1.00 | 1.00 | 1.00 | 0.435 | 0.415 | 0.425 | 0.424 | ok 8/9 | ok 9/10 | ok 10/10 | ok 11/13 | 1.382 | 1.070 | 2.629 | 1.120 |
| betaine | `R1_sediment__metabat2.bin.356` | QPJE01 | 679 | 814 | 1101 | 1156 | 172 | 148 | 129 | 126 | 0.00 | 1.00 | 1.00 | 1.00 | 0.401 | 0.415 | 0.434 | 0.430 | ok 5/6 | ok 5/6 | ok 5/6 | ok 5/6 | 1.034 | 1.082 | 2.564 | 2.592 |
| betaine | `R2_A_D1__maxbin2.Bin.018_sub` | M55B127 | 752 | 929 | 1311 | 1348 | 128 | 121 | 99 | 111 | 0.00 | 1.00 | 1.00 | 1.00 | 0.381 | 0.627 | 0.444 | 0.441 | ok 5/6 | ok 5/6 | ok 6/6 | ok 6/6 | 2.052 | **0** | 1.115 | 1.402 |
| betaine | `R2_B_D2__metabat2.bin.056` | M55B127 | 700 | 897 | 1278 | 1307 | 150 | 131 | 108 | 105 | 0.00 | 1.00 | 1.00 | 1.00 | 0.618 | 0.372 | 0.438 | 0.441 | ok 5/6 | ok 5/6 | ok 6/6 | ok 6/6 | **0** | 1.088 | 1.080 | 1.078 |
| betaine | `R2_B_D2__metabat2.bin.063` | M55B127 | 595 | 727 | 1053 | 1068 | 181 | 176 | 130 | 129 | 0.00 | 0.00 | 1.00 | 1.00 | 0.400 | 0.395 | 0.451 | 0.455 | ok 5/5 | ok 5/5 | ok 5/5 | ok 5/5 | 0.795 | 0.944 | 2.367 | 2.367 |
| betaine | `R2_C_D2__metabat2.bin.065` | Halanaerobium | 756 | 924 | 1275 | 1330 | 132 | 121 | 99 | 99 | 0.00 | 1.00 | 1.00 | 1.00 | 0.370 | 0.404 | 0.419 | 0.406 | ok 9/9 | ok 10/10 | ok 10/10 | ok 10/10 | 0.552 | 1.075 | 2.610 | 2.612 |
| betaine | `R2_sediment__metabat2.bin.141` | PWPR01 | 557 | 689 | 1021 | 1030 | 180 | 170 | 131 | 129 | 0.00 | 1.00 | 1.00 | 1.00 | 0.288 | 0.360 | 0.437 | 0.438 | ok 5/5 | ok 5/5 | ok 5/5 | ok 5/5 | 1.094 | 0.489 | 1.127 | 1.127 |
| betaine | `R2_sediment__metabat2.bin.162` | M55B127 | 681 | 866 | 1161 | 1192 | 169 | 153 | 123 | 120 | 0.00 | 1.00 | 1.00 | 1.00 | 0.383 | 0.383 | 0.447 | 0.432 | ok 5/5 | ok 5/6 | ok 5/6 | ok 5/6 | 0.737 | 0.580 | 1.098 | 2.617 |
| betaine | `R2_sediment__metabat2.bin.362_sub` | M55B127 | 711 | 890 | 1266 | 1276 | 163 | 140 | 143 | 142 | 0.00 | 1.00 | 1.00 | 1.00 | 0.359 | 0.405 | 0.433 | 0.436 | ok 2/2 | ok 2/2 | ok 1/2 | ok 1/2 | 1.400 | 1.878 | **0** | **0** |
| betaine | `R2_sediment__metabat2.bin.434_sub` | UBA8670 | 772 | 937 | 1349 | 1359 | 142 | 118 | 94 | 100 | 0.00 | 1.00 | 1.00 | 1.00 | 0.382 | 0.408 | 0.449 | 0.439 | ok 2/2 | ok 10/10 | ok 5/10 | ok 5/10 | 1.102 | 1.078 | 1.135 | 1.088 |
| betaine | `R2_sediment__metabat2.bin.435` | Halonatronomonas | 663 | 833 | 1250 | 1297 | 160 | 135 | 118 | 112 | 0.00 | 1.00 | 1.00 | 1.00 | 0.347 | 0.369 | 0.475 | 0.438 | ok 5/5 | ok 5/5 | ok 5/5 | ok 5/5 | 1.009 | 1.086 | 1.562 | 2.561 |
| mepn | `R1_A_D2__metabat2.bin.053` | Roseovarius | 1056 | 1315 | 1695 | 1740 | 54 | 35 | 31 | 33 | 1.00 | 1.00 | 1.00 | 1.00 | 0.434 | 0.464 | 0.482 | 0.480 | ok 16/16 | ok 16/16 | ok 16/16 | ok 16/16 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_B_D1__maxbin2.Bin.073` | Natrinema | 651 | 848 | 1153 | 1153 | 81 | 64 | 55 | 55 | 0.67 | 1.00 | 1.00 | 1.00 | 0.468 | 0.479 | 0.515 | 0.515 | ok 2/2 | ok 2/2 | ok 7/7 | ok 7/7 | 0.026 | 0.027 | 0.027 | 0.027 |
| mepn | `R1_B_D1__metabat2.bin.017_sub` | Roseovarius | 1103 | 1386 | 1758 | 1793 | 45 | 31 | 26 | 26 | 1.00 | 1.00 | 1.00 | 1.00 | 0.444 | 0.465 | 0.495 | 0.497 | ok 15/15 | ok 15/15 | ok 16/16 | ok 16/16 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_B_D1__metabat2.bin.040` | Roseovarius | 1018 | 1282 | 1676 | 1706 | 52 | 39 | 36 | 31 | 1.00 | 1.00 | 1.00 | 1.00 | 0.417 | 0.461 | 0.493 | 0.490 | ok 12/12 | ok 12/12 | ok 13/13 | ok 13/13 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_B_D2__metabat2.bin.065` | Roseovarius | 995 | 1271 | 1629 | 1667 | 62 | 51 | 47 | 43 | 1.00 | 1.00 | 1.00 | 1.00 | 0.450 | 0.473 | 0.496 | 0.493 | ok 15/15 | ok 15/15 | ok 16/16 | ok 16/16 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_B_H2O__concoct.bin.053` | Roseovarius | 1042 | 1276 | 1634 | 1666 | 48 | 39 | 33 | 31 | 1.00 | 1.00 | 1.00 | 1.00 | 0.415 | 0.441 | 0.476 | 0.473 | ok 12/12 | ok 12/12 | ok 13/13 | ok 13/13 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_C_D2__metabat2.bin.048` | Marivita | 1081 | 1368 | 1741 | 1761 | 58 | 41 | 35 | 34 | 1.00 | 1.00 | 1.00 | 1.00 | 0.434 | 0.463 | 0.476 | 0.472 | ok 13/13 | ok 15/15 | ok 15/15 | ok 16/16 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_sediment__concoct.bin.036` | CAZRZN01 | 1095 | 1415 | 1925 | 1928 | 54 | 49 | 44 | 40 | 1.00 | 1.00 | 1.00 | 1.00 | 0.415 | 0.458 | 0.500 | 0.500 | ok 15/15 | ok 15/15 | ok 16/16 | ok 16/16 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_sediment__concoct.bin.104` | Roseovarius | 1010 | 1251 | 1743 | 1746 | 63 | 54 | 55 | 53 | 1.00 | 1.00 | 1.00 | 1.00 | 0.466 | 0.486 | 0.532 | 0.530 | ok 3/6 | ok 3/6 | ok 7/7 | ok 7/7 | 0.042 | 0.042 | **0** | **0** |
| mepn | `R1_sediment__maxbin2.Bin.094_sub` | JAWLGC01 | 837 | 998 | 1340 | 1343 | 56 | 49 | 36 | 36 | 0.00 | 0.00 | 0.00 | 0.00 | 0.378 | 0.411 | 0.437 | 0.436 | ok 12/12 | ok 12/12 | ok 12/12 | ok 12/12 | 0.025 | 0.025 | 0.037 | 0.037 |
| mepn | `R1_sediment__metabat2.bin.121_sub` | Roseovarius | 1064 | 1372 | 1734 | 1772 | 56 | 46 | 36 | 36 | 1.00 | 1.00 | 1.00 | 1.00 | 0.421 | 0.454 | 0.478 | 0.480 | ok 13/13 | ok 14/14 | ok 15/15 | ok 15/15 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R1_sediment__metabat2.bin.142_sub` | RI-231 | 588 | 742 | 1227 | 1231 | 166 | 127 | 100 | 99 | 1.00 | 1.00 | 1.00 | 1.00 | 0.399 | 0.434 | 0.483 | 0.480 | ok 7/7 | ok 7/7 | **LOOP** 2/2 | **LOOP** 2/2 | 0.026 | 0.034 | 0.071 | 0.071 |
| mepn | `R1_sediment__metabat2.bin.258` | JACKII01 | 968 | 1292 | 1750 | 1763 | 62 | 46 | 38 | 38 | 1.00 | 1.00 | 1.00 | 1.00 | 0.429 | 0.461 | 0.492 | 0.482 | ok 10/10 | ok 9/10 | ok 10/10 | ok 10/10 | 0.048 | 0.048 | 0.049 | 0.049 |
| mepn | `R2A_sediment__concoct.bin.151` | WSZY01 | 777 | 916 | 1333 | 1337 | 77 | 62 | 51 | 51 | 1.00 | 1.00 | 1.00 | 1.00 | 0.419 | 0.434 | 0.483 | 0.479 | ok 8/8 | ok 8/8 | ok 10/10 | ok 10/10 | 0.044 | 0.034 | 0.045 | 0.045 |
| mepn | `R2A_sediment__concoct.bin.189_sub` | nan | 772 | 991 | 1397 | 1400 | 139 | 124 | 104 | 102 | 0.67 | 0.67 | 0.67 | 0.67 | 0.419 | 0.481 | 0.532 | 0.520 | **LOOP** 2/2 | **LOOP** 4/4 | ok 11/11 | ok 11/11 | 0.022 | 0.057 | 0.020 | 0.026 |
| mepn | `R2A_sediment__metabat2.bin.132` | CAMCWB01 | 1022 | 1318 | 1793 | 1798 | 71 | 61 | 48 | 46 | 1.00 | 1.00 | 1.00 | 1.00 | 0.451 | 0.470 | 0.497 | 0.500 | ok 12/12 | ok 13/14 | ok 18/18 | ok 18/18 | 0.023 | 0.035 | 0.048 | 0.048 |
| mepn | `R2_A_D1__concoct.bin.120` | Saliniramus | 1037 | 1300 | 1749 | 1788 | 55 | 43 | 39 | 29 | 1.00 | 1.00 | 1.00 | 1.00 | 0.414 | 0.432 | 0.472 | 0.470 | ok 12/13 | ok 13/13 | ok 13/13 | ok 13/13 | 0.048 | 0.048 | **0** | 0.050 |
| mepn | `R2_A_D2__metabat2.bin.055` | Haloplanus | 603 | 808 | 1102 | 1103 | 97 | 76 | 62 | 62 | 0.67 | 0.67 | 0.67 | 0.67 | 0.429 | 0.456 | 0.502 | 0.501 | ok 2/2 | ok 3/3 | ok 9/9 | ok 9/9 | 0.022 | 0.031 | 0.026 | 0.026 |
| mepn | `R2_A_H2O__metabat2.bin.026` | Halorubrum | 556 | 727 | 1067 | 1076 | 107 | 80 | 68 | 66 | 0.00 | 0.00 | 0.00 | 0.00 | 0.387 | 0.434 | 0.493 | 0.494 | ok 2/2 | ok 2/2 | ok 5/6 | ok 5/6 | 0.026 | 0.031 | 0.031 | 0.031 |
| mepn | `R2_B_D2__maxbin2.Bin.065` | Roseovarius | 1012 | 1265 | 1652 | 1680 | 59 | 43 | 36 | 35 | 1.00 | 1.00 | 1.00 | 1.00 | 0.421 | 0.449 | 0.493 | 0.496 | ok 12/12 | ok 12/12 | ok 13/13 | ok 13/13 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R2_B_D2__metabat2.bin.006_sub` | Roseovarius | 1058 | 1356 | 1858 | 1883 | 58 | 51 | 42 | 35 | 1.00 | 1.00 | 1.00 | 1.00 | 0.409 | 0.462 | 0.508 | 0.512 | ok 14/14 | ok 14/14 | ok 15/15 | ok 15/15 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R2_B_H2O__metabat2.bin.051` | Haloquadratum | 540 | 709 | 1019 | 1019 | 92 | 74 | 62 | 62 | 1.00 | 1.00 | 1.00 | 1.00 | 0.399 | 0.429 | 0.480 | 0.480 | ok 6/7 | ok 8/9 | ok 10/10 | ok 10/10 | 0.020 | 0.034 | 0.034 | 0.034 |
| mepn | `R2_C_D1__metabat2.bin.057_sub` | Desulfotignum | 837 | 1012 | 1503 | 1520 | 77 | 57 | 43 | 42 | 1.00 | 1.00 | 1.00 | 1.00 | 0.424 | 0.441 | 0.468 | 0.462 | ok 11/11 | ok 11/11 | ok 11/15 | ok 11/15 | 0.025 | 0.022 | 0.043 | 0.043 |
| mepn | `R2_C_D2__concoct.bin.063` | Roseibacterium | 1099 | 1423 | 1938 | 1959 | 62 | 48 | 39 | 34 | 1.00 | 1.00 | 1.00 | 1.00 | 0.396 | 0.441 | 0.497 | 0.497 | ok 14/15 | ok 14/15 | ok 16/16 | ok 16/16 | 0.048 | 0.048 | 0.050 | **0** |
| mepn | `R2_C_D2__concoct.bin.072` | JAVJPJ01 | 1060 | 1388 | 1876 | 1882 | 46 | 36 | 32 | 32 | 1.00 | 1.00 | 1.00 | 1.00 | 0.408 | 0.445 | 0.505 | 0.502 | ok 12/12 | ok 12/12 | ok 13/13 | ok 13/13 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R2_C_H2O__metabat2.bin.005` | Natronomonas | 586 | 795 | 1137 | 1137 | 88 | 69 | 57 | 57 | 0.67 | 1.00 | 1.00 | 1.00 | 0.452 | 0.470 | 0.527 | 0.527 | ok 4/5 | ok 5/5 | ok 4/5 | ok 4/5 | 0.020 | 0.028 | 0.035 | 0.035 |
| mepn | `R2_C_H2O__metabat2.bin.008` | Halonotius | 398 | 495 | 751 | 752 | 148 | 134 | 115 | 120 | 0.00 | 0.00 | 0.00 | 0.00 | 0.370 | 0.407 | 0.471 | 0.466 | ok 2/2 | ok 2/2 | ok 2/2 | ok 2/2 | 0.020 | 0.020 | 0.026 | 0.024 |
| mepn | `R2_restored_C_black__maxbin2.Bin.042` | Desulfotignum | 837 | 1035 | 1506 | 1520 | 64 | 53 | 47 | 49 | 1.00 | 1.00 | 1.00 | 1.00 | 0.424 | 0.404 | 0.460 | 0.451 | ok 13/13 | ok 13/13 | ok 17/17 | ok 17/17 | 0.039 | 0.039 | 0.045 | 0.045 |
| mepn | `R2_restored_C_black__metabat2.bin.011_sub` | Roseovarius | 984 | 1218 | 1698 | 1736 | 90 | 66 | 50 | 49 | 1.00 | 1.00 | 1.00 | 1.00 | 0.426 | 0.454 | 0.458 | 0.457 | ok 12/12 | ok 12/12 | ok 15/15 | ok 15/15 | 0.047 | 0.047 | 0.050 | 0.050 |
| mepn | `R2_restored_D_shore__concoct.bin.011_sub` | JASFMS01 | 991 | 1302 | 1763 | 1769 | 97 | 77 | 64 | 64 | 1.00 | 1.00 | 1.00 | 1.00 | 0.409 | 0.437 | 0.470 | 0.473 | ok 4/4 | ok 7/7 | ok 16/16 | ok 16/16 | **0** | 0.027 | 0.030 | 0.030 |
| mepn | `R2_sediment__concoct.bin.036_sub` | Roseovarius | 665 | 819 | 1250 | 1254 | 131 | 120 | 98 | 98 | 1.00 | 1.00 | 1.00 | 1.00 | 0.383 | 0.436 | 0.473 | 0.473 | **LOOP** 4/4 | **LOOP** 4/4 | ok 4/4 | ok 4/4 | 0.024 | 0.021 | 0.021 | 0.021 |
| mepn | `R2_sediment__concoct.bin.112` | CAJQQK01 | 979 | 1314 | 1777 | 1789 | 71 | 46 | 36 | 38 | 1.00 | 1.00 | 1.00 | 1.00 | 0.420 | 0.452 | 0.486 | 0.480 | ok 11/12 | ok 13/13 | ok 12/13 | ok 12/13 | 0.048 | 0.048 | 0.049 | 0.049 |
| mepn | `R2_sediment__concoct.bin.183_sub` | JBKDZR01 | 1041 | 1355 | 1832 | 1839 | 63 | 55 | 41 | 42 | 1.00 | 1.00 | 1.00 | 1.00 | 0.426 | 0.453 | 0.497 | 0.493 | ok 16/16 | ok 16/16 | ok 16/16 | ok 16/16 | 0.048 | 0.048 | 0.050 | 0.050 |
| mepn | `R2_sediment__metabat2.bin.112_sub` | Lentibacter | 853 | 1081 | 1491 | 1525 | 101 | 78 | 71 | 73 | 0.33 | 0.33 | 0.33 | 0.33 | 0.395 | 0.452 | 0.489 | 0.479 | ok 5/5 | ok 7/7 | ok 15/15 | ok 15/15 | 0.021 | 0.021 | 0.049 | 0.049 |
| mepn | `R2_sediment__metabat2.bin.128` | DATHSN01 | 782 | 1009 | 1572 | 1215 | 108 | 85 | 67 | 94 | 1.00 | 1.00 | 1.00 | 1.00 | 0.345 | 0.402 | 0.446 | 0.404 | ok 10/11 | ok 11/11 | ok 11/11 | ok 11/11 | 0.041 | 0.042 | 0.043 | 0.045 |
| mepn | `R2_sediment__metabat2.bin.154` | JAAAOG01 | 821 | 1031 | 1474 | 1482 | 122 | 97 | 64 | 63 | 1.00 | 1.00 | 1.00 | 1.00 | 0.379 | 0.411 | 0.446 | 0.449 | **LOOP** 2/3 | **LOOP** 3/4 | ok 13/13 | ok 13/13 | 0.031 | 0.071 | 0.028 | 0.028 |
| mepn | `R2_sediment__metabat2.bin.405` | Thermohalobaculum | 1006 | 1263 | 1683 | 1688 | 56 | 37 | 33 | 33 | 1.00 | 1.00 | 1.00 | 1.00 | 0.449 | 0.464 | 0.489 | 0.488 | ok 11/11 | ok 11/11 | ok 12/12 | ok 12/12 | 0.039 | 0.039 | 0.046 | 0.046 |
| mepn | `R2_sediment__metabat2.bin.414_sub` | DAINOF01 | 1032 | 1311 | 1871 | 1876 | 77 | 60 | 47 | 47 | 0.67 | 1.00 | 1.00 | 1.00 | 0.415 | 0.475 | 0.509 | 0.509 | ok 15/15 | ok 16/16 | ok 16/16 | ok 16/16 | 0.048 | 0.047 | 0.048 | 0.048 |
| mepn | `R2_sediment__metabat2.bin.434_sub` | UBA8670 | 772 | 937 | 1349 | 1359 | 138 | 119 | 95 | 98 | 1.00 | 1.00 | 1.00 | 1.00 | 0.354 | 0.405 | 0.446 | 0.442 | ok 2/2 | ok 10/10 | ok 10/10 | ok 10/10 | 0.023 | 0.021 | 0.038 | 0.038 |
| methanogen | `R1_A_D1__concoct.bin.007` | Methanolobus | 525 | 666 | 914 | 915 | 104 | 73 | 57 | 57 | 0.82 | 0.82 | 0.91 | 0.91 | 0.396 | 0.412 | 0.462 | 0.458 | ok 5/5 | ok 2/2 | ok 3/3 | ok 3/3 | 0.418 | 0.158 | 0.254 | 0.255 |
