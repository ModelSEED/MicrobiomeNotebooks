# Final genome-scale models of the 57 salterns guild species (2026 MAG set)

Built 2026-11-20/22. **57 models kept, 1 dropped**:

| Guild | Folder | Built | Kept |
|---|---|---:|---:|
| MePn consumers | `mepn_consumer/` | 39 | **38** |
| Betaine reducers | `betaine_reducer/` | 18 | **18** |
| Methanogen (*Methanolobus*) | `methanogen/` | 1 | **1** |

- *UBA8670* `R2_sediment__metabat2.bin.434_sub` is in both the MePn and the betaine guild. It has one model in each folder, built with each guild's protocol, so the 57 kept models cover 56 species.
- The community methane simulations and Figure 7 were **not** run; they are the next step.
- Nothing was written to KBase and nothing was committed.

## Files

| Path | What |
|---|---|
| `<guild>/<mag_id>.json` | the model (cobra JSON). Model notes `salterns_2026` list the gap-filled and hand-added reactions; reaction notes carry `gapfilled` (stage) or `hand_added`. |
| `reports/<guild>/<mag_id>.report.json` | the comparison build report (`build`) and the final-set repair report (`repair`): draft, ATP correction, each gap-fill, every QC value, and the evidence behind each gene-associated reaction |
| `models_summary.tsv` | one row per (species, guild) build: taxonomy, template and rule, size, gap-fill counts and lists, repairs, growth, obligacy, CH4/TMA, free ATP, ATP tests, blocked fraction, guild-pathway gene support, kept/dropped and why, md5 |
| `dropped/mepn_consumer/` | the dropped model and its report |
| `_repair_runs/A/`, `_repair_runs/logs/` | raw `repair_model.py` outputs and logs from TF (md5 list `_repair_runs/final_A.md5`, checked after copying) |
| `MANIFEST.md5` | md5 of every model and report above |

## How the models were made

### 1. Annotation: arm A (RAST), by the pre-registered rule

The four-arm comparison is in `../annotation_comparison/RESULTS.md`.
- No merged arm met all five criteria. All failed the +2 pp blocked-fraction limit, and C and D also failed the free-ATP criterion.
- The final models therefore use the RAST service functions (`../rast_service/functions/`) on the Prodigal genes of the KOfam tracks.
- Model gene ids are those Prodigal protein ids without the `<genome>|` prefix.

### 2. Template per genome (`template_rule` column)

| Rule | Template | Genomes |
|---|---|---:|
| Archaea | `ArchaeaTemplateV6.json` | 7 (all Halobacteriota) |
| Bacteria in Bacillota (and its GTDB splits), Actinomycetota, Chloroflexota or Patescibacteriota | `GramPosModelTemplateV7.json` | 19 (Bacillota 18, Actinomycetota 1) |
| Other Bacteria | `GramNegModelTemplateV7.json` | 31 (Pseudomonadota 26, Desulfobacterota 3, Acidobacteriota 1, Cloacimonadota 1) |

The counts are for the 57 kept models. The dropped *JASFMS01* (Pseudomonadota) also used GramNeg V7.
- The templates are the combined curated set of `../templates/` (see its `PROVENANCE.md`; md5s are in each report).
- They were chosen by file name, never by template id and never from ModelSEEDpy's bundled templates.
- The phyla come from `modelling_species.tsv` (GTDB R232).

### 3. Build (`../pipeline/build_model.py`)

- **Draft:** ModelSEEDpy `MSBuilder.build_base_model` logic (311f720), replicated step by step, with no RAST call inside the builder.
- **ATPM:** gapfill.ipynb cell 3.
- **ATP correction:** `MSATPCorrection` with Core-V6 and the 54 default ATP media (evaluate → determine → apply → expand → build tests). The tests go into every gap-fill.
- **Gap-filling:** `MSGapfill` with the template as database and `remove_unneeded_reactions=True`. Added reactions get bounds ±100 in the gap-filled direction, as MSGapfill integrates them.
- **After gap-filling:** every exchange is reset to (−1000, 1000), and every medium compound the model lacks gets an exchange.

**MePn consumers: the original protocol** (`cliff_mags/mac/salterns/gapfill.ipynb`)
- Cell 5 disables on the draft: rxn23850, rxn09193, rxn05315 and rxn05206 off; the `disable_lb` list made irreversible; H+ check on rxn15395 and rxn28338.
- Cell 72 additions to the in-memory gap-filling template:
  - `cpd25960_e`;
  - the hand-built C–P lyase `rxs43120` (CH4 released to e0) and MePn ABC transporter `rxs00001`;
  - the disable lists;
  - rxn26447's H+ fix, which was not needed: the curated rxn26447 is already balanced.
- Cells 74/76: rxn05602 (lactate transport) and the MePn module (`rxs00001`, rxn26456, rxn26457, `rxs43120`, rxn26447, rxn04651), added where absent.
  - The template lyase rxn26459, when annotated, is made irreversible, as in the scouting.
  - The hand-added `rxs00001` gets a GPR from the RAST PhnC/D/E genes, and `rxs43120` takes rxn26459's GPR.
- **Gap-fill 1** (scouting advice: first without MePn): the cell-5 simulation medium without MePn. Lactate 100, Pi 100, O2, NH3, sulfate, trace elements, thiamin; biomass ≥ 0.01.
- **Gap-fill 2** (cells 92–93): lactate 1 + MePn 1000, **no Pi**, biomass ≥ 0.02.

**Betaine reducers** (`correlational/betaine_gapfilling.ipynb`)
- Gap-fill on aerobic lactate (the cell-0 medium, without MePn per the scouting caution).
- Then `add_betaine_reductase_pathway()` (cell 2, same logic): rxn17220_c0 (betaine reductase), rxn07207_c0 (glycine reductase), TMA exchange, and a TMA transport where none exists.
- Then the guild check on lactate + betaine 0.1, the betaine simulation medium (`correlational/betaine_grdH_r232/README.md` §1).
- If the model cannot export TMA, a targeted gap-fill follows (TMA export ≥ 0.01). It was needed in 1 of 18: CAYRVO01 `R1_sediment__concoct.bin.052_sub`, which got the betaine:H+ symporter rxn05579. Every final model then exports TMA at the 0.1 betaine bound.

**Methanogen** (no original code exists)
- Archaea template.
- Gap-filled on a methanol + mono-/di-/trimethylamine medium: each at 100, CO2 100, H2S 1000 as the sulfur source, NH3, Pi, the cell-5 trace elements and thiamin, with no O2, nitrate or sulfate.
- Then a methane-export gap-fill (CH4 ≥ 0.01), because the biomass gap-fill did not need methanogenesis.
- **The manuscript states "dimethylamine + H2" for the original methanogen, but no code for that exists.** The model was also tested on DMA + H2 + CO2: it grows at 1.32 h⁻¹ and can export CH4 up to 300 with biomass unconstrained.

### 4. Final-set repairs (`../pipeline/repair_model.py`; applied to all 58 arm-A models after the decision)

These were added because QC of the comparison builds found three artefacts. They apply to every model alike. The comparison builds in `../annotation_comparison/runs/` are unchanged. The repair went through 4 versions on TF (`~/tf_runs/salterns_modelling/final_v1..v3` superseded, `final/` current); each version fixed what the previous one exposed. The final repair report of each model records its actions.

**1. Ammonia-oxidation block**

*Problem:* three MePn models made ATP from NH3 + O2 with no carbon (125–425 ATP). This is not ATP from nothing in the strict sense, because NH3 is in the free-ATP salts medium. But none of these organisms is a nitrifier. They oxidised NH3 to nitrite with nitrite reductases run backwards: rxn13647 (ubiquinol) or rxn03978 (cytochrome c), gap-filled in the energy-limited Pi-free stage. After that block, the re-gap-fill found the next route: hydroxylamine reductase rxn00109 run NH3 → NH2OH with NAD, plus rxn00570 NH2OH + O2 → nitrite. ModelSEEDpy's ATP-test media contain no NH3, which is why ATP correction passes these models.

*Fix:*
- The four reactions are restricted to their reductive direction in the model and in the gap-filling template.
- Every repair gap-fill also uses the three carbon-free media as ATP tests (ATPM ≤ 1e-5).

*Models affected:*
- Gap-filled N-oxidising reactions were restricted in 7 models (3 MePn, 4 betaine).
- The draft, RAST-annotated rxn00109 lost its oxidising direction in 23.

**1b. No gap-filled methane consumption**

*Problem:* in the energy-limited Pi-free stage, the gap-filler can make a model consume the CH4 that the template lyase rxn26459 releases in the cytosol.
- A methane monooxygenase was gap-filled in *Haloquadratum*.
- Reversed methyl-CoM reductase (anaerobic methane oxidation) appeared in a repair re-gap-fill.
- Either way, CH4 export is 0.

The original protocol's hand-built lyase releases CH4 outside the cell, and no annotation supports these routes.

*Fix:*
- Every template reaction with methane, except the transport rxn10471, may only produce methane.
- A gap-filled reaction loses its methane-consuming direction: rxn00843 was removed from *Haloquadratum*.

*Kept as annotated, flagged for Andrew:*
- **Two *Roseovarius* carry a RAST "Methane monooxygenase component A alpha chain (EC 1.14.13.25)" gene.** These are `R2_B_D2__metabat2.bin.006_sub` and `R2_restored_C_black__metabat2.bin.011_sub`, with reaction rxn00844.
  - In `R2_restored_C_black__metabat2.bin.011_sub` this reaction oxidises part of the lyase methane. CH4 export is 0.0245, against about 0.032 in the other *Roseovarius* at the same growth.
  - The annotation is probably a related di-iron monooxygenase. Decide whether to block it before simulating.
- The methanogen's Mcr (rxn03127) and the reversible lyase in the betaine model of *UBA8670* can also run in the methane-consuming direction. Neither medium contains methane.

**2. Growth restoration**

*Problem:* a stage whose model was below its minimum (biomass 0.01/0.02; CH4 or TMA export 0.01) was re-gap-filled, at most 2 passes, with the same medium, target and template, plus the free-ATP tests. The cause was either the blocks above or an integrated gap-fill solution that had lost a reaction.

*Re-gap-filled* (reactions added and stage):

| Model | Guild | Reactions added | Stage |
|---|---|---:|---|
| `R2_B_D2__metabat2.bin.056` (M55B127) | betaine | 1 (rxn00333, glycolate oxidase) | stage 1, growth 0 → 1.04 |
| `R1_A_D2__maxbin2.Bin.123_sub` (BM714) | betaine | 3 | stage 1 |
| `R1_C_D1__concoct.bin.031` (UBA8670) | betaine | 1 | stage 1 |
| `R1_sediment__metabat2.bin.356` (QPJE01) | betaine | 1 | stage 1 |
| `R2_sediment__metabat2.bin.141` (PWPR01) | betaine | 2 | stage 1 |
| `R2_B_H2O__metabat2.bin.051` (*Haloquadratum*) | MePn | 1 | stage 2 |
| `R2_sediment__concoct.bin.036_sub` (*Roseovarius*) | MePn | 11 | stage 2 |
| `R2_sediment__metabat2.bin.154` (JAAAOG01) | MePn | 25 | stage 2 |
| `R2A_sediment__concoct.bin.189_sub` (no genus) | MePn | 41 | stage 2 |

The last three were the NH3 models: without NH3 as an energy source, the Pi-free lactate-1 stage needed much more gap-filling to reach biomass 0.02.

**3. Betaine reductase GPR**

*Problem:* the templates carry the betaine reductase roles (ftr01398, ftr11139, ftr33907) but no complex. RAST therefore cannot associate genes with betaine reductase through the template, and the hand-added rxn17220_c0 had no GPR.

*Fix:* rxn17220_c0 gets the genes RAST annotates with those roles, noted `gpr_source = annotation-supported, not template-mapped`. All 18 betaine models have such genes. This does not change any flux.

**4. Methanogen Fpo test** (report only)

The template's F420H2:methanophenazine oxidoreductase rxn24614 is direction `<`. Freeing it does **not** make methanogenic growth possible (`fpo_freed_test` column). The model is unchanged.

## QC of the kept models (`models_summary.tsv`)

| | MePn consumers (38) | Betaine reducers (18) | Methanogen (1) |
|---|---|---|---|
| Growth on guild medium (1/h) | 0.020–0.048 (median 0.043), Pi-free lactate 1 + MePn | 0.53–2.05 (median 1.04), lactate + betaine 0.1 | 0.418, methanol + methylamines |
| Other growth | simulation medium 1.72–5.40 | aerobic lactate: all grow; **no model grows without O2** (as in the original aerobic protocol) | DMA + H2 (manuscript medium) 1.32 |
| Guild function at growth | CH4 export 0.013–0.048 at maximum growth (all > 0) | TMA export 0.10 in all (= betaine bound) | **CH4 0 at maximum growth** (see below) |
| `growth_triplet` | **0 obligate**: all grow with MePn blocked (≥ 1.72) and with Pi blocked (MePn as only P, ≥ 1.52) | – | – |
| ATP from carbon-free media (water; O2 + salts; nitrate + salts) | 0 in all | 0 in all | 0 |
| ATP tests passed | 30 of 38 pass all | 12 of 18 pass all | 5/5 |
| Blocked reactions (median) | 0.418 | 0.376 | 0.396 |
| Reactions / genes (median) | 1,165 / 923 | 960 / 608 | 718 / 509 |
| Gap-filled reactions (median, range) | 67.5 (45–180) | 154.5 (113–206) | 104 |
| Guild pathway gene-supported (RAST through the template) | C–P lyase all 3 reactions in 29/38. *Halonotius* `R2_C_H2O__metabat2.bin.008`, *Halorubrum* `R2_A_H2O__metabat2.bin.026` and JAWLGC01 `R1_sediment__maxbin2.Bin.094_sub` have none; 5 lack only PhnM (rxn26457); *Lentibacter* `R2_sediment__metabat2.bin.112_sub` has only PhnM. The MePn module is hand-added in every model, as in the original protocol | 0/18 through the template (no complex exists); 18/18 annotation-supported GPR on rxn17220 (repair 3) | 9/11 core reactions (methanol MtaA/MtaB missing) |

**ATP-test failures (reported, not dropped).** A model fails a test when it makes more than 1.2× the core model's ATP on that medium. No model makes ATP without carbon.

| Guild | Model | Failed test(s) |
|---|---|---|
| MePn | *Roseovarius* `R1_sediment__concoct.bin.104` | Etho.O2, Fum.O2, LLac.O2 |
| MePn | CAJQQK01, DATHSN01, *Haloquadratum*, *Natronomonas*, *Roseibacterium*, *Saliniramus* | Pyr.SO3 |
| MePn | JAAAOG01 | Etho.O2 |
| Betaine | DATFJG01 | Pyr.SO3 |
| Betaine | *Halanaerobium* `R1_sediment__metabat2.bin.206`, 2 M55B127, QPJE01, UBA8670 `R1_C_D1__concoct.bin.031` | Fum |

Pyr.SO3 failures come from a gap-filled assimilatory sulfite reductase used as an electron sink.

**The gap-filled reaction list** of every model is in `models_summary.tsv` (`gapfilled_reactions`, as `reaction:stage`) and in the model notes.

## Dropped (1)

**MePn consumer *JASFMS01* `R2_restored_D_shore__concoct.bin.011_sub`: no growth after gap-filling and repair.**
- The model cannot make ribose-5-phosphate, and so no nucleotides.
- Every gap-fill (2 in the build, 4 repair passes) returned a "solution" whose integrated model did not grow.
- A diagnostic gap-fill without ATP tests needs a single reaction, sedoheptulose-1,7-bisphosphatase (rxn01345). With it the model grows at 4.77 h⁻¹ with no ATP from carbon-free media, but it fails 2 of 4 ATP tests (Etho.O2, Fum.O2).
- So no growing solution that respects the ATP tests exists in the template. ModelSEEDpy then returns an infeasible solution rather than none.
- The original protocol dropped its 4 non-growing models. If Andrew wants this species, add rxn01345 and accept the two ATP-test failures.

## Caveats Andrew should see before simulating

1. **The methanogen is not methanogenic as built.**
   - *Methanolobus* grows (0.42 h⁻¹) as a methyl → acetate acetogen: CH4 is 0 at maximum growth, and with acetate export blocked growth is 0.
   - It can export CH4 (up to 100) only apart from growth, after the methane gap-fill added formate hydrogenlyase (rxn08518).
   - The arm-A draft has no working route that reduces methanophenazine: Fpo rxn24614 is `<` in the template, and the ATP correction filtered the `>` direction of Rnf (rxa45615).
   - In the KO arms (B–D), the KO-mapped methanol:CoM methyltransferase rxn10568 gives methanogenic growth on methanol (CH4 52–64 at maximum growth). Those models do not grow on any methylamine alone (`../annotation_comparison/RESULTS.md`, point 9).
   - This model needs curation before it is used for methane.
2. **MePn CH4 comes from rxn26459 (CH4 in the cytosol) or rxs43120 (CH4 outside), whichever the optimum picks.** Both give 1 CH4 per MePn. Two *Roseovarius* carry an annotated methane monooxygenase (repair 1b).
3. **Template direction items for the curation queue:**
   - rxn11678 (betaine reductase) and rxn11676 (glycine reductase) are written in the synthesis direction and set `>`. With KO evidence they form a null cycle with the hand-added rxn17220/rxn07207.
   - rxn24614 (Fpo) is `<`.
   - rxn13647, rxn03978, rxn00109 and rxn00570 are reversible, which opens NH3 oxidation in non-nitrifiers.
   - The betaine reductase roles have no complex.
4. **All guild media are aerobic, as in the original protocols.** No betaine model grows without O2 (Stickland fermenters).
5. **Gap-fill count convention.** `n_gapfilled` counts reactions added or reversed by the guild-medium gap-fills and the repairs. It excludes the ATP-correction core gap-fill (`atp_core_gapfill`) and the hand-added protocol reactions (`hand_added`).

## Software and provenance

- **Compute:** TF (`~/tf_runs/salterns_modelling`, env `~/venvs/salterns_modelling`), 14–18 parallel single-thread jobs at `nice 10`.
  - Python 3.12.3; cobra 0.29.1 (installed `--no-deps` because the local env runs it with pandas 3.0.3); optlang 1.9.1; Gurobi 13.0.1 with `Threads 1`, `TimeLimit 1200` (`gurobi.env`).
  - Full freeze: `../annotation_comparison/tf_env/pip_freeze.txt`.
- **ModelSEEDpy:** Andrew's fork, branch `sync-cshenry-main`, commit `311f72066736a44ea1c84799b595379a1bb14081`, installed on TF from a `git archive` of that commit (Python sources md5-identical to `~/Documents/ModelSEEDpy`).
- **ModelSEEDDatabase:** branch `tautomer-adenine-n1_20261113`, commit `41345b041f8c25d1e1aa8f476563165f7c759c15`. TF got a `git archive` of only the files `from_local()` reads. Neither repository was switched or modified.
- **Templates:** `../templates/` (md5s in `PROVENANCE.md` and in every report).
- **Code:**
  - `../pipeline/build_model.py`: v1 `7f8f7f65…` for 210 comparison builds; v2 `99652b0d…` for 22 re-runs; the current version `ea55c8ed…` only moves the QC into `qc_model()` (logic diff-checked identical) and was used by the repairs.
  - `repair_model.py`: `d0b5cf9e…`.
  - `finalize.py`: `e501e21f…`.
  - Each report records the md5s of the code that made it.
- **Re-running:**
  1. `pipeline/tf_setup.sh`, `pipeline/tf_sync.sh`.
  2. On TF: `bash code/run_tf.sh jobs.tsv 18`, then `bash code/run_repair_tf.sh A 14`.
  3. Locally: `pipeline/pull_tf.sh`, `pipeline/compare_arms.py`, copy back `final/A` (see `_repair_runs/final_A.md5`), `pipeline/finalize.py`.
  4. Every step skips finished work.

## Per model

"Without MePn / Pi-blocked" is `growth_triplet` growth with MePn uptake closed and with Pi uptake closed. "CH4 or TMA at growth" is the export at maximum growth on the guild medium.

| Guild | MAG | Genus | Phylum | Template | Reactions | Gap-filled | Growth (1/h) | Without MePn / Pi-blocked | Obligate | CH4 or TMA at growth | Free ATP | ATP tests | Pathway genes | Kept |
|---|---|---|---|---|---:|---:|---:|---|---|---:|---|---|---|---|
| betaine | `R1_A_D2__maxbin2.Bin.123_sub` | BM714 | Bacillota | GramPos V7 | 936 | 144 | 1.377 | – | – | TMA 0.10 | 0 | 2/2 | rxn17220: RAST roles | yes |
| betaine | `R1_sediment__concoct.bin.052_sub` | CAYRVO01 | Bacillota | GramPos V7 | 831 | 206 | 0.530 | – | – | TMA 0.10 | 0 | 7/7 | rxn17220: RAST roles | yes |
| betaine | `R1_sediment__metabat2.bin.186` | DAOWES01 | Cloacimonadota | GramNeg V7 | 785 | 137 | 1.062 | – | – | TMA 0.10 | 0 | 5/5 | rxn17220: RAST roles | yes |
| betaine | `R1_A_D2__metabat2.bin.020` | DATFJG01 | Bacillota | GramPos V7 | 1087 | 113 | 1.468 | – | – | TMA 0.10 | 0 | 5/6 | rxn17220: RAST roles | yes |
| betaine | `R1_B_D2__concoct.bin.112_sub` | Halanaerobium | Bacillota | GramPos V7 | 974 | 157 | 0.570 | – | – | TMA 0.10 | 0 | 7/7 | rxn17220: RAST roles | yes |
| betaine | `R1_sediment__metabat2.bin.206` | Halanaerobium | Bacillota | GramPos V7 | 1058 | 116 | 1.382 | – | – | TMA 0.10 | 0 | 8/9 | rxn17220: RAST roles | yes |
| betaine | `R2_C_D2__metabat2.bin.065` | Halanaerobium | Bacillota | GramPos V7 | 1020 | 132 | 0.552 | – | – | TMA 0.10 | 0 | 9/9 | rxn17220: RAST roles | yes |
| betaine | `R2_sediment__metabat2.bin.435` | Halonatronomonas | Bacillota | GramPos V7 | 929 | 160 | 1.009 | – | – | TMA 0.10 | 0 | 5/5 | rxn17220: RAST roles | yes |
| betaine | `R1_C_D1__maxbin2.Bin.003_sub` | M55B127 | Bacillota | GramPos V7 | 943 | 161 | 0.994 | – | – | TMA 0.10 | 0 | 5/5 | rxn17220: RAST roles | yes |
| betaine | `R2_A_D1__maxbin2.Bin.018_sub` | M55B127 | Bacillota | GramPos V7 | 988 | 128 | 2.052 | – | – | TMA 0.10 | 0 | 5/6 | rxn17220: RAST roles | yes |
| betaine | `R2_B_D2__metabat2.bin.056` | M55B127 | Bacillota | GramPos V7 | 956 | 151 | 1.043 | – | – | TMA 0.10 | 0 | 5/6 | rxn17220: RAST roles | yes |
| betaine | `R2_B_D2__metabat2.bin.063` | M55B127 | Bacillota | GramPos V7 | 875 | 181 | 0.795 | – | – | TMA 0.10 | 0 | 5/5 | rxn17220: RAST roles | yes |
| betaine | `R2_sediment__metabat2.bin.162` | M55B127 | Bacillota | GramPos V7 | 957 | 169 | 0.737 | – | – | TMA 0.10 | 0 | 5/5 | rxn17220: RAST roles | yes |
| betaine | `R2_sediment__metabat2.bin.362_sub` | M55B127 | Bacillota | GramPos V7 | 982 | 163 | 1.400 | – | – | TMA 0.10 | 0 | 2/2 | rxn17220: RAST roles | yes |
| betaine | `R2_sediment__metabat2.bin.141` | PWPR01 | Bacillota | GramPos V7 | 833 | 182 | 1.100 | – | – | TMA 0.10 | 0 | 5/5 | rxn17220: RAST roles | yes |
| betaine | `R1_sediment__metabat2.bin.356` | QPJE01 | Bacillota | GramPos V7 | 963 | 173 | 1.034 | – | – | TMA 0.10 | 0 | 5/6 | rxn17220: RAST roles | yes |
| betaine | `R1_C_D1__concoct.bin.031` | UBA8670 | Bacillota | GramPos V7 | 1005 | 152 | 0.588 | – | – | TMA 0.10 | 0 | 8/9 | rxn17220: RAST roles | yes |
| betaine | `R2_sediment__metabat2.bin.434_sub` | UBA8670 | Bacillota | GramPos V7 | 1033 | 142 | 1.102 | – | – | TMA 0.10 | 0 | 2/2 | rxn17220: RAST roles | yes |
| mepn | `R2_sediment__concoct.bin.112` | CAJQQK01 | Pseudomonadota | GramNeg V7 | 1150 | 71 | 0.048 | 5.35 / 4.78 | no | 0.0316 | 0 | 11/12 | 3/3 C–P lyase | yes |
| mepn | `R2A_sediment__metabat2.bin.132` | CAMCWB01 | Pseudomonadota | GramNeg V7 | 1206 | 71 | 0.023 | 2.89 / 2.30 | no | 0.0152 | 0 | 12/12 | 3/3 C–P lyase | yes |
| mepn | `R1_sediment__concoct.bin.036` | CAZRZN01 | Pseudomonadota | GramNeg V7 | 1274 | 54 | 0.048 | 5.40 / 4.82 | no | 0.0318 | 0 | 15/15 | 3/3 C–P lyase | yes |
| mepn | `R2_sediment__metabat2.bin.414_sub` | DAINOF01 | Pseudomonadota | GramNeg V7 | 1226 | 77 | 0.048 | 5.36 / 4.81 | no | 0.0318 | 0 | 15/15 | 2/3 C–P lyase | yes |
| mepn | `R2_sediment__metabat2.bin.128` | DATHSN01 | Actinomycetota | GramPos V7 | 988 | 108 | 0.041 | 5.01 / 4.07 | no | 0.0476 | 0 | 10/11 | 3/3 C–P lyase | yes |
| mepn | `R2_C_D1__metabat2.bin.057_sub` | Desulfotignum | Desulfobacterota | GramNeg V7 | 1037 | 77 | 0.025 | 2.78 / 2.51 | no | 0.0166 | 0 | 11/11 | 3/3 C–P lyase | yes |
| mepn | `R2_restored_C_black__maxbin2.Bin.042` | Desulfotignum | Desulfobacterota | GramNeg V7 | 1018 | 64 | 0.039 | 4.39 / 3.86 | no | 0.0255 | 0 | 13/13 | 3/3 C–P lyase | yes |
| mepn | `R2_C_H2O__metabat2.bin.008` | Halonotius | Halobacteriota | Archaea V6 | 645 | 148 | 0.020 | 2.43 / 2.03 | no | 0.0154 | 0 | 2/2 | 0/3 C–P lyase | yes |
| mepn | `R2_A_D2__metabat2.bin.055` | Haloplanus | Halobacteriota | Archaea V6 | 812 | 97 | 0.022 | 2.43 / 2.15 | no | 0.0164 | 0 | 2/2 | 2/3 C–P lyase | yes |
| mepn | `R2_B_H2O__metabat2.bin.051` | Haloquadratum | Halobacteriota | Archaea V6 | 738 | 92 | 0.029 | 3.33 / 2.82 | no | 0.0224 | 0 | 6/7 | 3/3 C–P lyase | yes |
| mepn | `R2_A_H2O__metabat2.bin.026` | Halorubrum | Halobacteriota | Archaea V6 | 762 | 107 | 0.026 | 3.16 / 2.63 | no | 0.0200 | 0 | 2/2 | 0/3 C–P lyase | yes |
| mepn | `R2_sediment__metabat2.bin.154` | JAAAOG01 | Pseudomonadota | GramNeg V7 | 1078 | 147 | 0.020 | 2.28 / 2.02 | no | 0.0134 | 0 | 2/3 | 3/3 C–P lyase | yes |
| mepn | `R1_sediment__metabat2.bin.258` | JACKII01 | Pseudomonadota | GramNeg V7 | 1134 | 62 | 0.048 | 5.38 / 4.80 | no | 0.0317 | 0 | 10/10 | 3/3 C–P lyase | yes |
| mepn | `R2_C_D2__concoct.bin.072` | JAVJPJ01 | Pseudomonadota | GramNeg V7 | 1216 | 46 | 0.048 | 5.37 / 4.79 | no | 0.0317 | 0 | 12/12 | 3/3 C–P lyase | yes |
| mepn | `R1_sediment__maxbin2.Bin.094_sub` | JAWLGC01 | Desulfobacterota | GramNeg V7 | 1006 | 56 | 0.025 | 3.06 / 2.52 | no | 0.0167 | 0 | 12/12 | 0/3 C–P lyase | yes |
| mepn | `R2_sediment__concoct.bin.183_sub` | JBKDZR01 | Pseudomonadota | GramNeg V7 | 1241 | 63 | 0.048 | 5.40 / 4.83 | no | 0.0319 | 0 | 16/16 | 3/3 C–P lyase | yes |
| mepn | `R2_sediment__metabat2.bin.112_sub` | Lentibacter | Pseudomonadota | GramNeg V7 | 1068 | 101 | 0.021 | 2.66 / 2.09 | no | 0.0138 | 0 | 5/5 | 1/3 C–P lyase | yes |
| mepn | `R1_C_D2__metabat2.bin.048` | Marivita | Pseudomonadota | GramNeg V7 | 1262 | 58 | 0.048 | 5.23 / 4.80 | no | 0.0317 | 0 | 13/13 | 3/3 C–P lyase | yes |
| mepn | `R1_B_D1__maxbin2.Bin.073` | Natrinema | Halobacteriota | Archaea V6 | 834 | 81 | 0.026 | 3.13 / 2.62 | no | 0.0201 | 0 | 2/2 | 2/3 C–P lyase | yes |
| mepn | `R2_C_H2O__metabat2.bin.005` | Natronomonas | Halobacteriota | Archaea V6 | 774 | 88 | 0.020 | 2.19 / 1.93 | no | 0.0155 | 0 | 4/5 | 2/3 C–P lyase | yes |
| mepn | `R1_sediment__metabat2.bin.142_sub` | RI-231 | Acidobacteriota | GramNeg V7 | 853 | 166 | 0.026 | 1.72 / 1.52 | no | 0.0172 | 0 | 7/7 | 3/3 C–P lyase | yes |
| mepn | `R2_C_D2__concoct.bin.063` | Roseibacterium | Pseudomonadota | GramNeg V7 | 1274 | 62 | 0.048 | 5.33 / 4.80 | no | 0.0317 | 0 | 14/15 | 3/3 C–P lyase | yes |
| mepn | `R1_A_D2__metabat2.bin.053` | Roseovarius | Pseudomonadota | GramNeg V7 | 1242 | 54 | 0.048 | 5.37 / 4.82 | no | 0.0318 | 0 | 16/16 | 3/3 C–P lyase | yes |
| mepn | `R1_B_D1__metabat2.bin.017_sub` | Roseovarius | Pseudomonadota | GramNeg V7 | 1274 | 45 | 0.048 | 5.37 / 4.79 | no | 0.0317 | 0 | 15/15 | 3/3 C–P lyase | yes |
| mepn | `R1_B_D1__metabat2.bin.040` | Roseovarius | Pseudomonadota | GramNeg V7 | 1180 | 52 | 0.048 | 5.34 / 4.79 | no | 0.0316 | 0 | 12/12 | 3/3 C–P lyase | yes |
| mepn | `R1_B_D2__metabat2.bin.065` | Roseovarius | Pseudomonadota | GramNeg V7 | 1180 | 62 | 0.048 | 5.37 / 4.80 | no | 0.0317 | 0 | 15/15 | 3/3 C–P lyase | yes |
| mepn | `R1_B_H2O__concoct.bin.053` | Roseovarius | Pseudomonadota | GramNeg V7 | 1207 | 48 | 0.048 | 5.38 / 4.82 | no | 0.0318 | 0 | 12/12 | 3/3 C–P lyase | yes |
| mepn | `R1_sediment__concoct.bin.104` | Roseovarius | Pseudomonadota | GramNeg V7 | 1185 | 63 | 0.042 | 4.68 / 4.21 | no | 0.0278 | 0 | 3/6 | 3/3 C–P lyase | yes |
| mepn | `R1_sediment__metabat2.bin.121_sub` | Roseovarius | Pseudomonadota | GramNeg V7 | 1236 | 56 | 0.048 | 5.35 / 4.81 | no | 0.0318 | 0 | 13/13 | 3/3 C–P lyase | yes |
| mepn | `R2_B_D2__maxbin2.Bin.065` | Roseovarius | Pseudomonadota | GramNeg V7 | 1185 | 59 | 0.048 | 5.38 / 4.80 | no | 0.0317 | 0 | 12/12 | 3/3 C–P lyase | yes |
| mepn | `R2_B_D2__metabat2.bin.006_sub` | Roseovarius | Pseudomonadota | GramNeg V7 | 1239 | 58 | 0.048 | 5.37 / 4.82 | no | 0.0318 | 0 | 14/14 | 3/3 C–P lyase | yes |
| mepn | `R2_restored_C_black__metabat2.bin.011_sub` | Roseovarius | Pseudomonadota | GramNeg V7 | 1197 | 90 | 0.047 | 5.15 / 4.66 | no | 0.0245 | 0 | 12/12 | 3/3 C–P lyase | yes |
| mepn | `R2_sediment__concoct.bin.036_sub` | Roseovarius | Pseudomonadota | GramNeg V7 | 907 | 142 | 0.020 | 2.50 / 2.00 | no | 0.0132 | 0 | 4/4 | 3/3 C–P lyase | yes |
| mepn | `R2_A_D1__concoct.bin.120` | Saliniramus | Pseudomonadota | GramNeg V7 | 1207 | 55 | 0.048 | 5.40 / 4.80 | no | 0.0318 | 0 | 12/13 | 3/3 C–P lyase | yes |
| mepn | `R2_sediment__metabat2.bin.405` | Thermohalobaculum | Pseudomonadota | GramNeg V7 | 1182 | 56 | 0.039 | 4.36 / 3.90 | no | 0.0258 | 0 | 11/11 | 3/3 C–P lyase | yes |
| mepn | `R2_sediment__metabat2.bin.434_sub` | UBA8670 | Bacillota | GramPos V7 | 1030 | 138 | 0.023 | 2.05 / 1.71 | no | 0.0268 | 0 | 2/2 | 3/3 C–P lyase | yes |
| mepn | `R2A_sediment__concoct.bin.151` | WSZY01 | Pseudomonadota | GramNeg V7 | 951 | 77 | 0.044 | 4.91 / 4.39 | no | 0.0290 | 0 | 8/8 | 3/3 C–P lyase | yes |
| mepn | `R2A_sediment__concoct.bin.189_sub` | (none) | Pseudomonadota | GramNeg V7 | 1051 | 180 | 0.020 | 2.17 / 2.00 | no | 0.0132 | 0 | 2/2 | 2/3 C–P lyase | yes |
| mepn | `R2_restored_D_shore__concoct.bin.011_sub` | JASFMS01 | Pseudomonadota | GramNeg V7 | 1266 | 152 | 0.000 | 0.00 / 0.00 | no | 0.0000 | 0 | 4/4 | 3/3 C–P lyase | **dropped** |
| methanogen | `R1_A_D1__concoct.bin.007` | Methanolobus | Halobacteriota | Archaea V6 | 718 | 104 | 0.418 | – | – | CH4 0.00 | 0 | 5/5 | 9/11 core | yes |

## 2026-11-24 updates

**Methanogen curated** (`methanogen/R1_A_D1__concoct.bin.007__curated.json`; `methanogen/CURATION.md`; `../pipeline/curate_methanogen.py`).
- **Why the arm-A build makes no CH4.** No route reduces methanophenazine:
  - Fpo was reverse-only, on the wrong F420 pool and moved no protons;
  - hdrDE was unmapped;
  - ATP correction had filtered the forward direction of Rnf.

  So heterodisulfide was never recycled and Mcr carried no flux. The build instead grew as an acetogen (CODH/ACS plus AcdAB).
- **The repair:** 19 gene-justified changes.
- **Validation:**
  - it is methanogenic on MMA, DMA and TMA, with textbook CH4/CO2/NH3 stoichiometry and no acetate;
  - it makes no free ATP.
- **Yield:** k = 126.5 mmol CH4 per g biomass formed, which is 7.9 g biomass per mol CH4. Published yields are 1.8–7.3, so the model sits at the efficient end. Across published ion couplings k ranges 119.5–159.1.
- **Methanol:** the genome has no MtaB or MtaC (confirmed by a six-assembly HMM search), so the model cannot use methanol. The methanol growth of the KOfam-arm builds (annotation_comparison runs B–D) is an artefact.
- **Yield file:** `methanogen/curation_k.json`.

**PrmA block** (`../pipeline/apply_prma_block.py`).
- **Change:** rxn00844 is closed in both *Roseovarius* models that carry the RAST "methane monooxygenase" gene (`R2_restored_C_black__metabat2.bin.011_sub`, `R2_B_D2__metabat2.bin.006_sub`).
- **Evidence:** the gene is propane 2-monooxygenase (PrmA). Its KOfam PrmA score is 1016 against 744 for MmoX. It shares 65 % identity with PrmA and 32–34 % with MmoX, and sits in a prmABCD operon. Details: CURATION.md.
- **Effect:** *Roseovarius*.2 now releases its MePn-derived methane.
- **MANIFEST.md5 is stale** for these two files; their new md5s are printed by the script.
