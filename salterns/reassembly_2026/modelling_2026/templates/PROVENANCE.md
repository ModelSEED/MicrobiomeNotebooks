# Salterns 2026 ModelSEED templates: provenance

2026-11-20. Built for the salterns MAG models (about 57 species: Gram-negative and Gram-positive bacteria, 7 archaea), to be built locally with ModelSEEDpy (`MSBuilder`, `MSATPCorrection`, `MSGapfill`) from RAST-annotated genomes.

**What this set is.** The newest upstream ModelSEEDTemplates files (`templates/v7.0` at ModelSEED/ModelSEEDTemplates `main` 2294299), plus every template curation Andrew Freiburger (`freiburgermsu`) made in his forks and submitted as PRs, except where one was superseded, was an explicitly rejected alternative, or was never meant to ship. Concretely: PR #14's branch `codif-curated` (b6fe86b), which already carries PRs #9, #11, #12 and #13, with PRs #15, #16 and #17 applied on top by their own apply scripts. ModelSEEDDatabase PR #325 is already in the templates through #14; PR #324 is not propagated (see below).

## Files

| File | Use as | md5 | Bytes | Reactions | Complexes |
|---|---|---|---|---|---|
| `GramNegModelTemplateV7.json` | genome-scale template, Gram-negative bacteria | `358a9119d63430c21e7be32a0040bbe5` | 14,604,789 | 8,583 | 3,297 |
| `GramPosModelTemplateV7.json` | genome-scale template, Gram-positive bacteria | `391f1ae6e6fee0d56925a627a026303a` | 14,606,765 | 8,583 | 3,297 |
| `ArchaeaTemplateV6.json` | genome-scale template, archaea | `7a47743e435dd1d16c745e26a8759bdb` | 24,468,731 | 8,599 | 3,323 |
| `Core-V6.json` | `template_core` for ATP correction, with every genome-scale template | `7ab51b2bc12136c066349c92c0e0bb07` | 523,552 | 252 | 468 |

- **The Archaea file is V6.** Upstream has no Archaea V7; `templates/v7.0/ArchaeaTemplateV6.json` was copied into the v7.0 folder by 2294299 ("Add Archaea template to v7.0 folder"). It is the newest Archaea template, and the one every curation below edits.
- **All three genome-scale files carry the same template id, `GramNegative.modeltemplate`**, and their biomass is named `GramNegativeBiomass` (an upstream quirk; the biomass compositions do differ). Choose the template by file, never by id.
- **Pass the templates explicitly.** Do not use `MSBuilder.auto_select_template`, `get_template` or ModelSEEDpy's bundled `modelseedpy/data/templates/template_gram_*.json`: those are a different, older lineage (the retired 2021 V3 release, per cshenry/ModelSEEDpy#109) and carry none of these curations.

```python
import json
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msbuilder import MSBuilder
T = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/modelling_2026/templates/"
load = lambda f: MSTemplateBuilder.from_dict(json.load(open(T + f))).build()
gramneg, grampos, archaea, core = (load(f) for f in ("GramNegModelTemplateV7.json",
    "GramPosModelTemplateV7.json", "ArchaeaTemplateV6.json", "Core-V6.json"))
builder = MSBuilder(genome, gramneg, ontology_term="RAST", template_core=core)
```

## Where the set lives in git

- **Branch `salterns-2026-combined`** in `~/Documents/ModelSEEDTemplates` (local only, not pushed), checked out in the new worktree `~/Documents/ModelSEEDTemplates-salterns2026`. Head **9696e5f**:

  | Commit | What |
  |---|---|
  | b6fe86b | base: `codif-curated`, PR #14 head (contains 21a8084 #11, 5516611 #12, fb58a76 replay of #9 + #13, a1a0986, 6fe204e, ed9b4ab, 5f53c0a, 2219d7e, b6fe86b) |
  | 06a1582 | replay of PR #15 (7fd1ab1) with its `apply_draft_loop_fix.py` |
  | e944dba | replay of PR #16 (0f6a5f6) with its `apply_thiolase_direction.py` |
  | 9696e5f | replay of PR #17 (e57ea75) with its `apply_proton_pyrophosphatase.py` |

- The four JSON files here are `git show salterns-2026-combined:templates/v7.0/<file>`, byte for byte. The PR #16 and #17 scripts also edit the v5.0 and v6.0 files on the branch, as the PRs do; those files are not part of this set.
- No existing branch, worktree or uncommitted file in Andrew's repositories was changed, and nothing was pushed.

## Curation inventory

Status is as of 2026-11-20, from `gh api -X GET search/issues` (author `freiburgermsu`, type `pr`) on ModelSEED/ModelSEEDTemplates, ModelSEED/ModelSEEDpy, cshenry/ModelSEEDpy, freiburgermsu/ModelSEEDpy and ModelSEED/ModelSEEDDatabase. Head commits come from the local branches and `git ls-remote` (`refs/pull/N/head`). "v7" means the four `templates/v7.0` files; GN/GP/Ar = GramNeg V7, GramPos V7, Archaea V6.

### ModelSEED/ModelSEEDTemplates

| # | Source branch | PR, status | Commit(s) | v7 files | Change | Included | Why |
|---|---|---|---|---|---|---|---|
| T1 | `remove-rxn39860` (on ModelSEED) | #9, open (2026-10-22) | e7fc410 | GN, GP, Ar (+ v5.0, v6.0) | Remove rxn39860, which with rxn07589 makes CO2 from nothing | **yes** | In the base: #14 replays it field by field in fb58a76 |
| T2 | `fix-irreversible-reaction-directions` (fork) | #11, open (2026-10-27) | 21a8084 | all four | rxn00006 (catalase, Core too), rxn10966, rxn00381, rxn00084, rxn08167 set to `>` | **yes** | Ancestor of #14 |
| T3 | `fix-opposite-template-directions` (fork) | #12, open, stacked on #11 | 5516611 | all four | rxn00206 `>`, rxn00380 `<`, rxn09566 `<`, rxn16828 `direction` `>`, rxn11937 `GapfillDirection` `<` | **yes** | Ancestor of #14 |
| T4 | `fix-reverse-penalty-sign` (fork) | #13, open, stacked on #11 | 09f6bab | all four | Negate the negative `reverse_penalty` where the database's graded dG favours forward: 1,234 per genome-scale template, 22 in Core | **yes** | In the base: replayed in fb58a76; every value equals 09f6bab's |
| T5 | `codif-combined` | part of #14 | fb58a76 | all four | Field-level replay of #9 and #13 onto #12 | **yes** | Ancestor of #14 |
| T6 | `codif-curated` stage A | #14, open (2026-11-01) | a1a0986 | all four | Electron counts: 1-electron cytochrome c and ferredoxin in 39 reactions (Core 3); NADPH restored in rxn04158 | **yes** | Base |
| T7 | stage B | #14 | 6fe204e | all four | Labels: cpd11620 charge 5 (all four), Mg-tetrapyrroles, A/AH2, flavodoxin, 118 acyl-ACPs; H+ of 63 reactions | **yes** | Base |
| T8 | stage C | #14 | ed9b4ab | GN, GP, Ar | rxn14415 `>`, rxn05902 `<`, rxn39679 `<`, rxn45223 `<` | **yes** | Base |
| T9 | stage D | #14 | 5f53c0a | all four | rxn10795, rxn00014, rxn00444, rxn09942, rxn20627 `>`; rxn13150 `<`; rxn00571 `<`, rxn00359 `>`, rxn00288 `<` (Core too); sul00002 `=` | **yes** | Base |
| T10 | stage E | #14 | 2219d7e | GN, GP, Ar | rxn00427, rxn02022, rxn02295, rxn02779, rxn11608, rxn08912, rxn08913, rxn09941, rxn22017, rxn48588 `>` | **yes** | Base |
| T11 | rxn10600 removal | #14 (head) | b6fe86b | GN, GP, Ar | Remove rxn10600, a cobalt chelatase without its electron acceptor | **yes** | Base |
| T12 | `fix-draft-atp-loops` (fork) | #15, open (2026-11-13), on #11 | 7fd1ab1 (1f84188, 627ffaa add validation only) | GN, GP, Ar | rxn23850 (GAPN), rxn05746, rxn08933, rxn08934, rxn08935, rxn08936 (maltodextrin glucosidases), rxn00014 set to `>` | **yes** | Applied by its own script; rxn00014 was already `>` from stage D, with the same four field values |
| T13 | `atp-cycle-directions` (fork) | #16, open (2026-11-18) | 0f6a5f6 | GN, GP, Ar (+ v5.0, v6.0) | rxn00598 (3-oxoadipyl-CoA thiolase) `<`, the thiolysis direction | **yes** | Applied by its own script |
| T14 | `hppase-pyrophosphate-proton-pump` (fork; clone `~/Documents/ModelSEEDTemplates-hppase`) | #17, open (2026-11-20) | e57ea75 | GN, GP, Ar (+ v5.0, v6.0) | Add rxn42569_c (H+-PPase, H2O + PPi -> 2 Pi + H+[e], `>`) and complex tcpx01431 for the orphan role ftr05360 | **yes** | Applied by its own script; entries identical to e57ea75's |
| T15 | `six-reaction-fix-nadx-variant` (fork) | none | 7d104ac | all four | #11's set with NadX (rxn05119, rxn05117) `>` instead of rxn08167 | no | Superseded by #11: a review found it blocks NAD synthesis in organisms without NadB |
| T16 | `negate-all-negative-reverse-penalties` (fork) | none; the alternative described in #13 | 6430128 | all four | Negate every negative reverse penalty (1,594 per genome-scale template, 26 in Core) | no | Conflicts with #13, which is the submitted form. The blanket form also flips the 80 reachable penalties whose reverse is favoured, and opened one new leak in #13's tests |
| T17 | `db-align-opposite-directions` (fork) | none ("Reference only", for issue #10) | 50ed6b4 | all four | 41 reactions set to the database direction | no | Reference only; #12 took the 5 that review judged template errors, and 24 of the 41 are database errors |
| T18 | `align-reversible-template-directions` (fork) | none ("Reference only") | 62fcc52 | all four | Every database-irreversible reaction restricted | no | Reference only |
| T19 | `codif-combined-aligned` (fork) | none ("Experiment") | 829013e | all four | Tier B and C direction alignment on codif-combined | no | Labelled an experiment |
| T20 | `freiburgermsu-patch-1`, `-patch-2`, `dev` | #4 closed unmerged; #5 merged 2023-10-11; #2, #3 merged 2023-09-05 | 65f9aff, 9849365, 22128ce | (2023 V5 layout) | V5 biomass/coefficient corrections; repository reorganisation | n/a | #2, #3 and #5 are already in `main`; #4 was superseded by #5 |

### ModelSEEDpy (freiburgermsu fork, cshenry/ModelSEEDpy, ModelSEED/ModelSEEDpy)

| # | Source branch | PR, status | Commit | Files | Change | Included | Why |
|---|---|---|---|---|---|---|---|
| P1 | `pr/gram-template-atp-directions` (worktree `ModelSEEDpy-gram-atp-directions`) | cshenry/ModelSEEDpy #109, open (2026-11-18) | 8506939 | bundled `template_gram_neg.json`, `template_gram_pos.json` | rxn00598 `<`; rxn05553, rxn08427, rxn08983, rxn08982 `>` | **covered**, no separate edit | rxn00598 is T13. rxn05553 is already `direction` `>`, `lower_bound` 0, `maxrevflux` 0 in all four v7 files, so `GapfillingPkg.convert_template_reaction` gives its candidate (0, 1000) even though `GapfillDirection` is `=`. rxn08427, rxn08982 and rxn08983 are not in any v7 file |
| P2 | `pr/mspredict-classifier` | cshenry #103, open | 3e7e6ca | bundled gram templates | Ships the bioseed files byte-identical | n/a | Not a curation |
| P3 | every other ModelSEEDpy PR (cshenry #82-#110, fork #19-#21, ModelSEED #143, #158) | various | | code | Gapfilling, ATP correction, thermodynamics, I/O | n/a | Code, no template content. The validation ran ModelSEEDpy 311f720, which carries the fork-merged #19 (gapfill reversals) and #110 (expansion-test filter) |

### ModelSEED/ModelSEEDDatabase

| # | Source branch | PR, status | Commit | Change | Included | Why |
|---|---|---|---|---|---|---|
| D1 | `ferredoxin-one-electron` | #325, open (2026-11-01) | 5389554a | cpd11620 charge 4 -> 5; Fd coefficients doubled in 43 reactions | **already in** (through #14) | 28 of the 43 reactions are in the genome-scale templates: 21 equal the corrected database form exactly and 7 up to the templates' own H+ coefficients. cpd11620 is +5 in all four files. The other 15 are not template reactions |
| D2 | `charge-imbalance-corrections` (worktree `ModelSEEDDatabase-e8446e3e`) | #324, open (2026-10-14) | e8446e3e | 33 reaction stoichiometries; 63 methyl-branched acyl-ACP charges -1 -> -2 | no | Database-side; never propagated into a template branch, and it does not map onto the templates (next section) |
| D3 | `tautomer-adenine-n1_20261113`, `thermo-impact_20261113`, `dgbyg-thermodynamics`; #326 (closed), #323 (open); the NAR 2026 manuscript and thermodynamics PRs | various | | pKa, protonation, free energies, manuscript | n/a | No template content |

### Not curations, or outside the forks and PRs

- **Issue #10** (ModelSEEDTemplates): the full direction audit. A reading, not an edit.
- **Project-local layers.** EmilyKin's member models use b6fe86b plus a local layer in that project (e.g. rxn48583 and rxn48586 `>`), and codiffusion uses its own curated v7 derivative. Neither is in a fork branch or a PR, so neither is included.

## Conflicts and how they were resolved

1. **#13 against its blanket alternative (T16).** Mutually exclusive treatments of the same 1,594 penalties. #13's rule is the submitted, documented one, and the one #14 replays: kept.
2. **#11 against the NadX variant (T15).** #11 replaced the variant after review: #11 kept.
3. **#14 stage D against issue #10 on rxn13150 and rxn00288.** #14 sets both to `<`. Issue #10 reads rxn13150 as mis-specified (remove or block it) and rxn00288 as reversible. #14 is the newer, submitted curation, with model evidence (rxn13150 `>` split water to O2 in 225 models; rxn00288 is an uphill step of a cross-member ATP cycle): **#14 kept**, disagreement recorded.
4. **#15 and #14 both set rxn00014.** Same values in all four fields: no conflict.
5. **#324 against #14's acyl-ACP labels.** For the 63 compounds, #324 sets charge -2 and keeps the database formula. #14 (6fe204e) also sets -2, but writes each acyl-ACP as ACP + acid - H2O, one H more (e.g. cpd11496: #14 C18H32N2O9PRS, #324 C18H31N2O9PRS), and adjusts the H+ of 60 template reactions to match. #14 is newer, and its labels balance the template reactions: **#14 kept**.
6. **#324 against #325 on rxn13974 and rxn14048.** #324 halves their ferredoxin to fit the old 4/6 convention; #325 and #14 keep them at 2 Fd with cpd11620 at +5, where they balance. **#325 / #14 kept** (newer, and the convention the templates use).
7. **rxn00598, `<` (#16) or the database's `=`.** #16 proposes `<` and documents why; `=` gave identical results in its tests. **#16 kept as submitted.**

## Why ModelSEEDDatabase #324 is not in the templates

The templates carry their own compound formulas and charges, which differ from the database's in many places. Each of #324's 33 corrected reactions was mapped into template compartments, oriented like the template entry, and balanced under the combined templates' own labels (`build/analysis/db_prs_vs_templates.py`, output `db_prs_vs_templates.tsv`):

- **2 are not template reactions:** rxn30129, and rxn10600 (removed by #14).
- **20 already balance in the templates.** Writing #324's form would unbalance them or bring in compounds the templates lack: rxn00689, rxn03121, rxn04158, rxn04160, rxn05266, rxn06031, rxn06335, rxn06586, rxn06621, rxn08437, rxn08480, rxn08550, rxn09507, rxn10002, rxn13690, rxn13727, rxn13838, rxn13974, rxn14048, rxn48582.
- **8 are imbalanced either way:** rxn00082, rxn00431, rxn01992, rxn03396, rxn04656, rxn04673, rxn43837, rxn45692.
- **3 would balance under #324's form:** rxn04046, rxn04219 and rxn41896. They are left as a follow-up, not edited here: no template branch made the change, and rxn04046's form swaps its product (cpd08369 -> cpd17535). rxn04219 is worth a look first: it is written with no electron acceptor (charge +2) and is reversible, so it can make electrons from protons in a gapfilled model.

## How the set was produced

`build/build_combined.sh` is the exact sequence; a rerun into a scratch worktree and branch on 2026-11-20 reproduced all four md5s (commit hashes differ, since they carry timestamps). In short:

```bash
REPO=~/Documents/ModelSEEDTemplates; WT=~/Documents/ModelSEEDTemplates-salterns2026
S=modelling_2026/templates/build/source_scripts           # verbatim copies, md5 below
git -C $REPO   show 7fd1ab1:validation/apply_draft_loop_fix.py         > $S/apply_draft_loop_fix.py
git -C $REPO   show 0f6a5f6:validation/apply_thiolase_direction.py     > $S/apply_thiolase_direction.py
git -C ~/Documents/ModelSEEDTemplates-hppase \
               show e57ea75:validation/apply_proton_pyrophosphatase.py > $S/apply_proton_pyrophosphatase.py
git -C $REPO worktree add -b salterns-2026-combined $WT b6fe86b       # codif-curated, PR #14
cd $WT
python $S/apply_draft_loop_fix.py --templates templates/v7.0           # PR #15
python $S/apply_draft_loop_fix.py --templates templates/v7.0 --check
git add templates/v7.0 && git commit                                   # 06a1582
python $S/apply_thiolase_direction.py --templates templates            # PR #16
python $S/apply_thiolase_direction.py --templates templates --check
git add templates && git commit                                        # e944dba
python $S/apply_proton_pyrophosphatase.py --templates templates        # PR #17
python $S/apply_proton_pyrophosphatase.py --templates templates --check
git add templates && git commit                                        # 9696e5f
for f in GramNegModelTemplateV7.json GramPosModelTemplateV7.json ArchaeaTemplateV6.json Core-V6.json; do
  git show salterns-2026-combined:templates/v7.0/$f > modelling_2026/templates/$f
done
```

`python` is `~/Documents/py_venv/bin/python`. Cherry-picking or merging #15, #16 and #17 is not possible: the Gram files are single-line JSON, so any two edits conflict. Each PR ships a script that edits only its fields and writes each file back in the serialization it already has, after checking that this reproduces the current bytes. Before any edit, a field-level check against b6fe86b found no conflicting field for any of the three.

Script md5s: `apply_draft_loop_fix.py` 6d59d7dcf87b918e4e8aabfb9fb85842, `apply_thiolase_direction.py` 72d204377938dcb5822bed40463e5dbe, `apply_proton_pyrophosphatase.py` 981fde3994d6aad357d3bbea7e56f9e8. Their JSON reports are in `build/logs/`.

**Combined against upstream 2294299**, per genome-scale file (the three are identical in this respect): 2 reactions removed (rxn39860, rxn10600), 1 added (rxn42569) and 1 complex added (tcpx01431); 1,363 reactions changed: `direction` 40, `GapfillDirection` 39, bounds and max fluxes 13-35, `reverse_penalty` 1,234, stoichiometry 102; 129 compounds relabelled (formula 121, charge 106) and 105 compartment-compound charges. Core-V6: 27 reactions (`reverse_penalty` 22, directions 4, stoichiometry 3) and cpd11620's charge.

## Validation

Scripts and outputs are in `build/validation/`: `check_templates.py` (template-level checks), `validate_build.py` (one MAG on one arm), `run_all.sh` (the six builds), `summarize.py`, `trace_draft_loop.py`, `results/*.json` and `logs/`.

**Environment.** `~/Documents/py_venv/bin/python` with the editable ModelSEEDpy in `~/Documents/ModelSEEDpy`, branch `sync-cshenry-main` at 311f720, no tracked-file changes (it was not switched); cobra 0.29.1; Gurobi 13.0.1, 2 threads per run. ModelSEEDpy's `reaction_expansion_test` also loads the local ModelSEEDDatabase checkout to order its candidates by reliability score; that checkout was on `tautomer-adenine-n1_20261113` (41345b04) and was the same for both arms. Run on this workstation rather than TF, because the task named this interpreter and editable install and TF has no ModelSEEDpy or Gurobi environment; the six runs took 3.3-6.6 min each.

### Templates

- **All eight files load** with `MSTemplateBuilder.from_dict(json).build()`: the four combined files and the four upstream 2294299 files.
- **Every curation is present, field by field.** For each included commit (parent -> commit), every field change in `templates/v7.0` was looked up in the combined files. All are there, and no later curation overrides an earlier one:

  | Curation | Field changes in the four v7 files | Present |
  |---|---|---|
  | #11 21a8084 | 72 | 72 |
  | #12 5516611 | 79 | 79 |
  | #13 09f6bab | 3,724 | 3,724 |
  | #9 e7fc410 | 3 | 3 |
  | #14 a1a0986 / 6fe204e / ed9b4ab / 5f53c0a / 2219d7e / b6fe86b | 120 / 1,186 / 48 / 138 / 120 / 3 | all |
  | #15 7fd1ab1 | 84 | 84 |
  | #16 0f6a5f6 | 18 | 18 |
  | #17 e57ea75 | 6 | 6 |

- **Nothing excluded leaked in.** Every field value of T15-T19 that also appears in the combined files is one an included curation sets (for example, #11's five reactions in the NadX variant).
- **DB #325:** cpd11620 is +5 in `compounds` and `compcompounds` of all four files.
- **Balance under the templates' own labels:** each genome-scale template goes from 237 element- or charge-imbalanced reactions upstream to 89; 148 are fixed and none becomes imbalanced. Core goes from 3 to 0.

### Draft models, ATP correction and gapfilling

Three MAGs (RAST genome JSONs in `jplfaria/Cliff MAGs/RASTified_MAG_json/`), each built on its template under both arms. CheckM and GTDB r232 from `MAG_QC/comparison_vs_2020.csv`:

- **Roseovarius** (Rhodobacteraceae), `Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs__`: 99.05% complete, 3.13% contaminated; GramNeg V7.
- **BM714** (Bacillota, Clostridia, Peptostreptococcales), `Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs__`: 100% / 1.06%; GramPos V7.
- **Methanolobus**, `Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_concoct_out.21.contigs__`: 98.37% / 2.61%; Archaea V6.

Pipeline: `MSBuilder.build_base_model(annotate_with_rast=False)` + `add_atpm`; `MSATPCorrection` with Core-V6 and ModelSEEDpy's 54 default ATP media (evaluate, determine, apply, expand to genome scale, build tests); `MSGapfill` with the ATP tests on complete medium (the 347 compounds of `ModelSEEDDatabase/Media/KBaseMedia.cpd` at 100), integrated with the gapfiller's reversals as `build_and_gapfill` does. RAST functions were split with `rast_client.split_annotation`.

Free-ATP test: maximise ATPM with every boundary reaction closed for uptake, except for the listed compounds. The test passes only if ATPM is 0. The three media carry no carbon: water + H+; + O2, NH3, phosphate, sulfate, K, Na, Mg, Ca, Cl; and the same with nitrate instead of O2. Values in the table are listed in that order.

| MAG | Arm | Draft reactions (gene-associated) | Free ATP, draft | Free ATP after ATP correction | Expansion-filtered | ATP test media selected (of 54) | Core gapfill | Gapfill new / reversed | Growth, complete medium | ATP tests passing, final | Free ATP, gapfilled | Final reactions |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Roseovarius | upstream | 1,201 (1,073) | **1000 / 1000 / 1000** | 0 / 0 / 0 | 7 | 16 | 0 | 47 / 0 | 21.57 | 17/17 | 0 / 0 / 0 | 1,248 |
| Roseovarius | combined | 1,201 (1,073) | **1000 / 1000 / 1000** | 0 / 0 / 0 | 3 | 16 | 0 | 47 / 0 | 17.16 | 17/17 | 0 / 0 / 0 | 1,248 |
| BM714 | upstream | 853 (750) | 0 / 0 / 0 | 0 / 0 / 0 | 1 | 5 | 0 | 94 / 0 | 8.67 | 6/6 | 0 / 0 / 0 | 947 |
| BM714 | combined | 853 (750) | 0 / 0 / 0 | 0 / 0 / 0 | 1 | 5 | 0 | 92 / 0 | 6.23 | 6/6 | 0 / 0 / 0 | 945 |
| Methanolobus | upstream | 645 (563) | **1000 / 1000 / 1000** | 0 / 0 / 0 | 3 | 4 | 8 | 72 / 0 | 3.96 | 5/5 | 0 / 0 / 0 | 725 |
| Methanolobus | combined | 645 (563) | 0 / 0 / 0 | 0 / 0 / 0 | 1 | 4 | 8 | 74 / 0 | 3.25 | 5/5 | 0 / 0 / 0 | 727 |

**No final model makes ATP from nothing.** After ATP correction, and again after gapfilling, all six models give ATPM 0 on all three media, with no biomass. Every ATP test passes, the empty-medium test included.

**Draft level:**
- **Upstream: 2 of 3 drafts leak.** Roseovarius and Methanolobus make ATP at the 1,000 bound through GAPN (rxn23850) run backwards, with GAPDH (rxn00781) and phosphoglycerate kinase (rxn01100). That is the loop PR #15 closes.
- **Combined: 1 of 3 drafts leaks.** Methanolobus is clean. The Roseovarius draft still leaks, through a second loop that the upstream draft also carries (`build/validation/logs/trace_draft_loop_Roseovarius.log`):
  - L-lactate:cytochrome c oxidoreductase (rxn00145, `=`) runs backwards and reduces pyruvate with ferrocytochrome c;
  - NAD-lactate dehydrogenase (rxn00499), NADH:ubiquinone oxidoreductase (rxn10122) and the bc1 complex (rxn13689) carry the electrons back, pumping protons for ATP synthase (rxn08173).
  - Knocking out rxn00145 alone removes it. This is the residual class PR #15 reports (95 of 1,181 drafts: one conversion with two reversible electron carriers). No curation covers it. ATP correction's expansion filter blocks `rxn00145 <` on both arms.
- **The expansion filter has less to do on the combined set.** On upstream it also had to block rxn14415, rxn15923, rxn11667 and rxn23850 (Roseovarius) and rxn23850 and rxn39860 (Methanolobus); the curations already fix all of these in the templates.

**Differences between the arms:**
- **The drafts hold the same reactions** except in Methanolobus, where rxn39860 is gone (#9) and rxn42569, the H+-PPase, is added with its gene (#17). Their bounds differ wherever a direction curation applies (for example, rxn23850 is forward-only on the combined set); that is what removes the GAPN leak above.
- **ATP yields on the aerobic and nitrate test media are lower on the combined set** (Roseovarius, maximum ATPM on the Ac.O2 test medium 9.5 -> 8.6, on Pyr.NO3 10.4 -> 8.7), consistent with #14's electron-count fixes (e.g. cytochrome c oxidase rxn14415, 3 -> 2 cyt c, which the upstream Roseovarius draft carries); this was not isolated by a knockout.
- **Growth on complete medium is 18-28% lower on the combined set.** The ATP-test media selected and the number of gapfilled reactions barely change, but the gapfill solutions differ (Roseovarius 18 of 47 reactions; Methanolobus 5 vs 7). On upstream, Methanolobus's solution used rxn00359 in reverse, an uphill step that #14 closes.

## Open items

- **Follow-ups that are not in any branch or PR:**
  - the three #324 reactions that would balance (rxn04046, rxn04219, rxn41896);
  - restricting rxn00145 (and its relatives in the two-carrier class) to `>`;
  - EmilyKin's rxn48583 and rxn48586 `>`;
  - #9's fuller fix (Co-precorrin-7 to the C44 convention);
  - #12's flagged database stoichiometries (rxn09566, rxn16828);
  - #12's rxn00148, rxn00157 and rxn00172, whose `direction` `<` and `GapfillDirection` `>` disagree, so gapfilling cannot add them.
- **89 reactions per genome-scale template stay imbalanced** under the templates' own labels (`check_templates.json`, `balance`).
- **None of PRs #9 and #11-#17 is merged upstream yet.** If maintainers change them before merging (for example, #16 to `=`), this set will differ from upstream's eventual v7.x.
