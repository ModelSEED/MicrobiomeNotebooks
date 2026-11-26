# Would a higher minimum contig length (or dropping MaxBin2) give better MAGs?

Local DAS Tool 1.1.7 re-runs on 8 of the 43 re-assembly units, 2026-11-11. Nothing was written to KBase.

## Answer

**No.** Each variant that raises the length floor or drops MaxBin2 loses MAGs, high-quality MAGs and original species. N50 rises, but much of that rise comes simply from deleting short contigs, and it stays well below the original MAGs' N50.

| projected for all 43 units (KBase CheckM) | MAGs ≥50/<10 | HQ ≥90/<5 | median N50 of passing bins (8 units) |
|---|---|---|---|
| V0 as run (current) | 1,324 | 323 | 5.2 kb |
| V1 MaxBin2 contigs ≥2,500 bp | ~1,150 (−13 %) | ~290 | 6.2 kb |
| V2 all binners ≥3,000 bp | ~1,055 (−20 %) | ~260 | 7.0 kb |
| V3 all binners ≥5,000 bp | ~710 (−47 %) | ~200 | 11.3 kb |
| V4 MetaBAT2 + CONCOCT only | ~1,085 (−18 %) | ~260 | 6.2 kb |
| V5 as run **plus** length-filtered copies as extra candidates | ~1,405 (+6 %) | ~350 | 5.3 kb |
| V6 as run **plus** ≥5,000 bp copies as extra candidates | ~1,365 (+3 %) | ~340 | 5.3 kb |

**Recommendation:**
- Keep the current binning (V0).
- Do not re-run MaxBin2 at 2,500 bp on KBase. That run would cost about 6 h median and up to 32 h per unit on KBase, plus a new extract, CheckM and GTDB-Tk. The expected result lies between V4 and V1, which both have fewer MAGs than now.
- The fragmentation is a property of the assembly, not of the binning floor:
  - 47–72 % of the bp in each assembly are in contigs under 2,500 bp.
  - Even at the 5 kb floor, a recovered original species reaches only 0.59× the original MAG's N50.
- The only gain comes from the ensemble V5, which leaves the floor alone and offers length-filtered bins as extra candidates:
  - The gain is small: +2 % expected MAGs in the per-sample units and +11 % in the one co-assembly tested. N50 does not change.
  - It would need KBase writes for every unit (a new BinnedContigs object, extract and CheckM), plus GTDB-Tk.
  - It is worth considering only for the 4 co-assemblies, and only if a KBase CheckM of the V5 bins confirms the gain.
- The 30 "lost" original MAGs are not a contig-length problem.
  - Every one of the 16 in scope is in the new assembly: 97.8–100 % ANI, with 79–99 % of the genome aligned.
  - They were either merged with other organisms into contaminated bins (CheckM contamination 16–22 %) or split across bins.
  - The fix for them is decontaminating or re-binning those bins, not a higher length floor.

## What was run

- **Units:** R1_A_D1, R2_B_D2, R2A_C_D1 and SF2_A_D1 (sediment cores); R1_B_H2O and SF2_B_H2O (water); R2_restored_D_shore (the original sample `R2_restored_DShore`); SF2_sediment (the smallest co-assembly, 1.66 Gbp).
- **DAS Tool command:** exactly as in `../../local/dastool_worker.py`: diamond, `--score_threshold 0.5 --duplicate_penalty 0.6 --megabin_penalty 0.5 --write_bin_evals --write_bins`, run under `nice -n 10`.
- **Reused single-copy-gene (SCG) calls:**
  - With `--proteins P --resume`, DAS Tool 1.1.7 only checks that P exists. It then reads `P.bacteria.scg` and `P.archaea.scg` and skips prodigal and diamond.
  - The production run's stored SCG tables were copied next to a one-line stub P, so every variant is a selection-only run on identical SCG calls. Each log says "Skipping SCG prediction" for both sets.
- **Assemblies:** each assembly was re-downloaded with `../../local/fetch_assembly.py`. The md5 was checked against the production run (`runs/<unit>/fetch.json`), and the FASTA was deleted after use.
- **Variants:** the inputs come from the stored `c2b/*.contigs2bin.tsv` and `<unit>.seqlength` files.
  - V0: the three binners as run.
  - V1: MaxBin2 bins with contigs <2,500 bp removed (a post-hoc stand-in for MaxBin2 at 2,500 bp).
  - V2 and V3: all three binners with contigs <3,000 or <5,000 bp removed.
  - V4: MetaBAT2 + CONCOCT only.
  - V5 and V6 (added): ensembles in which the bins as run compete with length-filtered copies of themselves. V5 adds MaxBin2 ≥2,500 bp and all binners ≥3,000 bp and ≥5,000 bp; V6 adds all binners ≥5,000 bp. DAS Tool picks a filtered copy only where it scores higher, and bin ids get the suffix `.geNNNN`.

**V0 reproduces production exactly.**
- In all 8 units the V0 contig2bin is identical to the stored `out/<unit>_DASTool_contig2bin.tsv`, contig for contig: 4–85 bins and 2,743–51,327 contig assignments.
- This holds both for the stored-SCG `--resume` path and for a fresh prodigal + diamond run (`V0fresh`, `-t 24`).
- The fresh SCG tables differ from the stored ones by at most 1 line per set in 4 units (the order of diamond ties), which has no effect on the result.

## Results (8 units)

These numbers are DAS Tool's own SCG estimates, so the variants are compared on the same basis. A bin passes at completeness ≥50 and redundancy <10. The expected KBase CheckM counts come from the calibration below.

| | V0 | V1 | V2 | V3 | V4 | V5 | V6 |
|---|---|---|---|---|---|---|---|
| DAS Tool bins | 310 | 261 | 233 | 156 | 248 | 327 | 317 |
| passing ≥50/<10 (DAS) | 262 | 237 | 213 | 149 | 219 | 278 | 272 |
| ≥90/<5 (DAS) | 76 | 62 | 50 | 38 | 56 | 78 | 78 |
| expected CheckM ≥50/<10 (actual V0: 219) | 225 | 191 | 171 | 114 | 180 | 236 | 231 |
| expected CheckM ≥90/<5 (actual V0: 44) | 49 | 40 | 33 | 24 | 37 | 50 | 50 |
| median N50 of passing bins (kb) | 5.2 | 6.2 | 7.0 | 11.3 | 6.2 | 5.3 | 5.3 |
| median contigs / size (Mb) | 500 / 2.83 | 384 / 2.54 | 322 / 2.45 | 181 / 2.18 | 394 / 2.61 | 471 / 2.66 | 487 / 2.83 |
| Mbp in passing bins | 720 | 595 | 520 | 336 | 568 | 744 | 743 |
| same-sample originals recovered (of 105; MQ+ of 95) | 75 (71) | 67 (63) | 61 (57) | 51 (48) | 62 (58) | 76 (72) | 76 (72) |
| originals gained / lost vs V0 | – | +1 / −9 | +1 / −15 | +1 / −25 | +1 / −14 | +1 / 0 | +1 / 0 |
| recovered originals: median % of genome covered; N50 vs original | 80 %; 0.38× | 78 %; 0.42× | 77 %; 0.48× | 68 %; 0.59× | 80 %; 0.44× | 77 %; 0.38× | 79 %; 0.38× |
| unit's current MAGs recovered (of 219) | 208 | 188 | 172 | 123 | 175 | 208 | 209 |
| V0 passing organisms dropped / new organisms | – | 38 / 13 | 68 / 19 | 126 / 13 | 53 / 10 | 2 / 18 | 0 / 10 |

- **Recovery rule:** ANI ≥95 and the larger aligned fraction ≥50 (skani 0.4.0), the same rule as `../../comparison/`. "Same-sample originals" are the original MAGs of that sample; for SF2_sediment they are the 28 originals of the six SF2 cores.
- **Higher floors only lose:**
  - Every variant that raises a floor loses originals and current MAGs, and gains at most one.
  - Losses are worst in the water units: R1_B_H2O recovers 6 originals at V0 but 2 at V1; SF2_B_H2O recovers 5 at V0 but 2–3 at V1, V2 and V4.
  - The co-assembly barely changes at V1 or V2.
- **V5 adds new organisms:** its 18 new organisms share almost no contigs with V0's passing bins (median 0 % of bp).
  - Only 2 of them match a current MAG.
  - 4 are <0.5 Mbp. Three are in SF2_B_H2O, at 0.18–0.29 Mbp with 61–71 % SCG completeness. At that size they are probably organelle or plasmid bins and would likely fail CheckM. The fourth is a 0.49 Mbp bin in SF2_sediment.
  - Most of the gain is in SF2_sediment (+10) and SF2_B_H2O (+5, including the three small bins). R1_A_D1, R2A_C_D1 and SF2_A_D1 gain nothing.

Per-unit numbers are in `results/variant_metrics.tsv` and `results/report_tables.md`.

## The 30 "lost" originals

"Lost" means the original MAG's species cluster has no re-assembly member: 58 MAGs in all. The 30 are the ones that are medium or high quality and have no re-assembly MAG above 80 % ANI. 23 of the 58, including 16 of the 30, belong to the samples of these 8 units.

- **All in the assembly:** skani of each original against its unit's whole assembly gives 97.8–100 % ANI, with a median of 88 % of the genome aligned (range 62–99 %; 79–99 % for the 16 of the 30).
- **Contaminated or split in the binning** (raw MetaBAT2/MaxBin2/CONCOCT bins, `results/lost_trace.tsv`). Of the 16 of the 30:
  - 8 sit mostly (≥50 %) in a chimeric raw bin, with DAS Tool redundancy ≥10 % (for example Kariarchaeum 94 % of the genome in a bin at 89/37, JBBDSH01, Halofilum, Guyparkeria).
  - 3 have a clean raw bin, and DAS Tool selected a bin of each: Psychroflexus, Vreelandella, and Fodinibius in R2_restored.
  - The other 5 are split across bins, so the best single bin covers 40–71 % of the genome.
- **Rejected by KBase CheckM for contamination:** the V0 DAS Tool output bins that do match lost originals failed at 15–22 %:
  - Psychroflexus: CheckM 84.5/20.4, while DAS Tool estimated 98/8.
  - Vreelandella 80.3/17.4, Halofilum 80.8/18.0, Fodinibius 82.7/16.4, D2472 61.9/22.1, Aquiluna 75.2/19.2.
- **Variants:**
  - V1, V2 and V4 recover no lost original that V0 does not already bin.
    - V1 and V2 lose the MaxBin2 bins that carry Psychroflexus, Vreelandella, Fodinibius, Aquiluna and D2472.
    - V4 keeps Psychroflexus, in a MetaBAT2 bin that covers 52 % of the genome, and loses the other four.
  - V3, V5 and V6 make one more lost original pass DAS Tool: JBBDSH01 at 67/0 in R2_restored_D_shore.
    - Removing its short contigs breaks it out of a chimeric MetaBAT2 bin.
    - Kariarchaeum also gets its own bin, but at 84/13 it still fails.
    - V5 and V6 keep every lost original that V0 binned.
  - Whether CheckM would pass any of these is untested.

## Calibration: DAS Tool SCG estimates vs KBase CheckM

Basis: all 1,792 production V0 bins of the 43 units (`results/calibration_bins.tsv`, `results/calibration.tsv`).

- **Completeness:** Pearson 0.77, Spearman 0.81. CheckM ≈ 1.02 × DAS − 6.6, with a median of 3 points below DAS Tool.
- **Redundancy vs contamination:** Pearson 0.68, Spearman 0.60. CheckM is a median 1.9 points higher, and DAS Tool badly underestimates contamination in chimeric bins (see above).
- **Pass rates:** among bins that pass DAS Tool (≥50/<10), 81 % pass CheckM (1,248 of 1,545); among bins that fail it, 31 % still pass CheckM (76 of 247). DAS Tool HQ → CheckM HQ is only 51 % (271 of 528).
- **The pass rate depends on DAS Tool completeness:** for bins that pass DAS Tool, 55 % pass CheckM at 50–70 % DAS completeness (381 bins), 84 % at 70–90 % (568) and 94 % at ≥90 % (596).
  - The projections therefore weight each passing bin by its completeness band, which matters because the new and surviving bins of the variants have low DAS completeness (median 59–75 for the new organisms).
  - A single flat 81 % rate (`exp_checkm_mq_two_rate`) would give every variant 2–9 more expected MAGs in the 8 units.
- **The model reproduces the 8 units' actual V0 counts:** expected 225 vs actual 219 at ≥50/<10, and 49 vs 44 at ≥90/<5 (two rates for HQ).
- **How the projections use it:** they are ratios to V0, so a constant bias cancels; a bias that changes with the variant does not. V5 picks filtered copies partly because they lower DAS Tool redundancy, so whether its bins really are less contaminated is untested.

## Post-hoc filtering is not a MaxBin2 re-run at 2,500 bp

Source: `results/maxbin2_posthoc.tsv`.

- **What V1 removes:**
  - MaxBin2 bins hold 47–69 % of their bp, and 42–67 % of their SCG hits, on contigs <2,500 bp.
  - V1 removes a per-bin median of 53–77 % of the bp but empties no bin, so it keeps MaxBin2's 1,000-bp partition of the long contigs.
- **What a re-run would change:** MaxBin2 sets its number of bins from marker-gene copy counts on contigs ≥ min_contig_length.
  - **Proxy:** the median copy number of DAS Tool's 51 bacterial SCGs drops from 20–362 on contigs ≥1,000 bp to 6–161 on contigs ≥2,500 bp, a ratio of 0.30–0.56.
  - **How good the proxy is:** at ≥1,000 bp it gives 0.48–1.02 of the actual MaxBin2 bin counts.
  - A re-run would therefore probably start from about half as many seeds and re-fit its EM, so its bins could differ from V1's. That cannot be measured locally: the multi-sample coverage is computed on KBase.
- **Bounds:**
  - V4 (no MaxBin2) is the floor for anything MaxBin2 at 2,500 bp could add.
  - **Replacement rate:** in the 8 units, a MetaBAT2/CONCOCT bin of the same organism replaces only 55 % (64 of 116) of the MaxBin2-derived passing bins.
  - **Projection:** 533 of the 1,324 current MAGs are MaxBin2 bins, so V4 projects to about 1,085 MAGs, matching the calibrated estimate of 1,083.
  - V1 already sits below V0, so a re-run is not expected to beat the current set.

## Cost

| step | time | measured where |
|---|---|---|
| download one assembly | 4–86 s | this run |
| DAS Tool selection, stored SCGs (5 variants concurrently) | 3–138 s per unit (V5: up to 151 s) | this run; the SCG tables exist for all 43 units in `work/dastool/` |
| fresh prodigal + diamond (`-t 24`) | 23–421 s per unit | this run |
| KBase MaxBin2, per unit | median 6.0 h, up to 32 h (co-assemblies) | `ledger.json`, submit to finish |
| KBase extract / CheckM, per unit | median 3.5 h / 4.5 h | `ledger.json` |
| GTDB-Tk R232, all MAGs (local) | about 2.5 h | `ledger.json` |

The ledger times are quantized by the watcher's 30-minute cycle.

A local-only variant such as V5 costs less than an hour of compute for all 43 units, but adopting it needs KBase writes per unit (a new BinnedContigs object, extract and CheckM), then GTDB-Tk and a re-run of `../../compare_vs_original.py`. A MaxBin2 re-run adds a KBase MaxBin2 job per unit before all of that.

## Caveats

- **Unit coverage:**
  - These are 8 of 43 units.
  - The co-assembly projection rests on SF2_sediment alone, the smallest co-assembly: 1.66 Gbp vs 1.84–2.24 Gbp for the other three.
  - The 4 co-assemblies hold 548 of the 1,324 MAGs and 189 of the 323 HQ MAGs.
- **N50 under a higher floor:** the N50 gains are partly mechanical, since the floor deletes short contigs. Read them together with bin size and genome coverage: the median share of an original covered by its matching bin falls from 80 % to 68 % at V3.
- **The assembly check:**
  - skani of an original against the whole assembly treats the assembly as one genome.
  - Its AF is the share of the original covered by any contigs, not by one bin.

## Files

- `scripts/run_unit.py`: fetches the assembly, builds the variant inputs, runs DAS Tool (`--only=V5,V6` and `--fresh-v0` options), runs the skani check against the assembly, then deletes the assembly.
- `scripts/prep_refs.py`: builds `results/originals.tsv` (lost58 and lost30 flags) and `results/current_mags.tsv`.
- `scripts/analyze.py`: metrics, the skani triangle over every bin and reference, recovery, organism turnover, calibration, the MaxBin2 post-hoc table, the 43-unit projection, costs and `results/report_tables.md`.
- `scripts/trace_lost.py`: lost originals vs every raw binner bin, written to `results/lost_trace.tsv`. The raw-bin FASTA and assemblies are deleted afterwards.
- `scripts/units.py`: the experiment units and the sample-name mapping.
- `runs/<unit>/`:
  - `fetch.json`: the md5s of the deleted assembly.
  - `run_checks.json`: the V0 reproduction check, the resume check, timings and filter counts.
  - `scg/`: the stored SCG tables and the stub.
  - `skani_originals_vs_assembly.tsv`.
  - `<variant>/{inputs,out}`: the bin FASTA, about 5 GB, ignored by git and regenerable. `V0fresh` keeps its SCG tables but not its bins.
- `results/`:
  - Variant results: `variant_totals.tsv`, `variant_metrics.tsv`, `bins.tsv`.
  - Recovery: `recovery_originals.tsv`, `recovery_current.tsv`, `organism_turnover.tsv`.
  - Lost originals: `lost_originals.tsv`, `lost_trace.tsv`.
  - Calibration: `calibration_bins.tsv`, `calibration.tsv`.
  - Projection and costs: `maxbin2_posthoc.tsv`, `extrapolation_43.tsv`, `costs.tsv`.
  - skani: `skani_all.tsv` (with its `skani_inputs.txt`).
  - Summaries: `summary.json`, `report_tables.md`.

**Reproduce:**

```bash
cd salterns/reassembly_2026/binning_experiments/dastool_min_contig
~/Documents/py_venv/bin/python -I scripts/prep_refs.py
for u in SF2_A_D1 R2A_C_D1 R1_A_D1 R2_B_D2 R1_B_H2O SF2_B_H2O R2_restored_D_shore SF2_sediment; do
  ~/Documents/py_venv/bin/python -I scripts/run_unit.py $u --fresh-v0; done
~/Documents/py_venv/bin/python -I scripts/analyze.py
~/Documents/py_venv/bin/python -I scripts/trace_lost.py
```

The reproduce commands above refer to the variant layout of the final script. In this run, V0–V4 (with `--fresh-v0`) and V5–V6 (with `--only=V5,V6`) ran in two passes.
