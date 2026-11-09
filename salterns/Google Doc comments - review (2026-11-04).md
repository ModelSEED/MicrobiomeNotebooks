# Google Doc comments, reviewed against the recent work (2026-11-04)

**Source:** the "Metabolic Manuscript" Google Doc (owner Susannah Tringe; last edited 2026-11-03), read through the Google Docs and Drive APIs on 2026-11-04. It has 237 comment threads, 70 open and 167 resolved, and this review covers the 70 open ones. Locations come from the comment markers in Google's own exports. Eleven open threads are no longer attached to any text, because the text they were on was deleted; they are listed at the end.

**Since the 2 July export** (`Metabolic Manuscript.docx`), the Google Doc has changed as follows:
- **Accepted suggestions:** on 27 August, A. Freiburger accepted 8 of Susannah's suggestions, all in the first Introduction paragraphs.
- **New suggestions:** on 3 November, Susannah added 37 suggestions. They include a new author (Kyle Hartman), renamed Methods headings and abbreviation expansions.
- **Comments:** 13 new threads (Susannah, 3 November), 1 new reply (Susannah, on the betaine-evidence thread), and 1 thread resolved (27 August).

All of this is in `Metabolic Manuscript - Google Doc mirror (2026-11-03).docx`, the base of the R232 review (see the notes file).

**Status labels:**
- **Addressed:** the R232 review (tracked changes) deals with it.
- **Partly addressed:** the review deals with part of it, or only in its own new text.
- **Answered by new analysis:** recent work answers the question, sometimes with a different outcome than the comment expected.
- **Open:** a decision or detail only the authors can supply.
- **Outdated:** no longer applies, or is already done in the text.

## Summary

| Status | Threads |
|---|---|
| Addressed | 16 |
| Partly addressed | 3 |
| Answered by new analysis | 11 |
| Answered in its thread | 1 |
| Open, for the authors | 36 |
| Outdated | 3 |
| **Total** | **70** |

- **DNA normalization is the main theme.** Several threads (Figure 7, Methods P38 and P40, Results P69, Figure S8) expected DNA normalization to sharpen the MPn signal or to show per-biomass activity. The recent work shows that the draft's good fits came from it:
  - The simulations fixed each model's biomass to relative abundance ÷ DNA, and two depth samples have almost no DNA (0.16 µg g⁻¹).
  - On a dimensionally consistent areal basis, the MPn models supply about half the measured methane (56 %). They do not track the differences among the six unrestored cores (ρ = −0.20), and neither do the betaine models (ρ = −0.09).
  - The revision rewrites the modeling Methods and Results on this basis and redraws Figure 7.
- **Specificity.** Susannah's 3 November reply on the betaine-evidence thread is confirmed. Everything abundant in the unrestored salterns correlates with methane, and no guild correlates more strongly than random sets of MAGs unless both axes are divided by DNA.
- **Other answers from the new analyses:**
  - **TMA:** no MAG has CutC. The TMA-forming enzymes present are betaine reductase and TMAO reductase.
  - **Substrate-addition simulation:** only three models make methane, and together they supply 19 % of the measured methane, not "more than enough".
  - **"Most reasonable media":** the betaine models need O₂.
  - **Network sub-communities:** not beyond chance.
  - ***Halomonas* sp. D:** this is isolate I1. No model of it is in the repository, and no MAG matches it.
- **Open for the authors (36):**
  - **Methods details** from Susannah's 3 November questions (substrate additions, culturing, metabolomics), best answered by Cliff and David.
  - **Decisions:** the central or graphical figure, naming the modeling method, and moving the lactate detail to the SI.
  - **References and housekeeping:**
    - KBase narratives to cite and make public;
    - Table S4;
    - READMEs and a Zenodo release;
    - abbreviations;
    - the bulk-density citation;
    - the Figure S7 title and wording;
    - the Zhou et al. year.
- **Changes I made to my own revision text because of these comments:**
  - The abstract now calls methylphosphonate "an important source of methane", not a "methanogenic substrate" (Susannah, 7 May).
  - The new MAG Methods paragraph cites Zhou et al. (2022) (Susannah, 7 May).
  - All new modeling text says "DNA concentration" (Clifton, 13 April).
- **Still in the authors' text** (accepted view of the revision):
  - **"Methanogenic substrate":** the Introduction's last paragraph uses it for methylphosphonate four times: "an important methanogenic substrate", "methanogenic substrate potential", "the hypothesized methanogenic substrates" and "the primary methanogenic substrate". The two uses in the Discussion (substrate additions, about the controls) are generic and fine.
  - **"Zhou et al. (2021)":** it appears six times in four places: Methods (Metagenomics), Discussion (Substrate additions, three times), Conclusions, and the Figure S9 caption. Reference 34 also says 2021.
  - **"DNA density":** the Figure S7 caption still uses it.
  - **"Trace consumption":** the abstract keeps the authors' phrase "trace consumption of methylphosphonic acid". Because the models release one methane per methylphosphonate, the methylphosphonate consumed equals the methane produced. The amount is small only relative to each model's phosphorus uptake, so the word may mislead.
- **Typo in one of Susannah's 3 November suggestions:** replacing "Measuring " with "The presence of" leaves "The presence oflactate" (Methods, lactate paragraph). The mirror reproduces it as suggested.

## Front matter

- **Open.** A. Freiburger, 2026-04-08, on the journal list ("Biogeochemistry - Springer Nature (first choice)").
  - **Comment:** "what is our central figure? … create a cartoon that depicts this system and our findings". Clifton linked an example paper.
  - **Context:** a draft graphical overview exists (`clauded_insights/conceptual_figure1.svg`, July 2026), but it is not in the manuscript. After the modeling correction, a cartoon should rest on the causal evidence: the substrate additions (Figure 4) and the *Halomonas* culturing (Figure 6). It should present the modeling as showing that MPn consumption can supply methane of the measured magnitude, not that it explains the differences among cores.

## Introduction

- **Open.** Susannah Tringe, 2026-05-07, on reference "(15)" after "…halophilic or halotolerant methanogenic taxa adapted for these resource conditions thrive, even in the presence of high sulfate".
  - **Comment:** "Need to confirm that this reference supports my statement about sulfate."
  - **Context:**
    - Reference 15 is McGenity & Sorokin (2018), "Methanogens and Methanogenesis in Hypersaline Environments", a book chapter. I have not checked its text.
    - The usual argument in this literature is that methylotrophic methanogens persist at high sulfate because sulfate reducers don't compete for methylated substrates.
    - The data here agree: across the 12 cores, methane rises with sulfate (ρ = 0.97; `clauded_insights/REASSESSMENT.md`).
- **Partly addressed.** Susannah Tringe, 2026-05-07, on "methanogens" in the same paragraph.
  - **Comment:** the field uses "methanogen" only for archaeal methanogens, not, for example, for methylphosphonate degraders.
  - **Context:**
    - The sentence it is on is already consistent ("These salt-adapted methanogens are archaeal methylotrophs").
    - The revision uses "methanogen" only for the *mcrA* MAG (*Methanolobus*), uses "MPn consumers" for C–P lyase carriers, and its abstract says "source of methane".
    - The Introduction's last paragraph still calls MPn a methanogenic substrate four times (see Summary).
- **Answered in its thread.** A. Freiburger, 2026-04-09, on the trehalose sentence.
  - **Comment:** "is this trivial or is trehalose a part of the discussion later?"
  - **Context:** Clifton replied that trehalose comes up in the metabolite results. Can be resolved.
- **Answered by new analysis.** A. Freiburger, 2026-04-09, on "…supporting that methylphosphonate is the methanogenic substrate…" (Introduction, last paragraph), with replies from Clifton (2026-04-13) and Susannah (2026-11-03).
  - **Comment:** "what was the evidence for betaine?" Susannah replied: "anything that is more abundant in the unrestored sites ends up being correlated with methane."
  - **Context:** the R232 analysis confirms her point.
    - **Without DNA normalization** (relative abundance against the measured flux, 12 cores), these all correlate with methane:
      - betaine reducers (*grdH*, ρ = 0.68);
      - MPn consumers (*phnJ*, 0.76);
      - the methanogen MAG (0.78);
      - 254 of 310 MAGs (FDR q < 0.05, either sign).

      No guild correlates more strongly than random sets of the same number of MAGs (Figure 9, Table 5).
    - **Metabolites:** measured betaine (ρ = 0.70) and TMA (0.86) correlate with methane too (10 cores).
    - **With DNA normalization:** MPn consumers stand out from random sets only when both axes are divided by DNA (Figure S8).
    - **Betaine models:** they don't track methane (ρ = −0.09). They are a weak test, though, because most of their fluxes sit at the assumed uptake bound.
- **Open.** A. Freiburger, 2026-04-09, on "…revealing that restoration efforts should…", with a reply from Susannah (2026-05-07).
  - **Comment:** "how would restoration efforts change knowing that methylphosphonate is the culprit?"
  - **Susannah's reply:** "I'm not sure I would state this so boldly…". MPn breakdown is aerobic and archaeal methanogenesis is anaerobic, so changes that add oxygen may not directly affect MPn breakdown.
  - **Context:**
    - The revision supports her caution: the modeling doesn't track the cores, and guild specificity appears only on the DNA-normalized basis.
    - Her oxygen-based framing is the defensible one.
    - The sentence is unfinished in the draft ("…restoration efforts should …").
    - 65 of 90 betaine-model simulations need O₂, so the anaerobic betaine route was simulated in aerobic models.
- **Open, resolved by her own suggestion.** Susannah Tringe, 2026-05-07, on "minerals".
  - **Comment:** "Methylphosphonates are also organic molecules, not minerals."
  - **Context:** her pending suggestion of the same day rewrites the sentence and removes "methylated minerals such as". Accepting it resolves the comment.
- **Open.** Susannah Tringe, 2026-05-07, on "among the most abundant of all studied".
  - **Comment:** "Most abundant genes in this environment overall? Or most abundant genes related to methane metabolism? Or more abundant relative to other environments?"
  - **Context:**
    - The wording ("…of all studied saline-hypersaline environments (34, 44)") reads as a comparison across environments. Saying so explicitly would answer her.
    - At MAG level, the 22 *phnJ* carriers hold 7.3 % of reads in the unrestored cores, against 0.10 % in the restored and 0.20 % in the reference sites (`correlational/r232/guild_membership_2020_vs_r232.md`).

## Methods

- **Open.** A. Freiburger, 2025-09-19, on the "Methods" heading, with replies from Clifton (2025-09-19) and A. Freiburger (2026-04-09).
  - **Comment:** a graphical figure that depicts the system, the sampling, and how these fed the substrate experiments and the models.
  - **Context:** same as the central-figure thread; a draft exists.
- **Outdated: already in the text.** A. Freiburger, 2026-04-21, Field sites and sampling.
  - **Comment:** "mention the triplicate sediment core sampling at each of these sites, in addition to deep measurements below".
  - **Context:** the anchored sentence already says "three 15 cm deep sediment cores were taken … from each of" the four sites, and the next paragraph gives the D1/D2 split. Can be resolved.
- **Open.** Susannah Tringe, 2026-11-03, Metagenomics, on the compound list ("trimethylamine, methylamine, methanol, glycine betaine, …").
  - **Comment:** "We might want to introduce the abbreviations for these compounds here."
  - **Context:** her suggestion in the first Results paragraph spells out DMS, DMSP and MMA. Defining all abbreviations here, at first use, would let later sections use them. Also harmonize "GB" (Figure 1) with "Bet." (Figure S5).
- **Open (methods detail).** Susannah Tringe, 2026-11-03, Metabolomics:
  - on "10 ml of sample": "Porewater?";
  - on "subset of 1-4 samples": "Do we not know how many? Or 1-4 samples per sample type?"
- **Answered by new analysis.** A. Freiburger, 2025-07-14, on the "Metabolic Modeling" heading, with his reply of 2025-08-10.
  - **Comment:** "discuss the TMA hypothesis as well". The reply adds that more methane after TMA addition "doesn't mean that this is the methane source".
  - **Context:**
    - No MAG carries CutC (choline TMA-lyase), and the four "TMA-route candidates" carry no TMA-forming enzyme.
    - The TMA-forming enzymes that are annotated are betaine reductase (*grdH*, 9 MAGs) and TMAO reductase (7 MAGs).
    - The only methanogen MAG, *Methanolobus*, has the full methylamine (TMA, DMA, MMA) methanogenesis machinery.
    - Of the seven measured metabolites tested, TMA correlates most strongly with methane (ρ = 0.86, 10 cores).
    - None of this is in the text yet. A short Discussion paragraph could draw on `correlational/r232/archaeal_mepn_and_tma_candidates.md`.
- **Open.** Andrew Freiburger, 2026-04-21, on "these MAGs in KBase".
  - **Comment:** "provide the specific Narrative here".
  - **Context:** Clifton's inventory (`data/Files Summary.xlsx`) lists these narratives:
    - 155805: draft models of the 310 MAGs;
    - 157061: archaeal MAG models, 30;
    - 104160: *Marivita* models;
    - 139158: *Halomonas* genomes, including the two isolates;
    - 134952: *Halomonas* South Bay MAGs.

    The MAG assemblies themselves are in workspaces 77580 ("Metabolism 310 Bins") and 163264 ("Unretained MAGs"), traced from KBase provenance.

    Your account can read 155805, 77580 and 163264. On 2026-11-04 it still could not read the other four, including 139158, where you have requested access. Make them public before citing. The revision now states that the three archaeal MPn consumers were not simulated.
- **Open.** Andrew Freiburger, 2026-04-21, on "…MAG-derived models were subjected to the following series of constraints and objectives…".
  - **Comment:** "may name this method: Community Assisted Modeling (I think?; CAM)".
  - **Context:** each MAG is simulated as a single-organism model with its growth fixed to its abundance, and the fluxes are summed afterwards. There is no exchange between members, so a "community" name may overstate the method.
- **Open.** Susannah Tringe, 2026-11-03, on "minimized to the consumption of lactate".
  - **Comment:** "Not sure what this means?"
  - **Context:** suggested wording: "Each model was then solved for the minimum lactate uptake that supports its fixed growth rate." In the code, the lactate exchange is pushed toward zero, so uptake is minimized.
- **Partly addressed.** Susannah Tringe, 2026-11-03, on "biomass compound".
  - **Comment:** "What does this mean?"
  - **Context:** my edit adds that the abundance was divided by the soil DNA concentration and averaged over the two depths. The term itself stays; suggest "the biomass reaction (growth rate) of each model was fixed to …".
- **Open.** Andrew Freiburger, 2026-04-21, on the lactate paragraph.
  - **Comment:** "this is somewhat extra detail that may be better placed in an SI section".
  - **Context:** the revision only corrects its figure citation (S9 → S10).
- **Open.** Susannah Tringe, 2026-11-03, on "hypersaline" in "This enzymatic capacity is prevalent in hypersaline archaea (55)" (same paragraph).
  - **Comment:** "halophilic?"
  - **Context:** agree. Organisms are halophilic, while "hypersaline" describes the habitat; "halophilic archaea" is the usual term.
- **Addressed, with a changed conclusion.** A. Freiburger, 2026-04-12, with two of his own replies, on the regression sentence ("Linear regression between the experimental and simulated methane emissions…"). The thread quotes an older wording: "Correlations … were conducted on DNA-normalized values…".
  - **Comment:** update the correlation figure, list the cumulative correlations (both about 95 %), and state the correlation type.
  - **Context:**
    - The revision names the methods (linear regression and Spearman correlation) and redraws Figure 7 on the consistent basis.
    - The ~95 % values came from dividing both axes by DNA and are not reported.
    - Guild abundance against methane is in Figure 9 and Table 5.
- **Addressed.** A. Freiburger, 2026-04-21, on "Pearson".
  - **Comment:** "should this be a spearman?"
  - **Context:** Methods and Results now report Pearson R² and Spearman ρ, with exact permutation p values.
- **Addressed.** Andrew Freiburger, 2026-04-21, on "…on DNA-normalized values…" (same paragraph).
  - **Comment:** "possibly not DNA normalized".
  - **Context:** confirmed and done. The revision converts each member's simulated flux to an areal emission from its biomass in each core, and compares that with the measured areal emissions.
- **Open.** A. Freiburger, 2026-04-21, and Susannah Tringe, 2026-11-03, on "The soil density was estimated from other studies of this same region to be 0.7 g soil cm⁻³" (two threads).
  - **Comments:** A. Freiburger supplied a citation (a South Bay Salt Pond brownbag talk, September 2024). Susannah asked whether this is a reference or our own observation, and whether it should say "sediment density".
  - **Context:** suggest "sediment bulk density" plus that citation. The simulated magnitudes in Figure 7 scale linearly with this 0.7 g cm⁻³, and with the assumed 10 % DNA content of biomass.
- **Open (methods detail).** Susannah Tringe, 2026-11-03: four threads in Substrate additions and one in Microbial culturing.
  - **Comments:**
    - **On "10 mL of substrate and 10 mL of sediment slurry were added to sterile 75 mL glass bottles":** does this repeat the next paragraph? As written, it sounds as if they were combined right after collection.
    - **On "another subset":** "For all the anaerobic incubations?"
    - **On "Glass bottles containing substrates":** individual bottles, each with 10 mL of substrate solution?
    - **On "A subset of sediment slurry was left to sit … with aerobic atmospheric conditions":** was this also in a glass bottle, holding all the slurry for the aerobic experiments? Was it then dispensed in 10 mL aliquots after equilibration, and with a pipette?
    - **On "septum bottles containing MM/W medium":** what size bottles, and how much medium?
  - **Context:** her wording suggestions in these paragraphs ("each glass bottle", "of slurry was", "microcosms") are in the revision as her tracked changes.
- **Addressed.** A. Freiburger, 2025-09-19, on the "Analyses" heading.
  - **Comment:** "I should include the correlational methods here".
  - **Context:**
    - The revision replaces P52 with the correlation design and adds a co-occurrence methods paragraph.
    - Susannah's suggestion renames the section "Omics data analyses and statistical methods".
- **Addressed.** Andrew Freiburger, 2026-04-22, on P52.
  - **Comment:** "elaborate the methods employed here".
  - **Context:** P52 is replaced.

## Results

- **Open.** Andrew Freiburger, 2025-09-18, Metabolomics, on "Metabolite profiles from all of the unrestored saltern samples clustered together…".
  - **Comment:** "consistent with the original community modeling". Clifton asked whether this means Andrew's work or Zhou et al.
  - **Context:** Andrew to clarify; no new analysis bears on this point.
- **Addressed.** A. Freiburger, 2026-04-09, on the "Metabolic modeling" heading.
  - **Comment:** "needs updating".
  - **Context:** the section is rewritten on the consistent areal basis (Figures 7 and 9, Table 5).
- **Answered by new analysis.** A. Freiburger, 2025-08-08, on the same heading, with his reply of 2025-08-10.
  - **Comment:** "add a section for the substrate addition experiments, where the models generated more than enough methane from MePn … only a few members were able to grow/generate methane". The reply notes that the addition sediment came from R2, with an average R2 emission of 1.3 µmol hr⁻¹.
  - **Context:**
    - In the stored substrate-addition simulation, only three Rhodobacteraceae models take up MPn: *Roseovarius*.9, *Roseovarius*.11 and *Albimonas*.1. That fits "only a few members", and the revision now says so.
    - On the consistent basis, that run supplies 19 % of the measured field methane, and only 5 % in the R2 cores. It is infeasible in cores R2_B and R2_C (`clauded_insights/modelling_consistent_basis_r232.json`, set "addition_no_ratio").
    - So the stored runs don't support "more than enough". The comparison with the 1.3 µmol hr⁻¹ was not redone.
- **Addressed, with a changed conclusion.** Andrew Freiburger, 2026-04-22, on the modeling equations (P69).
  - **Comment:** "these need to be updated with the results from the new DNA density weighted abundances".
  - **Context:**
    - The equations and statistics are replaced by fits on the consistent areal basis. For the MPn models, ρ = −0.20, R² = 0.08, and they supply 56 % of the measured methane.
    - The simulations used the older DNA-normalized abundances (the mean over depths of relative abundance ÷ DNA). Each model's methane is proportional to its abundance, so rescaling to each member's biomass in each core is exact and needed no re-simulation.
- **Open, with data available.** A. Freiburger, 2025-09-19, on *Roseovarius* and *Albimonas*.
  - **Comment:** include their total abundances in each sample. Clifton asked whether this means MAG abundances or IMG genus abundances.
  - **Context:**
    - Per-MAG abundances by core are in Table S5 and in `data/Cliff_310MAG_relabund.tsv`.
    - Suggest MAG abundances, consistent with the revision.
    - The three models whose ratio constraint does not bind (*Roseovarius*.9 and .11, *Albimonas*.1) produce 72 % of the simulated methane.
- **Addressed.** A. Freiburger, 2026-04-13, on the R² values in P78.
  - **Comment:** "include these figures in Figure 7, or possibly as their own figure". His reply: "relative abundance is better to be discussed".
  - **Context:** P78 is replaced. The new Figure 9 plots relative abundance against methane, Figure S8 gives the DNA-normalized version, and Table 5 summarizes both.

## Discussion

- **Answered by new analysis.** A. Freiburger, 2025-09-18, Culturing, on the *Halomonas* isolate sentence, with replies from Clifton.
  - **Comment:** which isolate, and did its model consume MPn? Clifton: "*Halomonas* sp. D … Do you all have a model for this organism? … I tried to enlist Jose's help modeling it several years ago."
  - **Context:**
    - *Halomonas* sp. D is isolate I1: GTDB R232 *Halomonas* sp032334135, NCBI GCA_032334135.1. It has a complete phnC–P cluster.
    - No model of it is in this repository or in narrative 155805, none of whose 1,482 objects is named for *Halomonas* or the isolates. Its genome is in narrative 139158, which is still not readable.
    - No MAG matches either isolate at species level, and neither Halomonadaceae MAG carries C–P lyase genes. The community models therefore contain no *Halomonas* MPn consumer.
    - A model could be built from the isolate genome.
- **Open.** A. Freiburger, 2025-09-21, Metabolic modeling, on "The binned reads also match to known and expected organisms, which suggests that the main actors are included…".
  - **Comment:** "verify this logic".
  - **Context:**
    - MAGs capture about 48 % of reads in the unrestored cores, but only 11–12 % in SF2 and R2A; unbinned reads make up 52, 88 and 89 %.
    - The methanogen MAG is at most 0.3 % of reads.
    - The logic holds reasonably for the unrestored salterns, but not for the restored and reference sites.

## Conclusions

- **Open; the revision argues against it.** A. Freiburger, 2025-09-19, on the "Conclusions" heading.
  - **Comment:** "try to make a compressive, definitely, statement regarding the methane source".
  - **Context:**
    - The modeling reproduces the magnitude of methane but not its pattern among cores.
    - Guild abundance is specific only on the DNA-normalized basis.
    - My comment on the Conclusions suggests "the most consistent candidate source", resting on Figures 4 and 6.
- **Open.** A. Freiburger, 2025-09-18, on "…metabolizing methylated compounds and acetoclastic and hydrogenotrophic archaeal methanogenesis".
  - **Comment:** "elaborate what evidence exists for these hypotheses".
  - **Context:**
    - The dereplicated MAGs include a single *mcrA* MAG, *Methanolobus* (methylotrophic), found only in the unrestored cores.
    - No acetoclastic or hydrogenotrophic methanogen MAG was recovered, so evidence for those routes would have to come from the IMG gene data (Figure 1) or the incubations.

## Data Availability

- **Open.** A. Freiburger, 2025-07-14, with a reply from Clifton.
  - **Comment:** the GitHub repository "needs to have READMEs". Clifton: make a release and link it to Zenodo.
  - **Context:**
    - The MicrobiomeNotebooks repository has READMEs for the R232 work (`correlational/r232/`, `clauded_insights/`, `MAG_QC/`), and the revision adds the repository to Data Availability.
    - The draft cites `HenryLabResearch/Salterns`, while the revision adds `ModelSEED/MicrobiomeNotebooks`. Reconcile the two before the Zenodo release.
- **Open.** Susannah Tringe, 2026-05-01, on "Table S4".
  - **Comment:** "Where can I find this table? Also, are there narratives and/or data in KBase that we should cite?"
  - **Context:**
    - Only the caption of Table S4 is in the document ("Metadata, including IMG/M and NCBI accessions, for the 24 metagenomes in this study"); the table itself is missing. The same is true of Table S3.
    - For the narratives, see the KBase entry under Methods.

## Figures and supplementary captions

- **Answered by new analysis, with a different outcome.** A. Freiburger, 2025-09-19, Figure 1 caption.
  - **Comment:** a family-level heatmap of abundance against methane, expecting the families whose models made methane to correlate strongly.
  - **Context:**
    - The guild- and MAG-level correlations (Figure 9, Table 5, `mag_methane_correlations_r232.csv`) show that most of the community correlates with methane because of the site gradient. A family heatmap would mostly show that gradient.
    - Rhodobacteraceae make up 20.6 % of binned reads in the unrestored cores (Table 3) and produce 99.6 % of the simulated methane, but their abundance correlation is not specific on the areal basis.
    - Either drop the heatmap or pair it with the random-set comparison.
- **Addressed.** Andrew Freiburger, 2026-04-22, Figure 3 caption.
  - **Comment:** "we'll perform a correlation with the MePn and betaine consuming relative abundances as well".
  - **Context:** Figure 9 and Table 5 do this, with gene-based guilds and a random-set comparison.
- **Addressed.** A. Freiburger, 2026-05-06, on the Figure 7 image, with his own reply.
  - **Comment:** clarify DNA normalization and update the axis units. The reply adds: make the lowest left axis label horizontal.
  - **Context:** Figure 7 is replaced by a single panel on the areal basis, with both axes in µmol CH₄ m⁻² d⁻¹.
- **Addressed.** Andrew Freiburger, 2026-04-22, Figure 7 caption.
  - **Comment:** "update the commentary here based on the results above".
  - **Context:** the caption is rewritten.
- **Partly addressed.** A. Freiburger, 2026-05-06, on "Figure" in the Figure 7 caption, with his own reply.
  - **Comment:** "confirm the data processing methods employed here". The reply adds: "consider also re-adding the metabolic pathways".
  - **Context:**
    - The new caption and the Methods state the biomass input (relative abundance ÷ DNA, averaged over the two depths) and the conversion to areal emissions. Both are reproduced in `clauded_insights/modelling_consistent_basis_r232.py`.
    - The pathway panels are not re-added.
- **Answered by new analysis.** A. Freiburger, 2025-08-08, Figure 7 caption.
  - **Comment:** "normalizing for the DNA means that biomass is held constant, so the universally high correlation must be because of greater per biomass metabolic activity".
  - **Context:** the premise doesn't hold.
    - Dividing both quantities by DNA gives them a shared 1/DNA factor. The two depth samples with 0.16 µg g⁻¹ make that factor dominant, and it produced the draft's R² ≈ 0.97.
    - The community-wide correlation with methane reflects the gradient between the unrestored salterns and the other sites (Table 5).
    - The revision deletes the Discussion sentence that attributed the correlation to metabolic activity. Suggest resolving.
- **Answered by new analysis.** A. Freiburger, 2026-05-06, Figure 8 caption.
  - **Comment:** "describe 1-2 biological insights about the particular subcommunities".
  - **Context:**
    - Networks fitted to permuted data are as large and as modular as the real one. The one feature beyond chance is that edges join related taxa, mostly the haloarchaea.
    - So the modules shouldn't be presented as functional sub-communities, and the revision's caption says so.
- **Addressed.** A. Freiburger, 2026-05-06, Figure 8 caption.
  - **Comment:** "Claude generated as a space-filler. subject to change".
  - **Context:** the old caption is deleted.
- **Open.** Susannah Tringe, 2026-11-03, Figure S5 caption.
  - **Comment:** "Bet. = glycine betaine … This is abbreviated as GB in Figure 1."
  - **Context:** harmonize the abbreviation.
- **Open.** Clifton, 2026-04-13, Figure S7 caption.
  - **Comment:** "Remove title and change 'DNA density' to 'DNA concentration'".
  - **Context:**
    - The revision's own text says "DNA concentration", but the Figure S7 image (with its title) and its caption are unchanged.
    - The figure's statistics are R² = 0.22 and Spearman ρ = 0.52, p = 0.08 (`modelling_numbers_r232.json`). The revised Results call the correlation weak and not significant.
- **Addressed.** Andrew Freiburger, 2026-04-22, on the Figure S8 image.
  - **Comment:** "the x-axis label should no-longer be normalized by DNA, since it is DNA weighted and has the same original abundance units".
  - **Context:**
    - The old Figure S8 is replaced.
    - The DNA-weighted relative abundance the comment describes is plotted in the new Figure 9, in relative-abundance units.
    - The new Figure S8 is DNA-normalized by design, with both axes divided by DNA as in the model inputs, and its caption says so.
- **Addressed.** Andrew Freiburger, 2026-04-22, Figure S8 caption.
  - **Comment:** "how were these organisms originally identified?"
  - **Context:** the caption is replaced. The guilds are defined by marker gene (*phnJ*, *grdH*) in the new Methods subsection and in Table 4.
- **Addressed.** Andrew Freiburger, 2026-04-22, on the stray caption fragment "…ne of the highest correlations…".
  - **Comment:** "fragment that could be integrated previously in the paper".
  - **Context:** the revision deletes it.
- **Answered by new analysis.** A. Freiburger, 2026-04-11, Figure S10 caption.
  - **Comment:** "are these actually the most reasonable media?"
  - **Context:**
    - For the betaine models, no. They were gap-filled and simulated on the aerobic lactate medium, and 65 of 90 betaine simulations are infeasible without O₂.
    - So an anaerobic, Stickland-type enzyme was simulated in aerobic models (`correlational/betaine_grdH_r232/`).

## References

- **Open; the revision's new text follows it.** Susannah Tringe, 2026-05-07, on the year of reference 34 (Zhou et al., ISME J 16:284–295).
  - **Comment:** "This was initially released in 2021 but official publication date appears to be January 2022."
  - **Context:**
    - Changing the year in the reference manager fixes reference 34.
    - The "Zhou et al. (year)" mentions in the text are typed separately from the numbered citation:
      - Susannah's pending suggestions change three of them in the Methods to 2022;
      - my new MAG paragraph cites 2022;
      - six "Zhou et al. (2021)" remain (see Summary).

## Threads whose anchor text is no longer in the document (11)

Google keeps these threads open, but the text they were attached to has been deleted, so its exports drop them. They are listed by date.

- **Open.** Clifton, 2022-08-24, on "85)".
  - **Comment:** "add other refs".
  - **Context:** probably the "(85)" in "only 3% of the TMA was converted to CH4, which is lower than other studies (85)" (Discussion, Substrate additions). Reference 85 is Oremland et al. (1989).
- **Outdated.** Clifton, 2023-10-20, on "Figure S4".
  - **Comment:** "Figure place holder for new metabolic modeling figure."
  - **Context:** Figure S4 is now the metabolite–gene heatmap, and the modeling figure is Figure 7. Can be resolved.
- **Open, with context.** A. Freiburger, 2025-07-13, with his reply of 2025-08-06, on "Betaine pathway".
  - **Comment:** "elaborate why Betaine is consumed at all". The reply adds: reference the first paper, which found correlations with betaine-consuming genes.
  - **Context:**
    - The nine betaine reducers carry the betaine-specific reductase subunit (*grdH*). Betaine reduction is a Stickland-type electron sink that releases TMA.
    - The carriers are Bacillota (*Halanaerobium*, M55B127, BM714, DATFJG01 and UBA8670), found mainly in the unrestored cores.
- **Answered by new analysis.** Andrew Freiburger, 2025-09-12, on ".".
  - **Comment:** "only a trace 0.0297 umol/hr of methylphosphonate must be consumed to produce the observed methane flux".
  - **Context:**
    - In every model, one methane is released per methylphosphonate consumed: exactly 1 : 1, which corrects the draft's "1.4 : 1". The MPn needed therefore equals the methane emitted.
    - On the consistent basis, the models supply 56 % of the measured flux.
    - I could not trace the 0.0297 figure. See also the note on "trace consumption" in the Summary.
- **Open.** Andrew Freiburger, 2025-09-12, on "…nerate the observed methane emissions".
  - **Comment:** "cite the specific figures or tables that contain these revelations".
  - **Context:** the modeling claims are rewritten. Cite Figure 7 for the modeling and Table 5 for abundance.
- **Answered by new analysis.** A. Freiburger, 2025-09-19, on "The correlational compar…", with his own reply.
  - **Comment:** "discuss the betaine pathways". The reply adds: "possibly discuss recent progress in metabolic modeling".
  - **Context:** as built, the betaine models cannot test the betaine route:
    - the betaine reductase reaction was added with no gene link;
    - 61 % of their fluxes sit at the uptake bound;
    - the simulations need O₂.
- **Outdated: already cited.** A. Freiburger, 2026-04-08, on "experimental flux of glycine betaine".
  - **Comment:** a link to PMC111389.
  - **Context:**
    - The link is reference 54 (Lai, Hong & Gunsalus 2000, J Bacteriol 182:5020–5024), which Methods P36 already cites for the betaine consumption rate of 0.1.
    - The paper reports a betaine transport Vmax of 8 nmol min⁻¹ (mg protein)⁻¹ in *Methanohalophilus portucalensis*. If protein is half the dry weight, that is about 0.25 mmol g DW⁻¹ h⁻¹, so 0.1 is of the same order.
    - How 0.1 was derived is not stated. That bound sets most of the betaine fluxes, so the derivation is worth adding.
- **Addressed.** A. Freiburger, 2026-04-08, on "experimental data".
  - **Comment:** "which specific data type?"
  - **Context:**
    - The phrase survives in Methods P36 ("0.07 yielded the best fit to the experimental data"), and the revision replaces it with the minimum MPn : Pi uptake ratio.
    - The "best fit" held only when both axes were divided by DNA.
- **Open.** A. Freiburger, 2026-04-09, on "microbial taxonomic and functional compositions".
  - **Comment:** "give notable examples of each".
  - **Context:** examples are now available. The unrestored salterns are dominated by Rhodobacteraceae (20.6 % of binned reads, mostly *Roseovarius*), Balneolaceae, Wenzhouxiangellaceae and haloarchaea (Haloferacaceae, 5.9 %; Table 3).
- **Open.** Andrew Freiburger, 2026-04-22, on the pathway descriptions ("…nate utilization pathway, which is a source …").
  - **Comment:** "add the beautified pathways above".
  - **Context:** the same pathway text is quoted in the Figure 7 data-processing thread. Whether to include a pathway figure is an author decision.
- **Answered by new analysis.** A. Freiburger, 2026-04-27, on "relative abundance of at least 1e-5 in 3 of the 4 experimental sites (R1, R2, R2A, SF2)".
  - **Comment:** "this is without adjusting for the non-binned samples".
  - **Context:**
    - Relative abundances are fractions of all reads, unbinned included (`correlational/r232/prefilter_sensitivity/abund_REPORT.md`).
    - That ≥ 3-of-4-sites filter belongs to the site-level correlation heatmap. It removes every MAG found only in the unrestored salterns, including *Methanolobus* and 8 of the 9 *grdH* MAGs.
    - The network (Figure 8) and Figure S11 instead keep MAGs present in at least 6 of the 12 cores.
