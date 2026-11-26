# Curation of the *Methanolobus* methanogen model (R1_A_D1__concoct.bin.007)

2026-11-24. Species *Methanolobus* sp020055655 (GTDB R232), representative MAG `R1_A_D1__concoct.bin.007` (CheckM2 99.91 % complete, 0.40 % contaminated). Nothing was committed, nothing was written to KBase, and neither the templates nor ModelSEEDpy were changed.

| File | What |
|---|---|
| `R1_A_D1__concoct.bin.007.json` | arm-A model of record (md5 `e8fd1042…`, as in `../MANIFEST.md5`), unchanged |
| **`R1_A_D1__concoct.bin.007__curated.json`** | curated model, md5 `893a75cdb6d83f6c79617fe9d27f1183` |
| `../reports/methanogen/R1_A_D1__concoct.bin.007__curated.validation.json` | every change (equation, bounds and GPR before and after) and every number in this file, for both the curated and the arm-A model |
| `../../pipeline/curate_methanogen.py` | the curation and validation (md5 `ce9d2200…`); run `~/Documents/py_venv/bin/python -I pipeline/curate_methanogen.py` (≈20 s, local, Gurobi) |

- Gene IDs are written `NODE_<n>_<gene>`; the model uses the full Prodigal IDs (`NODE_296_14` = `NODE_296_length_66342_cov_4.286295_14`).
- "K" values are the KOfam best KO with score/threshold (`../../annotations/kofam/`). This genome has no Sec read-through rows.
- Literature marked † was checked in this session (abstract or full text). Literature marked ‡ is cited from memory and not re-checked; see References.

---

## 1. Summary

- **Arm A grows as an acetogen because it cannot reduce methanophenazine (MP) from methyl oxidation.**
  - Mer/Mtd make F420H2 on coenzyme F420 (cpd00792). The template Fpo (rxn24614) runs on a different, disconnected compound pair (cpd28067/cpd27742), and only in reverse.
  - ATP correction removed Rnf's forward direction.
  - Without reduced MP, heterodisulfide cannot be recycled and Mcr carries no flux.
  - Every methyl group therefore ends as acetate: half through CODH/ACS (gene-supported), fed by a free, uphill F420H2 → ferredoxin step (rxn40355); half through Fae/Hps plus a gap-filled phosphoketolase/Pta. ATP comes from the ADP-forming acetyl-CoA synthetase (AcdAB): 0.25 ATP per methyl.
- **Arms B–D cannot grow on methylamines because they have no methylamine transport at all.**
  - Arm A's three transporters came only from its ATP-correction core gap-fill.
  - Their methanol growth is itself an artefact: the lumped reaction rxn10568 entered through the MtaA KO alone. They make CH4 by running the gene-less Frh backwards and exporting formate.
- **This genome cannot use methanol.** No MtaB (K04480) and no MtaC were found in six independent assemblies of the species, including the GTDB reference genome, by a 6-frame HMM search. The curated model does not grow on methanol.
- **The curated model is a methylamine methylotroph and makes CH4 only.**
  - It grows on MMA, DMA and TMA with acetate 0, both at maximum growth and at the substrate-minimising optimum.
  - At zero growth it reproduces the textbook stoichiometries exactly.
  - It gets no ATP from carbon-free media and does not grow on H2/CO2.
- **k = 126.5 mmol CH4 per g biomass formed for every methylamine and for the mixed medium** (7.91 g biomass per mol CH4). Across literature ion couplings it ranges 114–159 mmol g⁻¹. That is 1.1–1.55× the yields measured for *Methanosarcina barkeri* (5.1–7.3 g mol⁻¹) and 4.4× one *Methanococcoides* value (1.8 g mol⁻¹).
- ***Roseovarius*.2's "methane monooxygenase" is a propane 2-monooxygenase** (PrmA, SDIMO group 5: 65 % identity to PrmA, 32–34 % to MmoX). Block rxn00844 (§8).

---

## 2. Diagnosis

### 2.1 What arm A does on the methyl medium (FBA, maximum growth 0.418 h⁻¹)

| Flux (mmol gDW⁻¹ h⁻¹) | Reaction | What it does | Gene support |
|---:|---|---|---|
| 149.5 / 49.5 | rxn24608 + rxn24609 / rxn24610 + rxn24611 | MMA / DMA → methyl-CoM | yes (NODE_296 cluster) |
| −199 | rxn03020 (Mtr, reverse) | methyl-CoM → CH3-H4MPT, 2 Na⁺ in per methyl | yes (mtrA–H) |
| −99.7 | rxn15962 (ACS) | CH3-H4MPT + CO + CoA → acetyl-CoA | yes (cdhC) |
| −99.3 | rxn40505 (CODH) | CO2 + 2 Fd(red) → CO | yes (cdhAB) |
| +99.3 | rxn40355 | **F420H2 → 2 Fd(red), with no energy input (uphill)** | FpoF only |
| −99.3 | rxn03085 (Mer, reverse) | CH3-H4MPT → CH2=H4MPT + F420H2 | yes |
| −99.3 / 93.3 | rxn07847 (Fae, reverse) / rxn03643 + rxn03644 (Hps/Phi) | CH2=H4MPT → formaldehyde → F6P | yes (fae-hps) |
| 47.3 / −47.3 | rxn00548 / rxn00173 | **phosphoketolase F6P → acetyl-P; Pta acetyl-P → acetyl-CoA** | **none (gap-filled)** |
| −139.9 | rxn00172 (AcdAB) | acetyl-CoA + ADP + Pi → acetate + ATP | yes (NODE_1953_12, K24012) |
| −101 | rxn08173 (ATP synthase, hydrolysing) | pumps H⁺ out; the H⁺/Na⁺ antiporter then exports the Na⁺ that Mtr brought in | yes |
| 140.5 | acetate export | | |
| 0 | rxn03127 (Mcr), rxn03126 (Hdr) | | |

- Net stoichiometry: 4 CH3 + 2 CO2 → 3 acetate.
- With acetate export blocked, growth is 0 on every medium.
- The methane gap-fill (rxn08518, formate hydrogenlyase) made CH4 reachable only apart from growth, through formaldehyde → formate → H2 → Vht → Hdr.

**Every route to reduced methanophenazine (MPH2) is closed in arm A**, so Hdr cannot regenerate CoM-SH + CoB-SH after Mcr:

| Route to MPH2 | Arm A | Cause |
|---|---|---|
| F420H2 → Fpo | rxn24614, `<` only | (i) the template direction is `<`; (ii) **it runs on "Reduced-/Oxidized-Factor-F420" (cpd28067/cpd27742), a compound pair made only by the H4SPT copies of Mer/Mtd (rxn20187/rxn20214), which never connect to methyl-CoM**; (iii) no H⁺ is translocated. This is why the repair's `fpo_freed_test`, which changed only the direction, could not work |
| Fd(red) → Rnf | rxa45615, `<` only | the template has `=`; ATP correction's expansion filter removed `>` (`expansion_filtered` in the build report) |
| H2 → Vht | rxn15961 | no H2 in the methyl medium |

### 2.2 Arms B–D on methylamines

| | Arm A | Arms B, C, D |
|---|---|---|
| ATP-test media chosen by ATP correction | H2.CO2, Methanamine.H2, Dimethylamine.H2, Trimethylamine.H2 | B: H2.CO2; C, D: H2.CO2, Methanol.H2 |
| ATP-correction core gap-fill | rxn05759, rxn06299, and exchange + transport for MMA (rxn31759), DMA (rxn33011) and TMA (rxn09318) | rxn05759; in C and D also rxn24606 (MtaB/MtaC, no gene) |
| MMA / DMA / TMA transport | present (no GPR) | **absent**, so growth and CH4 are 0 on each methylamine |
| Methanol route | none | rxn10568, a lumped methanol:CoM methyltransferase (`gapfilling` type, no roles), from K14080 (MtaA-like) hits on NODE_485_13, NODE_1297_10, NODE_2869_3 (+ eggNOG NODE_427_14) |
| How CH4 is made on methanol | – | methanol → CH3-CoM (rxn10568); one in three methyls is oxidised; F420H2 → **gene-less Frh run backwards** (rxn06299) → H2 → Vht → MPH2 → Hdr; formyl-MFR is **hydrolysed to formate** (rxn07849) and exported instead of oxidised |

The methyltransferase reactions themselves (rxn24608–rxn24613) are present and gene-associated in every arm; only uptake is missing.

### 2.3 The methylotrophic pathway, step by step

Status codes: **ok** = present, gene-associated and correct; **dir** = wrong direction; **cof** = wrong cofactor or compound pool; **GPR** = genes present but missing from or wrong in the association; **absent** = no gene in the genome.

| Step | Genes (locus; KOfam score/threshold; RAST) | Arm-A reaction | Status | Curated |
|---|---|---|---|---|
| **Methanol** MtaB, MtaC | **none**: no K04480 hit in any frame; the K14081 hits (≤186 vs threshold 301) are the methylamine corrinoid proteins | none (no transport either) | **absent** | not added (§2.4) |
| MtaA-like (CmuA family) | NODE_485_13 (K14080 464/383; RAST "Methylcobalamin:CoM methyltransferase MA0779"), NODE_1297_10 (414/383), NODE_2869_3 (403/383), NODE_427_14 (eggNOG K14080) | none (template rxn24607 needs the MtaC role) | present, no partner | not added: MtaA cannot demethylate methanol without MtaB/MtaC |
| **MMA**: MtmB (Pyl) + MtmC | MtmB = NODE_296_14 + _15 and NODE_296_29 + _30. Each pair is one 458-aa protein that Prodigal split at the in-frame UAG (codon 202); K16176 304/260 and 370/260; RAST "… @ pyrrolysine-containing". MtmC = NODE_296_13, _28 (K16177 358/334) | rxn24608 `>` | **GPR** (pieces OR-ed) | GPR: pieces AND-ed |
| **DMA**: MtbB (Pyl) + MtbC | MtbB = NODE_296_26 + _27 (487 aa, UAG codon 376; K16178 667/497); MtbC NODE_296_25 (K16179 333/276) | rxn24610 `=` | **GPR** | GPR fixed |
| **TMA**: MttB (Pyl) + MttC | MttB = NODE_296_19 + _20 (496 aa, UAG codon 334; K14083 305/219); MttC NODE_296_21 (K14084 316/282) | rxn24612 `=` | **GPR** | GPR fixed |
| MtbA (methylcobamide:CoM) | NODE_296_33 (K14082 499/432), NODE_1297_13 (454/432) | rxn24609, rxn24611, rxn24613 | ok | – |
| Pyrrolysine | pylB NODE_296_11 (K16180), pylC _10 (K16181), pylD _9 (K16182), pylS _12 (K11627). tRNAscan did not call pylT; a 447-bp gap lies upstream of pylS | not modelled | – | – |
| RamA (reactivates corrinoids) | NODE_296_7, _8 (K16183 948/919, 1003/919) | not modelled | – | – |
| Non-Pyl MtmB homologue | NODE_296_31 (457 aa, K16176 568/260, no in-frame UAG) | – | substrate unknown | – |
| **Mcr** | mcrBDCGA NODE_601_17–21 (K00399 983/741, K00401, K00402, K03421, K03422) | rxn03127 `=` | **dir** (could consume CH4) | `>` only |
| **Mtr** (runs in reverse in the oxidative branch) | mtrEDCBAFGH NODE_85_97–104 (K00577–K00584, all above threshold) | rxn03020 `=`, 2 Na⁺ | ok | – |
| **Mer** | NODE_243_11 (K00320 419/220) | rxn03085 `=` | ok | – |
| **Mtd** (F420-dependent) | NODE_209_90 (K00319 417/331) | rxn03079 `=` | ok | – |
| **Mch** | NODE_2689_6 (K01499 418/287) | rxn02480 `=` | ok | – |
| **Ftr** | NODE_1502_7, NODE_2990_1 (K00672) | rxn02431 `<` (oxidative direction) | ok | – |
| **Fmd/Fwd** | Mo-Fmd: fmdA NODE_4070_15 (K00200 745/244), fmdB NODE_404_55, fmdC NODE_4070_16, fmdD NODE_404_56, fmdE NODE_2084_2, fmdF NODE_4070_14. W-Fwd: fwdB NODE_192_9, fwdD NODE_192_10, fwdF NODE_192_2, fwdG NODE_192_8 | rxn11938, **NAD⁺-linked**. The ferredoxin form rxn46184 is also present, but on the generic formyl-MFR pool, reachable only in the CO2-reducing direction (via rxn17445) | **cof** | MCUR_FWD (ferredoxin, specific pool) |
| **Fpo** | fpoA–D, H–N NODE_863_27–38; fpoMNO NODE_227_1–3; fpoF NODE_243_12 (K22158–K22170, all above threshold) | rxn24614, `<`, generic F420 pool, no H⁺. FpoK/FpoL (named "NADH-ubiquinone oxidoreductase chain K/L" by RAST) also carried **rxn10122, a ubiquinone complex I with 3.5 H⁺** | **dir + cof**; **GPR** | MCUR_FPO; rxn10122 removed |
| **HdrDE** (membrane) | hdrD NODE_664_33 (K08264 530/374), NODE_1064_33 (381/374); hdrE NODE_664_32 (K08265 376/212) | **not mapped**; rxn03126 carried hdrABC instead, with no H⁺ translocation | **GPR** | MCUR_HDRDE (HdrDE genes, 2 H⁺) |
| HdrA2B2C2 (soluble, electron-bifurcating) | hdrA NODE_192_3 (K03388 860/563), hdrB NODE_209_30 (K03389), hdrC NODE_209_29 (K03390) | mis-mapped to rxn03126 | **GPR** | genes left unassigned; reaction not added (§3.3) |
| **Vht** | vhtACG NODE_389_37–39 (K14068–K14070), vhtD NODE_389_36 | rxn15961 `=`, 2 H⁺ | ok | – |
| **Rnf** + MmcA | rnfB NODE_835_24, rnfA _25, rnfE _26, rnfG _27, rnfD _28 (eggNOG only), rnfC _29; multihaem cytochrome c NODE_835_30 (Pfam Cytochrom_NNT), the *M. acetivorans* rnf-mmcA arrangement | rxa45615, `<` only (after ATP correction), GPR OR-ed | **dir + GPR** | `>` only, GPR AND-ed with MmcA, ion coupling set (§4) |
| Ech | none (K14086–K14091) | – | absent | – |
| Frh | no FrhA/FrhG; FrhB-like single subunits NODE_209_24, NODE_1526_2 (K00441, Pfam FrhB_FdhB) | rxn06299 (no GPR, core gap-fill) | **absent** | removed |
| Fno (F420:NADP) | NODE_209_20 (RAST "NADPH-dependent F420 reductase", eggNOG K06988) | none; the template role name differs | **GPR** | rxn16596 added, `>` |
| **Na⁺/H⁺ antiporter** | Mrp: mrpA NODE_192_92, mrpCDEFG NODE_8305_4–8 (K05567–K05571). NhaP-type: NODE_835_14, NODE_1984_16 | rxn05209 `=` (NhaP genes only) | **GPR** | Mrp added to the GPR |
| **A1AO ATP synthase** | atpA NODE_2103_1, B NODE_863_1, C NODE_2103_3, D NODE_863_2, E NODE_2103_4, F NODE_2103_2, I NODE_2103_6, K NODE_2103_5 | rxn08173, GPR OR-ed | **GPR** | GPR AND-ed |
| Methylamine uptake | none assigned. NODE_296_22, an EamA-family "LicB" permease inside the mtt/mtb cluster, is a candidate MttP-type permease but unproven | rxn31759, rxn33011, rxn09318 (no GPR) | no GPR | unchanged |

### 2.4 Methanol: the genes are truly absent

A 6-frame search was run on six genomes of this species: stop-to-stop ORFs ≥ 50 aa, UAG read as Pyl, all KOfam profiles for K04480, K14081, K14080, K14083 and K16955; pyhmmer, E < 1e-5.

| Genome | Assembly | Completeness | K04480 (MtaB) | best K14081 (MtaC; threshold 301) | MttB family (K14083) |
|---|---|---|---|---|---|
| `R1_A_D1__concoct.bin.007` (the model) | rollout (metaSPAdes) | 99.91 % | **0 hits** | 186 (= MtmC/MttC) | 1 (Pyl MttB) |
| `R1_A_D2__concoct.bin.112` | rollout | 99.81 % | 0 | 184 | 1 |
| `R1_sediment__concoct.bin.101` | co-assembly (MEGAHIT) | 99.80 % | 0 | 184 | 1 |
| `R2_sediment__concoct.bin.079` | co-assembly | 99.80 % | 0 | 184 | 1 |
| 2020 original `…R1_A_D2…concoct_out.21` | 2020 pipeline | 98.37 % | 0 | – | 1 |
| **GCA_020055655.1**, the GTDB species representative | JGI MAG, same salt pond R1 (2014), 48 scaffolds † | 99.0 % / 99.9 % (CheckM / CheckM2, GTDB) † | **0** | 184 | 1. Also present: McrA 982, MtmB 709, MtbB 877 (Pyl joined) |

- No MtaB-family Pfam (PF12176) appears in eggNOG's annotations either.
- The full metagenome assemblies are on poplar and were not searched, so an unbinned *mtaCB* operon cannot be ruled out outright. Six assemblies at 98–99.9 % completeness make absence the strong reading.
- Every *Methanolobus* species whose description was checked does use methanol † (Shen et al. 2020 and others). But methylamine-only Methanosarcinaceae exist: *Methanohalophilus levihalophilus* grows on TMA and DMA, not methanol † (Katayama et al. 2014). This species would be another.
- **What-if** (on a copy, not the curated model; `methanol_what_if`): adding a lumped methanol:CoM methyltransferase gives growth 0.504 h⁻¹ on methanol and the same k (126.5), because methanol enters as methyl-CoM like the methylamines.

### 2.5 The acetogenic route and its genes

| Step | Reaction | Genes | After curation |
|---|---|---|---|
| CO2 → CO | rxn40505 (CODH) | cdhA NODE_192_71 (K00192 1261/423), cdhB NODE_192_72 | kept |
| CH3-H4MPT + CO + CoA → acetyl-CoA | rxn15962 (ACS, `<` = synthesis) | cdhC _73, cdhD _75, cdhE _76, cooC _74 | kept |
| acetyl-CoA + ADP → acetate + ATP | rxn00172 | AcdAB NODE_1953_12 (K24012 849/711) | kept |
| acetate ⇄ acetyl-CoA (AMP-forming) | rxn00175 | acsA NODE_1953_10, NODE_3428_5, NODE_1339_9 (K01895) | kept |
| acetate transport | rxn05488 | SatP family NODE_227_53 (K07034) | kept |
| F6P → acetyl-P; acetyl-P → acetyl-CoA | rxn00548, rxn00173 | **none** (no K01621, K01632, K00625, K13788 or K00925) | kept: rxn00548 is essential for biomass (probably as the E4P source). Neither carries catabolic flux at any curated optimum |
| F420H2 → Fd | rxn40355 | FpoF NODE_243_12 | restricted to Fd → F420 |

**The CODH/ACS + AcdAB shunt is gene-supported and was kept. After curation it gives no ATP.**
- With CH4 export blocked and no biomass, the best use of 1 MMA is 0.5 acetate + 1 H2, for 0 ATP. Methanogenesis gives 0.6875 ATP per MMA.
- **So methanogenesis wins outright.** Acetate is 0 at maximum growth and at the substrate-minimising optimum (FVA 0–2e-5 mmol per g biomass).
- Only with all ion couplings halved does the shunt compete; see §4.

---

## 3. Curation list

There are 19 logged changes (`changes` in the validation JSON):
- **7 reactions removed**, 4 of them with no replacement;
- **4 added**: 3 replacements and Fno;
- **3 re-directed**, Rnf also re-stoichiometried;
- **5 GPR-only fixes**.

Reactions 718 → 715; genes 509 → 522.

### 3.1 Replacements (new IDs, because the equation differs from the ModelSEED entry)

All three new reactions are element- and charge-balanced (checked by cobra).

**MCUR_FPO_c0** (ref. rxn24614; EC 1.5.98.3) replaces **rxn24614** (RXN-8106)
- Arm A: MP + Reduced-Factor-F420 (cpd28067) ← MPH2 + Oxidized-Factor-F420 (cpd27742)
- Curated: F420H2 (cpd00792) + MP + 2 H⁺[c] ⇌ F420 (cpd00649) + MPH2 + 2 H⁺[e]
- Bounds −1000, 1000
- GPR: fpoA–O + fpoF, as in rxn24614
- Evidence: Fpo genes (§2.3); 2 H⁺/2e⁻ (§4)

**MCUR_FWD_c0** (ref. ModelSEED rxn46184; EC 1.2.7.12) replaces **rxn11938**
- Arm A: formyl-MFR + NAD⁺ + H2O → CO2 + MFR + NADH
- Curated: formyl-MFR (cpd00735) + H2O + 2 Fd(ox) ⇌ CO2 + MFR (cpd00643) + 2 Fd(red) + H⁺
- Bounds −1000, 1000
- GPR: fmdA and fmdC and (fmdB or fwdB) and (fmdD or fwdD)
- Evidence: Fmd/Fwd genes; the enzyme is ferredoxin-dependent (Thauer et al. 2008 †); no NAD⁺ form is known

**MCUR_HDRDE_c0** (ref. rxn03126; EC 1.8.98.1) replaces **rxn03126**
- Arm A: H⁺ + CoM-S-S-CoB + MPH2 ⇌ CoM + CoB + MP, GPR hdrA or hdrB or hdrC
- Curated: CoM-S-S-CoB + MPH2 + 3 H⁺[c] → CoM + CoB (HTP) + MP + 2 H⁺[e]
- Bounds 0, 1000
- GPR: (hdrD NODE_664_33 or NODE_1064_33) and hdrE NODE_664_32
- Evidence: HdrDE genes; proton-translocating (Ide et al. 1999 †); coupling in §4

### 3.2 Removed without replacement

| Reaction | Equation | Added by | Why removed |
|---|---|---|---|
| rxn05759 | H2 + 2 Fd(ox) → 2 Fd(red) + 2 H⁺ | ATP-correction core gap-fill (H2.CO2 test) | No Ech or ferredoxin-reducing hydrogenase gene. With it, Rnf + Vht form an ATP-generating H2/Fd cycle, and arm A grows on H2/CO2 (0.53 h⁻¹) |
| rxn06299 | F420 + H2 → F420H2 | ATP-correction core gap-fill | No FrhA/FrhG; the two K00441 hits are single FrhB/FdhB/FpoF-type subunits |
| rxn08518 | formate + H⁺ → CO2 + H2 | methane gap-fill (gapfill_2) | No Fdh or Hyc genes; it only made CH4 possible apart from growth |
| rxn10122 | NADH + ubiquinone-8 + 4.5 H⁺[c] ⇌ NAD⁺ + ubiquinol-8 + 3.5 H⁺[e] | draft (GPR FpoK ∧ FpoL) | Its genes are Fpo subunits, now in MCUR_FPO, and methanogens have no ubiquinone. With the quinone-linked SDH rxn09272 it made succinate an energy source (8.1 ATP per 10 succinate, no methylamine needed). Growth on MMA drops 1.4 % |

### 3.3 Re-directed, re-stoichiometried or added

| Reaction | Equation | Before → after | Evidence |
|---|---|---|---|
| rxn03127 (Mcr) | CH3-S-CoM + CoB → CH4 + CoM-S-S-CoB + H⁺ | `=` → `>` | methane-forming direction (README caveat, repair 1b) |
| rxn40355 | F420H2 + 2 Fd(ox) ⇌ F420 + 2 Fd(red) + 2 H⁺ | `=` → `<` (Fd → F420 only) | F420H2 → Fd is uphill (E°′ F420 about −360 mV, Fd ≤ −420 mV); this step was arm A's free reductant for CODH |
| rxa45615 (Rnf) | 2 Fd(red) + MP + 2 H⁺[c] + 3 Na⁺[c] → 2 Fd(ox) + MPH2 + 3 Na⁺[e] | Equation: template's 3 Na⁺ + 2 H⁺ out → **3 Na⁺**. Direction: template `=`, arm A `<` after ATP correction → **`>` only**. GPR: OR → rnfA∧B∧C∧D∧E∧G∧mmcA | §4 |
| **rxn16596 (Fno), added** | NADP⁺ + F420H2 → NADPH + H⁺ + F420 | new, `>` | NODE_209_20 (RAST; eggNOG K06988); same equation as the template reaction. It changes no number in §5–6 |

**Fpo stays reversible.** It is treated as reversible, with the reverse direction costing Δp (2 H⁺ taken back in). With Rnf forward-only and no HdrA2B2C2, reverse Fpo does not open growth on H2/CO2 (0 h⁻¹, §5.2).

**HdrA2B2C2 (template rxn43076) was not added.** Its genes are present (§2.3). Tested on the model's F420 and Fd pools in the bifurcating direction, it gives a route from H2 to reduced ferredoxin: H2 → Vht → MPH2 → reverse Fpo → F420H2 → HdrA2B2C2 → Fd(red) + CoM-SH/CoB-SH.
- With it, the curated model **grew on H2/CO2 (0.239 h⁻¹)**.
- On the DMA + H2 medium it raised CH4 from 146 to 230 by reducing CO2.
- On methylamines alone it carried no flux (growth and k unchanged).
- *Methanosarcina* needs Ech to reduce ferredoxin with H2 for H2/CO2 growth (Meuer et al. 2002 ‡). This genome has no Ech, and growth of a *Methanolobus* on H2/CO2 is not reported in the sources checked. So the reaction is left out and its genes are recorded here.

### 3.4 GPR-only fixes (no flux change)

| Reaction | Before | After |
|---|---|---|
| rxn24608 (MtmB) | (296_14 or _15 or _29 or _30) and (296_13 or _28) | ((296_14 and _15) or (296_29 and _30)) and (296_13 or _28) |
| rxn24610 (MtbB) | 296_25 and (296_26 or _27) | (296_26 and _27) and 296_25 |
| rxn24612 (MttB) | (296_19 or _20) and 296_21 | (296_19 and _20) and 296_21 |
| rxn08173 (A1AO) | 8 subunits OR-ed | 8 subunits AND-ed |
| rxn05209 (Na⁺/H⁺) | NODE_835_14 or NODE_1984_16 | … or (mrpA and mrpC–G) |

Each amber join was checked on the contig: one in-frame TAG, same frame. The joined proteins are MtmB 458 aa (UAG at codon 202), MtbB 487 aa (codon 376) and MttB 496 aa (codon 334). The MtmB and MttB positions match the pyrrolysine positions of *M. barkeri* MtmB and MttB ‡.

---

## 4. Ion stoichiometry of the curated energy metabolism

| Step | Curated | Template / arm A | Measured | Published Methanosarcina models † |
|---|---|---|---|---|
| Fpo (F420H2 → MP) | 2 H⁺ / 2e⁻ | none (rxn24614) | 0.9 H⁺/2e⁻ in vesicles (Bäumer et al. 2000 †) | iAF692 1.8; iMB745 2 |
| HdrDE (MPH2 → CoM-S-S-CoB) | 2 H⁺ / 2e⁻ | 0 (rxn03126) | ≤ 0.9 H⁺/2e⁻ in vesicles (Ide et al. 1999 †) | iAF692 1.8 |
| Vht (H2 → MP) | 2 H⁺ / 2e⁻ (template, unchanged) | 2 | ≤ 0.9 H⁺/2e⁻ (Ide et al. 1999 †) | iAF692 1.8; iMB745 excludes Vht |
| Mtr (CH3-H4MPT → CH3-CoM) | 2 Na⁺ / CH3 (template, unchanged) | 2 | Na⁺-translocating; no number verified (Gottschalk & Thauer 2001 †) | iAF692 1.7; iMAC868 2 |
| Rnf (Fd → MP) | **3 Na⁺ / 2e⁻** | 3 Na⁺ + 2 H⁺ | Na⁺-translocating; **no measured ratio** (Schlegel et al. 2012 FEBS J †) | **iMB745 3 Na⁺**; iMAC868 4 Na⁺ |
| A1AO ATP synthase | 4 H⁺ / ATP (template, unchanged) | 4 | couples to H⁺ and Na⁺ (Schlegel et al. 2012 PNAS †); no c-ring count found | iAF692 4 H⁺; iMAC868 3 Na⁺ + 1 H⁺ |
| Na⁺/H⁺ antiporter | 1 : 1 (template, unchanged) | 1 : 1 | – | iMB745 1 : 1; iMAC868 2 : 1 |

**Why Rnf is 3 Na⁺ and not the template's 3 Na⁺ + 2 H⁺.**
- No paper measures Rnf's coupling. 3 Na⁺/2e⁻ is the value of the published *M. acetivorans* model iMB745 (Benedict et al. 2012 †).
- The template's extra 2 H⁺ makes the model inconsistent with its own Vht (2 H⁺/2e⁻). With 5 charges per 2e⁻, the model oxidises a methyl group completely to CO2 + 3 H2 with net ATP: 0.25 ATP per MMA with CH4 blocked (`rnf_sensitivity`). That reaction is endergonic (CH3OH + H2O → CO2 + 3 H2, ΔG°′ ≈ +23 kJ mol⁻¹ ‡).
- Such a route nets (n − 4) ions per methyl, so any n ≤ 4 removes it.

**The resulting budget per 4 MMA** (zero growth, maximum ATP; fluxes in the validation JSON):
- 1 methyl is oxidised: Mtr⁻¹ (−2 Na⁺), Mer⁻¹, Mtd, Mch, Ftr, Fwd.
- 3 methyls are reduced by Mcr.
- Fpo × 2 (+4 H⁺), Rnf × 1 (+3 Na⁺), HdrDE × 3 (+6 H⁺): 11 ions, 2.75 ATP.
- That is **0.917 ATP per CH4**. On DMA + H2 (methyl reduction only): Vht + HdrDE = 4 H⁺, **1.0 ATP per CH4**.
- For comparison, Thauer et al. 2008 † estimate about 1.5 ATP per CH4 for *M. barkeri*. The vesicle values (0.9 H⁺/2e⁻) would give about half the curated coupling.

**Sensitivity** (MMA; `rnf_sensitivity`, `coupling_sensitivity`; free ATP is 0 in every variant):

| Variant (Fpo, HdrDE, Vht H⁺; Mtr Na⁺; Rnf Na⁺) | ATP per MMA | ATP per MMA, CH4 export blocked | k (mmol CH4 g⁻¹) | Yield (g mol⁻¹ CH4) | Acetate at k optimum |
|---|---:|---:|---:|---:|---:|
| template Rnf (2, 2, 2; 2; 3 Na⁺ + 2 H⁺) | 0.8125 | **0.25 (artefact)** | 113.7 | 8.80 | 0 |
| **curated (2, 2, 2; 2; 3)** | **0.6875** | **0** | **126.5** | **7.91** | **0** |
| Rnf 2 Na⁺ (or 1 Na⁺; then Fd → F420 → Fpo takes over) | 0.625 | 0 | 134.8 | 7.42 | 0 |
| iAF692 couplings + iMB745 Rnf (1.8, 1.8, 1.8; 1.7; 3) | 0.644 | 0.075 | 133.3 | 7.50 | 0 |
| iMAC868-like Rnf (2, 2, 2; 2; 4) | 0.75 | 0 | 119.5 | 8.37 | 0 |
| half coupling, all × 0.5 (1, 1, 1; 1; 1.5), the vesicle-measured scale | 0.375 | 0.25 | 159.1 | 6.29 | **57.9 mmol g⁻¹** |

With the ion couplings halved, the ADP-forming acetyl-CoA synthetase's substrate-level ATP becomes competitive, and the substrate-minimising optimum co-produces acetate. Under every literature model coupling, acetate is 0.

---

## 5. Validation

The methyl medium uses the METHYL_MEDIUM minerals (`_MINERALS` in `build_model.py`): no O2, nitrate or sulfate; H2S as the sulfur source; CO2 100. Each substrate is given at 100 mmol gDW⁻¹ h⁻¹, objective bio1. Fluxes are in mmol gDW⁻¹ h⁻¹, growth in h⁻¹.

### 5.1 No ATP or electrons from nothing

- **`free_atp`** maximises ATPM_c0 on water + H⁺; + O2 and salts; + nitrate and salts. The curated model gives **0 / 0 / 0**; arm A also 0 / 0 / 0.
- **H2 + CO2 (ATP test).** The curated model gets **0 ATP**; arm A gets 1.0, through the gene-less H2 → Fd and H2 → F420 reactions.
- **Single non-methyl substrates** (1 mmol each, anaerobic, no biomass; `alt_substrates_ATP_per_mmol`):

| Substrate | Curated | Arm A |
|---|---|---|
| acetate | infeasible | infeasible |
| succinate | infeasible (after removing rxn10122) | 0.125 ATP |
| glycine betaine | infeasible | infeasible |
| H2 + CO2 | infeasible | 0.25 ATP |
| methanol | infeasible | infeasible |
| pyruvate | 1.56 ATP (4 pyruvate → 4 acetate + 3 CO2 + CH4, through PFOR porABCD + AcdAB, all gene-supported) | 1.25 ATP |
| formate | 0.56 ATP, a **residual gap-fill artefact**; see §9 | 0.25 ATP |

### 5.2 Growth and CH4 at maximum growth

| Medium | Arm A growth | Arm A CH4 / acetate | **Curated growth** | **Curated CH4 / acetate** | Curated CO2 out | CH4 per methyl taken up | Growth with acetate export blocked (A / curated) |
|---|---:|---:|---:|---:|---:|---:|---:|
| methanol only | 0 | 0 / 0 | **0** | **0 / 0** | 0 | – | 0 / 0 |
| MMA only | 0.210 | 0 / 70.6 | **0.504** | **63.7 / 0** | 13.0 | 0.637 | 0 / 0.504 |
| DMA only | 0.418 | 0 / 140.5 | **1.007** | **127.4 / 0** | 25.9 | 0.637 | 0 / 1.007 |
| TMA only | 0.418 | 0 / 140.5 | **1.511** | **191.1 / 0** | 38.9 | 0.637 | 0 / 1.511 |
| methanol + MMA + DMA + TMA | 0.418 | 0 / 140.5 | **3.022** | **382.3 / 0** | 77.7 | 0.637 | 0 / 3.022 |
| DMA + H2 + CO2 (manuscript medium) | 1.323 | 229.5 / 0 | **1.220** | **145.9 / 0** | H2 in 93.0, CO2 in 2.4 | 0.729 | – |
| H2 + CO2 only (check) | 0.529 | 51.8 / 0 | **0** | **0 / 0** | – | – | – |

- **Every curated growth is methanogenic.** Blocking acetate export changes nothing. Arm A's growth on DMA, TMA and the mixture is CO2-limited: its acetogenesis takes up CO2 at the 100 bound.
- **CH4 per methyl at maximum growth is 0.637**, against 0.75 for pure catabolism; the difference is methyl carbon built into biomass.
- **Manuscript medium.** The curated model uses H2 only to reduce methyl groups (Vht + HdrDE + Mcr): 146 CH4 from 200 methyls, and it cannot reduce CO2 with H2. Arm A's 229.5 CH4 there came from H2/CO2 methanogenesis through the gene-less rxn05759/rxn06299, with 87 CO2 taken up.

### 5.3 CH4 per substrate against the textbook stoichiometries

Zero growth, maximum ATP, substrate fixed at 1.

| Substrate | Textbook (per mol substrate) | **Curated CH4 / CO2 / NH3** | Curated ATP | Curated, CH4 export blocked | Arm A |
|---|---|---|---:|---|---|
| MMA | 4 MMA → 3 CH4 + CO2 + 4 NH3 (0.75 / 0.25 / 1) | **0.75 / 0.25 / 1.00** | 0.6875 | ATP 0 (0.5 acetate + 1 H2) | ATP 0.25; CH4 0; acetate 0.75; CO2 −0.5 |
| DMA | 2 DMA → 3 CH4 + CO2 + 2 NH3 (1.5 / 0.5 / 1) | **1.50 / 0.50 / 1.00** | 1.375 | ATP 0 | ATP 0.50; acetate 1.5 |
| TMA | 4 TMA → 9 CH4 + 3 CO2 + 4 NH3 (2.25 / 0.75 / 1) | **2.25 / 0.75 / 1.00** | 2.0625 | ATP 0 | ATP 0.75; acetate 2.25 |
| methanol | 4 CH3OH → 3 CH4 + CO2 + 2 H2O | not possible (no MtaB/MtaC); the what-if gives 0.75 / 0.25 | – | – | not possible |

### 5.4 ATP tests of the build (reported, not tuned to)

| Test medium (threshold = 1.2 × Core-V6) | Arm A ATP | Curated ATP | Curated passes |
|---|---:|---:|---|
| empty (1e-5) | 0 | 0 | yes |
| H2.CO2 (1.2) | 1.0 | **0** | yes |
| Methanamine.H2 (0.6) | 0.5 | **1.0** | **no** |
| Dimethylamine.H2 (1.2) | 1.0 | **2.0** | **no** |
| Trimethylamine.H2 (1.8) | 1.5 | **3.0** | **no** |

The three failures come from the thresholds, not from new energy leaks:
- Core-V6's own Hdr (rxn03126) translocates no protons, so the core reduces one methyl with H2 for 2 H⁺ (Vht only), i.e. 0.5 ATP. The thresholds are 1.2× that.
- The curated model gives 1.0 ATP per methyl + H2 (4 H⁺ at 4 H⁺/ATP), with the coupling of §4.
- The failures scale exactly with the methyl groups per substrate (1, 2, 3).

### 5.5 Aerobic community media and blocked reactions

- **Aerobic media.** Tested on the community MEDIA (aerobic lactate + MePn + sulfate) with MMA added (`aerobic_community_media`):
  - **No growth, because the model needs sulfide as its sulfur source and MEDIA supplies sulfate.** With H2S added it grows exactly as anaerobically (0.504 h⁻¹, CH4 63.7) and takes up **no O2**; it has no oxidase.
  - It cannot grow on lactate or MePn (no lactate exchange).
- **Blocked reactions.** Internal reactions go from 662 to 659; blocked reactions from 262 (0.396) to 263 (0.399). The ubiquinone reactions became dead once rxn10122 was removed. No curated reaction is blocked.

---

## 6. k: CH4 per gram of biomass formed

**Protocol.** bio1 is fixed at 1e-3 h⁻¹ and substrate uptake is minimised, as the community protocol minimises lactate. For single substrates the objective is that substrate's uptake. For the mixed medium two objectives were run: (a) total methyl-group uptake, 1·methanol + 1·MMA + 2·DMA + 3·TMA; (b) the plain sum of the four uptakes. FVA on CH4, acetate, formate, H2 and CO2 export was run at that optimum.

| Medium | Uptake per g biomass formed | **k (mmol CH4 g⁻¹)** | Yield (g mol⁻¹ CH4) = 1000/k | FVA CH4 | Acetate / formate / H2 export |
|---|---|---:|---:|---|---|
| methanol | – | infeasible (no methanol route); what-if 126.5 | – | – | – |
| MMA | 198.5 mmol MMA | **126.47** | **7.91** | 126.472–126.472 | 0 / 0 / 0 (≤ 4e-5) |
| DMA | 99.3 mmol DMA | **126.47** | 7.91 | same | 0 / 0 / 0 |
| TMA | 66.2 mmol TMA | **126.47** | 7.91 | same | 0 / 0 / 0 |
| mixed, objective (a) | 198.5 mmol methyl (degenerate; solver picked DMA) | **126.47** | 7.91 | same | 0 / 0 / 0 |
| mixed, objective (b) | 66.2 mmol TMA (picks TMA, 3 methyls per molecule) | **126.47** | 7.91 | same | 0 / 0 / 0 |

- **k is the same for every methylamine.** Each methyl group enters as methyl-CoM and demethylation is energy-neutral in the model, so k depends only on CH4 per methyl.
- Per g biomass formed, 198.5 methyl groups go to 126.5 CH4 + 25.7 CO2, and the rest is built into biomass.
- k spans 114–159 mmol g⁻¹ across the coupling variants of §4.
- **Arm A's k is 0**: at the same optimum it exports 335.7 mmol acetate per g and no CH4.

**Against measured yields (Y(CH4), g dry weight per mol CH4):**

| Organism, substrate | Y(CH4) | k equivalent (mmol g⁻¹) | Source |
|---|---:|---:|---|
| *M. barkeri*, methanol | 7.2 | 139 | Weimer & Zeikus 1978 † (abstract) |
| *M. barkeri*, methanol | 5.6 | 179 | Ranalli et al. 1986 † (4.2 g per mol methanol × 4/3) |
| *M. barkeri*, methanol | 6.0–7.3 | 137–167 | Mazumder et al. 1986 † (abstract, 0.14–0.17 g per g methanol) |
| *M. barkeri*, methanol (older) | 5.1; 4.4 | 196; 227 | Smith & Mah 1978; Stadtman, both as quoted by Ranalli et al. 1986 † |
| *M. barkeri*, MMA / DMA / TMA | 5.8 (the same for all three) | 172 | Hippe et al. 1979 † (abstract) |
| *Methanococcoides* NM1, TMA | 1.8 | 556 | Watkins et al. 2014 † (4.11 g per mol TMA × 4/9; dry weight from protein) |
| cytochrome-containing methanogens | "up to 7" | ≥ 143 | Thauer et al. 2008 † |
| model iMB745, methanol | 5.3 predicted / 6.9 measured | 189 / 145 | Benedict et al. 2012 † |
| **curated model, any methylamine** | **7.91** (sensitivity 6.3–8.8) | **126.5** (114–159) | this work |

- **The model's yield is at the top of the *Methanosarcina* range, or slightly above it.** It is 1.36× Hippe's methylamine value for *M. barkeri* and 1.1× Weimer & Zeikus' methanol value.
- **It is 4.4× the only measured value for a closer methylamine specialist**, *Methanococcoides* NM1 on TMA. That value is a protein-based estimate.
- Two features of the model push its yield up:
  - The biomass is the template's: GAM 40 mmol ATP gDW⁻¹ plus polymerisation costs, and **no NGAM** (ATPM lower bound 0). By contrast, iAF692 fixes NGAM at 0.8 and iMB745 uses GAM 65 and NGAM 2.5 mmol ATP gDW⁻¹ (h⁻¹) †.
  - The ion couplings in §4 are the models' 2 H⁺/2e⁻, above the 0.9 measured in vesicles.
- No yield for any *Methanolobus* was found.

**The model's maximum methanogenesis rate is set by the uptake bounds, not by physiology.**
- At maximum growth on MMA the uptake sits at its 100 mmol gDW⁻¹ h⁻¹ bound, giving CH4 63.7 mmol gDW⁻¹ h⁻¹ and µ = 0.50 h⁻¹. Both scale linearly with that bound.
- Nothing in the model (no kinetics, no enzyme capacity) limits the rate. The mixed medium's 3.0 h⁻¹ and 382 mmol CH4 gDW⁻¹ h⁻¹ only reflect four substrates each at 100.
- Only the ratio is meaningful: q_CH4 = k × µ = 126.5 × µ mmol CH4 gDW⁻¹ h⁻¹.

**Effect on the community evaluation.** `../../community_eval/methanogen_test.py` took k from the arm B–D methanol builds (196–406 mmol g⁻¹). Those builds make CH4 through the gene-less Frh and export formate (§2.2), and this genome cannot use methanol. The curated k (126.5, on methylamines) and methylamine-only substrate use should replace them there, and the methanol pool should not be credited to this species.

---

## 7. Glycine betaine as a methyl donor

**No evidence that this genome uses glycine betaine as a methanogenic substrate.**
- **What the pathway needs.** Methanogenic betaine use requires a non-pyrrolysine MttB-family betaine:corrinoid methyltransferase, MtgB, with its corrinoid protein MtgC and an MtaA-like CoM methyltransferase. This was shown for *Methanolobus vulcani* B1d (Creighbaum et al. 2019 †) and for the bacterial MtgB of *Desulfitobacterium* (Ticak et al. 2014 †). The *M. vulcani* type strain, *M. mangrovi* and *M. sediminis* do not use betaine (Ticak et al. 2015 † and species descriptions); 3 of 9 *Methanococcoides* strains do (Watkins et al. 2014 †).
- **MttB-family proteins in this genome.** The only one is the Pyl MttB (NODE_296_19 + _20, Pfam MTTB, in-frame UAG). The 6-frame search finds exactly one K14083 ORF in each of the six genomes (§2.4). There is no non-Pyl MttB homologue, so **no MtgB**. NODE_296_31, a 457-aa non-Pyl protein, is in the MtmB family (Pfam MtmB), not MttB.
- **MtaA-like genes.** These are present: NODE_485_13, NODE_1297_10, and NODE_2869_3 next to an OtaABC glycine betaine ABC transporter (NODE_2869_4–7). Their KO K14080 now covers both methanol- and betaine-specific CoM methyltransferases (EC 2.1.1.246, 2.1.1.377). Without an MtgB they cannot demethylate betaine.
- **Betaine transporters.** OtaABC (NODE_2869_4–6), a second OpuA-like ABC system (NODE_296_38–41) and BetT (NODE_843_24). These fit **betaine as a compatible solute** in a halophile.
- **Result.** The model has a betaine exchange but no betaine catabolism (betaine as sole substrate: infeasible). Nothing was added.

---

## 8. *Roseovarius*.2 (`R2_restored_C_black__metabat2.bin.011_sub`): the "methane monooxygenase" rxn00844

**Verdict: not a methane monooxygenase. The gene is the large subunit of a propane 2-monooxygenase (PrmA, SDIMO group 5). Recommend blocking rxn00844.**

The gene is `NODE_1821_length_19323_cov_3.402949_4` (554 aa). Together with NODE_1821_3 and _2 (+ NODE_10142_5) it forms rxn00844's GPR. The second *Roseovarius* with the same RAST call, `R2_B_D2__metabat2.bin.006_sub` (`NODE_33_length_127140_cov_4.273838_82`, 554 aa), carries the same operon.

| Evidence | Result |
|---|---|
| RAST | "Methane monooxygenase component A alpha chain (EC 1.14.13.25)"; the neighbours get the other sMMO names |
| KOfam best KO | **K18223 propane 2-monooxygenase large subunit, 1016/744** (bin.006: 1019/744) |
| KOfam K16157 (MmoX), same protein | **228 vs threshold 1169** (bin.006: 219) |
| other SDIMO alpha KOs | K21317 butane MO 309/939; K16242 phenol/toluene MO 175/449; K15760 toluene MO 154/301 |
| eggNOG | K18223 |
| Bakta | "Methane monooxygenase" (name only, no KO). Its neighbours are named "Propane monooxygenase small subunit PrmC" and "… reductase subunit PrmB" |
| Operon (bin.011 on the minus strand, bin.006 on the plus strand) | **prmA** (K18223), **prmB** reductase (K18225 502/401), **prmC** β (K18224 493/475), **prmD** coupling protein (K18226). Then an amidohydrolase (K07045), an MIP18-family protein, GroEL (K04077) and a secondary alcohol dehydrogenase Adh1 (K18382 493/424). **This α–reductase–β–coupling order is the group-5 order** (Yang, Kertesz & Coleman 2025 †). sMMO is mmoXYBZDC (α, β, coupling, γ, accessory, reductase). **There is no γ subunit** (MmoZ, K16159) |
| Smith–Waterman vs reference α subunits (BLOSUM62, gap −11/−1; UniProt FASTA in the session scratchpad, accessions below) | **PrmA *Gordonia* TY-5 Q768T5: 65.2 % identity over 532 aa. PrmA *Rhodococcus jostii* RHA1 Q0SJK9: 65.5 % over 533 aa.** ThmA (group 5, THF MO) 41 %; group-4 alkene MO (AmoC, EtnC) 37–39 %; group-6 Pr1MO 37.5 %; group-3 SmoX/BmoX 33–37 %; **MmoX *M. capsulatus* P22869 32.1 % and *M. trichosporium* P27353 33.8 %**; group-1/2 TMO/PhMO 27–30 % |
| Methanotrophy genes elsewhere | none: no pmoCAB (K10944–K10946) and no mmoXYBZCD (K16157–K16162) in either genome |
| Related reports | Several Roseobacter-clade proteins in GenBank carry the same "methane monooxygenase component A alpha chain" name: SMC10547.1 (*Roseovarius aestuarii*), ASP22755.1, QRF69149.1, GLO69945.1. They are 64–67 % identical to PrmA and 32–33 % to MmoX, i.e. misnamed group-5 enzymes. No genuine sMMO is reported in a Rhodobacteraceae † |

**Effect in the models** (Pi-free lactate 1 + MePn medium, maximum growth; these two models were not changed):

| Model | rxn00844 | Growth | CH4 export | rxn00844 flux | CH4 export FVA at 99.9 % growth |
|---|---|---:|---:|---:|---|
| `R2_restored_C_black__metabat2.bin.011_sub` | as built | 0.0466 | 0.0245 | 0.0063 | 0.0193–0.0261 |
| | blocked | 0.0464 | **0.0306** (= MePn uptake) | 0 | 0.0306–0.0306 |
| `R2_B_D2__metabat2.bin.006_sub` | as built | 0.0482 | 0.0318 | 0 | **0.0257**–0.0318 |
| | blocked | 0.0482 | 0.0318 | 0 | 0.0318–0.0318 |

- Blocking costs 0.4 % growth in bin.011_sub and nothing in bin.006_sub.
- It removes the degeneracy that lets either model oxidise its own MePn-derived CH4.
- It matches repair 1b, which already restricts every other template methane-consuming reaction.
- `../../simulations/` already runs an "mmo closed" sensitivity. With this evidence, closed should be the default.
- **Block rxn00844 (bounds 0, 0) in both models.**
- No source checked here shows a group-5 propane monooxygenase oxidising methane.

---

## 9. Unresolved

1. **Methanol.**
   - Six assemblies lack *mtaCB*. The full metagenome assemblies (on poplar) were not searched, so an unbinned, multi-copy operon is not strictly excluded.
   - If methanol use is wanted anyway, the what-if gives growth 0.504 h⁻¹ on methanol and the same k. It would be an assumption, not gene-supported.
2. **The ion couplings decide k** (§4, 114–159 mmol g⁻¹).
   - Rnf's coupling is unmeasured. The measured vesicle couplings are about half the model values; at that scale acetate co-production appears.
   - A1AO (4 H⁺/ATP), Mtr (2 Na⁺) and the 1:1 antiporter were left at template values. The antiporter treats ΔμNa⁺ and ΔμH⁺ as interchangeable.
3. **The biomass is the template's.** GAM is 40 mmol ATP gDW⁻¹ plus polymerisation, and there is no NGAM. This and the couplings put the model's yield at or above the top of measured *Methanosarcina* yields (§6). No growth-energetics fit to *Methanolobus* was attempted, and none is possible without its yield data.
4. **The gap-filled reactions of arm A were kept**, except the three removed above. Some are gene-less biosynthetic patches that stay essential:
   - phosphoketolase rxn00548 (E4P);
   - NAD⁺ formaldehyde dehydrogenase rxn00429 (formate supply);
   - a lactate shuttle (rxn00542, rxn00499, rxn00171): the only L-lactate source for F420 synthesis (gap-filled rxn10420).
   - Together with PFOR, these act as a **gene-less formate dehydrogenase**: external formate gives 0.56 ATP and 0.25 CH4 per formate.
   - Formate is in none of the media used, and each model is simulated alone. It matters only if formate is ever offered: keep formate uptake closed, or re-gap-fill F420 synthesis.
5. **Pyruvate** gives 1.56 ATP per pyruvate (4 pyruvate → 4 acetate + 3 CO2 + CH4) through gene-supported PFOR, AcdAB and CO2 reduction. Whether *Methanolobus* grows on pyruvate is unknown, so this was left as is.
6. **ATP tests.** 3 of 5 fail because of the Core-V6 thresholds (§5.4). The template defects found here should go into the template curation queue (the templates were not touched):
   - Core-V6 and Archaea V6 rxn03126 (no H⁺, HdrABC roles);
   - rxn24614 (the separate F420 pool);
   - rxn11938 (NAD⁺ Fwd);
   - the Rnf stoichiometry;
   - the K14080 → rxn10568 mapping;
   - FpoK/FpoL → rxn10122.
7. **Methylamine uptake has no GPR.** NODE_296_22 is a plausible permease but was not assigned.
8. **Sulfide and O2.**
   - The model needs H2S; on the aerobic sulfate-only community MEDIA it does not grow.
   - Being anaerobic and O2-insensitive in the model is not the same as being an obligate anaerobe in vivo. The model does not penalise O2.
9. **Not updated:** `../README.md` (caveat 1), `../models_summary.tsv`, `../MANIFEST.md5`, `../../community_eval/methanogen_test.py` and the two *Roseovarius* models. Point the community work at `R1_A_D1__concoct.bin.007__curated.json`.

---

## References

**† Checked in this session** (abstract or full text):
- Bäumer S et al. 2000 *J Biol Chem* 275:17968 (Fpo, 0.9 H⁺/2e⁻)
- Benedict MN et al. 2012 *J Bacteriol* 194:855 (iMB745)
- Coleman NV et al. 2006 *Environ Microbiol* 8:1228
- Creighbaum AJ et al. 2019 *Front Microbiol* 10:2572
- Feist AM et al. 2006 *Mol Syst Biol* 2:2006.0004 (iAF692; couplings read from the BiGG file)
- Gottschalk G, Thauer RK 2001 *Biochim Biophys Acta* 1505:28
- Hippe H et al. 1979 *PNAS* 76:494
- Ide T, Bäumer S, Deppenmeier U 1999 *J Bacteriol* 181:4076
- Katayama T et al. 2014 *IJSEM* 64:2089
- Mazumder TK et al. 1986 *Appl Environ Microbiol* 52:617
- Nazem-Bokaee H et al. 2016 *Microb Cell Fact* 15:10 (iMAC868)
- Ranalli G et al. 1986 *FEMS Microbiol Lett* 35:119
- Schlegel K et al. 2012 *FEBS J* 279:4444 (Rnf)
- Schlegel K et al. 2012 *PNAS* 109:947 (A1AO)
- Shen L et al. 2020 *IJSEM* 70:5586
- Thauer RK et al. 2008 *Nat Rev Microbiol* 6:579
- Ticak T et al. 2014 *PNAS* 111:E4668
- Ticak T et al. 2015 *Arch Microbiol* 197:197
- Watkins AJ et al. 2014 *Appl Environ Microbiol* 80:289
- Weimer PJ, Zeikus JG 1978 *Arch Microbiol* 119:49
- Yang, Kertesz & Coleman 2025 *Environ Microbiol* 27:e70050
- NCBI assembly GCA_020055655.1 (PRJNA365332)
- UniProt Q768T5, Q0SJK9, P22869, P27353, Q8KQF0, I4BTA3, O87082, Q00456, O87798, P19732, UPI00000B2FA5, Q6XBI6, Q9F3V6, Q08KF2

**‡ From memory, not re-checked:**
- Meuer J et al. 2002 *PNAS* 99:5632 (Ech required for H2/CO2 growth of *M. barkeri*)
- Thauer RK, Jungermann K, Decker K 1977 *Bacteriol Rev* 41:100 (ΔG°′ of CH3OH + H2O → CO2 + 3 H2)
- Jasso-Chávez R et al. 2013 *J Bacteriol* (Mrp in *M. acetivorans*)
- Kotani T et al. 2003 *J Bacteriol* 185:7120 (*Gordonia* TY-5 prmABCD; the source of UniProt Q768T5)
- Pyrrolysine at MtmB residue 202 and MttB residue 334 of *M. barkeri* (Hao et al. 2002 *Science*; Srinivasan et al. 2002 *Science*)
