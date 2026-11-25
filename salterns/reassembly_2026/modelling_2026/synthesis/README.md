# Synthesis of the hypothesized methane sources (2026-11-25)

One table and one figure (manuscript Figure 10) that set the hypothesized methane sources of the six unrestored cores side by side, plus an attribution of the measured CH4 among them.

| Script | Output |
|---|---|
| `pathway_summary.py` | `pathway_summary.json`, `pathway_abundance.tsv` (per sample, per 100 genome equivalents, reads and MAGs), `pathway_turnover.tsv` (per core: substrate pool turnover, d; doubling times, h), `attribution.tsv` (per core and scenario) |
| `figure10_pathways.py` | `figure10_pathways.png` / `.svg` |

**Pathways** (the manuscript's hypotheses):
- methylphosphonate via C–P lyase;
- glycine betaine → TMA (betaine reductase);
- choline / carnitine → TMA (CutC, CntA);
- methylamines → CH4;
- methanol → CH4;
- hydrogenotrophic and acetoclastic methanogenesis.

**Inputs:**
- `../../read_markers/results/`: `per_sample_KO.tsv`, `key_unrestored_cores.tsv`, `read_taxonomy.tsv`, `mcrA_taxonomy_by_sample.tsv`, and `methanol_taxonomy.tsv` (step 10d).
- `../simulations/specific_rates.tsv`, `per_gram_yields.tsv`, and `../models/methanogen/curation_k.json`.
- The porewater NMR (`data/Cliff_Sample_Metadata_BGC_NMR.csv`).

**Assumptions:**

| Topic | Assumption |
|---|---|
| *mcrA* groups | Pooled over the 12 samples. Methylotrophic: Methanosarcinaceae and the H2-dependent Methanomassiliicoccales. Hydrogenotrophic: Methanomicrobiales and Methanobacteriaceae. Acetoclastic: *Methanothrix*, *Methanocrinis*. ANME-1-related *Methanospirare* is excluded. |
| Methanol | Only the archaeal *mtaB* reads (21 %) count; the rest are bacterial, often distant. |
| Substrate pools | Concentration × porosity (0.74) × layer thickness. CH4 per molecule: TMA 2.25, DMA 1.5, MMA 0.75, methanol 0.75, formate 0.25, acetate 1, betaine 2.25 (through TMA). |
| Not quantified | Methylphosphonic acid was detected by GC-MS but not quantified. Choline, carnitine and H2 were not measured. |
| Doubling times | Per-gram yields × measured biomass. Hydrogenotrophs and acetoclasts have no model, so their doubling times use the methylotroph k and are lower bounds. |
| Attribution | Equal specific growth rates: share ∝ k × read-based biomass, and the methanogen share is split by *mcrA* reads. Sensitivities: MPn as the only P source; MPn consumers growing 10× faster. |

**Main numbers** (`pathway_summary.json`):
- Methylotrophic methanogens carry 91 % of the methanogen *mcrA* reads, hydrogenotrophic methanogens 9 % and acetoclastic methanogens 0.4 % (3 reads).
- **Substrate pool turnover needed:** methylamines 0.9–3.4 d; methanol 0.6–3.0 d; formate 4–12 h; acetate 8–335 d; betaine 61–168 d.
- **Attribution at equal growth rates:**

| Scenario | MPn consumers' share |
|---|---|
| MPn:Pi ratio 0.07 | 1–8 % (methylotrophic methanogenesis 84–90 %) |
| MPn:Pi ratio 0.33 | 3–24 % |
| MPn the only P source | 9–55 % |
| ratio 0.07, MPn consumers growing 10× faster | 6–45 % |

- **R2 cores:** at a one-day doubling time, growth-coupled methanogenesis supplies only 8–19 % of the measured CH4, so the required specific activity is 19–45 mmol CH4 g⁻¹ h⁻¹.
