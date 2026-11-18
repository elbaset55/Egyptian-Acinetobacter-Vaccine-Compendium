# Egyptian *Acinetobacter baumannii* Vaccine & Epidemiology Compendium

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24218732.svg)](https://doi.org/10.5281/zenodo.24218732)
[![PROSPERO](https://img.shields.io/badge/PROSPERO-CRD420261572661-0F766E.svg)](https://www.crd.york.ac.uk/prospero/display_record.php?ID=CRD420261572661)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Open-science data, PRISMA 2020 / SYRCLE audit matrices, 600-DPI publication figures, and reproducible Python 3.13 & R 4.6.2 analytical code accompanying:

> **Diab AA, Seddek MN (2026).** *Subunit and multi-epitope vaccines against carbapenem-resistant Acinetobacter baumannii: a systematic review and meta-analysis.*  
> **Affiliation:** Biotechnology Program, Department of Zoology, Faculty of Science, Benha University, Benha 13518, Egypt  
> **PROSPERO Registration:** `CRD420261572661` | **Zenodo Archive:** [https://doi.org/10.5281/zenodo.24218732](https://doi.org/10.5281/zenodo.24218732)

---

## Verified Three-Track Summary (`n = 55` Unique Peer-Reviewed Publications)

1. **Track 1 — Preclinical Murine Survival Meta-Analysis (`k = 22` Independent Experiments, `20` Publications, `48` Active Arms, `30` Platform Units, `N = 574` Mice):**
   - **Pooled Raw Survival:** `71.8%` (`280/390`) in vaccinated mice vs. `7.1%` (`13/184`) in mock controls.
   - **Multiplicative Risk Ratio (`DL` + `HKSJ`):** `RR_DL = 3.69` (`95% Wald CI: 2.61–5.22`; `HKSJ 95% CI: 2.56–5.33`; `95% PI: 2.56–5.33`; `z = 7.41, p < 0.0001`; `Q = 19.95, df = 21, p = 0.52, I² = 0.0%`).
   - **Additive Risk Difference (`RD`):** `+62.6` percentage points (`Wald 95% CI: 54.0–71.3`; `HKSJ 95% CI: 53.5–71.8`; `95% PI: +29.3 to +96.0 pp`; `Q_RD = 49.55, p = 0.0004, I²_RD = 57.6%`).
   - **Leave-One-Publication-Out Range (`20` iterations):** `RR = 3.52` to `6.34` (all `p < 0.0001`; lowest HKSJ lower bound: `2.42`).

2. **Track 2 — Methodological Quality Audit of *In Silico* MEVC Pipelines (`k = 22` Studies, 2018–2026):**
   - **Verified 7-Domain Rubric Score (`0–14`):** `4.91 ± 2.00 / 14` (`Median = 5.0`, `IQR = 3.25–6.00`, `Range = 2–8`).
   - `90.9%` (`20/22`) lacked strict β-barrel extracellular-loop confinement; `100%` (`22/22`) lacked empirical multi-locus regional HLA optimization; `90.9%` (`20/22`) omitted peptide-level short *k*-mer human proteome screening.

3. **Track 3 — Egyptian & MENA Clinical Molecular Epidemiology (`k = 17` Cohorts, `1,724` Analyzed Isolates/Episodes):**
   - **Strictly Unselected Carbapenem Resistance:** `90.7%` (`500/551`) across `k = 6` unselected Egyptian cohorts (`91.5%` [`550/601`] incl. outbreak MDR ICU surveillance; `91.4%` [`687/752`] across `k = 8` unselected Egyptian & MENA cohorts).
   - **Reference ISO 20776-1 Broth Microdilution Colistin Resistance:** `1.2%` (`2/172`, Türkiye) and `5.0%` (`2/40`), `49.0%` (`49/100`), and `52.9%` (`9/17`) across Egyptian cohorts (plus `40.9%` [`9/22`] in a companion Giza study).

---

## Quick-Start Reproduction

```bash
# 1. Exact Python 3.13 Track 1 DerSimonian-Laird, HKSJ, Mantel-Haenszel & Sensitivity Verification
python scripts/compute_track1_exact.py

# 2. Cross-Verified R 4.6.2 Replication across Tracks 1, 2, and 3
Rscript scripts/generate_figures_in_R.R

# 3. Regenerate All 600-DPI Publication Figures into figures/
python scripts/consolidate_and_generate_figures.py
```

---

## Repository Structure

- `data/Track1_22_experiments_verified.csv` — Master 22-experiment survival dataset (`ev/nv`, `ec/nc`, routes, doses, strains, `RR`, `RD`, weights).
- `data/Track1_48_arms_verified.csv` — Complete 48 active vaccine arms across the 22 experiments.
- `data/Track1_30_platform_units_verified.csv` — 30 experiment-by-platform subgroup units.
- `data/Track1_20_leave_one_out_verified.csv` — 20 leave-one-publication-out sensitivity iterations.
- `data/Track1_20_SYRCLE_verified.csv` — 20-publication SYRCLE risk-of-bias & ethics reporting audit.
- `data/Table_S2_77_excluded_studies_verified.csv` — All 77 excluded preclinical survival & candidate vaccine studies with reasons.
- `data/Track2_22_insilico_rubric_verified.csv` — Verified 7-domain rubric matrix across 22 *in silico* MEVC studies.
- `data/Track3_17_cohorts_verified.csv` — Verified 17-cohort Egyptian & MENA clinical epidemiology dataset.
- `data/meta_analysis_results.json` — Machine-precision JSON summary of all statistical models.
- `data/pubmed_search_summary.json` & `data/Q1..Q4_*.json` — Raw literature search and screening exports.
- `figures/Fig_SR1..Fig_SR10_*.png` — All 10 publication-grade 600-DPI figures (Figures 1–9 + 20-study SYRCLE matrix).
- `docs/PRISMA_2020_CHECKLIST.md`, `docs/Table_S1_PICOS_Framework.md`, & `docs/Supplementary_Material_Checklist_S1_Tables_S1_S2.pdf` — PRISMA 2020 checklist, three-track PICOS framework, and compiled Supplementary Material PDF.
- `scripts/compute_track1_exact.py`, `scripts/generate_figures_in_R.R`, `scripts/consolidate_and_generate_figures.py` — Standalone Python 3.13 & R 4.6.2 reproduction and 600-DPI figure scripts.
- `CITATION.cff` & `.zenodo.json` — Standardized bibliographic and Zenodo archival metadata.
