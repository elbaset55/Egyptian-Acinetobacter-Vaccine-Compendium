# Subunit, Outer Membrane Protein, and Multi-Epitope Vaccines Against Carbapenem-Resistant *Acinetobacter baumannii*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24214866.svg)](https://doi.org/10.5281/zenodo.24214866)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-blue.svg)](https://creativecommons.org/licenses/by/4.0/)
[![PRISMA 2020](https://img.shields.io/badge/PRISMA-2020_Compliant-green.svg)](http://www.prisma-statement.org/)
[![SYRCLE](https://img.shields.io/badge/RoB-SYRCLE_Audited-teal.svg)](https://www.radboudumc.nl/en/research/departments/syrcle)
[![R Core](https://img.shields.io/badge/R-v4.6.2-brightgreen.svg)](https://www.r-project.org/)
[![Target Journal](https://img.shields.io/badge/Target_Journal-Frontiers_in_Immunology-red.svg)](https://www.frontiersin.org/journals/immunology)

> **Official Open-Science Research Compendium** accompanying the systematic review and meta-analysis:  
> **"Subunit, Outer Membrane Protein, and Multi-Epitope Vaccines Against Carbapenem-Resistant *Acinetobacter baumannii*: A PRISMA 2020 Systematic Review, Preclinical Meta-Analysis, and Methodological Quality Appraisal"**  
> **Authors:** Abdelbaset A. Diab & M N Seddek  
> **Affiliation:** Department of Zoology, Faculty of Science, Benha University, Benha 13518, Egypt  
> **Permanent Zenodo DOI:** [10.5281/zenodo.24214866](https://doi.org/10.5281/zenodo.24214866) | **Zenodo Record:** [24214867](https://zenodo.org/records/24214867)

---

## 📌 Overview

Carbapenem-resistant *Acinetobacter baumannii* (CRAB) represents a critical WHO Priority 1 superbug characterized by extreme multi-drug resistance and hospital intensive care unit (ICU) mortality rates often exceeding 60%. This repository houses the complete, reproducible computational and empirical evidence corpus for our tri-track systematic review:

1. **Track 1 (Preclinical Murine Challenge Meta-Analysis):** Quantitative synthesis of $k = 25$ controlled animal challenge trials ($N = 510$ mice) evaluating recombinant outer membrane protein subunits (OmpA, BamA, OmpW, BauA) and chimeric constructs using DerSimonian–Laird (DL), Hartung–Knapp–Sidik–Jonkman (HKSJ), Generalized Linear Mixed Models (GLMM), and additive Risk Difference (RD) prediction intervals.
2. **Track 2 (Computational Quality Appraisal):** Objective 7-domain translational quality audit of $k = 22$ published *in silico* reverse vaccinology pipelines (2018–2026), exposing structural inclusion body aggregation risks, HLA-DRB1 restriction bias, and absence of human autoantigen screening.
3. **Track 3 (Regional ICU Molecular Epidemiology):** Synthesizing $k = 17$ clinical surveillance cohorts ($N = 1,813$ isolates) across Egypt and the MENA region, documenting dominant IC2/ST2 lineages and uncovering automated VITEK-2 diagnostic Very Major Errors (0.0% colistin resistance) contrasted against reference broth microdilution (BMD; up to 52.9%).

---

## 📂 Repository Structure

```
Egyptian-Acinetobacter-Vaccine-Compendium/
├── README.md                               # Project documentation & execution guide
├── LICENSE                                 # Creative Commons Attribution 4.0 International
├── data/
│   ├── R_meta_analysis_25_trials.csv       # Primary 25-trial preclinical challenge dataset
│   ├── R_leave_one_out_sensitivity.csv    # 25-iteration leave-one-out sensitivity analysis
│   ├── meta_analysis_results.json          # Pooled statistical effects, heterogeneity, & models
│   ├── pubmed_search_summary.json          # PubMed query metadata & hit statistics
│   ├── Q1_OMP_Subunit_InVivo.json          # Track 1 primary PubMed JSON corpus (46 PMIDs)
│   ├── Q2_MultiEpitope_InSilico.json       # Track 2 primary PubMed JSON corpus (56 PMIDs)
│   ├── Q3_Egypt_CRAB_Epidemiology.json     # Track 3 primary PubMed JSON corpus (74 PMIDs)
│   └── Q4_Targeted_Antigens_Extra.json     # Track 4 targeted high-recall PMIDs (43 unique PMIDs)
├── scripts/
│   └── generate_figures_in_R.R             # End-to-end reproducible R script generating Figs 1–9
└── docs/
    ├── PRISMA_2020_CHECKLIST.md            # Completed 27-item PRISMA 2020 compliance checklist
    └── Table_S1_PICOS_Framework.md         # Structured tri-track PICOS eligibility criteria
```

---

## 🚀 Quick Start & Reproducibility Guide

All statistical calculations, forest plots, funnel plots, bubble meta-regressions, and diagnostic error visualizations are 100% reproducible using **R (version $\ge 4.2.0$)**.

### 1. Prerequisites (R Libraries)
Open R or RStudio and install the required dependencies:
```R
install.packages(c("meta", "metafor", "ggplot2", "dplyr", "readr", "tidyr", "gridExtra"))
```

### 2. Run the Complete Meta-Analysis Pipeline
Clone or download the repository, set your working directory to the project root, and execute:
```R
source("scripts/generate_figures_in_R.R")
```
This generates all publication-grade figures (Figures 1–9) and outputs the quantitative model tables:
* **Pooled Risk Ratio (DL):** $\text{RR} = 14.38$ $[95\%\ \text{CI: } 8.59\text{--}24.05]$
* **Small-Sample Adjustment (HKSJ):** $\text{RR}_{\text{HKSJ}} = 14.38$ $[95\%\ \text{CI: } 8.36\text{--}24.72,\ p < 0.0001]$
* **Generalized Linear Mixed Model:** $\text{RR}_{\text{GLMM}} = 14.12$ $[95\%\ \text{CI: } 7.95\text{--}25.10]$
* **Pooled Risk Difference (Additive):** $\text{RD}_{\text{DL}} = +85.1\%$ $[95\%\ \text{CI: } 79.2\%\text{--}91.0\%]$
* **Additive 95% Prediction Interval:** $[+58.2\%,\ +100.0\%]$

---

## 🔬 Key Methodological Innovations

* **Resolution of the Scale-Dependent Heterogeneity Paradox:** Demonstrates mathematically why multiplicative estimators encounter artificial homogeneity ($I^2_{\text{RR}} = 0.0\%$) due to Haldane–Anscombe continuity inflation under zero-control events ($e_c = 0$ in 22/25 studies), while additive models ($I^2_{\text{RD}} = 46.7\%$) preserve true biological divergence.
* **The "Inclusion Body Trap" Resolution:** Highlights that full-length $\beta$-barrel expression produces denatured inclusion bodies eliciting non-neutralizing antibodies ("ELISA vanity"), whereas multivalent chimeric loop constructs (assembling ECLs of BamA, OmpA, LptD, OmpW) confer superior protection (92.1% survival, $-3.84\ \log_{10}\ \text{CFU/g}$ pulmonary clearance) by triggering an avidity-driven C1q hexamer nucleation threshold.
* **Biophysical Mechanism of VITEK-2 Colistin Failure:** Elucidates why amphipathic polycationic polymyxins adsorb to polystyrene automated well walls, causing Very Major Errors (0.0% reported colistin resistance) in clinical laboratories.

---

## 📜 Citation

If you use this dataset, code, or evidence synthesis in your research, please cite:

```bibtex
@article{DiabSeddek2026CRABVaccine,
  author    = {Abdelbaset A. Diab and M N Seddek},
  title     = {Subunit, Outer Membrane Protein, and Multi-Epitope Vaccines Against Carbapenem-Resistant Acinetobacter baumannii: A PRISMA 2020 Systematic Review, Preclinical Meta-Analysis, and Methodological Quality Appraisal},
  journal   = {Frontiers in Immunology},
  year      = {2026},
  note      = {Data and Code archived under Zenodo DOI: 10.5281/zenodo.24214866}
}
```

---

## 📄 License
This repository and all included data, scripts, and documentation are licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.
