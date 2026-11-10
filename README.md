# Macrophage Spatial Polarity Index (MSPI) in Breast Cancer

Analysis code for the manuscript *Macrophage spatial polarity to T cells predicts prognosis in young women with luminal breast cancer*.

## Overview

This repository contains the scripts that compute the spatial affinity scores and generate the figures and tables reported in the manuscript. The analysis quantifies a density-independent spatial affinity metric between macrophage subtypes (M1-like, M2-like) and T cells in multiplex immunohistochemistry (mIHC) breast cancer tissue, and evaluates its prognostic value as the Macrophage Spatial Polarity Index (MSPI).

## The affinity score

For every T cell (anchor), the scripts record the nearest-neighbour distance (NND) to the nearest M1-like (or M2-like) macrophage, and to the nearest cell of a random subsample of the total macrophage pool drawn to the same size (repeated 5 times). A nearer neighbour closer than 1 um is treated as the same position and skipped. Per patient, NNDs at or above a distance cap (set at the top of each affinity script) are excluded, separately for the subset and the pool, and

    delta = mean(subset NND) / (mean(pool NND) + mean(subset NND))

delta = 0.5 means no preference; values below 0.5 mean that the subset lies closer to T cells than the macrophage pool as a whole. Script 05 negates delta, so higher scores mean higher affinity, and defines MSPI = M2 score - M1 score.

## Repository structure

```
├── python/
│   └── 01_affinity_computation.py              # Validation cohort: affinity scores from cell data
├── R/
│   ├── 01_discovery_affinity_computation.R     # Discovery cohort: affinity scores from cell data
│   ├── 01b_discovery_scores_from_distances.R   # Discovery cohort: published scores from stored distances
│   ├── 02_Fig_S1_exploratory_density_vs_affinity.R
│   ├── 02b_Fig1c_exploratory_hugging_fractions.R
│   ├── 03_Fig1d_exploratory_ECDF_curves.R
│   ├── 04_Fig1fg_discovery_KM_survival.R
│   ├── 05_MSPI_integrated_analysis.R
│   ├── 06_clinical_subtype_figures.R
│   └── 07_Fig1a_patient_flow_and_tables.R
├── LICENSE
└── README.md
```

All R scripts are run from the repository root. They read inputs from `./data/` and write results to `./output/`.

## Script descriptions and output mapping

### 01 - Discovery cohort affinity scores (R)

**Script:** `R/01_discovery_affinity_computation.R`

Computes the M1- and M2-to-T cell affinity scores of the discovery cohort from the per-cell table (cell coordinates and marker calls for CD4, CD8, CD68 and CD163). T cells are CD4+ or CD8+ cells; M1-like cells are CD68+CD163-, M2-like cells CD68+CD163+, and the macrophage pool also includes CD68-CD163+ myeloid cells. Writes the per-cell distances and the per-patient scores. The random pool subsamples use a fixed seed, so reruns give identical results. The published scores came from an earlier run without a fixed seed; script 01b reproduces them exactly.

**Dependencies:** R, data.table, FNN

### 01b - Discovery cohort published scores (R)

**Script:** `R/01b_discovery_scores_from_distances.R`

Recomputes the per-patient discovery scores from the stored per-cell distance files of the published run and compares them with the published score files when these are present.

**Dependencies:** R, data.table

### 01 - Validation cohort affinity scores (Python)

**Script:** `python/01_affinity_computation.py`

Computes the M1- and M2-to-T cell affinity scores of the validation cohort from the per-cell table (cell coordinates and binary marker calls for CD3, CD56, NKp46, CD68, CD163 and CK) and the tissue table (core areas and core-to-patient map). Random subsamples use a fixed seed (42).

**Dependencies:** Python 3.8+, numpy, scipy, pandas, joblib, threadpoolctl, tqdm

---

### 02 - Supplementary Figure S1: Exploratory cohort density vs. affinity (R)

**Script:** `R/02_Fig_S1_exploratory_density_vs_affinity.R`

Demonstrates density dependence of contact fractions and density independence of the affinity metric in the exploratory cohort (n = 10 patients, 17 whole-tissue sections).

| Panel | Description |
|-------|-------------|
| **Suppl. Fig. S1a-d** | Contact fractions (T cells hugged by M1/M2) vs. cell densities |
| **Suppl. Fig. S1e-h** | Affinity scores vs. cell densities |

**Dependencies:** R, data.table, ggplot2, dplyr

---

### 02b - Figure 1c: Exploratory cohort hugging fractions (R)

**Script:** `R/02b_Fig1c_exploratory_hugging_fractions.R`

Computes macrophage-T cell contact fractions in the exploratory cohort: the proportion of M1 and M2 macrophages in contact with T cells, and the reciprocal proportion of T cells in contact with macrophages. The panels were assembled in Illustrator for the final figure.

| Panel | Description |
|-------|-------------|
| **Fig. 1c left** | Proportion of M1 and M2 macrophages adjacent to CD3+ T cells |
| **Fig. 1c right** | Proportion of CD3+ T cells adjacent to M1 and M2 macrophages |

**Dependencies:** R, data.table, ggplot2, dplyr, tidyverse, FNN, reshape2

---

### 03 - Figure 1d: Exploratory cohort ECDF curves (R)

**Script:** `R/03_Fig1d_exploratory_ECDF_curves.R`

Empirical cumulative distribution functions of nearest-neighbour distances from macrophage subtypes to T cells against the random macrophage pool (Kolmogorov-Smirnov test).

| Panel | Description |
|-------|-------------|
| **Fig. 1d upper** | M2-to-T cell distances vs. random macrophage pool |
| **Fig. 1d lower** | M1-to-T cell distances vs. random macrophage pool |

**Dependencies:** R, data.table, reshape2, dplyr, ggplot2

---

### 04 - Figure 1f-g: Discovery cohort Kaplan-Meier survival (R)

**Script:** `R/04_Fig1fg_discovery_KM_survival.R`

Kaplan-Meier analysis of M1- and M2-to-T cell affinity in the discovery cohort, dichotomised at optimal cutpoints (maximally selected rank statistics with bootstrap validation), restricted to patients with both scores.

| Panel | Description |
|-------|-------------|
| **Fig. 1f** | M1-to-T cell affinity, DDFS |
| **Fig. 1g** | M2-to-T cell affinity, DDFS |

**Dependencies:** R, foreign, data.table, survival, survminer, ggpubr, tidyverse

---

### 05 - MSPI integrated analysis (R), main analysis script

**Script:** `R/05_MSPI_integrated_analysis.R`

All MSPI computations across both cohorts.

| Panel | Description |
|-------|-------------|
| **Fig. 1h** | Additive M1 + M2 Cox model and M1 x M2 interaction test |
| **Fig. 2a** | MSPI scatter plot with qMSPI quadrants |
| **Fig. 2c** | C-index comparison of density and spatial models |
| **Fig. 2f** | Early luminal: qMSPI Kaplan-Meier curves, MSPI HRs, C-index |
| **Fig. 2g** | Young early luminal: qMSPI Kaplan-Meier curves, MSPI HRs, C-index |
| **Fig. 3a** | Validation cohort MSPI forest by ER status |
| **Fig. 3c** | Validation cohort ER+ age-stratified MSPI HRs and C-index |
| **Fig. 3d** | Pooled individual patient data analysis |
| **Suppl. Fig. S2** | Discovery cohort affinity scores vs. cell densities |
| CSV/XLSX tables | Univariable and multivariable Cox regression results |

**Dependencies:** R, foreign, data.table, dplyr, tidyr, ggplot2, ggpubr, survival, survminer, openxlsx, boot, gridExtra, scales, forcats, car, flextable, officer, survcomp (Bioconductor)

---

### 06 - Clinical and subtype figures (R)

**Script:** `R/06_clinical_subtype_figures.R`

Molecular subtype analyses, clinical associations and supplementary tables for both cohorts.

| Panel | Description |
|-------|-------------|
| **Fig. 2b** | Cell density and MSPI boxplots by molecular subtype |
| **Fig. 2d** | MSPI HR forest by molecular subtype |
| **Fig. 2e** | Clinical multivariable Cox forest, luminal subgroup |
| **Fig. 2h** | Clinical multivariable Cox forest, young early luminal |
| **Suppl. Table S4** | Discovery cohort: clinicopathological features by qMSPI quadrant |
| **Suppl. Table S5** | Discovery cohort: MSPI HRs by molecular subtype |
| **Suppl. Table S6** | Discovery cohort: clinical multivariable Cox, luminal |
| **Suppl. Table S7** | Discovery cohort: early luminal HRs and C-index |
| **Suppl. Table S8** | Discovery cohort: young early luminal HRs, C-index and clinical multivariable Cox |
| **Suppl. Table S9** | Validation cohort: clinicopathological features by qMSPI quadrant |
| **Suppl. Table S10** | Validation cohort: MSPI HRs by ER status |
| **Suppl. Table S11** | Validation cohort: ER+ age-stratified HRs and C-index |
| **Suppl. Table S12** | Pooled individual patient data analysis |

**Dependencies:** R, foreign, data.table, dplyr, tidyr, ggplot2, ggpubr, survival, survminer, scales, forcats, gridExtra, flextable, officer

---

### 07 - Patient flow diagram and clinicopathological tables (R)

**Script:** `R/07_Fig1a_patient_flow_and_tables.R`

| Panel | Description |
|-------|-------------|
| **Fig. 1a** | Patient flow and cohort structure |
| **Suppl. Table S2** | Discovery cohort clinicopathological characteristics |
| **Suppl. Table S3** | Validation cohort clinicopathological characteristics (all invasive) |
| Additional table | Validation cohort ER-positive clinicopathological characteristics |

**Dependencies:** R, grid, gridExtra, foreign, data.table, dplyr, openxlsx, flextable, officer

---

## Data availability

The clinical and spatial data used in this study are not publicly available due to patient privacy regulations. Requests for data access should be directed to the corresponding author.

## Input data files (not included)

Place the files in the folders below. The scripts accept the original file names listed here; the alternative names used in earlier versions of this repository are also recognised.

**`data/discovery/`**
- `BCY_workingr_832_updated_211024.sav` - clinical data (SPSS)
- `BCY_distances_mean_M1_to_CD3.txt`, `BCY_distances_mean_M2_to_CD3.txt` - per-patient affinity scores
- `BCY_immune_variables_with_CD3_for spatial.txt` - cell densities
- `Merged_BIN_for_spatial.txt`, `Tissue_category_areas.txt` - per-cell table and core-to-patient map (script 01)
- `nnALL_M1_to_CD3.txt`, `nnALL_M2_to_CD3.txt` - per-cell distances of the published run (script 01b)

**`data/validation/`**
- `NK_BRCA_for_affinity.txt` - clinical data and cell densities
- `BREAST_ENCIMA_affinity_CD3e_M1M2_myeloidpool_50um.txt` - per-patient affinity scores
- `Huvuddatabas.xlsx` - extended clinical variables (script 07)

**`data/exploratory/`**
- `nnALL_M1_to_CD3_MMP285.txt`, `nnALL_M2_to_CD3_MMP285.txt` - per-cell distances
- `ia_areas_11062025_images.tsv` - tissue areas
- `hmf.tsv`, `cd3ki67.tsv` - cell-level contact data

## Execution order

1. Affinity scores: `R/01_discovery_affinity_computation.R` (or `R/01b_discovery_scores_from_distances.R` for the published discovery scores) and `python/01_affinity_computation.py` for the validation cohort.
2. `R/02` to `R/07`: each script loads its own inputs, so they can run in any order. The numbers follow the order of results in the manuscript.

## Software environment

Tested with R 4.5.1 (CRAN build for macOS arm64) and Python 3.11. PDF output uses cairo_pdf with a quartz fallback.
