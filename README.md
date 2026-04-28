# Macrophage Spatial Polarity Index (MSPI) in Breast Cancer

Analysis code for the manuscript: *Macrophage-T cell spatial polarity predicts prognosis in luminal breast cancer*.

## Overview

This repository contains the scripts that generate all figures and tables reported in the manuscript. The analysis quantifies a density-independent spatial affinity metric between macrophage subtypes (M1, M2) and T cells in multiplex immunohistochemistry (mIHC) breast cancer tissue sections, and evaluates its prognostic value as the Macrophage Spatial Polarity Index (MSPI).

## Repository structure

```
├── python/
│   └── 01_affinity_computation.py       # Core affinity metric computation
├── R/
│   ├── 02_Fig_S1_exploratory_density_vs_affinity.R
│   ├── 02b_Fig1c_exploratory_hugging_fractions.R
│   ├── 03_Fig1d_exploratory_ECDF_curves.R
│   ├── 04_Fig1fg_discovery_KM_survival.R
│   ├── 05_MSPI_integrated_analysis.R
│   ├── 06_clinical_subtype_figures.R
│   └── 07_Fig1a_patient_flow_and_tables.R
└── README.md
```

## Script descriptions and output mapping

### 01 - Affinity computation (Python)

**Script:** `python/01_affinity_computation.py`

Computes the density-independent spatial affinity score for each tumor. For every T cell (anchor), the script measures:
- 1-nearest-neighbor (1-NN) distance to the nearest macrophage subtype (M1 or M2)
- 1-NN distance to a size-matched random subsample from the total macrophage pool (repeated 5 times, averaged)

The ratio (subtype distance / reference distance) yields the affinity score. Values below 0.5 indicate that the subtype is closer to T cells than expected by chance.

**Dependencies:** Python 3.8+, numpy, scipy (cKDTree), pandas, joblib

**Input:** Cell coordinate files (x, y, phenotype) per tissue core/section

**Output:** Per-patient affinity scores (M1-to-T cell and M2-to-T cell distance files for each cohort)

---

### 02 - Supplementary Figure S1: Exploratory cohort density vs. affinity (R)

**Script:** `R/02_Fig_S1_exploratory_density_vs_affinity.R`

Demonstrates density dependence of contact fractions and density independence of the affinity metric in the exploratory cohort (n = 10 patients, 17 whole-tissue sections).

**Outputs:**

| Panel | Description |
|-------|-------------|
| **Suppl. Fig. S1a-d** | Contact fractions (T cells hugged by M1/M2) vs. cell densities - showing strong density dependence |
| **Suppl. Fig. S1e-h** | Affinity scores vs. cell densities - showing density independence |

**Dependencies:** R 4.0+, data.table, ggplot2, dplyr

---

### 02b - Figure 1c: Exploratory cohort hugging fractions (R)

**Script:** `R/02b_Fig1c_exploratory_hugging_fractions.R`

Computes macrophage-T cell spatial adjacency (contact) fractions in the exploratory cohort. Generates pie/bar charts showing the proportion of macrophages adjacent to T cells, stratified by M1 and M2 subtype, and the reciprocal (T cells adjacent to macrophages).

**Outputs:**

| Panel | Description |
|-------|-------------|
| **Fig. 1c left** | Proportion of M1 and M2 macrophages adjacent to CD3+ T cells |
| **Fig. 1c right** | Proportion of CD3+ T cells adjacent to M1 and M2 macrophages |

**Note:** This script was developed early in the project and reads from raw cell-level data files. The output panels were assembled manually in Illustrator for the final figure.

**Dependencies:** R 4.0+, data.table, ggplot2, dplyr, tidyverse, FNN, reshape2

---

### 03 - Figure 1d: Exploratory cohort ECDF curves (R)

**Script:** `R/03_Fig1d_exploratory_ECDF_curves.R`

Generates empirical cumulative distribution function (ECDF) curves comparing nearest-neighbor distances from macrophage subtypes to T cells versus the random macrophage pool in the exploratory cohort.

**Outputs:**

| Panel | Description |
|-------|-------------|
| **Fig. 1d upper** | M2-to-T cell distances vs. random macrophage pool (KS test) |
| **Fig. 1d lower** | M1-to-T cell distances vs. random macrophage pool (KS test) |

**Dependencies:** R 4.0+, data.table, reshape2, dplyr, ggplot2

---

### 04 - Figure 1f-g: Discovery cohort Kaplan-Meier survival (R)

**Script:** `R/04_Fig1fg_discovery_KM_survival.R`

Kaplan-Meier survival analysis for M1-to-T cell and M2-to-T cell affinity in the discovery cohort. Uses optimal cutpoints determined by maximally selected rank statistics with bootstrap validation. Analyses are restricted to the MSPI-available subset (patients with both M1 and M2 affinity scores).

**Outputs:**

| Panel | Description |
|-------|-------------|
| **Fig. 1f** | KM curves for M1-to-T cell affinity (optimal cutpoint), DDFS |
| **Fig. 1g** | KM curves for M2-to-T cell affinity (optimal cutpoint), DDFS |

**Dependencies:** R 4.0+, foreign, data.table, survival, survminer, ggpubr, tidyverse

---

### 05 - MSPI integrated analysis (R) - main analysis script

**Script:** `R/05_MSPI_integrated_analysis.R`

Core analysis script performing all MSPI-related computations across both cohorts. Score convention: distances are inverted (negated) so that higher values indicate higher affinity. MSPI = M2_score - M1_score (high = suppression dominant = adverse prognosis).

**Outputs:**

| Panel | Description |
|-------|-------------|
| **Fig. 1i** | M1 vs M2 interaction forest plot (Task 3: interaction analysis) |
| **Fig. 2a** | MSPI scatter plot with qMSPI quadrants (Task 2: 4-quadrant analysis) |
| **Fig. 2c** | C-index comparison across density and spatial metrics (Task 5) |
| **Fig. 2f** | KM curves for early luminal with qMSPI quadrants (Task 1) |
| **Fig. 2g** | KM curves for young early luminal with qMSPI quadrants (Task 1) |
| **Fig. 3a** | Validation cohort MSPI forest by ER status (Task 1) |
| **Fig. 3c** | Age-stratified MSPI HR + C-index in validation cohort (Task 1) |
| **Fig. 3d** | Pooled IPD meta-analysis forest (Task 4) |
| **Suppl. Fig. S2** | Discovery cohort density independence scatter plots (affinity vs. cell densities) |
| Numerous CSV/XLSX tables | Univariable and multivariable Cox regression results |

**Dependencies:** R 4.0+, foreign, data.table, dplyr, tidyr, ggplot2, ggpubr, survival, survminer, openxlsx, boot, gridExtra, scales, forcats, car, flextable, officer, survcomp (Bioconductor)

---

### 06 - Clinical and subtype figures (R)

**Script:** `R/06_clinical_subtype_figures.R`

Generates molecular subtype-specific analyses, clinical association tables, and supplementary materials for both cohorts.

**Outputs:**

| Panel | Description |
|-------|-------------|
| **Fig. 2b** | Cell density and MSPI boxplots by molecular subtype |
| **Fig. 2d** | MSPI HR forest by molecular subtype |
| **Fig. 2e** | Clinical multivariable Cox forest - luminal subgroup |
| **Fig. 2h** | Clinical multivariable Cox forest - young early luminal |
| **Suppl. S21** | Discovery cohort MSPI clinical association table (qMSPI quadrants) |
| **Suppl. S22** | Discovery cohort MSPI subtype-specific forest data table |
| **Suppl. S23** | Discovery cohort luminal clinical MV forest data table |
| **Suppl. S24** | Discovery cohort early luminal mini forest and C-index table |
| **Suppl. S25** | Discovery cohort young early luminal mini forest and C-index table |
| **Suppl. S31** | Validation cohort MSPI clinical association table (qMSPI quadrants) |
| **Suppl. S32** | Validation cohort MSPI forest data table |
| **Suppl. S33** | Validation cohort ER+ age-stratified forest and C-index table |
| **Suppl. S34** | Pooled IPD forest data table |

**Dependencies:** R 4.0+, foreign, data.table, dplyr, tidyr, ggplot2, ggpubr, survival, survminer, rstatix, gtsummary, labelled, scales, forcats, gridExtra, flextable, officer, patchwork

---

### 07 - Patient flow diagram and clinicopathological tables (R)

**Script:** `R/07_Fig1a_patient_flow_and_tables.R`

Generates the CONSORT-style patient flow diagram and clinicopathological characterization tables for the MSPI-available subsets of both cohorts.

**Outputs:**

| Panel | Description |
|-------|-------------|
| **Fig. 1a** | Patient flow and cohort structure diagram |
| **Suppl. S12** | Discovery cohort clinicopathological characteristics (MSPI-available subset) |
| **Suppl. S13** | Validation cohort clinicopathological characteristics (all invasive, MSPI-available subset) |
| Additional table | Validation cohort ER-positive clinicopathological characteristics |

**Dependencies:** R 4.0+, grid, gridExtra, foreign, data.table, dplyr, openxlsx, flextable, officer

---

## Data availability

The clinical and spatial data used in this study are not publicly available due to patient privacy regulations. Requests for data access should be directed to the corresponding author.

## Input data files (not included)

The scripts expect the following input files. Paths are configurable at the top of each script.

**Discovery cohort:**
- Clinical data (SPSS .sav format)
- M1-to-T cell affinity scores (tab-delimited text)
- M2-to-T cell affinity scores (tab-delimited text)
- Cell density data (tab-delimited text)

**Validation cohort:**
- Clinical data (tab-delimited text)
- M1 and M2 affinity scores (tab-delimited text)
- Extended clinical variables (Excel .xlsx)

**Exploratory cohort:**
- M1-to-T cell nearest-neighbor distances (tab-delimited text)
- M2-to-T cell nearest-neighbor distances (tab-delimited text)
- Tissue area measurements (tab-delimited text)
- Cell-level contact fraction data (tab-delimited text)

## Execution order

1. **`01_affinity_computation.py`** - Run first to compute affinity scores from raw cell coordinate data (produces input files for downstream R scripts)
2. **`02` through `07`** - R scripts can be run in any order as they independently load required data files. The numbered prefixes reflect the order in which results appear in the manuscript, not execution dependencies.

## Software environment

- **Python:** 3.8+ with numpy, scipy, pandas, joblib
- **R:** 4.0+ with CRAN packages listed per script above
- **Bioconductor:** survcomp (for C-index computation in script 05)
- **System:** Scripts were developed on macOS. PDF output uses cairo_pdf with quartz fallback.
