# Macrophage Spatial Polarity Index (MSPI) in Breast Cancer

Analysis code for the manuscript *Macrophage spatial polarity to T cells predicts prognosis in young women with luminal breast cancer*.

Version 2.0.0 accompanies the revised manuscript. In the discovery cohort, T cells are now CD8+ cells and CD4+ cells with a lymphocyte-like nucleus (see below and `CHANGELOG.md`). Version 1.1.0 holds the code of the original submission.

## Overview

This repository contains the scripts that compute the spatial affinity scores and generate the figures and tables reported in the manuscript. The analysis quantifies a density-independent spatial affinity metric between macrophage subtypes (M1-like, M2-like) and T cells in multiplex immunohistochemistry (mIHC) breast cancer tissue, and evaluates its prognostic value as the Macrophage Spatial Polarity Index (MSPI).

## The affinity score

For every T cell (anchor), the scripts record the nearest-neighbour distance (NND) to the nearest M1-like (or M2-like) macrophage, and to the nearest cell of a random subsample of the total macrophage pool drawn to the same size (repeated 5 times). A nearer neighbour closer than 1 um is treated as the same position and skipped. Per patient, NNDs at or above a distance cap are excluded, separately for the subset and the pool, and

    delta = mean(subset NND) / (mean(pool NND) + mean(subset NND))

delta = 0.5 means no preference, and values below 0.5 mean that the subset lies closer to T cells than the macrophage pool as a whole. Script 05 negates delta, so higher scores mean higher affinity, and defines MSPI = M2 score - M1 score.

## T cells in the discovery cohort

The discovery cohort was stained for CD4, CD8, CD68 and CD163 without CD3. Because CD4 is also expressed by monocytes, macrophages and dendritic cells, a CD4+ cell counts as a T cell only when its nucleus is classified as lymphocyte-like (script 00b). The classifier is gradient boosting of decision trees (xgboost) on nucleus area, axis ratio, compactness and minor axis, cell area and the mean, minimum, maximum, standard deviation and total nuclear DAPI intensity from the inForm segmentation. It is trained on CD8+ cells without CD4, CD68 or CD163 signal against CD68+ or CD163+ cells without CD4 or CD8 signal in 70% of the TMA cores (tree depth and number of rounds chosen by early stopping on 15% of these cores) and reaches an area under the curve of 0.80 in the other 30%. T cells are all CD8+ cells and the CD4+ cells with a predicted probability above 0.5. The same definition is used for the affinity scores (python/01) and the T-cell densities (R/01c). The exploratory and validation cohorts use CD3.

## Repository structure

```
├── python/
│   ├── 01_discovery_affinity_nucleus.py       # Discovery cohort: affinity scores
│   ├── 01_affinity_computation.py             # Validation cohort: affinity scores
│   ├── 09_Fig2_compose.py                     # Figure 2: data panels onto the figure page
│   └── 12_Fig1_Fig3_compose.py                # Figures 1 and 3: data panels onto the figure pages
├── R/
│   ├── 00a_discovery_nucleus_features.R       # Discovery cohort: nucleus measures for every cell
│   ├── 00b_discovery_nucleus_classifier.R     # Discovery cohort: nucleus classifier
│   ├── 01c_discovery_T_cell_density_nucleus.R # Discovery cohort: T-cell densities
│   ├── 02b_Fig1c_exploratory_hugging_fractions.R
│   ├── 03_Fig1d_exploratory_ECDF_curves.R
│   ├── 04_Fig1fg_discovery_KM_survival.R
│   ├── 05_MSPI_integrated_analysis.R
│   ├── 06_clinical_subtype_figures.R
│   ├── 07_Tables_S2_S3_clinicopathological.R
│   ├── 08_Fig2_draw.R
│   ├── 10_Fig1_draw.R
│   ├── 11_Fig3d_draw.R
│   ├── 13_Supplementary_Figures_S1_S2.R
│   ├── 14_parity_sensitivity.R
│   ├── 15_Supplementary_Figure_S4_qMSPI.R
│   └── 16_Supplementary_Figure_S3_nucleus_classifier.R
├── assets/                                    # Figure pages without the data panels (scripts 09 and 12)
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
└── README.md
```

All scripts are run from the repository root. They read inputs from `./data/` and write results to `./output/`.

## Script descriptions and output mapping

### 00a, 00b - Discovery cohort nucleus classifier (R)

**Scripts:** `R/00a_discovery_nucleus_features.R`, `R/00b_discovery_nucleus_classifier.R`

00a attaches the inForm nucleus measures (area, shape, DAPI) to every analysis cell. 00b fits the nucleus classifier described above (fixed seed) and writes the probability of a lymphocyte-like nucleus for every cell.

**Dependencies:** R, data.table, xgboost (3.2)

### 01 - Discovery cohort affinity scores (Python)

**Script:** `python/01_discovery_affinity_nucleus.py`

Computes the M1- and M2-to-T cell affinity scores of the discovery cohort. T cells as above. M1-like cells are CD68+CD163-, M2-like cells CD68+CD163+, and the macrophage pool also includes CD68-CD163+ myeloid cells, in each case without the T cells. The random pool subsamples are seeded per core, so reruns give identical scores.

**Dependencies:** Python 3.8+, numpy, scipy, pandas

### 01c - Discovery cohort T-cell densities (R)

**Script:** `R/01c_discovery_T_cell_density_nucleus.R`

T-cell density (cells/mm²) per patient in all tissue, tumour and stroma, with the T cells defined above.

**Dependencies:** R, data.table

### 01 - Validation cohort affinity scores (Python)

**Script:** `python/01_affinity_computation.py`

Computes the M1- and M2-to-T cell affinity scores of the validation cohort from the per-cell table (cell coordinates and binary marker calls for CD3, CD56, NKp46, CD68, CD163 and CK) and the tissue table (core areas and core-to-patient map). Random subsamples use a fixed seed (42).

**Dependencies:** Python 3.8+, numpy, scipy, pandas, joblib, threadpoolctl, tqdm

---

### 02b - Figure 1c: Exploratory cohort hugging fractions (R)

**Script:** `R/02b_Fig1c_exploratory_hugging_fractions.R`

Computes macrophage-T cell contact fractions in the exploratory cohort: the proportion of M1 and M2 macrophages in contact with T cells, and the reciprocal proportion of T cells in contact with macrophages.

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

### 04 - Discovery cohort Kaplan-Meier survival (R)

**Script:** `R/04_Fig1fg_discovery_KM_survival.R`

Kaplan-Meier analysis of M1- and M2-to-T cell affinity in the discovery cohort, dichotomised at optimal cutpoints (maximally selected rank statistics with bootstrap validation), restricted to patients with both scores. The panels of Figure 1f-g are drawn by script 10 with the same computation.

**Dependencies:** R, foreign, data.table, survival, survminer, ggpubr, tidyverse

---

### 05 - MSPI integrated analysis (R), main analysis script

**Script:** `R/05_MSPI_integrated_analysis.R`

All MSPI computations across both cohorts: MSPI and qMSPI, univariable and density-adjusted Cox models by subgroup, C-indices, clinical multivariable models and the pooled individual patient data analysis. Its tables are the source of Figures 1-3 (drawn by scripts 08, 10 and 11) and of several supplementary tables.

**Dependencies:** R, foreign, data.table, dplyr, tidyr, ggplot2, ggpubr, survival, survminer, openxlsx, boot, gridExtra, scales, forcats, car, flextable, officer, survcomp (Bioconductor)

---

### 06 - Clinical and subtype analyses (R)

**Script:** `R/06_clinical_subtype_figures.R`

Molecular subtype analyses, clinical associations, the additive M1 + M2 Cox model with the interaction test, and supplementary tables for both cohorts.

| Output | Description |
|-------|-------------|
| **Fig. 1h** | Additive M1 + M2 Cox model and M1 x M2 interaction test (drawn by script 10) |
| **Fig. 2b, 2d, 2e, 2h** | Subtype boxplots and Dunn's tests, MSPI HR by subtype, clinical multivariable models (drawn by script 08) |
| **Suppl. Table S4** | Discovery cohort: clinicopathological features by qMSPI quadrant |
| **Suppl. Table S5** | Discovery cohort: MSPI HRs by molecular subtype |
| **Suppl. Tables S6-S8** | Discovery cohort: clinical multivariable Cox, early and young early luminal HRs and C-index |
| **Suppl. Tables S9-S11** | Validation cohort: qMSPI quadrants, MSPI HRs by ER status, ER+ age-stratified HRs and C-index |
| **Suppl. Table S12** | Pooled individual patient data analysis |

**Dependencies:** R, foreign, data.table, dplyr, tidyr, ggplot2, ggpubr, survival, survminer, scales, forcats, gridExtra, flextable, officer, rstatix, car

---

### 07 - Clinicopathological tables (R)

**Script:** `R/07_Tables_S2_S3_clinicopathological.R`

| Output | Description |
|-------|-------------|
| **Suppl. Table S2** | Discovery cohort clinicopathological characteristics |
| **Suppl. Table S3** | Validation cohort clinicopathological characteristics (all invasive) |
| Additional table | Validation cohort ER-positive clinicopathological characteristics |

The discovery numbers of Figure 1a are drawn by script 10.

**Dependencies:** R, foreign, data.table, dplyr, openxlsx, flextable, officer

---

### 08, 10, 11 and 09, 12 - Figures 1, 2 and 3 (R and Python)

**Scripts:** `R/08_Fig2_draw.R`, `R/10_Fig1_draw.R`, `R/11_Fig3d_draw.R`, `python/09_Fig2_compose.py`, `python/12_Fig1_Fig3_compose.py`

The data panels of the figures, drawn from the outputs of scripts 05 and 06 (and the Kaplan-Meier, boxplot and pie computations they repeat) at the positions, scales, line widths, marker sizes, colours and fonts of the figure layout. Scripts 09 and 12 place them onto the figure pages in `assets/`, which hold everything that does not depend on the data (panel letters, labels, images and schematics, panels of the exploratory and validation cohorts).

| Output | Description |
|-------|-------------|
| **Fig. 1** | Panel a (discovery patient numbers), f and g (Kaplan-Meier by M1- and M2-to-T affinity), h (additive model) |
| **Fig. 2** | All data panels a-h (f and g: Kaplan-Meier curves by MSPI quartile, cut within each subgroup) |
| **Fig. 3** | Panel d (pooled analysis) |

**Dependencies:** R, foreign, data.table, survival, survminer, car, rstatix, grid, systemfonts, labeling; Python 3.8+, PyMuPDF

---

### 13 - Supplementary Figures S1 and S2 (R)

**Script:** `R/13_Supplementary_Figures_S1_S2.R`

| Output | Description |
|-------|-------------|
| **Suppl. Fig. S1** | Exploratory cohort: contact fractions and affinity scores against cell densities |
| **Suppl. Fig. S2** | Discovery cohort: affinity scores against cell densities |

Both figures are drawn at the text width of the supplement (152 mm), each group of four panels as a 2 x 2 block.

**Dependencies:** R, data.table, dplyr, ggplot2, patchwork, foreign

---

### 14 - Parity sensitivity analysis (R)

**Script:** `R/14_parity_sensitivity.R`

Parity in the discovery cohort, MSPI in parous and nulliparous patients, and the clinical multivariable models of the luminal and young early luminal subgroups without and with parity.

**Dependencies:** R, data.table, survival, foreign

---

### 15 - Supplementary Figure S4 (R)

**Script:** `R/15_Supplementary_Figure_S4_qMSPI.R`

| Output | Description |
|-------|-------------|
| **Suppl. Fig. S4** | Discovery cohort: Kaplan-Meier curves of DDFS by qMSPI phenotype in early-stage and young early-stage luminal patients (quadrants by the subgroup medians, log-rank p for Q4 vs Q1) |

**Dependencies:** R, data.table, survival, ggplot2, patchwork, foreign

---

### 16 - Supplementary Figure S3 (R)

**Script:** `R/16_Supplementary_Figure_S3_nucleus_classifier.R`

| Output | Description |
|-------|-------------|
| **Suppl. Fig. S3** | Discovery cohort: nucleus-supported identification of CD4+ T cells. Nucleus measures of the classifier with DAPI intensity and axis ratio in the stroma by cell group, schematic of the classifier (training on CD8+ T cells against CD68+ or CD163+ cells, application to CD4+ cells), lymphocyte-like nuclei by cell group, and CD4+ cells counted as T cells or excluded |

Reads the outputs of `R/00a` and `R/00b` in `./output/discovery_nucleus/`. Numbers behind the figure: `S3_numbers.tsv` and `S3_totals.tsv`.

**Dependencies:** R, data.table, ggplot2, patchwork, ggtext, foreign

---

## Data availability

The clinical and spatial data used in this study are not publicly available due to patient privacy regulations. Requests for data access should be directed to the corresponding author.

## Input data files (not included)

**`data/discovery/`**
- `BCY_workingr_832_updated_211024.sav` - clinical data (SPSS)
- `Merged_BIN_for_spatial.txt` - per-cell table (coordinates and marker calls)
- `Tissue_category_areas.txt` - tissue areas per core and compartment
- `BCY_immune_variables_with_CD3_for spatial.txt` - cell densities
- the inForm cell segmentation export (`BCY_merged_cell_seg_data.txt`, `BCY_merged_additional_cell_seg_data.txt`), in the folder named by the environment variable `MSPI_INFORM_DIR` (default `./data/discovery/inform/`)
- after step 1 below: `discovery_distances_mean_M1_to_CD3_nucleus.txt`, `discovery_distances_mean_M2_to_CD3_nucleus.txt` and `discovery_immune_variables_with_CD3_nucleus.txt`, copied from `./output/discovery_affinity_nucleus/`

**`data/validation/`**
- `NK_BRCA_for_affinity.txt` - clinical data and cell densities
- `BREAST_ENCIMA_affinity_CD3e_M1M2_myeloidpool_50um.txt` - per-patient affinity scores
- `Huvuddatabas.xlsx` - extended clinical variables (script 07)

**`data/exploratory/`**
- `nnALL_M1_to_CD3_MMP285.txt`, `nnALL_M2_to_CD3_MMP285.txt` - per-cell distances
- `ia_areas_11062025_images.tsv` - tissue areas and cell counts
- `hmf.tsv`, `cd3ki67.tsv` - cell-level contact data

## Execution order

1. Discovery cohort: `R/00a`, `R/00b`, `python/01_discovery_affinity_nucleus.py` and `R/01c` (about 5 minutes), then copy the three files named above to `./data/discovery/`. Validation cohort: `python/01_affinity_computation.py`.
2. `R/02b` to `R/07` and `R/13` to `R/16`: each script loads its own inputs, so they can run in any order (script 05 takes about 7 minutes).
3. Figures: `R/08_Fig2_draw.R`, `R/10_Fig1_draw.R` and `R/11_Fig3d_draw.R`, then `python/09_Fig2_compose.py` and `python/12_Fig1_Fig3_compose.py` (Figures 1-3 in `./output/fig1/`, `./output/fig2/` and `./output/fig3/`).

## Software environment

Tested with R 4.5.1 (CRAN build for macOS, run with `LANG=en_US.UTF-8`) and Python 3.11. PDF output uses cairo_pdf with a quartz fallback.
