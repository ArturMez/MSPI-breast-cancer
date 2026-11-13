# Changelog

## 2.0.0 (revised manuscript)

- Discovery cohort T cells: CD8+ cells and CD4+ cells with a lymphocyte-like nucleus, instead of all CD4+ or CD8+
  cells, because CD4 is also expressed by monocytes, macrophages and dendritic cells. New scripts `R/00a` and `R/00b`
  (nucleus measures and nucleus classifier, gradient boosting with xgboost), `python/01_discovery_affinity_nucleus.py` (affinity scores) and
  `R/01c` (T-cell densities). Scripts 04-07 read the new scores and densities and are otherwise unchanged.
- New scripts for the figures of the revised manuscript: `R/08`, `R/10`, `R/11` with `python/09` and `python/12`
  (Figures 1, 2 and 3), `R/13` (Supplementary Figures S1 and S2, densities in cells/mm²) and `R/14` (parity
  sensitivity analysis).
- Figure 2f and g show Kaplan-Meier curves by MSPI quartile (`R/08`). The qMSPI curves of the submitted panels are
  Supplementary Figure S4 (new script `R/15`).
- New script `R/16`: Supplementary Figure S3, the nucleus classifier that decides which CD4+ cells count as T cells
  (nucleus measures, training and application, lymphocyte-like nuclei by cell group, CD4+ cells kept or excluded).
- Supplementary Figures S1 to S4 have lowercase panel letters of one size, as in their legends. S1 and S2 (`R/13`) are
  drawn at the text width of the supplement with printed font sizes, each group of four panels as a 2 x 2 block.
- `R/07_Fig1a_patient_flow_and_tables.R` is now `R/07_Tables_S2_S3_clinicopathological.R`: the tables only, without
  the schematic flow diagram, which printed fixed numbers (the discovery numbers of Figure 1a are drawn by `R/10`).
- Removed: `R/01_discovery_affinity_computation.R` and `R/01b_discovery_scores_from_distances.R` (discovery scores
  with the previous T-cell definition) and `R/02_Fig_S1_exploratory_density_vs_affinity.R` (replaced by `R/13`).
  They remain available in version 1.1.0.

## 1.1.0

- README with script roles, score formula, input files and table numbering. Discovery affinity scripts 01 and 01b.

## 1.0.0

- Analysis scripts of the original submission.
