################################################################################
##
##  Figure S1 - Exploratory Cohort: Density vs Contact & Affinity
##
##  Purpose:
##    Supplementary Figure S1 demonstrates that:
##    - Upper panel: contact (hugging) fractions are correlated with cell
##      densities, confirming density dependence of direct cell-cell contact
##    - Lower panel: the affinity metric is NOT correlated with cell densities,
##      confirming density independence of the spatial affinity approach
##
##  Inputs:
##    - ia_areas_exploratory_images.tsv  (per-sample tissue area + cell counts +
##      hugging counts)
##    - nnALL_M1_to_CD3_exploratory.txt  (M1-to-T cell NND distances)
##    - nnALL_M2_to_CD3_exploratory.txt  (M2-to-T cell NND distances)
##
##  Output:
##    - Individual panel PDFs (for Illustrator assembly)
##    - Combined multi-panel PDF
##    - Summary CSV with all Spearman correlations
##
##  Score convention: affinity metric (delta) = ratio, lower = closer = higher
##    affinity. Density = cell count / AREA.ALL (cells per unit area).
##
##  Design: v6 style (matching discovery cohort survival and MSPI scripts)
##
##  2025-06 Original
##  2026-03-25 Rewritten: v6 design, combined panels, summary CSV
##
################################################################################

## =========================================================================
## 0. CONFIGURATION
## =========================================================================

## Configure paths - adjust these to your local setup
base_dir   <- "./data"
expl_dir   <- file.path(base_dir, "exploratory")

## ia_areas: tissue area + cell counts for exploratory cohort
areas_file <- file.path(expl_dir, "ia_areas_exploratory_images.tsv")

output_dir <- file.path(
  "./output",
  paste0(format(Sys.Date(), "%Y%m%d"), "_Fig_S1_exploratory")
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

cat("Output directory:", output_dir, "\n")

## =========================================================================
## 1. LIBRARIES
## =========================================================================

pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggpubr",
  "scales", "gridExtra", "cowplot"
)

for (p in pkgs) {
  if (!suppressWarnings(require(p, character.only = TRUE, quietly = TRUE))) {
    install.packages(p, repos = "https://cran.rstudio.com/")
    library(p, character.only = TRUE)
  }
}

## =========================================================================
## 2. PLOT SETTINGS (v6 style)
## =========================================================================

BASE_SIZE  <- 14
FONT_FAM   <- "Arial"

## Palette
COL_M1     <- "#5B9BD5"   # blue
COL_M2     <- "#E07B54"   # coral/red
COL_CD3    <- "#6aad8e"   # green (matching old script's "darkgreen")
COL_GREY   <- "#888888"

## Detect quartz fallback for macOS
USE_QUARTZ <- (.Platform$OS.type != "windows") &&
  tryCatch({ cairo_pdf(tempfile()); dev.off(); FALSE },
           error = function(e) TRUE)

theme_set(
  theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
    theme(
      text            = element_text(family = FONT_FAM),
      plot.title      = element_text(size = BASE_SIZE, face = "bold"),
      plot.subtitle   = element_text(size = BASE_SIZE - 2, colour = "grey30"),
      axis.title      = element_text(size = BASE_SIZE - 1),
      axis.text       = element_text(size = BASE_SIZE - 2),
      legend.text     = element_text(size = BASE_SIZE - 2),
      legend.title    = element_text(size = BASE_SIZE - 1),
      strip.text      = element_text(size = BASE_SIZE - 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin     = margin(8, 10, 8, 8)
    )
)

## Helper: save single panel PDF
save_panel <- function(p, fname, w = 3.5, h = 3.5) {
  fpath <- file.path(plots_dir, paste0(fname, ".pdf"))
  if (USE_QUARTZ) {
    quartz(type = "pdf", file = fpath, width = w, height = h,
           family = FONT_FAM, bg = "white")
  } else {
    cairo_pdf(fpath, width = w, height = h, bg = "white")
  }
  print(p)
  dev.off()
  cat("  Saved:", fpath, "\n")
}

## Helper: format p-value nicely
fmt_p <- function(pv) {
  if (pv < 0.001) return(formatC(pv, format = "e", digits = 1))
  if (pv < 0.01)  return(formatC(pv, format = "f", digits = 3))
  formatC(pv, format = "f", digits = 2)
}


## =========================================================================
## 3. LOAD DATA
## =========================================================================

cat("\n=== Loading exploratory cohort data ===\n")

## --- 3a. Tissue area + cell counts (ia_areas) ----------------------------

stopifnot(file.exists(areas_file))
dataTissue <- read.table(areas_file, sep = "\t", header = TRUE, fill = TRUE)

cat("  Loaded ia_areas:", nrow(dataTissue), "samples\n")
cat("  Columns:", paste(colnames(dataTissue), collapse = ", "), "\n")

## Standardize column names (TSV has spaces/dots)
colnames(dataTissue) <- gsub("[. ]+", "_", colnames(dataTissue))
colnames(dataTissue) <- gsub("-", "_", colnames(dataTissue))

## Rename Name -> core_ID
names(dataTissue)[names(dataTissue) == "Name"] <- "core_ID"

## Compute densities (cells / tissue area)
dataTissue$d_CD3      <- (dataTissue$CD3 + dataTissue$HCD3_BY_M1 +
                           dataTissue$HCD3_BY_M2) / dataTissue$AREA_ALL
dataTissue$d_CD3_by_M1 <- dataTissue$HCD3_BY_M1 / dataTissue$AREA_ALL
dataTissue$d_CD3_by_M2 <- dataTissue$HCD3_BY_M2 / dataTissue$AREA_ALL
dataTissue$d_HM1       <- dataTissue$HM1 / dataTissue$AREA_ALL
dataTissue$d_HM2       <- dataTissue$HM2 / dataTissue$AREA_ALL
dataTissue$d_M1        <- dataTissue$M1 / dataTissue$AREA_ALL
dataTissue$d_M2        <- dataTissue$M2 / dataTissue$AREA_ALL

cat("  Density range d_CD3:", round(range(dataTissue$d_CD3, na.rm = TRUE), 6), "\n")
cat("  Density range d_M1:",  round(range(dataTissue$d_M1,  na.rm = TRUE), 6), "\n")
cat("  Density range d_M2:",  round(range(dataTissue$d_M2,  na.rm = TRUE), 6), "\n")

## --- 3b. NND distance files (for affinity scores) -----------------------

nn_m1_file <- file.path(expl_dir, "nnALL_M1_to_CD3_exploratory.txt")
nn_m2_file <- file.path(expl_dir, "nnALL_M2_to_CD3_exploratory.txt")
stopifnot(file.exists(nn_m1_file), file.exists(nn_m2_file))

nn_m1 <- data.table::fread(nn_m1_file, sep = "\t", header = TRUE, fill = TRUE)
nn_m2 <- data.table::fread(nn_m2_file, sep = "\t", header = TRUE, fill = TRUE)

cat("  M1 NND rows:", nrow(nn_m1), "\n")
cat("  M2 NND rows:", nrow(nn_m2), "\n")


## =========================================================================
## 4. COMPUTE AFFINITY SCORES (delta) FROM NND FILES
## =========================================================================

cat("\n=== Computing affinity scores ===\n")

## For each sample: delta = mean(neg) / (mean(POS) + mean(neg))
## POS = distance to M1/M2, neg = distance to random MF pool
## Lower delta = macrophage subtype closer to T cells = higher affinity

compute_delta <- function(nn_df, max_dist = 25) {
  nn_df$POS <- ifelse(nn_df$POS < max_dist, nn_df$POS, NA)
  nn_df$neg <- ifelse(nn_df$neg < max_dist, nn_df$neg, NA)
  nn_df <- na.omit(as.data.frame(nn_df))

  aff <- nn_df %>%
    dplyr::group_by(core_ID) %>%
    dplyr::summarise(
      delta = mean(neg, na.rm = TRUE) /
        (mean(POS, na.rm = TRUE) + mean(neg, na.rm = TRUE)),
      .groups = "drop"
    )
  return(aff)
}

delta_m1 <- compute_delta(nn_m1)
names(delta_m1)[2] <- "delta_M1"

delta_m2 <- compute_delta(nn_m2)
names(delta_m2)[2] <- "delta_M2"


## =========================================================================
## 5. MERGE ALL DATA
## =========================================================================

sample_df <- dataTissue %>%
  dplyr::left_join(delta_m1, by = "core_ID") %>%
  dplyr::left_join(delta_m2, by = "core_ID")

cat("  Merged samples:", nrow(sample_df), "\n")
cat("  With M1 affinity:", sum(!is.na(sample_df$delta_M1)), "\n")
cat("  With M2 affinity:", sum(!is.na(sample_df$delta_M2)), "\n")

## Save summary table
write.csv(sample_df, file.path(output_dir, "exploratory_sample_summary.csv"),
          row.names = FALSE)
cat("  Saved sample summary CSV\n")


## =========================================================================
## 6. CORRELATION SCATTER PLOT FUNCTION
## =========================================================================

make_scatter <- function(df, x_var, y_var, x_lab, y_lab,
                         point_col = COL_GREY, title_text = NULL) {

  df_plot <- df[!is.na(df[[x_var]]) & !is.na(df[[y_var]]), ]

  sp <- cor.test(df_plot[[x_var]], df_plot[[y_var]], method = "spearman")
  rho_val  <- round(sp$estimate, 2)
  p_val    <- fmt_p(sp$p.value)

  sub_text <- paste0("Spearman \u03C1 = ", rho_val, ", p = ", p_val)

  p <- ggplot(df_plot, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(colour = point_col, alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7,
                colour = "firebrick", fill = "firebrick", alpha = 0.15) +
    ## Fewer x-axis ticks, plain labels (no scientific notation),
    ## each number as a single text element for clean Illustrator editing
    scale_x_continuous(n.breaks = 4,
                       labels = function(x) formatC(x, format = "g", digits = 3)) +
    scale_y_continuous(n.breaks = 5,
                       labels = function(x) formatC(x, format = "g", digits = 3)) +
    labs(
      x        = x_lab,
      y        = y_lab,
      title    = title_text,
      subtitle = sub_text
    )

  return(list(plot = p, rho = sp$estimate, p = sp$p.value,
              n = nrow(df_plot), x = x_var, y = y_var))
}


## =========================================================================
## 7. UPPER PANEL: CONTACT DENSITIES vs CELL DENSITIES
##    (Density-DEPENDENT - expected to be correlated)
## =========================================================================

cat("\n=== Upper panel: contact densities vs cell densities ===\n")

upper_plots <- list()

## 7a. T cells hugged by M1 (density) vs CD3 density
upper_plots[["dCD3byM1_vs_dCD3"]] <- make_scatter(
  sample_df, "d_CD3", "d_CD3_by_M1",
  x_lab = "CD3 density (cells/area)",
  y_lab = "T cells hugged by M1 (density)",
  point_col = COL_M1,
  title_text = "Contact by M1 vs T cell density"
)

## 7b. T cells hugged by M1 (density) vs M1 density
upper_plots[["dCD3byM1_vs_dM1"]] <- make_scatter(
  sample_df, "d_M1", "d_CD3_by_M1",
  x_lab = "M1 density (cells/area)",
  y_lab = "T cells hugged by M1 (density)",
  point_col = COL_M1,
  title_text = "Contact by M1 vs M1 density"
)

## 7c. T cells hugged by M2 (density) vs CD3 density
upper_plots[["dCD3byM2_vs_dCD3"]] <- make_scatter(
  sample_df, "d_CD3", "d_CD3_by_M2",
  x_lab = "CD3 density (cells/area)",
  y_lab = "T cells hugged by M2 (density)",
  point_col = COL_M2,
  title_text = "Contact by M2 vs T cell density"
)

## 7d. T cells hugged by M2 (density) vs M2 density
upper_plots[["dCD3byM2_vs_dM2"]] <- make_scatter(
  sample_df, "d_M2", "d_CD3_by_M2",
  x_lab = "M2 density (cells/area)",
  y_lab = "T cells hugged by M2 (density)",
  point_col = COL_M2,
  title_text = "Contact by M2 vs M2 density"
)

## 7e. Hugging M1 macrophages (density) vs M1 density
upper_plots[["dHM1_vs_dM1"]] <- make_scatter(
  sample_df, "d_M1", "d_HM1",
  x_lab = "M1 density (cells/area)",
  y_lab = "Hugging M1 density",
  point_col = COL_M1,
  title_text = "Hugging M1 vs M1 density"
)

## 7f. Hugging M2 macrophages (density) vs M2 density
upper_plots[["dHM2_vs_dM2"]] <- make_scatter(
  sample_df, "d_M2", "d_HM2",
  x_lab = "M2 density (cells/area)",
  y_lab = "Hugging M2 density",
  point_col = COL_M2,
  title_text = "Hugging M2 vs M2 density"
)


## =========================================================================
## 8. LOWER PANEL: AFFINITY SCORES vs CELL DENSITIES
##    (Density-INDEPENDENT - expected NOT to be correlated)
## =========================================================================

cat("\n=== Lower panel: affinity scores vs densities ===\n")

lower_plots <- list()

## 8a. M1 affinity vs CD3 density
lower_plots[["affM1_vs_dCD3"]] <- make_scatter(
  sample_df, "d_CD3", "delta_M1",
  x_lab = "CD3 density (cells/area)",
  y_lab = "M1-to-T cell affinity",
  point_col = COL_M1,
  title_text = "M1 affinity vs T cell density"
)

## 8b. M1 affinity vs M1 density
lower_plots[["affM1_vs_dM1"]] <- make_scatter(
  sample_df, "d_M1", "delta_M1",
  x_lab = "M1 density (cells/area)",
  y_lab = "M1-to-T cell affinity",
  point_col = COL_M1,
  title_text = "M1 affinity vs M1 density"
)

## 8c. M2 affinity vs CD3 density
lower_plots[["affM2_vs_dCD3"]] <- make_scatter(
  sample_df, "d_CD3", "delta_M2",
  x_lab = "CD3 density (cells/area)",
  y_lab = "M2-to-T cell affinity",
  point_col = COL_M2,
  title_text = "M2 affinity vs T cell density"
)

## 8d. M2 affinity vs M2 density
lower_plots[["affM2_vs_dM2"]] <- make_scatter(
  sample_df, "d_M2", "delta_M2",
  x_lab = "M2 density (cells/area)",
  y_lab = "M2-to-T cell affinity",
  point_col = COL_M2,
  title_text = "M2 affinity vs M2 density"
)


## =========================================================================
## 9. SAVE INDIVIDUAL PANELS
## =========================================================================

cat("\n=== Saving individual panel PDFs ===\n")

for (nm in names(upper_plots)) {
  save_panel(upper_plots[[nm]]$plot, paste0("S1_upper_", nm))
}

for (nm in names(lower_plots)) {
  save_panel(lower_plots[[nm]]$plot, paste0("S1_lower_", nm))
}


## =========================================================================
## 10. COMBINED MULTI-PANEL FIGURE
## =========================================================================

cat("\n=== Assembling combined panel ===\n")

## Upper row: M1 contact vs CD3, M1 contact vs M1, M2 contact vs CD3, M2 contact vs M2
## Lower row: M1 affinity vs CD3, M1 affinity vs M1, M2 affinity vs CD3, M2 affinity vs M2

upper_row <- cowplot::plot_grid(
  upper_plots[["dCD3byM1_vs_dCD3"]]$plot + ggtitle(NULL),
  upper_plots[["dCD3byM1_vs_dM1"]]$plot  + ggtitle(NULL),
  upper_plots[["dCD3byM2_vs_dCD3"]]$plot + ggtitle(NULL),
  upper_plots[["dCD3byM2_vs_dM2"]]$plot  + ggtitle(NULL),
  ncol = 4, labels = c("a", "b", "c", "d"),
  label_size = BASE_SIZE + 2, label_fontfamily = FONT_FAM
)

lower_row <- cowplot::plot_grid(
  lower_plots[["affM1_vs_dCD3"]]$plot + ggtitle(NULL),
  lower_plots[["affM1_vs_dM1"]]$plot  + ggtitle(NULL),
  lower_plots[["affM2_vs_dCD3"]]$plot + ggtitle(NULL),
  lower_plots[["affM2_vs_dM2"]]$plot  + ggtitle(NULL),
  ncol = 4, labels = c("e", "f", "g", "h"),
  label_size = BASE_SIZE + 2, label_fontfamily = FONT_FAM
)

## Row labels
upper_label <- cowplot::ggdraw() +
  cowplot::draw_label("Contact fractions vs cell density (density-dependent)",
                      fontface = "bold", size = BASE_SIZE,
                      fontfamily = FONT_FAM, hjust = 0, x = 0.02)

lower_label <- cowplot::ggdraw() +
  cowplot::draw_label("Affinity score vs cell density (density-independent)",
                      fontface = "bold", size = BASE_SIZE,
                      fontfamily = FONT_FAM, hjust = 0, x = 0.02)

combined <- cowplot::plot_grid(
  upper_label, upper_row,
  lower_label, lower_row,
  ncol = 1, rel_heights = c(0.06, 1, 0.06, 1)
)

## Save combined figure
fpath <- file.path(plots_dir, "Fig_S1_combined.pdf")
if (USE_QUARTZ) {
  quartz(type = "pdf", file = fpath, width = 14, height = 8,
         family = FONT_FAM, bg = "white")
} else {
  cairo_pdf(fpath, width = 14, height = 8, bg = "white")
}
print(combined)
dev.off()
cat("  Saved:", fpath, "\n")


## =========================================================================
## 11. SUMMARY TABLE
## =========================================================================

cat("\n=== Correlation summary ===\n")

all_results <- c(upper_plots, lower_plots)

summary_df <- data.frame(
  Panel    = names(all_results),
  x_var    = sapply(all_results, function(r) r$x),
  y_var    = sapply(all_results, function(r) r$y),
  n        = sapply(all_results, function(r) r$n),
  Spearman_rho = sapply(all_results, function(r) round(r$rho, 3)),
  p_value  = sapply(all_results, function(r) signif(r$p, 3)),
  stringsAsFactors = FALSE
)
rownames(summary_df) <- NULL

print(summary_df)

write.csv(summary_df, file.path(output_dir, "S1_correlation_summary.csv"),
          row.names = FALSE)
cat("  Saved: S1_correlation_summary.csv\n")


cat("\n============================================================\n")
cat("  Done. All outputs in:", output_dir, "\n")
cat("============================================================\n")
