################################################################################
##
##  Discovery and Validation Cohort MSPI — SUBTYPE FOREST PLOTS, BOXPLOTS & CLINICAL ASSOCIATIONS
##  VERSION 2
##
##  Purpose:
##    1. Forest plot of MSPI HR per SD by molecular subtype (Discovery cohort)
##    2. Boxplots of MSPI by molecular subtype (Discovery cohort)
##    3. Clinical associations: MSPI vs clinicopathological features (Discovery cohort)
##    4. Cell density boxplots by subtype (Discovery cohort)
##    5. Validation cohort: forest (ER+/ER-), boxplots, clinical associations
##
##  v2 changes:
##    - CHANGE 1: Clinical tables split by 4-quadrant MSPI groups (Q1 Best, Q2 Mixed,
##                Q3 Interm, Q4 Worst) instead of median split, for both Discovery cohort and Validation cohort.
##                Keep median-split tables as secondary output for backward compatibility.
##    - CHANGE 2: All tables output as .docx (primary) via flextable + officer,
##                with .pdf fallback via gridExtra. Consistent Arial styling throughout.
##    - CHANGE 3: Ensure all forest plots have matching CSVs + docx/PDF tables.
##    - CHANGE 4: Supplementary tables - Table S2 (Discovery cohort clinicopathological
##                summary) and Table S3 (Validation cohort clinicopathological summary),
##                matching format of existing manual Table S1 (exploratory cohort).
##
##  Discovery cohort inputs:
##    - discovery_workingr_832_updated_211024.sav (SPSS clinical data)
##    - discovery_distances_mean_M1_to_CD3.txt
##    - discovery_distances_mean_M2_to_CD3.txt
##    - discovery_immune_variables_with_CD3_for spatial.txt
##
##  Validation cohort inputs:
##    - NK_BRCA_for_affinity.txt (clinical data)
##    - BREAST_validation_affinity_CD3e_M1M2_myeloidpool_50um.txt (affinity scores)
##
##  Cell densities: CD3, M1, M2 (raw, no log transform) in both cohorts
##
##  Output: PDF figures + summary CSV + forest table PDFs in output folder
##
##  Score convention: INVERTED (negated distance) — higher = more affinity
##  MSPI = M2_score - M1_score  (high = suppression dominant = adverse)
##
##  2026-03-18  Artur Mezheyeuski
##  2026-03-20  Added Validation cohort analyses
##  2026-03-25  v2: quadrant tables, forest table PDFs, enhanced outputs
##
################################################################################

## ═══════════════════════════════════════════════════════════════════════════
## 0. CONFIGURATION
## ═══════════════════════════════════════════════════════════════════════════

disc_base <- "./data/discovery"

## Input files may carry the names used in this repository or their original names
## (as shared on request); the first file found is used.
find_input <- function(dir, ...) {
  for (f in c(...)) {
    p <- file.path(dir, f)
    if (file.exists(p)) return(p)
  }
  stop("Cannot find any of: ", paste(c(...), collapse = ", "), " in ", dir)
}

output_dir <- file.path(
  "./output",
  paste0(format(Sys.Date(), "%Y%m%d"), "_MSPI_clinical")
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", output_dir, "\n")

## ═══════════════════════════════════════════════════════════════════════════
## 1. LIBRARIES
## ═══════════════════════════════════════════════════════════════════════════

pkgs <- c(
  "foreign", "data.table", "dplyr", "tidyr", "ggplot2", "ggpubr",
  "survival", "survminer", "rstatix", "gtsummary", "labelled",
  "scales", "forcats", "gridExtra", "flextable", "officer", "patchwork"
)

for (p in pkgs) {
  if (!suppressWarnings(require(p, character.only = TRUE, quietly = TRUE))) {
    install.packages(p, repos = "https://cran.rstudio.com/")
    library(p, character.only = TRUE)
  }
}

## ═══════════════════════════════════════════════════════════════════════════
## 2. LOAD AND MERGE DATA
## ═══════════════════════════════════════════════════════════════════════════

cat("Loading Discovery cohort clinical data...\n")

disc_spss <- find_input(disc_base, "discovery_workingr_832_updated_211024.sav",
                        "Matching to clinical data/discovery_workingr_832_updated_211024.sav",
                        "BCY_workingr_832_updated_211024.sav",
                        "Matching to clinical data/BCY_workingr_832_updated_211024.sav")
stopifnot(file.exists(disc_spss))

disc_raw <- foreign::read.spss(
  disc_spss,
  use.value.labels = TRUE,
  to.data.frame    = TRUE,
  max.value.labels = Inf,
  trim.factor.names = TRUE,
  trim_values       = TRUE,
  reencode          = NA
)

## Survival endpoint: DDFS
disc_raw$status_num <- ifelse(
  as.character(disc_raw$DDFScens2021) == "Distant recurrence" |
    as.character(disc_raw$BCcens2021) == "Breast cancer death",
  1L, 0L
)
disc_raw$DDFStime2021 <- as.numeric(disc_raw$DDFStime2021)

## Rename data column from input file to generic name
## The original clinical file calls the subtype variable Subtype_BCY21
for (old in c("Subtype_BCY21", "Subtype_Discovery21")) {
  if (!"Subtype" %in% names(disc_raw) && old %in% names(disc_raw)) {
    names(disc_raw)[names(disc_raw) == old] <- "Subtype"
  }
}

## Remove unclassified
disc_raw <- disc_raw[as.character(disc_raw$Subtype) != "Unclassified", ]
disc_raw$Subtype <- droplevels(disc_raw$Subtype)

cat("  Clinical rows after cleanup:", nrow(disc_raw), "\n")

## Spatial distance files
find_file <- function(...) {
  fnames <- c(...)
  candidates <- c(
    file.path(disc_base, fnames),
    file.path(disc_base, "Discovery cohort", fnames),
    file.path(disc_base, "analysis", fnames)
  )
  for (f in candidates) if (file.exists(f)) return(f)
  stop(paste("Cannot find:", paste(fnames, collapse = " / ")))
}

m1_raw  <- data.table::fread(find_file("discovery_distances_mean_M1_to_CD3.txt", "BCY_distances_mean_M1_to_CD3.txt"), data.table = FALSE)
m2_raw  <- data.table::fread(find_file("discovery_distances_mean_M2_to_CD3.txt", "BCY_distances_mean_M2_to_CD3.txt"), data.table = FALSE)
cd3_raw <- data.table::fread(find_file("discovery_immune_variables_with_CD3_for spatial.txt", "BCY_immune_variables_with_CD3_for spatial.txt"),
                              data.table = FALSE)

colnames(m1_raw)[colnames(m1_raw) == "delta_CD3"] <- "M1_to_CD3_distance"
colnames(m2_raw)[colnames(m2_raw) == "delta_CD3"] <- "M2_to_CD3_distance"
cd3_raw <- cd3_raw[, 1:2]
colnames(cd3_raw) <- c("AnnUnit", "CD3_density")

## Merge
disc <- disc_raw %>%
  dplyr::left_join(m1_raw[, c("AnnUnit", "M1_to_CD3_distance")], by = "AnnUnit") %>%
  dplyr::left_join(m2_raw[, c("AnnUnit", "M2_to_CD3_distance")], by = "AnnUnit") %>%
  dplyr::left_join(cd3_raw, by = "AnnUnit")

## Score inversion (higher = more affinity)
disc$M1_score <- -disc$M1_to_CD3_distance
disc$M2_score <- -disc$M2_to_CD3_distance

## MSPI = M2 - M1 (high = suppression dominant)
disc$MSPI <- disc$M2_score - disc$M1_score

## Filter to analysis-ready (both scores + valid survival)
disc_s4 <- disc %>%
  dplyr::filter(
    !is.na(M1_score), !is.na(M2_score),
    !is.na(DDFStime2021), DDFStime2021 > 0
  )

cat("  Analysis-ready (S4): n =", nrow(disc_s4),
    "| DDFS events:", sum(disc_s4$status_num), "\n")

## Cell densities: M1, M2, CD3
if ("M1" %in% names(disc_s4) && is.numeric(disc_s4[["M1"]])) {
  disc_s4$M1_density <- disc_s4[["M1"]]
  cat("  M1 density: using 'M1' column\n")
} else {
  cat("  WARNING: no M1 density column found\n")
  disc_s4$M1_density <- NA
}

if ("M2" %in% names(disc_s4) && is.numeric(disc_s4[["M2"]])) {
  disc_s4$M2_density <- disc_s4[["M2"]]
  cat("  M2 density: using 'M2' column\n")
} else {
  cat("  WARNING: no M2 density column found\n")
  disc_s4$M2_density <- NA
}

cat("  CD3 density available:", sum(!is.na(disc_s4$CD3_density) & disc_s4$CD3_density > 0), "cases\n")
cat("  M1 density available:",  sum(!is.na(disc_s4$M1_density) & disc_s4$M1_density > 0), "cases\n")
cat("  M2 density available:",  sum(!is.na(disc_s4$M2_density) & disc_s4$M2_density > 0), "cases\n")


## ═══════════════════════════════════════════════════════════════════════════
## 3. MOLECULAR SUBTYPE DEFINITIONS
## ═══════════════════════════════════════════════════════════════════════════

## Merged groups for forest plot
disc_s4$Subtype_merged <- dplyr::case_when(
  as.character(disc_s4$Subtype) %in% c("Lum A", "Lum B") ~ "Luminal",
  as.character(disc_s4$Subtype) %in% c("Her2 exp", "Lum Her2") ~ "HER2+",
  as.character(disc_s4$Subtype) == "Triple neg" ~ "TNBC",
  TRUE ~ NA_character_
)
disc_s4$Subtype_merged <- factor(disc_s4$Subtype_merged,
                                 levels = c("Luminal", "HER2+", "TNBC"))

## Five-class subtype for boxplots
disc_s4$Subtype5 <- factor(
  as.character(disc_s4$Subtype),
  levels = c("Lum A", "Lum B", "Lum Her2", "Her2 exp", "Triple neg")
)

cat("  Subtype distribution (merged):\n")
print(table(disc_s4$Subtype_merged, useNA = "ifany"))
cat("\n  Subtype distribution (5-class):\n")
print(table(disc_s4$Subtype5, useNA = "ifany"))


## ═══════════════════════════════════════════════════════════════════════════
## 4. HELPER FUNCTIONS: TABLE OUTPUT (docx + PDF)
## ═══════════════════════════════════════════════════════════════════════════

## ── 4a. Style a flextable consistently (Arial, thin borders, header shading) ──
style_flextable <- function(ft) {
  bdr <- officer::fp_border(color = "black", width = 0.5)
  ft %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 10, part = "body") %>%
    flextable::fontsize(size = 10, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "#D5E8F0", part = "header") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = bdr, part = "header") %>%
    flextable::hline_bottom(border = bdr, part = "header") %>%
    flextable::hline_bottom(border = bdr, part = "body") %>%
    flextable::set_table_properties(width = 1, layout = "autofit") %>%
    flextable::padding(padding.top = 2, padding.bottom = 2, part = "all")
}

## ── 4b. Save forest data as docx table (primary) + PDF fallback ───────────
save_forest_table <- function(forest_df, out_base, title = NULL) {

  ## Prepare data frame for table display
  tab_df <- forest_df %>%
    select(Subgroup, n, events, HR, lower, upper, p, HR_adj, lower_adj, upper_adj, p_adj) %>%
    mutate(
      `UV HR (95% CI)` = ifelse(
        is.na(HR), "-",
        sprintf("%.2f (%.2f-%.2f)", HR, lower, upper)
      ),
      `UV p` = ifelse(
        is.na(p), "-",
        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
      ),
      `MV HR (95% CI)` = ifelse(
        is.na(HR_adj), "-",
        sprintf("%.2f (%.2f-%.2f)", HR_adj, lower_adj, upper_adj)
      ),
      `MV p` = ifelse(
        is.na(p_adj), "-",
        ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))
      )
    ) %>%
    select(Subgroup, n, events, `UV HR (95% CI)`, `UV p`, `MV HR (95% CI)`, `MV p`) %>%
    rename(N = n, Events = events)

  ## Build flextable
  ft <- flextable::flextable(tab_df) %>%
    style_flextable() %>%
    flextable::align(j = 2:7, align = "center", part = "body")

  ## Save as .docx (primary)
  docx_path <- paste0(out_base, ".docx")
  tryCatch({
    doc <- officer::read_docx()
    if (!is.null(title)) {
      doc <- doc %>%
        officer::body_add_par(title, style = "Normal") %>%
        officer::body_add_par("", style = "Normal")
    }
    doc <- doc %>% flextable::body_add_flextable(ft)
    print(doc, target = docx_path)
    cat("  Saved:", basename(docx_path), "\n")
  }, error = function(e) {
    cat("  WARNING: docx export failed:", conditionMessage(e), "\n")
  })

  ## Also save as PDF (fallback/backup for manuscript assembly)
  pdf_path <- paste0(out_base, ".pdf")
  n_rows <- nrow(tab_df)
  core_bg <- rep(c("white", "#f5f5f5"), length.out = n_rows)
  g <- gridExtra::tableGrob(
    tab_df, rows = NULL,
    theme = gridExtra::ttheme_default(
      base_size = 9,
      base_family = "Arial",
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = core_bg)
      ),
      colhead = list(
        fg_params = list(col = "black", fontface = "bold"),
        bg_params = list(fill = "#D5E8F0")
      )
    )
  )

  cairo_pdf(pdf_path, width = 10, height = max(2.5, 1.0 + 0.4 * n_rows))
  grid::grid.newpage()
  if (!is.null(title)) {
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(2, 1,
        heights = grid::unit(c(0.08, 0.92), "npc"))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(title,
      gp = grid::gpar(fontsize = 12, fontface = "bold", fontfamily = "Arial"))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    grid::grid.draw(g)
    grid::popViewport()
    grid::popViewport()
  } else {
    grid::grid.draw(g)
  }
  dev.off()
  cat("  Saved:", basename(pdf_path), "\n")

  invisible(ft)
}

## ── 4c. Save a general data frame as a styled docx table ──────────────────
## Used for supplementary tables and any other tabular output
save_table_docx <- function(tab_df, out_path, title = NULL,
                            col_widths = NULL, align_cols = NULL) {
  ft <- flextable::flextable(tab_df) %>%
    style_flextable()

  if (!is.null(align_cols)) {
    ft <- ft %>% flextable::align(j = align_cols, align = "center", part = "body")
  }
  if (!is.null(col_widths)) {
    ft <- ft %>% flextable::width(width = col_widths)
  }

  doc <- officer::read_docx()
  if (!is.null(title)) {
    doc <- doc %>%
      officer::body_add_par(title, style = "Normal") %>%
      officer::body_add_par("", style = "Normal")
  }
  doc <- doc %>% flextable::body_add_flextable(ft)
  print(doc, target = out_path)
  cat("  Saved:", basename(out_path), "\n")
  invisible(ft)
}


## ═══════════════════════════════════════════════════════════════════════════
## 5. FOREST PLOT — MSPI HR PER SD BY MOLECULAR SUBTYPE
## ═══════════════════════════════════════════════════════════════════════════

cat("\n=== PART 1: Molecular subtype MSPI forest plot ===\n")

## Function: compute MSPI Cox HR per SD for a subset
compute_mspi_hr <- function(d, label) {
  n  <- nrow(d)
  ev <- sum(d$status_num)
  mspi_sd <- sd(d$MSPI, na.rm = TRUE)

  if (n < 15 || ev < 5 || is.na(mspi_sd) || mspi_sd < 0.001) {
    return(data.frame(
      Subgroup = label, n = n, events = ev,
      HR = NA, lower = NA, upper = NA, p = NA,
      HR_adj = NA, lower_adj = NA, upper_adj = NA, p_adj = NA,
      stringsAsFactors = FALSE
    ))
  }

  d$MSPI_sd <- d$MSPI / mspi_sd

  ## UV
  cx_uv <- coxph(Surv(DDFStime2021, status_num) ~ MSPI_sd, data = d)
  s_uv  <- summary(cx_uv)

  ## MV (density-adjusted): raw M1 + M2 + CD3 densities
  has_m1  <- sum(!is.na(d$M1_density))  > 10
  has_m2  <- sum(!is.na(d$M2_density))  > 10
  has_cd3 <- sum(!is.na(d$CD3_density)) > 10

  if (has_cd3 && has_m1 && has_m2) {
    d_mv <- d %>% filter(!is.na(CD3_density), !is.na(M1_density), !is.na(M2_density))
    cx_mv <- coxph(Surv(DDFStime2021, status_num) ~ MSPI_sd + CD3_density + M1_density + M2_density, data = d_mv)
  } else if (has_cd3) {
    d_mv <- d %>% filter(!is.na(CD3_density))
    cx_mv <- coxph(Surv(DDFStime2021, status_num) ~ MSPI_sd + CD3_density, data = d_mv)
    cat("    WARNING:", label, "- using CD3-only MV (M1/M2 density missing)\n")
  } else {
    cx_mv <- cx_uv
    cat("    WARNING:", label, "- no density columns; MV = UV fallback\n")
  }
  s_mv <- summary(cx_mv)
  idx_mv <- which(rownames(s_mv$coefficients) == "MSPI_sd")

  data.frame(
    Subgroup  = label,
    n         = n,
    events    = ev,
    HR        = s_uv$conf.int[1, 1],
    lower     = s_uv$conf.int[1, 3],
    upper     = s_uv$conf.int[1, 4],
    p         = s_uv$coefficients[1, 5],
    HR_adj    = s_mv$conf.int[idx_mv, 1],
    lower_adj = s_mv$conf.int[idx_mv, 3],
    upper_adj = s_mv$conf.int[idx_mv, 4],
    p_adj     = s_mv$coefficients[idx_mv, 5],
    stringsAsFactors = FALSE
  )
}

## Compute for each subgroup
forest_data <- rbind(
  compute_mspi_hr(disc_s4, "All-comers"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype_merged == "Luminal"), "Luminal (A+B)"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype5 == "Lum A"), "  Lum A"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype5 == "Lum B"), "  Lum B"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype_merged == "HER2+"), "HER2+ (all)"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype5 == "Her2 exp"), "  HER2-enriched"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype5 == "Lum Her2"), "  Luminal-HER2"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype_merged == "TNBC"), "TNBC")
)

## Print summary
cat("\nMSPI HR per SD by molecular subtype (DDFS):\n")
cat(sprintf("%-20s  n=%4s  ev=%3s  UV: HR=%5s (%5s-%5s) p=%8s  |  MV: HR=%5s p=%8s\n",
            "Subgroup", "n", "ev", "HR", "lo", "hi", "p", "HR_adj", "p_adj"))
for (i in seq_len(nrow(forest_data))) {
  r <- forest_data[i, ]
  cat(sprintf("%-20s  n=%4d  ev=%3d  UV: HR=%5.2f (%5.2f-%5.2f) p=%8.2e  |  MV: HR=%5.2f p=%8.2e\n",
              r$Subgroup, r$n, r$events,
              ifelse(is.na(r$HR), NA, r$HR),
              ifelse(is.na(r$lower), NA, r$lower),
              ifelse(is.na(r$upper), NA, r$upper),
              ifelse(is.na(r$p), NA, r$p),
              ifelse(is.na(r$HR_adj), NA, r$HR_adj),
              ifelse(is.na(r$p_adj), NA, r$p_adj)))
}

## Save CSV
write.csv(forest_data, file.path(output_dir, "discovery_MSPI_subtype_forest_data.csv"),
          row.names = FALSE)

## Save forest table PDF (CHANGE 2)
save_forest_table(forest_data,
                      file.path(output_dir, "discovery_MSPI_subtype_forest_data_TABLE"),
                      "Discovery cohort MSPI by Molecular Subtype - Forest Data")

## ── Forest plot (UV) ────────────────────────────────────────────────────

forest_data$Subgroup <- factor(forest_data$Subgroup,
                                levels = rev(forest_data$Subgroup))

## Color by significance
forest_data$sig <- ifelse(is.na(forest_data$p), "NS",
                          ifelse(forest_data$p < 0.05, "Significant", "NS"))

## Format labels
forest_data$label_text <- sprintf(
  "%.2f (%.2f-%.2f)%s",
  forest_data$HR, forest_data$lower, forest_data$upper,
  ifelse(forest_data$p < 0.001, "***",
         ifelse(forest_data$p < 0.01, "**",
                ifelse(forest_data$p < 0.05, "*", "")))
)
forest_data$label_text[is.na(forest_data$HR)] <- "—"

## ── v4-matched forest theme ──────────────────────────────────────────────
## Matching style from Script_Integrated_MSPI_Analysis_v4 T4/T6 forests:
## base_size=14, shape=18, size=6, linewidth=0.9, grey65 vline, no y-grid

forest_theme <- theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 15, face = "bold"),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    axis.text.y        = element_text(size = 12)
  )

## Helper: compute forest height
fh <- function(n_rows) max(3.5, 1.5 + 0.7 * n_rows)

## ── Forest plot (UV) — full subtypes ─────────────────────────────────────

p_forest_uv <- ggplot(forest_data %>% filter(!is.na(HR)),
                       aes(x = HR, y = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.22, linewidth = 0.9, colour = "#4d4d4d") +
  geom_point(aes(colour = sig), size = 6, shape = 18) +
  geom_text(aes(label = label_text),
            hjust = -0.08, vjust = -0.6, size = 3.5, family = "Arial") +
  scale_colour_manual(values = c("Significant" = "#c0392b", "NS" = "grey60"),
                       guide = "none") +
  scale_x_continuous(breaks = seq(0.5, 2.5, 0.25), limits = c(0.5, 2.5)) +
  labs(title = "MSPI and DDFS by molecular subtype",
       subtitle = "Discovery cohort discovery cohort - univariable HR per SD",
       x = "Hazard Ratio per SD (95% CI)", y = NULL) +
  forest_theme

cairo_pdf(file.path(output_dir, "discovery_MSPI_subtype_forest_UV.pdf"),
          width = 6.5, height = fh(nrow(forest_data %>% filter(!is.na(HR)))))
print(p_forest_uv)
dev.off()
cat("\n  Saved: discovery_MSPI_subtype_forest_UV.pdf\n")

## ── Forest plot (MV density-adjusted) — full subtypes ────────────────────

forest_data$label_text_adj <- sprintf(
  "%.2f (%.2f-%.2f)%s",
  forest_data$HR_adj, forest_data$lower_adj, forest_data$upper_adj,
  ifelse(forest_data$p_adj < 0.001, "***",
         ifelse(forest_data$p_adj < 0.01, "**",
                ifelse(forest_data$p_adj < 0.05, "*", "")))
)
forest_data$label_text_adj[is.na(forest_data$HR_adj)] <- "—"

forest_data$sig_adj <- ifelse(is.na(forest_data$p_adj), "NS",
                               ifelse(forest_data$p_adj < 0.05, "Significant", "NS"))

p_forest_mv <- ggplot(forest_data %>% filter(!is.na(HR_adj)),
                       aes(x = HR_adj, y = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
  geom_errorbarh(aes(xmin = lower_adj, xmax = upper_adj),
                 height = 0.22, linewidth = 0.9, colour = "#4d4d4d") +
  geom_point(aes(colour = sig_adj), size = 6, shape = 18) +
  geom_text(aes(label = label_text_adj),
            hjust = -0.08, vjust = -0.6, size = 3.5, family = "Arial") +
  scale_colour_manual(values = c("Significant" = "#c0392b", "NS" = "grey60"),
                       guide = "none") +
  scale_x_continuous(breaks = seq(0.5, 2.5, 0.25), limits = c(0.5, 2.5)) +
  labs(title = "MSPI and DDFS by molecular subtype",
       subtitle = "Discovery cohort discovery cohort - density-adjusted HR per SD",
       x = "Hazard Ratio per SD (95% CI)", y = NULL) +
  forest_theme

cairo_pdf(file.path(output_dir, "discovery_MSPI_subtype_forest_MV.pdf"),
          width = 6.5, height = fh(nrow(forest_data %>% filter(!is.na(HR_adj)))))
print(p_forest_mv)
dev.off()
cat("  Saved: discovery_MSPI_subtype_forest_MV.pdf\n")

## ── Combined UV + MV — full subtypes ─────────────────────────────────────
## UV on top (first in factor levels), Density-adjusted below

forest_long <- rbind(
  forest_data %>%
    select(Subgroup, HR, lower, upper, p, sig) %>%
    mutate(Model = "Univariable"),
  forest_data %>%
    select(Subgroup, HR_adj, lower_adj, upper_adj, p_adj, sig_adj) %>%
    rename(HR = HR_adj, lower = lower_adj, upper = upper_adj,
           p = p_adj, sig = sig_adj) %>%
    mutate(Model = "Density-adjusted")
) %>%
  filter(!is.na(HR))

forest_long$Model <- factor(forest_long$Model,
                             levels = c("Density-adjusted", "Univariable"))

p_forest_combined <- ggplot(forest_long,
                             aes(x = HR, y = Subgroup, colour = sig, shape = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.2, linewidth = 0.8, colour = "#4d4d4d",
                 position = position_dodge(width = 0.7)) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  scale_colour_manual(values = c("Significant" = "#c0392b", "NS" = "grey55"),
                       guide = "none") +
  scale_shape_manual(values = c("Univariable" = 16, "Density-adjusted" = 18)) +
  scale_x_continuous(breaks = seq(0.5, 2.5, 0.5), limits = c(0.4, 2.6)) +
  labs(title = "MSPI and DDFS by molecular subtype",
       subtitle = "Discovery cohort discovery cohort",
       x = "Hazard Ratio per SD (95% CI)", y = NULL, shape = NULL) +
  forest_theme +
  theme(legend.position = "bottom")

cairo_pdf(file.path(output_dir, "discovery_MSPI_subtype_forest_combined.pdf"),
          width = 6.5, height = fh(length(unique(forest_long$Subgroup))))
print(p_forest_combined)
dev.off()
cat("  Saved: discovery_MSPI_subtype_forest_combined.pdf\n")

## ── MERGED-ONLY forest plots (no Lum A/B or HER2 split) ─────────────────

forest_merged <- rbind(
  compute_mspi_hr(disc_s4, "All-comers"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype_merged == "Luminal"), "Luminal"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype_merged == "HER2+"), "HER2+"),
  compute_mspi_hr(disc_s4 %>% filter(Subtype_merged == "TNBC"), "TNBC")
)

forest_merged$Subgroup <- factor(forest_merged$Subgroup,
                                  levels = rev(forest_merged$Subgroup))
forest_merged$sig <- ifelse(is.na(forest_merged$p), "NS",
                            ifelse(forest_merged$p < 0.05, "Significant", "NS"))
forest_merged$sig_adj <- ifelse(is.na(forest_merged$p_adj), "NS",
                                ifelse(forest_merged$p_adj < 0.05, "Significant", "NS"))
forest_merged$label_text <- ifelse(
  is.na(forest_merged$HR), "—",
  sprintf("%.2f (%.2f-%.2f)%s", forest_merged$HR, forest_merged$lower, forest_merged$upper,
          ifelse(forest_merged$p < 0.001, "***", ifelse(forest_merged$p < 0.01, "**",
          ifelse(forest_merged$p < 0.05, "*", "")))))
forest_merged$label_text_adj <- ifelse(
  is.na(forest_merged$HR_adj), "—",
  sprintf("%.2f (%.2f-%.2f)%s", forest_merged$HR_adj, forest_merged$lower_adj, forest_merged$upper_adj,
          ifelse(forest_merged$p_adj < 0.001, "***", ifelse(forest_merged$p_adj < 0.01, "**",
          ifelse(forest_merged$p_adj < 0.05, "*", "")))))

## UV merged-only
p_forest_uv_m <- ggplot(forest_merged %>% filter(!is.na(HR)),
                          aes(x = HR, y = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.22, linewidth = 0.9, colour = "#4d4d4d") +
  geom_point(aes(colour = sig), size = 6, shape = 18) +
  geom_text(aes(label = label_text),
            hjust = -0.08, vjust = -0.6, size = 3.5, family = "Arial") +
  scale_colour_manual(values = c("Significant" = "#c0392b", "NS" = "grey60"), guide = "none") +
  scale_x_continuous(breaks = seq(0.5, 2.5, 0.25), limits = c(0.5, 2.5)) +
  labs(title = "MSPI and DDFS by molecular subtype (merged)",
       subtitle = "Discovery cohort discovery cohort - univariable HR per SD",
       x = "Hazard Ratio per SD (95% CI)", y = NULL) +
  forest_theme

cairo_pdf(file.path(output_dir, "discovery_MSPI_subtype_forest_UV_merged.pdf"),
          width = 6.5, height = fh(nrow(forest_merged %>% filter(!is.na(HR)))))
print(p_forest_uv_m)
dev.off()
cat("  Saved: discovery_MSPI_subtype_forest_UV_merged.pdf\n")

## MV merged-only
p_forest_mv_m <- ggplot(forest_merged %>% filter(!is.na(HR_adj)),
                          aes(x = HR_adj, y = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
  geom_errorbarh(aes(xmin = lower_adj, xmax = upper_adj),
                 height = 0.22, linewidth = 0.9, colour = "#4d4d4d") +
  geom_point(aes(colour = sig_adj), size = 6, shape = 18) +
  geom_text(aes(label = label_text_adj),
            hjust = -0.08, vjust = -0.6, size = 3.5, family = "Arial") +
  scale_colour_manual(values = c("Significant" = "#c0392b", "NS" = "grey60"), guide = "none") +
  scale_x_continuous(breaks = seq(0.5, 2.5, 0.25), limits = c(0.5, 2.5)) +
  labs(title = "MSPI and DDFS by molecular subtype (merged)",
       subtitle = "Discovery cohort discovery cohort - density-adjusted HR per SD",
       x = "Hazard Ratio per SD (95% CI)", y = NULL) +
  forest_theme

cairo_pdf(file.path(output_dir, "discovery_MSPI_subtype_forest_MV_merged.pdf"),
          width = 6.5, height = fh(nrow(forest_merged %>% filter(!is.na(HR_adj)))))
print(p_forest_mv_m)
dev.off()
cat("  Saved: discovery_MSPI_subtype_forest_MV_merged.pdf\n")

## Combined UV + MV merged-only — T5-style (one row per model, no dodging)
forest_merged_long <- rbind(
  forest_merged %>%
    select(Subgroup, HR, lower, upper, p, sig) %>%
    mutate(Model = "Univariable"),
  forest_merged %>%
    select(Subgroup, HR_adj, lower_adj, upper_adj, p_adj, sig_adj) %>%
    rename(HR = HR_adj, lower = lower_adj, upper = upper_adj, p = p_adj, sig = sig_adj) %>%
    mutate(Model = "Density-adjusted")
) %>%
  filter(!is.na(HR))

## Build a combined label: "Subgroup — Model" so each row is unique
forest_merged_long$Label <- paste0(forest_merged_long$Subgroup, " - ", forest_merged_long$Model)

## Numeric y positions: tight within subgroup, bigger gap between subgroups
## All-comers at top → highest y value; TNBC at bottom → lowest
sg_order <- c("All-comers", "Luminal", "HER2+", "TNBC")
within_gap  <- 0.6
between_gap <- 1.0

ypos   <- c()
ylabs  <- c()
ymap   <- list()
y_cur  <- length(sg_order) * (within_gap + between_gap)
for (i in seq_along(sg_order)) {
  sg <- sg_order[i]
  if (i > 1) y_cur <- y_cur - between_gap
  ## UV on top (higher y), MV just below
  ymap[[paste0(sg, " - Univariable")]]      <- y_cur
  ymap[[paste0(sg, " - Density-adjusted")]] <- y_cur - within_gap
  ypos  <- c(ypos, y_cur, y_cur - within_gap)
  ylabs <- c(ylabs, paste0(sg, " - Univariable"), paste0(sg, " - Density-adjusted"))
  y_cur <- y_cur - within_gap
}

forest_merged_long$ynum <- sapply(as.character(forest_merged_long$Label),
                                   function(lb) ymap[[lb]])

## Colour: significant → subgroup colour, NS → grey
sg_cols <- c(
  "All-comers" = "#a67c52",
  "Luminal"    = "#d4b86a",
  "HER2+"      = "#5b7fa0",
  "TNBC"       = "#c47a7a"
)

forest_merged_long$point_col <- ifelse(
  forest_merged_long$sig == "Significant",
  sg_cols[as.character(forest_merged_long$Subgroup)],
  "grey55"
)

p_forest_comb_m <- ggplot(forest_merged_long,
                          aes(x = HR, y = ynum, xmin = lower, xmax = upper,
                              shape = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(height = 0.2, linewidth = 0.9, colour = "#4d4d4d") +
  geom_point(size = 5, colour = forest_merged_long$point_col) +
  geom_text(aes(label = sprintf("HR=%.2f (p=%s)",
                                HR, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))),
            hjust = -0.08, size = 3.5, family = "Arial", colour = forest_merged_long$point_col,
            show.legend = FALSE) +
  scale_y_continuous(breaks = ypos, labels = ylabs, expand = expansion(mult = 0.08)) +
  scale_shape_manual(values = c("Univariable" = 16, "Density-adjusted" = 18),
                     guide = guide_legend(title = NULL)) +
  labs(title = "MSPI and DDFS by molecular subtype (merged)",
       subtitle = "Discovery cohort discovery cohort",
       x = "Hazard Ratio per SD (95% CI)", y = NULL) +
  theme_bw(base_size = 14, base_family = "Arial") +
  xlim(min(forest_merged_long$lower, na.rm = TRUE) - 0.05,
       max(forest_merged_long$upper, na.rm = TRUE) + 0.6) +
  theme(legend.position = "bottom")

cairo_pdf(file.path(output_dir, "discovery_MSPI_subtype_forest_combined_merged.pdf"),
          width = 6.5, height = 5.0)
print(p_forest_comb_m)
dev.off()
cat("  Saved: discovery_MSPI_subtype_forest_combined_merged.pdf\n")

write.csv(forest_merged, file.path(output_dir, "discovery_MSPI_subtype_forest_data_merged.csv"),
          row.names = FALSE)

## Save forest merged table PDF (CHANGE 2)
save_forest_table(forest_merged,
                      file.path(output_dir, "discovery_MSPI_subtype_forest_data_merged_TABLE"),
                      "Discovery cohort MSPI by Molecular Subtype (Merged) - Forest Data")


## ═══════════════════════════════════════════════════════════════════════════
## 5b. M1/M2 INTERACTION FOREST — MSPI JUSTIFICATION
##     Adapted from Alfonso's "Script for interactions Figure 2 and 4.R"
##     Restyled to match discovery_MSPI_subtype_forest_combined_merged.pdf design
## ═══════════════════════════════════════════════════════════════════════════

cat("\n=== PART 1b: M1/M2 interaction forest (MSPI justification) ===\n")

## ── 5b-1. Compute M1 and M2 additive Cox HRs per SD for each subgroup ───

compute_m1m2_additive_hr <- function(d, label, time_var, status_var) {
  n  <- nrow(d)
  ev <- sum(d[[status_var]], na.rm = TRUE)
  if (n < 15 || ev < 5) {
    return(data.frame(Subgroup = label, n = n, events = ev,
      M1_HR = NA, M1_lo = NA, M1_hi = NA, M1_p = NA,
      M2_HR = NA, M2_lo = NA, M2_hi = NA, M2_p = NA,
      LRT_p = NA, EPP = NA, stringsAsFactors = FALSE))
  }

  m1_sd <- sd(d$M1_score, na.rm = TRUE)
  m2_sd <- sd(d$M2_score, na.rm = TRUE)
  if (is.na(m1_sd) || m1_sd < 0.001 || is.na(m2_sd) || m2_sd < 0.001) {
    return(data.frame(Subgroup = label, n = n, events = ev,
      M1_HR = NA, M1_lo = NA, M1_hi = NA, M1_p = NA,
      M2_HR = NA, M2_lo = NA, M2_hi = NA, M2_p = NA,
      LRT_p = NA, EPP = NA, stringsAsFactors = FALSE))
  }

  d$M1_sd <- d$M1_score / m1_sd
  d$M2_sd <- d$M2_score / m2_sd

  surv_obj <- Surv(d[[time_var]], d[[status_var]])

  ## Additive model: M1 + M2 as independent terms
  cx_add <- coxph(surv_obj ~ M1_sd + M2_sd, data = d)
  s_add  <- summary(cx_add)
  i_m1 <- which(rownames(s_add$coefficients) == "M1_sd")
  i_m2 <- which(rownames(s_add$coefficients) == "M2_sd")

  ## Interaction model: M1 + M2 + M1:M2
  cx_int <- coxph(surv_obj ~ M1_sd * M2_sd, data = d)
  lrt <- anova(cx_add, cx_int, test = "LRT")
  lrt_p <- lrt[["Pr(>|Chi|)"]][2]

  epp <- ev / 3  # 3 parameters in additive model

  data.frame(
    Subgroup = label, n = n, events = ev,
    M1_HR = s_add$conf.int[i_m1, 1],
    M1_lo = s_add$conf.int[i_m1, 3],
    M1_hi = s_add$conf.int[i_m1, 4],
    M1_p  = s_add$coefficients[i_m1, 5],
    M2_HR = s_add$conf.int[i_m2, 1],
    M2_lo = s_add$conf.int[i_m2, 3],
    M2_hi = s_add$conf.int[i_m2, 4],
    M2_p  = s_add$coefficients[i_m2, 5],
    LRT_p = lrt_p,
    EPP   = epp,
    stringsAsFactors = FALSE
  )
}

## ── 5b-2. Compute for Discovery cohort subgroups (same as merged forest) ──────────────

int_bcy <- rbind(
  compute_m1m2_additive_hr(disc_s4, "All-comers",
    "DDFStime2021", "status_num"),
  compute_m1m2_additive_hr(disc_s4 %>% filter(Subtype_merged == "Luminal"), "Luminal",
    "DDFStime2021", "status_num"),
  compute_m1m2_additive_hr(disc_s4 %>% filter(Subtype_merged == "HER2+"), "HER2+",
    "DDFStime2021", "status_num"),
  compute_m1m2_additive_hr(disc_s4 %>% filter(Subtype_merged == "TNBC"), "TNBC",
    "DDFStime2021", "status_num")
)

cat("\nM1/M2 additive HRs per SD (Discovery cohort, DDFS):\n")
print(int_bcy)

## Save CSV + docx
write.csv(int_bcy, file.path(output_dir, "discovery_M1M2_interaction_data.csv"),
          row.names = FALSE)
if (requireNamespace("flextable", quietly = TRUE)) {
  ft_int <- flextable::flextable(int_bcy)
  ft_int <- ft_int %>%
    flextable::set_header_labels(
      Subgroup = "Subgroup", n = "n", events = "Events",
      M1_HR = "M1 HR", M1_lo = "M1 95% lo", M1_hi = "M1 95% hi", M1_p = "M1 p",
      M2_HR = "M2 HR", M2_lo = "M2 95% lo", M2_hi = "M2 95% hi", M2_p = "M2 p",
      LRT_p = "LRT p", EPP = "EPP") %>%
    flextable::autofit() %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::bold(part = "header")
  doc_int <- officer::read_docx() %>%
    officer::body_add_par("Discovery cohort M1/M2 Additive & Interaction Test", style = "Normal") %>%
    flextable::body_add_flextable(ft_int)
  print(doc_int, target = file.path(output_dir, "discovery_M1M2_interaction_data.docx"))
}

## ── 5b-3. Reshape to long format for forest plot ─────────────────────────

sg_order <- c("All-comers", "Luminal", "HER2+", "TNBC")
int_bcy$Subgroup <- factor(int_bcy$Subgroup, levels = sg_order)

int_long <- rbind(
  int_bcy %>%
    transmute(Subgroup, Metric = "M1 affinity",
              HR = M1_HR, lower = M1_lo, upper = M1_hi, p = M1_p,
              LRT_p, EPP),
  int_bcy %>%
    transmute(Subgroup, Metric = "M2 affinity",
              HR = M2_HR, lower = M2_lo, upper = M2_hi, p = M2_p,
              LRT_p, EPP)
) %>%
  filter(!is.na(HR))

## Subgroup colors (matching combined_merged forest)
sg_cols <- c(
  "All-comers" = "#a67c52",
  "Luminal"    = "#d4b86a",
  "HER2+"      = "#5b7fa0",
  "TNBC"       = "#c47a7a"
)

## Numeric y positions: M1 above M2 within each subgroup, gap between groups
within_gap  <- 0.6
between_gap <- 1.0

ypos_map <- list()
ylabs_vec <- c()
ypos_vec  <- c()
y_cur <- length(sg_order) * (within_gap + between_gap)
for (i in seq_along(sg_order)) {
  sg <- sg_order[i]
  if (i > 1) y_cur <- y_cur - between_gap
  ypos_map[[paste0(sg, "_M1 affinity")]] <- y_cur
  ypos_map[[paste0(sg, "_M2 affinity")]] <- y_cur - within_gap
  ypos_vec <- c(ypos_vec, y_cur, y_cur - within_gap)
  ylabs_vec <- c(ylabs_vec,
                 paste0(sg, " - M1 affinity"),
                 paste0(sg, " - M2 affinity"))
  y_cur <- y_cur - within_gap
}

int_long$ynum <- sapply(seq_len(nrow(int_long)), function(j) {
  ypos_map[[paste0(as.character(int_long$Subgroup[j]), "_", int_long$Metric[j])]]
})

int_long$sig <- ifelse(int_long$p < 0.05, "Significant", "NS")

## Point colour: significant → metric colour (M1 blue, M2 red), NS → grey
col_M1 <- "#1a7abf"
col_M2 <- "#c0392b"
int_long$point_col <- ifelse(int_long$sig == "NS", "grey55",
  ifelse(int_long$Metric == "M1 affinity", col_M1, col_M2))

int_long$text_col <- int_long$point_col

## HR label
int_long$hr_label <- sprintf("HR=%.2f (p=%s)",
  int_long$HR,
  ifelse(int_long$p < 0.001, "<0.001", sprintf("%.3f", int_long$p)))

## ── 5b-4. FOREST PANEL ──────────────────────────────────────────────────

p_int_forest <- ggplot(int_long,
    aes(x = HR, y = ynum, xmin = lower, xmax = upper,
        shape = Metric)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(height = 0.2, linewidth = 0.9, colour = "#4d4d4d") +
  geom_point(size = 5, colour = int_long$point_col) +
  geom_text(aes(label = hr_label),
            hjust = -0.08, size = 3.5, family = "Arial",
            colour = int_long$text_col, show.legend = FALSE) +
  scale_y_continuous(breaks = ypos_vec, labels = ylabs_vec,
                     expand = expansion(mult = 0.08)) +
  scale_shape_manual(values = c("M1 affinity" = 16, "M2 affinity" = 18),
                     guide = guide_legend(title = NULL)) +
  scale_x_continuous(
    breaks = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
    limits = c(min(int_long$lower, na.rm = TRUE) - 0.05,
               max(int_long$upper, na.rm = TRUE) + 0.8)) +
  labs(title = "M1 and M2 affinity - additive model",
       subtitle = "Discovery cohort discovery cohort (DDFS)",
       x = "Hazard Ratio per SD (95% CI)", y = NULL) +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 15, face = "bold"),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    legend.position    = "bottom"
  )

## ── 5b-5. LRT PANEL ─────────────────────────────────────────────────────

int_lrt <- int_bcy %>%
  filter(!is.na(LRT_p)) %>%
  mutate(
    y_num     = sapply(as.character(Subgroup), function(sg)
      mean(c(ypos_map[[paste0(sg, "_M1 affinity")]],
             ypos_map[[paste0(sg, "_M2 affinity")]]))),
    lrt_sig   = LRT_p < 0.05,
    lrt_label = ifelse(LRT_p < 0.001, "p < 0.001",
                       sprintf("p = %.2f", LRT_p)),
    lrt_col   = ifelse(lrt_sig, "#e6550d", "#2ca25f")
  )

p_int_lrt <- ggplot(int_lrt, aes(x = LRT_p, y = y_num)) +
  geom_vline(xintercept = 0.05, colour = "#e6550d",
             linewidth = 0.7, linetype = "dashed") +
  annotate("text", x = 0.07,
           y = max(int_lrt$y_num, na.rm = TRUE) + 0.8,
           label = "p = 0.05", colour = "#e6550d",
           size = 3.0, hjust = 0, fontface = "italic") +
  geom_segment(aes(x = 0, xend = LRT_p, y = y_num, yend = y_num),
               colour = "grey55", linewidth = 0.65) +
  geom_point(aes(colour = lrt_sig, size = EPP), shape = 16) +
  geom_text(aes(label = lrt_label),
            hjust = -0.25, size = 3.5, colour = "grey20", family = "Arial") +
  scale_colour_manual(
    values = c("FALSE" = "#2ca25f", "TRUE" = "#e6550d"),
    labels = c("FALSE" = "Additive (p > 0.05)",
               "TRUE"  = "Significant interaction"),
    guide = guide_legend(title = NULL, override.aes = list(size = 4))
  ) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  scale_x_continuous(limits = c(0, 1.3),
    breaks = c(0.25, 0.50, 0.75, 1.00)) +
  scale_y_continuous(
    limits = range(ypos_vec) + c(-0.5, 0.5),
    expand = c(0, 0)) +
  labs(x = "LRT p-value (M1 x M2 interaction)",
       y = NULL,
       title = "M1 x M2 interaction test",
       subtitle = "All p > 0.05: additive contributions") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(size = 15, face = "bold"),
    plot.subtitle      = element_text(size = 10, colour = "grey50"),
    legend.position    = "bottom"
  )

## ── 5b-6. COMBINED FIGURE (patchwork) ────────────────────────────────────

p_int_combined <- p_int_forest + p_int_lrt +
  patchwork::plot_layout(widths = c(3.0, 0.8)) +
  patchwork::plot_annotation(
    title = "M1 and M2 spatial affinity act independently - justification for MSPI composite",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13,
                                colour = "grey10", family = "Arial",
                                margin = margin(b = 6)),
      plot.margin = margin(6, 6, 4, 6)
    )
  )

cairo_pdf(file.path(output_dir, "discovery_M1M2_interaction_forest.pdf"),
          width = 12, height = 5.0)
print(p_int_combined)
dev.off()
cat("  Saved: discovery_M1M2_interaction_forest.pdf\n")

## PNG version
ggsave(file.path(output_dir, "discovery_M1M2_interaction_forest.png"),
       p_int_combined, width = 12, height = 5.0, dpi = 300, bg = "white")
cat("  Saved: discovery_M1M2_interaction_forest.png\n")


## ═══════════════════════════════════════════════════════════════════════════
## 6. BOXPLOTS — MSPI BY MOLECULAR SUBTYPE
## ═══════════════════════════════════════════════════════════════════════════

cat("\n=== PART 2: MSPI boxplots by molecular subtype ===\n")

## ── 6a. Five-class subtype boxplot ──────────────────────────────────────

sub5_colors <- c(
  "Lum A"      = "#4e79a7",
  "Lum B"      = "#f28e2b",
  "Lum Her2"   = "#e15759",
  "Her2 exp"   = "#76b7b2",
  "Triple neg" = "#59a14f"
)

## Kruskal-Wallis + pairwise Dunn tests
kw_test <- disc_s4 %>%
  filter(!is.na(Subtype5)) %>%
  kruskal_test(MSPI ~ Subtype5)

dunn_test_res <- disc_s4 %>%
  filter(!is.na(Subtype5)) %>%
  dunn_test(MSPI ~ Subtype5, p.adjust.method = "BH")

cat("  Kruskal-Wallis p =", kw_test$p, "\n")
cat("  Pairwise Dunn tests (BH-adjusted, significant pairs):\n")
print(dunn_test_res %>% filter(p.adj < 0.05) %>% select(group1, group2, p.adj))

p_box5 <- ggplot(disc_s4 %>% filter(!is.na(Subtype5)),
                  aes(x = Subtype5, y = MSPI, fill = Subtype5)) +
  geom_jitter(width = 0.2, alpha = 0.15, size = 0.8, colour = "grey40") +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.80) +
  scale_fill_manual(values = sub5_colors, guide = "none") +
  labs(
    title = "MSPI distribution by molecular subtype",
    subtitle = sprintf("Discovery cohort discovery cohort (n = %d) | Kruskal-Wallis p = %.3g",
                        sum(!is.na(disc_s4$Subtype5)), kw_test$p),
    x = NULL, y = "MSPI (M2 score - M1 score)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 11)
  )

cairo_pdf(file.path(output_dir, "discovery_MSPI_boxplot_5subtypes.pdf"),
          width = 7, height = 5)
print(p_box5)
dev.off()
cat("  Saved: discovery_MSPI_boxplot_5subtypes.pdf\n")

## ── 6b. Three-class subtype boxplot ─────────────────────────────────────

sub3_colors <- c("Luminal" = "#4e79a7", "HER2+" = "#f28e2b", "TNBC" = "#e15759")

kw_test_3 <- disc_s4 %>%
  filter(!is.na(Subtype_merged)) %>%
  kruskal_test(MSPI ~ Subtype_merged)

p_box3 <- ggplot(disc_s4 %>% filter(!is.na(Subtype_merged)),
                  aes(x = Subtype_merged, y = MSPI, fill = Subtype_merged)) +
  geom_jitter(width = 0.2, alpha = 0.15, size = 0.8, colour = "grey40") +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.80) +
  scale_fill_manual(values = sub3_colors, guide = "none") +
  labs(
    title = "MSPI by molecular subtype",
    subtitle = sprintf("Discovery cohort (n = %d) | Kruskal-Wallis p = %.3g",
                        sum(!is.na(disc_s4$Subtype_merged)), kw_test_3$p),
    x = NULL, y = "MSPI (M2 score - M1 score)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    axis.text.x   = element_text(size = 12)
  )

cairo_pdf(file.path(output_dir, "discovery_MSPI_boxplot_3subtypes.pdf"),
          width = 6, height = 5)
print(p_box3)
dev.off()
cat("  Saved: discovery_MSPI_boxplot_3subtypes.pdf\n")

## ── 6c. M1 and M2 affinity scores by subtype ────────────────────────────

disc_score_long <- disc_s4 %>%
  filter(!is.na(Subtype5)) %>%
  select(Subtype5, M1_score, M2_score) %>%
  tidyr::pivot_longer(cols = c(M1_score, M2_score),
                      names_to = "Marker", values_to = "Score") %>%
  mutate(Marker = ifelse(Marker == "M1_score", "M1 affinity", "M2 affinity"))

p_m1m2 <- ggplot(disc_score_long, aes(x = Subtype5, y = Score, fill = Subtype5)) +
  geom_jitter(width = 0.2, alpha = 0.10, size = 0.5, colour = "grey40") +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.80) +
  facet_wrap(~ Marker, scales = "free_y") +
  scale_fill_manual(values = sub5_colors, guide = "none") +
  labs(
    title = "M1 and M2 affinity scores by molecular subtype",
    subtitle = "Discovery cohort discovery cohort",
    x = NULL, y = "Affinity score (inverted, higher = closer)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    strip.text    = element_text(size = 12, face = "bold"),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 10)
  )

cairo_pdf(file.path(output_dir, "discovery_MSPI_boxplot_M1_M2_by_subtype.pdf"),
          width = 10, height = 5)
print(p_m1m2)
dev.off()
cat("  Saved: discovery_MSPI_boxplot_M1_M2_by_subtype.pdf\n")

## ── 6d. Cell density boxplots by subtype ────────────────────────────────

## Prepare long data frame for densities
disc_dens_long <- disc_s4 %>%
  filter(!is.na(Subtype_merged), !is.na(CD3_density), !is.na(M1_density), !is.na(M2_density)) %>%
  select(Subtype_merged, CD3_density, M1_density, M2_density) %>%
  tidyr::pivot_longer(cols = c(CD3_density, M1_density, M2_density),
                      names_to = "Cell_type", values_to = "Density") %>%
  mutate(Cell_type = recode(Cell_type,
                             CD3_density = "CD3+ T cells",
                             M1_density = "M1 macrophages",
                             M2_density = "M2 macrophages"))

p_dens <- ggplot(disc_dens_long, aes(x = Subtype_merged, y = Density, fill = Subtype_merged)) +
  geom_jitter(width = 0.2, alpha = 0.10, size = 0.5, colour = "grey40") +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.80) +
  facet_wrap(~ Cell_type, scales = "free_y") +
  scale_fill_manual(values = sub3_colors, guide = "none") +
  scale_y_sqrt(breaks = c(1, 10, 100, 1000)) +
  labs(
    title = "Cell densities by molecular subtype",
    subtitle = "Discovery cohort discovery cohort",
    x = NULL, y = "Cell density (counts per mm2, sqrt scale)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    strip.text    = element_text(size = 12, face = "bold"),
    axis.text.x   = element_text(size = 11)
  )

cairo_pdf(file.path(output_dir, "discovery_Cell_density_boxplot_3subtypes.pdf"),
          width = 10, height = 5)
print(p_dens)
dev.off()
cat("  Saved: discovery_Cell_density_boxplot_3subtypes.pdf\n")

## ── 6e. Clinical associations: MSPI vs Age, Grade, Stage, etc. ──────────

cat("\n  Building Discovery cohort clinical association plots...\n")

## Initialize storage for summary
clin_summary <- data.frame()

## List of clinical variables to plot
clin_vars_plot <- c("Age_cat_2", "Grade", "Stage_TNM_44", "Subtype_merged", "Syst_treat_NE", "diadtm_final")
clin_vars_plot <- clin_vars_plot[clin_vars_plot %in% names(disc_s4)]

clin_plots <- list()
clin_plots_nd <- list()

for (var in clin_vars_plot) {
  disc_clin_filt <- disc_s4 %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(group = factor(.data[[var]]))

  if (nrow(disc_clin_filt) < 20) next

  ## Compute p-value (Kruskal-Wallis)
  wt <- tryCatch(
    kruskal_test(MSPI ~ group, data = disc_clin_filt),
    error = function(e) list(p = NA)
  )

  ## Base plot
  p_base <- ggplot(disc_clin_filt, aes(x = group, y = MSPI)) +
    labs(
      title = sprintf("MSPI by %s", var),
      subtitle = sprintf("KW p = %.3g", wt$p),
      x = NULL, y = "MSPI"
    ) +
    theme_bw(base_size = 11, base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )

  ## With dots
  p_with_dots <- p_base +
    geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7)

  clin_plots[[var]] <- p_with_dots

  ## Without dots
  p_no_dots <- p_base +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7)

  clin_plots_nd[[var]] <- p_no_dots

  ## Add to summary
  clin_summary <- rbind(clin_summary,
    data.frame(Variable = var, p_value = wt$p, stringsAsFactors = FALSE))
}

## Save individual PDFs
for (var in names(clin_plots)) {
  fname <- sprintf("discovery_MSPI_clinical_%s.pdf", var)
  cairo_pdf(file.path(output_dir, fname), width = 6, height = 4.5)
  print(clin_plots[[var]])
  dev.off()
}
cat("  Saved individual clinical boxplots\n")

## Combined panel
if (length(clin_plots) >= 4) {
  ncol_grid <- min(3, length(clin_plots))
  nrow_grid <- ceiling(length(clin_plots) / ncol_grid)

  cairo_pdf(file.path(output_dir, "discovery_MSPI_clinical_associations_panel.pdf"),
            width = 5.5 * ncol_grid, height = 4.5 * nrow_grid)
  do.call(gridExtra::grid.arrange, c(clin_plots, ncol = ncol_grid))
  dev.off()
  cat("  Saved: discovery_MSPI_clinical_associations_panel.pdf\n")
}

## Combined panel - no dots version
if (length(clin_plots_nd) >= 4) {
  ncol_grid <- min(3, length(clin_plots_nd))
  nrow_grid <- ceiling(length(clin_plots_nd) / ncol_grid)

  cairo_pdf(file.path(output_dir, "discovery_MSPI_clinical_associations_panel_nodots.pdf"),
            width = 5.5 * ncol_grid, height = 4.5 * nrow_grid)
  do.call(gridExtra::grid.arrange, c(clin_plots_nd, ncol = ncol_grid))
  dev.off()
  cat("  Saved: discovery_MSPI_clinical_associations_panel_nodots.pdf\n")
}

## ── 6f. Correlation: MSPI vs continuous variables ───────────────────────

cat("\n  Correlations (Spearman) — MSPI vs cell densities:\n")

cor_vars <- c("CD3_density", "M1_density", "M2_density")
if ("Ki67_num" %in% names(disc_s4)) cor_vars <- c(cor_vars, "Ki67_num")

for (cv in cor_vars) {
  if (cv %in% names(disc_s4) && sum(!is.na(disc_s4[[cv]])) > 20) {
    ct <- cor.test(disc_s4$MSPI, disc_s4[[cv]], method = "spearman",
                   use = "complete.obs")
    cat(sprintf("    MSPI vs %s: rho = %.3f, p = %.3g\n",
                cv, ct$estimate, ct$p.value))
  }
}

## ── 6g. gtsummary table: CHANGE 1 - quadrant split (primary) ────────────

cat("\n  Building Discovery cohort gtsummary table (MSPI quadrant split)...\n")

## Compute quadrant cutpoints using inverted scores (higher = more affinity)
m1c <- median(disc_s4$M1_score, na.rm = TRUE)
m2c <- median(disc_s4$M2_score, na.rm = TRUE)

disc_s4$MSPI_quad <- factor(
  ifelse(disc_s4$M1_score >= m1c & disc_s4$M2_score <  m2c, "Q1 Best",
  ifelse(disc_s4$M1_score >= m1c & disc_s4$M2_score >= m2c, "Q2 Mixed",
  ifelse(disc_s4$M1_score <  m1c & disc_s4$M2_score <  m2c, "Q3 Interm",
         "Q4 Worst"))),
  levels = c("Q1 Best", "Q2 Mixed", "Q3 Interm", "Q4 Worst")
)

cat("  MSPI quadrant distribution (Discovery cohort):\n")
print(table(disc_s4$MSPI_quad, useNA = "ifany"))

## Select variables for the table
tbl_vars <- c("Age_cat_2", "Grade", "Stage_TNM_44",
              "Subtype_merged", "Syst_treat_NE", "diadtm_final")
if ("ER10" %in% names(disc_s4))              tbl_vars <- c(tbl_vars, "ER10")
if ("HER2_IHC_dicotome" %in% names(disc_s4)) tbl_vars <- c(tbl_vars, "HER2_IHC_dicotome")

tbl_vars <- tbl_vars[tbl_vars %in% names(disc_s4)]

## Build quadrant table
tbl_labels <- list(
  Age_cat_2       ~ "Age",
  Grade           ~ "Histological grade",
  Stage_TNM_44    ~ "TNM stage",
  Subtype_merged  ~ "Molecular subtype",
  Syst_treat_NE   ~ "Systemic treatment",
  diadtm_final    ~ "Diagnosis period"
)
tbl_labels <- tbl_labels[sapply(tbl_labels, function(f) as.character(f[[2]]) %in% tbl_vars)]

tbl_clin_quad <- disc_s4 %>%
  select(all_of(c("MSPI_quad", tbl_vars))) %>%
  tbl_summary(
    by = MSPI_quad,
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = tbl_labels
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct")) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

cat("
Discovery cohort clinical table by MSPI quadrant:\n")
print(tbl_clin_quad)

tbl_csv_quad <- as.data.frame(tbl_clin_quad$table_body)
write.csv(tbl_csv_quad, file.path(output_dir, "discovery_MSPI_clinical_table_quadrants.csv"),
          row.names = FALSE)

if (requireNamespace("flextable", quietly = TRUE)) {
  ft_quad <- tbl_clin_quad %>% as_flex_table() %>% style_flextable()
  doc <- officer::read_docx() %>%
    officer::body_add_par("Discovery cohort Discovery Cohort - Clinicopathological Features by MSPI Quadrant",
                          style = "Normal") %>%
    officer::body_add_par("", style = "Normal") %>%
    flextable::body_add_flextable(ft_quad)
  print(doc, target = file.path(output_dir, "discovery_MSPI_clinical_table_quadrants.docx"))
  cat("  Saved: discovery_MSPI_clinical_table_quadrants.csv / .docx\n")
}

## ── 6h. gtsummary table (SECONDARY): median split for backward compatibility ─

cat("\n  Building Discovery cohort gtsummary table (MSPI median split - secondary)...\n")

disc_s4$MSPI_group <- ifelse(
  disc_s4$MSPI >= median(disc_s4$MSPI, na.rm = TRUE),
  "High MSPI", "Low MSPI"
)
disc_s4$MSPI_group <- factor(disc_s4$MSPI_group, levels = c("Low MSPI", "High MSPI"))

tbl_clin_median <- disc_s4 %>%
  select(all_of(c("MSPI_group", tbl_vars))) %>%
  tbl_summary(
    by = MSPI_group,
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = tbl_labels
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct")) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

tbl_csv_median <- as.data.frame(tbl_clin_median$table_body)
write.csv(tbl_csv_median, file.path(output_dir, "discovery_MSPI_clinical_table_median_split.csv"),
          row.names = FALSE)

if (requireNamespace("flextable", quietly = TRUE)) {
  ft_med <- tbl_clin_median %>% as_flex_table() %>% style_flextable()
  doc <- officer::read_docx() %>%
    officer::body_add_par("Discovery cohort Discovery Cohort - Clinicopathological Features by MSPI Median Split",
                          style = "Normal") %>%
    officer::body_add_par("", style = "Normal") %>%
    flextable::body_add_flextable(ft_med)
  print(doc, target = file.path(output_dir, "discovery_MSPI_clinical_table_median_split.docx"))
  cat("  Saved: discovery_MSPI_clinical_table_median_split.csv / .docx (backward compatibility)\n")
}

## Save clinical association summary
write.csv(clin_summary, file.path(output_dir, "discovery_MSPI_clinical_summary.csv"),
          row.names = FALSE)


## ═══════════════════════════════════════════════════════════════════════════
## 7. Validation cohort VALIDATION COHORT — LOAD AND ANALYSE
## ═══════════════════════════════════════════════════════════════════════════

cat("\n", strrep("═", 70), "\n")
cat("PART 3: Validation cohort VALIDATION COHORT\n")
cat(strrep("═", 70), "\n\n")

## ── 7a. Load Validation cohort data ─────────────────────────────────────────────────

val_base <- "./data/validation"

val_clin <- data.table::fread(file.path(val_base, "NK_BRCA_for_affinity.txt"),
                               data.table = FALSE)
val_aff  <- data.table::fread(
  find_input(val_base, "BREAST_validation_affinity_CD3e_M1M2_myeloidpool_50um.txt",
             "BREAST_ENCIMA_affinity_CD3e_M1M2_myeloidpool_50um.txt"),
  data.table = FALSE
)
val <- dplyr::left_join(val_clin, val_aff, by = "Case")

## Score inversion
val$M1_score <- -val$M1_to_CD3e
val$M2_score <- -val$M2_to_CD3e
val$MSPI     <- val$M2_score - val$M1_score

## Endpoints
val$RFS_months <- suppressWarnings(as.numeric(val$RFS)) / 30.4375
val$status_rfs <- ifelse(val$RFS.event == "Yes", 1L, 0L)
val$OS_months  <- suppressWarnings(as.numeric(val$Time_Diagnosis_Last_followup..OS.)) / 30.4375
val$status_os  <- ifelse(val$Event_last_followup == "Dead", 1L, 0L)

## Cell densities
## The validation CD3 density column is named with a Cyrillic "С" ("СD3"); match on the ending so this works in any locale
cd3_col_val <- grep("D3$", names(val), value = TRUE)[1]
val$CD3_density <- val[[cd3_col_val]]
val$M1_density  <- val[["M1"]]
val$M2_density  <- val[["M2"]]
cat("  Validation cohort CD3 density column:", cd3_col_val, "\n")

## Clinical variables
val$Age_bin <- factor(
  ifelse(val$Age == "<=50", "<=50",
  ifelse(val$Age %in% c(">50-<65", ">=65"), ">50", NA)),
  levels = c(">50", "<=50")
)
val$Grade_val <- factor(
  ifelse(val$Diff_grade %in% c("1", "2"), "1-2",
  ifelse(val$Diff_grade == "3", "3", NA)),
  levels = c("1-2", "3")
)
val$Node_bin <- factor(
  ifelse(val$pN_stage == "NO axil metastasis", "N0",
  ifelse(val$pN_stage == "axil metastasis", "N+", NA)),
  levels = c("N0", "N+")
)

## ER subsets (exclude DCIS)
val_erp <- val %>%
  dplyr::filter(
    Tumor_type_code == "BRCA_ER_pos",
    !is.na(M1_score), !is.na(M2_score),
    !is.na(OS_months), OS_months > 0
  )

val_ern <- val %>%
  dplyr::filter(
    Tumor_type_code == "BRCA_ER_neg",
    !is.na(M1_score), !is.na(M2_score),
    !is.na(OS_months), OS_months > 0
  )

## Age subsets
val_young <- val_erp %>% filter(Age == "<=50")
val_old   <- val_erp %>% filter(Age %in% c(">50-<65", ">=65"))

cat("  Validation cohort ER+ (all): n =", nrow(val_erp), "| OS events:", sum(val_erp$status_os), "\n")
cat("  Validation cohort ER+ Young (<=50y): n =", nrow(val_young), "| OS events:", sum(val_young$status_os), "\n")
cat("  Validation cohort ER+ Old (>50y): n =", nrow(val_old), "| OS events:", sum(val_old$status_os), "\n")
cat("  Validation cohort ER-negative: n =", nrow(val_ern), "| OS events:", sum(val_ern$status_os), "\n")

## ── 7b. Helper: compute MSPI HR per SD for Validation cohort (generic) ──────────────

compute_mspi_hr_gen <- function(d, label, time_col, event_col) {
  n  <- nrow(d)
  ev <- sum(d[[event_col]])
  mspi_sd <- sd(d$MSPI, na.rm = TRUE)

  if (n < 15 || ev < 5 || is.na(mspi_sd) || mspi_sd < 0.001) {
    return(data.frame(
      Subgroup = label, n = n, events = ev,
      HR = NA, lower = NA, upper = NA, p = NA,
      HR_adj = NA, lower_adj = NA, upper_adj = NA, p_adj = NA,
      stringsAsFactors = FALSE
    ))
  }

  d$MSPI_sd <- d$MSPI / mspi_sd

  ## UV
  surv_formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ MSPI_sd"))
  cx_uv <- coxph(surv_formula, data = d)
  s_uv  <- summary(cx_uv)

  ## MV
  has_m1  <- sum(!is.na(d$M1_density))  > 10
  has_m2  <- sum(!is.na(d$M2_density))  > 10
  has_cd3 <- sum(!is.na(d$CD3_density)) > 10

  if (has_cd3 && has_m1 && has_m2) {
    d_mv <- d %>% filter(!is.na(CD3_density), !is.na(M1_density), !is.na(M2_density))
    mv_formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ MSPI_sd + CD3_density + M1_density + M2_density"))
    cx_mv <- coxph(mv_formula, data = d_mv)
  } else if (has_cd3) {
    d_mv <- d %>% filter(!is.na(CD3_density))
    mv_formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ MSPI_sd + CD3_density"))
    cx_mv <- coxph(mv_formula, data = d_mv)
    cat("    WARNING:", label, "- using CD3-only MV (M1/M2 density missing)\n")
  } else {
    cx_mv <- cx_uv
    cat("    WARNING:", label, "- no density columns; MV = UV fallback\n")
  }
  s_mv <- summary(cx_mv)
  idx_mv <- which(rownames(s_mv$coefficients) == "MSPI_sd")

  data.frame(
    Subgroup  = label,
    n         = n,
    events    = ev,
    HR        = s_uv$conf.int[1, 1],
    lower     = s_uv$conf.int[1, 3],
    upper     = s_uv$conf.int[1, 4],
    p         = s_uv$coefficients[1, 5],
    HR_adj    = s_mv$conf.int[idx_mv, 1],
    lower_adj = s_mv$conf.int[idx_mv, 3],
    upper_adj = s_mv$conf.int[idx_mv, 4],
    p_adj     = s_mv$coefficients[idx_mv, 5],
    stringsAsFactors = FALSE
  )
}

## ── 7c. Validation cohort forest plots (OS endpoint) ────────────────────────────────

cat("\n  Computing Validation cohort MSPI HRs per SD...\n")

val_forest <- rbind(
  compute_mspi_hr_gen(val_erp,   "ER+ (all)", "OS_months", "status_os"),
  compute_mspi_hr_gen(val_young, "ER+ Young (<=50y)", "OS_months", "status_os"),
  compute_mspi_hr_gen(val_old,   "ER+ Old (>50y)", "OS_months", "status_os"),
  compute_mspi_hr_gen(val_ern,   "ER-negative", "OS_months", "status_os")
)

## Print summary
cat("
Validation Cohort MSPI HR per SD by ER status/age (OS):\n")
cat(sprintf("%-22s  n=%4s  ev=%3s  UV: HR=%5s (%5s-%5s) p=%8s  |  MV: HR=%5s p=%8s\n",
            "Subgroup", "n", "ev", "HR", "lo", "hi", "p", "HR_adj", "p_adj"))
for (i in seq_len(nrow(val_forest))) {
  r <- val_forest[i, ]
  cat(sprintf("%-22s  n=%4d  ev=%3d  UV: HR=%5.2f (%5.2f-%5.2f) p=%8.2e  |  MV: HR=%5.2f p=%8.2e\n",
              r$Subgroup, r$n, r$events,
              ifelse(is.na(r$HR), NA, r$HR),
              ifelse(is.na(r$lower), NA, r$lower),
              ifelse(is.na(r$upper), NA, r$upper),
              ifelse(is.na(r$p), NA, r$p),
              ifelse(is.na(r$HR_adj), NA, r$HR_adj),
              ifelse(is.na(r$p_adj), NA, r$p_adj)))
}

write.csv(val_forest, file.path(output_dir, "VAL_MSPI_forest_data.csv"),
          row.names = FALSE)

## Save forest table PDF (CHANGE 2)
save_forest_table(val_forest,
                      file.path(output_dir, "VAL_MSPI_forest_data_TABLE"),
                      "Validation Cohort MSPI by ER Status - Forest Data")

## Formatting
val_forest$sig     <- ifelse(is.na(val_forest$p), "NS",
                             ifelse(val_forest$p < 0.05, "Significant", "NS"))
val_forest$sig_adj <- ifelse(is.na(val_forest$p_adj), "NS",
                             ifelse(val_forest$p_adj < 0.05, "Significant", "NS"))

## ── Helper: build Discovery cohort-style combined forest ──────────────────────────────

build_combined_forest <- function(fdat, sg_order, sg_cols,
                                  plot_title, plot_subtitle, out_name) {

  fdat_long <- rbind(
    fdat %>% select(Subgroup, HR, lower, upper, p, sig) %>%
      mutate(Model = "Univariable"),
    fdat %>% select(Subgroup, HR_adj, lower_adj, upper_adj, p_adj, sig_adj) %>%
      rename(HR = HR_adj, lower = lower_adj, upper = upper_adj,
             p = p_adj, sig = sig_adj) %>%
      mutate(Model = "Density-adjusted")
  ) %>%
    filter(!is.na(HR))

  fdat_long$Subgroup <- as.character(fdat_long$Subgroup)
  fdat_long$Label <- paste0(fdat_long$Subgroup, " - ", fdat_long$Model)

  ## Y-positioning
  within_gap  <- 0.6
  between_gap <- 1.0

  ypos   <- c()
  ylabs  <- c()
  ymap   <- list()
  y_cur  <- length(sg_order) * (within_gap + between_gap)
  for (i in seq_along(sg_order)) {
    sg <- sg_order[i]
    if (i > 1) y_cur <- y_cur - between_gap
    ymap[[paste0(sg, " - Univariable")]]      <- y_cur
    ymap[[paste0(sg, " - Density-adjusted")]] <- y_cur - within_gap
    ypos  <- c(ypos, y_cur, y_cur - within_gap)
    ylabs <- c(ylabs, paste0(sg, " - Univariable"),
               paste0(sg, " - Density-adjusted"))
    y_cur <- y_cur - within_gap
  }

  fdat_long$ynum <- sapply(fdat_long$Label, function(lb) {
    if (!is.null(ymap[[lb]])) ymap[[lb]] else NA
  })
  fdat_long <- fdat_long %>% filter(!is.na(ynum))

  ## Colours
  fdat_long$point_col <- ifelse(
    fdat_long$sig == "Significant",
    sg_cols[fdat_long$Subgroup],
    "grey55"
  )

  p <- ggplot(fdat_long,
              aes(x = HR, y = ynum, xmin = lower, xmax = upper,
                  shape = Model)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(height = 0.2, linewidth = 0.9, colour = "#4d4d4d") +
    geom_point(size = 5, colour = fdat_long$point_col) +
    geom_text(aes(label = sprintf("HR=%.2f (p=%s)",
                                  HR, ifelse(p < 0.001, "<0.001",
                                             sprintf("%.3f", p)))),
              hjust = -0.08, size = 3.5, family = "Arial",
              colour = fdat_long$point_col, show.legend = FALSE) +
    scale_y_continuous(breaks = ypos, labels = ylabs,
                       expand = expansion(mult = 0.08)) +
    scale_shape_manual(values = c("Univariable" = 16, "Density-adjusted" = 18),
                       guide = guide_legend(title = NULL)) +
    labs(title = plot_title, subtitle = plot_subtitle,
         x = "Hazard Ratio per SD (95% CI)", y = NULL) +
    theme_bw(base_size = 14, base_family = "Arial") +
    xlim(min(fdat_long$lower, na.rm = TRUE) - 0.05,
         max(fdat_long$upper, na.rm = TRUE) + 0.6) +
    theme(legend.position = "bottom")

  cairo_pdf(file.path(output_dir, out_name),
            width = 6.5, height = 5.0)
  print(p)
  dev.off()
  cat("  Saved:", out_name, "\n")
  invisible(p)
}

## ── Detailed forest: ER+, ER+ Young, ER+ Old, ER- ─────────────────────────

sg_order_enc <- c("ER+ (all)", "ER+ Young (<=50y)", "ER+ Old (>50y)", "ER-negative")
sg_cols_enc  <- c(
  "ER+ (all)"          = "#a67c52",
  "ER+ Young (<=50y)"  = "#d4b86a",
  "ER+ Old (>50y)"     = "#5b7fa0",
  "ER-negative"        = "#c47a7a"
)

build_combined_forest(
  val_forest, sg_order_enc, sg_cols_enc,
  "MSPI and overall survival by ER status",
  "Validation Cohort validation cohort",
  "VAL_MSPI_forest_combined.pdf"
)

## ── 7c2. Simplified forest (All invasive / ER+ / ER-) ────────────────────

val_invasive <- rbind(val_erp, val_ern)

val_forest_s <- rbind(
  compute_mspi_hr_gen(val_invasive, "All invasive", "OS_months", "status_os"),
  compute_mspi_hr_gen(val_erp,      "ER-positive",  "OS_months", "status_os"),
  compute_mspi_hr_gen(val_ern,      "ER-negative",  "OS_months", "status_os")
)

write.csv(val_forest_s, file.path(output_dir, "VAL_MSPI_forest_data_simple.csv"),
          row.names = FALSE)

## Save simplified forest table PDF (CHANGE 2)
save_forest_table(val_forest_s,
                      file.path(output_dir, "VAL_MSPI_forest_data_simple_TABLE"),
                      "Validation Cohort MSPI - Simplified ER Status Analysis")

val_forest_s$sig     <- ifelse(is.na(val_forest_s$p), "NS",
                               ifelse(val_forest_s$p < 0.05, "Significant", "NS"))
val_forest_s$sig_adj <- ifelse(is.na(val_forest_s$p_adj), "NS",
                               ifelse(val_forest_s$p_adj < 0.05, "Significant", "NS"))

sg_order_s <- c("All invasive", "ER-positive", "ER-negative")
sg_cols_s  <- c("All invasive" = "#a67c52", "ER-positive" = "#5b7fa0", "ER-negative" = "#c47a7a")

build_combined_forest(
  val_forest_s, sg_order_s, sg_cols_s,
  "MSPI and overall survival by ER status",
  "Validation Cohort validation cohort",
  "VAL_MSPI_forest_combined_simple.pdf"
)

## ── 7d. Validation cohort MSPI boxplots (ER+ vs ER-) ───────────────────────────────

val_all <- rbind(
  val_erp %>% mutate(ER_status = "ER+"),
  val_ern %>% mutate(ER_status = "ER-")
) %>% mutate(ER_status = factor(ER_status, levels = c("ER+", "ER-")))

er_colors <- c("ER+" = "#4e79a7", "ER-" = "#c47a7a")

## Wilcoxon test
val_wt_mspi <- wilcox.test(MSPI ~ ER_status, data = val_all)
val_wt_mspi_p <- val_wt_mspi$p.value
val_wt_mspi_star <- ifelse(val_wt_mspi_p < 0.0001, "****",
                    ifelse(val_wt_mspi_p < 0.001, "***",
                    ifelse(val_wt_mspi_p < 0.01, "**",
                    ifelse(val_wt_mspi_p < 0.05, "*", "ns"))))

val_wt_bracket <- data.frame(
  group1 = "ER+", group2 = "ER-",
  y.position = max(val_all$MSPI, na.rm = TRUE) * 1.05,
  p.adj.signif = val_wt_mspi_star,
  stringsAsFactors = FALSE
)

p_enc_mspi_box <- ggplot(val_all, aes(x = ER_status, y = MSPI, fill = ER_status)) +
  geom_jitter(width = 0.2, alpha = 0.15, size = 0.8, colour = "grey40") +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.80) +
  scale_fill_manual(values = er_colors, guide = "none") +
  stat_pvalue_manual(
    val_wt_bracket,
    label = "p.adj.signif", size = 4, bracket.size = 0.5, tip.length = 0.01
  ) +
  labs(
    title = "MSPI by ER status",
    subtitle = sprintf("Validation Cohort (n = %d) | Wilcoxon p = %.3g",
                        nrow(val_all), val_wt_mspi_p),
    x = NULL, y = "MSPI (M2 score - M1 score)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    axis.text.x   = element_text(size = 12)
  )

cairo_pdf(file.path(output_dir, "VAL_MSPI_boxplot_ER.pdf"),
          width = 5.5, height = 5)
print(p_enc_mspi_box)
dev.off()
cat("  Saved: VAL_MSPI_boxplot_ER.pdf\n")

## ── 7d-nodots. Validation cohort MSPI boxplot (no individual dots) ─────────────────

bp_enc <- boxplot.stats(val_all$MSPI)
ypad_enc <- diff(range(bp_enc$stats)) * 0.1

p_enc_mspi_box_nd <- ggplot(val_all, aes(x = ER_status, y = MSPI, fill = ER_status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.80) +
  coord_cartesian(ylim = c(bp_enc$stats[1] - ypad_enc, bp_enc$stats[5] + ypad_enc * 3)) +
  scale_fill_manual(values = er_colors, guide = "none") +
  stat_pvalue_manual(
    val_wt_bracket,
    label = "p.adj.signif", size = 4, bracket.size = 0.5, tip.length = 0.01
  ) +
  labs(
    title = "MSPI by ER status",
    subtitle = sprintf("Validation Cohort (n = %d) | Wilcoxon p = %.3g",
                        nrow(val_all), val_wt_mspi_p),
    x = NULL, y = "MSPI (M2 score - M1 score)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    axis.text.x   = element_text(size = 12)
  )

cairo_pdf(file.path(output_dir, "VAL_MSPI_boxplot_ER_nodots.pdf"),
          width = 5.5, height = 5)
print(p_enc_mspi_box_nd)
dev.off()
cat("  Saved: VAL_MSPI_boxplot_ER_nodots.pdf\n")

## ── 7d2. M1/M2 affinity boxplots by ER status ──────────────────────────

val_score_long <- val_all %>%
  select(Case, ER_status, M1_score, M2_score) %>%
  tidyr::pivot_longer(cols = c(M1_score, M2_score),
                      names_to = "Marker", values_to = "Score") %>%
  mutate(Marker = ifelse(Marker == "M1_score", "M1 affinity", "M2 affinity"))

p_enc_m1m2 <- ggplot(val_score_long,
                      aes(x = ER_status, y = Score, fill = ER_status)) +
  geom_jitter(width = 0.2, alpha = 0.10, size = 0.5, colour = "grey40") +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.80) +
  facet_wrap(~ Marker, scales = "free_y") +
  scale_fill_manual(values = er_colors, guide = "none") +
  labs(
    title = "M1 and M2 affinity scores by ER status",
    subtitle = "Validation Cohort validation cohort",
    x = NULL, y = "Affinity score (inverted, higher = closer)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    strip.text    = element_text(size = 12, face = "bold"),
    axis.text.x   = element_text(size = 11)
  )

cairo_pdf(file.path(output_dir, "VAL_MSPI_boxplot_M1_M2_ER.pdf"),
          width = 9, height = 5)
print(p_enc_m1m2)
dev.off()
cat("  Saved: VAL_MSPI_boxplot_M1_M2_ER.pdf\n")

## ── 7d2-nodots. M1/M2 affinity boxplots (no individual dots) ───────────

p_enc_m1m2_nd <- ggplot(val_score_long,
                         aes(x = ER_status, y = Score, fill = ER_status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.80) +
  facet_wrap(~ Marker, scales = "free_y") +
  scale_fill_manual(values = er_colors, guide = "none") +
  labs(
    title = "M1 and M2 affinity scores by ER status",
    subtitle = "Validation Cohort validation cohort",
    x = NULL, y = "Affinity score (inverted, higher = closer)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    strip.text    = element_text(size = 12, face = "bold"),
    axis.text.x   = element_text(size = 11)
  )

cairo_pdf(file.path(output_dir, "VAL_MSPI_boxplot_M1_M2_ER_nodots.pdf"),
          width = 9, height = 5)
print(p_enc_m1m2_nd)
dev.off()
cat("  Saved: VAL_MSPI_boxplot_M1_M2_ER_nodots.pdf\n")

## ── 7e. Validation cohort cell density boxplots (ER+ vs ER-) ───────────────────────

val_dens_long <- val_all %>%
  filter(!is.na(CD3_density), !is.na(M1_density), !is.na(M2_density)) %>%
  select(ER_status, CD3_density, M1_density, M2_density) %>%
  tidyr::pivot_longer(cols = c(CD3_density, M1_density, M2_density),
                      names_to = "Cell_type", values_to = "Density") %>%
  mutate(Cell_type = recode(Cell_type,
                             CD3_density = "CD3+ T cells",
                             M1_density = "M1 macrophages",
                             M2_density = "M2 macrophages"))

p_enc_dens <- ggplot(val_dens_long,
                      aes(x = ER_status, y = Density, fill = ER_status)) +
  geom_jitter(width = 0.2, alpha = 0.10, size = 0.5, colour = "grey40") +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.80) +
  facet_wrap(~ Cell_type, scales = "free_y") +
  scale_fill_manual(values = er_colors, guide = "none") +
  scale_y_sqrt(breaks = c(1, 10, 100, 1000)) +
  labs(
    title = "Cell densities by ER status",
    subtitle = "Validation Cohort validation cohort",
    x = NULL, y = "Cell density (counts per mm2, sqrt scale)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    strip.text    = element_text(size = 12, face = "bold"),
    axis.text.x   = element_text(size = 11)
  )

cairo_pdf(file.path(output_dir, "VAL_Cell_density_boxplot_ER.pdf"),
          width = 10, height = 5)
print(p_enc_dens)
dev.off()
cat("  Saved: VAL_Cell_density_boxplot_ER.pdf\n")

## ── 7f. Validation cohort clinical associations: MSPI vs Age, Grade, Node ──────────

cat("\n  Building Validation cohort clinical association plots...\n")

val_clin_summary <- data.frame()
val_clin_plots <- list()
val_clin_plots_nd <- list()

val_clin_vars <- c("Age_bin", "Grade_val", "Node_bin")
val_clin_vars <- val_clin_vars[val_clin_vars %in% names(val_erp)]

for (var in val_clin_vars) {
  val_clin_filt <- val_erp %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(group = factor(.data[[var]]))

  if (nrow(val_clin_filt) < 20) next

  wt <- tryCatch(
    kruskal_test(MSPI ~ group, data = val_clin_filt),
    error = function(e) list(p = NA)
  )

  p_base <- ggplot(val_clin_filt, aes(x = group, y = MSPI)) +
    labs(
      title = sprintf("MSPI by %s", var),
      subtitle = sprintf("KW p = %.3g", wt$p),
      x = NULL, y = "MSPI"
    ) +
    theme_bw(base_size = 11, base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )

  p_with_dots <- p_base +
    geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7)

  val_clin_plots[[var]] <- p_with_dots

  p_no_dots <- p_base +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7)

  val_clin_plots_nd[[var]] <- p_no_dots

  val_clin_summary <- rbind(val_clin_summary,
    data.frame(Variable = var, p_value = wt$p, stringsAsFactors = FALSE))
}

for (nm in names(val_clin_plots)) {
  fname <- paste0("VAL_MSPI_clinical_", nm, ".pdf")
  cairo_pdf(file.path(output_dir, fname), width = 6, height = 4.5)
  print(val_clin_plots[[nm]])
  dev.off()
}

for (nm in names(val_clin_plots_nd)) {
  fname <- paste0("VAL_MSPI_clinical_", nm, "_nodots.pdf")
  cairo_pdf(file.path(output_dir, fname), width = 6, height = 4.5)
  print(val_clin_plots_nd[[nm]])
  dev.off()
}

if (length(val_clin_plots) >= 2) {
  ncol_grid <- min(3, length(val_clin_plots))
  nrow_grid <- ceiling(length(val_clin_plots) / ncol_grid)

  cairo_pdf(file.path(output_dir, "VAL_MSPI_clinical_associations_panel.pdf"),
            width = 5.5 * ncol_grid, height = 4.5 * nrow_grid)
  do.call(gridExtra::grid.arrange, c(val_clin_plots, ncol = ncol_grid))
  dev.off()
  cat("  Saved: VAL_MSPI_clinical_associations_panel.pdf\n")
}

## ── 7g. Validation cohort correlations ──────────────────────────────────────────────

cat("\n  Correlations (Spearman) - Validation cohort ER+ MSPI vs cell densities:\n")

val_cor_vars <- c("CD3_density", "M1_density", "M2_density")
for (cv in val_cor_vars) {
  if (cv %in% names(val_erp) && sum(!is.na(val_erp[[cv]])) > 20) {
    ct <- cor.test(val_erp$MSPI, val_erp[[cv]], method = "spearman",
                   use = "complete.obs")
    cat(sprintf("    MSPI vs %s: rho = %.3f, p = %.3g\n",
                cv, ct$estimate, ct$p.value))
  }
}

## ── 7h. gtsummary table: CHANGE 1 - quadrant split (primary) ─────────────

cat("\n  Building Validation cohort gtsummary table (MSPI quadrant split)...\n")

## Compute quadrant cutpoints
m1c_enc <- median(val_erp$M1_score, na.rm = TRUE)
m2c_enc <- median(val_erp$M2_score, na.rm = TRUE)

val_erp$MSPI_quad <- factor(
  ifelse(val_erp$M1_score >= m1c_enc & val_erp$M2_score <  m2c_enc, "Q1 Best",
  ifelse(val_erp$M1_score >= m1c_enc & val_erp$M2_score >= m2c_enc, "Q2 Mixed",
  ifelse(val_erp$M1_score <  m1c_enc & val_erp$M2_score <  m2c_enc, "Q3 Interm",
         "Q4 Worst"))),
  levels = c("Q1 Best", "Q2 Mixed", "Q3 Interm", "Q4 Worst")
)

cat("  Validation cohort ER+ quadrant distribution:\n")
print(table(val_erp$MSPI_quad, useNA = "ifany"))

val_tbl_vars <- c("Age_bin", "Grade_val", "Node_bin")
val_tbl_vars <- val_tbl_vars[val_tbl_vars %in% names(val_erp)]

val_tbl_labels <- list(
  Age_bin   ~ "Age",
  Grade_val ~ "Histological grade",
  Node_bin  ~ "Nodal status"
)
val_tbl_labels <- val_tbl_labels[sapply(val_tbl_labels, function(f) as.character(f[[2]]) %in% val_tbl_vars)]

val_tbl_clin_quad <- val_erp %>%
  select(all_of(c("MSPI_quad", val_tbl_vars))) %>%
  tbl_summary(
    by = MSPI_quad,
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = val_tbl_labels
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct")) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

cat("
Validation Cohort ER+ clinical table by MSPI quadrant:\n")
print(val_tbl_clin_quad)

val_tbl_csv_quad <- as.data.frame(val_tbl_clin_quad$table_body)
write.csv(val_tbl_csv_quad, file.path(output_dir, "VAL_MSPI_clinical_table_quadrants.csv"),
          row.names = FALSE)

if (requireNamespace("flextable", quietly = TRUE)) {
  ft_eq <- val_tbl_clin_quad %>% as_flex_table() %>% style_flextable()
  doc <- officer::read_docx() %>%
    officer::body_add_par("Validation Cohort Validation Cohort (ER+) - Clinicopathological Features by MSPI Quadrant",
                          style = "Normal") %>%
    officer::body_add_par("", style = "Normal") %>%
    flextable::body_add_flextable(ft_eq)
  print(doc, target = file.path(output_dir, "VAL_MSPI_clinical_table_quadrants.docx"))
  cat("  Saved: VAL_MSPI_clinical_table_quadrants.csv / .docx\n")
}

## ── 7h2. gtsummary table: All Invasive (ER+ + ER-) by MSPI quadrant ─────────

cat("\n  Building Validation cohort All Invasive gtsummary table (MSPI quadrant split)...\n")

## val_invasive already created above (line ~1392) as rbind(val_erp, val_ern)
## Add ER status column for the all-invasive table
val_invasive$ER_status <- factor(
  ifelse(val_invasive$Tumor_type_code == "BRCA_ER_pos", "ER+", "ER-"),
  levels = c("ER+", "ER-")
)

## Compute quadrant cutpoints on the ALL-INVASIVE population (own medians)
m1c_inv <- median(val_invasive$M1_score, na.rm = TRUE)
m2c_inv <- median(val_invasive$M2_score, na.rm = TRUE)

val_invasive$MSPI_quad <- factor(
  ifelse(val_invasive$M1_score >= m1c_inv & val_invasive$M2_score <  m2c_inv, "Q1 Best",
  ifelse(val_invasive$M1_score >= m1c_inv & val_invasive$M2_score >= m2c_inv, "Q2 Mixed",
  ifelse(val_invasive$M1_score <  m1c_inv & val_invasive$M2_score <  m2c_inv, "Q3 Interm",
         "Q4 Worst"))),
  levels = c("Q1 Best", "Q2 Mixed", "Q3 Interm", "Q4 Worst")
)

cat("  Validation cohort All Invasive quadrant distribution:\n")
print(table(val_invasive$MSPI_quad, useNA = "ifany"))

inv_tbl_vars <- c("ER_status", "Age_bin", "Grade_val", "Node_bin")
inv_tbl_vars <- inv_tbl_vars[inv_tbl_vars %in% names(val_invasive)]

inv_tbl_labels <- list(
  ER_status ~ "ER status",
  Age_bin   ~ "Age",
  Grade_val ~ "Histological grade",
  Node_bin  ~ "Nodal status"
)
inv_tbl_labels <- inv_tbl_labels[sapply(inv_tbl_labels, function(f) as.character(f[[2]]) %in% inv_tbl_vars)]

val_tbl_clin_quad_inv <- val_invasive %>%
  select(all_of(c("MSPI_quad", inv_tbl_vars))) %>%
  tbl_summary(
    by = MSPI_quad,
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = inv_tbl_labels
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct")) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

cat("
Validation Cohort All Invasive clinical table by MSPI quadrant:\n")
print(val_tbl_clin_quad_inv)

val_tbl_csv_quad_inv <- as.data.frame(val_tbl_clin_quad_inv$table_body)
write.csv(val_tbl_csv_quad_inv, file.path(output_dir, "VAL_AllInv_MSPI_clinical_table_quadrants.csv"),
          row.names = FALSE)

if (requireNamespace("flextable", quietly = TRUE)) {
  ft_eqi <- val_tbl_clin_quad_inv %>% as_flex_table() %>% style_flextable()
  doc <- officer::read_docx() %>%
    officer::body_add_par("Validation Cohort Validation Cohort (All Invasive) - Clinicopathological Features by MSPI Quadrant",
                          style = "Normal") %>%
    officer::body_add_par("", style = "Normal") %>%
    flextable::body_add_flextable(ft_eqi)
  print(doc, target = file.path(output_dir, "VAL_AllInv_MSPI_clinical_table_quadrants.docx"))
  cat("  Saved: VAL_AllInv_MSPI_clinical_table_quadrants.csv / .docx\n")
}

## ── 7h3. Median split table for All Invasive ────────────────────────────────

val_invasive$MSPI <- val_invasive$M2_score - val_invasive$M1_score
val_invasive$MSPI_group <- factor(
  ifelse(val_invasive$MSPI >= median(val_invasive$MSPI, na.rm = TRUE),
         "High MSPI", "Low MSPI"),
  levels = c("Low MSPI", "High MSPI")
)

val_tbl_clin_median_inv <- val_invasive %>%
  select(all_of(c("MSPI_group", inv_tbl_vars))) %>%
  tbl_summary(
    by = MSPI_group,
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = inv_tbl_labels
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct")) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

val_tbl_csv_median_inv <- as.data.frame(val_tbl_clin_median_inv$table_body)
write.csv(val_tbl_csv_median_inv, file.path(output_dir, "VAL_AllInv_MSPI_clinical_table_median_split.csv"),
          row.names = FALSE)

if (requireNamespace("flextable", quietly = TRUE)) {
  ft_emi <- val_tbl_clin_median_inv %>% as_flex_table() %>% style_flextable()
  doc <- officer::read_docx() %>%
    officer::body_add_par("Validation Cohort Validation Cohort (All Invasive) - Clinicopathological Features by MSPI Median Split",
                          style = "Normal") %>%
    officer::body_add_par("", style = "Normal") %>%
    flextable::body_add_flextable(ft_emi)
  print(doc, target = file.path(output_dir, "VAL_AllInv_MSPI_clinical_table_median_split.docx"))
  cat("  Saved: VAL_AllInv_MSPI_clinical_table_median_split.csv / .docx\n")
}

## ── 7i. gtsummary table (SECONDARY): median split for backward compatibility ─

cat("\n  Building Validation cohort gtsummary table (MSPI median split - secondary)...\n")

val_erp$MSPI_group <- ifelse(
  val_erp$MSPI >= median(val_erp$MSPI, na.rm = TRUE),
  "High MSPI", "Low MSPI"
)
val_erp$MSPI_group <- factor(val_erp$MSPI_group, levels = c("Low MSPI", "High MSPI"))

val_tbl_clin_median <- val_erp %>%
  select(all_of(c("MSPI_group", val_tbl_vars))) %>%
  tbl_summary(
    by = MSPI_group,
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = val_tbl_labels
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct")) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

val_tbl_csv_median <- as.data.frame(val_tbl_clin_median$table_body)
write.csv(val_tbl_csv_median, file.path(output_dir, "VAL_MSPI_clinical_table_median_split.csv"),
          row.names = FALSE)

if (requireNamespace("flextable", quietly = TRUE)) {
  ft_em <- val_tbl_clin_median %>% as_flex_table() %>% style_flextable()
  doc <- officer::read_docx() %>%
    officer::body_add_par("Validation Cohort Validation Cohort (ER+) - Clinicopathological Features by MSPI Median Split",
                          style = "Normal") %>%
    officer::body_add_par("", style = "Normal") %>%
    flextable::body_add_flextable(ft_em)
  print(doc, target = file.path(output_dir, "VAL_MSPI_clinical_table_median_split.docx"))
  cat("  Saved: VAL_MSPI_clinical_table_median_split.csv / .docx (backward compatibility)\n")
}

write.csv(val_clin_summary, file.path(output_dir, "VAL_MSPI_clinical_summary.csv"),
          row.names = FALSE)


## ═══════════════════════════════════════════════════════════════════════════
## 8. SUPPLEMENTARY TABLES (S2 Discovery, S3 Validation)
##    Matching format of manually created Table S1 (Exploratory) and Table S2
##    (Discovery) in Figures/Supplementary/.
##    Table S1 = patient-level listing (Exploratory, n=17) - requires manual
##    clinical data not in the analysis pipeline; keep existing manual version.
##    Table S2 = aggregated clinicopathological summary (Discovery/Discovery cohort)
##    Table S3 = aggregated clinicopathological summary (Validation/Validation cohort)
## ═══════════════════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("SUPPLEMENTARY TABLES\n")
cat(strrep("=", 70), "\n\n")

## ── 8a. Table S2: Discovery cohort (Discovery cohort) clinicopathological summary ──────

cat("  Building Table S2 (Discovery cohort summary)...\n")

## Use disc_raw (full cohort before S4 filtering) for the demographic table
## Ensure factor levels
disc_tbl <- disc_raw
disc_tbl$Age_cat_2 <- factor(as.character(disc_tbl$Age_cat_2))
disc_tbl$Stage_TNM_44 <- factor(as.character(disc_tbl$Stage_TNM_44))
disc_tbl$Subtype <- factor(as.character(disc_tbl$Subtype))
if ("Syst_treat_NE" %in% names(disc_tbl)) {
  disc_tbl$Syst_treat_NE <- factor(as.character(disc_tbl$Syst_treat_NE))
}
if ("diadtm_final" %in% names(disc_tbl)) {
  disc_tbl$diadtm_final <- factor(as.character(disc_tbl$diadtm_final))
}
if ("Grade" %in% names(disc_tbl)) {
  disc_tbl$Grade <- factor(as.character(disc_tbl$Grade))
}

## Build two-column summary: Characteristic | Cases (percentage), N=total
s2_vars <- c("Age_cat_2", "Stage_TNM_44", "Subtype")
if ("Grade" %in% names(disc_tbl))          s2_vars <- c(s2_vars, "Grade")
if ("Syst_treat_NE" %in% names(disc_tbl))  s2_vars <- c(s2_vars, "Syst_treat_NE")
if ("diadtm_final" %in% names(disc_tbl))   s2_vars <- c(s2_vars, "diadtm_final")
s2_vars <- s2_vars[s2_vars %in% names(disc_tbl)]

s2_labels <- list(
  Age_cat_2      ~ "Age",
  Stage_TNM_44   ~ "TNM stage",
  Subtype  ~ "Molecular subtype",
  Grade          ~ "Histological grade",
  Syst_treat_NE  ~ "Systemic treatment",
  diadtm_final   ~ "Diagnosis period"
)
s2_labels <- s2_labels[sapply(s2_labels, function(f) as.character(f[[2]]) %in% s2_vars)]

tbl_s2 <- disc_tbl %>%
  select(all_of(s2_vars)) %>%
  tbl_summary(
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = s2_labels
  ) %>%
  bold_labels()

cat("\nTable S2 - Discovery cohort (Discovery cohort):\n")
print(tbl_s2)

## Save as CSV
tbl_s2_csv <- as.data.frame(tbl_s2$table_body)
write.csv(tbl_s2_csv, file.path(output_dir, "Table_S2_Discovery_clinicopath.csv"),
          row.names = FALSE)

## Save as docx
tryCatch({
  tbl_s2 %>%
    as_flex_table() %>%
    style_flextable() %>%
    {
      doc <- officer::read_docx()
      doc <- doc %>%
        officer::body_add_par(
          "Table S2. Distribution of clinicopathological parameters in tumors within discovery cohort.",
          style = "Normal") %>%
        officer::body_add_par("", style = "Normal") %>%
        flextable::body_add_flextable(.)
      print(doc, target = file.path(output_dir, "Table_S2_Discovery_clinicopath.docx"))
    }
  cat("  Saved: Table_S2_Discovery_clinicopath.docx\n")
}, error = function(e) {
  cat("  WARNING: Table S2 docx failed:", conditionMessage(e), "\n")
})


## ── 8b. Table S3: Validation cohort (Validation cohort) clinicopathological summary ──

cat("\n  Building Table S3 (Validation cohort summary)...\n")

## Use the full Validation cohort dataset (all tumor types, not just ER+ subset)
## Add ER status column for the summary
val_tbl <- val
val_tbl$ER_status <- factor(
  ifelse(val_tbl$Tumor_type_code == "BRCA_ER_pos", "ER-positive",
  ifelse(val_tbl$Tumor_type_code == "BRCA_ER_neg", "ER-negative",
  ifelse(grepl("DCIS", val_tbl$Tumor_type_code), "DCIS", "Other"))),
  levels = c("ER-positive", "ER-negative", "DCIS", "Other")
)

## Prepare clinical variables
val_tbl$Age_cat <- factor(val_tbl$Age)
if ("Diff_grade" %in% names(val_tbl)) {
  val_tbl$Grade_enc2 <- factor(as.character(val_tbl$Diff_grade))
}
if ("pN_stage" %in% names(val_tbl)) {
  val_tbl$Node_status <- factor(as.character(val_tbl$pN_stage))
}
if ("pT_stage" %in% names(val_tbl)) {
  val_tbl$T_stage <- factor(as.character(val_tbl$pT_stage))
}

## Filter to MSPI-available (both M1 and M2 affinity scores present)
val_tbl_s3 <- val_tbl %>%
  filter(!is.na(M1_to_CD3e) & !is.na(M2_to_CD3e))

s3_vars <- c("ER_status", "Age_cat")
if ("Grade_enc2" %in% names(val_tbl_s3))  s3_vars <- c(s3_vars, "Grade_enc2")
if ("Node_status" %in% names(val_tbl_s3)) s3_vars <- c(s3_vars, "Node_status")
if ("T_stage" %in% names(val_tbl_s3))     s3_vars <- c(s3_vars, "T_stage")
s3_vars <- s3_vars[s3_vars %in% names(val_tbl_s3)]

s3_labels <- list(
  ER_status   ~ "ER status",
  Age_cat     ~ "Age group",
  Grade_enc2  ~ "Histological grade",
  Node_status ~ "Nodal status",
  T_stage     ~ "T stage"
)
s3_labels <- s3_labels[sapply(s3_labels, function(f) as.character(f[[2]]) %in% s3_vars)]

tbl_s3 <- val_tbl_s3 %>%
  select(all_of(s3_vars)) %>%
  tbl_summary(
    type = list(everything() ~ "categorical"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    label = s3_labels
  ) %>%
  bold_labels()

cat("\nTable S3 - Validation cohort (Validation cohort):\n")
print(tbl_s3)

## Save as CSV
tbl_s3_csv <- as.data.frame(tbl_s3$table_body)
write.csv(tbl_s3_csv, file.path(output_dir, "Table_S3_VAL_clinicopath.csv"),
          row.names = FALSE)

## Save as docx
tryCatch({
  tbl_s3 %>%
    as_flex_table() %>%
    style_flextable() %>%
    {
      doc <- officer::read_docx()
      doc <- doc %>%
        officer::body_add_par(
          "Table S3. Distribution of clinicopathological parameters in the validation cohort.",
          style = "Normal") %>%
        officer::body_add_par("", style = "Normal") %>%
        flextable::body_add_flextable(.)
      print(doc, target = file.path(output_dir, "Table_S3_VAL_clinicopath.docx"))
    }
  cat("  Saved: Table_S3_VAL_clinicopath.docx\n")
}, error = function(e) {
  cat("  WARNING: Table S3 docx failed:", conditionMessage(e), "\n")
})


## ═══════════════════════════════════════════════════════════════════════════
## 9. SUMMARY AND SAVE
## ═══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("================================================================\n")
cat("  ALL OUTPUTS SAVED TO:\n")
cat("  ", output_dir, "\n")
cat("================================================================\n")
cat("\n  Discovery cohort:\n")
cat("    Forest plots: discovery_MSPI_subtype_forest_{UV,MV,combined,combined_merged}.pdf\n")
cat("    Interaction:  discovery_M1M2_interaction_forest.pdf/.png + _data.csv/.docx\n")
cat("    Forest data:  discovery_MSPI_subtype_forest_data{,_merged}.csv + _TABLE.docx/.pdf\n")
cat("    Boxplots:     discovery_MSPI_boxplot_{5subtypes,3subtypes,M1_M2_by_subtype}.pdf\n")
cat("    Densities:    discovery_Cell_density_boxplot_3subtypes.pdf\n")
cat("    Clinical:     discovery_MSPI_clinical_*.pdf (individual + panel)\n")
cat("    Tables:       discovery_MSPI_clinical_table_quadrants.csv/.docx (v2 primary)\n")
cat("                  discovery_MSPI_clinical_table_median_split.csv/.docx (secondary)\n")
cat("    Supp Table:   Table_S2_Discovery_clinicopath.csv/.docx\n")
cat("\n  Validation cohort:\n")
cat("    Forest plots: VAL_MSPI_forest_combined{,_simple}.pdf\n")
cat("    Forest data:  VAL_MSPI_forest_data{,_simple}.csv + _TABLE.docx/.pdf\n")
cat("    Boxplots:     VAL_MSPI_boxplot_{ER,M1_M2_ER}.pdf\n")
cat("    Densities:    VAL_Cell_density_boxplot_ER.pdf\n")
cat("    Clinical:     VAL_MSPI_clinical_*.pdf (individual + panel)\n")
cat("    Tables:       VAL_MSPI_clinical_table_quadrants.csv/.docx (v2 primary)\n")
cat("                  VAL_MSPI_clinical_table_median_split.csv/.docx (secondary)\n")
cat("    Supp Table:   Table_S3_VAL_clinicopath.csv/.docx\n")
cat("\n  ALL tables saved as .docx (primary) with .csv backup.\n")
cat("  Forest tables also have .pdf backup for figure assembly.\n")
cat("================================================================\n")
