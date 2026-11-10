################################################################################
##
##  INTEGRATED SPATIAL IMMUNE ANALYSIS
##  Discovery cohort + Validation cohort ER+
##
##  VERSION v8 (2026-03-25)
##
##  Changes from v7:
##    1. Fixed paths (case-sensitive: ./data/)
##    2. All CSV tables now also saved as styled .docx (flextable + officer)
##    3. Added flextable, officer to library list
##
################################################################################

## ═══════════════════════════════════════════════════════════════════════════
## STEP 0 — OUTPUT DIRECTORIES AND LOG SETUP
## ═══════════════════════════════════════════════════════════════════════════

output_date <- format(Sys.Date(), "%Y%m%d")

disc_base <- "./data/discovery"



val_base <- "./data/validation"

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
  paste0(output_date, "_Integrated_MSPI_Analysis")
)

for (sub in c("plots", "tables", "validation")) {
  dir.create(file.path(output_dir, sub), showWarnings = FALSE, recursive = TRUE)
}

plots_dir      <- file.path(output_dir, "plots")
tables_dir     <- file.path(output_dir, "tables")
validation_dir <- file.path(output_dir, "validation")
mlog           <- file.path(output_dir, "analysis_log.txt")

writeLines(
  c(
    "================================================================",
    "INTEGRATED SPATIAL IMMUNE ANALYSIS — Discovery + Validation ER+",
    paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "MOD: scores INVERTED (negated) — higher value = higher affinity",
    "MSPI = M2_score - M1_score  (high = suppression dominant = BAD)",
    "4-Quadrant cutpoints: individual score medians (pre-specified)",
    "Discovery endpoint: DDFS | Validation endpoint: OS",
    "Added: Discovery cohort subtype-specific analyses and dedicated Q4 KM outputs",
    "================================================================",
    ""
  ),
  mlog
)

## ═══════════════════════════════════════════════════════════════════════════
## STEP 1 — LIBRARIES
## ═══════════════════════════════════════════════════════════════════════════

cat("STEP 1: Loading libraries...\n")

pkgs <- c(
  "foreign", "data.table", "dplyr", "tidyr", "ggplot2", "ggpubr",
  "survival", "survminer", "openxlsx", "boot", "gridExtra", "scales", "forcats",
  "car",   # v4: for Type II Anova chi-sq in pie charts
  "flextable", "officer"   # v8: for styled .docx table output
)

for (p in pkgs) {
  if (!suppressWarnings(require(p, character.only = TRUE, quietly = TRUE))) {
    install.packages(p, repos = "https://cran.rstudio.com/")
    library(p, character.only = TRUE)
  }
}

if (!suppressWarnings(require("survcomp", character.only = TRUE, quietly = TRUE))) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.rstudio.com/")
  }
  BiocManager::install("survcomp", update = FALSE, ask = FALSE)
  library(survcomp)
}

## v4: Ensure Arial is available via systemfonts/ragg for cairo_pdf
if (requireNamespace("systemfonts", quietly = TRUE)) {
  fam_check <- systemfonts::system_fonts()
  if (any(grepl("Arial", fam_check$family, ignore.case = TRUE))) {
    cat("  Arial found via systemfonts.\n")
  } else {
    cat("  WARNING: Arial not found in system fonts. Falling back to sans-serif.\n")
    FONT_FAM <<- "sans"
  }
} else {
  cat("  systemfonts not available; cairo_pdf should find Arial on macOS.\n")
}

## v6: PDF device selection — try cairo_pdf first, fall back to quartz on macOS
USE_QUARTZ <- FALSE
tryCatch({
  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 1, height = 1)
  plot.new()
  dev.off()
  unlink(tf)
  cat("  cairo_pdf device: OK\n")
}, warning = function(w) {
  if (grepl("cairo", w$message, ignore.case = TRUE)) {
    USE_QUARTZ <<- TRUE
  }
}, error = function(e) {
  USE_QUARTZ <<- TRUE
})
if (USE_QUARTZ) {
  cat("  cairo_pdf unavailable, using quartz PDF device (macOS native, vector output)\n")
}

cat("STEP 1: Complete — libraries loaded\n\n")

## ── v8: docx table helpers ──────────────────────────────────────────────
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

save_table_docx <- function(tab_df, out_path, title = NULL) {
  ft <- flextable::flextable(tab_df) %>% style_flextable()
  doc <- officer::read_docx()
  if (!is.null(title) && nzchar(title)) {
    doc <- doc %>%
      officer::body_add_par(title, style = "Normal") %>%
      officer::body_add_par("", style = "Normal")
  }
  doc <- doc %>% flextable::body_add_flextable(ft)
  print(doc, target = out_path)
}

## ═══════════════════════════════════════════════════════════════════════════
## STEP 2 — LOAD DISCOVERY COHORT CLINICAL DATA (SPSS)
## ═══════════════════════════════════════════════════════════════════════════

cat("STEP 2: Loading Discovery cohort clinical data (SPSS)...\n")

disc_spss_candidates <- c(
  file.path(disc_base, "Discovery_workingr_832_updated_211024.sav"),
  file.path(disc_base, "Matching to clinical data", "Discovery_workingr_832_updated_211024.sav"),
  file.path(disc_base, "BCY_workingr_832_updated_211024.sav"),
  file.path(disc_base, "Matching to clinical data", "BCY_workingr_832_updated_211024.sav")
)

disc_spss <- NA
for (p in disc_spss_candidates) {
  if (file.exists(p)) {
    disc_spss <- p
    break
  }
}
if (is.na(disc_spss)) stop("Cannot find Discovery_workingr_832_updated_211024.sav")

disc_raw <- foreign::read.spss(
  disc_spss,
  use.value.labels = TRUE,
  to.data.frame = TRUE,
  max.value.labels = Inf,
  trim.factor.names = TRUE,
  trim_values = TRUE,
  reencode = NA
)

cat("  Clinical rows:", nrow(disc_raw), "\n")
## The original clinical file calls the subtype variable Subtype_BCY21
if (!"Subtype_Discovery21" %in% names(disc_raw) && "Subtype_BCY21" %in% names(disc_raw)) {
  names(disc_raw)[names(disc_raw) == "Subtype_BCY21"] <- "Subtype_Discovery21"
}

disc_raw$status_num <- ifelse(
  as.character(disc_raw$DDFScens2021) == "Distant recurrence" |
    as.character(disc_raw$BCcens2021) == "Breast cancer death",
  1L, 0L
)
disc_raw$DDFStime2021 <- as.numeric(disc_raw$DDFStime2021)

disc_raw <- disc_raw[as.character(disc_raw$Subtype_Discovery21) != "Unclassified", ]
disc_raw$Subtype_Discovery21 <- droplevels(disc_raw$Subtype_Discovery21)

cat(
  "  After removing unclassified: n=", nrow(disc_raw),
  " | DDFS events:", sum(disc_raw$status_num, na.rm = TRUE), "\n\n"
)
cat("STEP 2: Discovery cohort clinical loaded\n\n", file = mlog, append = TRUE)

## ═══════════════════════════════════════════════════════════════════════════
## STEP 3 — LOAD DISCOVERY COHORT SPATIAL AND DENSITY DATA
## ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: Loading Discovery cohort spatial distance and density files...\n")

find_disc <- function(...) {
  fnames <- c(...)
  candidates <- c(
    file.path(disc_base, fnames),
    file.path(disc_base, "Discovery cohort", fnames),
    file.path(disc_base, "analysis", fnames)
  )
  for (f in candidates) {
    if (file.exists(f)) {
      cat("  Found:", f, "\n")
      return(f)
    }
  }
  stop(paste("Cannot find:", paste(fnames, collapse = " / ")))
}

m1_raw  <- data.table::fread(find_disc("Discovery_distances_mean_M1_to_CD3.txt", "BCY_distances_mean_M1_to_CD3.txt"), data.table = FALSE)
m2_raw  <- data.table::fread(find_disc("Discovery_distances_mean_M2_to_CD3.txt", "BCY_distances_mean_M2_to_CD3.txt"), data.table = FALSE)
cd3_raw <- data.table::fread(find_disc("Discovery_immune_variables_with_CD3_for spatial.txt", "BCY_immune_variables_with_CD3_for spatial.txt"), data.table = FALSE)

colnames(m1_raw)[colnames(m1_raw) == "delta_CD3"] <- "M1_to_CD3_distance"
colnames(m2_raw)[colnames(m2_raw) == "delta_CD3"] <- "M2_to_CD3_distance"
cd3_raw <- cd3_raw[, 1:2]
colnames(cd3_raw) <- c("AnnUnit", "CD3_density")

cat("  M1 rows:", nrow(m1_raw), "| M2 rows:", nrow(m2_raw), "\n\n")
cat("STEP 3: Discovery cohort spatial data loaded\n\n", file = mlog, append = TRUE)

## ═══════════════════════════════════════════════════════════════════════════
## STEP 4 — MERGE DISCOVERY COHORT AND DEFINE COHORTS
## ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Merging Discovery cohort data and defining analysis cohorts...\n")

disc <- disc_raw %>%
  dplyr::left_join(m1_raw[, c("AnnUnit", "M1_to_CD3_distance")], by = "AnnUnit") %>%
  dplyr::left_join(m2_raw[, c("AnnUnit", "M2_to_CD3_distance")], by = "AnnUnit") %>%
  dplyr::left_join(cd3_raw, by = "AnnUnit")

disc$M1_score <- -disc$M1_to_CD3_distance
disc$M2_score <- -disc$M2_to_CD3_distance

cat(
  "  Post-merge: n=", nrow(disc),
  "| M1_score non-NA:", sum(!is.na(disc$M1_score)),
  "| M2_score non-NA:", sum(!is.na(disc$M2_score)),
  "| Both:", sum(!is.na(disc$M1_score) & !is.na(disc$M2_score)), "\n"
)

disc_young <- disc %>%
  dplyr::filter(
    Stage_TNM_44 %in% c("Stage I", "Stage IIa", "Stage IIb"),
    Subtype_Discovery21 %in% c("Lum A", "Lum B"),
    Age_cat_2 == "<40",
    !is.na(M1_score), !is.na(M2_score),
    !is.na(DDFStime2021), DDFStime2021 > 0
  )

disc_early <- disc %>%
  dplyr::filter(
    Stage_TNM_44 %in% c("Stage I", "Stage IIa", "Stage IIb"),
    Subtype_Discovery21 %in% c("Lum A", "Lum B"),
    !is.na(M1_score), !is.na(M2_score),
    !is.na(DDFStime2021), DDFStime2021 > 0
  )

disc_s4 <- disc %>%
  dplyr::filter(
    !is.na(M1_score), !is.na(M2_score),
    !is.na(DDFStime2021), DDFStime2021 > 0
  )

for (nm in c("Young Luminal", "Early Luminal", "S4")) {
  d <- switch(
    nm,
    "Young Luminal" = disc_young,
    "Early Luminal" = disc_early,
    "S4" = disc_s4
  )
  cat(sprintf(
    "  Discovery cohort %-18s n=%4d | DDFS ev=%3d (%.1f%%)\n",
    nm, nrow(d), sum(d$status_num), 100 * mean(d$status_num)
  ))
}

## Added: molecular subtype subsets within Discovery cohort S4
## Individual subtypes (kept for per-subtype output CSVs)
disc_subtypes <- c("Lum A", "Lum B", "Lum Her2", "Her2 exp", "Triple neg")
disc_subtype_list <- list()

for (st in disc_subtypes) {
  ds <- disc_s4 %>%
    dplyr::filter(
      as.character(Subtype_Discovery21) == st,
      !is.na(M1_score), !is.na(M2_score),
      !is.na(DDFStime2021), DDFStime2021 > 0
    )
  disc_subtype_list[[st]] <- ds
  cat(sprintf(
    "  Discovery cohort subtype %-10s n=%4d | DDFS ev=%3d (%.1f%%)\n",
    st, nrow(ds), sum(ds$status_num), ifelse(nrow(ds) > 0, 100 * mean(ds$status_num), NA_real_)
  ))
}

## Merged molecular subtype groups (for manuscript-level analyses)
cat("\n  Building merged molecular subtype groups:\n")
disc_merged_list <- list()

## Luminal = Lum A + Lum B
disc_merged_list[["Luminal"]] <- disc_s4 %>%
  dplyr::filter(
    as.character(Subtype_Discovery21) %in% c("Lum A", "Lum B"),
    !is.na(M1_score), !is.na(M2_score),
    !is.na(DDFStime2021), DDFStime2021 > 0
  )

## HER2+ = Her2 exp + Lum Her2
disc_merged_list[["HER2+"]] <- disc_s4 %>%
  dplyr::filter(
    as.character(Subtype_Discovery21) %in% c("Her2 exp", "Lum Her2"),
    !is.na(M1_score), !is.na(M2_score),
    !is.na(DDFStime2021), DDFStime2021 > 0
  )

## TNBC = Triple neg (unchanged)
disc_merged_list[["TNBC"]] <- disc_s4 %>%
  dplyr::filter(
    as.character(Subtype_Discovery21) == "Triple neg",
    !is.na(M1_score), !is.na(M2_score),
    !is.na(DDFStime2021), DDFStime2021 > 0
  )

for (mg in names(disc_merged_list)) {
  ds <- disc_merged_list[[mg]]
  cat(sprintf(
    "  Discovery cohort merged %-10s n=%4d | DDFS ev=%3d (%.1f%%)\n",
    mg, nrow(ds), sum(ds$status_num), ifelse(nrow(ds) > 0, 100 * mean(ds$status_num), NA_real_)
  ))
}

## Prepare Discovery cohort clinical variables for Task 6 (Clinical MV)
## Using original column names: Age_cat_2, Grade, Stage_TNM_44, Syst_treat_NE, diadtm_final
cat("  Preparing Discovery cohort clinical variables for multivariable analysis...\n")
cat("  Using raw Discovery cohort columns: Age_cat_2, Grade, Stage_TNM_44, Syst_treat_NE, diadtm_final\n")

## Ensure each is a factor (they may already be from SPSS import)
for (v in c("Age_cat_2", "Grade", "Stage_TNM_44", "Syst_treat_NE", "diadtm_final")) {
  if (v %in% names(disc_s4)) {
    disc_s4[[v]] <- factor(as.character(disc_s4[[v]]))
    cat(sprintf("    %s: %d non-NA, levels: %s\n",
                v, sum(!is.na(disc_s4[[v]])),
                paste(levels(disc_s4[[v]]), collapse = ", ")))
  } else {
    cat(sprintf("    %s: NOT FOUND in dataset\n", v))
  }
}

## Propagate to subsets (disc_young, disc_early)
for (nm in c("disc_young", "disc_early")) {
  d <- get(nm)
  idx <- match(d$AnnUnit, disc_s4$AnnUnit)
  for (v in c("Age_cat_2", "Grade", "Stage_TNM_44", "Syst_treat_NE", "diadtm_final")) {
    if (v %in% names(disc_s4)) d[[v]] <- disc_s4[[v]][idx]
  }
  assign(nm, d)
}

## Propagate to merged subtype subsets (disc_merged_list)
for (mg in names(disc_merged_list)) {
  ds <- disc_merged_list[[mg]]
  idx <- match(ds$AnnUnit, disc_s4$AnnUnit)
  for (v in c("Age_cat_2", "Grade", "Stage_TNM_44", "Syst_treat_NE", "diadtm_final")) {
    if (v %in% names(disc_s4)) ds[[v]] <- disc_s4[[v]][idx]
  }
  disc_merged_list[[mg]] <- ds
}
cat("  Clinical variables propagated to disc_young, disc_early, and disc_merged_list\n")

cat("STEP 4: Discovery cohort cohorts defined\n\n", file = mlog, append = TRUE)

## ═══════════════════════════════════════════════════════════════════════════
## STEP 5 — LOAD VALIDATION COHORT DATA
## ═══════════════════════════════════════════════════════════════════════════

cat("STEP 5: Loading Validation cohort clinical and affinity data...\n")

val_clin <- data.table::fread(file.path(val_base, "NK_BRCA_for_affinity.txt"), data.table = FALSE)
val_aff  <- data.table::fread(
  find_input(val_base, "BREAST_validation_affinity_CD3e_M1M2_myeloidpool_50um.txt",
             "BREAST_ENCIMA_affinity_CD3e_M1M2_myeloidpool_50um.txt"),
  data.table = FALSE
)
val <- dplyr::left_join(val_clin, val_aff, by = "Case")

val$M1_score <- -val$M1_to_CD3e
val$M2_score <- -val$M2_to_CD3e

val$RFS_months <- suppressWarnings(as.numeric(val$RFS)) / 30.4375
val$status_rfs <- ifelse(val$RFS.event == "Yes", 1L, 0L)
val$OS_months  <- suppressWarnings(as.numeric(val$Time_Diagnosis_Last_followup..OS.)) / 30.4375
val$status_os  <- ifelse(val$Event_last_followup == "Dead", 1L, 0L)

## The validation CD3 density column is named with a Cyrillic "С" ("СD3"); match on the ending so this works in any locale
cd3_col <- grep("D3$", names(val), value = TRUE)[1]
val$CD3_density <- val[[cd3_col]]
cat("  CD3 density column identified as:", cd3_col, "\n")

val$Grade <- factor(
  ifelse(val$Diff_grade %in% c("1", "2", "3"), val$Diff_grade, NA),
  levels = c("1", "2", "3")
)
val$Node_pos <- ifelse(
  val$pN_stage == "axil metastasis", 1L,
  ifelse(val$pN_stage == "NO axil metastasis", 0L, NA_integer_)
)

val_erp <- val %>%
  dplyr::filter(
    Tumor_type_code == "BRCA_ER_pos",
    !is.na(M1_score), !is.na(M2_score),
    !is.na(OS_months), OS_months > 0
  )

cat(sprintf("  Validation cohort ER+ (both scores, valid OS): n=%d\n", nrow(val_erp)))
cat(sprintf(
  "  OS ev=%d (%.1f%%)\n",
  sum(val_erp$status_os), 100 * mean(val_erp$status_os)
))

## Validation cohort ER+ age subgroups
val_young <- val_erp %>% dplyr::filter(Age == "<=50")
val_oldb   <- val_erp %>% dplyr::filter(Age %in% c(">50-<65", ">=65"))
cat(sprintf("  Validation cohort ER+ Young (<=50y):  n=%d  OS ev=%d\n",
            nrow(val_young), sum(val_young$status_os)))
cat(sprintf("  Validation cohort ER+ Old   (>50y):   n=%d  OS ev=%d\n",
            nrow(val_oldb), sum(val_oldb$status_os)))

## Validation cohort ER+ node-negative subsets (N0 only)
val_erp_N0 <- val_erp %>% dplyr::filter(Node_pos == 0L)
val_young_N0 <- val_erp_N0 %>% dplyr::filter(Age == "<=50")
cat(sprintf("  Validation cohort ER+ N0:             n=%d  OS ev=%d\n",
            nrow(val_erp_N0), sum(val_erp_N0$status_os)))
cat(sprintf("  Validation cohort ER+ Young N0:       n=%d  OS ev=%d\n",
            nrow(val_young_N0), sum(val_young_N0$status_os)))

## Validation cohort ER-negative subset (exploratory)
val_ern <- val %>%
  dplyr::filter(
    Tumor_type_code == "BRCA_ER_neg",
    !is.na(M1_score), !is.na(M2_score),
    !is.na(OS_months), OS_months > 0
  )
cat(sprintf("  Validation cohort ER- (both scores, valid OS): n=%d  OS ev=%d (%.1f%%)\n",
            nrow(val_ern), sum(val_ern$status_os), 100 * mean(val_ern$status_os)))

## Validation cohort DCIS subset (exploratory)
val_dcis <- val %>%
  dplyr::filter(
    Tumor_type_code == "DCIS",
    !is.na(M1_score), !is.na(M2_score),
    !is.na(OS_months), OS_months > 0
  )
cat(sprintf("  Validation cohort DCIS (both scores, valid OS): n=%d  OS ev=%d (%.1f%%)\n",
            nrow(val_dcis), sum(val_dcis$status_os), 100 * mean(val_dcis$status_os)))

## Validation cohort all-invasive (ER+ + ER-, excluding DCIS) — for pan-cohort OS analyses
val_inv <- val %>%
  dplyr::filter(
    Tumor_type_code %in% c("BRCA_ER_pos", "BRCA_ER_neg"),
    !is.na(M1_score), !is.na(M2_score),
    !is.na(OS_months), OS_months > 0
  )
cat(sprintf("  Validation cohort all invasive (ER+ + ER-, valid OS): n=%d  OS ev=%d (%.1f%%)\n",
            nrow(val_inv), sum(val_inv$status_os), 100 * mean(val_inv$status_os)))

## Validation cohort all-invasive age subgroups
val_inv_young <- val_inv %>% dplyr::filter(Age == "<=50")
val_inv_old   <- val_inv %>% dplyr::filter(Age %in% c(">50-<65", ">=65"))
cat(sprintf("  Validation cohort All Inv Young (<=50y): n=%d  OS ev=%d\n",
            nrow(val_inv_young), sum(val_inv_young$status_os)))
cat(sprintf("  Validation cohort All Inv Old   (>50y):  n=%d  OS ev=%d\n",
            nrow(val_inv_old), sum(val_inv_old$status_os)))

## Prepare Validation cohort clinical variables for Task 6 (Clinical MV)
cat("  Preparing Validation cohort clinical variables...\n")

val$Age_bin <- factor(
  ifelse(val$Age == "<=50", "<=50",
  ifelse(val$Age %in% c(">50-<65", ">=65"), ">50", NA)),
  levels = c(">50", "<=50")
)

val$Grade_enc <- factor(
  ifelse(val$Diff_grade %in% c("1", "2"), "1-2",
  ifelse(val$Diff_grade == "3", "3", NA)),
  levels = c("1-2", "3")
)

val$Node_bin <- factor(
  ifelse(val$pN_stage == "NO axil metastasis", "N0",
  ifelse(val$pN_stage == "axil metastasis", "N+", NA)),
  levels = c("N0", "N+")
)

cat(sprintf("  Validation cohort clinical completeness (all, n=%d):\n", nrow(val)))
cat(sprintf("    Age_bin:   %d non-NA\n", sum(!is.na(val$Age_bin))))
cat(sprintf("    Grade_enc: %d non-NA\n", sum(!is.na(val$Grade_enc))))
cat(sprintf("    Node_bin:  %d non-NA\n", sum(!is.na(val$Node_bin))))

## Propagate to Validation cohort subsets
for (nm in c("val_erp", "val_young", "val_oldb", "val_erp_N0", "val_young_N0",
              "val_ern", "val_dcis",
              "val_inv", "val_inv_young", "val_inv_old")) {
  if (exists(nm)) {
    d <- get(nm)
    idx <- match(d$Case, val$Case)
    for (v in c("Age_bin", "Grade_enc", "Node_bin")) {
      d[[v]] <- val[[v]][idx]
    }
    assign(nm, d)
  }
}

cat("STEP 5: Validation cohort data loaded (ER+, age subgroups, N0 subsets, ER-, DCIS)\n\n", file = mlog, append = TRUE)

## ═══════════════════════════════════════════════════════════════════════════
## STEP 6 — DOUBLE-CHECK
## ═══════════════════════════════════════════════════════════════════════════

cat("\nSTEP 6 [DOUBLE-CHECK]: Distributions, MSPI variance, M1-M2 correlation\n")
cat("STEP 6 double-check\n", file = mlog, append = TRUE)

dbl_check_ok <- TRUE

check_scores <- function(dat, label) {
  m1 <- dat$M1_score
  m2 <- dat$M2_score
  mspi <- m2 - m1
  cat(sprintf("  %s (n=%d):\n", label, nrow(dat)))
  cat(sprintf(
    "    M1_score : mean=%.4f  SD=%.4f  range [%.4f, %.4f]\n",
    mean(m1, na.rm = TRUE), sd(m1, na.rm = TRUE), min(m1, na.rm = TRUE), max(m1, na.rm = TRUE)
  ))
  cat(sprintf(
    "    M2_score : mean=%.4f  SD=%.4f  range [%.4f, %.4f]\n",
    mean(m2, na.rm = TRUE), sd(m2, na.rm = TRUE), min(m2, na.rm = TRUE), max(m2, na.rm = TRUE)
  ))
  cat(sprintf(
    "    MSPI     : mean=%.4f  SD=%.4f  range [%.4f, %.4f]\n",
    mean(mspi, na.rm = TRUE), sd(mspi, na.rm = TRUE), min(mspi, na.rm = TRUE), max(mspi, na.rm = TRUE)
  ))
  ct <- cor.test(m1, m2, method = "spearman", exact = FALSE)
  rho_flag <- if (abs(ct$estimate) > 0.6) {
    dbl_check_ok <<- FALSE
    "HIGH COLLINEARITY - MSPI results may be unreliable"
  } else "collinearity acceptable"
  cat(sprintf(
    "    M1-M2 Spearman rho=%.3f  p=%.4f  %s\n",
    ct$estimate, ct$p.value, rho_flag
  ))
  var_flag <- if (sd(mspi, na.rm = TRUE) < 0.01) {
    dbl_check_ok <<- FALSE
    "NEAR-ZERO VARIANCE"
  } else "adequate variance"
  cat(sprintf("    MSPI SD=%.4f  %s\n", sd(mspi, na.rm = TRUE), var_flag))
  invisible(list(rho = ct$estimate, mspi_sd = sd(mspi, na.rm = TRUE)))
}

chk_by <- check_scores(disc_young, "Discovery Young Luminal")
chk_be <- check_scores(disc_early, "Discovery Early Luminal")
chk_bs <- check_scores(disc_s4,    "Discovery S4 All-Comers")
chk_en <- check_scores(val_erp,   "Validation cohort ER+")

for (st in names(disc_subtype_list)) {
  ds <- disc_subtype_list[[st]]
  if (nrow(ds) >= 5) {
    check_scores(ds, paste0("Discovery cohort subtype ", st))
  }
}

cat(sprintf(
  "\n  Overall double-check: %s\n\n",
  if (dbl_check_ok) "PASSED" else "WARNING - review before interpretation"
))
cat(sprintf(
  "STEP 6: %s\n\n",
  if (dbl_check_ok) "PASSED" else "WARNING"
), file = mlog, append = TRUE)

## ── Global plot settings ─────────────────────────────────────────────────

BASE_SIZE  <- 14          # default ggplot base_size (was 11)
KM_BASE    <- 13          # base_size for KM plots
FONT_FAM   <- "Arial"     # universal font

## ── Muted / old-journal colour palette ──────────────────────────────────
## Pastel-toned versions of the original colours; gentler on the eye.

COL_GREEN   <- "#6aad8e"   # was #1a6b3c  (Q1 best / low MSPI / M1 protective)
COL_RED     <- "#c47a7a"   # was #d7191c  (Q4 worst / high MSPI / M2 risk)
COL_AMBER   <- "#d4b86a"   # was #e6b800  (Q2 mixed)
COL_GREY    <- "#a0a0a0"   # was #888888  (Q3 intermediate / neutral)
COL_BLUE    <- "#7ba3c4"   # was #4a82b0  (general accent, pooled)
COL_NAVY    <- "#5b7fa0"   # was #1d5c8a  (forest plot emphasis)
COL_DKGREEN <- "#5a8a6e"   # was #2e6b3e  (forest individual cohorts)
COL_BROWN   <- "#a67c52"   # was #7a2e0e  (density + MSPI model)
COL_DKGREY  <- "#4d4d4d"   # non-Q4 group

## Named vectors for quick re-use
PAL_4Q <- c(Q1_Best = COL_GREEN, Q2_Mixed = COL_AMBER,
            Q3_Interm = COL_GREY, Q4_Worst = COL_RED)
PAL_2G <- c(COL_GREEN, COL_RED)       # low vs high
PAL_Q4_REST <- c(COL_DKGREY, COL_RED) # non-Q4 vs Q4
PAL_Q12_Q34 <- c(COL_GREEN, COL_RED)  # Q1+Q2 (M1-hi) vs Q3+Q4 (M1-lo)

## Set ggplot2 default theme once: Arial font, larger text
theme_set(
  theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
    theme(
      text            = element_text(family = FONT_FAM),
      plot.title      = element_text(size = BASE_SIZE + 2, face = "bold"),
      plot.subtitle   = element_text(size = BASE_SIZE),
      axis.title      = element_text(size = BASE_SIZE),
      axis.text       = element_text(size = BASE_SIZE - 1),
      legend.text     = element_text(size = BASE_SIZE - 1),
      legend.title    = element_text(size = BASE_SIZE),
      strip.text      = element_text(size = BASE_SIZE)
    )
)

## ── Shared helpers ────────────────────────────────────────────────────────

## PDF-only, vector output.  Default width shrunk so that
## elements (text, ticks, legends) appear larger relative to the page.
## v6: uses quartz PDF on macOS when cairo_pdf DLL is unavailable.
save_pdf <- function(fpath, w, h) {
  if (USE_QUARTZ) {
    quartz(type = "pdf", file = fpath, width = w, height = h,
           family = FONT_FAM, bg = "white")
  } else {
    cairo_pdf(fpath, width = w, height = h, bg = "white")
  }
}

save_plot <- function(p, fname, w = 6, h = 5) {
  fpath <- file.path(plots_dir, paste0(fname, ".pdf"))
  save_pdf(fpath, w, h)
  print(p)
  dev.off()
}

## KM plots — combine plot + risk table into one panel, then save.
save_km <- function(km, fname, w = 5.5, h = 6) {
  g <- ggpubr::ggarrange(
    km$plot   + theme(text = element_text(family = FONT_FAM)),
    km$table  + theme(text = element_text(family = FONT_FAM)),
    ncol = 1, heights = c(0.75, 0.25), align = "v"
  )
  fpath <- file.path(plots_dir, paste0(fname, ".pdf"))
  save_pdf(fpath, w, h)
  print(g)
  dev.off()
}

## ── Mini forest: UV + MV (density-adjusted), compact vertical ───────────
##    Produces a narrow horizontal strip for placement above KM plots.
##    rows = two bars (Univariable / Adjusted).
##    Returns the ggplot object AND saves PDF.

make_mini_forest <- function(uv_cx, mv_cx, mv_var = "MSPI_sd",
                             cohort_label, endpoint_label, file_prefix) {

  extract <- function(cx, lbl) {
    s <- summary(cx)
    i <- which(rownames(s$coefficients) == mv_var)
    if (length(i) == 0) i <- 1
    data.frame(
      Model = lbl,
      HR    = s$conf.int[i, 1],
      CI_lo = s$conf.int[i, 3],
      CI_hi = s$conf.int[i, 4],
      p     = s$coefficients[i, 5],
      stringsAsFactors = FALSE
    )
  }

  rows <- rbind(extract(uv_cx, "Univariable"),
                extract(mv_cx, "Density-adjusted"))
  rows$Model <- factor(rows$Model, levels = c("Density-adjusted", "Univariable"))
  rows$sig <- ifelse(rows$p < 0.05, "Significant", "NS")

  p <- ggplot(rows, aes(x = HR, y = Model, shape = Model, colour = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.22,
                   linewidth = 0.85, colour = COL_DKGREY) +
    geom_point(size = 6) +
    geom_text(aes(label = sprintf("HR %.2f (%.2f-%.2f) p=%s",
                                  HR, CI_lo, CI_hi, fmt_p(p))),
              hjust = -0.08, size = 3.5, family = FONT_FAM, colour = "black") +
    scale_shape_manual(values = c("Univariable" = 16, "Density-adjusted" = 18),
                       guide = "none") +
    scale_colour_manual(values = c("Significant" = "#c0392b", "NS" = "grey55"),
                        guide = "none") +
    scale_x_log10() +
    labs(x = "HR per SD", y = NULL,
         title = paste0(cohort_label, " | ", endpoint_label)) +
    theme_bw(base_size = 12, base_family = FONT_FAM) +
    theme(plot.margin = margin(2, 60, 2, 2),
          axis.text.y = element_text(size = 11),
          plot.title  = element_text(size = 11, face = "bold"),
          panel.grid.major.y = element_blank())

  save_plot(p, paste0(file_prefix, "_mini_forest"), w = 6, h = 1.8)

  ## v8: save mini-forest data as CSV + docx
  rows$p_fmt <- sapply(rows$p, fmt_p)
  write.csv(rows,
    file.path(tables_dir, paste0(file_prefix, "_mini_forest_", output_date, ".csv")),
    row.names = FALSE)
  save_table_docx(rows,
    file.path(tables_dir, paste0(file_prefix, "_mini_forest_", output_date, ".docx")),
    title = paste0("Mini Forest - ", cohort_label, " | ", endpoint_label))

  invisible(p)
}

## ── Pie chart: chi-square proportions in MV Cox ─────────────────────────
##    Uses Type II Anova (car::Anova) to get per-variable chi-sq.
##    Each slice = proportion of total chi-sq; labelled with p-value.

make_mv_pie <- function(cx, var_labels = NULL,
                        cohort_label, file_prefix,
                        palette = NULL) {

  ## Get chi-sq from Type II Anova
  if (!requireNamespace("car", quietly = TRUE)) {
    install.packages("car", repos = "https://cran.rstudio.com/")
  }
  a2 <- car::Anova(cx, type = "II")
  chisq <- a2[, grep("Chisq|LR", colnames(a2))[1]]
  pvals <- a2[, grep("Pr|p", colnames(a2))[1]]
  vars  <- rownames(a2)

  total_chi <- sum(chisq, na.rm = TRUE)
  pct       <- 100 * chisq / total_chi

  ## Build labels
  labels_nice <- if (!is.null(var_labels)) {
    sapply(vars, function(v) {
      ## strip factor level suffixes for lookup
      base_v <- sub("(High|Yes|III|II|N\\+|<=50|3|_sd)$", "", v)
      if (base_v %in% names(var_labels)) var_labels[[base_v]]
      else if (v %in% names(var_labels)) var_labels[[v]]
      else v
    })
  } else vars

  df <- data.frame(
    Variable = labels_nice,
    chisq    = chisq,
    pct      = pct,
    p        = pvals,
    stringsAsFactors = FALSE
  )
  df <- df[order(-df$pct), ]
  df$Variable <- factor(df$Variable, levels = df$Variable)

  if (is.null(palette)) {
    n <- nrow(df)
    palette <- if (n <= 3) c(COL_BLUE, COL_RED, COL_AMBER)[1:n]
               else if (n == 4) c(COL_BLUE, COL_RED, COL_AMBER, COL_GREEN)
               else scales::hue_pal(l = 75)(n)  # light/pastel fallback
  }

  ## Build p-value label
  df$plabel <- sapply(df$p, function(x) {
    if (is.na(x)) "NS"
    else if (x < 0.001) sprintf("p=%.1e", x)
    else sprintf("p=%.3f", x)
  })

  ## Pie
  df$pos <- cumsum(df$pct) - df$pct / 2  # midpoint for label placement

  p <- ggplot(df, aes(x = "", y = pct, fill = Variable)) +
    geom_bar(stat = "identity", width = 1, colour = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    geom_text(aes(y = pos, label = plabel), size = 3.5,
              colour = "black", fontface = "bold", family = FONT_FAM) +
    scale_fill_manual(values = palette) +
    labs(title = paste0("adjusted: ", cohort_label),
         fill = NULL) +
    theme_void(base_size = 12, base_family = FONT_FAM) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 10),
          legend.position = "right")

  save_plot(p, paste0(file_prefix, "_mv_pie"), w = 4.5, h = 3.5)
  invisible(p)
}

fmt_p <- function(p) {
  sapply(p, function(x) {
    if (is.null(x) || is.na(x)) "NA" else format(signif(x, 3), scientific = TRUE)
  })
}

fmt_hr <- function(R, s) {
  r <- R[[s]]
  if (is.null(r) || is.na(r$HR)) {
    "—"
  } else {
    sprintf("%.3f (%.3f-%.3f)", r$HR, r$CI_lo, r$CI_hi)
  }
}

safe_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("\\+", "pos", x)
  x <- gsub("[^A-Za-z0-9_]", "", x)
  x
}

################################################################################
##
##  TASK 1: MACROPHAGE SPATIAL POLARITY INDEX (MSPI)
##
################################################################################

cat("\n", strrep("═", 70), "\n")
cat("TASK 1: MACROPHAGE SPATIAL POLARITY INDEX (MSPI)\n")
cat(strrep("═", 70), "\n\n")
cat("TASK 1 — MSPI\n", file = mlog, append = TRUE)

run_task1 <- function(dat, time_col, event_col,
                      m1_dens_col, m2_dens_col, cd3_dens_col,
                      cohort_label, endpoint_label, file_prefix,
                      n_boot = 200, seed = 12345) {
  
  cat(sprintf("\n─── TASK 1 | %-22s | %s ───\n", cohort_label, endpoint_label))
  
  d <- dat[, c(time_col, event_col, "M1_score", "M2_score",
               m1_dens_col, m2_dens_col, cd3_dens_col), drop = FALSE]
  names(d) <- c("TIME", "EVENT", "M1s", "M2s", "dm1", "dm2", "dcd3")
  d <- d[!is.na(d$M1s) & !is.na(d$M2s) & !is.na(d$TIME) & d$TIME > 0, ]
  d$EVENT <- as.integer(d$EVENT)
  d$MSPI  <- d$M2s - d$M1s
  
  n       <- nrow(d)
  ev      <- sum(d$EVENT)
  mspi_sd <- sd(d$MSPI, na.rm = TRUE)
  mspi_md <- median(d$MSPI, na.rm = TRUE)
  
  R <- list(
    n = n, ev = ev, mspi_sd = mspi_sd, mspi_med = mspi_md,
    cohort = cohort_label, ep = endpoint_label, skipped = FALSE
  )
  
  cat(sprintf(
    "  n=%d | events=%d (%.1f%%) | MSPI median=%.4f SD=%.4f\n",
    n, ev, 100 * ev / n, mspi_md, mspi_sd
  ))
  
  if (n < 10 || ev < 5 || is.na(mspi_sd) || mspi_sd == 0) {
    cat("  SKIPPED: insufficient sample size, events, or MSPI variance\n")
    R$skipped <- TRUE
    return(invisible(R))
  }
  
  ## ── 1a. Continuous Cox ───────────────────────────────────────────────
  d$MSPI_sd <- d$MSPI / mspi_sd
  cx1 <- coxph(Surv(TIME, EVENT) ~ MSPI_sd, data = d)
  s1  <- summary(cx1)
  ph1 <- cox.zph(cx1)
  
  R$cont <- list(
    HR = s1$conf.int[1, 1],
    CI_lo = s1$conf.int[1, 3],
    CI_hi = s1$conf.int[1, 4],
    p = s1$coefficients[1, 5],
    PH_p = ph1$table["GLOBAL", "p"]
  )
  
  cat(sprintf(
    "  1a Continuous: HR=%.3f (%.3f-%.3f)  p=%s  PH p=%.3f\n",
    R$cont$HR, R$cont$CI_lo, R$cont$CI_hi, fmt_p(R$cont$p), R$cont$PH_p
  ))
  
  ## ── 1a-extra. Quartilised MSPI analysis (Q1-Q4) ─────────────────────
  qbreaks <- unique(quantile(d$MSPI, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 2))
  
  if (length(qbreaks) >= 5) {
    d$MSPI_Q <- cut(
      d$MSPI,
      breaks = qbreaks,
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  } else {
    d$MSPI_Q <- suppressWarnings(dplyr::ntile(d$MSPI, 4))
    d$MSPI_Q <- factor(paste0("Q", d$MSPI_Q), levels = c("Q1", "Q2", "Q3", "Q4"))
  }
  
  d$MSPI_Q <- factor(d$MSPI_Q, levels = c("Q1", "Q2", "Q3", "Q4"))
  
  quart_summary <- d %>%
    dplyr::group_by(MSPI_Q) %>%
    dplyr::summarise(
      n = dplyr::n(),
      Events = sum(EVENT, na.rm = TRUE),
      Event_pct = round(100 * mean(EVENT, na.rm = TRUE), 1),
      MSPI_min = round(min(MSPI, na.rm = TRUE), 4),
      MSPI_median = round(median(MSPI, na.rm = TRUE), 4),
      MSPI_max = round(max(MSPI, na.rm = TRUE), 4),
      .groups = "drop"
    )
  
  write.csv(
    quart_summary,
    file.path(tables_dir, paste0(file_prefix, "_T1_MSPI_quartiles_summary_", output_date, ".csv")),
    row.names = FALSE
  )

  save_table_docx(
    quart_summary,
    file.path(tables_dir, paste0(file_prefix, "_T1_MSPI_quartiles_summary_", output_date, ".docx")),
    title = "MSPI Quartiles Summary"
  )

  if (length(unique(na.omit(d$MSPI_Q))) >= 2 && sum(d$EVENT) >= 5) {
    fit_q <- survfit(Surv(TIME, EVENT) ~ MSPI_Q, data = d)
    km_q  <- ggsurvplot(
      fit_q, data = d, pval = TRUE, pval.method = TRUE, conf.int = FALSE,
      risk.table = TRUE, xlab = paste0(endpoint_label, " (months)"),
      title = paste0("TASK 1 -- MSPI quartiles | ", cohort_label),
      subtitle = "Quartilised MSPI: Q1-Q4",
      legend.labs = paste0(levels(d$MSPI_Q), " (n=", as.integer(table(d$MSPI_Q)[levels(d$MSPI_Q)]), ")"),
      palette = c(COL_GREEN, "#91bfdb", "#e6c87a", COL_RED),
      ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM),
      risk.table.height = 0.25
    )
    save_km(km_q, paste0(file_prefix, "_T1_MSPI_quartiles_KM"))
  }
  
  quart_cox <- NULL
  if (length(unique(na.omit(d$MSPI_Q))) >= 2 && sum(d$EVENT) >= 5) {
    d$MSPI_Q_ref <- relevel(d$MSPI_Q, ref = "Q1")
    cx_q <- coxph(Surv(TIME, EVENT) ~ MSPI_Q_ref, data = d)
    s_q  <- summary(cx_q)
    
    quart_cox <- data.frame(
      Comparison = paste0(rownames(s_q$conf.int), " vs Q1"),
      HR = round(s_q$conf.int[, 1], 3),
      CI_lo = round(s_q$conf.int[, 3], 3),
      CI_hi = round(s_q$conf.int[, 4], 3),
      p = fmt_p(s_q$coefficients[, 5])
    )
    
    write.csv(
      quart_cox,
      file.path(tables_dir, paste0(file_prefix, "_T1_MSPI_quartiles_cox_", output_date, ".csv")),
      row.names = FALSE
    )

    save_table_docx(
      quart_cox,
      file.path(tables_dir, paste0(file_prefix, "_T1_MSPI_quartiles_cox_", output_date, ".docx")),
      title = "MSPI Quartiles Cox Regression"
    )
  }

  ## ── 1b. Median split ─────────────────────────────────────────────────
  d$grp_med <- factor(
    ifelse(d$MSPI >= mspi_md, "High_MSPI", "Low_MSPI"),
    levels = c("Low_MSPI", "High_MSPI")
  )
  cx2 <- coxph(Surv(TIME, EVENT) ~ grp_med, data = d)
  s2  <- summary(cx2)
  ph2 <- cox.zph(cx2)
  
  nh <- sum(d$grp_med == "High_MSPI")
  eh <- sum(d$EVENT[d$grp_med == "High_MSPI"])
  nl <- sum(d$grp_med == "Low_MSPI")
  el <- sum(d$EVENT[d$grp_med == "Low_MSPI"])
  
  R$med <- list(
    HR = s2$conf.int[1, 1],
    CI_lo = s2$conf.int[1, 3],
    CI_hi = s2$conf.int[1, 4],
    p = s2$coefficients[1, 5],
    PH_p = ph2$table["GLOBAL", "p"],
    cut = mspi_md, n_hi = nh, ev_hi = eh, n_lo = nl, ev_lo = el
  )
  
  fit_m <- survfit(Surv(TIME, EVENT) ~ grp_med, data = d)
  km_m  <- ggsurvplot(
    fit_m, data = d, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
    risk.table = TRUE, xlab = paste0(endpoint_label, " (months)"),
    title = paste0("TASK 1 -- MSPI Median Split | ", cohort_label),
    subtitle = paste0("MSPI = M2_score - M1_score | cut=", round(mspi_md, 4)),
    legend.labs = c(
      paste0("Low MSPI (n=", nl, ")"),
      paste0("High MSPI (n=", nh, ")")
    ),
    palette = PAL_2G,
    ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM), risk.table.height = 0.25
  )
  save_km(km_m, paste0(file_prefix, "_T1_MSPI_median"))
  
  ## ── 1c. Optimal cutpoint + bootstrap ─────────────────────────────────
  d_opt <- data.frame(TIME = d$TIME, EVENT = d$EVENT, MSPI = d$MSPI)
  
  cp <- tryCatch(
    surv_cutpoint(d_opt, time = "TIME", event = "EVENT", variables = "MSPI"),
    error = function(e) {
      cat("  surv_cutpoint failed:", e$message, "\n")
      NULL
    }
  )
  
  if (!is.null(cp)) {
    opt_cut <- cp$cutpoint$cutpoint
    d_cat   <- surv_categorize(cp, labels = c("Low_MSPI", "High_MSPI"))
    d_cat$grp <- factor(d_cat$MSPI, levels = c("Low_MSPI", "High_MSPI"))
    
    cx3 <- coxph(Surv(TIME, EVENT) ~ grp, data = d_cat)
    s3  <- summary(cx3)
    ph3 <- cox.zph(cx3)
    
    nh3 <- sum(d_cat$grp == "High_MSPI")
    eh3 <- sum(d_cat$EVENT[d_cat$grp == "High_MSPI"])
    nl3 <- sum(d_cat$grp == "Low_MSPI")
    el3 <- sum(d_cat$EVENT[d_cat$grp == "Low_MSPI"])
    
    R$opt <- list(
      HR = s3$conf.int[1, 1],
      CI_lo = s3$conf.int[1, 3],
      CI_hi = s3$conf.int[1, 4],
      p = s3$coefficients[1, 5],
      PH_p = ph3$table["GLOBAL", "p"],
      cut = opt_cut, n_hi = nh3, ev_hi = eh3, n_lo = nl3, ev_lo = el3,
      unstable = (nh3 < 10 | eh3 < 3)
    )
    
    fit_o <- survfit(Surv(TIME, EVENT) ~ grp, data = d_cat)
    km_o  <- ggsurvplot(
      fit_o, data = d_cat, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
      risk.table = TRUE, xlab = paste0(endpoint_label, " (months)"),
      title = paste0("TASK 1 -- MSPI Optimal Cutpoint | ", cohort_label),
      subtitle = paste0("cut=", round(opt_cut, 4)),
      legend.labs = c(
        paste0("Low MSPI (n=", nl3, ")"),
        paste0("High MSPI (n=", nh3, ")")
      ),
      palette = PAL_2G,
      ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM), risk.table.height = 0.25
    )
    save_km(km_o, paste0(file_prefix, "_T1_MSPI_optimal"))
    
    set.seed(seed)
    find_cut_b <- function(dat, idx) {
      db <- dat[idx, ]
      if (sum(db$EVENT) < 3) return(NA_real_)
      tryCatch({
        cs <- unique(quantile(db$MSPI, probs = seq(0.10, 0.90, 0.05), na.rm = TRUE))
        bp <- Inf
        bc <- NA_real_
        for (cv in cs) {
          g <- as.integer(db$MSPI >= cv)
          if (length(unique(g)) < 2) next
          ft <- survdiff(Surv(TIME, EVENT) ~ g, data = db)
          pv <- 1 - pchisq(ft$chisq, df = 1)
          if (!is.na(pv) && pv < bp) {
            bp <- pv
            bc <- cv
          }
        }
        bc
      }, error = function(e) NA_real_)
    }
    
    bo    <- boot(data = d_opt, statistic = function(dat, i) find_cut_b(dat, i), R = n_boot)
    bcuts <- bo$t[!is.na(bo$t)]
    
    if (length(bcuts) >= 10) {
      R$boot <- list(
        orig = opt_cut,
        mean = mean(bcuts),
        CI_lo = quantile(bcuts, 0.025),
        CI_hi = quantile(bcuts, 0.975),
        SD = sd(bcuts),
        n_ok = length(bcuts)
      )
      
      p_boot <- ggplot(data.frame(cut = bcuts), aes(x = cut)) +
        geom_histogram(bins = 25, fill = COL_BLUE, alpha = 0.75, colour = "white") +
        geom_vline(xintercept = opt_cut, colour = "red", linewidth = 1.2, linetype = "dashed") +
        geom_vline(xintercept = mean(bcuts), colour = "darkgreen", linewidth = 1) +
        labs(
          x = "Bootstrap MSPI cutpoint",
          y = "Frequency",
          title = paste0("TASK 1 -- Bootstrap Stability | ", cohort_label),
          subtitle = paste0(length(bcuts), " valid iters | ", endpoint_label),
          caption = paste0("Red=original(", round(opt_cut, 4), ")  Green=mean(", round(mean(bcuts), 4), ")")
        ) +
        theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM)
      
      save_plot(p_boot, paste0(file_prefix, "_T1_bootstrap"), w = 6, h = 4)
      
      boot_df <- data.frame(
        Metric = c("Original", "Boot_mean", "CI_lo_2.5", "CI_hi_97.5", "SD", "n_valid"),
        Value = round(c(
          opt_cut, mean(bcuts), R$boot$CI_lo, R$boot$CI_hi,
          R$boot$SD, length(bcuts)
        ), 4)
      )
      write.csv(
        boot_df,
        file.path(validation_dir, paste0(file_prefix, "_T1_bootstrap_", output_date, ".csv")),
        row.names = FALSE
      )

      save_table_docx(
        boot_df,
        file.path(validation_dir, paste0(file_prefix, "_T1_bootstrap_", output_date, ".docx")),
        title = "Bootstrap Validation"
      )
    } else {
      R$boot <- list(orig = opt_cut, mean = NA, CI_lo = NA, CI_hi = NA, SD = NA, n_ok = length(bcuts))
    }
  } else {
    R$opt <- list(HR = NA, CI_lo = NA, CI_hi = NA, p = NA, PH_p = NA, cut = NA, unstable = TRUE)
    R$boot <- list(orig = NA, mean = NA, CI_lo = NA, CI_hi = NA, SD = NA, n_ok = NA)
  }
  
  ## ── 1d. Multivariate Cox (density-adjusted, raw densities — no log) ────
  d_mv <- d
  d_mvc <- d_mv[complete.cases(d_mv[, c("MSPI_sd", "dm1", "dm2", "dcd3")]), ]

  if (nrow(d_mvc) > 10 && sum(d_mvc$EVENT) >= 5) {
    cx4 <- coxph(Surv(TIME, EVENT) ~ MSPI_sd + dm1 + dm2 + dcd3, data = d_mvc)
    s4  <- summary(cx4)
    ph4 <- cox.zph(cx4)
    
    R$mv <- list(
      HR = s4$conf.int["MSPI_sd", 1],
      CI_lo = s4$conf.int["MSPI_sd", 3],
      CI_hi = s4$conf.int["MSPI_sd", 4],
      p = s4$coefficients["MSPI_sd", 5],
      PH_p = ph4$table["GLOBAL", "p"],
      n = nrow(d_mvc), ev = sum(d_mvc$EVENT)
    )
  } else {
    R$mv <- list(HR = NA, CI_lo = NA, CI_hi = NA, p = NA, PH_p = NA, n = NA, ev = NA)
  }
  
  task1_tab <- data.frame(
    Cohort = cohort_label,
    Endpoint = endpoint_label,
    Approach = c("Continuous(per SD)", "Median split", "Optimal cutpoint", "MV Continuous(adj)"),
    n = c(n, n, n, R$mv$n),
    Events = c(ev, ev, ev, R$mv$ev),
    Cut = c(NA, round(mspi_md, 4), round(R$opt$cut, 4), NA),
    HR = round(c(R$cont$HR, R$med$HR, R$opt$HR, R$mv$HR), 3),
    CI_lo = round(c(R$cont$CI_lo, R$med$CI_lo, R$opt$CI_lo, R$mv$CI_lo), 3),
    CI_hi = round(c(R$cont$CI_hi, R$med$CI_hi, R$opt$CI_hi, R$mv$CI_hi), 3),
    p = c(fmt_p(R$cont$p), fmt_p(R$med$p), fmt_p(R$opt$p), fmt_p(R$mv$p)),
    PH_p = c(fmt_p(R$cont$PH_p), fmt_p(R$med$PH_p), fmt_p(R$opt$PH_p), fmt_p(R$mv$PH_p))
  )
  
  write.csv(
    task1_tab,
    file.path(tables_dir, paste0(file_prefix, "_T1_", output_date, ".csv")),
    row.names = FALSE
  )

  save_table_docx(
    task1_tab,
    file.path(tables_dir, paste0(file_prefix, "_T1_", output_date, ".docx")),
    title = "Task 1 - MSPI Survival Analysis"
  )

  R$quartiles <- quart_summary
  R$quart_cox <- quart_cox

  ## v4: Store Cox model objects for mini-forest and pie chart generation
  R$model_uv <- if (exists("cx1")) cx1 else NULL
  R$model_mv <- if (exists("cx4")) cx4 else NULL
  R$data_mv  <- if (exists("d_mvc")) d_mvc else NULL

  cat(sprintf(
    "  → Cont HR=%s p=%s | MV HR=%s p=%s\n",
    fmt_hr(R, "cont"), fmt_p(R$cont$p), fmt_hr(R, "mv"), fmt_p(R$mv$p)
  ), file = mlog, append = TRUE)

  invisible(R)
}




################################################################################
##
##  TASK 2: 4-QUADRANT SPATIAL IMMUNE PHENOTYPING
##
################################################################################

cat("\n", strrep("═", 70), "\n")
cat("TASK 2: 4-QUADRANT SPATIAL IMMUNE PHENOTYPING\n")
cat(strrep("═", 70), "\n\n")
cat("TASK 2 — 4-Quadrant\n", file = mlog, append = TRUE)

run_task2 <- function(dat, time_col, event_col,
                      cohort_label, endpoint_label, file_prefix) {
  
  cat(sprintf("\n─── TASK 2 | %-22s | %s ───\n", cohort_label, endpoint_label))
  
  d <- dat[!is.na(dat$M1_score) & !is.na(dat$M2_score) &
             !is.na(dat[[time_col]]) & dat[[time_col]] > 0, ]
  d$TIME  <- d[[time_col]]
  d$EVENT <- as.integer(d[[event_col]])
  n <- nrow(d)
  ev <- sum(d$EVENT)
  
  if (n < 10 || ev < 5) {
    cat("  SKIPPED: insufficient data for quadrant analysis\n")
    return(invisible(NULL))
  }
  
  m1c <- median(d$M1_score, na.rm = TRUE)
  m2c <- median(d$M2_score, na.rm = TRUE)
  
  d$Quad <- factor(
    ifelse(d$M1_score >= m1c & d$M2_score <  m2c, "Q1_Best",
           ifelse(d$M1_score >= m1c & d$M2_score >= m2c, "Q2_Mixed",
                  ifelse(d$M1_score <  m1c & d$M2_score <  m2c, "Q3_Interm",
                         "Q4_Worst"))),
    levels = c("Q1_Best", "Q2_Mixed", "Q3_Interm", "Q4_Worst")
  )
  
  qt  <- table(d$Quad)
  evt <- tapply(d$EVENT, d$Quad, sum)
  
  summ <- data.frame(
    Quadrant = names(qt),
    Definition = c("M1-close/M2-far", "M1-close/M2-close", "M1-far/M2-far", "M1-far/M2-close"),
    Phenotype = c("Activation dominant", "Competing signals", "Immune disengagement", "Suppression dominant"),
    n = as.integer(qt),
    Events = as.integer(evt),
    Event_pct = round(100 * evt / qt, 1)
  )
  
  write.csv(
    summ,
    file.path(tables_dir, paste0(file_prefix, "_T2_summary_", output_date, ".csv")),
    row.names = FALSE
  )

  save_table_docx(
    summ,
    file.path(tables_dir, paste0(file_prefix, "_T2_summary_", output_date, ".docx")),
    title = "Task 2 - Quadrant Summary"
  )

  ## 2a. 4-group KM
  fit4 <- survfit(Surv(TIME, EVENT) ~ Quad, data = d)
  km4  <- ggsurvplot(
    fit4, data = d, pval = TRUE, pval.method = TRUE, conf.int = FALSE,
    risk.table = TRUE, xlab = paste0(endpoint_label, " (months)"),
    title = paste0("TASK 2 -- 4-Quadrant Survival | ", cohort_label),
    subtitle = paste0("M1 cut=", round(m1c, 4), "  M2 cut=", round(m2c, 4)),
    legend.labs = c(
      "Q1 Best: M1-close/M2-far",
      "Q2 Mixed: both-close",
      "Q3 Intermediate: both-far",
      "Q4 Worst: M2-close/M1-far"
    ),
    palette = unname(PAL_4Q),
    ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM), risk.table.height = 0.30
  )
  save_km(km4, paste0(file_prefix, "_T2_4quad_KM"))
  
  ## 2b. Pairwise HRs vs Q1
  d$Q_ref <- relevel(d$Quad, ref = "Q1_Best")
  cx_all <- coxph(Surv(TIME, EVENT) ~ Q_ref, data = d)
  s_all  <- summary(cx_all)
  
  pw_tab <- data.frame(
    Comparison = paste0(rownames(s_all$conf.int), " vs Q1_Best"),
    HR = round(s_all$conf.int[, 1], 3),
    CI_lo = round(s_all$conf.int[, 3], 3),
    CI_hi = round(s_all$conf.int[, 4], 3),
    p = fmt_p(s_all$coefficients[, 5])
  )
  
  write.csv(
    pw_tab,
    file.path(tables_dir, paste0(file_prefix, "_T2_pairwise_", output_date, ".csv")),
    row.names = FALSE
  )

  save_table_docx(
    pw_tab,
    file.path(tables_dir, paste0(file_prefix, "_T2_pairwise_", output_date, ".docx")),
    title = "Task 2 - Pairwise HRs vs Q1"
  )

  ## Dedicated Q4 metric analyses
  d$Q4_vs_Q1 <- ifelse(d$Quad %in% c("Q1_Best", "Q4_Worst"), as.character(d$Quad), NA)
  d_q41 <- d[!is.na(d$Q4_vs_Q1), ]
  q4_q1_res <- NULL
  
  if (nrow(d_q41) >= 10 && length(unique(d_q41$Q4_vs_Q1)) == 2 && sum(d_q41$EVENT) >= 5) {
    d_q41$Q4_vs_Q1 <- factor(d_q41$Q4_vs_Q1, levels = c("Q1_Best", "Q4_Worst"))
    fit_q41 <- survfit(Surv(TIME, EVENT) ~ Q4_vs_Q1, data = d_q41)
    km_q41  <- ggsurvplot(
      fit_q41, data = d_q41, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
      risk.table = TRUE, xlab = paste0(endpoint_label, " (months)"),
      title = paste0("TASK 2 -- Q4 vs Q1 KM | ", cohort_label),
      subtitle = "Primary worst-vs-best quadrant contrast",
      legend.labs = c(
        paste0("Q1 Best (n=", sum(d_q41$Q4_vs_Q1 == "Q1_Best"), ")"),
        paste0("Q4 Worst (n=", sum(d_q41$Q4_vs_Q1 == "Q4_Worst"), ")")
      ),
      palette = PAL_2G,
      ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM), risk.table.height = 0.25
    )
    save_km(km_q41, paste0(file_prefix, "_T2_Q4_vs_Q1_KM"))
    
    cx_q41 <- coxph(Surv(TIME, EVENT) ~ Q4_vs_Q1, data = d_q41)
    s_q41  <- summary(cx_q41)
    q4_q1_res <- data.frame(
      Contrast = "Q4_Worst vs Q1_Best",
      n = nrow(d_q41),
      Events = sum(d_q41$EVENT),
      HR = round(s_q41$conf.int[1, 1], 3),
      CI_lo = round(s_q41$conf.int[1, 3], 3),
      CI_hi = round(s_q41$conf.int[1, 4], 3),
      p = fmt_p(s_q41$coefficients[1, 5])
    )
    write.csv(
      q4_q1_res,
      file.path(tables_dir, paste0(file_prefix, "_T2_Q4_vs_Q1_", output_date, ".csv")),
      row.names = FALSE
    )

    save_table_docx(
      q4_q1_res,
      file.path(tables_dir, paste0(file_prefix, "_T2_Q4_vs_Q1_", output_date, ".docx")),
      title = "Task 2 - Q4 vs Q1"
    )
  }

  d$Q4_flag <- factor(ifelse(d$Quad == "Q4_Worst", "Q4_Worst", "Non_Q4"),
                      levels = c("Non_Q4", "Q4_Worst"))
  q4_rest_res <- NULL
  
  if (length(unique(d$Q4_flag)) == 2 && sum(d$EVENT) >= 5) {
    fit_q4rest <- survfit(Surv(TIME, EVENT) ~ Q4_flag, data = d)
    km_q4rest  <- ggsurvplot(
      fit_q4rest, data = d, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
      risk.table = TRUE, xlab = paste0(endpoint_label, " (months)"),
      title = paste0("TASK 2 -- Q4 vs non-Q4 KM | ", cohort_label),
      subtitle = "Q4 worst quadrant versus all other patients",
      legend.labs = c(
        paste0("Non-Q4 (n=", sum(d$Q4_flag == "Non_Q4"), ")"),
        paste0("Q4 Worst (n=", sum(d$Q4_flag == "Q4_Worst"), ")")
      ),
      palette = PAL_Q4_REST,
      ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM), risk.table.height = 0.25
    )
    save_km(km_q4rest, paste0(file_prefix, "_T2_Q4_vs_nonQ4_KM"))
    
    cx_q4rest <- coxph(Surv(TIME, EVENT) ~ Q4_flag, data = d)
    s_q4rest  <- summary(cx_q4rest)
    
    q4_rest_res <- data.frame(
      Contrast = "Q4_Worst vs Non_Q4",
      n = nrow(d),
      Events = sum(d$EVENT),
      HR = round(s_q4rest$conf.int[1, 1], 3),
      CI_lo = round(s_q4rest$conf.int[1, 3], 3),
      CI_hi = round(s_q4rest$conf.int[1, 4], 3),
      p = fmt_p(s_q4rest$coefficients[1, 5])
    )
    write.csv(
      q4_rest_res,
      file.path(tables_dir, paste0(file_prefix, "_T2_Q4_vs_nonQ4_", output_date, ".csv")),
      row.names = FALSE
    )

    save_table_docx(
      q4_rest_res,
      file.path(tables_dir, paste0(file_prefix, "_T2_Q4_vs_nonQ4_", output_date, ".docx")),
      title = "Task 2 - Q4 vs Non-Q4"
    )
  }

  ## 2d. Q1+Q2 vs Q3+Q4 (M1-high vs M1-low based on quad groups)
  d$Q12_flag <- factor(
    ifelse(d$Quad %in% c("Q1_Best", "Q2_Mixed"), "Q12_M1hi", "Q34_M1lo"),
    levels = c("Q12_M1hi", "Q34_M1lo")
  )
  q12_q34_res <- NULL

  if (length(unique(d$Q12_flag)) == 2 && sum(d$EVENT) >= 5) {
    fit_q12q34 <- survfit(Surv(TIME, EVENT) ~ Q12_flag, data = d)
    km_q12q34  <- ggsurvplot(
      fit_q12q34, data = d, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
      risk.table = TRUE, xlab = paste0(endpoint_label, " (months)"),
      title = paste0("TASK 2 -- Q1+Q2 vs Q3+Q4 KM | ", cohort_label),
      subtitle = "M1-high affinity (Q1+Q2) versus M1-low affinity (Q3+Q4)",
      legend.labs = c(
        paste0("Q1+Q2 M1-hi (n=", sum(d$Q12_flag == "Q12_M1hi"), ")"),
        paste0("Q3+Q4 M1-lo (n=", sum(d$Q12_flag == "Q34_M1lo"), ")")
      ),
      palette = PAL_Q12_Q34,
      ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM), risk.table.height = 0.25
    )
    save_km(km_q12q34, paste0(file_prefix, "_T2_Q12_vs_Q34_KM"))

    cx_q12q34 <- coxph(Surv(TIME, EVENT) ~ Q12_flag, data = d)
    s_q12q34  <- summary(cx_q12q34)

    q12_q34_res <- data.frame(
      Contrast = "Q34_M1lo vs Q12_M1hi",
      n = nrow(d),
      Events = sum(d$EVENT),
      HR = round(s_q12q34$conf.int[1, 1], 3),
      CI_lo = round(s_q12q34$conf.int[1, 3], 3),
      CI_hi = round(s_q12q34$conf.int[1, 4], 3),
      p = fmt_p(s_q12q34$coefficients[1, 5])
    )
    write.csv(
      q12_q34_res,
      file.path(tables_dir, paste0(file_prefix, "_T2_Q12_vs_Q34_", output_date, ".csv")),
      row.names = FALSE
    )

    save_table_docx(
      q12_q34_res,
      file.path(tables_dir, paste0(file_prefix, "_T2_Q12_vs_Q34_", output_date, ".docx")),
      title = "Task 2 - Q1+Q2 vs Q3+Q4"
    )
  }

  ## Scatter
  p_sc <- ggplot(d, aes(x = M1_score, y = M2_score, colour = Quad, shape = factor(EVENT))) +
    geom_point(alpha = 0.55, size = 2.5) +
    geom_vline(xintercept = m1c, linetype = "dashed", linewidth = 0.6, colour = "grey40") +
    geom_hline(yintercept = m2c, linetype = "dashed", linewidth = 0.6, colour = "grey40") +
    scale_colour_manual(
      values = PAL_4Q,
      labels = c("Q1 Best", "Q2 Mixed", "Q3 Intermediate", "Q4 Worst")
    ) +
    scale_shape_manual(values = c("0" = 1, "1" = 16), labels = c("Censored", "Event")) +
    labs(
      x = "M1 affinity score",
      y = "M2 affinity score",
      title = paste0("TASK 2 -- Spatial Immune Landscape | ", cohort_label),
      subtitle = "Dashed lines = cohort medians | Filled circles = events",
      colour = "Quadrant", shape = "Status"
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM)
  save_plot(p_sc, paste0(file_prefix, "_T2_scatter"), w = 6, h = 5)
  
  invisible(list(
    m1c = m1c, m2c = m2c, pairwise = pw_tab,
    q4_q1 = q4_q1_res, q4_rest = q4_rest_res, q12_q34 = q12_q34_res
  ))
}

################################################################################
##
##  TASK 3: M1 × M2 INTERACTION IN JOINT COX MODEL
##
################################################################################

cat("\n", strrep("═", 70), "\n")
cat("TASK 3: M1 × M2 INTERACTION IN JOINT COX MODEL\n")
cat(strrep("═", 70), "\n\n")
cat("TASK 3 — M1×M2 interaction\n", file = mlog, append = TRUE)

run_task3 <- function(dat, time_col, event_col,
                      cohort_label, endpoint_label, file_prefix) {
  
  cat(sprintf("\n─── TASK 3 | %-22s | %s ───\n", cohort_label, endpoint_label))
  
  d <- dat[!is.na(dat$M1_score) & !is.na(dat$M2_score) &
             !is.na(dat[[time_col]]) & dat[[time_col]] > 0, ]
  d$TIME  <- d[[time_col]]
  d$EVENT <- as.integer(d[[event_col]])
  n <- nrow(d)
  ev <- sum(d$EVENT)
  
  if (n < 15 || ev < 10) {
    cat("  SKIPPED: insufficient sample size or events\n")
    return(invisible(NULL))
  }
  
  d$M1_sd <- d$M1_score / sd(d$M1_score, na.rm = TRUE)
  d$M2_sd <- d$M2_score / sd(d$M2_score, na.rm = TRUE)
  
  epp <- ev / 3
  
  cx_add <- coxph(Surv(TIME, EVENT) ~ M1_sd + M2_sd, data = d)
  s_add  <- summary(cx_add)
  
  cx_int <- coxph(Surv(TIME, EVENT) ~ M1_sd * M2_sd, data = d)
  s_int  <- summary(cx_int)
  ph_int <- cox.zph(cx_int)
  
  lrt   <- anova(cx_add, cx_int)
  lrt_pcol <- grep("Pr|P\\(", colnames(lrt), value = TRUE)
  lrt_p <- if (length(lrt_pcol) > 0) suppressWarnings(as.numeric(lrt[2, lrt_pcol[1]])) else NA_real_
  
  m1_hi <- quantile(d$M1_sd, 0.75)
  m1_lo <- quantile(d$M1_sd, 0.25)
  m2_hi <- quantile(d$M2_sd, 0.75)
  m2_lo <- quantile(d$M2_sd, 0.25)
  
  pg <- data.frame(
    Scenario = c("Q1 Best\nM1-hi/M2-lo", "Q2 Mixed\nM1-hi/M2-hi",
                 "Q3 Interm\nM1-lo/M2-lo", "Q4 Worst\nM1-lo/M2-hi"),
    M1_sd = c(m1_hi, m1_hi, m1_lo, m1_lo),
    M2_sd = c(m2_lo, m2_hi, m2_lo, m2_hi),
    Fill = c(COL_GREEN, COL_AMBER, COL_GREY, COL_RED)
  )
  pg$HR <- exp(predict(cx_int, newdata = pg, type = "lp"))
  pg$Scenario <- factor(pg$Scenario, levels = pg$Scenario)
  
  p_bar <- ggplot(pg, aes(x = Scenario, y = HR, fill = Scenario)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_text(aes(label = round(HR, 2)), vjust = -0.4, fontface = "bold", size = 4.5, family = FONT_FAM) +
    scale_fill_manual(values = setNames(pg$Fill, pg$Scenario)) +
    labs(
      x = NULL, y = "Predicted HR (vs cohort mean)",
      title = paste0("TASK 3 -- Predicted HR by Spatial Phenotype | ", cohort_label),
      subtitle = paste0(
        "LRT p=", round(lrt_p, 4),
        "  Interaction HR=", round(s_int$conf.int["M1_sd:M2_sd", 1], 3),
        "  p=", fmt_p(s_int$coefficients["M1_sd:M2_sd", 5])
      ),
      caption = "HR at 25th and 75th percentiles of each score from interaction Cox model"
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
    theme(legend.position = "none", axis.text.x = element_text(size = 12))
  save_plot(p_bar, paste0(file_prefix, "_T3_interaction_HR"), w = 6, h = 4.5)
  
  d$M1_tert <- cut(
    d$M1_score,
    breaks = quantile(d$M1_score, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("M1 T1 (high affinity)", "M1 T2 (mid)", "M1 T3 (low affinity)"),
    include.lowest = TRUE
  )
  
  p_loe <- ggplot(
    d[!is.na(d$M1_tert), ],
    aes(x = M2_score, y = EVENT, colour = M1_tert, group = M1_tert)
  ) +
    stat_smooth(method = "loess", se = TRUE, linewidth = 1.1) +
    geom_vline(xintercept = median(d$M2_score, na.rm = TRUE), linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c(
      "M1 T1 (high affinity)" = COL_GREEN,
      "M1 T2 (mid)" = COL_AMBER,
      "M1 T3 (low affinity)" = COL_RED
    )) +
    labs(
      x = "M2 affinity score",
      y = "Smoothed event probability (LOESS)",
      title = paste0("TASK 3 -- M2 Effect by M1 Tertile | ", cohort_label),
      subtitle = "Diverging lines suggest M1 modifies the M2-outcome relationship",
      colour = "M1 tertile"
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM)
  save_plot(p_loe, paste0(file_prefix, "_T3_LOESS"), w = 6, h = 4.5)
  
  res_tab <- data.frame(
    Cohort = cohort_label, Endpoint = endpoint_label, EPP = round(epp, 1),
    Model = c("Additive-M1", "Additive-M2", "Int-M1", "Int-M2", "Int-M1:M2"),
    HR = round(c(
      s_add$conf.int[1, 1], s_add$conf.int[2, 1],
      s_int$conf.int[1, 1], s_int$conf.int[2, 1], s_int$conf.int[3, 1]
    ), 3),
    CI_lo = round(c(
      s_add$conf.int[1, 3], s_add$conf.int[2, 3],
      s_int$conf.int[1, 3], s_int$conf.int[2, 3], s_int$conf.int[3, 3]
    ), 3),
    CI_hi = round(c(
      s_add$conf.int[1, 4], s_add$conf.int[2, 4],
      s_int$conf.int[1, 4], s_int$conf.int[2, 4], s_int$conf.int[3, 4]
    ), 3),
    p = c(
      fmt_p(s_add$coefficients[1, 5]), fmt_p(s_add$coefficients[2, 5]),
      fmt_p(s_int$coefficients[1, 5]), fmt_p(s_int$coefficients[2, 5]),
      fmt_p(s_int$coefficients[3, 5])
    ),
    LRT_p = c(NA, NA, NA, NA, round(lrt_p, 4))
  )
  
  write.csv(
    res_tab,
    file.path(tables_dir, paste0(file_prefix, "_T3_", output_date, ".csv")),
    row.names = FALSE
  )

  save_table_docx(
    res_tab,
    file.path(tables_dir, paste0(file_prefix, "_T3_", output_date, ".docx")),
    title = "Task 3 - M1 x M2 Interaction"
  )

  invisible(list(
    lrt_p = lrt_p,
    int_HR = s_int$conf.int["M1_sd:M2_sd", 1],
    int_p = s_int$coefficients["M1_sd:M2_sd", 5],
    n = n, ev = ev, epp = epp
  ))
}

################################################################################
##
##  TASK 4: CROSS-COHORT POOLED INDIVIDUAL PATIENT DATA ANALYSIS
##
################################################################################

cat("\n", strrep("═", 70), "\n")
cat("TASK 4: CROSS-COHORT POOLED INDIVIDUAL PATIENT DATA ANALYSIS\n")
cat(strrep("═", 70), "\n\n")
cat("TASK 4 — Pooled IPD\n", file = mlog, append = TRUE)

prep_pool <- function(dat, time_col, event_col, label) {
  d <- dat[!is.na(dat$M1_score) & !is.na(dat$M2_score) &
             !is.na(dat[[time_col]]) & dat[[time_col]] > 0, ]
  d$TIME  <- d[[time_col]]
  d$EVENT <- as.integer(d[[event_col]])
  d$M1_sd   <- as.numeric(scale(d$M1_score))
  d$M2_sd   <- as.numeric(scale(d$M2_score))
  d$MSPI_sd <- as.numeric(scale(d$M2_score - d$M1_score))
  d$Cohort  <- label
  d[, c("TIME", "EVENT", "M1_sd", "M2_sd", "MSPI_sd", "Cohort")]
}

run_task4 <- function(pool, pool_label, file_prefix) {

  if (nrow(pool) < 20 || sum(pool$EVENT) < 10) {
    cat("  TASK 4 skipped for", pool_label, ": insufficient pooled data\n")
    return(invisible(NULL))
  }

  cx_m <- coxph(Surv(TIME, EVENT) ~ MSPI_sd + strata(Cohort), data = pool)
  s_m  <- summary(cx_m)
  ph_m <- cox.zph(cx_m)

  cx_j <- coxph(Surv(TIME, EVENT) ~ M1_sd + M2_sd + strata(Cohort), data = pool)
  s_j  <- summary(cx_j)

  cx_het <- tryCatch(
    coxph(Surv(TIME, EVENT) ~ MSPI_sd + MSPI_sd:Cohort + strata(Cohort), data = pool),
    error = function(e) NULL
  )

  if (!is.null(cx_het)) {
    lrt_h <- tryCatch(anova(cx_m, cx_het), error = function(e) NULL)
    hpcol <- if (!is.null(lrt_h)) grep("Pr|P\\(", colnames(lrt_h), value = TRUE) else character(0)
    het_p <- if (length(hpcol) > 0) suppressWarnings(as.numeric(lrt_h[2, hpcol[1]])) else NA_real_
  } else {
    het_p <- NA_real_
  }

  forest <- do.call(rbind, lapply(levels(factor(pool$Cohort)), function(coh) {
    d <- pool[pool$Cohort == coh, ]
    if (sum(d$EVENT) < 5) return(NULL)
    cx <- coxph(Surv(TIME, EVENT) ~ MSPI_sd, data = d)
    s  <- summary(cx)
    data.frame(
      Label = coh,
      HR = s$conf.int[1, 1],
      CI_lo = s$conf.int[1, 3],
      CI_hi = s$conf.int[1, 4],
      p = s$coefficients[1, 5],
      n = nrow(d),
      events = sum(d$EVENT)
    )
  }))

  forest <- rbind(
    forest,
    data.frame(
      Label = "Pooled (stratified Cox)",
      HR = s_m$conf.int[1, 1],
      CI_lo = s_m$conf.int[1, 3],
      CI_hi = s_m$conf.int[1, 4],
      p = s_m$coefficients[1, 5],
      n = nrow(pool),
      events = sum(pool$EVENT)
    )
  )

  forest$Label <- factor(forest$Label, levels = rev(forest$Label))
  is_pooled <- forest$Label == levels(forest$Label)[1]

  ## v4: Nicer forest — coloured diamonds, muted tones, narrower
  p_f <- ggplot(forest, aes(x = HR, y = Label, xmin = CI_lo, xmax = CI_hi)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
    geom_errorbarh(height = 0.22, linewidth = 0.9, colour = COL_DKGREY) +
    geom_point(aes(size = events, colour = is_pooled),
               shape = 18) +
    geom_text(
      aes(label = sprintf("HR=%.2f (%.2f-%.2f)  n=%d ev=%d  p=%s",
                          HR, CI_lo, CI_hi, n, events, fmt_p(p))),
      x = max(forest$CI_hi) + 0.05, hjust = 0, size = 3.5, family = FONT_FAM
    ) +
    scale_size_continuous(range = c(5, 9), guide = "none") +
    scale_colour_manual(values = c("FALSE" = COL_DKGREEN, "TRUE" = COL_NAVY), guide = "none") +
    labs(
      x = "HR for MSPI per within-cohort SD",
      title = paste0("Forest Plot | ", pool_label),
      subtitle = "Diamond = pooled stratified Cox | Scores standardised within each cohort"
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
    theme(panel.grid.major.y = element_blank()) +
    xlim(min(forest$CI_lo) - 0.05, max(forest$CI_hi) + 2.0)

  save_plot(p_f, paste0(file_prefix, "_T4_forest"), w = 6.5, h = 3.5)

  write.csv(
    forest,
    file.path(tables_dir, paste0(file_prefix, "_T4_forest_", output_date, ".csv")),
    row.names = FALSE
  )

  save_table_docx(
    forest,
    file.path(tables_dir, paste0(file_prefix, "_T4_forest_", output_date, ".docx")),
    title = "Task 4 - Pooled IPD Forest"
  )

  invisible(list(
    HR_mspi = s_m$conf.int[1, 1], p_mspi = s_m$coefficients[1, 5],
    HR_M1 = s_j$conf.int["M1_sd", 1], p_M1 = s_j$coefficients["M1_sd", 5],
    HR_M2 = s_j$conf.int["M2_sd", 1], p_M2 = s_j$coefficients["M2_sd", 5],
    het_p = het_p
  ))
}




################################################################################
##
##  TASK 5: INCREMENTAL PROGNOSTIC VALUE — AFFINITY vs DENSITY
##
################################################################################

cat("\n", strrep("═", 70), "\n")
cat("TASK 5: INCREMENTAL PROGNOSTIC VALUE — AFFINITY vs DENSITY\n")
cat(strrep("═", 70), "\n\n")
cat("TASK 5 — Incremental C-index\n", file = mlog, append = TRUE)

run_task5 <- function(dat, time_col, event_col,
                      m1_dens_col, m2_dens_col, cd3_dens_col,
                      has_clinical = FALSE,
                      cohort_label, endpoint_label, file_prefix,
                      n_boot = 200, seed = 12345) {
  
  d <- dat[!is.na(dat$M1_score) & !is.na(dat$M2_score) &
             !is.na(dat[[time_col]]) & dat[[time_col]] > 0, ]
  d$TIME  <- d[[time_col]]
  d$EVENT <- as.integer(d[[event_col]])
  d$M1_sd   <- as.numeric(scale(d$M1_score))
  d$M2_sd   <- as.numeric(scale(d$M2_score))
  d$MSPI_sd <- as.numeric(scale(d$M2_score - d$M1_score))
  d$dens_m1  <- d[[m1_dens_col]]
  d$dens_m2  <- d[[m2_dens_col]]
  d$dens_cd3 <- d[[cd3_dens_col]]
  
  base_v <- c("TIME", "EVENT", "dens_m1", "dens_m2", "dens_cd3", "M1_sd", "M2_sd", "MSPI_sd")
  
  if (has_clinical) {
    d$Grade_n <- suppressWarnings(as.integer(as.character(d$Grade)))
    d$Node_n  <- as.integer(d$Node_pos)
    if (!all(is.na(d$Grade_n))) base_v <- c(base_v, "Grade_n", "Node_n")
  }
  
  d_c <- d[complete.cases(d[, base_v]), ]
  n   <- nrow(d_c)
  ev  <- sum(d_c$EVENT)
  
  if (n < 20 || ev < 10) {
    cat("  TASK 5 skipped for", cohort_label, ": insufficient data\n")
    return(invisible(NULL))
  }
  
  cx_dens    <- coxph(Surv(TIME, EVENT) ~ dens_m1 + dens_m2 + dens_cd3, data = d_c, x = TRUE, y = TRUE)
  cx_m1      <- coxph(Surv(TIME, EVENT) ~ M1_sd, data = d_c, x = TRUE, y = TRUE)
  cx_m2      <- coxph(Surv(TIME, EVENT) ~ M2_sd, data = d_c, x = TRUE, y = TRUE)
  cx_mspi_o  <- coxph(Surv(TIME, EVENT) ~ MSPI_sd, data = d_c, x = TRUE, y = TRUE)
  cx_dm      <- coxph(Surv(TIME, EVENT) ~ dens_m1 + dens_m2 + dens_cd3 + MSPI_sd, data = d_c, x = TRUE, y = TRUE)
  
  get_ci <- function(cx, dat, B = n_boot, s = seed) {
    c0 <- concordance(cx)$concordance
    set.seed(s)
    boots <- replicate(B, {
      idx <- sample(seq_len(nrow(dat)), replace = TRUE)
      db  <- dat[idx, ]
      tryCatch(concordance(update(cx, data = db))$concordance, error = function(e) NA_real_)
    })
    boots <- boots[!is.na(boots)]
    list(
      c = round(c0, 4),
      lo = round(quantile(boots, 0.025), 4),
      hi = round(quantile(boots, 0.975), 4),
      n = length(boots)
    )
  }
  
  ci_d  <- get_ci(cx_dens,   d_c)
  ci_m1 <- get_ci(cx_m1,     d_c)
  ci_m2 <- get_ci(cx_m2,     d_c)
  ci_mo <- get_ci(cx_mspi_o, d_c)
  ci_dm <- get_ci(cx_dm,     d_c)
  
  models <- c("Density only", "M1 affinity", "M2 affinity", "MSPI only", "Density + MSPI")
  c_vec  <- c(ci_d$c, ci_m1$c, ci_m2$c, ci_mo$c, ci_dm$c)
  lo_v   <- c(ci_d$lo, ci_m1$lo, ci_m2$lo, ci_mo$lo, ci_dm$lo)
  hi_v   <- c(ci_d$hi, ci_m1$hi, ci_m2$hi, ci_mo$hi, ci_dm$hi)
  aic_v  <- c(AIC(cx_dens), AIC(cx_m1), AIC(cx_m2), AIC(cx_mspi_o), AIC(cx_dm))
  delta_c <- ci_dm$c - ci_d$c
  
  lrt_t5    <- anova(cx_dens, cx_dm)
  lrt_pcol5 <- grep("Pr|P\\(", colnames(lrt_t5), value = TRUE)
  lrt_p     <- if (length(lrt_pcol5) > 0) suppressWarnings(as.numeric(lrt_t5[2, lrt_pcol5[1]])) else NA_real_
  
  df_c <- data.frame(
    Model = factor(models, levels = models),
    C = c_vec,
    CI_lo = lo_v,
    CI_hi = hi_v,
    AIC = round(aic_v, 1),
    stringsAsFactors = FALSE
  )
  
  col_map <- c(
    "Density only"   = COL_GREY,
    "M1 affinity"    = COL_BLUE,
    "M2 affinity"    = COL_RED,
    "MSPI only"      = COL_NAVY,
    "Density + MSPI" = COL_DKGREEN
  )
  
  shape_map <- c(
    "Density only"   = 18,
    "M1 affinity"    = 16,
    "M2 affinity"    = 16,
    "MSPI only"      = 16,
    "Density + MSPI" = 18
  )
  
  ## Main vertical plot
  p_c <- ggplot(df_c, aes(x = Model, y = C, ymin = CI_lo, ymax = CI_hi,
                          colour = Model, shape = Model)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_errorbar(width = 0.18, linewidth = 0.9) +
    geom_point(size = 5) +
    geom_text(
      aes(label = sprintf("C=%.3f\nAIC=%.0f", C, AIC)),
      vjust = -0.7, size = 3.6, family = FONT_FAM, lineheight = 0.95
    ) +
    scale_colour_manual(values = col_map, guide = "none") +
    scale_shape_manual(values = shape_map, guide = "none") +
    labs(
      x = NULL,
      y = "C-index (Harrell) with bootstrap 95% CI",
      title = paste0("TASK 5 -- Incremental Prognostic Value | ", cohort_label),
      subtitle = paste0("LRT (density vs density+MSPI): p=", round(lrt_p, 4), "  |  ΔC=", round(delta_c, 4)),
      caption = "Bootstrap 200 iter. | AIC lower = better"
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    ylim(min(lo_v) - 0.03, max(hi_v) + 0.08)
  
  save_plot(p_c, paste0(file_prefix, "_T5_cindex"), w = 6.2, h = 5.2)
  
  ## Additional compact plot: MSPI only vs Density + MSPI
  df_c_small <- data.frame(
    Model = factor(c("MSPI only", "Density + MSPI"),
                   levels = c("MSPI only", "Density + MSPI")),
    C     = c(ci_mo$c, ci_dm$c),
    CI_lo = c(ci_mo$lo, ci_dm$lo),
    CI_hi = c(ci_mo$hi, ci_dm$hi),
    AIC   = round(c(AIC(cx_mspi_o), AIC(cx_dm)), 1),
    stringsAsFactors = FALSE
  )
  
  small_col_map <- c(
    "MSPI only"      = COL_NAVY,
    "Density + MSPI" = COL_DKGREEN
  )
  
  small_shape_map <- c(
    "MSPI only"      = 16,
    "Density + MSPI" = 18
  )
  
  p_c_small <- ggplot(df_c_small, aes(x = Model, y = C, ymin = CI_lo, ymax = CI_hi,
                                      colour = Model, shape = Model)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_errorbar(width = 0.14, linewidth = 0.9) +
    geom_point(size = 5.5) +
    geom_text(
      aes(label = sprintf("C=%.3f\nAIC=%.0f", C, AIC)),
      vjust = -0.7, size = 3.8, family = FONT_FAM, lineheight = 0.95
    ) +
    scale_colour_manual(values = small_col_map, guide = "none") +
    scale_shape_manual(values = small_shape_map, guide = "none") +
    labs(
      x = NULL,
      y = "C-index (Harrell) with bootstrap 95% CI",
      title = paste0("TASK 5 -- MSPI-focused Comparison | ", cohort_label),
      subtitle = paste0("MSPI only vs Density + MSPI  |  ΔC=", round(delta_c, 4)),
      caption = "Bootstrap 200 iter. | AIC lower = better"
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
    theme(
      axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    ylim(min(df_c_small$CI_lo) - 0.03, max(df_c_small$CI_hi) + 0.08)
  
  save_plot(p_c_small, paste0(file_prefix, "_T5_cindex_MSPIfocused"), w = 2.0, h = 5.0)

  ## v8: save MSPIfocused c-index data as CSV + docx
  write.csv(
    df_c_small,
    file.path(tables_dir, paste0(file_prefix, "_T5_cindex_MSPIfocused_", output_date, ".csv")),
    row.names = FALSE
  )
  save_table_docx(
    df_c_small,
    file.path(tables_dir, paste0(file_prefix, "_T5_cindex_MSPIfocused_", output_date, ".docx")),
    title = paste0("Task 5 - MSPI-focused C-index | ", cohort_label)
  )

  write.csv(
    df_c,
    file.path(tables_dir, paste0(file_prefix, "_T5_cindex_", output_date, ".csv")),
    row.names = FALSE
  )

  save_table_docx(
    df_c,
    file.path(tables_dir, paste0(file_prefix, "_T5_cindex_", output_date, ".docx")),
    title = "Task 5 - C-index Comparison"
  )

  invisible(list(
    c_dens = ci_d$c, c_m1 = ci_m1$c, c_m2 = ci_m2$c,
    c_mspi = ci_mo$c, c_dens_mspi = ci_dm$c,
    delta_c = delta_c, lrt_p = lrt_p
  ))
}





## ═══════════════════════════════════════════════════════════════════════════
## MAIN RUNS
## ═══════════════════════════════════════════════════════════════════════════

## Main Discovery/Validation cohorts
t1_by <- run_task1(disc_young, "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                   "Discovery Young Luminal", "DDFS", "Discovery_Young")
t1_be <- run_task1(disc_early, "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                   "Discovery Early Luminal", "DDFS", "Discovery_Early")
t1_bs <- run_task1(disc_s4,    "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                   "Discovery S4", "DDFS", "Discovery_S4")
## Validation cohort Task 1 — all 7 subgroups (OS only)
t1_val_inv   <- run_task1(val_inv,       "OS_months", "status_os", "M1", "M2", "CD3_density",
                           "Validation cohort All Invasive", "OS", "VAL_AllInv")
t1_val_inv_y <- tryCatch(
  run_task1(val_inv_young, "OS_months", "status_os", "M1", "M2", "CD3_density",
            "Validation cohort All Inv Young", "OS", "VAL_AllInv_Young"),
  error = function(e) { cat("  T1 All Inv Young error:", e$message, "\n"); NULL })
t1_val_inv_o <- tryCatch(
  run_task1(val_inv_old,   "OS_months", "status_os", "M1", "M2", "CD3_density",
            "Validation cohort All Inv Old", "OS", "VAL_AllInv_Old"),
  error = function(e) { cat("  T1 All Inv Old error:", e$message, "\n"); NULL })
t1_val_erp   <- run_task1(val_erp,       "OS_months", "status_os", "M1", "M2", "CD3_density",
                           "Validation cohort ER+", "OS", "VAL_ERpos")
t1_val_ey    <- tryCatch(
  run_task1(val_young,     "OS_months", "status_os", "M1", "M2", "CD3_density",
            "Validation cohort ER+ Young", "OS", "VAL_ERpos_Young"),
  error = function(e) { cat("  T1 ER+ Young error:", e$message, "\n"); NULL })
t1_val_erp_N0 <- tryCatch(
  run_task1(val_erp_N0,   "OS_months", "status_os", "M1", "M2", "CD3_density",
            "Validation cohort ER+ N0", "OS", "VAL_ERpos_N0"),
  error = function(e) { cat("  T1 ER+ N0 error:", e$message, "\n"); NULL })
t1_val_ey_N0  <- tryCatch(
  run_task1(val_young_N0,  "OS_months", "status_os", "M1", "M2", "CD3_density",
            "Validation cohort ER+ Young N0", "OS", "VAL_ERpos_Young_N0"),
  error = function(e) { cat("  T1 ER+ Young N0 error:", e$message, "\n"); NULL })
t1_val_ern   <- tryCatch(
  run_task1(val_ern,       "OS_months", "status_os", "M1", "M2", "CD3_density",
            "Validation cohort ER-", "OS", "VAL_ERneg"),
  error = function(e) { cat("  T1 ER-neg error:", e$message, "\n"); NULL })
t1_val_dcis  <- tryCatch(
  run_task1(val_dcis,      "OS_months", "status_os", "M1", "M2", "CD3_density",
            "Validation cohort DCIS", "OS", "VAL_DCIS"),
  error = function(e) { cat("  T1 DCIS error:", e$message, "\n"); NULL })

## ── v4: Mini forest plots (UV + MV density-adjusted) ────────────────────
## Each produces a compact 2-bar horizontal strip for placement above KMs.

cat("\n─── v4: Generating mini forest plots + pie charts ───\n")

mv_var_labels <- c(
  "MSPI_sd"   = "MSPI",
  "dm1"       = "M1 density",
  "dm2"       = "M2 density",
  "dcd3"      = "T cell\ndensity"
)

for (info in list(
  list(r = t1_by, lab = "Discovery Young Luminal", ep = "DDFS", pfx = "Discovery_Young"),
  list(r = t1_be, lab = "Discovery Early Luminal", ep = "DDFS", pfx = "Discovery_Early"),
  list(r = t1_bs, lab = "Discovery S4",            ep = "DDFS", pfx = "Discovery_S4"),
  list(r = t1_val_inv,   lab = "Validation cohort All Invasive",  ep = "OS", pfx = "VAL_AllInv"),
  list(r = t1_val_inv_y, lab = "Validation cohort All Inv Young", ep = "OS", pfx = "VAL_AllInv_Young"),
  list(r = t1_val_inv_o, lab = "Validation cohort All Inv Old",   ep = "OS", pfx = "VAL_AllInv_Old"),
  list(r = t1_val_erp,   lab = "Validation cohort ER+",           ep = "OS", pfx = "VAL_ERpos"),
  list(r = t1_val_ey,    lab = "Validation cohort ER+ Young",     ep = "OS", pfx = "VAL_ERpos_Young"),
  list(r = t1_val_erp_N0, lab = "Validation cohort ER+ N0",       ep = "OS", pfx = "VAL_ERpos_N0"),
  list(r = t1_val_ey_N0,  lab = "Validation cohort ER+ Young N0", ep = "OS", pfx = "VAL_ERpos_Young_N0"),
  list(r = t1_val_ern,   lab = "Validation cohort ER-",           ep = "OS", pfx = "VAL_ERneg"),
  list(r = t1_val_dcis,  lab = "Validation cohort DCIS",          ep = "OS", pfx = "VAL_DCIS")
)) {
  R <- info$r
  if (!is.null(R) && !isTRUE(R$skipped) && !is.null(R$model_uv) && !is.null(R$model_mv)) {
    tryCatch({
      make_mini_forest(R$model_uv, R$model_mv, "MSPI_sd",
                       info$lab, info$ep, info$pfx)
      cat("  Mini forest:", info$lab, info$ep, "\n")
    }, error = function(e) cat("  Mini forest ERROR:", info$lab, e$message, "\n"))

    tryCatch({
      make_mv_pie(R$model_mv, mv_var_labels, info$lab, info$pfx)
      cat("  Pie chart:", info$lab, info$ep, "\n")
    }, error = function(e) cat("  Pie chart ERROR:", info$lab, e$message, "\n"))
  }
}

cat("─── v4: Mini forests + pies complete ───\n\n")

t2_by <- run_task2(disc_young, "DDFStime2021", "status_num", "Discovery Young Luminal", "DDFS", "Discovery_Young")
t2_be <- run_task2(disc_early, "DDFStime2021", "status_num", "Discovery Early Luminal", "DDFS", "Discovery_Early")
t2_bs <- run_task2(disc_s4,    "DDFStime2021", "status_num", "Discovery S4", "DDFS", "Discovery_S4")
## Validation cohort Task 2 — all 7 subgroups (OS only)
t2_val_inv   <- run_task2(val_inv,       "OS_months", "status_os", "Validation cohort All Invasive",  "OS", "VAL_AllInv")
t2_val_inv_y <- tryCatch(
  run_task2(val_inv_young, "OS_months", "status_os", "Validation cohort All Inv Young", "OS", "VAL_AllInv_Young"),
  error = function(e) { cat("  T2 All Inv Young error:", e$message, "\n"); NULL })
t2_val_inv_o <- tryCatch(
  run_task2(val_inv_old,   "OS_months", "status_os", "Validation cohort All Inv Old",   "OS", "VAL_AllInv_Old"),
  error = function(e) { cat("  T2 All Inv Old error:", e$message, "\n"); NULL })
t2_val_erp   <- run_task2(val_erp,       "OS_months", "status_os", "Validation cohort ER+",           "OS", "VAL_ERpos")
t2_val_ey    <- tryCatch(
  run_task2(val_young,     "OS_months", "status_os", "Validation cohort ER+ Young",     "OS", "VAL_ERpos_Young"),
  error = function(e) { cat("  T2 ER+ Young error:", e$message, "\n"); NULL })
t2_val_erp_N0 <- tryCatch(
  run_task2(val_erp_N0,   "OS_months", "status_os", "Validation cohort ER+ N0",        "OS", "VAL_ERpos_N0"),
  error = function(e) { cat("  T2 ER+ N0 error:", e$message, "\n"); NULL })
t2_val_ey_N0  <- tryCatch(
  run_task2(val_young_N0,  "OS_months", "status_os", "Validation cohort ER+ Young N0", "OS", "VAL_ERpos_Young_N0"),
  error = function(e) { cat("  T2 ER+ Young N0 error:", e$message, "\n"); NULL })
t2_val_ern   <- tryCatch(
  run_task2(val_ern,       "OS_months", "status_os", "Validation cohort ER-",           "OS", "VAL_ERneg"),
  error = function(e) { cat("  T2 ER-neg error:", e$message, "\n"); NULL })
t2_val_dcis  <- tryCatch(
  run_task2(val_dcis,      "OS_months", "status_os", "Validation cohort DCIS",          "OS", "VAL_DCIS"),
  error = function(e) { cat("  T2 DCIS error:", e$message, "\n"); NULL })

t3_by <- run_task3(disc_young, "DDFStime2021", "status_num", "Discovery Young Luminal", "DDFS", "Discovery_Young")
t3_be <- run_task3(disc_early, "DDFStime2021", "status_num", "Discovery Early Luminal", "DDFS", "Discovery_Early")
t3_bs <- run_task3(disc_s4,    "DDFStime2021", "status_num", "Discovery S4", "DDFS", "Discovery_S4")
t3_enc <- run_task3(val_erp, "OS_months", "status_os", "Validation cohort ER+", "OS", "VAL_ERpos")
t3_val_N0 <- tryCatch(
  run_task3(val_erp_N0, "OS_months", "status_os", "Validation cohort ER+ N0", "OS", "VAL_ERpos_N0"),
  error = function(e) { cat("  T3 ER+ N0 error:", e$message, "\n"); NULL })

pool_early <- prep_pool(disc_early, "DDFStime2021", "status_num", "Discovery_EarlyLuminal")
pool_young <- prep_pool(disc_young, "DDFStime2021", "status_num", "Discovery_YoungLuminal")
pool_enc   <- prep_pool(val_erp,   "OS_months", "status_os", "Validation cohort_ER_pos")

pool_A <- rbind(pool_early, pool_enc)
pool_A$Cohort <- factor(pool_A$Cohort)

pool_val_young <- prep_pool(val_young, "OS_months", "status_os", "Validation cohort_ER_pos_Young")
pool_B <- rbind(pool_young, pool_val_young)
pool_B$Cohort <- factor(pool_B$Cohort)

t4_A <- run_task4(pool_A, "Pool A: Discovery Early Luminal + Validation cohort ER+", "Pool_A")
t4_B <- run_task4(pool_B, "Pool B: Discovery Young Luminal + Validation cohort ER+ Young <=50y", "Pool_B")

## Pooled IPD — N0 subsets
pool_val_N0       <- prep_pool(val_erp_N0,   "OS_months", "status_os", "Validation cohort_ER_pos_N0")
pool_A_N0         <- rbind(pool_early, pool_val_N0)
pool_A_N0$Cohort  <- factor(pool_A_N0$Cohort)
pool_val_young_N0 <- prep_pool(val_young_N0, "OS_months", "status_os", "Validation cohort_ER_pos_Young_N0")
pool_B_N0         <- rbind(pool_young, pool_val_young_N0)
pool_B_N0$Cohort  <- factor(pool_B_N0$Cohort)
t4_A_N0 <- tryCatch(
  run_task4(pool_A_N0, "Pool A N0: Discovery Early Luminal + Validation cohort ER+ N0", "Pool_A_N0"),
  error = function(e) { cat("  T4 Pool A N0 error:", e$message, "\n"); NULL })
t4_B_N0 <- tryCatch(
  run_task4(pool_B_N0, "Pool B N0: Discovery Young Luminal + Validation cohort ER+ Young N0", "Pool_B_N0"),
  error = function(e) { cat("  T4 Pool B N0 error:", e$message, "\n"); NULL })

t5_by <- run_task5(disc_young, "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                   FALSE, "Discovery Young Luminal", "DDFS", "Discovery_Young")
t5_be <- run_task5(disc_early, "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                   FALSE, "Discovery Early Luminal", "DDFS", "Discovery_Early")
t5_bs <- run_task5(disc_s4,    "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                   FALSE, "Discovery S4", "DDFS", "Discovery_S4")
## Validation cohort Task 5 — all 7 subgroups (OS only)
t5_val_inv   <- run_task5(val_inv,   "OS_months", "status_os", "M1", "M2", "CD3_density",
                           TRUE, "Validation cohort All Invasive", "OS", "VAL_AllInv")
t5_val_inv_y <- tryCatch(
  run_task5(val_inv_young, "OS_months", "status_os", "M1", "M2", "CD3_density",
            TRUE, "Validation cohort All Inv Young", "OS", "VAL_AllInv_Young"),
  error = function(e) { cat("  T5 All Inv Young error:", e$message, "\n"); NULL })
t5_val_inv_o <- tryCatch(
  run_task5(val_inv_old,   "OS_months", "status_os", "M1", "M2", "CD3_density",
            TRUE, "Validation cohort All Inv Old", "OS", "VAL_AllInv_Old"),
  error = function(e) { cat("  T5 All Inv Old error:", e$message, "\n"); NULL })
t5_val_erp   <- run_task5(val_erp,   "OS_months", "status_os", "M1", "M2", "CD3_density",
                           TRUE, "Validation cohort ER+", "OS", "VAL_ERpos")
t5_val_ey    <- tryCatch(
  run_task5(val_young,     "OS_months", "status_os", "M1", "M2", "CD3_density",
            TRUE, "Validation cohort ER+ Young", "OS", "VAL_ERpos_Young"),
  error = function(e) { cat("  T5 ER+ Young error:", e$message, "\n"); NULL })
t5_val_erp_N0 <- tryCatch(
  run_task5(val_erp_N0,   "OS_months", "status_os", "M1", "M2", "CD3_density",
            TRUE, "Validation cohort ER+ N0", "OS", "VAL_ERpos_N0"),
  error = function(e) { cat("  T5 ER+ N0 error:", e$message, "\n"); NULL })
t5_val_ey_N0  <- tryCatch(
  run_task5(val_young_N0,  "OS_months", "status_os", "M1", "M2", "CD3_density",
            TRUE, "Validation cohort ER+ Young N0", "OS", "VAL_ERpos_Young_N0"),
  error = function(e) { cat("  T5 ER+ Young N0 error:", e$message, "\n"); NULL })




## ═══════════════════════════════════════════════════════════════════════════
## v6: Validation cohort ER+ N03 — AGE-STRATIFIED ANALYSIS (3 groups)
## Same logic as previous N0 block, but WITHOUT N0 prefiltering:
## uses full val_erp and labels outputs as N03
## Runs T1 + T5 per age group, then produces COMBINED mini forest and
## COMBINED C-index MSPIfocused plots (all 3 age groups on one graph each)
## ═══════════════════════════════════════════════════════════════════════════

cat("\n", strrep("─", 60), "\n")
cat("v6: Validation cohort ER+ N03 age-stratified analysis (<=50, >50-<65, >=65)\n")
cat(strrep("─", 60), "\n")

val_N03_age_groups <- list(
  list(filt = "<=50",    lab = "Age <=50",  pfx = "VAL_ERpos_N03_Age50"),
  list(filt = ">50-<65", lab = "Age 51-64", pfx = "VAL_ERpos_N03_Age51_64"),
  list(filt = ">=65",    lab = "Age >=65",  pfx = "VAL_ERpos_N03_Age65")
)

t1_N03_age <- list()
t5_N03_age <- list()
forest_rows <- data.frame()
cindex_rows <- data.frame()

for (ag in val_N03_age_groups) {
  sub <- val_erp %>% dplyr::filter(Age == ag$filt)
  n_sub <- nrow(sub)
  ev_sub <- sum(sub$status_os)
  cat(sprintf("  %s: n=%d  OS events=%d\n", ag$lab, n_sub, ev_sub))
  
  if (n_sub < 10 || ev_sub < 5) {
    cat("    Skipping - too few patients/events\n")
    t1_N03_age[[ag$pfx]] <- NULL
    t5_N03_age[[ag$pfx]] <- NULL
    next
  }
  
  ## Task 1 — individual per-age outputs (tables, KMs) + collect UV/MV for forest
  t1_res <- tryCatch(
    run_task1(sub, "OS_months", "status_os", "M1", "M2", "CD3_density",
              paste0("Validation cohort ER+ N03 ", ag$lab), "OS", ag$pfx),
    error = function(e) { cat("    T1 error:", e$message, "\n"); NULL })
  t1_N03_age[[ag$pfx]] <- t1_res
  
  ## Extract HR/CI/p from UV and MV models for combined mini forest
  if (!is.null(t1_res) && !isTRUE(t1_res$skipped) &&
      !is.null(t1_res$model_uv) && !is.null(t1_res$model_mv)) {
    for (mod_info in list(
      list(cx = t1_res$model_uv, mtype = "Univariable"),
      list(cx = t1_res$model_mv, mtype = "Density-adjusted")
    )) {
      s <- summary(mod_info$cx)
      i <- which(rownames(s$coefficients) == "MSPI_sd")
      if (length(i) == 0) i <- 1
      forest_rows <- rbind(forest_rows, data.frame(
        AgeGroup = ag$lab,
        AgeN     = sprintf("%s (n=%d, %d ev)", ag$lab, n_sub, ev_sub),
        Model    = mod_info$mtype,
        HR       = s$conf.int[i, 1],
        CI_lo    = s$conf.int[i, 3],
        CI_hi    = s$conf.int[i, 4],
        p        = s$coefficients[i, 5],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  ## Task 5 — individual per-age outputs + collect C-index for combined plot
  t5_res <- tryCatch(
    run_task5(sub, "OS_months", "status_os", "M1", "M2", "CD3_density",
              TRUE, paste0("Validation cohort ER+ N03 ", ag$lab), "OS", ag$pfx),
    error = function(e) { cat("    T5 error:", e$message, "\n"); NULL })
  t5_N03_age[[ag$pfx]] <- t5_res
  
  if (!is.null(t5_res)) {
    cindex_rows <- rbind(cindex_rows, data.frame(
      AgeGroup = ag$lab,
      AgeN     = sprintf("%s (n=%d)", ag$lab, n_sub),
      Model    = c("MSPI only", "Density + MSPI"),
      C        = c(t5_res$c_mspi, t5_res$c_dens_mspi),
      stringsAsFactors = FALSE
    ))
  }
}

## ── COMBINED MINI FOREST (all age groups on one plot) ────────────────────
if (nrow(forest_rows) > 0) {
  tryCatch({
    forest_rows$AgeN <- factor(forest_rows$AgeN, levels = rev(unique(forest_rows$AgeN)))
    forest_rows$Model <- factor(forest_rows$Model,
                                levels = c("Density-adjusted", "Univariable"))
    forest_rows$sig <- ifelse(forest_rows$p < 0.05, "Significant", "NS")
    
    p_cf <- ggplot(forest_rows, aes(x = HR, y = AgeN, shape = Model, colour = sig)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
      geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.22,
                     linewidth = 0.85, colour = COL_DKGREY,
                     position = position_dodge(width = 0.55)) +
      geom_point(size = 5, position = position_dodge(width = 0.55)) +
      geom_text(aes(label = sprintf("%.2f (%.2f-%.2f) p=%s",
                                    HR, CI_lo, CI_hi, fmt_p(p))),
                hjust = -0.06, size = 3.2, family = FONT_FAM, colour = "black",
                position = position_dodge(width = 0.55)) +
      scale_shape_manual(values = c("Univariable" = 16, "Density-adjusted" = 18)) +
      scale_colour_manual(values = c("Significant" = "#c0392b", "NS" = "grey55"),
                          guide = "none") +
      scale_x_log10() +
      labs(x = "HR per SD (MSPI)", y = NULL, shape = NULL,
           title = "Validation cohort ER+ N03 | MSPI by age group | OS") +
      theme_bw(base_size = 12, base_family = FONT_FAM) +
      theme(plot.margin = margin(4, 80, 4, 4),
            axis.text.y  = element_text(size = 11),
            plot.title   = element_text(size = 11, face = "bold"),
            legend.position = "bottom",
            panel.grid.major.y = element_blank())
    
    n_ages <- length(unique(forest_rows$AgeN))
    save_plot(p_cf, "VAL_ERpos_N03_AgeStrat_mini_forest",
              w = 7, h = 1.2 + n_ages * 0.8)
    cat("  Combined mini forest: OK\n")
  }, error = function(e) cat("  Combined mini forest error:", e$message, "\n"))
}

## ── COMBINED C-INDEX MSPIfocused (all age groups on one plot) ────────────
if (nrow(cindex_rows) > 0) {
  tryCatch({
    ci_detail <- data.frame()
    
    for (ag in val_N03_age_groups) {
      sub <- val_erp %>% dplyr::filter(Age == ag$filt)
      if (nrow(sub) < 10 || sum(sub$status_os) < 5) next
      
      d_c <- sub[!is.na(sub$M1_score) & !is.na(sub$M2_score) &
                   !is.na(sub$OS_months) & sub$OS_months > 0, ]
      d_c$TIME  <- d_c$OS_months
      d_c$EVENT <- as.integer(d_c$status_os)
      mspi_sd <- sd(d_c$M2_score - d_c$M1_score, na.rm = TRUE)
      d_c$MSPI_sd   <- (d_c$M2_score - d_c$M1_score) / mspi_sd
      d_c$dens_m1   <- d_c$M1
      d_c$dens_m2   <- d_c$M2
      d_c$dens_cd3  <- d_c$CD3_density
      
      d_c <- d_c[complete.cases(d_c[, c("MSPI_sd", "dens_m1", "dens_m2", "dens_cd3",
                                        "TIME", "EVENT")]), ]
      
      cx_mo <- coxph(Surv(TIME, EVENT) ~ MSPI_sd, data = d_c, x = TRUE, y = TRUE)
      cx_dm <- coxph(Surv(TIME, EVENT) ~ dens_m1 + dens_m2 + dens_cd3 + MSPI_sd,
                     data = d_c, x = TRUE, y = TRUE)
      
      get_ci_b <- function(cx, dat, B = 200, s = 12345) {
        c0 <- concordance(cx)$concordance
        set.seed(s)
        boots <- replicate(B, {
          idx <- sample(seq_len(nrow(dat)), replace = TRUE)
          tryCatch(concordance(update(cx, data = dat[idx, ]))$concordance,
                   error = function(e) NA_real_)
        })
        boots <- boots[!is.na(boots)]
        list(c = round(c0, 4),
             lo = round(quantile(boots, 0.025), 4),
             hi = round(quantile(boots, 0.975), 4))
      }
      
      ci_mo_b <- get_ci_b(cx_mo, d_c)
      ci_dm_b <- get_ci_b(cx_dm, d_c)
      n_sub <- nrow(d_c)
      
      ci_detail <- rbind(ci_detail, data.frame(
        AgeGroup = sprintf("%s (n=%d)", ag$lab, n_sub),
        Model    = c("MSPI only", "Density + MSPI"),
        C        = c(ci_mo_b$c, ci_dm_b$c),
        CI_lo    = c(ci_mo_b$lo, ci_dm_b$lo),
        CI_hi    = c(ci_mo_b$hi, ci_dm_b$hi),
        AIC      = round(c(AIC(cx_mo), AIC(cx_dm)), 1),
        stringsAsFactors = FALSE
      ))
    }
    
    if (nrow(ci_detail) > 0) {
      ci_detail$AgeGroup <- factor(ci_detail$AgeGroup,
                                   levels = unique(ci_detail$AgeGroup))
      ci_detail$Model <- factor(ci_detail$Model,
                                levels = c("MSPI only", "Density + MSPI"))
      
      p_cc <- ggplot(ci_detail, aes(x = AgeGroup, y = C, ymin = CI_lo, ymax = CI_hi,
                                    colour = Model, shape = Model)) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
        geom_errorbar(width = 0.2, linewidth = 0.9,
                      position = position_dodge(width = 0.4)) +
        geom_point(size = 5, position = position_dodge(width = 0.4)) +
        geom_text(aes(label = sprintf("C=%.3f\nAIC=%.0f", C, AIC)),
                  vjust = -0.8, size = 3.4, family = FONT_FAM, lineheight = 0.95,
                  position = position_dodge(width = 0.4), colour = "black") +
        scale_colour_manual(values = c("MSPI only" = COL_NAVY,
                                       "Density + MSPI" = COL_DKGREEN)) +
        scale_shape_manual(values = c("MSPI only" = 16, "Density + MSPI" = 18)) +
        labs(x = NULL, y = "C-index (Harrell) with bootstrap 95% CI",
             colour = NULL, shape = NULL,
             title = "Validation cohort ER+ N03 | MSPI-focused C-index by age group | OS",
             caption = "Bootstrap 200 iter. | AIC lower = better") +
        theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
        theme(legend.position = "bottom",
              axis.text.x = element_text(size = 11),
              panel.grid.major.x = element_blank()) +
        ylim(min(ci_detail$CI_lo) - 0.03, max(ci_detail$CI_hi) + 0.1)
      
      save_plot(p_cc, "VAL_ERpos_N03_AgeStrat_cindex_MSPIfocused",
                w = 5, h = 5.5)
      cat("  Combined C-index: OK\n")
      
      write.csv(ci_detail,
                file.path(tables_dir,
                          paste0("VAL_ERpos_N03_AgeStrat_cindex_", output_date, ".csv")),
                row.names = FALSE)

      save_table_docx(ci_detail,
                      file.path(tables_dir,
                                paste0("VAL_ERpos_N03_AgeStrat_cindex_", output_date, ".docx")),
                      title = "Validation cohort ER+ Age-Stratified C-index")
    }
  }, error = function(e) cat("  Combined C-index error:", e$message, "\n"))
}

## Save combined forest data
if (nrow(forest_rows) > 0) {
  write.csv(forest_rows,
            file.path(tables_dir,
                      paste0("VAL_ERpos_N03_AgeStrat_forest_", output_date, ".csv")),
            row.names = FALSE)

  save_table_docx(forest_rows,
                  file.path(tables_dir,
                            paste0("VAL_ERpos_N03_AgeStrat_forest_", output_date, ".docx")),
                  title = "Validation cohort ER+ Age-Stratified Forest")
}

cat("v6: ER+ N03 age-stratified analysis complete\n\n")





t5_val_ern   <- tryCatch(
  run_task5(val_ern,       "OS_months", "status_os", "M1", "M2", "CD3_density",
            FALSE, "Validation cohort ER-", "OS", "VAL_ERneg"),
  error = function(e) { cat("  T5 ER-neg error:", e$message, "\n"); NULL })
t5_val_dcis  <- tryCatch(
  run_task5(val_dcis,      "OS_months", "status_os", "M1", "M2", "CD3_density",
            FALSE, "Validation cohort DCIS", "OS", "VAL_DCIS"),
  error = function(e) { cat("  T5 DCIS error:", e$message, "\n"); NULL })

## (Validation cohort ER-neg and DCIS analyses now integrated into main T1/T2/T5/T6 blocks above)

## ═══════════════════════════════════════════════════════════════════════════
## ADDED: Discovery MOLECULAR SUBTYPE ANALYSES (individual subtypes)
## ═══════════════════════════════════════════════════════════════════════════

cat("\n", strrep("═", 70), "\n")
cat("ADDED MODULE: Discovery MOLECULAR SUBTYPE ANALYSES (individual)\n")
cat(strrep("═", 70), "\n\n")
cat("Discovery individual subtype analyses\n", file = mlog, append = TRUE)

subtype_master <- list()

for (st in names(disc_subtype_list)) {
  ds <- disc_subtype_list[[st]]
  prefix <- paste0("Discovery_Subtype_", safe_name(st))
  label  <- paste0("Discovery subtype: ", st)

  cat("\nRunning individual subtype analyses for:", st, "\n")

  r1 <- run_task1(ds, "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                  label, "DDFS", prefix)
  r2 <- run_task2(ds, "DDFStime2021", "status_num",
                  label, "DDFS", prefix)

  ## Individual M1 and M2 affinity Cox (continuous per SD + density-adjusted)
  d_ind <- ds[!is.na(ds$M1_score) & !is.na(ds$M2_score) &
                !is.na(ds$DDFStime2021) & ds$DDFStime2021 > 0 &
                !is.na(ds$CD3_density), ]
  d_ind$TIME  <- d_ind$DDFStime2021
  d_ind$EVENT <- as.integer(d_ind$status_num)
  d_ind$M1_sd <- as.numeric(scale(d_ind$M1_score))
  d_ind$M2_sd <- as.numeric(scale(d_ind$M2_score))
  m2_uv <- m2_mv <- m1_uv <- m1_mv <- NULL
  if (nrow(d_ind) >= 20 && sum(d_ind$EVENT) >= 5) {
    tryCatch({
      cx_m2u <- coxph(Surv(TIME, EVENT) ~ M2_sd, data = d_ind)
      s2 <- summary(cx_m2u)
      m2_uv <- list(HR = s2$conf.int[1,1], CI_lo = s2$conf.int[1,3],
                     CI_hi = s2$conf.int[1,4], p = s2$coefficients[1,5])
    }, error = function(e) NULL)
    tryCatch({
      cx_m2m <- coxph(Surv(TIME, EVENT) ~ M2_sd + CD3_density + M1 + M2, data = d_ind)
      s2m <- summary(cx_m2m)
      i <- which(rownames(s2m$coefficients) == "M2_sd")
      m2_mv <- list(HR = s2m$conf.int[i,1], CI_lo = s2m$conf.int[i,3],
                     CI_hi = s2m$conf.int[i,4], p = s2m$coefficients[i,5])
    }, error = function(e) NULL)
    tryCatch({
      cx_m1u <- coxph(Surv(TIME, EVENT) ~ M1_sd, data = d_ind)
      s1 <- summary(cx_m1u)
      m1_uv <- list(HR = s1$conf.int[1,1], CI_lo = s1$conf.int[1,3],
                     CI_hi = s1$conf.int[1,4], p = s1$coefficients[1,5])
    }, error = function(e) NULL)
    tryCatch({
      cx_m1m <- coxph(Surv(TIME, EVENT) ~ M1_sd + CD3_density + M1 + M2, data = d_ind)
      s1m <- summary(cx_m1m)
      i <- which(rownames(s1m$coefficients) == "M1_sd")
      m1_mv <- list(HR = s1m$conf.int[i,1], CI_lo = s1m$conf.int[i,3],
                     CI_hi = s1m$conf.int[i,4], p = s1m$coefficients[i,5])
    }, error = function(e) NULL)
  }

  subtype_master[[st]] <- list(task1 = r1, task2 = r2,
                                m2_uv = m2_uv, m2_mv = m2_mv,
                                m1_uv = m1_uv, m1_mv = m1_mv,
                                n_ind = nrow(d_ind), ev_ind = sum(d_ind$EVENT))
}

## ═══════════════════════════════════════════════════════════════════════════
## ADDED: Discovery MERGED MOLECULAR SUBTYPE ANALYSES (Luminal, HER2+, TNBC)
## ═══════════════════════════════════════════════════════════════════════════

cat("\n", strrep("═", 70), "\n")
cat("ADDED MODULE: Discovery MERGED MOLECULAR SUBTYPE ANALYSES\n")
cat(strrep("═", 70), "\n\n")
cat("Discovery merged subtype analyses\n", file = mlog, append = TRUE)

merged_master <- list()

for (mg in names(disc_merged_list)) {
  ds <- disc_merged_list[[mg]]
  prefix <- paste0("Discovery_Merged_", safe_name(mg))
  label  <- paste0("Discovery ", mg)

  cat("\nRunning merged subtype analyses for:", mg,
      "(n=", nrow(ds), ", events=", sum(ds$status_num), ")\n")

  r1 <- run_task1(ds, "DDFStime2021", "status_num", "M1", "M2", "CD3_density",
                  label, "DDFS", prefix)
  r2 <- run_task2(ds, "DDFStime2021", "status_num",
                  label, "DDFS", prefix)

  ## Individual M1 and M2 affinity Cox for merged subtypes
  d_ind <- ds[!is.na(ds$M1_score) & !is.na(ds$M2_score) &
                !is.na(ds$DDFStime2021) & ds$DDFStime2021 > 0 &
                !is.na(ds$CD3_density), ]
  d_ind$TIME  <- d_ind$DDFStime2021
  d_ind$EVENT <- as.integer(d_ind$status_num)
  d_ind$M1_sd <- as.numeric(scale(d_ind$M1_score))
  d_ind$M2_sd <- as.numeric(scale(d_ind$M2_score))
  m2_uv <- m2_mv <- m1_uv <- m1_mv <- NULL
  if (nrow(d_ind) >= 20 && sum(d_ind$EVENT) >= 5) {
    tryCatch({
      cx <- coxph(Surv(TIME, EVENT) ~ M2_sd, data = d_ind)
      s <- summary(cx)
      m2_uv <- list(HR = s$conf.int[1,1], CI_lo = s$conf.int[1,3],
                     CI_hi = s$conf.int[1,4], p = s$coefficients[1,5])
    }, error = function(e) NULL)
    tryCatch({
      cx <- coxph(Surv(TIME, EVENT) ~ M2_sd + CD3_density + M1 + M2, data = d_ind)
      s <- summary(cx)
      i <- which(rownames(s$coefficients) == "M2_sd")
      m2_mv <- list(HR = s$conf.int[i,1], CI_lo = s$conf.int[i,3],
                     CI_hi = s$conf.int[i,4], p = s$coefficients[i,5])
    }, error = function(e) NULL)
    tryCatch({
      cx <- coxph(Surv(TIME, EVENT) ~ M1_sd, data = d_ind)
      s <- summary(cx)
      m1_uv <- list(HR = s$conf.int[1,1], CI_lo = s$conf.int[1,3],
                     CI_hi = s$conf.int[1,4], p = s$coefficients[1,5])
    }, error = function(e) NULL)
    tryCatch({
      cx <- coxph(Surv(TIME, EVENT) ~ M1_sd + CD3_density + M1 + M2, data = d_ind)
      s <- summary(cx)
      i <- which(rownames(s$coefficients) == "M1_sd")
      m1_mv <- list(HR = s$conf.int[i,1], CI_lo = s$conf.int[i,3],
                     CI_hi = s$conf.int[i,4], p = s$coefficients[i,5])
    }, error = function(e) NULL)
  }

  merged_master[[mg]] <- list(task1 = r1, task2 = r2,
                               m2_uv = m2_uv, m2_mv = m2_mv,
                               m1_uv = m1_uv, m1_mv = m1_mv,
                               n_ind = nrow(d_ind), ev_ind = sum(d_ind$EVENT))

  ## v4: Mini forests + pies for merged subtype MSPI
  if (!is.null(r1) && !isTRUE(r1$skipped) && !is.null(r1$model_uv) && !is.null(r1$model_mv)) {
    tryCatch({
      make_mini_forest(r1$model_uv, r1$model_mv, "MSPI_sd",
                       paste0("Discovery ", mg), "DDFS", paste0("Discovery_Merged_", safe_name(mg)))
      make_mv_pie(r1$model_mv, mv_var_labels, paste0("Discovery ", mg),
                  paste0("Discovery_Merged_", safe_name(mg)))
    }, error = function(e) cat("  Subtype mini-forest/pie error:", mg, e$message, "\n"))
  }

  ## (Task 6 for merged subtypes is handled below in the dedicated T6 section
  ##  via t6_bl — no duplicate call here to avoid empty/redundant PDFs)
}

################################################################################
##
##  TASK 6: CLINICAL MULTIVARIABLE COX REGRESSION
##
################################################################################

cat("\n", strrep("═", 70), "\n")
cat("TASK 6: CLINICAL MULTIVARIABLE COX REGRESSION\n")
cat(strrep("═", 70), "\n\n")
cat("TASK 6 — Clinical MV\n", file = mlog, append = TRUE)

## run_task6: Fits clinical MV Cox models with MSPI (or MSPI_Q4) + clinical covariates
## Returns a list with model summaries suitable for forest plot and manuscript table

run_task6 <- function(dat, time_col, event_col,
                      clinical_vars,        # named list: varname -> label
                      cohort_label, endpoint_label, file_prefix) {

  d <- dat[!is.na(dat$M1_score) & !is.na(dat$M2_score) &
             !is.na(dat[[time_col]]) & dat[[time_col]] > 0, ]
  d$TIME  <- d[[time_col]]
  d$EVENT <- as.integer(d[[event_col]])
  d$MSPI_sd <- as.numeric(scale(d$M2_score - d$M1_score))

  ## Median split MSPI for dichotomized model
  mspi_raw <- d$M2_score - d$M1_score
  d$MSPI_high <- factor(
    ifelse(mspi_raw >= median(mspi_raw, na.rm = TRUE), "High", "Low"),
    levels = c("Low", "High")
  )

  ## Build formula with available clinical variables
  clin_available <- character(0)
  for (v in names(clinical_vars)) {
    if (v %in% names(d) && sum(!is.na(d[[v]])) >= 20) {
      clin_available <- c(clin_available, v)
    } else {
      cat(sprintf("  %s: variable '%s' (%s) not available or too few non-NA, skipping\n",
                  cohort_label, v, clinical_vars[[v]]))
    }
  }

  if (length(clin_available) == 0) {
    cat("  No clinical variables available for", cohort_label, "- skipping Task 6\n")
    return(invisible(NULL))
  }

  ## Complete cases on all needed variables
  keep_vars <- c("TIME", "EVENT", "MSPI_sd", "MSPI_high", clin_available)
  d_cc <- d[complete.cases(d[, keep_vars]), ]
  n <- nrow(d_cc)
  ev <- sum(d_cc$EVENT)

  cat(sprintf("  %s: n=%d events=%d | Clinical vars: %s\n",
              cohort_label, n, ev, paste(clin_available, collapse = ", ")))

  if (n < 20 || ev < 10) {
    cat("  Insufficient data for clinical MV, skipping\n")
    return(invisible(NULL))
  }

  ## Model 1: MSPI continuous + clinical
  fml_cont <- as.formula(paste0(
    "Surv(TIME, EVENT) ~ MSPI_sd + ",
    paste(clin_available, collapse = " + ")
  ))
  cx_cont <- tryCatch(coxph(fml_cont, data = d_cc), error = function(e) NULL)

  ## Model 2: MSPI dichotomized (median) + clinical
  fml_dich <- as.formula(paste0(
    "Surv(TIME, EVENT) ~ MSPI_high + ",
    paste(clin_available, collapse = " + ")
  ))
  cx_dich <- tryCatch(coxph(fml_dich, data = d_cc), error = function(e) NULL)

  ## Extract results
  extract_cox <- function(cx) {
    if (is.null(cx)) return(NULL)
    s <- summary(cx)
    coefs <- s$coefficients
    ci <- s$conf.int
    vars <- rownames(coefs)

    data.frame(
      Variable = vars,
      HR = round(coefs[, "exp(coef)"], 3),
      CI_lo = round(ci[, "lower .95"], 3),
      CI_hi = round(ci[, "upper .95"], 3),
      p = coefs[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }

  res_cont <- extract_cox(cx_cont)
  res_dich <- extract_cox(cx_dich)

  ## Save CSV outputs
  if (!is.null(res_cont)) {
    write.csv(res_cont,
      file.path(tables_dir, paste0(file_prefix, "_T6_clinical_MV_continuous_", output_date, ".csv")),
      row.names = FALSE)

    save_table_docx(res_cont,
      file.path(tables_dir, paste0(file_prefix, "_T6_clinical_MV_continuous_", output_date, ".docx")),
      title = "Task 6 - Clinical MV Continuous")
  }
  if (!is.null(res_dich)) {
    write.csv(res_dich,
      file.path(tables_dir, paste0(file_prefix, "_T6_clinical_MV_dichotomized_", output_date, ".csv")),
      row.names = FALSE)

    save_table_docx(res_dich,
      file.path(tables_dir, paste0(file_prefix, "_T6_clinical_MV_dichotomized_", output_date, ".docx")),
      title = "Task 6 - Clinical MV Dichotomized")
  }

  ## Forest plot for continuous model
  ## Standardised layout: fixed row height (0.55 in), fixed total width (7 in),
  ## fixed x-axis range and y-axis label width so Discovery and Validation plots are
  ## visually comparable regardless of number of covariates.
  ## All levels of multi-level factors are shown (against their reference).
  if (!is.null(res_cont) && nrow(res_cont) >= 2) {
    tryCatch({
      pd <- res_cont

      ## ── Label mapping ──────────────────────────────────────────────────
      ## Multi-level factors (e.g. Stage_TNM_44IIa, diadtm_final1990-1999)
      ## produce many expanded coefficient rows. Use PREFIX matching: for
      ## each variable name from coxph output, find the clinical_vars key
      ## that is the longest prefix, then format as "Label: level".
      map_label <- function(v) {
        if (v == "MSPI_sd") return("MSPI (per SD)")
        if (v == "MSPI_highHigh") return("MSPI High vs Low")
        ## Try exact match (variable is not a factor, or has only 2 levels
        ## and coxph used the variable name directly)
        if (v %in% names(clinical_vars)) return(clinical_vars[[v]])
        ## Prefix match: find the longest clinical_vars key that is a prefix
        matches <- names(clinical_vars)[sapply(names(clinical_vars),
                                               function(k) startsWith(v, k))]
        if (length(matches) > 0) {
          best <- matches[which.max(nchar(matches))]
          suffix <- substring(v, nchar(best) + 1)
          return(paste0(clinical_vars[[best]], ": ", suffix))
        }
        return(v)
      }
      pd$label <- sapply(pd$Variable, map_label, USE.NAMES = FALSE)

      ## Remove rows with non-finite or extreme HR/CI (from singular
      ## factor levels with 0 events), which would collapse the log10 axis
      ok <- is.finite(pd$HR) & is.finite(pd$CI_lo) & is.finite(pd$CI_hi) &
            pd$HR > 0 & pd$CI_lo > 0 & pd$CI_hi < 1e4
      if (sum(!ok) > 0) {
        cat(sprintf("  Forest: dropping %d rows with extreme/non-finite HR: %s\n",
                    sum(!ok), paste(pd$Variable[!ok], collapse = ", ")))
        pd <- pd[ok, ]
      }

      if (nrow(pd) < 2) {
        cat("  Forest: fewer than 2 plottable rows after filtering, skipping\n")
      } else {

      pd$label <- factor(pd$label, levels = rev(pd$label))

      pd$is_mspi <- grepl("MSPI", pd$Variable)
      pd$point_col <- ifelse(pd$is_mspi, COL_NAVY, COL_AMBER)
      pd$p_fmt <- sapply(pd$p, fmt_p)

      ## Fixed x-axis limits (log10 scale) to keep forest area the same width
      ## Floor at 0.05 to prevent log10(0) = -Inf from collapsing the axis
      x_lo <- max(0.05, min(0.3, min(pd$CI_lo, na.rm = TRUE) * 0.8))
      x_hi <- min(50, max(3.0, max(pd$CI_hi, na.rm = TRUE) * 1.2))
      ## Padding for HR text annotation
      x_text <- x_hi * 1.15

      p_forest <- ggplot(pd, aes(x = HR, y = label)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey65") +
        geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.22,
                       linewidth = 0.85, colour = COL_DKGREY) +
        geom_point(size = 6, shape = 18, colour = pd$point_col) +
        geom_text(aes(x = x_text,
                      label = sprintf("HR=%.2f (p=%s)", HR, p_fmt)),
                  hjust = 0, size = 3.3, family = FONT_FAM,
                  colour = pd$point_col, show.legend = FALSE) +
        scale_x_log10(limits = c(x_lo, x_text * 1.8),
                      breaks = c(0.5, 1, 2, 5),
                      oob = scales::squish) +
        labs(
          title = paste0("Clinical MV Cox - ", cohort_label, " (", endpoint_label, ")"),
          subtitle = paste0("n=", n, ", events=", ev),
          x = "Hazard Ratio (95% CI)", y = ""
        ) +
        theme_bw(base_size = BASE_SIZE, base_family = FONT_FAM) +
        theme(
          plot.title         = element_text(face = "bold", size = BASE_SIZE + 1),
          plot.subtitle      = element_text(size = 11, colour = "grey40"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y        = element_text(size = 12),
          plot.margin        = margin(t = 8, r = 8, b = 8, l = 8)
        )

      ## Fixed height: 1.6 in header/footer + 0.55 in per row
      plot_h <- 1.6 + 0.55 * nrow(pd)
      save_plot(p_forest, paste0(file_prefix, "_T6_forest"), w = 9, h = plot_h)

      ## v8: save T6 forest data as CSV + docx
      pd_out <- pd[, c("Variable", "label", "HR", "CI_lo", "CI_hi", "p")]
      pd_out$p_fmt <- sapply(pd_out$p, fmt_p)
      write.csv(pd_out,
        file.path(tables_dir, paste0(file_prefix, "_T6_forest_", output_date, ".csv")),
        row.names = FALSE)
      save_table_docx(pd_out,
        file.path(tables_dir, paste0(file_prefix, "_T6_forest_", output_date, ".docx")),
        title = paste0("Task 6 - Clinical MV Forest | ", cohort_label))
      } ## end nrow check
    }, error = function(e) cat("  Forest plot error:", e$message, "\n"))
  }

  cat(sprintf("  Task 6 complete for %s\n", cohort_label))

  list(
    n = n, ev = ev,
    cont = res_cont,
    dich = res_dich,
    clin_vars = clin_available,
    model_cont = cx_cont,
    model_dich = cx_dich
  )
}

## ── Task 6 runs ────────────────────────────────────────────────────────────

## Discovery cohort clinical variables - using original column names
disc_clin_vars <- c(
  "Age_cat_2"       = "Age",
  "Grade"           = "Grade",
  "Stage_TNM_44"    = "Stage",
  "Syst_treat_NE"   = "Systemic treatment",
  "diadtm_final"    = "Diagnosis period"
)

## Validation cohort clinical variables
val_clin_vars <- c(
  "Age_bin"    = "Age",
  "Grade_enc"  = "Grade",
  "Node_bin"   = "Nodal status"
)

## Discovery S4 (all-comers)
t6_bs <- run_task6(disc_s4, "DDFStime2021", "status_num",
                   disc_clin_vars, "Discovery All-Comers", "DDFS", "Discovery_S4")

## Discovery Luminal (Lum A + Lum B)
t6_bl <- run_task6(disc_merged_list[["Luminal"]], "DDFStime2021", "status_num",
                   disc_clin_vars, "Discovery Luminal", "DDFS", "Discovery_Merged_Luminal")

## Discovery HER2+
t6_bh <- tryCatch(
  run_task6(disc_merged_list[["HER2+"]], "DDFStime2021", "status_num",
            disc_clin_vars, "Discovery HER2+", "DDFS", "Discovery_Merged_HER2pos"),
  error = function(e) { cat("  T6 HER2+ error:", e$message, "\n"); NULL })

## Discovery TNBC
t6_bt <- tryCatch(
  run_task6(disc_merged_list[["TNBC"]], "DDFStime2021", "status_num",
            disc_clin_vars, "Discovery TNBC", "DDFS", "Discovery_Merged_TNBC"),
  error = function(e) { cat("  T6 TNBC error:", e$message, "\n"); NULL })

## Discovery Early Luminal
t6_be <- run_task6(disc_early, "DDFStime2021", "status_num",
                   disc_clin_vars, "Discovery Early Luminal", "DDFS", "Discovery_Early")

## Discovery Young Luminal
t6_by <- run_task6(disc_young, "DDFStime2021", "status_num",
                   disc_clin_vars, "Discovery Young Luminal", "DDFS", "Discovery_Young")

## Validation cohort Task 6 — all 7 subgroups (OS only)
t6_val_inv   <- run_task6(val_inv,   "OS_months", "status_os",
                           val_clin_vars, "Validation cohort All Invasive", "OS", "VAL_AllInv")
t6_val_inv_y <- tryCatch(
  run_task6(val_inv_young, "OS_months", "status_os",
            val_clin_vars, "Validation cohort All Inv Young", "OS", "VAL_AllInv_Young"),
  error = function(e) { cat("  T6 All Inv Young error:", e$message, "\n"); NULL })
t6_val_inv_o <- tryCatch(
  run_task6(val_inv_old,   "OS_months", "status_os",
            val_clin_vars, "Validation cohort All Inv Old", "OS", "VAL_AllInv_Old"),
  error = function(e) { cat("  T6 All Inv Old error:", e$message, "\n"); NULL })
t6_val_erp   <- run_task6(val_erp,   "OS_months", "status_os",
                           val_clin_vars, "Validation cohort ER+", "OS", "VAL_ERpos")
t6_val_ey    <- tryCatch(
  run_task6(val_young,     "OS_months", "status_os",
            val_clin_vars, "Validation cohort ER+ Young", "OS", "VAL_ERpos_Young"),
  error = function(e) { cat("  T6 ER+ Young error:", e$message, "\n"); NULL })
t6_val_erp_N0 <- tryCatch(
  run_task6(val_erp_N0,   "OS_months", "status_os",
            val_clin_vars, "Validation cohort ER+ N0", "OS", "VAL_ERpos_N0"),
  error = function(e) { cat("  T6 ER+ N0 error:", e$message, "\n"); NULL })
t6_val_ey_N0  <- tryCatch(
  run_task6(val_young_N0,  "OS_months", "status_os",
            val_clin_vars, "Validation cohort ER+ Young N0", "OS", "VAL_ERpos_Young_N0"),
  error = function(e) { cat("  T6 ER+ Young N0 error:", e$message, "\n"); NULL })
t6_val_ern   <- tryCatch(
  run_task6(val_ern,       "OS_months", "status_os",
            val_clin_vars, "Validation cohort ER-", "OS", "VAL_ERneg"),
  error = function(e) { cat("  T6 ER-neg error:", e$message, "\n"); NULL })
t6_val_dcis  <- tryCatch(
  run_task6(val_dcis,      "OS_months", "status_os",
            val_clin_vars, "Validation cohort DCIS", "OS", "VAL_DCIS"),
  error = function(e) { cat("  T6 DCIS error:", e$message, "\n"); NULL })

## Task 6 summary table
t6_summary <- do.call(rbind, lapply(
  list(
    list(r = t6_bs, label = "Discovery All-Comers (DDFS)"),
    list(r = t6_bl, label = "Discovery Luminal (DDFS)"),
    list(r = t6_bh, label = "Discovery HER2+ (DDFS)"),
    list(r = t6_bt, label = "Discovery TNBC (DDFS)"),
    list(r = t6_be, label = "Discovery Early Luminal (DDFS)"),
    list(r = t6_by, label = "Discovery Young Luminal (DDFS)"),
    list(r = t6_val_inv,   label = "Validation cohort All Invasive (OS)"),
    list(r = t6_val_inv_y, label = "Validation cohort All Inv Young (OS)"),
    list(r = t6_val_inv_o, label = "Validation cohort All Inv Old (OS)"),
    list(r = t6_val_erp,   label = "Validation cohort ER+ (OS)"),
    list(r = t6_val_ey,    label = "Validation cohort ER+ Young (OS)"),
    list(r = t6_val_erp_N0, label = "Validation cohort ER+ N0 (OS)"),
    list(r = t6_val_ey_N0,  label = "Validation cohort ER+ Young N0 (OS)"),
    list(r = t6_val_ern,   label = "Validation cohort ER- (OS)"),
    list(r = t6_val_dcis,  label = "Validation cohort DCIS (OS)")
  ),
  function(x) {
    r <- x$r
    if (is.null(r) || is.null(r$cont)) {
      return(data.frame(Cohort = x$label, n = NA, Events = NA,
                        MSPI_HR_cont = NA, MSPI_p_cont = NA,
                        MSPI_HR_dich = NA, MSPI_p_dich = NA,
                        stringsAsFactors = FALSE))
    }
    ## Extract MSPI row from results
    mspi_cont <- r$cont[r$cont$Variable == "MSPI_sd", ]
    mspi_dich <- if (!is.null(r$dich)) {
      r$dich[grep("MSPI_high", r$dich$Variable), ]
    } else NULL

    data.frame(
      Cohort = x$label,
      n = r$n, Events = r$ev,
      MSPI_HR_cont = if (nrow(mspi_cont) > 0) mspi_cont$HR[1] else NA,
      MSPI_p_cont  = if (nrow(mspi_cont) > 0) fmt_p(mspi_cont$p[1]) else NA,
      MSPI_HR_dich = if (!is.null(mspi_dich) && nrow(mspi_dich) > 0) mspi_dich$HR[1] else NA,
      MSPI_p_dich  = if (!is.null(mspi_dich) && nrow(mspi_dich) > 0) fmt_p(mspi_dich$p[1]) else NA,
      stringsAsFactors = FALSE
    )
  }
))

write.csv(t6_summary,
  file.path(tables_dir, paste0("MASTER_T6_Clinical_MV_", output_date, ".csv")),
  row.names = FALSE)

save_table_docx(t6_summary,
  file.path(tables_dir, paste0("MASTER_T6_Clinical_MV_", output_date, ".docx")),
  title = "Task 6 - Clinical MV Master Summary")

cat("\nTask 6 MASTER summary written.\n")

## ═══════════════════════════════════════════════════════════════════════════
## MASTER SUMMARY TABLES
## ═══════════════════════════════════════════════════════════════════════════

cat("\n", strrep("═", 70), "\n")
cat("BUILDING MASTER SUMMARY TABLES AND EXCEL WORKBOOK\n")
cat(strrep("═", 70), "\n")

t1_list <- list(t1_by, t1_be, t1_bs,
                t1_val_inv, t1_val_inv_y, t1_val_inv_o,
                t1_val_erp, t1_val_ey,
                t1_val_erp_N0, t1_val_ey_N0,
                t1_val_ern, t1_val_dcis)
t1_labs <- c("Discovery Young Luminal", "Discovery Early Luminal", "Discovery S4",
             "Validation cohort All Invasive(OS)", "Validation cohort All Inv Young(OS)",
             "Validation cohort All Inv Old(OS)",
             "Validation cohort ER+(OS)", "Validation cohort ER+ Young(OS)",
             "Validation cohort ER+ N0(OS)", "Validation cohort ER+ Young N0(OS)",
             "Validation cohort ER-(OS)", "Validation cohort DCIS(OS)")

get_v <- function(R, nm, fld, dg = 3) {
  if (is.null(R)) return(NA)
  r <- R[[nm]]
  if (is.null(r)) return(NA)
  v <- r[[fld]]
  if (is.null(v) || is.na(v)) return(NA)
  round(v, dg)
}

master_t1 <- data.frame(
  Cohort = t1_labs,
  n = sapply(t1_list, function(r) if (is.null(r)) NA else r$n),
  Events = sapply(t1_list, function(r) if (is.null(r)) NA else r$ev),
  Cont_HR = sapply(t1_list, function(r) get_v(r, "cont", "HR")),
  Cont_p = sapply(t1_list, function(r) fmt_p(get_v(r, "cont", "p"))),
  Med_HR = sapply(t1_list, function(r) get_v(r, "med", "HR")),
  Med_p = sapply(t1_list, function(r) fmt_p(get_v(r, "med", "p"))),
  Opt_HR = sapply(t1_list, function(r) get_v(r, "opt", "HR")),
  Opt_p = sapply(t1_list, function(r) fmt_p(get_v(r, "opt", "p"))),
  MV_HR = sapply(t1_list, function(r) get_v(r, "mv", "HR")),
  MV_p = sapply(t1_list, function(r) fmt_p(get_v(r, "mv", "p"))),
  stringsAsFactors = FALSE
)

write.csv(
  master_t1,
  file.path(tables_dir, paste0("MASTER_T1_MSPI_", output_date, ".csv")),
  row.names = FALSE
)

save_table_docx(
  master_t1,
  file.path(tables_dir, paste0("MASTER_T1_MSPI_", output_date, ".docx")),
  title = "Master Table - MSPI Analysis (T1)"
)

## Subtype summary tables
subtype_t1_master <- do.call(rbind, lapply(names(subtype_master), function(st) {
  r <- subtype_master[[st]]$task1
  data.frame(
    Subtype = st,
    n = ifelse(is.null(r), NA, r$n),
    Events = ifelse(is.null(r), NA, r$ev),
    Cont_HR = ifelse(is.null(r) || is.null(r$cont), NA, round(r$cont$HR, 3)),
    Cont_p = ifelse(is.null(r) || is.null(r$cont), NA, fmt_p(r$cont$p)),
    Med_HR = ifelse(is.null(r) || is.null(r$med), NA, round(r$med$HR, 3)),
    Med_p = ifelse(is.null(r) || is.null(r$med), NA, fmt_p(r$med$p)),
    Opt_HR = ifelse(is.null(r) || is.null(r$opt), NA, round(r$opt$HR, 3)),
    Opt_p = ifelse(is.null(r) || is.null(r$opt), NA, fmt_p(r$opt$p)),
    MV_HR = ifelse(is.null(r) || is.null(r$mv), NA, round(r$mv$HR, 3)),
    MV_p = ifelse(is.null(r) || is.null(r$mv), NA, fmt_p(r$mv$p)),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  subtype_t1_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Subtypes_T1_MSPI_", output_date, ".csv")),
  row.names = FALSE
)

save_table_docx(
  subtype_t1_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Subtypes_T1_MSPI_", output_date, ".docx")),
  title = "Discovery Subtypes - MSPI Analysis (T1)"
)

subtype_t1_quart_master <- do.call(rbind, lapply(names(subtype_master), function(st) {
  r <- subtype_master[[st]]$task1
  qc <- if (!is.null(r) && !is.null(r$quart_cox)) r$quart_cox else NULL
  
  out <- data.frame(
    Subtype = st,
    Q2_vs_Q1_HR = NA, Q2_vs_Q1_p = NA,
    Q3_vs_Q1_HR = NA, Q3_vs_Q1_p = NA,
    Q4_vs_Q1_HR = NA, Q4_vs_Q1_p = NA,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(qc) && nrow(qc) > 0) {
    if (any(grepl("Q2", qc$Comparison))) {
      i <- grep("Q2", qc$Comparison)[1]
      out$Q2_vs_Q1_HR <- qc$HR[i]
      out$Q2_vs_Q1_p  <- qc$p[i]
    }
    if (any(grepl("Q3", qc$Comparison))) {
      i <- grep("Q3", qc$Comparison)[1]
      out$Q3_vs_Q1_HR <- qc$HR[i]
      out$Q3_vs_Q1_p  <- qc$p[i]
    }
    if (any(grepl("Q4", qc$Comparison))) {
      i <- grep("Q4", qc$Comparison)[1]
      out$Q4_vs_Q1_HR <- qc$HR[i]
      out$Q4_vs_Q1_p  <- qc$p[i]
    }
  }
  
  out
}))

## Build subtype T2 Q4 master table from task2 results
subtype_t2_master <- do.call(rbind, lapply(names(subtype_master), function(st) {
  r2 <- subtype_master[[st]]$task2
  out <- data.frame(
    Subtype = st,
    Q4_vs_Q1_HR = NA, Q4_vs_Q1_p = NA,
    Q4_vs_nonQ4_HR = NA, Q4_vs_nonQ4_p = NA,
    Q12_vs_Q34_HR = NA, Q12_vs_Q34_p = NA,
    stringsAsFactors = FALSE
  )
  if (!is.null(r2)) {
    if (!is.null(r2$q4_q1)) {
      out$Q4_vs_Q1_HR <- r2$q4_q1$HR[1]
      out$Q4_vs_Q1_p  <- r2$q4_q1$p[1]
    }
    if (!is.null(r2$q4_rest)) {
      out$Q4_vs_nonQ4_HR <- r2$q4_rest$HR[1]
      out$Q4_vs_nonQ4_p  <- r2$q4_rest$p[1]
    }
    if (!is.null(r2$q12_q34)) {
      out$Q12_vs_Q34_HR <- r2$q12_q34$HR[1]
      out$Q12_vs_Q34_p  <- r2$q12_q34$p[1]
    }
  }
  out
}))

write.csv(
  subtype_t2_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Subtypes_T2_Q4_", output_date, ".csv")),
  row.names = FALSE
)

save_table_docx(
  subtype_t2_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Subtypes_T2_Q4_", output_date, ".docx")),
  title = "Discovery Merged Subtypes - Q4 Analysis (T2)"
)

## Individual M1/M2 affinity per subtype master table
subtype_indiv_master <- do.call(rbind, lapply(names(subtype_master), function(st) {
  r <- subtype_master[[st]]
  data.frame(
    Subtype = st,
    n = ifelse(is.null(r$n_ind), NA, r$n_ind),
    Events = ifelse(is.null(r$ev_ind), NA, r$ev_ind),
    M2_UV_HR = ifelse(is.null(r$m2_uv), NA, round(r$m2_uv$HR, 3)),
    M2_UV_CI = ifelse(is.null(r$m2_uv), NA, sprintf("%.3f-%.3f", r$m2_uv$CI_lo, r$m2_uv$CI_hi)),
    M2_UV_p = ifelse(is.null(r$m2_uv), NA, fmt_p(r$m2_uv$p)),
    M2_MV_HR = ifelse(is.null(r$m2_mv), NA, round(r$m2_mv$HR, 3)),
    M2_MV_p = ifelse(is.null(r$m2_mv), NA, fmt_p(r$m2_mv$p)),
    M1_UV_HR = ifelse(is.null(r$m1_uv), NA, round(r$m1_uv$HR, 3)),
    M1_UV_CI = ifelse(is.null(r$m1_uv), NA, sprintf("%.3f-%.3f", r$m1_uv$CI_lo, r$m1_uv$CI_hi)),
    M1_UV_p = ifelse(is.null(r$m1_uv), NA, fmt_p(r$m1_uv$p)),
    M1_MV_HR = ifelse(is.null(r$m1_mv), NA, round(r$m1_mv$HR, 3)),
    M1_MV_p = ifelse(is.null(r$m1_mv), NA, fmt_p(r$m1_mv$p)),
    stringsAsFactors = FALSE
  )
}))

write.csv(subtype_indiv_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Subtypes_Individual_M1M2_", output_date, ".csv")),
  row.names = FALSE)

save_table_docx(subtype_indiv_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Subtypes_Individual_M1M2_", output_date, ".docx")),
  title = "Discovery Subtypes - Individual M1/M2")

## Same for merged subtypes
merged_indiv_master <- do.call(rbind, lapply(names(merged_master), function(mg) {
  r <- merged_master[[mg]]
  data.frame(
    Subtype = mg,
    n = ifelse(is.null(r$n_ind), NA, r$n_ind),
    Events = ifelse(is.null(r$ev_ind), NA, r$ev_ind),
    M2_UV_HR = ifelse(is.null(r$m2_uv), NA, round(r$m2_uv$HR, 3)),
    M2_UV_CI = ifelse(is.null(r$m2_uv), NA, sprintf("%.3f-%.3f", r$m2_uv$CI_lo, r$m2_uv$CI_hi)),
    M2_UV_p = ifelse(is.null(r$m2_uv), NA, fmt_p(r$m2_uv$p)),
    M2_MV_HR = ifelse(is.null(r$m2_mv), NA, round(r$m2_mv$HR, 3)),
    M2_MV_p = ifelse(is.null(r$m2_mv), NA, fmt_p(r$m2_mv$p)),
    M1_UV_HR = ifelse(is.null(r$m1_uv), NA, round(r$m1_uv$HR, 3)),
    M1_UV_CI = ifelse(is.null(r$m1_uv), NA, sprintf("%.3f-%.3f", r$m1_uv$CI_lo, r$m1_uv$CI_hi)),
    M1_UV_p = ifelse(is.null(r$m1_uv), NA, fmt_p(r$m1_uv$p)),
    M1_MV_HR = ifelse(is.null(r$m1_mv), NA, round(r$m1_mv$HR, 3)),
    M1_MV_p = ifelse(is.null(r$m1_mv), NA, fmt_p(r$m1_mv$p)),
    stringsAsFactors = FALSE
  )
}))

write.csv(merged_indiv_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Merged_Individual_M1M2_", output_date, ".csv")),
  row.names = FALSE)

save_table_docx(merged_indiv_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Merged_Individual_M1M2_", output_date, ".docx")),
  title = "Discovery Merged Subtypes - Individual M1/M2")

cat("\nIndividual M1/M2 affinity per-subtype MASTER tables written.\n")

## ---- Merged molecular subtype MASTER tables ----

build_t1_row <- function(mg, master_obj) {
  r <- master_obj[[mg]]$task1
  data.frame(
    Subtype = mg,
    n = ifelse(is.null(r), NA, r$n),
    Events = ifelse(is.null(r), NA, r$ev),
    Cont_HR = ifelse(is.null(r) || is.null(r$cont), NA, round(r$cont$HR, 3)),
    Cont_p = ifelse(is.null(r) || is.null(r$cont), NA, fmt_p(r$cont$p)),
    Med_HR = ifelse(is.null(r) || is.null(r$med), NA, round(r$med$HR, 3)),
    Med_p = ifelse(is.null(r) || is.null(r$med), NA, fmt_p(r$med$p)),
    Opt_HR = ifelse(is.null(r) || is.null(r$opt), NA, round(r$opt$HR, 3)),
    Opt_p = ifelse(is.null(r) || is.null(r$opt), NA, fmt_p(r$opt$p)),
    MV_HR = ifelse(is.null(r) || is.null(r$mv), NA, round(r$mv$HR, 3)),
    MV_p = ifelse(is.null(r) || is.null(r$mv), NA, fmt_p(r$mv$p)),
    stringsAsFactors = FALSE
  )
}

merged_t1_master <- do.call(rbind, lapply(names(merged_master), function(mg) {
  build_t1_row(mg, merged_master)
}))

write.csv(
  merged_t1_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Merged_T1_MSPI_", output_date, ".csv")),
  row.names = FALSE
)

save_table_docx(
  merged_t1_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Merged_T1_MSPI_", output_date, ".docx")),
  title = "Discovery Merged Subtypes - MSPI Analysis (T1)"
)

## Merged subtype T2 Q4 table
merged_t2_master <- do.call(rbind, lapply(names(merged_master), function(mg) {
  r2 <- merged_master[[mg]]$task2
  out <- data.frame(
    Subtype = mg,
    Q4_vs_Q1_HR = NA, Q4_vs_Q1_p = NA,
    Q4_vs_nonQ4_HR = NA, Q4_vs_nonQ4_p = NA,
    Q12_vs_Q34_HR = NA, Q12_vs_Q34_p = NA,
    stringsAsFactors = FALSE
  )
  if (!is.null(r2)) {
    if (!is.null(r2$q4_q1)) {
      out$Q4_vs_Q1_HR <- r2$q4_q1$HR[1]
      out$Q4_vs_Q1_p  <- r2$q4_q1$p[1]
    }
    if (!is.null(r2$q4_rest)) {
      out$Q4_vs_nonQ4_HR <- r2$q4_rest$HR[1]
      out$Q4_vs_nonQ4_p  <- r2$q4_rest$p[1]
    }
    if (!is.null(r2$q12_q34)) {
      out$Q12_vs_Q34_HR <- r2$q12_q34$HR[1]
      out$Q12_vs_Q34_p  <- r2$q12_q34$p[1]
    }
  }
  out
}))

write.csv(
  merged_t2_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Merged_T2_Q4_", output_date, ".csv")),
  row.names = FALSE
)

save_table_docx(
  merged_t2_master,
  file.path(tables_dir, paste0("MASTER_Discovery_Merged_T2_Q4_", output_date, ".docx")),
  title = "Discovery Merged Subtypes - Q4 Master (T2)"
)

cat("\nMerged molecular subtype MASTER tables written.\n")

## Excel workbook
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Master_T1_MSPI")
openxlsx::writeData(wb, "Master_T1_MSPI", master_t1)

openxlsx::addWorksheet(wb, "Discovery_Subtypes_T1")
openxlsx::writeData(wb, "Discovery_Subtypes_T1", subtype_t1_master)

openxlsx::addWorksheet(wb, "Discovery_Merged_T1")
openxlsx::writeData(wb, "Discovery_Merged_T1", merged_t1_master)

openxlsx::addWorksheet(wb, "Discovery_Merged_T2_Q4")
openxlsx::writeData(wb, "Discovery_Merged_T2_Q4", merged_t2_master)

openxlsx::addWorksheet(wb, "Discovery_Subtypes_T1_Quart")
openxlsx::writeData(wb, "Discovery_Subtypes_T1_Quart", subtype_t1_quart_master)

openxlsx::addWorksheet(wb, "Discovery_Subtypes_Indiv_M1M2")
openxlsx::writeData(wb, "Discovery_Subtypes_Indiv_M1M2", subtype_indiv_master)

openxlsx::addWorksheet(wb, "Discovery_Merged_Indiv_M1M2")
openxlsx::writeData(wb, "Discovery_Merged_Indiv_M1M2", merged_indiv_master)

openxlsx::addWorksheet(wb, "T6_Clinical_MV_Summary")
openxlsx::writeData(wb, "T6_Clinical_MV_Summary", t6_summary)

csv_files <- list.files(tables_dir, pattern = "\\.csv$", full.names = TRUE)
for (f in csv_files) {
  sn <- substr(tools::file_path_sans_ext(basename(f)), 1, 31)
  tryCatch({
    if (!(sn %in% names(wb))) {
      openxlsx::addWorksheet(wb, sn)
      openxlsx::writeData(wb, sn, read.csv(f, stringsAsFactors = FALSE))
    }
  }, error = function(e) NULL)
}

openxlsx::saveWorkbook(
  wb,
  file.path(tables_dir, paste0("Integrated_Analysis_", output_date, ".xlsx")),
  overwrite = TRUE
)

## ═══════════════════════════════════════════════════════════════════════════
## COMPREHENSIVE TEXT SUMMARY
## ═══════════════════════════════════════════════════════════════════════════

smry <- file.path(output_dir, paste0("Integrated_Analysis_Summary_", output_date, ".txt"))
sink(smry)

cat("================================================================\n")
cat("INTEGRATED SPATIAL IMMUNE ANALYSIS — Discovery + Validation ER+\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("SCORE CONVENTION:\n")
cat("  Higher score = higher affinity (closer macrophage-T-cell proximity)\n")
cat("  MSPI = M2_score - M1_score\n")
cat("  High MSPI = suppression dominant = worse\n\n")

cat("MAIN COHORTS — TASK 1 MSPI SUMMARY:\n")
print(master_t1)

cat("
Discovery INDIVIDUAL SUBTYPES — TASK 1 SUMMARY:\n")
print(subtype_t1_master)

cat("
Discovery MERGED MOLECULAR SUBTYPES (Luminal, HER2+, TNBC) — TASK 1 SUMMARY:\n")
print(merged_t1_master)

cat("
Discovery MERGED MOLECULAR SUBTYPES — Q4 SUMMARY:\n")
print(merged_t2_master)

cat("
Discovery INDIVIDUAL SUBTYPES — TASK 1 MSPI QUARTILE SUMMARY:\n")
print(subtype_t1_quart_master)

## (Validation cohort ER-neg and DCIS results now included in MASTER T1 table above)

cat("\nINDIVIDUAL M1/M2 BY SUBTYPE:\n")
print(subtype_indiv_master)
cat("\nINDIVIDUAL M1/M2 BY MERGED SUBTYPE:\n")
print(merged_indiv_master)

cat("\nTASK 6 — CLINICAL MV SUMMARY:\n")
print(t6_summary)

cat("\nTASK 4 — POOLED IPD:\n")
if (!is.null(t4_A)) {
  cat(sprintf(
    "  Pool A  MSPI HR=%.3f p=%s | M1 HR=%.3f p=%s | M2 HR=%.3f p=%s\n",
    t4_A$HR_mspi, fmt_p(t4_A$p_mspi), t4_A$HR_M1, fmt_p(t4_A$p_M1),
    t4_A$HR_M2, fmt_p(t4_A$p_M2)
  ))
}
if (!is.null(t4_B)) {
  cat(sprintf(
    "  Pool B  MSPI HR=%.3f p=%s | M1 HR=%.3f p=%s | M2 HR=%.3f p=%s\n",
    t4_B$HR_mspi, fmt_p(t4_B$p_mspi), t4_B$HR_M1, fmt_p(t4_B$p_M1),
    t4_B$HR_M2, fmt_p(t4_B$p_M2)
  ))
}
if (!is.null(t4_A_N0)) {
  cat(sprintf(
    "  Pool A N0  MSPI HR=%.3f p=%s | M1 HR=%.3f p=%s | M2 HR=%.3f p=%s\n",
    t4_A_N0$HR_mspi, fmt_p(t4_A_N0$p_mspi), t4_A_N0$HR_M1, fmt_p(t4_A_N0$p_M1),
    t4_A_N0$HR_M2, fmt_p(t4_A_N0$p_M2)
  ))
}
if (!is.null(t4_B_N0)) {
  cat(sprintf(
    "  Pool B N0  MSPI HR=%.3f p=%s | M1 HR=%.3f p=%s | M2 HR=%.3f p=%s\n",
    t4_B_N0$HR_mspi, fmt_p(t4_B_N0$p_mspi), t4_B_N0$HR_M1, fmt_p(t4_B_N0$p_M1),
    t4_B_N0$HR_M2, fmt_p(t4_B_N0$p_M2)
  ))
}

cat("\nTASK 5 — INCREMENTAL VALUE:\n")
for (nm in c("Discovery Young", "Discovery Early", "Discovery S4",
              "Val All Inv", "Val All Inv Young", "Val All Inv Old",
              "Val ER+", "Val ER+ Young",
              "Val ER+ N0", "Val ER+ Young N0",
              "Val ER-", "Val DCIS")) {
  r <- switch(
    nm,
    "Discovery Young"        = t5_by,
    "Discovery Early"        = t5_be,
    "Discovery S4"           = t5_bs,
    "Val All Inv"      = t5_val_inv,
    "Val All Inv Young"= t5_val_inv_y,
    "Val All Inv Old"  = t5_val_inv_o,
    "Val ER+"          = t5_val_erp,
    "Val ER+ Young"    = t5_val_ey,
    "Val ER+ N0"       = t5_val_erp_N0,
    "Val ER+ Young N0" = t5_val_ey_N0,
    "Val ER-"          = t5_val_ern,
    "Val DCIS"         = t5_val_dcis
  )
  if (!is.null(r)) {
    cat(sprintf(
      "  %-20s C_dens=%.3f  C_aff=%.3f  C_full=%.3f  ΔC=%.4f  LRT p=%s\n",
      nm, r$c_dens, r$c_aff, r$c_full, r$delta_c, fmt_p(r$lrt_p)
    ))
  }
}

sink()

cat("\n", strrep("=", 70), "\n")
cat("INTEGRATED ANALYSIS — COMPLETE\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output folder:", output_dir, "\n")
cat(strrep("=", 70), "\n")

cat(paste0("Completed: ", format(Sys.time()), "\n"), file = mlog, append = TRUE)