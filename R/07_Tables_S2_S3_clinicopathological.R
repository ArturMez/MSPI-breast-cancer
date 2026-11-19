########################
## Supplementary Tables S2-S3 (clinicopathological characteristics of the MSPI analysis sets)
##
## Run from the repository root. Inputs (original names; the repository names in brackets also work):
##   ./data/discovery/   BCY_workingr_832_updated_211024.sav [cohort_clinical.sav, clinical_data.sav]
##                       discovery_distances_mean_M1_to_CD3_nucleus.txt, discovery_distances_mean_M2_to_CD3_nucleus.txt
##   ./data/validation/  NK_BRCA_for_affinity.txt [clinical.txt]
##                       BREAST_ENCIMA_affinity_CD3e_M1M2_myeloidpool_50um.txt [affinity_scores.txt]
##                       Huvuddatabas.xlsx [clinical_extra.xlsx]
## Outputs: ./output/tables_S2_S3/ (the discovery numbers of Figure 1a are counted by script 10)
## In the code, bcy_* objects hold the discovery cohort and enc_* objects the validation cohort.
##
## Version 2 (nucleus-supported T cells): reads the discovery affinity scores of
## python/01_discovery_affinity_nucleus.py
## (discovery_distances_mean_M1/M2_to_CD3_nucleus.txt in
## ./data/discovery/). The input names of version 1 are not accepted, so its inputs cannot be used by mistake.
## Everything else is unchanged from version 1.
########################


## ── Paths ──────────────────────────────────────────────────────────────────
proj_root  <- "./data"
output_dir <- "./output/tables_S2_S3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Data dir:", proj_root, "\n")
cat("Output dir:", output_dir, "\n")

## Input files may carry their original names or the names used in this repository;
## the first file found is used.
find_input <- function(dir, ...) {
  for (f in c(...)) {
    p <- file.path(dir, f)
    if (file.exists(p)) return(p)
  }
  stop("Cannot find any of: ", paste(c(...), collapse = ", "), " in ", dir)
}

## ════════════════════════════════════════════════════════════════════════════
## CLINICOPATHOLOGICAL TABLES FOR MSPI-AVAILABLE SUBSETS
## Produces .docx tables (styled via flextable + officer)
## for BCY (discovery) and ENCIMA (validation)
## ════════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════\n")
cat("Clinicopathological Tables\n")
cat("═══════════════════════════════════════════════\n\n")

## ── Libraries for data loading and table generation ────────────────────────
pkgs2 <- c("foreign", "data.table", "dplyr", "openxlsx", "flextable", "officer")
for (p in pkgs2) {
  if (!suppressWarnings(require(p, character.only = TRUE, quietly = TRUE))) {
    install.packages(p, repos = "https://cran.rstudio.com/")
    library(p, character.only = TRUE)
  }
}

## ── Shared flextable styling ───────────────────────────────────────────────
style_ft <- function(ft) {
  bdr <- officer::fp_border(color = "black", width = 0.5)
  ft %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 10, part = "body") %>%
    flextable::fontsize(size = 10, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "#D5E8F0", part = "header") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::align(align = "left", j = 1, part = "body") %>%
    flextable::align(align = "center", j = 2, part = "body") %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = bdr, part = "header") %>%
    flextable::hline_bottom(border = bdr, part = "header") %>%
    flextable::hline_bottom(border = bdr, part = "body") %>%
    flextable::set_table_properties(width = 1, layout = "autofit") %>%
    flextable::padding(padding.top = 2, padding.bottom = 2, part = "all")
}

## helper: format count (%)
fmt_n <- function(n, total) {
  sprintf("%d (%.1f)", n, 100 * n / total)
}

save_clinpath_docx <- function(tab_df, out_path, title) {
  ft <- flextable::flextable(tab_df) %>% style_ft()
  ## Bold the category header rows (where Value column is empty)
  bold_rows <- which(tab_df[[2]] == "")
  if (length(bold_rows) > 0) {
    ft <- ft %>% flextable::bold(i = bold_rows, j = 1, part = "body")
  }
  ## Indent sub-rows (where Value column is non-empty)
  indent_rows <- which(tab_df[[2]] != "")
  if (length(indent_rows) > 0) {
    ft <- ft %>%
      flextable::padding(i = indent_rows, j = 1, padding.left = 15, part = "body")
  }
  doc <- officer::read_docx()
  doc <- doc %>%
    officer::body_add_par(title, style = "heading 1") %>%
    officer::body_add_par("", style = "Normal") %>%
    flextable::body_add_flextable(ft)
  print(doc, target = out_path)
  cat("  Saved:", out_path, "\n")
}


## ══════════════════════════════════════════════════════════════════════════
## BCY (Discovery Cohort)
## ══════════════════════════════════════════════════════════════════════════

cat("Loading discovery cohort data...\n")

bcy_base <- file.path(proj_root, "discovery")

bcy_spss_candidates <- c(
  file.path(bcy_base, "BCY_workingr_832_updated_211024.sav"),
  file.path(bcy_base, "Matching to clinical data", "BCY_workingr_832_updated_211024.sav"),
  file.path(bcy_base, "cohort_clinical.sav"),
  file.path(bcy_base, "clinical_data.sav")
)
bcy_spss <- NA
for (p in bcy_spss_candidates) {
  if (file.exists(p)) { bcy_spss <- p; break }
}
if (is.na(bcy_spss)) stop("Cannot find the discovery cohort SPSS file")

bcy_raw <- foreign::read.spss(
  bcy_spss, use.value.labels = TRUE, to.data.frame = TRUE,
  max.value.labels = Inf, trim.factor.names = TRUE,
  trim_values = TRUE, reencode = NA
)

## Load affinity distances
m1_raw <- data.table::fread(find_input(bcy_base, "discovery_distances_mean_M1_to_CD3_nucleus.txt"),
                            data.table = FALSE)
m2_raw <- data.table::fread(find_input(bcy_base, "discovery_distances_mean_M2_to_CD3_nucleus.txt"),
                            data.table = FALSE)
colnames(m1_raw)[colnames(m1_raw) == "delta_CD3"] <- "M1_to_CD3_distance"
colnames(m2_raw)[colnames(m2_raw) == "delta_CD3"] <- "M2_to_CD3_distance"

## Merge
bcy <- bcy_raw %>%
  dplyr::left_join(m1_raw[, c("AnnUnit", "M1_to_CD3_distance")], by = "AnnUnit") %>%
  dplyr::left_join(m2_raw[, c("AnnUnit", "M2_to_CD3_distance")], by = "AnnUnit")

## Remove unclassified
bcy <- bcy[as.character(bcy$Subtype_BCY21) != "Unclassified", ]
bcy$Subtype_BCY21 <- droplevels(bcy$Subtype_BCY21)

## Define subsets
bcy$has_M1 <- !is.na(bcy$M1_to_CD3_distance)
bcy$has_M2 <- !is.na(bcy$M2_to_CD3_distance)
bcy$has_MSPI <- bcy$has_M1 & bcy$has_M2

## Survival filter (valid DDFS)
bcy$status_num <- ifelse(
  as.character(bcy$DDFScens2021) == "Distant recurrence" |
    as.character(bcy$BCcens2021) == "Breast cancer death",
  1L, 0L
)
bcy$DDFStime2021 <- as.numeric(bcy$DDFStime2021)
bcy$valid_surv <- !is.na(bcy$DDFStime2021) & bcy$DDFStime2021 > 0

cat(sprintf("  Discovery total (excl. unclassified): n=%d\n", nrow(bcy)))
cat(sprintf("  M1 available: n=%d\n", sum(bcy$has_M1)))
cat(sprintf("  M2 available: n=%d\n", sum(bcy$has_M2)))
cat(sprintf("  MSPI available: n=%d\n", sum(bcy$has_MSPI)))
cat(sprintf("  MSPI + valid DDFS: n=%d\n", sum(bcy$has_MSPI & bcy$valid_surv)))

## MSPI analysis subset
bcy_mspi <- bcy %>% dplyr::filter(has_MSPI, valid_surv)
N_bcy <- nrow(bcy_mspi)

cat(sprintf("  Discovery MSPI analysis set: n=%d, events=%d\n", N_bcy, sum(bcy_mspi$status_num)))

## Build table rows
bcy_age_num <- suppressWarnings(as.numeric(as.character(bcy_mspi$age)))
if (all(is.na(bcy_age_num))) {
  ## age might be stored as a labelled variable - try Age_cat_2
  bcy_age_num <- NULL
}

bcy_rows <- data.frame(
  Characteristics = character(),
  Value = character(),
  stringsAsFactors = FALSE
)

## Age
bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Age", Value = ""))
if (!is.null(bcy_age_num) && !all(is.na(bcy_age_num))) {
  bcy_rows <- rbind(bcy_rows,
    data.frame(Characteristics = "Mean age, years",
               Value = sprintf("%.1f", mean(bcy_age_num, na.rm = TRUE))))
}
## FIX: SPSS Age_cat_2 labels are "<40" and "\u226540" (unicode \u2265), NOT ">=40".
## Previous code matched ">=40" which silently returned 0 -> table showed
## "Age \u226540 = 0 (0.0%)". Match against the actual unicode character.
age_cat <- as.character(bcy_mspi$Age_cat_2)
n_young <- sum(age_cat == "<40", na.rm = TRUE)
n_old   <- sum(age_cat == "\u226540", na.rm = TRUE)
bcy_rows <- rbind(bcy_rows,
  data.frame(Characteristics = "<40", Value = fmt_n(n_young, N_bcy)),
  data.frame(Characteristics = "\u226540", Value = fmt_n(n_old, N_bcy))
)
stopifnot(n_young + n_old == N_bcy)  # safeguard: catch future mis-coding

## Stage
bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Stage", Value = ""))
stage <- as.character(bcy_mspi$Stage_TNM_44)
for (s in c("Stage I", "Stage IIa", "Stage IIb")) {
  n_s <- sum(stage == s, na.rm = TRUE)
  bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = s, Value = fmt_n(n_s, N_bcy)))
}
n_iii_iv <- sum(stage %in% c("Stage IIIa", "Stage IIIb", "Stage IIIc", "Stage IV",
                              "Stage III", "Stage III-IV"), na.rm = TRUE)
bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Stage III-IV", Value = fmt_n(n_iii_iv, N_bcy)))
## FIX: Stage_TNM_44 SPSS label for the missing level is "Stage missing"
## (level 5), not "missing". Previous code checked stage=="missing" and
## therefore never counted the 1 missing case.
n_stage_miss <- sum(is.na(stage) | stage == "" |
                    stage == "missing" | stage == "Stage missing")
if (n_stage_miss > 0) {
  bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Stage missing", Value = fmt_n(n_stage_miss, N_bcy)))
}

## Grade
## FIX: The BCY SAV does not have a "Grade" column. Histological grade is
## stored in `Elston_final`, with SPSS labels "Elston I", "Elston II",
## "Elston III", and "Missing". Previous code referenced bcy_mspi$Grade,
## which evaluated to NULL and produced an empty Grade row.
bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Grade", Value = ""))
grade_col <- if ("Elston_final" %in% names(bcy_mspi)) {
  "Elston_final"
} else if ("Elston_final_NE" %in% names(bcy_mspi)) {
  "Elston_final_NE"
} else if ("Grade" %in% names(bcy_mspi)) {
  "Grade"
} else {
  NA_character_
}
if (!is.na(grade_col)) {
  grade <- as.character(bcy_mspi[[grade_col]])
  ## Map the various possible spellings to "Grade 1/2/3"
  for (lbl_pair in list(
        c("Grade 1", c("Elston I", "I", "Grade I", "1")),
        c("Grade 2", c("Elston II", "II", "Grade II", "2")),
        c("Grade 3", c("Elston III", "III", "Grade III", "3")))) {
    out_lbl <- lbl_pair[[1]]
    matches <- lbl_pair[[2]]
    n_g <- sum(grade %in% matches, na.rm = TRUE)
    if (n_g > 0) bcy_rows <- rbind(bcy_rows,
      data.frame(Characteristics = out_lbl, Value = fmt_n(n_g, N_bcy)))
  }
  ## Anything not matched above is treated as missing
  matched_set <- c("Elston I", "I", "Grade I", "1",
                   "Elston II", "II", "Grade II", "2",
                   "Elston III", "III", "Grade III", "3")
  n_grade_miss <- sum(is.na(grade) | grade == "" | !(grade %in% matched_set), na.rm = TRUE)
  if (n_grade_miss > 0) {
    bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Grade missing", Value = fmt_n(n_grade_miss, N_bcy)))
  }
} else {
  bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Grade missing", Value = fmt_n(N_bcy, N_bcy)))
  warning("Discovery cohort: no grade column found (Elston_final / Grade) - Table S2 grade row will be all-missing")
}

## Molecular Subtype
bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Molecular subtype", Value = ""))
subtype <- as.character(bcy_mspi$Subtype_BCY21)
for (st in c("Lum A", "Lum B", "Lum Her2", "Her2 exp", "Triple neg")) {
  n_st <- sum(subtype == st, na.rm = TRUE)
  bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = st, Value = fmt_n(n_st, N_bcy)))
}

## Systemic treatment
bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Systemic treatment", Value = ""))
treat <- as.character(bcy_mspi$Syst_treat_NE)
treat_levels <- sort(unique(treat[!is.na(treat) & treat != ""]))
for (tl in treat_levels) {
  n_tl <- sum(treat == tl, na.rm = TRUE)
  bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = tl, Value = fmt_n(n_tl, N_bcy)))
}
n_treat_miss <- sum(is.na(treat) | treat == "")
if (n_treat_miss > 0) {
  bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Missing", Value = fmt_n(n_treat_miss, N_bcy)))
}

## Diagnosis period
bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = "Diagnosis period", Value = ""))
diag_p <- as.character(bcy_mspi$diadtm_final)
diag_levels <- sort(unique(diag_p[!is.na(diag_p) & diag_p != ""]))
for (dl in diag_levels) {
  n_dl <- sum(diag_p == dl, na.rm = TRUE)
  bcy_rows <- rbind(bcy_rows, data.frame(Characteristics = dl, Value = fmt_n(n_dl, N_bcy)))
}

## Set column name
colnames(bcy_rows) <- c("Characteristics", sprintf("Cases (percentage)\nN = %d", N_bcy))

## Save
bcy_table_path <- file.path(output_dir, "Table_S2_discovery_clinicopathological.docx")
save_clinpath_docx(bcy_rows, bcy_table_path,
  sprintf("Table S2. Distribution of clinicopathological parameters in tumors within the discovery cohort (MSPI-available subset, N = %d).", N_bcy))


## ══════════════════════════════════════════════════════════════════════════
## ENCIMA (Validation Cohort)
## ══════════════════════════════════════════════════════════════════════════

cat("\nLoading validation cohort data...\n")

enc_base <- file.path(proj_root, "validation")

enc_clin <- data.table::fread(find_input(enc_base, "NK_BRCA_for_affinity.txt", "clinical.txt"),
                              data.table = FALSE)
enc_aff <- data.table::fread(
  find_input(enc_base, "BREAST_ENCIMA_affinity_CD3e_M1M2_myeloidpool_50um.txt",
             "BREAST_validation_affinity_CD3e_M1M2_myeloidpool_50um.txt", "affinity_scores.txt"),
  data.table = FALSE
)
enc <- dplyr::left_join(enc_clin, enc_aff, by = "Case")

## Load Huvuddatabas for extra clinical variables
enc_hd <- openxlsx::read.xlsx(find_input(enc_base, "Huvuddatabas.xlsx", "clinical_extra.xlsx"))
## Rename PAD column to avoid issues with special characters
## openxlsx may read it as "PAD+AR" or "PAD.AR" depending on version
pad_col <- grep("^PAD", names(enc_hd), value = TRUE)[1]
if (!is.na(pad_col)) {
  names(enc_hd)[names(enc_hd) == pad_col] <- "PAD_AR"
  cat("  Huvuddatabas PAD column found as:", pad_col, "-> renamed to PAD_AR\n")
}
## Similarly fix DIAG YEAR column
diag_col <- grep("^DIAG.*YEAR", names(enc_hd), value = TRUE)[1]
if (!is.na(diag_col)) {
  names(enc_hd)[names(enc_hd) == diag_col] <- "DIAG_YEAR"
  cat("  Huvuddatabas DIAG YEAR column found as:", diag_col, "-> renamed to DIAG_YEAR\n")
}

## Extract case number for merge with Huvuddatabas
enc$case_num <- as.numeric(sub(".*# (\\d+)_.*", "\\1", enc$Case))
## FIX: Huvuddatabas has 124 duplicated PAD+AR rows (bilateral / metachronous tumors
## recorded as separate rows). A naive left_join inflates enc by ~9 rows and creates
## 7 duplicated ER+ patients. Deduplicate Huvuddatabas on PAD_AR before the join.
n_hd_before <- nrow(enc_hd)
enc_hd <- enc_hd %>% dplyr::distinct(PAD_AR, .keep_all = TRUE)
cat(sprintf("  Huvuddatabas dedup on PAD_AR: %d -> %d rows\n",
            n_hd_before, nrow(enc_hd)))
n_enc_before <- nrow(enc)
enc <- dplyr::left_join(enc, enc_hd, by = c("case_num" = "PAD_AR"))
stopifnot(nrow(enc) == n_enc_before)  # safeguard against future one-to-many duplication
cat(sprintf("  Huvuddatabas matched: %d / %d cases\n",
            sum(!is.na(enc$AGEYEARS)), nrow(enc)))

## Define MSPI availability
enc$has_M1 <- !is.na(enc$M1_to_CD3e)
enc$has_M2 <- !is.na(enc$M2_to_CD3e)
enc$has_MSPI <- enc$has_M1 & enc$has_M2

## Survival
enc$OS_days <- suppressWarnings(as.numeric(enc$Time_Diagnosis_Last_followup..OS.))
enc$status_os <- ifelse(enc$Event_last_followup == "Dead", 1L, 0L)
enc$valid_surv <- !is.na(enc$OS_days) & enc$OS_days > 0

cat(sprintf("  Validation total: n=%d\n", nrow(enc)))
cat(sprintf("  M1 available: n=%d\n", sum(enc$has_M1)))
cat(sprintf("  M2 available: n=%d\n", sum(enc$has_M2)))
cat(sprintf("  MSPI available: n=%d\n", sum(enc$has_MSPI)))

## ----- All invasive MSPI subset -----
## FIX: Drop "BRCA" (unspecified) tumor type to match the analysis pipeline
## (Script_Integrated_MSPI_Analysis_v8) and Figure 1 patient flow. The Tumor_type_code
## "BRCA" is the unspecified-subtype label (n=9, 6 with both M1+M2 affinity); it was
## previously lumped into ER-negative, which inflated the count by 6.
enc_inv <- enc %>%
  dplyr::filter(
    Tumor_type_code %in% c("BRCA_ER_pos", "BRCA_ER_neg"),
    has_MSPI, valid_surv
  )
N_enc_inv <- nrow(enc_inv)
cat(sprintf("  Validation all invasive MSPI + valid OS: n=%d, events=%d\n",
            N_enc_inv, sum(enc_inv$status_os)))

## Build table for ALL invasive ENCIMA
enc_rows <- data.frame(Characteristics = character(), Value = character(),
                       stringsAsFactors = FALSE)

## Age
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Age", Value = ""))
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Mean age, years",
             Value = sprintf("%.1f", mean(enc_inv$AGEYEARS, na.rm = TRUE))))
age_cat <- as.character(enc_inv$Age)
for (ac in c("<=50", ">50-<65", ">=65")) {
  n_ac <- sum(age_cat == ac, na.rm = TRUE)
  enc_rows <- rbind(enc_rows, data.frame(Characteristics = ac, Value = fmt_n(n_ac, N_enc_inv)))
}

## ER status
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "ER status", Value = ""))
n_erp <- sum(enc_inv$Tumor_type_code == "BRCA_ER_pos", na.rm = TRUE)
n_ern <- sum(enc_inv$Tumor_type_code == "BRCA_ER_neg", na.rm = TRUE)
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "ER-positive", Value = fmt_n(n_erp, N_enc_inv)),
  data.frame(Characteristics = "ER-negative", Value = fmt_n(n_ern, N_enc_inv))
)

## Grade (from Diff_grade in NK_BRCA, more complete than GRADE in Huvuddatabas)
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Grade", Value = ""))
dg <- enc_inv$Diff_grade
for (g in c("1", "2", "3")) {
  n_g <- sum(dg == g, na.rm = TRUE)
  if (n_g > 0) enc_rows <- rbind(enc_rows,
    data.frame(Characteristics = paste0("Grade ", g), Value = fmt_n(n_g, N_enc_inv)))
}
n_g_miss <- sum(dg == "missing" | is.na(dg), na.rm = TRUE)
if (n_g_miss > 0) {
  enc_rows <- rbind(enc_rows,
    data.frame(Characteristics = "Grade missing", Value = fmt_n(n_g_miss, N_enc_inv)))
}

## Nodal status (from NK_BRCA pN_stage)
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Nodal status", Value = ""))
pn <- enc_inv$pN_stage
n_n0 <- sum(pn == "NO axil metastasis", na.rm = TRUE)
n_n1 <- sum(pn == "axil metastasis", na.rm = TRUE)
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Node-negative (N0)", Value = fmt_n(n_n0, N_enc_inv)),
  data.frame(Characteristics = "Node-positive", Value = fmt_n(n_n1, N_enc_inv))
)

## Tumor size (from Huvuddatabas SIZE, mm)
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Tumor size", Value = ""))
sz <- enc_inv$SIZE
sz_clean <- sz[!is.na(sz) & sz < 900]  # exclude 999 (multifocal) and 88 (missing)
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Mean size, mm",
             Value = sprintf("%.1f", mean(sz_clean, na.rm = TRUE))))
n_le10 <- sum(sz_clean <= 10, na.rm = TRUE)
n_11_20 <- sum(sz_clean > 10 & sz_clean <= 20, na.rm = TRUE)
n_gt20 <- sum(sz_clean > 20, na.rm = TRUE)
n_sz_miss <- sum(is.na(sz) | sz >= 900)  # 999 = multifocal, 88 = missing
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "\u226410 mm", Value = fmt_n(n_le10, N_enc_inv)),
  data.frame(Characteristics = "11-20 mm", Value = fmt_n(n_11_20, N_enc_inv)),
  data.frame(Characteristics = ">20 mm", Value = fmt_n(n_gt20, N_enc_inv))
)
if (n_sz_miss > 0) {
  enc_rows <- rbind(enc_rows,
    data.frame(Characteristics = "Size missing/multifocal", Value = fmt_n(n_sz_miss, N_enc_inv)))
}

## Surgery type (from Huvuddatabas TREAT: 1=mastectomy, 2=sector, 4=no surgery)
## NOTE - a small number of records carry an undocumented TREAT = 3 (and some
## have NA). They were previously dropped silently, so Mastectomy + Breast-
## conserving did not sum to N. We now report them under "Unknown/other" so
## the column adds up to N.
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Surgery", Value = ""))
tr <- enc_inv$TREAT
n_mast   <- sum(tr == 1, na.rm = TRUE)
n_sect   <- sum(tr == 2, na.rm = TRUE)
n_nosurg <- sum(tr == 4, na.rm = TRUE)
n_other  <- sum(is.na(tr) | !(tr %in% c(1, 2, 4)))
if (n_mast > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Mastectomy", Value = fmt_n(n_mast, N_enc_inv)))
if (n_sect > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Breast-conserving", Value = fmt_n(n_sect, N_enc_inv)))
if (n_nosurg > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "No surgery", Value = fmt_n(n_nosurg, N_enc_inv)))
if (n_other > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Unknown/other", Value = fmt_n(n_other, N_enc_inv)))

## Radiotherapy (from Huvuddatabas RT: 1=yes, 2=no)
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Radiotherapy", Value = ""))
rt <- enc_inv$RT
n_rt_yes <- sum(rt == 1, na.rm = TRUE)
n_rt_no  <- sum(rt == 2, na.rm = TRUE)
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_rt_yes, N_enc_inv)),
  data.frame(Characteristics = "No", Value = fmt_n(n_rt_no, N_enc_inv))
)

## Chemotherapy (from Huvuddatabas CHEMOTHERAPY: 1=yes, 2=no)
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Chemotherapy", Value = ""))
chemo <- enc_inv$CHEMOTHERAPY
n_ch_yes <- sum(chemo == 1, na.rm = TRUE)
n_ch_no  <- sum(chemo == 2, na.rm = TRUE)
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_ch_yes, N_enc_inv)),
  data.frame(Characteristics = "No", Value = fmt_n(n_ch_no, N_enc_inv))
)

## Endocrine therapy (from Huvuddatabas ENDOCRINETHERAPY: 1=yes, 2=no, 4=?)
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Endocrine therapy", Value = ""))
endo <- enc_inv$ENDOCRINETHERAPY
n_en_yes <- sum(endo == 1, na.rm = TRUE)
n_en_no  <- sum(endo == 2, na.rm = TRUE)
n_en_oth <- sum(!is.na(endo) & !(endo %in% c(1, 2)))
enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_en_yes, N_enc_inv)),
  data.frame(Characteristics = "No", Value = fmt_n(n_en_no, N_enc_inv))
)
if (n_en_oth > 0) {
  enc_rows <- rbind(enc_rows,
    data.frame(Characteristics = "Unknown", Value = fmt_n(n_en_oth, N_enc_inv)))
}

## Diagnosis period
enc_rows <- rbind(enc_rows, data.frame(Characteristics = "Diagnosis period", Value = ""))
dy <- enc_inv$DIAG_YEAR
## Convert 2-digit years: 85-99 -> 1985-1999, 2002-2004 stay as-is
dy_full <- ifelse(!is.na(dy) & dy < 100, dy + 1900, dy)
n_85_89 <- sum(dy_full >= 1985 & dy_full <= 1989, na.rm = TRUE)
n_90_94 <- sum(dy_full >= 1990 & dy_full <= 1994, na.rm = TRUE)
n_95_99 <- sum(dy_full >= 1995 & dy_full <= 1999, na.rm = TRUE)
n_00_04 <- sum(dy_full >= 2000 & dy_full <= 2004, na.rm = TRUE)
if (n_85_89 > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "1985-1989", Value = fmt_n(n_85_89, N_enc_inv)))
if (n_90_94 > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "1990-1994", Value = fmt_n(n_90_94, N_enc_inv)))
if (n_95_99 > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "1995-1999", Value = fmt_n(n_95_99, N_enc_inv)))
if (n_00_04 > 0) enc_rows <- rbind(enc_rows,
  data.frame(Characteristics = "2000-2004", Value = fmt_n(n_00_04, N_enc_inv)))

## Set column name
colnames(enc_rows) <- c("Characteristics", sprintf("Cases (percentage)\nN = %d", N_enc_inv))

## Save
enc_table_path <- file.path(output_dir, "Table_S3_validation_clinicopathological.docx")
save_clinpath_docx(enc_rows, enc_table_path,
  sprintf("Table S3. Distribution of clinicopathological parameters in tumors within the validation cohort (MSPI-available subset, all invasive, N = %d).", N_enc_inv))


## ----- ER-positive MSPI subset (primary analysis population) -----
enc_erp <- enc %>%
  dplyr::filter(
    Tumor_type_code == "BRCA_ER_pos",
    has_MSPI, valid_surv
  )
N_enc_erp <- nrow(enc_erp)
cat(sprintf("  Validation ER+ MSPI + valid OS: n=%d, events=%d\n",
            N_enc_erp, sum(enc_erp$status_os)))

## Build ER+ table (similar structure, minus ER status row)
erp_rows <- data.frame(Characteristics = character(), Value = character(),
                       stringsAsFactors = FALSE)

## Age
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Age", Value = ""))
erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Mean age, years",
             Value = sprintf("%.1f", mean(enc_erp$AGEYEARS, na.rm = TRUE))))
age_cat <- as.character(enc_erp$Age)
for (ac in c("<=50", ">50-<65", ">=65")) {
  n_ac <- sum(age_cat == ac, na.rm = TRUE)
  erp_rows <- rbind(erp_rows, data.frame(Characteristics = ac, Value = fmt_n(n_ac, N_enc_erp)))
}

## Grade
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Grade", Value = ""))
dg <- enc_erp$Diff_grade
for (g in c("1", "2", "3")) {
  n_g <- sum(dg == g, na.rm = TRUE)
  if (n_g > 0) erp_rows <- rbind(erp_rows,
    data.frame(Characteristics = paste0("Grade ", g), Value = fmt_n(n_g, N_enc_erp)))
}
n_g_miss <- sum(dg == "missing" | is.na(dg), na.rm = TRUE)
if (n_g_miss > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Grade missing", Value = fmt_n(n_g_miss, N_enc_erp)))

## Nodal status
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Nodal status", Value = ""))
pn <- enc_erp$pN_stage
n_n0 <- sum(pn == "NO axil metastasis", na.rm = TRUE)
n_n1 <- sum(pn == "axil metastasis", na.rm = TRUE)
erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Node-negative (N0)", Value = fmt_n(n_n0, N_enc_erp)),
  data.frame(Characteristics = "Node-positive", Value = fmt_n(n_n1, N_enc_erp))
)

## Tumor size
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Tumor size", Value = ""))
sz <- enc_erp$SIZE
sz_clean <- sz[!is.na(sz) & sz < 900]
erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Mean size, mm",
             Value = sprintf("%.1f", mean(sz_clean, na.rm = TRUE))))
n_le10 <- sum(sz_clean <= 10, na.rm = TRUE)
n_11_20 <- sum(sz_clean > 10 & sz_clean <= 20, na.rm = TRUE)
n_gt20 <- sum(sz_clean > 20, na.rm = TRUE)
n_sz_miss <- sum(is.na(sz) | sz >= 900)
erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "\u226410 mm", Value = fmt_n(n_le10, N_enc_erp)),
  data.frame(Characteristics = "11-20 mm", Value = fmt_n(n_11_20, N_enc_erp)),
  data.frame(Characteristics = ">20 mm", Value = fmt_n(n_gt20, N_enc_erp))
)
if (n_sz_miss > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Size missing/multifocal", Value = fmt_n(n_sz_miss, N_enc_erp)))

## Surgery (TREAT: 1=mastectomy, 2=sector, 4=no surgery; anything else -> Unknown/other)
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Surgery", Value = ""))
tr <- enc_erp$TREAT
n_mast   <- sum(tr == 1, na.rm = TRUE)
n_sect   <- sum(tr == 2, na.rm = TRUE)
n_nosurg <- sum(tr == 4, na.rm = TRUE)
n_other  <- sum(is.na(tr) | !(tr %in% c(1, 2, 4)))
if (n_mast > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Mastectomy", Value = fmt_n(n_mast, N_enc_erp)))
if (n_sect > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Breast-conserving", Value = fmt_n(n_sect, N_enc_erp)))
if (n_nosurg > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "No surgery", Value = fmt_n(n_nosurg, N_enc_erp)))
if (n_other > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Unknown/other", Value = fmt_n(n_other, N_enc_erp)))

## Radiotherapy
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Radiotherapy", Value = ""))
rt <- enc_erp$RT
erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(sum(rt == 1, na.rm = TRUE), N_enc_erp)),
  data.frame(Characteristics = "No", Value = fmt_n(sum(rt == 2, na.rm = TRUE), N_enc_erp))
)

## Chemotherapy
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Chemotherapy", Value = ""))
chemo <- enc_erp$CHEMOTHERAPY
erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(sum(chemo == 1, na.rm = TRUE), N_enc_erp)),
  data.frame(Characteristics = "No", Value = fmt_n(sum(chemo == 2, na.rm = TRUE), N_enc_erp))
)

## Endocrine therapy
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Endocrine therapy", Value = ""))
endo <- enc_erp$ENDOCRINETHERAPY
n_en_yes <- sum(endo == 1, na.rm = TRUE)
n_en_no  <- sum(endo == 2, na.rm = TRUE)
n_en_oth <- sum(!is.na(endo) & !(endo %in% c(1, 2)))
erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_en_yes, N_enc_erp)),
  data.frame(Characteristics = "No", Value = fmt_n(n_en_no, N_enc_erp))
)
if (n_en_oth > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "Unknown", Value = fmt_n(n_en_oth, N_enc_erp)))

## Diagnosis period
erp_rows <- rbind(erp_rows, data.frame(Characteristics = "Diagnosis period", Value = ""))
dy <- enc_erp$DIAG_YEAR
dy_full <- ifelse(!is.na(dy) & dy < 100, dy + 1900, dy)
n_85_89 <- sum(dy_full >= 1985 & dy_full <= 1989, na.rm = TRUE)
n_90_94 <- sum(dy_full >= 1990 & dy_full <= 1994, na.rm = TRUE)
n_95_99 <- sum(dy_full >= 1995 & dy_full <= 1999, na.rm = TRUE)
n_00_04 <- sum(dy_full >= 2000 & dy_full <= 2004, na.rm = TRUE)
if (n_85_89 > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "1985-1989", Value = fmt_n(n_85_89, N_enc_erp)))
if (n_90_94 > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "1990-1994", Value = fmt_n(n_90_94, N_enc_erp)))
if (n_95_99 > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "1995-1999", Value = fmt_n(n_95_99, N_enc_erp)))
if (n_00_04 > 0) erp_rows <- rbind(erp_rows,
  data.frame(Characteristics = "2000-2004", Value = fmt_n(n_00_04, N_enc_erp)))

## Set column name
colnames(erp_rows) <- c("Characteristics", sprintf("Cases (percentage)\nN = %d", N_enc_erp))

## Save
erp_table_path <- file.path(output_dir, "Table_validation_ERpos_clinicopathological.docx")
save_clinpath_docx(erp_rows, erp_table_path,
  sprintf("Table S4. Distribution of clinicopathological parameters in ER-positive tumors within the validation cohort (MSPI-available subset, N = %d).", N_enc_erp))


cat("\n═══════════════════════════════════════════════\n")
cat("PART 2: Complete — all clinicopathological tables saved\n")
cat("═══════════════════════════════════════════════\n")
