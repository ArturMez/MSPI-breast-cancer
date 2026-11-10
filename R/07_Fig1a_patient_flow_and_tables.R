########################
## Figure 1a (patient flow) and Supplementary Tables S2-S3 (clinicopathological characteristics)
##
## Run from the repository root. Inputs (original names; the repository names in brackets also work):
##   ./data/discovery/   BCY_workingr_832_updated_211024.sav [cohort_clinical.sav, clinical_data.sav]
##                       BCY_distances_mean_M1_to_CD3.txt, BCY_distances_mean_M2_to_CD3.txt
##                       [M1_to_CD3_distances.txt, M2_to_CD3_distances.txt]
##   ./data/validation/  NK_BRCA_for_affinity.txt [clinical.txt]
##                       BREAST_ENCIMA_affinity_CD3e_M1M2_myeloidpool_50um.txt [affinity_scores.txt]
##                       Huvuddatabas.xlsx [clinical_extra.xlsx]
## Outputs: ./output/fig1a_and_tables/
## In the code, bcy_* objects hold the discovery cohort and enc_* objects the validation cohort.
########################

library(grid)
library(gridExtra)

## ── Paths ──────────────────────────────────────────────────────────────────
proj_root  <- "./data"
output_dir <- "./output/fig1a_and_tables"
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

## ── Quartz/cairo fallback ──────────────────────────────────────────────────
FONT_FAM <- "Arial"
USE_QUARTZ <- FALSE
if (.Platform$OS.type != "windows") {
  USE_QUARTZ <- tryCatch({
    tf <- tempfile(fileext = ".pdf")
    cairo_pdf(tf, width = 1, height = 1)
    plot.new()
    dev.off()
    unlink(tf)
    FALSE  # cairo works
  }, warning = function(w) TRUE, error = function(e) TRUE)
}

save_pdf <- function(fpath, w, h) {
  ## Ensure parent directory exists
  dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
  if (USE_QUARTZ) {
    quartz(type = "pdf", file = fpath, width = w, height = h,
           family = FONT_FAM, bg = "white")
  } else {
    tryCatch(
      cairo_pdf(fpath, width = w, height = h, bg = "white"),
      error = function(e) {
        message("cairo_pdf failed, falling back to quartz: ", e$message)
        quartz(type = "pdf", file = fpath, width = w, height = h,
               family = FONT_FAM, bg = "white")
      }
    )
  }
}

## ── Colour palette ─────────────────────────────────────────────────────────
COL_DISC   <- "#7ba3c4"   # Discovery cohort (blue)
COL_VAL    <- "#c47a7a"   # Validation cohort (red-ish)
COL_EXPL   <- "#a0a0a0"   # Exploratory cohort (grey)
COL_POOL   <- "#d4b86a"   # Pooled IPD (amber)
COL_BOX    <- "#f5f5f0"   # Light fill for boxes
COL_DISC_L <- "#e8f0f7"   # Light fill for discovery boxes
COL_VAL_L  <- "#fae8e8"   # Light fill for validation boxes
COL_POOL_L <- "#faf3e0"   # Light fill for pooled box

## ── Helper: draw a rounded box with text ───────────────────────────────────
draw_box <- function(x, y, w, h, labels, fill = COL_BOX, border = "grey40",
                     cex = 0.65, font = 1, line_spacing = 1.15) {
  grid.roundrect(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(w, "npc"), height = unit(h, "npc"),
    r = unit(0.008, "npc"),
    gp = gpar(fill = fill, col = border, lwd = 1.2)
  )
  n <- length(labels)
  for (i in seq_along(labels)) {
    lbl <- labels[[i]]
    f <- if (is.list(lbl)) lbl$font else font
    cx <- if (is.list(lbl)) lbl$cex else cex
    txt <- if (is.list(lbl)) lbl$text else lbl
    yoff <- y + h * 0.5 * (1 - (2 * i - 1) / n) * 0.85
    grid.text(txt, x = unit(x, "npc"), y = unit(yoff, "npc"),
              gp = gpar(fontsize = cx * 12, fontfamily = FONT_FAM, fontface = f))
  }
}

## ── Helper: draw an arrow ──────────────────────────────────────────────────
draw_arrow <- function(x0, y0, x1, y1, col = "grey40") {
  grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
             arrow = arrow(length = unit(0.012, "npc"), type = "closed"),
             gp = gpar(col = col, fill = col, lwd = 1.2))
}

draw_line <- function(x0, y0, x1, y1, col = "grey40", lty = 1) {
  grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
             gp = gpar(col = col, lwd = 1.2, lty = lty))
}


## ════════════════════════════════════════════════════════════════════════════
## DRAW THE FIGURE
## ════════════════════════════════════════════════════════════════════════════

fpath <- file.path(output_dir, "Fig1a_patient_flow.pdf")
save_pdf(fpath, w = 12, h = 8)
grid.newpage()

## Background
grid.rect(gp = gpar(fill = "white", col = NA))

## ── Column headers ─────────────────────────────────────────────────────────
grid.text("DISCOVERY COHORT", x = 0.25, y = 0.96,
          gp = gpar(fontsize = 14, fontface = "bold", fontfamily = FONT_FAM, col = "#4a6e91"))
grid.text("VALIDATION COHORT", x = 0.75, y = 0.96,
          gp = gpar(fontsize = 14, fontface = "bold", fontfamily = FONT_FAM, col = "#91504a"))

## ── Exploratory cohort (top bar) ───────────────────────────────────────────
draw_box(0.50, 0.935, 0.55, 0.025,
         list(list(text = "EXPLORATORY COHORT  n = 10 patients (17 whole-tissue sections)", font = "bold", cex = 0.6)),
         fill = "#e8e8e8", border = COL_EXPL)

## ══════════════════════════════════════════════════════════════════════════
## DISCOVERY SIDE (left)
## ══════════════════════════════════════════════════════════════════════════

## Registry box
draw_box(0.25, 0.87, 0.38, 0.055,
         list(list(text = "Swedish population-based registry", font = "bold", cex = 0.65),
              list(text = "Women diagnosed 1992-2005", font = 1, cex = 0.55),
              list(text = "Age-stratified sampling (enriched <40 y)", font = 1, cex = 0.55)),
         fill = COL_DISC_L, border = COL_DISC)

draw_arrow(0.25, 0.843, 0.25, 0.82)

## MSPI available (both M1 and M2 affinity scores + valid DDFS)
draw_box(0.25, 0.795, 0.38, 0.04,
         list(list(text = "MSPI available", font = "bold", cex = 0.6),
              list(text = "n = 687, 189 DDFS events", font = 1, cex = 0.55)),
         fill = COL_DISC_L, border = COL_DISC)

## Arrows to molecular subtypes
draw_line(0.25, 0.775, 0.25, 0.685)  # vertical
draw_line(0.10, 0.685, 0.40, 0.685)  # horizontal
draw_arrow(0.10, 0.685, 0.10, 0.665)
draw_arrow(0.25, 0.685, 0.25, 0.665)
draw_arrow(0.40, 0.685, 0.40, 0.665)

## Molecular subtypes
draw_box(0.10, 0.64, 0.15, 0.04,
         list(list(text = "Luminal A+B", font = "bold", cex = 0.55),
              list(text = "n = 415, 103 events", font = 1, cex = 0.50)),
         fill = COL_DISC_L, border = COL_DISC)

draw_box(0.25, 0.64, 0.12, 0.04,
         list(list(text = "HER2+", font = "bold", cex = 0.55),
              list(text = "n = 114, 36 events", font = 1, cex = 0.50)),
         fill = COL_DISC_L, border = COL_DISC)

draw_box(0.40, 0.64, 0.14, 0.04,
         list(list(text = "Triple-negative", font = "bold", cex = 0.55),
              list(text = "n = 158, 50 events", font = 1, cex = 0.50)),
         fill = COL_DISC_L, border = COL_DISC)

## Arrows from Luminal to subgroups
draw_line(0.10, 0.62, 0.10, 0.60)
draw_line(0.10, 0.60, 0.20, 0.60)
draw_arrow(0.10, 0.60, 0.10, 0.58)
draw_arrow(0.20, 0.60, 0.20, 0.58)

## Early-stage luminal
draw_box(0.10, 0.55, 0.16, 0.045,
         list(list(text = "Early-stage luminal", font = "bold", cex = 0.55),
              list(text = "Stage I-IIb", font = "italic", cex = 0.50),
              list(text = "n = 360, 77 events", font = 1, cex = 0.50)),
         fill = COL_DISC_L, border = COL_DISC)

## Young luminal
draw_box(0.235, 0.55, 0.14, 0.045,
         list(list(text = "Young luminal", font = "bold", cex = 0.55),
              list(text = "Age <40", font = "italic", cex = 0.50),
              list(text = "n = 152, 51 events", font = 1, cex = 0.50)),
         fill = COL_DISC_L, border = COL_DISC)

## Endpoint box
draw_box(0.17, 0.48, 0.28, 0.03,
         list(list(text = "Endpoint: DDFS (distant disease-free survival)", font = "bold.italic", cex = 0.55)),
         fill = "#dbe8f4", border = COL_DISC)


## ══════════════════════════════════════════════════════════════════════════
## VALIDATION SIDE (right)
## ══════════════════════════════════════════════════════════════════════════

## Registry box
draw_box(0.75, 0.87, 0.38, 0.055,
         list(list(text = "Uppsala / Vasteras population-based", font = "bold", cex = 0.65),
              list(text = "Women diagnosed 1985-2004", font = 1, cex = 0.55),
              list(text = "TMA collection, Opal multiplex IF", font = 1, cex = 0.55)),
         fill = COL_VAL_L, border = COL_VAL)

draw_arrow(0.75, 0.843, 0.75, 0.82)

## Full cohort MSPI available
draw_box(0.75, 0.795, 0.38, 0.04,
         list(list(text = "Full cohort, MSPI available", font = "bold", cex = 0.6),
              list(text = "n = 305, 132 OS events", font = 1, cex = 0.55)),
         fill = COL_VAL_L, border = COL_VAL)

## Arrows to ER subgroups
draw_line(0.75, 0.775, 0.75, 0.755)
draw_line(0.65, 0.755, 0.85, 0.755)
draw_arrow(0.65, 0.755, 0.65, 0.735)
draw_arrow(0.85, 0.755, 0.85, 0.735)

## ER-positive
draw_box(0.65, 0.71, 0.18, 0.04,
         list(list(text = "ER-positive", font = "bold", cex = 0.55),
              list(text = "n = 247, 105 OS events", font = 1, cex = 0.50)),
         fill = COL_VAL_L, border = COL_VAL)

## ER-negative
draw_box(0.85, 0.71, 0.18, 0.04,
         list(list(text = "ER-negative", font = "bold", cex = 0.55),
              list(text = "n = 58, 27 OS events", font = 1, cex = 0.50)),
         fill = COL_VAL_L, border = COL_VAL)

## Arrow to ER+ N0
draw_arrow(0.65, 0.69, 0.65, 0.665)

## ER+ N0
draw_box(0.65, 0.64, 0.21, 0.04,
         list(list(text = "ER+ node-negative (N0)", font = "bold", cex = 0.55),
              list(text = "n = 206, 89 OS events", font = 1, cex = 0.50)),
         fill = COL_VAL_L, border = COL_VAL)

## Arrows to age subgroups
draw_line(0.65, 0.62, 0.65, 0.605)
draw_line(0.555, 0.605, 0.745, 0.605)
draw_arrow(0.555, 0.605, 0.555, 0.585)
draw_arrow(0.65, 0.605, 0.65, 0.585)
draw_arrow(0.745, 0.605, 0.745, 0.585)

## Age <=50
draw_box(0.555, 0.555, 0.13, 0.045,
         list(list(text = "Age <=50 y", font = "bold", cex = 0.50),
              list(text = "n = 43", font = 1, cex = 0.48),
              list(text = "11 OS events", font = 1, cex = 0.48)),
         fill = COL_VAL_L, border = COL_VAL)

## Age 51-64
draw_box(0.65, 0.555, 0.13, 0.045,
         list(list(text = "Age 51-64 y", font = "bold", cex = 0.50),
              list(text = "n = 76", font = 1, cex = 0.48),
              list(text = "19 OS events", font = 1, cex = 0.48)),
         fill = COL_VAL_L, border = COL_VAL)

## Age >=65
draw_box(0.745, 0.555, 0.13, 0.045,
         list(list(text = "Age >=65 y", font = "bold", cex = 0.50),
              list(text = "n = 87", font = 1, cex = 0.48),
              list(text = "59 OS events", font = 1, cex = 0.48)),
         fill = COL_VAL_L, border = COL_VAL)

## Endpoint box
draw_box(0.75, 0.48, 0.28, 0.03,
         list(list(text = "Endpoint: OS (overall survival)", font = "bold.italic", cex = 0.55)),
         fill = "#fae0e0", border = COL_VAL)


## ══════════════════════════════════════════════════════════════════════════
## POOLED IPD (bottom)
## ══════════════════════════════════════════════════════════════════════════

## Connecting arrows from Discovery and Validation to pooled box
draw_line(0.17, 0.465, 0.17, 0.42, lty = 2)
draw_line(0.75, 0.465, 0.75, 0.42, lty = 2)
draw_line(0.17, 0.42, 0.75, 0.42, lty = 2)
draw_arrow(0.46, 0.42, 0.46, 0.40)

## Pooled IPD box
draw_box(0.50, 0.355, 0.70, 0.07,
         list(list(text = "Pooled Individual Patient Data Meta-analysis", font = "bold", cex = 0.7),
              list(text = "Stratified Cox regression, within-cohort standardization of MSPI", font = "italic", cex = 0.55),
              list(text = "Discovery young luminal + Validation ER+ young N0 (n = 195, 62 events)", font = 1, cex = 0.55)),
         fill = COL_POOL_L, border = COL_POOL)

dev.off()
cat("Saved:", fpath, "\n")


## ════════════════════════════════════════════════════════════════════════════
## PART 2 — CLINICOPATHOLOGICAL TABLES FOR MSPI-AVAILABLE SUBSETS
## Produces .docx tables (styled via flextable + officer)
## for BCY (discovery) and ENCIMA (validation)
## ════════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════\n")
cat("PART 2: Clinicopathological Tables\n")
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
m1_raw <- data.table::fread(find_input(bcy_base, "BCY_distances_mean_M1_to_CD3.txt", "M1_to_CD3_distances.txt"),
                            data.table = FALSE)
m2_raw <- data.table::fread(find_input(bcy_base, "BCY_distances_mean_M2_to_CD3.txt", "M2_to_CD3_distances.txt"),
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
