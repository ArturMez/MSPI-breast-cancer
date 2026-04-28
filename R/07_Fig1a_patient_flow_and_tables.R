########################
## Figure 1a - Patient Flow and Cohort Structure
## Corrected version with verified numbers from analysis outputs
########################

library(grid)
library(gridExtra)

## ── Output ─────────────────────────────────────────────────────────────────
## Set project root and output directories
## Modify these paths as needed for your local environment
proj_root <- "./data"  # Set this to your project data directory
output_dir <- file.path(proj_root, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Project root:", proj_root, "\n")
cat("Output dir:", output_dir, "\n")

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
COL_EXPL   <- "#a0a0a0"   # Exploratory cohort (grey) — for reference only
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

fpath <- file.path(output_dir, "Fig_patient_flow.pdf")
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
## for the Discovery and Validation cohorts
## ════════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════\n")
cat("PART 2: Clinicopathological Tables (Discovery & Validation Cohorts)\n")
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
## Discovery Cohort
## ══════════════════════════════════════════════════════════════════════════

cat("Loading discovery cohort data...\n")

disc_base <- file.path(proj_root, "discovery")

disc_spss_candidates <- c(
  file.path(disc_base, "cohort_clinical.sav"),
  file.path(disc_base, "clinical_data.sav")
)
disc_spss <- NA
for (p in disc_spss_candidates) {
  if (file.exists(p)) { disc_spss <- p; break }
}
if (is.na(disc_spss)) stop("Cannot find discovery cohort SPSS file")

disc_raw <- foreign::read.spss(
  disc_spss, use.value.labels = TRUE, to.data.frame = TRUE,
  max.value.labels = Inf, trim.factor.names = TRUE,
  trim_values = TRUE, reencode = NA
)

## Load affinity distances
m1_raw <- data.table::fread(file.path(disc_base, "M1_to_CD3_distances.txt"),
                            data.table = FALSE)
m2_raw <- data.table::fread(file.path(disc_base, "M2_to_CD3_distances.txt"),
                            data.table = FALSE)
colnames(m1_raw)[colnames(m1_raw) == "delta_CD3"] <- "M1_to_CD3_distance"
colnames(m2_raw)[colnames(m2_raw) == "delta_CD3"] <- "M2_to_CD3_distance"

## Merge
disc <- disc_raw %>%
  dplyr::left_join(m1_raw[, c("patient_id", "M1_to_CD3_distance")], by = "patient_id") %>%
  dplyr::left_join(m2_raw[, c("patient_id", "M2_to_CD3_distance")], by = "patient_id")

## Remove unclassified
disc <- disc[as.character(disc$Subtype) != "Unclassified", ]
disc$Subtype <- droplevels(disc$Subtype)

## Define subsets
disc$has_M1 <- !is.na(disc$M1_to_CD3_distance)
disc$has_M2 <- !is.na(disc$M2_to_CD3_distance)
disc$has_MSPI <- disc$has_M1 & disc$has_M2

## Survival filter (valid event time)
disc$status_num <- ifelse(
  as.character(disc$Event) == "Distant recurrence" |
    as.character(disc$Event) == "Death",
  1L, 0L
)
disc$time_to_event <- as.numeric(disc$time_to_event)
disc$valid_surv <- !is.na(disc$time_to_event) & disc$time_to_event > 0

cat(sprintf("  Discovery cohort total (excl. unclassified): n=%d\n", nrow(disc)))
cat(sprintf("  M1 available: n=%d\n", sum(disc$has_M1)))
cat(sprintf("  M2 available: n=%d\n", sum(disc$has_M2)))
cat(sprintf("  MSPI available: n=%d\n", sum(disc$has_MSPI)))
cat(sprintf("  MSPI + valid survival: n=%d\n", sum(disc$has_MSPI & disc$valid_surv)))

## MSPI analysis subset
disc_mspi <- disc %>% dplyr::filter(has_MSPI, valid_surv)
N_disc <- nrow(disc_mspi)

cat(sprintf("  Discovery cohort MSPI analysis set: n=%d, events=%d\n", N_disc, sum(disc_mspi$status_num)))

## Build table rows
disc_age_num <- suppressWarnings(as.numeric(as.character(disc_mspi$age)))
if (all(is.na(disc_age_num))) {
  ## age might be stored as a labelled variable - try age_cat
  disc_age_num <- NULL
}

disc_rows <- data.frame(
  Characteristics = character(),
  Value = character(),
  stringsAsFactors = FALSE
)

## Age
disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Age", Value = ""))
if (!is.null(disc_age_num) && !all(is.na(disc_age_num))) {
  disc_rows <- rbind(disc_rows,
    data.frame(Characteristics = "Mean age, years",
               Value = sprintf("%.1f", mean(disc_age_num, na.rm = TRUE))))
}
age_cat <- as.character(disc_mspi$age_cat)
n_young <- sum(age_cat == "<40", na.rm = TRUE)
n_old   <- sum(age_cat == ">=40", na.rm = TRUE)
disc_rows <- rbind(disc_rows,
  data.frame(Characteristics = "<40", Value = fmt_n(n_young, N_disc)),
  data.frame(Characteristics = "\u226540", Value = fmt_n(n_old, N_disc))
)

## Stage
disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Stage", Value = ""))
stage <- as.character(disc_mspi$stage)
for (s in c("Stage I", "Stage IIa", "Stage IIb")) {
  n_s <- sum(stage == s, na.rm = TRUE)
  disc_rows <- rbind(disc_rows, data.frame(Characteristics = s, Value = fmt_n(n_s, N_disc)))
}
n_iii_iv <- sum(stage %in% c("Stage IIIa", "Stage IIIb", "Stage IIIc", "Stage IV",
                              "Stage III", "Stage III-IV"), na.rm = TRUE)
disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Stage III-IV", Value = fmt_n(n_iii_iv, N_disc)))
n_stage_miss <- sum(is.na(stage) | stage == "" | stage == "missing")
if (n_stage_miss > 0) {
  disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Stage missing", Value = fmt_n(n_stage_miss, N_disc)))
}

## Grade
disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Grade", Value = ""))
grade <- as.character(disc_mspi$grade)
for (g in c("Grade I", "Grade II", "Grade III")) {
  n_g <- sum(grade == g, na.rm = TRUE)
  if (n_g > 0) disc_rows <- rbind(disc_rows, data.frame(Characteristics = g, Value = fmt_n(n_g, N_disc)))
}
n_grade_miss <- sum(is.na(grade) | grade == "" | grepl("[Mm]issing", grade))
if (n_grade_miss > 0) {
  disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Grade missing", Value = fmt_n(n_grade_miss, N_disc)))
}

## Molecular Subtype
disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Molecular subtype", Value = ""))
subtype <- as.character(disc_mspi$Subtype)
for (st in c("Lum A", "Lum B", "Lum Her2", "Her2 exp", "Triple neg")) {
  n_st <- sum(subtype == st, na.rm = TRUE)
  disc_rows <- rbind(disc_rows, data.frame(Characteristics = st, Value = fmt_n(n_st, N_disc)))
}

## Systemic treatment
disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Systemic treatment", Value = ""))
treat <- as.character(disc_mspi$systemic_treatment)
treat_levels <- sort(unique(treat[!is.na(treat) & treat != ""]))
for (tl in treat_levels) {
  n_tl <- sum(treat == tl, na.rm = TRUE)
  disc_rows <- rbind(disc_rows, data.frame(Characteristics = tl, Value = fmt_n(n_tl, N_disc)))
}
n_treat_miss <- sum(is.na(treat) | treat == "")
if (n_treat_miss > 0) {
  disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Missing", Value = fmt_n(n_treat_miss, N_disc)))
}

## Diagnosis period
disc_rows <- rbind(disc_rows, data.frame(Characteristics = "Diagnosis period", Value = ""))
diag_p <- as.character(disc_mspi$diagnosis_year)
diag_levels <- sort(unique(diag_p[!is.na(diag_p) & diag_p != ""]))
for (dl in diag_levels) {
  n_dl <- sum(diag_p == dl, na.rm = TRUE)
  disc_rows <- rbind(disc_rows, data.frame(Characteristics = dl, Value = fmt_n(n_dl, N_disc)))
}

## Set column name
colnames(disc_rows) <- c("Characteristics", sprintf("Cases (percentage)\nN = %d", N_disc))

## Save
disc_table_path <- file.path(output_dir, "Table_discovery_clinicopathological.docx")
save_clinpath_docx(disc_rows, disc_table_path,
  sprintf("Table. Distribution of clinicopathological parameters in discovery cohort (MSPI-available subset, N = %d).", N_disc))


## ══════════════════════════════════════════════════════════════════════════
## Validation Cohort
## ══════════════════════════════════════════════════════════════════════════

cat("\nLoading validation cohort data...\n")

val_base <- file.path(proj_root, "validation")

val_clin <- data.table::fread(file.path(val_base, "clinical.txt"),
                              data.table = FALSE)
val_aff <- data.table::fread(
  file.path(val_base, "affinity_scores.txt"),
  data.table = FALSE
)
val <- dplyr::left_join(val_clin, val_aff, by = "case_id")

## Load additional clinical variables
val_extra <- openxlsx::read.xlsx(file.path(val_base, "clinical_extra.xlsx"))
## Rename columns to avoid issues with special characters
patient_col <- grep("^patient_id", names(val_extra), value = TRUE)[1]
if (!is.na(patient_col) && patient_col != "patient_id") {
  names(val_extra)[names(val_extra) == patient_col] <- "patient_id"
  cat("  Clinical data patient column found as:", patient_col, "-> renamed to patient_id\n")
}
## Similarly fix diagnosis year column
diag_col <- grep("^diag", names(val_extra), value = TRUE)[1]
if (!is.na(diag_col) && diag_col != "diagnosis_year") {
  names(val_extra)[names(val_extra) == diag_col] <- "diagnosis_year"
  cat("  Clinical data diagnosis year column found as:", diag_col, "-> renamed to diagnosis_year\n")
}

## Extract patient identifier for merge
val$patient_id <- as.numeric(sub(".*_(\\d+)_.*", "\\1", val$case_id))
val <- dplyr::left_join(val, val_extra, by = "patient_id")
cat(sprintf("  Clinical data matched: %d / %d cases\n",
            sum(!is.na(val$age)), nrow(val)))

## Define MSPI availability
val$has_M1 <- !is.na(val$M1_to_CD3)
val$has_M2 <- !is.na(val$M2_to_CD3)
val$has_MSPI <- val$has_M1 & val$has_M2

## Survival
val$os_days <- suppressWarnings(as.numeric(val$time_to_followup_os))
val$status_os <- ifelse(val$vital_status == "Dead", 1L, 0L)
val$valid_surv <- !is.na(val$os_days) & val$os_days > 0

cat(sprintf("  Validation cohort total: n=%d\n", nrow(val)))
cat(sprintf("  M1 available: n=%d\n", sum(val$has_M1)))
cat(sprintf("  M2 available: n=%d\n", sum(val$has_M2)))
cat(sprintf("  MSPI available: n=%d\n", sum(val$has_MSPI)))

## ----- All invasive MSPI subset -----
val_inv <- val %>%
  dplyr::filter(
    tumor_type %in% c("invasive_ER_pos", "invasive_ER_neg", "invasive"),
    has_MSPI, valid_surv
  )
N_val_inv <- nrow(val_inv)
cat(sprintf("  Validation cohort all invasive MSPI + valid OS: n=%d, events=%d\n",
            N_val_inv, sum(val_inv$status_os)))

## Build table for ALL invasive validation
val_rows <- data.frame(Characteristics = character(), Value = character(),
                       stringsAsFactors = FALSE)

## Age
val_rows <- rbind(val_rows, data.frame(Characteristics = "Age", Value = ""))
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Mean age, years",
             Value = sprintf("%.1f", mean(val_inv$age, na.rm = TRUE))))
age_cat <- as.character(val_inv$age_cat)
for (ac in c("<=50", ">50-<65", ">=65")) {
  n_ac <- sum(age_cat == ac, na.rm = TRUE)
  val_rows <- rbind(val_rows, data.frame(Characteristics = ac, Value = fmt_n(n_ac, N_val_inv)))
}

## ER status
val_rows <- rbind(val_rows, data.frame(Characteristics = "ER status", Value = ""))
n_erp <- sum(val_inv$tumor_type == "invasive_ER_pos", na.rm = TRUE)
n_ern <- sum(val_inv$tumor_type %in% c("invasive_ER_neg", "BRCA"), na.rm = TRUE)
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "ER-positive", Value = fmt_n(n_erp, N_val_inv)),
  data.frame(Characteristics = "ER-negative", Value = fmt_n(n_ern, N_val_inv))
)

## Grade (from grade in NK_BRCA, more complete than GRADE in Huvuddatabas)
val_rows <- rbind(val_rows, data.frame(Characteristics = "Grade", Value = ""))
dg <- val_inv$grade
for (g in c("1", "2", "3")) {
  n_g <- sum(dg == g, na.rm = TRUE)
  if (n_g > 0) val_rows <- rbind(val_rows,
    data.frame(Characteristics = paste0("Grade ", g), Value = fmt_n(n_g, N_val_inv)))
}
n_g_miss <- sum(dg == "missing" | is.na(dg), na.rm = TRUE)
if (n_g_miss > 0) {
  val_rows <- rbind(val_rows,
    data.frame(Characteristics = "Grade missing", Value = fmt_n(n_g_miss, N_val_inv)))
}

## Nodal status (from NK_BRCA node_status)
val_rows <- rbind(val_rows, data.frame(Characteristics = "Nodal status", Value = ""))
pn <- val_inv$node_status
n_n0 <- sum(pn == "N0", na.rm = TRUE)
n_n1 <- sum(pn == "N+", na.rm = TRUE)
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Node-negative (N0)", Value = fmt_n(n_n0, N_val_inv)),
  data.frame(Characteristics = "Node-positive", Value = fmt_n(n_n1, N_val_inv))
)

## Tumor size (from Huvuddatabas SIZE, mm)
val_rows <- rbind(val_rows, data.frame(Characteristics = "Tumor size", Value = ""))
sz <- val_inv$tumor_size_mm
sz_clean <- sz[!is.na(sz) & sz < 900]  # exclude 999 (multifocal) and 88 (missing)
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Mean size, mm",
             Value = sprintf("%.1f", mean(sz_clean, na.rm = TRUE))))
n_le10 <- sum(sz_clean <= 10, na.rm = TRUE)
n_11_20 <- sum(sz_clean > 10 & sz_clean <= 20, na.rm = TRUE)
n_gt20 <- sum(sz_clean > 20, na.rm = TRUE)
n_sz_miss <- sum(is.na(sz) | sz >= 900)  # 999 = multifocal, 88 = missing
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "\u226410 mm", Value = fmt_n(n_le10, N_val_inv)),
  data.frame(Characteristics = "11-20 mm", Value = fmt_n(n_11_20, N_val_inv)),
  data.frame(Characteristics = ">20 mm", Value = fmt_n(n_gt20, N_val_inv))
)
if (n_sz_miss > 0) {
  val_rows <- rbind(val_rows,
    data.frame(Characteristics = "Size missing/multifocal", Value = fmt_n(n_sz_miss, N_val_inv)))
}

## Surgery type (from Huvuddatabas TREAT: 1=mastectomy, 2=sector, 4=no surgery)
val_rows <- rbind(val_rows, data.frame(Characteristics = "Surgery", Value = ""))
tr <- val_inv$surgery_type
n_mast <- sum(tr == "mastectomy", na.rm = TRUE)
n_sect <- sum(tr == "breast_conserving", na.rm = TRUE)
n_nosurg <- sum(tr == "none", na.rm = TRUE)
if (n_mast > 0) val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Mastectomy", Value = fmt_n(n_mast, N_val_inv)))
if (n_sect > 0) val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Breast-conserving", Value = fmt_n(n_sect, N_val_inv)))
if (n_nosurg > 0) val_rows <- rbind(val_rows,
  data.frame(Characteristics = "No surgery", Value = fmt_n(n_nosurg, N_val_inv)))

## Radiotherapy (from Huvuddatabas RT: 1=yes, 2=no)
val_rows <- rbind(val_rows, data.frame(Characteristics = "Radiotherapy", Value = ""))
rt <- val_inv$radiotherapy
n_rt_yes <- sum(rt == "yes", na.rm = TRUE)
n_rt_no  <- sum(rt == "no", na.rm = TRUE)
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_rt_yes, N_val_inv)),
  data.frame(Characteristics = "No", Value = fmt_n(n_rt_no, N_val_inv))
)

## Chemotherapy (from Huvuddatabas CHEMOTHERAPY: 1=yes, 2=no)
val_rows <- rbind(val_rows, data.frame(Characteristics = "Chemotherapy", Value = ""))
chemo <- val_inv$chemotherapy
n_ch_yes <- sum(chemo == "yes", na.rm = TRUE)
n_ch_no  <- sum(chemo == "no", na.rm = TRUE)
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_ch_yes, N_val_inv)),
  data.frame(Characteristics = "No", Value = fmt_n(n_ch_no, N_val_inv))
)

## Endocrine therapy (from Huvuddatabas ENDOCRINETHERAPY: 1=yes, 2=no, 4=?)
val_rows <- rbind(val_rows, data.frame(Characteristics = "Endocrine therapy", Value = ""))
endo <- val_inv$endocrine_therapy
n_en_yes <- sum(endo == "yes", na.rm = TRUE)
n_en_no  <- sum(endo == "no", na.rm = TRUE)
n_en_oth <- sum(!is.na(endo) & !(endo %in% c("yes", "no")))
val_rows <- rbind(val_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_en_yes, N_val_inv)),
  data.frame(Characteristics = "No", Value = fmt_n(n_en_no, N_val_inv))
)
if (n_en_oth > 0) {
  val_rows <- rbind(val_rows,
    data.frame(Characteristics = "Unknown", Value = fmt_n(n_en_oth, N_val_inv)))
}

## Diagnosis period
val_rows <- rbind(val_rows, data.frame(Characteristics = "Diagnosis period", Value = ""))
dy <- val_inv$diagnosis_year
## Convert 2-digit years: 85-99 -> 1985-1999, 2002-2004 stay as-is
dy_full <- ifelse(!is.na(dy) & dy < 100, dy + 1900, dy)
n_85_89 <- sum(dy_full >= 1985 & dy_full <= 1989, na.rm = TRUE)
n_90_94 <- sum(dy_full >= 1990 & dy_full <= 1994, na.rm = TRUE)
n_95_99 <- sum(dy_full >= 1995 & dy_full <= 1999, na.rm = TRUE)
n_00_04 <- sum(dy_full >= 2000 & dy_full <= 2004, na.rm = TRUE)
if (n_85_89 > 0) val_rows <- rbind(val_rows,
  data.frame(Characteristics = "1985-1989", Value = fmt_n(n_85_89, N_val_inv)))
if (n_90_94 > 0) val_rows <- rbind(val_rows,
  data.frame(Characteristics = "1990-1994", Value = fmt_n(n_90_94, N_val_inv)))
if (n_95_99 > 0) val_rows <- rbind(val_rows,
  data.frame(Characteristics = "1995-1999", Value = fmt_n(n_95_99, N_val_inv)))
if (n_00_04 > 0) val_rows <- rbind(val_rows,
  data.frame(Characteristics = "2000-2004", Value = fmt_n(n_00_04, N_val_inv)))

## Set column name
colnames(val_rows) <- c("Characteristics", sprintf("Cases (percentage)\nN = %d", N_val_inv))

## Save
val_table_path <- file.path(output_dir, "Table_validation_clinicopathological.docx")
save_clinpath_docx(val_rows, val_table_path,
  sprintf("Table S3. Distribution of clinicopathological parameters in tumors within the validation cohort (MSPI-available subset, all invasive, N = %d).", N_val_inv))


## ----- ER-positive MSPI subset (primary analysis population) -----
val_erp <- val %>%
  dplyr::filter(
    tumor_type == "invasive_ER_pos",
    has_MSPI, valid_surv
  )
N_val_erp <- nrow(val_erp)
cat(sprintf("  Validation cohort ER+ MSPI + valid OS: n=%d, events=%d\n",
            N_val_erp, sum(val_erp$status_os)))

## Build ER+ table (similar structure, minus ER status row)
val_erp_rows <- data.frame(Characteristics = character(), Value = character(),
                       stringsAsFactors = FALSE)

## Age
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Age", Value = ""))
val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Mean age, years",
             Value = sprintf("%.1f", mean(val_erp$AGEYEARS, na.rm = TRUE))))
age_cat <- as.character(val_erp$Age)
for (ac in c("<=50", ">50-<65", ">=65")) {
  n_ac <- sum(age_cat == ac, na.rm = TRUE)
  val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = ac, Value = fmt_n(n_ac, N_val_erp)))
}

## Grade
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Grade", Value = ""))
dg <- val_erp$grade
for (g in c("1", "2", "3")) {
  n_g <- sum(dg == g, na.rm = TRUE)
  if (n_g > 0) val_erp_rows <- rbind(val_erp_rows,
    data.frame(Characteristics = paste0("Grade ", g), Value = fmt_n(n_g, N_val_erp)))
}
n_g_miss <- sum(dg == "missing" | is.na(dg), na.rm = TRUE)
if (n_g_miss > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Grade missing", Value = fmt_n(n_g_miss, N_val_erp)))

## Nodal status
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Nodal status", Value = ""))
pn <- val_erp$node_status
n_n0 <- sum(pn == "N0", na.rm = TRUE)
n_n1 <- sum(pn == "N+", na.rm = TRUE)
val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Node-negative (N0)", Value = fmt_n(n_n0, N_val_erp)),
  data.frame(Characteristics = "Node-positive", Value = fmt_n(n_n1, N_val_erp))
)

## Tumor size
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Tumor size", Value = ""))
sz <- val_erp$tumor_size_mm
sz_clean <- sz[!is.na(sz) & sz < 900]
val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Mean size, mm",
             Value = sprintf("%.1f", mean(sz_clean, na.rm = TRUE))))
n_le10 <- sum(sz_clean <= 10, na.rm = TRUE)
n_11_20 <- sum(sz_clean > 10 & sz_clean <= 20, na.rm = TRUE)
n_gt20 <- sum(sz_clean > 20, na.rm = TRUE)
n_sz_miss <- sum(is.na(sz) | sz >= 900)
val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "\u226410 mm", Value = fmt_n(n_le10, N_val_erp)),
  data.frame(Characteristics = "11-20 mm", Value = fmt_n(n_11_20, N_val_erp)),
  data.frame(Characteristics = ">20 mm", Value = fmt_n(n_gt20, N_val_erp))
)
if (n_sz_miss > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Size missing/multifocal", Value = fmt_n(n_sz_miss, N_val_erp)))

## Surgery
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Surgery", Value = ""))
tr <- val_erp$surgery_type
n_mast <- sum(tr == "mastectomy", na.rm = TRUE)
n_sect <- sum(tr == "breast_conserving", na.rm = TRUE)
if (n_mast > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Mastectomy", Value = fmt_n(n_mast, N_val_erp)))
if (n_sect > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Breast-conserving", Value = fmt_n(n_sect, N_val_erp)))

## Radiotherapy
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Radiotherapy", Value = ""))
rt <- val_erp$radiotherapy
val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(sum(rt == "yes", na.rm = TRUE), N_val_erp)),
  data.frame(Characteristics = "No", Value = fmt_n(sum(rt == "no", na.rm = TRUE), N_val_erp))
)

## Chemotherapy
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Chemotherapy", Value = ""))
chemo <- val_erp$chemotherapy
val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(sum(chemo == "yes", na.rm = TRUE), N_val_erp)),
  data.frame(Characteristics = "No", Value = fmt_n(sum(chemo == "no", na.rm = TRUE), N_val_erp))
)

## Endocrine therapy
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Endocrine therapy", Value = ""))
endo <- val_erp$endocrine_therapy
n_en_yes <- sum(endo == "yes", na.rm = TRUE)
n_en_no  <- sum(endo == "no", na.rm = TRUE)
n_en_oth <- sum(!is.na(endo) & !(endo %in% c("yes", "no")))
val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Yes", Value = fmt_n(n_en_yes, N_val_erp)),
  data.frame(Characteristics = "No", Value = fmt_n(n_en_no, N_val_erp))
)
if (n_en_oth > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "Unknown", Value = fmt_n(n_en_oth, N_val_erp)))

## Diagnosis period
val_erp_rows <- rbind(val_erp_rows, data.frame(Characteristics = "Diagnosis period", Value = ""))
dy <- val_erp$diagnosis_year
dy_full <- ifelse(!is.na(dy) & dy < 100, dy + 1900, dy)
n_85_89 <- sum(dy_full >= 1985 & dy_full <= 1989, na.rm = TRUE)
n_90_94 <- sum(dy_full >= 1990 & dy_full <= 1994, na.rm = TRUE)
n_95_99 <- sum(dy_full >= 1995 & dy_full <= 1999, na.rm = TRUE)
n_00_04 <- sum(dy_full >= 2000 & dy_full <= 2004, na.rm = TRUE)
if (n_85_89 > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "1985-1989", Value = fmt_n(n_85_89, N_val_erp)))
if (n_90_94 > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "1990-1994", Value = fmt_n(n_90_94, N_val_erp)))
if (n_95_99 > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "1995-1999", Value = fmt_n(n_95_99, N_val_erp)))
if (n_00_04 > 0) val_erp_rows <- rbind(val_erp_rows,
  data.frame(Characteristics = "2000-2004", Value = fmt_n(n_00_04, N_val_erp)))

## Set column name
colnames(val_erp_rows) <- c("Characteristics", sprintf("Cases (percentage)\nN = %d", N_val_erp))

## Save
erp_table_path <- file.path(output_dir, "Table_validation_ERpos_clinicopathological.docx")
save_clinpath_docx(val_erp_rows, erp_table_path,
  sprintf("Table S4. Distribution of clinicopathological parameters in ER-positive tumors within the validation cohort (MSPI-available subset, N = %d).", N_val_erp))


cat("\n═══════════════════════════════════════════════\n")
cat("PART 2: Complete — all clinicopathological tables saved\n")
cat("═══════════════════════════════════════════════\n")
