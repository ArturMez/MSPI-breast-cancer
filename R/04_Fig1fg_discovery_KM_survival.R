########################
## Figure 1f-g: Discovery cohort - KM survival curves by affinity
##
## Performs survival analysis of M1 and M2 macrophage-to-T cell affinity
## in the discovery cohort using optimal cutpoints and KM estimation.
## Generates publication-quality forest plots with risk tables.
########################

## ── Libraries ──────────────────────────────────────────────────────────────
require(gdata)
require(xlsx)
require(ggplot2)
library(reshape)
library(survival)
library(tidyverse)
library(schoRsch)
library(DescTools)
library(foreign)
library(survminer)
library(ggpubr)

## ── Paths ──────────────────────────────────────────────────────────────────
## Adjust these to point to your local discovery cohort data and output directories
data_dir   <- "./data/discovery/"

## Input files may carry the names used in this repository or their original names
## (as shared on request); the first file found is used.
find_input <- function(dir, ...) {
  for (f in c(...)) {
    p <- file.path(dir, f)
    if (file.exists(p)) return(p)
  }
  stop("Cannot find any of: ", paste(c(...), collapse = ", "), " in ", dir)
}
output_dir <- "./output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

## ── v6-style plot settings ─────────────────────────────────────────────────
BASE_SIZE  <- 14
KM_BASE    <- 13
FONT_FAM   <- "Arial"

COL_GREEN  <- "#6aad8e"
COL_RED    <- "#c47a7a"
PAL_2G     <- c(COL_GREEN, COL_RED)

## Detect quartz fallback for macOS when cairo_pdf is unavailable
USE_QUARTZ <- (.Platform$OS.type != "windows") &&
  tryCatch({ cairo_pdf(tempfile()); dev.off(); FALSE },
           error = function(e) TRUE)

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

## ── v6-style helpers ───────────────────────────────────────────────────────
save_pdf <- function(fpath, w, h) {
  if (USE_QUARTZ) {
    quartz(type = "pdf", file = fpath, width = w, height = h,
           family = FONT_FAM, bg = "white")
  } else {
    cairo_pdf(fpath, width = w, height = h, bg = "white")
  }
}

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
  cat("  Saved:", fpath, "\n")
}

## ── Load clinical data ─────────────────────────────────────────────────────
## Expected: SPSS file with clinical variables, endpoint data
dat <- NULL
disc_spss <- find_input(data_dir, "discovery_clinical.sav", "discovery_clinical_updated.sav",
                        "BCY_workingr_832_updated_211024.sav")
dat <- read.spss(disc_spss,
                 use.value.labels = T, to.data.frame = T,
                 max.value.labels = Inf, trim.factor.names = T,
                 trim_values = TRUE, reencode = NA)
## The original clinical file calls the subtype variable Subtype_BCY21
if (!"Subtype_Discovery21" %in% names(dat) && "Subtype_BCY21" %in% names(dat)) {
  names(dat)[names(dat) == "Subtype_BCY21"] <- "Subtype_Discovery21"
}

## ── Load distance and density data ─────────────────────────────────────────
## Affinity scores and CD3 density computed from spatial analysis
distance_M1 <- data.table::fread(find_input(data_dir, "discovery_distances_mean_M1_to_CD3.txt", "BCY_distances_mean_M1_to_CD3.txt"),
                                  sep = "\t", header = T, fill = T)
distance_M2 <- data.table::fread(find_input(data_dir, "discovery_distances_mean_M2_to_CD3.txt", "BCY_distances_mean_M2_to_CD3.txt"),
                                  sep = "\t", header = T, fill = T)

CD3_density_data <- data.table::fread(find_input(data_dir, "discovery_immune_variables_with_CD3.txt",
                                                  "BCY_immune_variables_with_CD3_for spatial.txt"),
                                       sep = "\t", header = T, fill = T)
CD3_density_data <- CD3_density_data[, c(1, 2)]

## ── Filter to MSPI-available subset (both M1 and M2 affinity available) ───
## This ensures all analyses use the same patient set
mspi_ids <- intersect(distance_M1$AnnUnit, distance_M2$AnnUnit)
cat("MSPI-available filter: M1 n=", nrow(distance_M1),
    ", M2 n=", nrow(distance_M2),
    ", both n=", length(mspi_ids), "\n")
distance_M1 <- distance_M1[distance_M1$AnnUnit %in% mspi_ids, ]
distance_M2 <- distance_M2[distance_M2$AnnUnit %in% mspi_ids, ]
cat("After filter: M1 n=", nrow(distance_M1),
    ", M2 n=", nrow(distance_M2), "\n")

## ── Recode inflammatory variable ───────────────────────────────────────────
dat$Re_Inflammatory_1 <- as.character(dat$Re_Inflammatory_1)
dat$Re_Inflammatory_1 <- ifelse((dat$Re_Inflammatory_1 == "Missing"), NA, dat$Re_Inflammatory_1)
dat$Re_Inflammatory_1 <- ifelse((dat$Re_Inflammatory_1 == "None or mild"), "1", dat$Re_Inflammatory_1)
dat$Re_Inflammatory_1 <- ifelse((dat$Re_Inflammatory_1 == "Moderate"), "2", dat$Re_Inflammatory_1)
dat$Re_Inflammatory_1 <- ifelse((dat$Re_Inflammatory_1 == "Strong"), "3", dat$Re_Inflammatory_1)
dat$Re_Inflammatory_1 <- as.factor(dat$Re_Inflammatory_1)
dat$Re_Inflammatory_1 <- relevel(dat$Re_Inflammatory_1, ref = "1")

quant_cut <- function(x, n) {
  qs <- quantile(x, 1:(n - 1) / n)
  brks <- c(-Inf, qs, Inf)
  cut(x, breaks = brks, labels = FALSE)
}

## ── TIL variable ───────────────────────────────────────────────────────────
dat$TIL <- dat$CD4 + dat$CD8 + dat$CD20
xx <- NULL
xx <- dat[, c(603, 637)]
xx <- na.omit(xx)
xx$TIL <- quant_cut(xx$TIL, 3)
dat <- left_join(dat, xx, by = "AnnUnit")

xx <- NULL
xx <- dat[, c(603, 590, 637)]
xx <- na.omit(xx)
xx$TIL.4 <- quant_cut((xx$CD8_T / (xx$CD8 + 1)), 3)
xx <- xx[, c(1, 4)]
dat <- left_join(dat, xx, by = "AnnUnit")

markers <- colnames(dat[, c("Re_Inflammatory_1", "TIL.y", "CD4_cat", "CD8_cat",
                             "CD20_cat", "M1_cat", "M2_cat", "Myeloid_cat",
                             "SIA_cat", "TIL.4")])

## ── Define survival endpoint (DDFS) ───────────────────────────────────────
dat$OS <- as.numeric(dat$futime2021)
dat$status <- (dat$BCcens2021)
dat$SurvObj <- with(dat, Surv(dat$OS, status == "Breast cancer death"))

dat$EVENT_DR_or_BC_death <- as.character(dat$BCcens2021)
dat$EVENT_DR_or_BC_death <- ifelse((dat$DDFScens2021 == "Distant recurrence"),
                                    "DDFScens2021", as.character(dat$BCcens2021))
dat$status <- dat$EVENT_DR_or_BC_death
dat$SurvObj <- with(dat, Surv(dat$DDFStime2021, status == "DDFScens2021"))

## ── Exclude unclassified ───────────────────────────────────────────────────
dat <- dat[which(dat$Subtype_Discovery21 != "Unclassified"), ]

## ── Save intermediate table ────────────────────────────────────────────────
write.table(dat, file = file.path(output_dir, "discovery_temp.txt"),
            append = FALSE, sep = "\t", row.names = F, col.names = TRUE)


############################################################################
## ── ANALYSIS FUNCTION: optimal cutpoint KM for a given affinity score ────
############################################################################

run_affinity_km <- function(dat_base, distance_data, affinity_label) {

  cat("\n========================================\n")
  cat("  Running:", affinity_label, "\n")
  cat("========================================\n")

  ## Join distance data
  dat_work <- left_join(dat_base, distance_data, by = "AnnUnit")
  dat_work <- left_join(dat_work, CD3_density_data, by = "AnnUnit")

  ## Reference dataset (full cohort, no subgroup filtering)
  dat_ref <- dat_work

  ## Correlation test
  df_corr <- dat_ref %>%
    dplyr::select(any_of(c("M1", "M2")), delta_CD3) %>%
    tidyr::drop_na()

  if ("M1" %in% colnames(df_corr) && affinity_label == "M1-to-T cell affinity") {
    spearman <- cor.test(df_corr$M1, df_corr$delta_CD3, method = "spearman")
    cat("  Spearman (M1 vs delta_CD3): rho =", round(spearman$estimate, 3),
        ", p =", format.pval(spearman$p.value, digits = 3), "\n")
  }
  if ("M2" %in% colnames(df_corr) && affinity_label == "M2-to-T cell affinity") {
    spearman <- cor.test(df_corr$M2, df_corr$delta_CD3, method = "spearman")
    cat("  Spearman (M2 vs delta_CD3): rho =", round(spearman$estimate, 3),
        ", p =", format.pval(spearman$p.value, digits = 3), "\n")
  }

  ## Prepare for survival analysis
  data <- dat_ref[!is.na(dat_ref$delta_CD3), ]
  d <- data

  ## Recode status for surv_cutpoint
  d$status.x <- ifelse((d$status == "Censored"), 0, d$status)
  d$status.x <- ifelse((d$status == "DDFScens2021"), 1, d$status.x)
  d$status <- as.numeric(d$status.x)

  ## Optimal cutpoint
  res.cut <- surv_cutpoint(d, time = "DDFStime2021", event = "status",
                            variables = c("delta_CD3"))
  opt_cut <- res.cut$cutpoint$cutpoint
  cat("  Optimal cutpoint:", opt_cut, "\n")

  res.cat_2 <- surv_categorize(res.cut)

  ## high delta_CD3 = farther = low affinity; low delta_CD3 = closer = high affinity
  ## XXX = 1 means low delta_CD3 (high affinity), XXX = 0 means high delta_CD3 (low affinity)
  res.cat_2$XXX <- ifelse((res.cat_2$delta_CD3 == "high"), 0, 1)

  n_hi <- sum(res.cat_2$XXX == 1)
  n_lo <- sum(res.cat_2$XXX == 0)

  res.cat_2$SurvObj <- with(res.cat_2, Surv(res.cat_2$DDFStime2021, status == "1"))
  cx <- coxph(SurvObj ~ XXX, data = res.cat_2)
  p_value_optimal <- summary(cx)$coefficients[5]

  cat("  Cox HR:", round(summary(cx)$conf.int[1, 1], 3),
      " (", round(summary(cx)$conf.int[1, 3], 3), "-",
      round(summary(cx)$conf.int[1, 4], 3), ")\n")
  cat("  p-value:", format.pval(p_value_optimal, digits = 3), "\n")

  ## KM plot - v6 style
  fit <- survfit(SurvObj ~ XXX, data = res.cat_2, conf.type = "log-log")

  ## Determine group labels based on affinity type
  if (grepl("M2", affinity_label)) {
    legend_labs <- c(
      paste0("Low M2 affinity (n=", n_lo, ")"),
      paste0("High M2 affinity (n=", n_hi, ")")
    )
    plot_title <- "M2-to-T Cell Affinity | Discovery Cohort"
    file_prefix <- "disc_M2_affinity"
  } else {
    legend_labs <- c(
      paste0("Low M1 affinity (n=", n_lo, ")"),
      paste0("High M1 affinity (n=", n_hi, ")")
    )
    plot_title <- "M1-to-T Cell Affinity | Discovery Cohort"
    file_prefix <- "disc_M1_affinity"
  }

  km <- ggsurvplot(
    fit, data = res.cat_2,
    pval = TRUE, pval.method = TRUE, conf.int = TRUE,
    risk.table = TRUE,
    xlab = "DDFS (months)",
    title = plot_title,
    subtitle = paste0("Optimal cutpoint: ", round(opt_cut, 4)),
    legend.labs = legend_labs,
    palette = PAL_2G,
    ggtheme = theme_classic2(base_size = KM_BASE, base_family = FONT_FAM),
    risk.table.height = 0.25
  )

  save_km(km, paste0(file_prefix, "_optimal_cutpoint_KM"))

  ## Print summary
  cat("\n  --- Summary for", affinity_label, "---\n")
  print(summary(cx))
  cat("\n")

  invisible(list(cox = cx, fit = fit, cutpoint = res.cat_2))
}


############################################################################
## ── Run both analyses ────────────────────────────────────────────────────
############################################################################

cat("\n\n################################################################\n")
cat("## Discovery Cohort Survival Analysis - M1 and M2 Affinity\n")
cat("################################################################\n")

## M1-to-T cell affinity
res_M1 <- run_affinity_km(dat, distance_M1, "M1-to-T cell affinity")

## M2-to-T cell affinity
res_M2 <- run_affinity_km(dat, distance_M2, "M2-to-T cell affinity")

cat("\n\nDone. Outputs saved to:", output_dir, "\n")
cat("KM plots saved to:", plots_dir, "\n")
