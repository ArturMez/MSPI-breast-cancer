##############################################################################
## Figure 1d: Exploratory cohort - ECDF Proximity Curves (M1 & M2 to CD3)
##
## Generates empirical cumulative distribution functions (ECDFs) showing
## nearest-neighbor distances from macrophage subtypes (M1, M2) to CD3+ T cells.
## Includes KS test statistics for density-independent comparison.
## Pre-computed ECDF at 500 evaluation points for efficient PDF rendering.
##############################################################################

library(data.table)
library(reshape2)
library(dplyr)
library(ggplot2)

options(digits = 8)

FONT_FAM <- "Helvetica"

## ── Paths ────────────────────────────────────────────────────────────────
## Adjust these to point to your local exploratory cohort data and output directories
data_dir <- "./data/exploratory/"
out_dir  <- "./output/"

## ── Colour scheme ────────────────────────────────────────────────────────
COL_NEG  <- "darkolivegreen4"
COL_M2   <- "coral1"
COL_M1   <- "deepskyblue2"

## ── Files to process ─────────────────────────────────────────────────────
## Nearest-neighbor distance files for M1 and M2 macrophages to CD3+ T cells
runs <- list(
  list(file  = "nnALL_M2_to_CD3_exploratory.txt",
       label = "M2",
       col   = COL_M2),
  list(file  = "nnALL_M1_to_CD3_exploratory.txt",
       label = "M1",
       col   = COL_M1)
)

## ── Helper: format p-value ───────────────────────────────────────────────
fmt_p <- function(pv) {
  if (pv < 0.0001) return("< 0.0001")
  formatC(pv, format = "f", digits = 4)
}

## ── Helper: pre-compute ECDF as a light data.frame ───────────────────────
## Evaluates the empirical CDF at n_pts evenly spaced x values.
## Returns a tidy data.frame ready for geom_line.
ecdf_line <- function(values, group_label, n_pts = 500) {
  fn <- ecdf(values)
  xs <- seq(0, 50, length.out = n_pts)
  data.frame(x = xs, y = fn(xs), Group = group_label,
             stringsAsFactors = FALSE)
}

## ══════════════════════════════════════════════════════════════════════════
## Main loop
## ══════════════════════════════════════════════════════════════════════════
for (run in runs) {

  cat("\n=== Processing", run$label, "===\n")

  ## ── 1. Read & cap at 50 um ─────────────────────────────────────────────
  nn <- data.table::fread(file.path(data_dir, run$file),
                          sep = "\t", header = TRUE, fill = TRUE)

  nn$POS <- ifelse(nn$POS < 50, nn$POS, NA)
  nn$neg <- ifelse(nn$neg < 50, nn$neg, NA)
  nn <- na.omit(as.data.frame(nn))

  ## ── 2. Compute ECDF + KS test on full data ────────────────────────────
  vals_neg <- nn$neg
  vals_pos <- nn$POS

  ks <- ks.test(vals_neg, vals_pos)
  ks_D <- round(ks$statistic, 4)
  ks_p <- fmt_p(ks$p.value)

  n1 <- as.numeric(length(vals_neg))
  n2 <- as.numeric(length(vals_pos))
  critical_D <- round(1.36 * sqrt((n1 + n2) / (n1 * n2)), 4)

  med_neg <- quantile(vals_neg, 0.50)
  med_pos <- quantile(vals_pos, 0.50)

  cat("  n neg:", n1, "  n POS:", n2, "\n")
  cat("  Median neg:", round(med_neg, 4), "\n")
  cat("  Median POS:", round(med_pos, 4), "\n")
  cat("  KS D =", ks_D, " p =", ks_p, "\n")

  ## ── 3. Pre-compute lightweight ECDF lines ──────────────────────────────
  nice_neg <- "Random M\u03c6 pool"
  nice_pos <- paste0(run$label, " macrophages")

  ecdf_df <- rbind(ecdf_line(vals_neg, nice_neg),
                   ecdf_line(vals_pos, nice_pos))
  ecdf_df$Group <- factor(ecdf_df$Group, levels = c(nice_neg, nice_pos))

  ## ── 4. Plot ────────────────────────────────────────────────────────────
  p <- ggplot(ecdf_df, aes(x = x, y = y, colour = Group)) +
    geom_line(linewidth = 0.7) +
    coord_cartesian(xlim = c(0, 50), ylim = c(0, 1.00)) +

    ## median dashed lines
    geom_vline(xintercept = med_neg, linetype = "dashed",
               colour = COL_NEG, linewidth = 0.4) +
    geom_vline(xintercept = med_pos, linetype = "dashed",
               colour = run$col, linewidth = 0.4) +

    ## median value annotations (offset vertically to avoid overlap)
    annotate("text", x = med_neg, y = 0.56,
             label = paste0(round(med_neg, 2), " \u00b5m"),
             hjust = -0.08, colour = COL_NEG,
             size = 3.2, family = FONT_FAM) +
    annotate("text", x = med_pos, y = 0.44,
             label = paste0(round(med_pos, 2), " \u00b5m"),
             hjust = -0.08, colour = run$col,
             size = 3.2, family = FONT_FAM) +

    ## KS test annotation
    annotate("text", x = 2, y = 0.78,
             label = paste0("KS test\nD = ", ks_D,
                            "\nCritical D = ", critical_D,
                            "\np = ", ks_p),
             hjust = 0, colour = "black",
             size = 3.5, family = FONT_FAM) +

    scale_colour_manual(values = setNames(c(COL_NEG, run$col),
                                          c(nice_neg, nice_pos))) +
    labs(x = "Nearest-neighbor distance to CD3+ T cell (\u00b5m)",
         y = "Cumulative fraction",
         title = paste0("Exploratory cohort - ",
                        run$label, " to CD3 proximity"),
         colour = NULL) +
    theme_minimal(base_size = 13, base_family = FONT_FAM) +
    theme(
      legend.position   = c(0.75, 0.20),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.text       = element_text(size = 11),
      plot.title        = element_text(face = "bold", size = 13),
      axis.title        = element_text(size = 12),
      axis.text         = element_text(size = 11),
      plot.margin       = margin(10, 12, 10, 10)
    )

  ## ── 5. Save PDF (base pdf device - best Illustrator compatibility) ─────
  fname <- paste0("exploratory_ECDF_", run$label, "_to_CD3.pdf")
  ggsave(file.path(out_dir, fname), p,
         width = 3.0, height = 3.0, device = pdf, useDingbats = FALSE)
  cat("  Saved:", fname, "\n")
}

cat("\nDone.\n")
