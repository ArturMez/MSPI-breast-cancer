########################
## 11 - Figure 3, panel d (pooled individual patient data analysis)
##
## Draws the marks and HR labels of Figure 3d at the page positions, scales, line widths, marker sizes, colours and
## fonts of the figure layout (measured from the vector content of the figure page). The frame, grid lines, axis and
## row labels and the other panels come from assets/Fig3_template.pdf in script 12. Panels a-c use the validation
## cohort only.
##
## Numbers: MSPI HR per within-cohort SD for the discovery young early-stage luminal patients, the validation
##          ER-positive patients aged <= 50 years and the pooled stratified Cox model, from Pool_B_T4_forest_<date>.csv
##          (script 05, task 4).
## Inputs : ./output/<date>_Integrated_MSPI_Analysis/tables/
## Outputs: ./output/fig3/Fig3_panels.pdf   one transparent A4 page (drawn 16 times larger), placed onto the
##                                           template by script 12
##          ./output/fig3/Fig3_numbers.tsv   every number drawn
## Usage  : Rscript R/11_Fig3d_draw.R
########################

suppressMessages({ library(data.table); library(grid); library(systemfonts) })
out_dir <- "./output/fig3"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
latest <- function(pat) { f <- sort(list.files("./output", pattern = pat, full.names = TRUE)); f[length(f)] }
TAB <- file.path(latest("_Integrated_MSPI_Analysis$"), "tables")
pf <- list.files(TAB, pattern = "^Pool_B_T4_forest_[0-9]{8}\\.csv$", full.names = TRUE)[1]
cat("pooled forest table:", pf, "\n")
fo <- fread(pf)
stopifnot(nrow(fo) == 3, fo$Label[3] == "Pooled (stratified Cox)")

## ------------------------------------------------------------------ page, units and drawing helpers (as script 08)
S <- 16; PW <- 595.28; PH <- 841.89
X  <- function(x) unit(x * S, "bigpts")
Y  <- function(y) unit((PH - y) * S, "bigpts")
LW <- function(pt) pt * S * 96 / 72
seg <- function(x0, y0, x1, y1, col = "black", w = 0.5)
  grid.segments(X(x0), Y(y0), X(x1), Y(y1), gp = gpar(col = col, lwd = LW(w), lineend = "butt"))
## Text at a baseline origin with the published tracking: the labels are Arial with every advance shortened by
## 0.168 pt (Illustrator tracking -42/1000 em at 4 pt), measured from the published glyph positions. Each part is
## list(text, fontface); the letters are placed one by one.
tracked <- function(parts, x, y, size = 4, col = "black", track = -0.168) {
  for (pt in parts) for (ch in strsplit(pt[[1]], "")[[1]]) {
    if (ch != " ") grid.text(ch, X(x), Y(y), hjust = 0, vjust = 0, gp = gpar(fontsize = size * S, col = col, fontface = pt[[2]], fontfamily = "Arial"))
    adv <- systemfonts::shape_string(ch, family = "Arial", size = 1000, weight = if (pt[[2]] %in% c(2, 4)) "bold" else "normal",
                                     italic = pt[[2]] %in% c(3, 4))$metrics$width / 1000 * size    # (integer pixels at 1000 pt)
    x <- x + adv + track
  }
}

## ------------------------------------------------------------------ panel d
## HR axis of the published panel: x = 294.12 + 14.13 * (HR - 1), frame from HR 1.0 to 3.70; rows at y 36.72
## (discovery), 47.46 (validation), 58.21 (pooled). Error bars 0.5 pt #4D4D4D with caps 2.37 pt high, cohort
## estimates as circles (3.12 and 1.92 pt, sized by events as published), pooled estimate as a 12.0 x 4.26 pt diamond.
TX <- "#4D4D4D"; RED <- "#C47A7A"; POOL <- "#C0392B"
hx <- function(v) 294.12 + 14.13 * (v - 1)
if (min(hx(fo$CI_lo)) < 294.08 || max(hx(fo$CI_hi)) > 332.24) stop("a confidence interval is outside the published axis")
ROW <- c(36.72, 47.46, 58.21); BASE <- c(37.47, 48.22, 58.97)

cairo_pdf(file.path(out_dir, "Fig3_panels.pdf"), width = PW * S / 72, height = PH * S / 72, family = "Arial", bg = "transparent")
grid.newpage()
for (i in 1:3) {
  seg(hx(fo$CI_lo[i]), ROW[i], hx(fo$CI_hi[i]), ROW[i], TX)
  for (e in c(fo$CI_lo[i], fo$CI_hi[i])) seg(hx(e), ROW[i] - 1.185, hx(e), ROW[i] + 1.185, TX)
}
for (i in 1:2) grid.circle(X(hx(fo$HR[i])), Y(ROW[i]), r = unit(c(3.12, 1.92)[i] / 2 * S, "bigpts"), gp = gpar(fill = RED, col = NA))
grid.polygon(X(hx(fo$HR[3]) + c(-6, 0, 6, 0)), Y(ROW[3] + c(0, -2.13, 0, 2.13)), gp = gpar(fill = POOL, col = NA))
nums <- list()
for (i in 1:3) {
  b <- sprintf("HR=%.2f", fo$HR[i])
  r <- sprintf(" (%.2f-%.2f)  n=%d ev=%d  p=%.2e", fo$CI_lo[i], fo$CI_hi[i], fo$n[i], fo$events[i], fo$p[i])
  tracked(list(list(b, 4), list(r, 3)), 335.77, BASE[i], 4, TX)
  nums[[i]] <- data.table(panel = "d", row = fo$Label[i], label = paste0(b, r))
}
invisible(dev.off())
fwrite(rbindlist(nums), file.path(out_dir, "Fig3_numbers.tsv"), sep = "\t")
cat("written:", file.path(out_dir, "Fig3_panels.pdf"), "and Fig3_numbers.tsv\n")
