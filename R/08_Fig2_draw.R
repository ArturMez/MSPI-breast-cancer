########################
## 08 - Figure 2, data panels
##
## Draws every data element of Figure 2 at the page positions, scales, line widths, marker sizes, colours,
## opacities and fonts of the figure layout (measured from the vector content of the figure page). Fixed
## elements (panel letters, cohort label, cell schematic) come from assets/Fig2_template.pdf in script 09.
##
## Numbers: C-indices with bootstrap CI (c, f, g), MSPI HR per SD by subtype (d), clinical multivariable Cox
## (e, h) and MSPI HR per SD UV and density-adjusted (f, g) are read from the tables of scripts 05 and 06.
## The scatter (a), boxplots and Dunn's tests with BH correction (b), Kaplan-Meier curves by MSPI quartile with
## the 4th vs 1st log-rank p (f, g) and the Type II Anova pies (f, g) are computed here with the data preparation
## of script 05. MSPI quartiles are cut within each subgroup as in script 05 (revision, 2026-11-15: the qMSPI
## curves of the submitted panels f and g are in Supplementary Figure S4, script 15).
##
## Inputs : ./data/discovery/ (SAV, affinity scores, densities), ./output/<date>_Integrated_MSPI_Analysis/,
##          ./output/<date>_MSPI_clinical/
## Outputs: ./output/fig2/Fig2_panels.pdf   one transparent A4 page (drawn 16 times larger), placed onto the
##                                           template by script 09
##          ./output/fig2/Fig2_numbers.tsv   every number drawn in the figure
## Usage  : Rscript R/08_Fig2_draw.R
########################

suppressMessages({
  library(foreign); library(data.table); library(survival); library(car); library(rstatix); library(grid)
})
disc_dir <- "./data/discovery"
out_dir  <- "./output/fig2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
latest <- function(pat) { f <- sort(list.files("./output", pattern = pat, full.names = TRUE)); f[length(f)] }
TAB <- file.path(latest("_Integrated_MSPI_Analysis$"), "tables")
PLT <- file.path(latest("_Integrated_MSPI_Analysis$"), "plots")
CLI <- latest("_MSPI_clinical$")
rt  <- function(stem) fread(list.files(TAB, pattern = paste0("^", stem, "_[0-9]{8}\\.csv$"), full.names = TRUE)[1])
inp <- c(m1 = "discovery_distances_mean_M1_to_CD3_nucleus.txt", m2 = "discovery_distances_mean_M2_to_CD3_nucleus.txt",
         dens = "discovery_immune_variables_with_CD3_nucleus.txt")
cat("tables:", TAB, "\n")

## ------------------------------------------------------------------ data (as script 05, steps 2 to 4)
sav <- list.files(disc_dir, pattern = "workingr_832_updated_211024\\.sav$", full.names = TRUE)[1]
cl <- suppressWarnings(read.spss(sav, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
                                 trim.factor.names = TRUE, trim_values = TRUE, reencode = NA))
if (!"Subtype_BCY21" %in% names(cl)) names(cl)[names(cl) == "Subtype_Discovery21"] <- "Subtype_BCY21"
cl$status_num <- ifelse(as.character(cl$DDFScens2021) == "Distant recurrence" | as.character(cl$BCcens2021) == "Breast cancer death", 1L, 0L)
cl$DDFStime2021 <- as.numeric(cl$DDFStime2021)
cl <- cl[as.character(cl$Subtype_BCY21) != "Unclassified", ]
cl <- as.data.table(cl)[, .(AnnUnit, status_num, DDFStime2021, Subtype = as.character(Subtype_BCY21),
                            Stage = as.character(Stage_TNM_44), Age = as.character(Age_cat_2), M1d = M1, M2d = M2)]
m1 <- fread(file.path(disc_dir, inp[["m1"]])); m2 <- fread(file.path(disc_dir, inp[["m2"]]))
cd <- fread(file.path(disc_dir, inp[["dens"]]))[, 1:2]
setnames(cd, c("AnnUnit", "Td"))
d <- merge(cl, m1[, .(AnnUnit, M1s = -delta_CD3)], by = "AnnUnit", all.x = TRUE)
d <- merge(d, m2[, .(AnnUnit, M2s = -delta_CD3)], by = "AnnUnit", all.x = TRUE)
d <- merge(d, cd, by = "AnnUnit", all.x = TRUE)
s4 <- d[!is.na(M1s) & !is.na(M2s) & !is.na(DDFStime2021) & DDFStime2021 > 0]
s4[, MSPI := M2s - M1s]
s4[, Sub3 := factor(fcase(Subtype %in% c("Lum A", "Lum B"), "L", Subtype %in% c("Lum Her2", "Her2 exp"), "H",
                          Subtype == "Triple neg", "T"), levels = c("L", "H", "T"))]
early <- s4[Stage %in% c("Stage I", "Stage IIa", "Stage IIb") & Subtype %in% c("Lum A", "Lum B")]
young <- early[Age == "<40"]
cat("S4", nrow(s4), sum(s4$status_num), "| early", nrow(early), sum(early$status_num), "| young", nrow(young), sum(young$status_num), "\n")

## ------------------------------------------------------------------ page, units and drawing helpers
S <- 16; PW <- 595.28; PH <- 841.89                         # A4 in pt, drawn S times larger
X  <- function(x) unit(x * S, "bigpts")
Y  <- function(y) unit((PH - y) * S, "bigpts")
LW <- function(pt) pt * S * 96 / 72                          # grid lwd (1/96 inch) for a width in pt
lin <- function(v, lo, hi, p0, p1) p0 + (v - lo) / (hi - lo) * (p1 - p0)
seg <- function(x0, y0, x1, y1, col = "black", w = 0.5, lty = "solid", alpha = 1)
  grid.segments(X(x0), Y(y0), X(x1), Y(y1), gp = gpar(col = col, lwd = LW(w), lty = lty, alpha = alpha, lineend = "butt"))
box <- function(x0, y0, x1, y1, col = NA, fill = NA, w = 0.25, lty = "solid", alpha = 1)
  grid.rect(X(x0), Y(y1), width = unit((x1 - x0) * S, "bigpts"), height = unit((y1 - y0) * S, "bigpts"), just = c("left", "bottom"),
            gp = gpar(col = col, fill = fill, lwd = LW(w), lty = lty, alpha = alpha, linejoin = "mitre"))
ellipses <- function(cx, cy, w, h, fill = NA, col = NA, lw = 0.25, alpha = 1, n = 28) {
  if (!length(cx)) return(invisible())
  t <- seq(0, 2 * pi, length.out = n + 1)[-1]
  k <- length(cx)
  grid.polygon(X(rep(cx, each = n) + rep(w / 2 * cos(t), k)), Y(rep(cy, each = n) + rep(h / 2 * sin(t), k)), id.lengths = rep(n, k),
               gp = gpar(fill = fill, col = col, lwd = LW(lw), alpha = alpha))
}
diamonds <- function(cx, cy, s, fill) {
  if (!length(cx)) return(invisible())
  k <- length(cx); dx <- c(0, s / 2, 0, -s / 2); dy <- c(-s / 2, 0, s / 2, 0)
  grid.polygon(X(rep(cx, each = 4) + rep(dx, k)), Y(rep(cy, each = 4) + rep(dy, k)), id.lengths = rep(4, k),
               gp = gpar(fill = fill, col = NA))
}
## Text at a baseline origin (x, y) as in the published PDF; hjust 0 left, 0.5 centred, 1 right.
## The published labels are set narrower than plain Arial (Illustrator tracking and scaling). Each fixed label is
## compressed along its baseline to its published ink width (WT, keyed by label and rounded x), and numbers by
## 0.965, through a grid group transform, so that height and font size stay as published.
WT <- c("qMSPI subgroups@35" = 41.04, "M1-to-T affinity score@39" = 35.28, "M2-to-T affinity score@23" = 35.28,
        "Quadrant@83" = 17.27, "Q1: M1-close/M2-far@92" = 37.82, "Q2: both close@92" = 25.70, "Q3: both far@92" = 21.05,
        "Q4: M2-close/M1-far@92" = 37.87, "Status@83" = 11.54, "Censored@92" = 17.41, "Event@92" = 10.14,
        "T cell@165" = 8.72, "M1@186" = 5.39, "M2@203" = 5.39, "MSPI@228" = 9.28, "Dunn's test with BH correction@176" = 55.07,
        "C-index@143" = 12.44, "with bootstrap 95% CI@148" = 35.26, "MSPI@356" = 9.68, "MSPI+density@356" = 24.44,
        "Full cohort@271" = 16.99, "Luminal@274" = 13.40, "HER2+@274" = 12.70, "TNBC@276" = 10.74,
        "Hazard Ratio per SD (95% CI)@304" = 61.57, "Hazard Ratio (95% CI)@466" = 36.87, "Luminal (DDFS)@405" = 28.20,
        "Hazard Ratio (95% CI)@464" = 36.87, "Young Early Luminal (DDFS)@392" = 50.67,
        "MSPI (per SD)@454" = 25.27, "Age ≥ 40@454" = 15.58, "Stage IIa@454" = 15.46, "Stage IIb@454" = 15.46,
        "Stage III-IV@454" = 19.48, "Chemo- and Endocrine treatment@454" = 57.45, "Endocrine treatment@454" = 35.06,
        "No systemic treatment@454" = 39.36, "Diagnosis period: 1998-2002@454" = 47.67, "Diagnosis period: 2003-2005@454" = 47.67,
        "MSPI (per SD)@445" = 23.94, "Stage IIa@445" = 14.64, "Stage IIb@445" = 14.64, "Chemo- and Endocrine treatment@445" = 54.15,
        "Endocrine treatment@445" = 32.92, "No systemic treatment@445" = 36.35, "Diagnosis period: 1998-2002@445" = 46.66,
        "Diagnosis period: 2003-2005@445" = 46.66, "Early Luminal (DDFS)@60" = 47.76, "Young Early Luminal (DDFS)@238" = 63.33,
        "Survival probability@26" = 33.56, "Survival probability@208" = 33.56, "DDFS (years)@69" = 30.76, "DDFS (years)@253" = 30.76,
        "qMSPI@64" = 15.28, "qMSPI@248" = 15.28, "Q1@49" = 5.44, "Q2@49" = 5.44, "Q3@49" = 5.44, "Q4@49" = 5.44,
        "Q1@232" = 5.44, "Q2@232" = 5.44, "Q3@232" = 5.44, "Q4@232" = 5.44, "UV@90" = 5.58, "MV@90" = 6.02, "UV@275" = 5.58,
        "MV@275" = 6.02, "MSPI@93" = 12.22, "MSPI@277" = 12.23, "HR per SD@106" = 24.17, "HR per SD@291" = 24.17,
        "MSPI@177" = 9.84, "d M1@177" = 8.95, "d M2@177" = 8.95, "d T cell@177" = 13.03, "MSPI@358" = 9.84, "d M1@358" = 8.95,
        "d M2@358" = 8.95, "d T cell@358" = 13.03, "C-index with bootstrap 95% CI@187" = 55.1, "C-index with bootstrap 95% CI@369" = 55.1)
ink_w <- function(label, size, face) {
  m <- systemfonts::shape_string(label, family = "Arial", size = size * S, weight = if (face == 2) "bold" else "normal", italic = face == 3)$metrics
  (m$width - m$left_bearing - m$right_bearing) / S
}
gid <- 0
txt <- function(label, x, y, size = 4, col = "black", face = 1, hjust = 0, rot = 0, family = "Arial") {
  label <- as.character(label)
  gp <- gpar(fontsize = size * S, col = col, fontface = face, fontfamily = family)
  key <- paste0(label, "@", round(x))
  sx <- if (!is.na(WT[key])) WT[[key]] / ink_w(label, size, face) else if (grepl("^[-0-9.=pn, e]+$", label)) .965 else 1
  if (abs(sx - 1) < .005 || !(rot %in% c(0, 90))) {
    grid.text(label, X(x), Y(y), hjust = hjust, vjust = 0, rot = rot, gp = gp); return(invisible())
  }
  wd <- convertWidth(grobWidth(textGrob(label, gp = gp)), "bigpts", TRUE)          # advance width, device units
  gid <<- gid + 1; nm <- paste0("t", gid); hh <- size * S * 2
  if (rot == 0) {
    x0 <- x * S - hjust * wd * sx
    pushViewport(viewport(x = unit(x0, "bigpts"), y = Y(y), width = unit(wd, "bigpts"), height = unit(hh, "bigpts"), just = c("left", "bottom")))
    grid.define(textGrob(label, x = 0, y = 0, hjust = 0, vjust = 0, gp = gp), name = nm); popViewport()
    pushViewport(viewport(x = unit(x0, "bigpts"), y = Y(y), width = unit(wd * sx, "bigpts"), height = unit(hh, "bigpts"), just = c("left", "bottom")))
  } else {
    pushViewport(viewport(x = X(x), y = Y(y), width = unit(hh, "bigpts"), height = unit(wd, "bigpts"), just = c("left", "bottom")))
    grid.define(textGrob(label, x = 0, y = 0, hjust = 0, vjust = 0, rot = 90, gp = gp), name = nm); popViewport()
    pushViewport(viewport(x = X(x), y = Y(y), width = unit(hh, "bigpts"), height = unit(wd * sx, "bigpts"), just = c("left", "bottom")))
  }
  grid.use(nm); popViewport()
}
TX <- "#4D4D4D"; INK <- "#333333"; GRID <- "#EBEBEB"; BAR <- "#4D4D4D"
QC <- c(Q1 = "#6AAD8E", Q2 = "#D4B86A", Q3 = "#A0A0A0", Q4 = "#C47A7A")
MQC <- c(Q1 = "#5B7FA0", Q2 = "#8EABBF", Q3 = "#C3ABAD", Q4 = "#C47A7A")   # MSPI quartiles, low (blue) to high (red)
nums <- list()
addn <- function(panel, what, value) nums[[length(nums) + 1]] <<- data.table(panel = panel, quantity = what, value = value)
p_txt <- function(p) if (p < .001) sprintf("%.1e", p) else sprintf("%#.2g", p)
ci_txt <- function(h, l, u, p) sprintf("%.2f (%.2f-%.2f) p=%.2g", h, l, u, p)
quad <- function(x, c1, c2) factor(fifelse(x$M1s >= c1 & x$M2s < c2, "Q1", fifelse(x$M1s >= c1 & x$M2s >= c2, "Q2",
                                  fifelse(x$M1s < c1 & x$M2s < c2, "Q3", "Q4"))), levels = names(QC))
## ggplot's default breaks: extended breaks over the expanded range, minor breaks halfway
maj_br <- function(lim) { b <- scales::extended_breaks()(lim); b[b >= lim[1] - 1e-9 & b <= lim[2] + 1e-9] }
min_br <- function(lim) { b <- scales::regular_minor_breaks()(scales::extended_breaks()(lim), lim, 2); b[b >= lim[1] - 1e-9 & b <= lim[2] + 1e-9] }

cairo_pdf(file.path(out_dir, "Fig2_panels.pdf"), width = PW * S / 72, height = PH * S / 72, family = "Arial", bg = "transparent")
grid.newpage()

## ------------------------------------------------------------------ (a) qMSPI scatter, full cohort
A <- c(33.07, 80.23, 79.82, 130.27)
c1 <- median(s4$M1s); c2 <- median(s4$M2s)
s4[, Quad := quad(s4, c1, c2)]
addn("a", "medians of the M1-to-T / M2-to-T affinity (0.5 = no preference)", sprintf("%.3f / %.3f", -c1, -c2))
addn("a", "qMSPI n per quadrant, full cohort (Q1/Q2/Q3/Q4)", paste(table(s4$Quad), collapse = "/"))
xl <- range(s4$M1s) + c(-1, 1) * .05 * diff(range(s4$M1s)); yl <- range(s4$M2s) + c(-1, 1) * .05 * diff(range(s4$M2s))
ax <- function(v) lin(v, xl[1], xl[2], A[1], A[3]); ay <- function(v) lin(v, yl[1], yl[2], A[4], A[2])
for (g in sort(c(maj_br(xl), min_br(xl)))) seg(ax(g), A[2], ax(g), A[4], GRID, .25)
for (g in sort(c(maj_br(yl), min_br(yl)))) seg(A[1], ay(g), A[3], ay(g), GRID, .25)
seg(ax(c1), A[2], ax(c1), A[4], "#666666", .25, "44"); seg(A[1], ay(c2), A[3], ay(c2), "#666666", .25, "44")
ev <- s4$status_num == 1
for (q in names(QC)) {
  i <- s4$Quad == q
  ellipses(ax(s4$M1s[i & !ev]), ay(s4$M2s[i & !ev]), 1.19, 1.57, col = QC[[q]], lw = .25, alpha = .8)
  ellipses(ax(s4$M1s[i & ev]), ay(s4$M2s[i & ev]), 1.19, 1.57, fill = QC[[q]], alpha = .8)
}
box(A[1], A[2], A[3], A[4], INK, w = .25)
## ticks at the major breaks; the 0.5 tick (no preference, score -0.5) is always drawn and is the only label
for (b in unique(c(maj_br(xl), if (xl[1] < -.5 && xl[2] > -.5) -.5))) seg(ax(b), A[4], ax(b), A[4] + .67, INK, .25)
for (b in unique(c(maj_br(yl), if (yl[1] < -.5 && yl[2] > -.5) -.5))) seg(A[1] - .67, ay(b), A[1], ay(b), INK, .25)
if (xl[1] < -.5 && xl[2] > -.5) txt("0.5", ax(-.5), 135.52, 4.02, hjust = .5)
if (yl[1] < -.5 && yl[2] > -.5) txt("0.5", 29.87, ay(-.5) + 1.91, 4.02, hjust = 1)
txt("qMSPI subgroups", 34.54, 77.07, 5.03, "#3A3A3A", 2)
txt("M1-to-T affinity score", 39.19, 140.66, 4, "#3A3A3A")
txt("M2-to-T affinity score", 22.79, 127.99, 4, "#3A3A3A", rot = 90)
txt("Quadrant", 82.70, 81.95, 4.02)
ky <- 88.53 + 6.513 * 0:3
ellipses(rep(86.0, 4), ky, 2.38, 3.14, fill = QC, col = QC, lw = .5, alpha = .8)
for (k in 1:4) txt(c("Q1: M1-close/M2-far", "Q2: both close", "Q3: both far", "Q4: M2-close/M1-far")[k], 91.95, 89.90 + 6.513 * (k - 1),
                   4.02, face = if (k %in% c(1, 4)) 2 else 1)
txt("Status", 82.70, 118.48, 4.02)
ellipses(86.0, 123.50, 2.38, 3.13, col = "black", lw = .5, alpha = .8)
ellipses(86.0, 130.02, 2.38, 3.13, fill = "black", alpha = .8)
txt("Censored", 91.95, 124.87, 4.02); txt("Event", 91.95, 131.38, 4.02)

## ------------------------------------------------------------------ (b) densities and MSPI by subtype
bx <- s4[!is.na(Sub3) & !is.na(Td) & !is.na(M1d) & !is.na(M2d)]
dunn <- function(v, g) { r <- as.data.table(dunn_test(data.frame(v = v, g = g), v ~ g, p.adjust.method = "BH"))
  r[, lab := fcase(p.adj <= 1e-4, "****", p.adj <= 1e-3, "***", p.adj <= .01, "**", p.adj <= .05, "*", default = "ns")]; r }
bstats <- function(y) { q <- quantile(y, c(.25, .5, .75)); i <- q[3] - q[1]
  c(lo = min(y[y >= q[1] - 1.5 * i]), q1 = q[[1]], med = q[[2]], q3 = q[[3]], hi = max(y[y <= q[3] + 1.5 * i])) }
cat_x <- function(r, k, add = .6, n = 3) lin(k, 1 - add, n + add, r[1], r[3])
## square-root scales of the published panel (limits in sqrt(cells/mm2) and ticks in cells/mm2, from its geometry);
## kept for both versions so that the panels can be compared, and widened only if a whisker would not fit
BP <- list("T cell" = list(r = c(161.17, 36.94, 177.42, 60.83), v = "Td", fill = "#E0CE62", lim = c(-5.41, 114.28), tk = seq(0, 12500, 2500), lab = 164.69),
           M1 = list(r = c(179.68, 36.94, 195.93, 60.83), v = "M1d", fill = "#CCD5DE", lim = c(-2.02, 42.65), tk = seq(0, 1500, 500), lab = 185.70),
           M2 = list(r = c(197.73, 36.94, 213.98, 60.83), v = "M2d", fill = "#C3ABAD", lim = c(-1.63, 34.47), tk = seq(0, 900, 300), lab = 203.01))
bracket <- function(r, g1, g2, lab, yb) {
  x1 <- cat_x(r, match(g1, c("L", "H", "T"))); x2 <- cat_x(r, match(g2, c("L", "H", "T")))
  seg(x1, yb, x2, yb, "black", .25); seg(x1, yb, x1, yb + .7, "black", .25); seg(x2, yb, x2, yb + .7, "black", .25)
  if (lab == "ns") txt("ns", (x1 + x2) / 2, yb - .30, 3, hjust = .5) else txt(lab, (x1 + x2) / 2, yb + 1.3, 4, hjust = .5)
}
BY <- c("L-T" = 39.65, "H-T" = 42.66, "L-H" = 45.67)
for (nm in names(BP)) {
  P <- BP[[nm]]; r <- P$r; y <- sqrt(bx[[P$v]])
  st <- lapply(levels(bx$Sub3), function(g) bstats(y[bx$Sub3 == g]))
  lim <- c(P$lim[1], max(P$lim[2], max(sapply(st, `[[`, "hi")) * 1.05))
  by <- function(v) lin(v, lim[1], lim[2], r[4], r[2])
  for (k in 1:3) {
    s <- st[[k]]; xc <- cat_x(r, k); hw <- 2.79 / 2
    seg(xc, by(s[["hi"]]), xc, by(s[["q3"]]), INK, .25); seg(xc, by(s[["q1"]]), xc, by(s[["lo"]]), INK, .25)
    box(xc - hw, by(s[["q3"]]), xc + hw, by(s[["q1"]]), NA, P$fill, 0, alpha = .8)       # fill at 80%, as published
    box(xc - hw, by(s[["q3"]]), xc + hw, by(s[["q1"]]), INK, NA, .25)
    seg(xc - hw, by(s[["med"]]), xc + hw, by(s[["med"]]), INK, .25)
    seg(xc, r[4], xc, r[4] + .28, INK, .25)
    txt(c("L", "H", "T")[k], xc, 64.94, 4, TX, hjust = .5)
  }
  for (b in sqrt(P$tk)) if (b <= lim[2]) seg(r[1] - .28, by(b), r[1], by(b), INK, .25)
  box(r[1], r[2], r[3], r[4], INK, w = .25)
  dn <- dunn(y, bx$Sub3)
  for (i in seq_len(nrow(dn))) {
    key <- paste(dn$group1[i], dn$group2[i], sep = "-")
    bracket(r, dn$group1[i], dn$group2[i], dn$lab[i], BY[[key]])
    addn("b", sprintf("Dunn (BH) %s, %s vs %s", nm, dn$group1[i], dn$group2[i]), sprintf("%.2g %s", dn$p.adj[i], dn$lab[i]))
  }
  txt(nm, P$lab, 35.46, 4, TX)
  for (g in levels(bx$Sub3)) addn("b", sprintf("median density %s, %s (cells/mm2)", nm, g), sprintf("%.0f", median(bx[[P$v]][bx$Sub3 == g])))
}
r <- c(219.35, 36.72, 247.34, 60.91)
st <- lapply(levels(bx$Sub3), function(g) bstats(bx$MSPI[bx$Sub3 == g]))
lo_w <- min(sapply(st, `[[`, "lo")); hi_w <- max(sapply(st, `[[`, "hi"))
lim <- c(min(-.0548, lo_w - .0111), max(.1011, hi_w + .0288))
by <- function(v) lin(v, lim[1], lim[2], r[4], r[2])
mfill <- c("#E0C580", "#5B7FA0", "#5A8A6E")
for (k in 1:3) {
  s <- st[[k]]; xc <- cat_x(r, k); hw <- 4.37 / 2
  seg(xc, by(s[["hi"]]), xc, by(s[["q3"]]), INK, .25); seg(xc, by(s[["q1"]]), xc, by(s[["lo"]]), INK, .25)
  box(xc - hw, by(s[["q3"]]), xc + hw, by(s[["q1"]]), NA, mfill[k], 0, alpha = .8)
  box(xc - hw, by(s[["q3"]]), xc + hw, by(s[["q1"]]), INK, NA, .25)
  seg(xc - hw, by(s[["med"]]), xc + hw, by(s[["med"]]), INK, .25)
  seg(xc, r[4], xc, r[4] + .27, INK, .25)
  txt(c("L", "H", "T")[k], xc, 64.94, 4, TX, hjust = .5)
}
for (b in seq(ceiling(lim[1] / .05) * .05, lim[2], .05)) seg(r[1] - .27, by(b), r[1], by(b), INK, .25)
box(r[1], r[2], r[3], r[4], INK, w = .25)
dm <- dunn(bx$MSPI, bx$Sub3)
for (i in seq_len(nrow(dm))) {
  addn("b", sprintf("Dunn (BH) MSPI, %s vs %s", dm$group1[i], dm$group2[i]), sprintf("%.2g %s", dm$p.adj[i], dm$lab[i]))
  if (dm$lab[i] != "ns") bracket(r, dm$group1[i], dm$group2[i], dm$lab[i], BY[[paste(dm$group1[i], dm$group2[i], sep = "-")]])
}
txt("MSPI", 227.97, 35.46, 4, TX)
txt("Dunn's test with BH correction", 175.77, 68.99, 4)

## ------------------------------------------------------------------ (c) C-index, five models, full cohort
ci <- rt("Discovery_S4_T5_cindex")
LC <- c("Density only", "M1-to-T score", "M2-to-T score", "MSPI only", "Density + MSPI")
ci[, k := match(Model, c("Density only", "M1 affinity", "M2 affinity", "MSPI only", "Density + MSPI"))]
setorder(ci, k)
for (i in seq_len(nrow(ci))) addn("c", paste("C-index, full cohort,", LC[ci$k[i]]), sprintf("%.3f (%.3f-%.3f)", ci$C[i], ci$CI_lo[i], ci$CI_hi[i]))
r <- c(160.22, 79.43, 246.60, 125.59)
lim <- c(min(.5, min(ci$CI_lo)) - .021, max(ci$CI_hi) + .0436)
cy <- function(v) lin(v, lim[1], lim[2], r[4], r[2])
cx <- 167.57 + 17.92 * (ci$k - 1)
br <- seq(ceiling(lim[1] / .05) * .05, lim[2], .05)
for (b in br[abs(br - .5) > 1e-9]) seg(r[1], cy(b), r[3], cy(b), GRID, .5)
if (lim[1] < .5) seg(r[1], cy(.5), r[3], cy(.5), "#999999", .5, "42")
ccol <- c("#A0A0A0", "#7BA3C4", "#C47A7A", "#5B7FA0", "#5A8A6E")
for (i in seq_len(nrow(ci))) {
  col <- ccol[ci$k[i]]
  seg(cx[i], cy(ci$CI_lo[i]), cx[i], cy(ci$CI_hi[i]), col, .5)
  seg(cx[i] - 1.7, cy(ci$CI_lo[i]), cx[i] + 1.7, cy(ci$CI_lo[i]), col, .5); seg(cx[i] - 1.7, cy(ci$CI_hi[i]), cx[i] + 1.7, cy(ci$CI_hi[i]), col, .5)
  if (ci$k[i] %in% c(1, 5)) diamonds(cx[i], cy(ci$C[i]), 2.94, col) else ellipses(cx[i], cy(ci$C[i]), 2.94, 2.94, fill = col)
}
box(r[1], r[2], r[3], r[4], INK, w = .25)
for (b in br) { seg(r[1] - .88, cy(b), r[1], cy(b), INK, .5); txt(sprintf("%.2f", b), 158.29, cy(b) + .93, 4, TX, hjust = 1) }
for (x in 167.57 + 17.92 * 0:4) seg(x, r[4], x, r[4] + .88, INK, .5)
lo_xy <- list(c(149.35, 139.32), c(166.06, 140.44), c(183.98, 140.44), c(204.88, 137.66), c(218.34, 141.81))
for (k in 1:5) txt(LC[k], lo_xy[[k]][1], lo_xy[[k]][2], 4, if (k %in% 2:3) "#3A3A3A" else TX, rot = 24.8)
txt("C-index", 143.49, 108.74, 3.56, rot = 90); txt("with bootstrap 95% CI", 147.76, 120.15, 3.56, rot = 90)

## ------------------------------------------------------------------ (d) MSPI HR per SD by subtype
fd <- fread(file.path(CLI, "discovery_MSPI_subtype_forest_data_merged.csv"))
GL <- c("Full cohort", "Luminal", "HER2+", "TNBC")
fd[, k := match(Subgroup, c("All-comers", "Luminal", "HER2+", "TNBC"))]
setorder(fd, k)
for (i in seq_len(nrow(fd))) {
  addn("d", paste0(GL[fd$k[i]], ": n, events"), sprintf("%d, %d", fd$n[i], fd$events[i]))
  addn("d", paste0(GL[fd$k[i]], ": MSPI HR per SD, UV"), ci_txt(fd$HR[i], fd$lower[i], fd$upper[i], fd$p[i]))
  addn("d", paste0(GL[fd$k[i]], ": MSPI HR per SD, density-adjusted"), ci_txt(fd$HR_adj[i], fd$lower_adj[i], fd$upper_adj[i], fd$p_adj[i]))
}
r <- c(290.98, 38.29, 382.82, 128.29)
lo0 <- min(c(fd$lower, fd$lower_adj)) - .05; hi0 <- max(c(fd$upper, fd$upper_adj)) + .6
lim <- c(min(.4716, lo0 - .05 * (hi0 - lo0)), max(1.9229, max(c(fd$upper, fd$upper_adj)) + .178))   # published edges at least
dx <- function(v) lin(v, lim[1], lim[2], r[1], r[3])
rows <- list(c(45.88, 54.19), c(68.05, 76.36), c(90.21, 98.53), c(112.38, 120.69))
dcol <- c("#A67C52", "#D4B86A", "#8C8C8C", "#8C8C8C")
brd <- seq(ceiling(lim[1] / .5) * .5, lim[2], .5)
for (b in brd[abs(brd - 1) > 1e-9]) seg(dx(b), r[2], dx(b), r[4], GRID, .5)
seg(dx(1), r[2], dx(1), r[4], "#999999", .5, "22")
for (i in seq_len(nrow(fd))) {
  yy <- rows[[fd$k[i]]]
  for (j in 1:2) {
    h <- if (j == 1) fd$HR[i] else fd$HR_adj[i]; l <- if (j == 1) fd$lower[i] else fd$lower_adj[i]; u <- if (j == 1) fd$upper[i] else fd$upper_adj[i]
    seg(dx(l), yy[j], dx(u), yy[j], BAR, .5)
    seg(dx(l), yy[j] - 1.385, dx(l), yy[j] + 1.385, BAR, .5); seg(dx(u), yy[j] - 1.385, dx(u), yy[j] + 1.385, BAR, .5)
    if (j == 1) ellipses(dx(h), yy[j], 4.73, 4.73, fill = dcol[fd$k[i]]) else diamonds(dx(h), yy[j], 4.73, dcol[fd$k[i]])
  }
}
box(349.12, 114.51, 381.81, 127.34, NA, "white", 0)
box(349.12, 114.50, 381.81, 127.34, "#999999", NA, .5, "22")
ellipses(352.77, 118.56, 3.8, 3.8, fill = "#8C8C8C"); diamonds(352.77, 124.07, 3.8, "#8C8C8C")
txt("MSPI", 355.98, 119.70, 4, TX); txt("MSPI+density", 355.98, 124.90, 4, TX)
box(r[1], r[2], r[3], r[4], INK, w = .25)
for (yy in unlist(rows)) seg(r[1] - 1.43, yy, r[1], yy, INK, .5)
for (b in brd) { seg(dx(b), r[4], dx(b), r[4] + 1.42, INK, .5); txt(sprintf("%.1f", b), dx(b), 134.13, 4, hjust = .5) }
txt("Full cohort", 271.23, 50.92, 4, TX); txt("Luminal", 273.55, 73.39, 4, TX)
txt("HER2+", 274.25, 95.87, 4, TX); txt("TNBC", 276.20, 118.34, 4, TX)
txt("Hazard Ratio per SD (95% CI)", 304.07, 140.88, 5, TX)

## ------------------------------------------------------------------ (e, h) clinical multivariable Cox
nev <- function(pdf) { t <- paste(system2("pdftotext", c(shQuote(pdf), "-"), stdout = TRUE), collapse = " ")
  m <- regmatches(t, regexec("n=([0-9]+), events=([0-9]+)", t))[[1]]; as.integer(m[2:3]) }
relab <- function(x) { x <- sub("^Age: ≥40$", "Age ≥ 40", x); x <- sub("^Stage: ", "", x)
  x <- sub("^Systemic treatment: Chemotherapy and Endocrine treatment$", "Chemo- and Endocrine treatment", x)
  sub("^Systemic treatment: ", "", x) }
forest_clin <- function(stem, panel, title, G) {
  ft <- rt(stem); ne <- nev(file.path(PLT, paste0(stem, ".pdf")))
  ft[, lab := relab(label)]
  for (i in seq_len(nrow(ft))) addn(panel, paste0(title, ": ", ft$lab[i]), ci_txt(ft$HR[i], ft$CI_lo[i], ft$CI_hi[i], ft$p[i]))
  addn(panel, paste0(title, ": n, events"), sprintf("%d, %d", ne[1], ne[2]))
  fx <- function(v) G$x1 + G$sl * log10(v)
  for (b in c(.5, 2, 5)) seg(fx(b), G$top, fx(b), G$bot, GRID, .5)
  seg(fx(1), G$top, fx(1), G$bot, "#A6A6A6", .5, "22")
  for (i in seq_len(nrow(ft))) {
    yy <- G$y1 + G$dy * (i - 1)
    seg(fx(ft$CI_lo[i]), yy, fx(ft$CI_hi[i]), yy, BAR, .5)
    seg(fx(ft$CI_lo[i]), yy - G$cap / 2, fx(ft$CI_lo[i]), yy + G$cap / 2, BAR, .5)
    seg(fx(ft$CI_hi[i]), yy - G$cap / 2, fx(ft$CI_hi[i]), yy + G$cap / 2, BAR, .5)
    diamonds(fx(ft$HR[i]), yy, G$dm, if (ft$Variable[i] == "MSPI_sd") "#5B7FA0" else "#D4B86A")
    txt(ft$lab[i], G$labx, yy + G$off, 4, TX, hjust = 1)
  }
  for (b in c(.5, 1, 2, 5)) txt(sprintf("%.1f", b), fx(b), G$xlab_y, 4, hjust = .5)
  txt("Hazard Ratio (95% CI)", G$at[1], G$at[2], 4, TX)
  txt(title, G$tt[1], G$tt[2], 4, TX, 3)
  txt(sprintf("n=%d, events=%d", ne[1], ne[2]), G$st[1], G$st[2], 4, TX, 3, hjust = .5)
}
forest_clin("Discovery_Merged_Luminal_T6_forest", "e", "Luminal (DDFS)",
            list(x1 = 477.71, sl = 29.09, top = 38.66, bot = 129.49, y1 = 44.00, dy = 8.906, cap = 1.96, dm = 2.97, labx = 454.05,
                 off = 1.13, xlab_y = 134.88, at = c(465.55, 140.98), tt = c(404.92, 32.59), st = c(472.47, 36.49)))
forest_clin("Discovery_Young_T6_forest", "h", "Young Early Luminal (DDFS)",
            list(x1 = 476.75, sl = 28.77, top = 184.45, bot = 265.30, y1 = 190.37, dy = 9.86, cap = 2.17, dm = 3.26, labx = 444.52,
                 off = 1.03, xlab_y = 270.45, at = c(463.76, 276.55), tt = c(392.26, 171.86), st = c(463.29, 182.07)))

## ------------------------------------------------------------------ (f, g) early and young early luminal
panel_fg <- function(x, stem, panel, title, G) {
  x <- copy(x)
  qb <- unique(quantile(x$MSPI, probs = seq(0, 1, .25), type = 2)); stopifnot(length(qb) == 5)
  x[, Quad := cut(MSPI, breaks = qb, include.lowest = TRUE, labels = names(MQC))]   # MSPI quartile, as in script 05
  fit <- survfit(Surv(DDFStime2021, status_num) ~ Quad, data = x)
  p41 <- 1 - pchisq(survdiff(Surv(DDFStime2021, status_num) ~ Quad, data = x[Quad %in% c("Q1", "Q4")])$chisq, 1)
  p4g <- 1 - pchisq(survdiff(Surv(DDFStime2021, status_num) ~ Quad, data = x)$chisq, 3)
  ptr <- summary(coxph(Surv(DDFStime2021, status_num) ~ as.integer(Quad), data = x))$coefficients[1, 5]
  addn(panel, "n, events", sprintf("%d, %d", nrow(x), sum(x$status_num)))
  addn(panel, "MSPI quartiles, n (1st/2nd/3rd/4th)", paste(table(x$Quad), collapse = "/"))
  addn(panel, "MSPI quartiles, events (1st/2nd/3rd/4th)", paste(x[, sum(status_num), keyby = Quad]$V1, collapse = "/"))
  addn(panel, "MSPI quartiles, log-rank 4th vs 1st (shown)", p_txt(p41))
  addn(panel, "MSPI quartiles, log-rank, four groups (not shown)", p_txt(p4g))
  addn(panel, "MSPI quartiles, trend (Cox, quartile 1-4, not shown)", p_txt(ptr))
  ## Kaplan-Meier
  K <- G$km; kx <- function(t) lin(t, -1.5, 31.5, K[1], K[3]); ky <- function(s) lin(s, -.05, 1.05, K[4], K[2])
  st <- data.table(time = fit$time, surv = fit$surv, n.censor = fit$n.censor, Quad = rep(sub("Quad=", "", names(fit$strata)), fit$strata))
  for (q in names(MQC)) {
    s <- st[Quad == q]
    tt <- c(0, rep(s$time, each = 2))                             # step path: (0,1), (t1,1), (t1,s1), (t2,s1), ...
    ss <- c(1, as.vector(rbind(c(1, head(s$surv, -1)), s$surv)))
    grid.lines(X(kx(tt)), Y(ky(ss)), gp = gpar(col = MQC[[q]], lwd = LW(1), lineend = "butt", linejoin = "mitre"))
  }
  for (q in names(MQC)) { s <- st[Quad == q & n.censor > 0]
    if (nrow(s)) grid.text("+", X(kx(s$time)), Y(ky(s$surv) - .44), hjust = .5, vjust = .5,
                           gp = gpar(fontsize = G$plus * S, fontfamily = "Helvetica", col = MQC[[q]])) }
  seg(K[1], K[2], K[1], K[4], "black", .5); seg(K[1], K[4], K[3], K[4], "black", .5)
  for (b in seq(0, 1, .25)) { seg(K[1] - 1.32, ky(b), K[1], ky(b), "black", .5); txt(sprintf("%.2f", b), G$ylab_r, ky(b) + 1.39, 4, hjust = 1) }
  for (b in c(0, 10, 20, 30)) { seg(kx(b), K[4], kx(b), K[4] + 1.32, "black", .5); txt(b, kx(b), G$xlab_y, 4, hjust = .5) }
  ## legend in the empty lower left corner of the plot, title above the rows: "MSPI quartile" is longer than the
  ## submitted "qMSPI" and would cross the young early luminal curves at the submitted place
  stopifnot(min(st$surv) > .45)                                    # the corner stays free of curves
  Ty <- K[4] - 31.22; ly <- Ty + 5.2 + 5.94 * 0:3
  txt("MSPI", K[1] + 3.2, Ty, 5, face = 2)
  wm <- convertWidth(grobWidth(textGrob("MSPI", gp = gpar(fontsize = 5 * S, fontface = 2, fontfamily = "Arial"))), "bigpts", TRUE) / S
  txt("quartile", K[1] + 3.2 + wm + 1.4, Ty, 5, face = 3)
  ellipses(rep(K[1] + 5.16, 4), ly, 4.0, 5.33, fill = MQC, col = MQC, lw = .25, alpha = .8)
  for (k in 1:4) txt(c("1st", "2nd", "3rd", "4th")[k], K[1] + 9.63, ly[k] + 1.03, 4)
  bx <- K[1] + 21.2
  seg(bx, ly[1], bx, ly[4], BAR, .5); seg(bx - 3.79, ly[1], bx, ly[1], BAR, .5); seg(bx - 3.79, ly[4], bx, ly[4], BAR, .5)
  txt(paste("p =", p_txt(p41)), bx + 2.6, mean(ly) + 1.4, 4, face = 3)
  txt(title, G$tt[1], G$tt[2], 5, TX, 3)
  txt("Survival probability", G$sp[1], G$sp[2], 4, rot = 90)
  txt("DDFS (years)", G$dd[1], G$dd[2], 5)
  ## mini forest (log scale)
  mf <- rt(paste0(stem, "_mini_forest"))
  for (i in seq_len(nrow(mf))) addn(panel, paste("MSPI HR per SD,", mf$Model[i]), ci_txt(mf$HR[i], mf$CI_lo[i], mf$CI_hi[i], mf$p[i]))
  M <- G$mf; U <- log(max(mf$CI_hi)); R <- U / (1 - 2 * .045); ml <- c(-.045 * R, U + .045 * R)
  mx <- function(v) lin(log(v), ml[1], ml[2], M[1], M[3])
  maj <- labeling::extended(exp(ml[1]), exp(ml[2]), m = 3); maj <- maj[maj >= exp(ml[1]) & maj <= exp(ml[2])]
  stp <- diff(maj)[1]; mnr <- if (stp < 1) head(maj, -1) + stp / 2 else numeric(0)
  box(M[1], M[2], M[3], M[4], NA, "white", 0)
  for (b in maj[maj > 1]) seg(mx(b), M[2], mx(b), M[4], GRID, .5)
  seg(mx(1), M[2], mx(1), M[4], "#C47A7A", .5, "22")
  ry <- c(249.41, 254.96)
  for (j in 1:2) {
    seg(mx(mf$CI_lo[j]), ry[j], mx(mf$CI_hi[j]), ry[j], BAR, .5)
    seg(mx(mf$CI_lo[j]), ry[j] - .61, mx(mf$CI_lo[j]), ry[j] + .61, BAR, .5); seg(mx(mf$CI_hi[j]), ry[j] - .61, mx(mf$CI_hi[j]), ry[j] + .61, BAR, .5)
    if (j == 1) ellipses(mx(mf$HR[j]), ry[j], 3, 3, fill = "#C47A7A") else diamonds(mx(mf$HR[j]), ry[j], 3, "#C47A7A")
  }
  box(M[1], M[2], M[3], M[4], INK, NA, .25)
  for (yy in ry) seg(M[1] - .47, yy, M[1], yy, INK, .5)
  for (b in c(maj, mnr)) seg(mx(b), M[4], mx(b), M[4] + .47, INK, .5)
  for (b in maj) txt(sprintf("%.1f", b), mx(b), 262.48, 4, hjust = .5)
  txt("UV", G$uv[1], G$uv[2], 4, TX, hjust = 1); txt("MV", G$uv[1], G$uv[3], 4, TX, hjust = 1)
  txt("MSPI", G$mt[1], G$mt[2], 5, face = 2); txt("HR per SD", G$mt[1] + 13.61, G$mt[2], 5, face = 3)
  ## pie: Type II Anova of the density-adjusted Cox model
  xx <- copy(x)[!is.na(Td) & !is.na(M1d) & !is.na(M2d)]; xx[, MSPI_sd := MSPI / sd(MSPI)]
  a2 <- car::Anova(coxph(Surv(DDFStime2021, status_num) ~ MSPI_sd + M1d + M2d + Td, data = xx), type = "II")
  pp <- data.table(v = c("MSPI", "d M1", "d M2", "d T cell"), chisq = a2[, 1], p = a2[, 3])
  pp[, pct := 100 * chisq / sum(chisq)]
  for (i in seq_len(nrow(pp))) addn(panel, paste0("pie, ", pp$v[i], " (share of chi-square, p)"), sprintf("%.1f%%, p=%.2g", pp$pct[i], pp$p[i]))
  pc <- c(MSPI = "#8EABBF", "d M1" = "#CCD5DE", "d M2" = "#C3ABAD", "d T cell" = "#E0CE62")
  pp <- pp[match(c("d T cell", "d M2", "d M1", "MSPI"), v)]            # clockwise from 12 o'clock, as published
  pp[, end := cumsum(pct) / 100 * 2 * pi]; pp[, start := end - pct / 100 * 2 * pi]
  ctr <- G$pie; R0 <- 12.5
  for (i in seq_len(nrow(pp))) if (pp$pct[i] > 0) {
    th <- seq(pp$start[i], pp$end[i], length.out = max(3, ceiling(pp$pct[i] * 2)))
    grid.polygon(X(c(ctr[1], ctr[1] + R0 * sin(th))), Y(c(ctr[2], ctr[2] - R0 * cos(th))), gp = gpar(fill = pc[[pp$v[i]]], col = NA))
  }
  pp[, plab := fifelse(p < .001, sprintf("p=%.1e", p), sprintf("p=%.3f", p))]
  pp[, mid := (start + end) / 2]; pp[, w := nchar(plab) * 2.03]
  ## label places as in the published pies: the MSPI p-value low in the centre when the MSPI slice covers the
  ## bottom, otherwise along its middle; the other slices left-aligned from 0.4 R towards their middle
  pp[, `:=`(hj = 0, cxl = ctr[1] + .4 * R0 * sin(mid), cyl = ctr[2] - .75 * R0 * cos(mid))]
  im <- which(pp$v == "MSPI")
  if (pp$start[im] < pi - .5 && pp$end[im] > pi + .5) pp[im, `:=`(hj = .5, cxl = ctr[1] - .7, cyl = ctr[2] + 3.7)] else
    pp[im, `:=`(hj = .5, cxl = ctr[1] + .4 * R0 * sin(mid), cyl = ctr[2] - .4 * R0 * cos(mid))]
  pp[sin(mid) < 0 & v != "MSPI", hj := 1]
  ins <- which(pp$pct >= 4)
  for (i in ins) for (j in ins) if (i < j && abs(pp$cyl[i] - pp$cyl[j]) < 4.3 &&
                                    abs((pp$cxl[i] + (.5 - pp$hj[i]) * pp$w[i]) - (pp$cxl[j] + (.5 - pp$hj[j]) * pp$w[j])) < (pp$w[i] + pp$w[j]) / 2) {
    lo <- if (pp$cyl[i] > pp$cyl[j]) i else j; hi <- setdiff(c(i, j), lo); pp$cyl[lo] <- pp$cyl[hi] + 4.3 }
  for (i in seq_len(nrow(pp))) {
    if (pp$pct[i] >= 4) txt(pp$plab[i], pp$cxl[i], pp$cyl[i] + 1.4, 4, TX, 3, hjust = pp$hj[i]) else {
      e <- c(ctr[1] + (R0 + .4) * sin(pp$mid[i]), ctr[2] - (R0 + .4) * cos(pp$mid[i]))
      lx <- e[1] + 4.0; ly <- e[2] - .5
      seg(e[1] + .2, e[2] - .1, lx - .45, ly - .9, pc[[pp$v[i]]], .5)
      txt(pp$plab[i], lx, ly, 4, TX, 3)
    }
  }
  for (k in 1:4) { box(G$sq[1], G$sq[2] + 4.45 * (k - 1), G$sq[1] + 4.09, G$sq[2] + 4.45 * (k - 1) + 4.09, NA, pc[[k]], 0)
    txt(names(pc)[k], G$sql[1], G$sql[2] + 4.495 * (k - 1), 4, TX) }
  ## C-index, MSPI only and density + MSPI
  cf <- rt(paste0(stem, "_T5_cindex_MSPIfocused"))
  for (i in seq_len(nrow(cf))) addn(panel, paste("C-index,", cf$Model[i]), sprintf("%.3f (%.3f-%.3f)", cf$C[i], cf$CI_lo[i], cf$CI_hi[i]))
  Cb <- G$cb; R1 <- (max(cf$CI_hi) + .08) - (min(cf$CI_lo) - .03)
  cl <- c(min(cf$CI_lo) - .03 - .05 * R1, max(cf$CI_hi) + .043)
  cyy <- function(v) lin(v, cl[1], cl[2], Cb[4], Cb[2]); cxx <- lin(c(1, 2), .4, 2.6, Cb[1], Cb[3])
  brc <- seq(ceiling(cl[1] / .1) * .1, cl[2], .1)
  for (b in brc) seg(Cb[1], cyy(b), Cb[3], cyy(b), GRID, .5)
  cc <- c("#5B7FA0", "#5A8A6E")
  for (j in 1:2) {
    i <- match(c("MSPI only", "Density + MSPI")[j], cf$Model)
    seg(cxx[j], cyy(cf$CI_lo[i]), cxx[j], cyy(cf$CI_hi[i]), cc[j], .5)
    seg(cxx[j] - .835, cyy(cf$CI_lo[i]), cxx[j] + .835, cyy(cf$CI_lo[i]), cc[j], .5)
    seg(cxx[j] - .835, cyy(cf$CI_hi[i]), cxx[j] + .835, cyy(cf$CI_hi[i]), cc[j], .5)
    if (j == 1) ellipses(cxx[j], cyy(cf$C[i]), 3.63, 3.63, fill = cc[j]) else diamonds(cxx[j], cyy(cf$C[i]), 3.63, cc[j])
  }
  box(Cb[1], Cb[2], Cb[3], Cb[4], INK, NA, .25, "44")
  for (b in brc) { seg(Cb[3], cyy(b), Cb[3] + 1.0, cyy(b), INK, .5); txt(sprintf("%.1f", b), Cb[3] + 1.66, cyy(b) + 1.04, 4, TX) }
  for (xx2 in cxx) seg(xx2, Cb[4], xx2, Cb[4] + 1.0, INK, .5)
  txt("C-index with bootstrap 95% CI", G$ct[1], G$ct[2], 4.02, rot = 90)
  txt("MSPI only", G$xl1[1], G$xl1[2], 4, TX, rot = 15); txt("Density + MSPI", G$xl2[1], G$xl2[2], 4, TX, rot = 15)
}
panel_fg(early, "Discovery_Early", "f", "Early Luminal (DDFS)",
         list(km = c(38.95, 180.43, 136.61, 264.22), plus = 5.69, ylab_r = 35.83, xlab_y = 270.91,
              tt = c(60.46, 177.50), sp = c(25.78, 240.70), dd = c(69.20, 278.58), mf = c(93.05, 246.08, 137.50, 258.29),
              uv = c(89.59, 251.12, 256.67), mt = c(92.58, 244.08), pie = c(147.50, 191.33), sq = c(171.35, 182.16),
              sql = c(177.16, 184.99), cb = c(149.19, 205.45, 175.42, 264.72), ct = c(186.98, 260.30),
              xl1 = c(140.39, 274.47), xl2 = c(147.81, 276.47)))
panel_fg(young, "Discovery_Young", "g", "Young Early Luminal (DDFS)",
         list(km = c(222.21, 180.01, 320.36, 264.22), plus = 5.72, ylab_r = 218.83, xlab_y = 270.44,
              tt = c(237.78, 177.94), sp = c(208.18, 242.81), dd = c(252.59, 278.15), mf = c(276.68, 246.08, 321.13, 258.29),
              uv = c(274.66, 251.45, 256.99), mt = c(277.22, 243.93), pie = c(330.59, 190.02), sq = c(351.92, 181.31),
              sql = c(357.73, 184.14), cb = c(332.01, 205.03, 358.24, 264.17), ct = c(369.48, 259.78),
              xl1 = c(323.20, 273.91), xl2 = c(330.63, 276.92)))
invisible(dev.off())
fwrite(rbindlist(nums), file.path(out_dir, "Fig2_numbers.tsv"), sep = "\t")
cat("written:", file.path(out_dir, "Fig2_panels.pdf"), "and Fig2_numbers.tsv\n")
