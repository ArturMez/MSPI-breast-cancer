########################
## 16 - Supplementary Figure S3: nucleus-supported identification of CD4+ T cells in the discovery cohort
##
## The nucleus classifier of script 00b decides which CD4+ cells count as T cells (CD4 is also expressed by monocytes,
## macrophages and dendritic cells). Cells of the patients of the MSPI analysis (both affinity scores, DDFS time > 0,
## classified subtype, as script 05). Cell groups as in python/01: macrophage marker = CD68 > 5 or CD163+,
## lymphocyte-like nucleus = predicted probability above 0.5.
##   a  schematic of a CD8+ T-cell nucleus and a macrophage nucleus with the classifier measures, and DAPI intensity
##      (scaled within each core) and nucleus axis ratio in the stroma by cell group (median, IQR, 10th-90th percentile)
##   b  schematic of the classifier: training on CD8+ T cells against CD68+ or CD163+ cells without CD4 or CD8, test
##      in the held-out cores, application to CD4+ cells without CD8
##   c  proportion of cells with a lymphocyte-like nucleus by cell group, in all tissue compartments
##   d  CD4+ cells without CD8 counted as T cells (lymphocyte-like nucleus) or excluded (non-lymphocyte-like)
##
## Inputs : ./data/discovery/ (SAV, affinity scores), ./output/discovery_nucleus/ (scripts 00a and 00b:
##          cells_with_nuclei.tsv, cells_for_affinity_nucleus.tsv, nucleus_classifier_auc.tsv)
## Outputs: ./output/supplementary_figures/Supplementary_Figure_S3.pdf/.png (152 x 200 mm, the text width of the
##          supplement, so the font sizes set here are the printed sizes), S3_numbers.tsv, S3_totals.tsv
## Usage  : Rscript R/16_Supplementary_Figure_S3_nucleus_classifier.R
########################

suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(foreign); library(ggtext)})
disc <- "./data/discovery"
nuc <- "./output/discovery_nucleus"
OUT <- "./output/supplementary_figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- palette and theme (the colours of Figures 1 and 2)
col_T <- "#D4B86A"; col_mac <- "#A0A0A0"; col_excl <- "#8C8C8C"; col_tag <- "#CCAE8F"
ink <- "#333333"; ink2 <- "#4D4D4D"; ref <- "grey55"
nucT <- "#4F72B3"; nucT2 <- "#3F5F9E"; nucM <- "#B3C3DD"; nucM2 <- "#7F98C4"; nucC <- "#819ACB"; nucC2 <- "#5F7DB8"
pt <- function(x) x / ggplot2::.pt          # printed size in points to the mm size of annotate() and geom_text()
th <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2, size = 7.5),
        axis.line = element_line(linewidth = 0.3, colour = ink), axis.ticks = element_line(linewidth = 0.3, colour = ink),
        plot.title = element_text(size = 8.5, colour = ink, margin = margin(b = 5)),
        plot.tag = element_text(face = "bold", size = 12, colour = "black"),
        legend.key.size = unit(3.2, "mm"), legend.text = element_text(size = 7.5, colour = ink2),
        strip.background = element_blank(), strip.text = element_text(face = "italic", size = 8, colour = ink2))
grid_x <- theme(panel.grid.major.x = element_line(colour = "#EBEBEB", linewidth = 0.3))
blank_axes <- theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                    axis.title = element_blank())
glev <- c("CD8+", "CD8+ marker+", "CD4+", "CD4+ marker+", "Macrophages")   # internal keys (numbers files)
# labels the reader sees: each group named by its markers, in two lines (markdown, drawn by ggtext)
glab <- c("CD8+" = "CD8<sup>+</sup><br>no CD68 or CD163", "CD8+ marker+" = "CD8<sup>+</sup><br>with CD68 or CD163",
          "CD4+" = "CD4<sup>+</sup><br>no CD68 or CD163", "CD4+ marker+" = "CD4<sup>+</sup><br>with CD68 or CD163",
          "Macrophages" = "Macrophages<br>CD68<sup>+</sup>, no CD4 or CD8")
y_groups <- scale_y_discrete(labels = function(x) unname(glab[x]))
# groups drawn in the boxplots and in the lymphocyte-like panel. The numbers files keep all five groups.
gshow <- c("CD8+", "CD4+", "CD4+ marker+", "Macrophages")
y_md <- theme(axis.text.y = element_markdown(hjust = 1, halign = 1, lineheight = 1.15, size = 7.5, colour = ink2))
grp5 <- function(dt) dt[, g := fcase(
  CD8 == 1 & CD4 == 0 & !(CD68 > 5 | CD163 == 1), "CD8+", CD8 == 1 & CD4 == 0, "CD8+ marker+",
  CD4 == 1 & CD8 == 0 & !(CD68 > 5 | CD163 == 1), "CD4+", CD4 == 1 & CD8 == 0, "CD4+ marker+",
  CD4 == 0 & CD8 == 0 & CD68 > 5, "Macrophages", default = NA_character_)]

## ---- analysis set: both affinity scores, DDFS time > 0, classified subtype (682 patients)
cl <- suppressWarnings(read.spss(list.files(disc, pattern = "workingr_832_updated_211024\\.sav$", full.names = TRUE)[1],
                                 use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
                                 trim.factor.names = TRUE, trim_values = TRUE, reencode = NA))
if (!"Subtype_BCY21" %in% names(cl)) names(cl)[names(cl) == "Subtype_Discovery21"] <- "Subtype_BCY21"
cl <- as.data.table(cl[as.character(cl$Subtype_BCY21) != "Unclassified", c("AnnUnit", "DDFStime2021")])
cl[, `:=`(AnnUnit = as.character(AnnUnit), DDFStime2021 = as.numeric(DDFStime2021))]
rdsc <- function(f) fread(file.path(disc, f))[, .(AnnUnit = as.character(AnnUnit), d = delta_CD3)]
sc <- merge(setnames(rdsc("discovery_distances_mean_M1_to_CD3_nucleus.txt"), "d", "M1d"),
            setnames(rdsc("discovery_distances_mean_M2_to_CD3_nucleus.txt"), "d", "M2d"), by = "AnnUnit")
aset <- merge(cl, sc, by = "AnnUnit")[!is.na(M1d) & !is.na(M2d) & !is.na(DDFStime2021) & DDFStime2021 > 0, AnnUnit]
cat("patients:", length(aset), "\n")

## ---- A: nucleus features in the stroma
nb <- fread(file.path(nuc, "cells_with_nuclei.tsv"),
            select = c("Core.ID", "Sample.Name", "Tissue.Category", "NucArea", "NucAxisRatio", "NucDAPI", "CD4", "CD8",
                       "CD68", "CD163"), colClasses = list(character = "Sample.Name"))
nb <- nb[!is.na(NucArea) & NucArea > 0]
nb[, dapi_z := (NucDAPI - median(NucDAPI)) / mad(NucDAPI), by = Core.ID]   # all cells of the core, as in 00b
nb <- nb[Sample.Name %in% aset & Tissue.Category == "Stroma"]
grp5(nb)
nb <- nb[!is.na(g)]
qq <- function(x, p) as.numeric(quantile(x, p))
bstat <- rbindlist(lapply(c("dapi_z", "NucAxisRatio"), function(v)
  nb[, .(feature = v, n = .N, lo = qq(get(v), .1), q1 = qq(get(v), .25), md = qq(get(v), .5), q3 = qq(get(v), .75),
         hi = qq(get(v), .9)), by = g]))
rm(nb); invisible(gc())

## ---- C and D: classifier calls by group and CD4+ cells kept or excluded
pl <- fread(file.path(nuc, "cells_for_affinity_nucleus.tsv"),
            select = c("Sample.Name", "CD4", "CD8", "CD68", "CD163", "p_lymph"), colClasses = list(character = "Sample.Name"))
pl <- pl[Sample.Name %in% aset]
pl[, lymph := !is.na(p_lymph) & p_lymph > 0.5]
grp5(pl)
lstat <- pl[!is.na(g), .(n = .N, lymph = sum(lymph)), by = g][, pct := 100 * lymph / n]
tot <- pl[, .(patients = uniqueN(Sample.Name), T_original = sum(CD4 == 1 | CD8 == 1),
              T_nucleus = sum(CD8 == 1 | (CD4 == 1 & lymph)), CD4_only = sum(CD4 == 1 & CD8 == 0),
              CD4_only_excluded = sum(CD4 == 1 & CD8 == 0 & !lymph))]
rm(pl); invisible(gc())
auc <- fread(file.path(nuc, "nucleus_classifier_auc.tsv"))$AUC_heldout_cores
stopifnot(length(auc) == 1)
fwrite(merge(bstat, lstat[, .(g, n_all = n, lymph, pct_lymph = pct)], by = "g", all.x = TRUE),
       file.path(OUT, "S3_numbers.tsv"), sep = "\t")
fwrite(tot, file.path(OUT, "S3_totals.tsv"), sep = "\t")

## ---- drawing helpers
ell <- function(cx, cy, area, ratio, n = 240) {
  a <- sqrt(area * ratio / pi); b <- sqrt(area / (pi * ratio)); t <- seq(0, 2 * pi, length.out = n)
  data.table(x = cx + a * cos(t), y = cy + b * sin(t))
}
circ <- function(cx, cy, r, n = 240) { t <- seq(0, 2 * pi, length.out = n); data.table(x = cx + r * cos(t), y = cy + r * sin(t)) }
spots <- function(cx, cy, dx, dy, r) rbindlist(lapply(seq_along(dx), function(i) cbind(circ(cx + dx[i], cy + dy[i], r[i], 60), id = i)))
# a small nucleus icon: dashed cell outline, nucleus, chromatin spots; s scales the icon
icon <- function(p, cx, cy, kind, s = 1) {
  pars <- list(T = list(ratio = 1.12, f = nucT, f2 = nucT2, cell = 4.0),
               M = list(ratio = 1.55, f = nucM, f2 = nucM2, cell = 4.6),
               C = list(ratio = 1.32, f = nucC, f2 = nucC2, cell = 4.3))[[kind]]
  sp <- spots(cx, cy, s * c(-1.1, 0.9, 0.2, -0.3), s * c(0.8, 0.6, -1.0, -0.1), s * c(0.42, 0.38, 0.4, 0.34))
  p + geom_polygon(data = circ(cx, cy, pars$cell * s), aes(x, y), fill = NA, colour = "grey55", linewidth = 0.25,
                   linetype = "22", inherit.aes = FALSE) +
    geom_polygon(data = ell(cx, cy, 22 * s^2, pars$ratio), aes(x, y), fill = pars$f, colour = ink, linewidth = 0.25,
                 inherit.aes = FALSE) +
    geom_polygon(data = sp, aes(x, y, group = id), fill = pars$f2, colour = NA, inherit.aes = FALSE)
}
txt <- function(x, y, label, size = pt(7), hjust = 0.5, parse = FALSE, colour = ink2, fontface = "plain")
  annotate("text", x = x, y = y, label = label, size = size, hjust = hjust, parse = parse, colour = colour,
           family = "Arial", fontface = fontface)
# text with markdown (italic group names in a darker ink), as the classifier measures of the nucleus schematic
rtxt <- function(x, y, label, size = pt(7), hjust = 0.5)
  annotate("richtext", x = x, y = y, label = label, size = size, hjust = hjust, vjust = 0.5, colour = ink2,
           family = "Arial", fill = NA, label.colour = NA, label.padding = unit(0, "pt"))
box <- function(x0, x1, y0, y1, colour = ink, fill = "white", lw = 0.35)
  annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, colour = colour, fill = fill, linewidth = lw)
arr <- function(x0, y0, x1, y1)
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1, linewidth = 0.35, colour = ink,
           arrow = arrow(length = unit(1.1, "mm"), type = "closed"))

## ---- B: training and application of the classifier (schematic). Units are 1 mm at the printed width (152 mm): the
## boxes are sized to their text, so the panel must be drawn at that scale (row height of at least 50 mm).
panel_training <- function(tag = "b") {
  s7 <- pt(7); s75 <- pt(7.5); s85 <- pt(8.5)
  p <- ggplot() + coord_fixed(xlim = c(-1, 147), ylim = c(0, 46), expand = FALSE, clip = "off")
  # training: two classes, model, test
  p <- p + txt(0, 43.6, "Training", size = s85, hjust = 0, colour = ink, fontface = "bold") +
    box(0, 40, 25.5, 39.5, colour = col_T, lw = 0.5) + box(0, 40, 9.5, 23.5, colour = col_mac, lw = 0.5)
  p <- icon(p, 5.4, 32.5, "T", s = 0.8)
  p <- icon(p, 5.4, 16.5, "M", s = 0.8)
  p <- p + txt(10.6, 34.6, "bold('CD8'^'+'~'T cells')", size = s75, hjust = 0, parse = TRUE, colour = ink) +
    txt(10.6, 30.4, "no CD4, CD68 or CD163", size = s7, hjust = 0) +
    txt(10.6, 18.6, "bold('CD68'^'+'~'or CD163'^'+'~'cells')", size = s75, hjust = 0, parse = TRUE, colour = ink) +
    txt(10.6, 14.4, "no CD4 or CD8", size = s7, hjust = 0) +
    txt(20, 5.9, "150,000 cells per class,", size = s7) + txt(20, 2.6, "70% of the TMA cores", size = s7) +
    arr(40, 32.5, 45, 28.5) + arr(40, 16.5, 45, 20.5) +
    box(45.5, 85, 9.5, 39.5, colour = ink, fill = "#F7F5EF") +
    txt(65.25, 35.4, "Gradient boosting", size = s75, colour = ink, fontface = "bold") +
    txt(65.25, 32.0, "decision trees, xgboost", size = s7) +
    rtxt(65.25, 27.6, "<span style='color:#333333'><i>Size and shape</i></span>: nucleus area,", size = s7) +
    txt(65.25, 24.4, "axis ratio, compactness,", size = s7) +
    txt(65.25, 21.2, "minor axis, cell area", size = s7) +
    rtxt(65.25, 16.8, "<span style='color:#333333'><i>DAPI (core-scaled)</i></span>: mean,", size = s7) +
    txt(65.25, 13.6, "minimum, maximum, SD, total", size = s7) +
    txt(65.25, 5.9, "Tested in the other 30%", size = s7) +
    txt(65.25, 2.6, sprintf("of the cores: AUC %.2f", auc), size = s7)
  # application to CD4+ cells
  p <- p + annotate("segment", x = 87.5, xend = 87.5, y = 1, yend = 45.5, linewidth = 0.3, colour = "grey70",
                    linetype = "22") +
    txt(90, 43.6, "Application", size = s85, hjust = 0, colour = ink, fontface = "bold") +
    box(90, 116, 17, 32, colour = ink2, lw = 0.4)
  p <- icon(p, 95.3, 24.5, "C", s = 0.8)
  p <- p + txt(100.4, 26.6, "bold('CD4'^'+'~'cells')", size = s75, hjust = 0, parse = TRUE, colour = ink) +
    txt(100.4, 22.3, "no CD8", size = s7, hjust = 0) +
    arr(116, 27.5, 121, 32) + arr(116, 21.5, 121, 17) +
    box(121.5, 146, 25.5, 39.5, colour = col_T, fill = "#F3EBD3", lw = 0.5) +
    box(121.5, 146, 9.5, 23.5, colour = col_excl, fill = "#E6E6E6", lw = 0.5) +
    txt(133.75, 36.0, "italic(p) > 0.5", size = s75, parse = TRUE, colour = ink) +
    txt(133.75, 32.5, "lymphocyte-like", size = s7) + txt(133.75, 29.0, "counted as T cell", size = s7) +
    txt(133.75, 20.0, "italic(p) <= 0.5", size = s75, parse = TRUE, colour = ink) +
    txt(133.75, 16.5, "non-lymphocyte-like", size = s7) + txt(133.75, 13.0, "excluded", size = s7) +
    txt(118, 5.9, "italic(p)*','~'predicted probability of'", size = s7, parse = TRUE) +
    txt(118, 2.6, "a lymphocyte-like nucleus", size = s7)
  p + labs(tag = tag) + th + blank_axes + theme(plot.margin = margin(0, 0, 0, 0))
}

## ---- A: nucleus schematic
panel_nucleus <- function(tag = "a") {
  cyN <- 1.2; aM <- sqrt(22 * 1.55 / pi); bM <- sqrt(22 / (pi * 1.55))
  spT <- spots(0, cyN, c(-1.1, 0.9, 0.2, -0.3, 1.3), c(0.8, 0.9, -1.2, -0.1, -0.5), c(0.42, 0.38, 0.4, 0.34, 0.3))
  spM <- spots(10.5, cyN, c(-1.9, 1.6, 0.4, -0.6), c(0.4, -0.5, 0.9, -0.8), c(0.36, 0.32, 0.3, 0.26))
  # one text block, so the line spacing stays the same whatever the panel scale
  # the same words and line breaks as the model box of the classifier schematic, group names in italics
  metrics <- paste("<b>Classifier measures</b>",
                   "<span style='color:#333333'><i>Size and shape</i></span>: nucleus area,",
                   "axis ratio, compactness,", "minor axis, cell area", "(dashed outline)",
                   "<span style='color:#333333'><i>DAPI (core-scaled)</i></span>: mean,",
                   "minimum, maximum, SD, total", sep = "<br>")
  ggplot() +
    geom_polygon(data = circ(0, cyN, 4.0), aes(x, y), fill = NA, colour = "grey55", linewidth = 0.3, linetype = "22") +
    geom_polygon(data = circ(10.5, cyN, 4.6), aes(x, y), fill = NA, colour = "grey55", linewidth = 0.3, linetype = "22") +
    geom_polygon(data = ell(0, cyN, 22, 1.12), aes(x, y), fill = nucT, colour = ink, linewidth = 0.3) +
    geom_polygon(data = spT, aes(x, y, group = id), fill = nucT2, colour = NA) +
    geom_polygon(data = ell(10.5, cyN, 22, 1.55), aes(x, y), fill = nucM, colour = ink, linewidth = 0.3) +
    geom_polygon(data = spM, aes(x, y, group = id), fill = nucM2, colour = NA) +
    annotate("segment", x = 10.5 - aM, xend = 10.5 + aM, y = cyN, yend = cyN, linewidth = 0.25, colour = ink,
             arrow = arrow(ends = "both", length = unit(0.8, "mm"), type = "open")) +
    annotate("segment", x = 10.5, xend = 10.5, y = cyN - bM, yend = cyN + bM, linewidth = 0.25, colour = ink,
             arrow = arrow(ends = "both", length = unit(0.8, "mm"), type = "open")) +
    annotate("text", x = 0, y = 6.6, label = "CD8^'+'~'T cell'", parse = TRUE, size = pt(7.5), colour = ink,
             family = "Arial") +
    annotate("text", x = 10.5, y = 6.6, label = "Macrophage", size = pt(7.5), colour = ink, family = "Arial") +
    annotate("text", x = c(-0.4, 10.9), y = -4.4, label = c("Brighter, rounder", "Paler, more oval"), size = pt(7),
             colour = ink2, family = "Arial") +
    annotate("richtext", x = -5.4, y = -6.1, label = metrics, hjust = 0, vjust = 1, size = pt(7), colour = ink2,
             family = "Arial", fill = NA, label.colour = NA, label.padding = unit(0, "pt"), lineheight = 1.2) +
    coord_fixed(xlim = c(-5.5, 16), ylim = c(-17.4, 7.8), expand = FALSE, clip = "off") +
    labs(title = "Nucleus, schematic", tag = tag) + th + blank_axes
}

## ---- A: DAPI and axis ratio in the stroma by cell group
panel_box <- function(v, xlab, title, lim, show_y = TRUE) {
  s <- bstat[feature == v & g %in% gshow][, g := factor(g, levels = rev(gshow))][, lin := fifelse(g == "Macrophages", "M", "T")]
  p <- ggplot(s, aes(y = g)) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.3, colour = ink) +
    geom_crossbar(aes(x = md, xmin = q1, xmax = q3, fill = lin), width = 0.6, linewidth = 0.3, colour = ink,
                  middle.linewidth = 0.6, orientation = "y") +
    scale_fill_manual(values = c(T = col_T, M = col_mac), guide = "none") +
    y_groups + coord_cartesian(xlim = lim) +
    labs(x = xlab, y = NULL, title = title) + th + grid_x + y_md
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

## ---- C: lymphocyte-like nuclei by cell group
panel_lymph <- function(tag = "c", show_y = TRUE, title = "Nucleus classifier") {
  d <- lstat[g %in% gshow][, g := factor(g, levels = rev(gshow))][, lin := fifelse(g == "Macrophages", "M", "T")]
  p <- ggplot(d, aes(pct, g, fill = lin)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.0f", pct)), hjust = -0.3, size = pt(7.5), colour = ink2, family = "Arial") +
    scale_fill_manual(values = c(T = col_T, M = col_mac), guide = "none") +
    y_groups +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
    annotate("text", x = 98, y = 0.62, label = sprintf("AUC %.2f", auc), hjust = 1, vjust = 0, size = pt(7),
             colour = ink2, family = "Arial") +
    labs(x = "Lymphocyte-like nucleus (%)", y = NULL, title = title, tag = tag) + th + grid_x + y_md
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

## ---- D: CD4+ cells kept as T cells or excluded
panel_counts <- function(tag = "d") {
  d <- rbind(lstat[g %in% c("CD4+", "CD4+ marker+"), .(g, cls = "Lymphocyte-like nucleus, counted as T cell", cells = lymph, n)],
             lstat[g %in% c("CD4+", "CD4+ marker+"), .(g, cls = "Non-lymphocyte-like nucleus, excluded", cells = n - lymph, n)])
  d[, `:=`(g = factor(g, levels = c("CD4+ marker+", "CD4+")),
           cls = factor(cls, levels = c("Non-lymphocyte-like nucleus, excluded", "Lymphocyte-like nucleus, counted as T cell")),
           lab = sprintf("%s\n(%.0f%%)", formatC(cells, format = "d", big.mark = ","), 100 * cells / n))]
  ggplot(d, aes(cells / 1000, g, fill = cls)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab, colour = cls), position = position_stack(vjust = 0.5), size = pt(7), family = "Arial",
              lineheight = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = c("Lymphocyte-like nucleus, counted as T cell" = col_T,
                                 "Non-lymphocyte-like nucleus, excluded" = col_excl), name = NULL,
                      breaks = c("Lymphocyte-like nucleus, counted as T cell", "Non-lymphocyte-like nucleus, excluded")) +
    guides(fill = guide_legend(nrow = 2)) +
    scale_colour_manual(values = c("Lymphocyte-like nucleus, counted as T cell" = ink,
                                   "Non-lymphocyte-like nucleus, excluded" = "white")) +
    y_groups +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Cells (thousands)", y = NULL, title = "CD4<sup>+</sup> cells kept or excluded", tag = tag) +
    th + grid_x + y_md +
    theme(legend.position = "bottom", legend.justification = "left", legend.location = "plot",
          legend.margin = margin(0, 0, 0, 0),
          plot.title = ggtext::element_markdown(size = 8.5, colour = ink, margin = margin(b = 5)))
}

## ---- the figure
W <- 152; H <- 200                                   # mm, printed size
cohort_label <- theme(plot.title = element_textbox(family = "Arial", size = 7.5, colour = "white", fill = col_tag,
                                                   width = NULL, r = unit(2, "pt"), padding = margin(2, 6, 1.5, 6),
                                                   margin = margin(0, 0, 3, 0)))
row1 <- (panel_nucleus("a") | panel_box("dapi_z", "Core-scaled intensity", "DAPI", c(-1.6, 3.1)) |
           panel_box("NucAxisRatio", "Length to width", "Axis ratio", c(1, 2.1), show_y = FALSE)) +
  plot_layout(widths = c(44, 31, 30))
row2 <- panel_training("b")
row3 <- (panel_lymph("c", title = "Lymphocyte-like nuclei") | panel_counts("d")) +
  plot_layout(widths = c(36, 54))
fig <- wrap_elements(full = row1) / wrap_elements(full = row2) / wrap_elements(full = row3) +
  plot_layout(heights = c(73, 51, 68)) +
  plot_annotation(title = "DISCOVERY COHORT", theme = cohort_label + theme(plot.margin = margin(4, 4, 4, 4)))
f <- file.path(OUT, "Supplementary_Figure_S3")
ggsave(paste0(f, ".pdf"), fig, width = W, height = H, units = "mm", device = cairo_pdf)
ggsave(paste0(f, ".png"), fig, width = W, height = H, units = "mm", dpi = 300, bg = "white")
print(lstat[, .(g, n, lymph, pct = round(pct, 1))]); print(tot)
cat("written:", paste0(f, ".pdf"), "\n")
