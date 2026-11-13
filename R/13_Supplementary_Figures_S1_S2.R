########################
## 13 - Supplementary Figures S1 and S2: cell density against contact fractions and affinity scores
##
## S1, exploratory cohort (17 whole sections): densities of T cells, M1 and M2 macrophages and of T cells in
##     contact with M1 or M2 macrophages (counts / tissue area, cells/mm2), against the per-section affinity
##     computed from the nearest-neighbour distance files (distances of 25 um or more set to missing,
##     mean(pool distance) / (mean(subset distance) + mean(pool distance))). Spearman correlations.
## S2, discovery cohort: per-patient M1-to-T and M2-to-T affinity scores (delta, 0.5 = no preference;
##     python/01_discovery_affinity_nucleus.py) against T-cell density (R/01c) and M1 and M2 densities (SAV),
##     in every patient with the score and complete data. Spearman correlations.
##
## Inputs : ./data/exploratory/ia_areas_exploratory_images.tsv (or ia_areas_11062025_images.tsv)
##          ./data/exploratory/nnALL_M1_to_CD3_exploratory.txt, nnALL_M2_to_CD3_exploratory.txt (or *_MMP285.txt)
##          ./data/discovery/ SAV, discovery_distances_mean_M1/M2_to_CD3_nucleus.txt,
##          discovery_immune_variables_with_CD3_nucleus.txt
## Outputs: ./output/supplementary_figures/Supplementary_Figure_S1.pdf/.png, Supplementary_Figure_S2.pdf/.png,
##          S1_spearman.tsv, S2_spearman.tsv
## Both figures are drawn at the text width of the supplement (152 mm), so the font sizes set here are the printed
## sizes (axis text 7.5 pt, axis titles 8 pt, as Supplementary Figures S3 and S4), with lowercase panel letters as in
## the legends. Each group of four panels is a 2 x 2 block: M1 above M2, T-cell density left, macrophage density right.
## Usage  : Rscript R/13_Supplementary_Figures_S1_S2.R
########################

suppressMessages({library(data.table); library(dplyr); library(ggplot2); library(patchwork); library(foreign)})

expl_dir <- "./data/exploratory"
disc_dir <- "./data/discovery"
OUT <- "./output/supplementary_figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
find_input <- function(dir, ...) {
  for (f in c(...)) { p <- file.path(dir, f); if (file.exists(p)) return(p) }
  stop("Cannot find any of: ", paste(c(...), collapse = ", "), " in ", dir)
}
M1C <- "#7BA3C4"; M2C <- "#C47A7A"; INK <- "#333333"; AX <- "#4D4D4D"; GRID <- "#EBEBEB"

theme_paper <- function() {
  theme_minimal(base_size = 8, base_family = "Arial") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          axis.text = element_text(colour = AX, size = 7.5), axis.title = element_text(colour = INK, size = 8),
          plot.title = element_text(colour = INK, size = 8.5, face = "bold"),
          plot.subtitle = element_text(colour = AX, size = 7.5),
          axis.line = element_line(colour = INK, linewidth = 0.3), axis.ticks = element_line(colour = INK, linewidth = 0.3),
          plot.tag = element_text(face = "bold", size = 12, colour = "black"), plot.margin = margin(4, 10, 6, 4))
}

## p-value as plotmath text: p < 0.001 written as a x 10^-b with `sig` significant digits
p_expr <- function(p, sig = 2) {
  if (p == 0) return("italic(p) < 0.001")
  if (p < 0.001) {
    e <- floor(log10(p)); m <- signif(p / 10^e, sig)
    m <- formatC(m, format = "f", digits = sig - 1)
    return(sprintf("italic(p) == %s %%*%% 10^%d", m, e))
  }
  sprintf("italic(p) == '%s'", if (sig >= 3 || p < 0.1) sprintf("%.3f", p) else sprintf("%.2f", p))
}

scatter <- function(d, x, y, xl, yl, col, sig = 2) {
  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "spearman"))
  rho <- round(unname(ct$estimate), 2)
  rho_txt <- if (rho == 0) "0.00" else formatC(rho, format = "f", digits = 2)   # no "-0.00"
  sub <- parse(text = sprintf("'Spearman'~rho == '%s'*','~%s", rho_txt, p_expr(ct$p.value, sig)))
  g <- ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK, fill = "#BDBDBD", alpha = 0.35,
                linewidth = 0.5) +
    geom_point(colour = col, alpha = if (nrow(d) > 100) 0.35 else 0.85, size = if (nrow(d) > 100) 0.8 else 1.6,
               stroke = 0) +
    labs(x = xl, y = yl, title = sub) + theme_paper() +
    theme(plot.title = element_text(colour = AX, size = 7.5, face = "plain"))
  list(plot = g, rho = unname(ct$estimate), p = ct$p.value, n = sum(complete.cases(d[, c(x, y)])))
}
xT <- "T cell density (cells/mm²)"; xM1 <- "M1 density (cells/mm²)"; xM2 <- "M2 density (cells/mm²)"

## ---------------- Supplementary Figure S1: exploratory cohort
areas <- read.table(find_input(expl_dir, "ia_areas_exploratory_images.tsv", "ia_areas_11062025_images.tsv"),
                    sep = "\t", header = TRUE, fill = TRUE)
names(areas) <- gsub("-", "_", gsub("[. ]+", "_", names(areas)))
names(areas)[names(areas) == "Name"] <- "core_ID"
um2 <- areas$AREA_ALL                                   # tissue area in um2, so x 1e6 gives cells/mm2
areas <- areas %>% mutate(d_CD3 = (CD3 + HCD3_BY_M1 + HCD3_BY_M2) / um2 * 1e6,
                          d_CD3_by_M1 = HCD3_BY_M1 / um2 * 1e6, d_CD3_by_M2 = HCD3_BY_M2 / um2 * 1e6,
                          d_M1 = M1 / um2 * 1e6, d_M2 = M2 / um2 * 1e6)
delta <- function(...) {
  nn <- fread(find_input(expl_dir, ...), sep = "\t", header = TRUE, fill = TRUE)
  nn$POS <- ifelse(nn$POS < 25, nn$POS, NA); nn$neg <- ifelse(nn$neg < 25, nn$neg, NA)
  nn <- na.omit(as.data.frame(nn))
  nn %>% group_by(core_ID) %>% summarise(delta = mean(neg) / (mean(POS) + mean(neg)), .groups = "drop")
}
s1 <- areas %>%
  left_join(rename(delta("nnALL_M1_to_CD3_exploratory.txt", "nnALL_M1_to_CD3_MMP285.txt"), delta_M1 = delta), by = "core_ID") %>%
  left_join(rename(delta("nnALL_M2_to_CD3_exploratory.txt", "nnALL_M2_to_CD3_MMP285.txt"), delta_M2 = delta), by = "core_ID")
stopifnot(nrow(s1) == 17)
p1 <- list(
  scatter(s1, "d_CD3", "d_CD3_by_M1", xT,  "T cells hugged by M1\n(cells/mm²)", M1C),
  scatter(s1, "d_M1",  "d_CD3_by_M1", xM1, "T cells hugged by M1\n(cells/mm²)", M1C),
  scatter(s1, "d_CD3", "d_CD3_by_M2", xT,  "T cells hugged by M2\n(cells/mm²)", M2C),
  scatter(s1, "d_M2",  "d_CD3_by_M2", xM2, "T cells hugged by M2\n(cells/mm²)", M2C),
  scatter(s1, "d_CD3", "delta_M1", xT,  "M1-to-T cell affinity", M1C),
  scatter(s1, "d_M1",  "delta_M1", xM1, "M1-to-T cell affinity", M1C),
  scatter(s1, "d_CD3", "delta_M2", xT,  "M2-to-T cell affinity", M2C),
  scatter(s1, "d_M2",  "delta_M2", xM2, "M2-to-T cell affinity", M2C))
row_title <- function(t) wrap_elements(grid::textGrob(t, x = 0, hjust = 0, gp = grid::gpar(fontsize = 8.5,
                                       fontface = "bold", fontfamily = "Arial", col = INK)))
P1 <- lapply(seq_along(p1), function(i) p1[[i]]$plot + labs(tag = letters[i]))
fig1 <- (row_title("Contact fractions vs cell density (density-dependent)") / (P1[[1]] | P1[[2]]) / (P1[[3]] | P1[[4]]) /
         row_title("Affinity score vs cell density (density-independent)") / (P1[[5]] | P1[[6]]) / (P1[[7]] | P1[[8]])) +
  plot_layout(heights = c(0.12, 1, 1, 0.12, 1, 1))
ggsave(file.path(OUT, "Supplementary_Figure_S1.pdf"), fig1, width = 152, height = 190, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Supplementary_Figure_S1.png"), fig1, width = 152, height = 190, units = "mm", dpi = 300, bg = "white")
r1 <- data.frame(panel = letters[1:8], rho = sapply(p1, `[[`, "rho"), p = sapply(p1, `[[`, "p"), n = sapply(p1, `[[`, "n"))
print(r1)
fwrite(r1, file.path(OUT, "S1_spearman.tsv"), sep = "\t")

## ---------------- Supplementary Figure S2: discovery cohort
sav <- list.files(disc_dir, pattern = "workingr_832_updated_211024\\.sav$", full.names = TRUE)[1]
sav <- suppressWarnings(read.spss(sav, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
                                  trim.factor.names = TRUE, trim_values = TRUE, reencode = NA))
if (!"Subtype_BCY21" %in% names(sav)) names(sav)[names(sav) == "Subtype_Discovery21"] <- "Subtype_BCY21"
sav <- sav[which(sav$Subtype_BCY21 != "Unclassified"), ]
base <- data.frame(AnnUnit = as.character(sav$AnnUnit), M1 = sav$M1, M2 = sav$M2)
cd3 <- fread(file.path(disc_dir, "discovery_immune_variables_with_CD3_nucleus.txt"))[, 1:2]
setnames(cd3, c("AnnUnit", "CD3")); cd3$AnnUnit <- as.character(cd3$AnnUnit)
sel <- function(f, dens) {
  sc <- fread(file.path(disc_dir, f))[, .(AnnUnit = as.character(AnnUnit), score = delta_CD3)]
  d <- base %>% left_join(as.data.frame(sc), by = "AnnUnit") %>% left_join(as.data.frame(cd3), by = "AnnUnit")
  d <- d[!is.na(d$score), c("AnnUnit", dens, "CD3", "score")]
  d[complete.cases(d), ]
}
s2m1 <- sel("discovery_distances_mean_M1_to_CD3_nucleus.txt", "M1")
s2m2 <- sel("discovery_distances_mean_M2_to_CD3_nucleus.txt", "M2")
cat("patients: M1-to-T", nrow(s2m1), " M2-to-T", nrow(s2m2), "\n")
p2 <- list(
  scatter(s2m1, "CD3", "score", xT,  "M1-to-T cell affinity score", M1C, sig = 3),
  scatter(s2m1, "M1",  "score", xM1, "M1-to-T cell affinity score", M1C, sig = 3),
  scatter(s2m2, "CD3", "score", xT,  "M2-to-T cell affinity score", M2C, sig = 3),
  scatter(s2m2, "M2",  "score", xM2, "M2-to-T cell affinity score", M2C, sig = 3))
P2 <- lapply(seq_along(p2), function(i) p2[[i]]$plot + labs(tag = letters[i]))
fig2 <- (P2[[1]] | P2[[2]]) / (P2[[3]] | P2[[4]])
ggsave(file.path(OUT, "Supplementary_Figure_S2.pdf"), fig2, width = 152, height = 100, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Supplementary_Figure_S2.png"), fig2, width = 152, height = 100, units = "mm", dpi = 300, bg = "white")
r2 <- data.frame(panel = letters[1:4], rho = sapply(p2, `[[`, "rho"), p = sapply(p2, `[[`, "p"), n = sapply(p2, `[[`, "n"))
print(r2)
fwrite(r2, file.path(OUT, "S2_spearman.tsv"), sep = "\t")
cat("written:", OUT, "\n")
