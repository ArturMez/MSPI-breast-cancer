########################
## 08 - Figure 2: statistics of its data panels, and the panels as plain plots
##
## Computes the numbers of Figure 2 that do not come from the tables of scripts 05 and 06, with the data preparation
## of script 05 (steps 2 to 4):
##   (a)    medians of the M1-to-T and M2-to-T affinity scores and the qMSPI quadrant counts, full cohort
##   (b)    T cell, M1-like and M2-like densities and MSPI by subtype (luminal, HER2-positive, triple-negative) with
##          Dunn's tests and Benjamini-Hochberg correction
##   (f, g) early-stage (stage I-IIb) and young (< 40 years) early-stage luminal patients: Kaplan-Meier curves by
##          MSPI quartile (cut within each subgroup, as in script 05) with the 4th vs 1st log-rank p, the four-group
##          log-rank p and the trend p, and the Type II Anova shares of the density-adjusted Cox model (pies)
## The C-indices (c, f, g), the MSPI hazard ratios by subtype (d), the clinical multivariable models (e, h) and the
## MSPI hazard ratios of the mini forests (f, g) are in the tables of scripts 05 and 06. The figure page was composed
## from these outputs by separate layout scripts, which are not part of this repository.
##
## Inputs : ./data/discovery/ (SAV, affinity scores, densities)
## Outputs: ./output/fig2_statistics/Fig2_statistics.tsv   panel, quantity, value
##          ./output/fig2_statistics/Fig2a_qMSPI_scatter.pdf
##          ./output/fig2_statistics/Fig2b_densities_MSPI_by_subtype.pdf
##          ./output/fig2_statistics/Fig2f_early_luminal_MSPI_quartiles.pdf (Kaplan-Meier curves and pie)
##          ./output/fig2_statistics/Fig2g_young_early_luminal_MSPI_quartiles.pdf (Kaplan-Meier curves and pie)
## Usage  : Rscript R/08_Fig2_statistics.R
########################

suppressMessages({
  library(foreign); library(data.table); library(survival); library(car); library(rstatix)
})
disc_dir <- "./data/discovery"
out_dir  <- "./output/fig2_statistics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inp <- c(m1 = "discovery_distances_mean_M1_to_CD3_nucleus.txt", m2 = "discovery_distances_mean_M2_to_CD3_nucleus.txt",
         dens = "discovery_immune_variables_with_CD3_nucleus.txt")

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

nums <- list()
addn <- function(panel, what, value) nums[[length(nums) + 1]] <<- data.table(panel = panel, quantity = what, value = value)
p_txt <- function(p) if (p < .001) sprintf("%.1e", p) else sprintf("%#.2g", p)
QC  <- c(Q1 = "#6AAD8E", Q2 = "#D4B86A", Q3 = "#A0A0A0", Q4 = "#C47A7A")   # qMSPI quadrants
MQC <- c(Q1 = "#5B7FA0", Q2 = "#8EABBF", Q3 = "#C3ABAD", Q4 = "#C47A7A")   # MSPI quartiles, low (blue) to high (red)
quad <- function(x, c1, c2) factor(fifelse(x$M1s >= c1 & x$M2s < c2, "Q1", fifelse(x$M1s >= c1 & x$M2s >= c2, "Q2",
                                  fifelse(x$M1s < c1 & x$M2s < c2, "Q3", "Q4"))), levels = names(QC))
GL <- c(L = "Luminal", H = "HER2+", T = "TNBC")

## ------------------------------------------------------------------ (a) qMSPI quadrants, full cohort
c1 <- median(s4$M1s); c2 <- median(s4$M2s)
s4[, Quad := quad(s4, c1, c2)]
addn("a", "medians of the M1-to-T / M2-to-T affinity (0.5 = no preference)", sprintf("%.3f / %.3f", -c1, -c2))
addn("a", "qMSPI n per quadrant, full cohort (Q1/Q2/Q3/Q4)", paste(table(s4$Quad), collapse = "/"))
pdf(file.path(out_dir, "Fig2a_qMSPI_scatter.pdf"), width = 5, height = 5)
plot(s4$M1s, s4$M2s, pch = 16, cex = .5, col = QC[as.character(s4$Quad)], xlab = "M1-to-T affinity score",
     ylab = "M2-to-T affinity score", main = sprintf("qMSPI quadrants, discovery cohort (n = %d)", nrow(s4)))
abline(v = c1, h = c2, lty = 2)
legend("topright", legend = c("Q1: M1-close/M2-far", "Q2: both close", "Q3: both far", "Q4: M2-close/M1-far"),
       col = QC, pch = 16, cex = .7, bty = "n")
invisible(dev.off())

## ------------------------------------------------------------------ (b) densities and MSPI by subtype
bx <- s4[!is.na(Sub3) & !is.na(Td) & !is.na(M1d) & !is.na(M2d)]
dunn <- function(v, g) { r <- as.data.table(dunn_test(data.frame(v = v, g = g), v ~ g, p.adjust.method = "BH"))
  r[, lab := fcase(p.adj <= 1e-4, "****", p.adj <= 1e-3, "***", p.adj <= .01, "**", p.adj <= .05, "*", default = "ns")]; r }
pdf(file.path(out_dir, "Fig2b_densities_MSPI_by_subtype.pdf"), width = 11, height = 3.6)
par(mfrow = c(1, 4), mar = c(5.5, 4.2, 3.2, 1))
for (nm in c("T cell", "M1", "M2")) {
  v <- c("T cell" = "Td", M1 = "M1d", M2 = "M2d")[[nm]]
  y <- sqrt(bx[[v]])                                            # square-root scale, as in the figure
  dn <- dunn(y, bx$Sub3)
  for (i in seq_len(nrow(dn)))
    addn("b", sprintf("Dunn (BH) %s, %s vs %s", nm, dn$group1[i], dn$group2[i]), sprintf("%.2g %s", dn$p.adj[i], dn$lab[i]))
  for (g in levels(bx$Sub3))
    addn("b", sprintf("median density %s, %s (cells/mm2)", nm, g), sprintf("%.0f", median(bx[[v]][bx$Sub3 == g])))
  boxplot(y ~ bx$Sub3, names = GL, outline = FALSE, ylab = "Density (square root of cells/mm2)", xlab = "",
          main = nm)
  mtext(paste(sprintf("%s vs %s %s", dn$group1, dn$group2, dn$lab), collapse = ", "), side = 1, line = 3.5, cex = .7)
}
dm <- dunn(bx$MSPI, bx$Sub3)
for (i in seq_len(nrow(dm)))
  addn("b", sprintf("Dunn (BH) MSPI, %s vs %s", dm$group1[i], dm$group2[i]), sprintf("%.2g %s", dm$p.adj[i], dm$lab[i]))
boxplot(bx$MSPI ~ bx$Sub3, names = GL, outline = FALSE, ylab = "MSPI", xlab = "", main = "MSPI")
mtext(paste(sprintf("%s vs %s %s", dm$group1, dm$group2, dm$lab), collapse = ", "), side = 1, line = 3.5, cex = .7)
invisible(dev.off())

## ------------------------------------------------------------------ (f, g) early and young early luminal
panel_fg <- function(x, panel, title, file) {
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
  addn(panel, "MSPI quartiles, log-rank 4th vs 1st", p_txt(p41))
  addn(panel, "MSPI quartiles, log-rank, four groups", p_txt(p4g))
  addn(panel, "MSPI quartiles, trend (Cox, quartile 1-4)", p_txt(ptr))
  ## pie: Type II Anova of the density-adjusted Cox model (share of the chi-square of each term)
  xx <- copy(x)[!is.na(Td) & !is.na(M1d) & !is.na(M2d)]; xx[, MSPI_sd := MSPI / sd(MSPI)]
  a2 <- car::Anova(coxph(Surv(DDFStime2021, status_num) ~ MSPI_sd + M1d + M2d + Td, data = xx), type = "II")
  pp <- data.table(v = c("MSPI", "d M1", "d M2", "d T cell"), chisq = a2[, 1], p = a2[, 3])
  pp[, pct := 100 * chisq / sum(chisq)]
  for (i in seq_len(nrow(pp))) addn(panel, paste0("pie, ", pp$v[i], " (share of chi-square, p)"), sprintf("%.1f%%, p=%.2g", pp$pct[i], pp$p[i]))
  pdf(file, width = 9, height = 4.2)
  par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 1))
  plot(fit, col = MQC, lwd = 1.5, mark.time = TRUE, xlab = "DDFS (years)", ylab = "Survival probability",
       main = sprintf("%s, n = %d", title, nrow(x)))
  legend("bottomleft", legend = sprintf("MSPI %s quartile", c("1st", "2nd", "3rd", "4th")), col = MQC, lwd = 1.5, bty = "n", cex = .8)
  legend("topright", legend = sprintf("4th vs 1st, log-rank p = %s", p_txt(p41)), bty = "n", cex = .8)
  pie(pp$pct, labels = sprintf("%s %.1f%%", pp$v, pp$pct), col = c("#8EABBF", "#CCD5DE", "#C3ABAD", "#E0CE62"),
      clockwise = TRUE, main = "Share of chi-square\n(density-adjusted Cox model)", cex.main = .9)
  invisible(dev.off())
}
panel_fg(early, "f", "Early luminal (DDFS)", file.path(out_dir, "Fig2f_early_luminal_MSPI_quartiles.pdf"))
panel_fg(young, "g", "Young early luminal (DDFS)", file.path(out_dir, "Fig2g_young_early_luminal_MSPI_quartiles.pdf"))

fwrite(rbindlist(nums), file.path(out_dir, "Fig2_statistics.tsv"), sep = "\t")
cat("written", out_dir, "\n")
