########################
## 10 - Figure 1, data panels a (discovery numbers), f, g and h
##
## Draws the data-dependent elements of Figure 1 at the page positions, scales, line widths, marker sizes, colours,
## opacities and fonts of the figure layout (measured from the vector content of the figure page). Everything else
## (panels b-e and i, the boxes and labels of panel a, the axes, labels and legend of panel h) comes from
## assets/Fig1_template.pdf in script 12.
##
## Panel a : discovery patient numbers of the MSPI analysis set (patients with both affinity scores and DDFS
##           time > 0, as in script 05) by merged subtype, early-stage luminal (stage I-IIb) and age < 40 / >= 40
## Panels f, g: Kaplan-Meier curves by the optimal cutpoint of the M1-to-T and M2-to-T affinity score, computed
##           exactly as in script 04 (surv_cutpoint, log-log confidence bands, survminer's log-rank p label)
## Panel h : M1 and M2 affinity HR per SD from the additive Cox model and the LRT p of the M1 x M2 interaction,
##           all-comers row of discovery_M1M2_interaction_data.csv (script 06)
##
## Inputs : ./data/discovery/ (SAV, affinity scores), ./output/<date>_MSPI_clinical/
## Outputs: ./output/fig1/Fig1_panels.pdf   one transparent A4 page (drawn 16 times larger), placed onto the
##                                           template by script 12
##          ./output/fig1/Fig1_numbers.tsv   every number drawn
## Usage  : Rscript R/10_Fig1_draw.R
########################

suppressMessages({
  library(foreign); library(data.table); library(survival); library(survminer); library(grid)
})
disc_dir <- "./data/discovery"
out_dir  <- "./output/fig1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
latest <- function(pat) { f <- sort(list.files("./output", pattern = pat, full.names = TRUE)); f[length(f)] }
CLI <- latest("_MSPI_clinical$")
inp <- c(m1 = "discovery_distances_mean_M1_to_CD3_nucleus.txt", m2 = "discovery_distances_mean_M2_to_CD3_nucleus.txt")
cat("clinical tables:", CLI, "\n")

## ------------------------------------------------------------------ data
sav <- list.files(disc_dir, pattern = "workingr_832_updated_211024\\.sav$", full.names = TRUE)[1]
cl <- suppressWarnings(read.spss(sav, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
                                 trim.factor.names = TRUE, trim_values = TRUE, reencode = NA))
if (!"Subtype_BCY21" %in% names(cl)) names(cl)[names(cl) == "Subtype_Discovery21"] <- "Subtype_BCY21"
cl <- cl[which(as.character(cl$Subtype_BCY21) != "Unclassified"), ]
m1 <- fread(file.path(disc_dir, inp[["m1"]])); m2 <- fread(file.path(disc_dir, inp[["m2"]]))

## panel a: the MSPI analysis set of script 05 (both scores, DDFS time > 0)
st5 <- ifelse(as.character(cl$DDFScens2021) == "Distant recurrence" | as.character(cl$BCcens2021) == "Breast cancer death", 1L, 0L)
a <- data.table(AnnUnit = cl$AnnUnit, ev = st5, time = as.numeric(cl$DDFStime2021), Subtype = as.character(cl$Subtype_BCY21),
                Stage = as.character(cl$Stage_TNM_44), Age = as.character(cl$Age_cat_2))
a <- a[AnnUnit %in% m1[!is.na(delta_CD3), AnnUnit] & AnnUnit %in% m2[!is.na(delta_CD3), AnnUnit] & !is.na(time) & time > 0]
lum <- a[Subtype %in% c("Lum A", "Lum B")]; her <- a[Subtype %in% c("Lum Her2", "Her2 exp")]; tnbc <- a[Subtype == "Triple neg"]
early <- lum[Stage %in% c("Stage I", "Stage IIa", "Stage IIb")]
young <- early[Age == "<40"]; old <- early[Age != "<40"]
stopifnot(nrow(lum) + nrow(her) + nrow(tnbc) == nrow(a), nrow(young) + nrow(old) == nrow(early))
cnt <- function(x) c(n = nrow(x), e = sum(x$ev))
A_N <- list(all = cnt(a), lum = cnt(lum), her = cnt(her), tnbc = cnt(tnbc), early = cnt(early), young = cnt(young), old = cnt(old))
print(sapply(A_N, identity))

## panels f and g: script 04 (patients with both scores in the score files, the 04 status coding)
cl$EVENT04 <- ifelse(cl$DDFScens2021 == "Distant recurrence", "DDFScens2021", as.character(cl$BCcens2021))
ids <- intersect(m1$AnnUnit, m2$AnnUnit)
km_fit <- function(dist) {
  d <- merge(data.frame(AnnUnit = cl$AnnUnit, DDFStime2021 = as.numeric(cl$DDFStime2021), EVENT = cl$EVENT04),
             as.data.frame(dist[dist$AnnUnit %in% ids, ]), by = "AnnUnit", all.x = TRUE)
  d <- d[!is.na(d$delta_CD3), ]
  st <- ifelse(d$EVENT == "Censored", 0, d$EVENT); st <- ifelse(d$EVENT == "DDFScens2021", 1, st)
  d$status <- suppressWarnings(as.numeric(st))
  rc <- surv_cutpoint(d, time = "DDFStime2021", event = "status", variables = "delta_CD3")
  cc <- surv_categorize(rc)
  cc$XXX <- ifelse(cc$delta_CD3 == "high", 0, 1)          # 1 = low delta = closer = high affinity
  fit <- survfit(Surv(DDFStime2021, status == 1) ~ XXX, data = cc, conf.type = "log-log")
  lr  <- survdiff(Surv(DDFStime2021, status == 1) ~ XXX, data = cc)
  p   <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  list(fit = fit, p = p, cut = rc$cutpoint$cutpoint, n = sum(fit$n), n_lo = fit$n[1], n_hi = fit$n[2],
       tmax = max(cc$DDFStime2021, na.rm = TRUE))
}
KM <- list(f = km_fit(m1), g = km_fit(m2))
## survminer's p label (surv_pvalue): "p < 0.0001", otherwise p to 2 significant digits
p_lab <- function(p) if (p < 1e-4) "p < 0.0001" else paste("p =", signif(p, 2))

## panel h: additive model, all-comers (script 06)
H <- fread(file.path(CLI, "discovery_M1M2_interaction_data.csv"))[Subgroup == "All-comers"]
stopifnot(nrow(H) == 1)

## ------------------------------------------------------------------ page, units and drawing helpers (as script 08)
S <- 16; PW <- 595.28; PH <- 841.89                         # A4 in pt, drawn S times larger
X  <- function(x) unit(x * S, "bigpts")
Y  <- function(y) unit((PH - y) * S, "bigpts")
LW <- function(pt) pt * S * 96 / 72                          # grid lwd (1/96 inch) for a width in pt
seg <- function(x0, y0, x1, y1, col = "black", w = 0.5, cap = "butt")
  grid.segments(X(x0), Y(y0), X(x1), Y(y1), gp = gpar(col = col, lwd = LW(w), lineend = cap))
circ <- function(cx, cy, dia, fill) grid.circle(X(cx), Y(cy), r = unit(dia / 2 * S, "bigpts"), gp = gpar(fill = fill, col = NA))
## text at a baseline origin (x, y) as in the published PDF, optionally scaled along the baseline by sx
gid <- 0
txt <- function(label, x, y, size = 4, col = "black", face = 1, hjust = 0, rot = 0, sx = 1) {
  gp <- gpar(fontsize = size * S, col = col, fontface = face, fontfamily = "Arial")
  if (abs(sx - 1) < 1e-3) { grid.text(label, X(x), Y(y), hjust = hjust, vjust = 0, rot = rot, gp = gp); return(invisible()) }
  wd <- convertWidth(grobWidth(textGrob(label, gp = gp)), "bigpts", TRUE); hh <- size * S * 2
  gid <<- gid + 1; nm <- paste0("t", gid); x0 <- x * S - hjust * wd * sx
  pushViewport(viewport(x = unit(x0, "bigpts"), y = Y(y), width = unit(wd, "bigpts"), height = unit(hh, "bigpts"), just = c("left", "bottom")))
  grid.define(textGrob(label, x = 0, y = 0, hjust = 0, vjust = 0, gp = gp), name = nm); popViewport()
  pushViewport(viewport(x = unit(x0, "bigpts"), y = Y(y), width = unit(wd * sx, "bigpts"), height = unit(hh, "bigpts"), just = c("left", "bottom")))
  grid.use(nm); popViewport()
}
nums <- list()
addn <- function(panel, what, value) nums[[length(nums) + 1]] <<- data.table(panel = panel, quantity = what, value = as.character(value))

cairo_pdf(file.path(out_dir, "Fig1_panels.pdf"), width = PW * S / 72, height = PH * S / 72, family = "Arial", bg = "transparent")
grid.newpage()

## ------------------------------------------------------------------ (a) discovery numbers, centred where the published ones are
## centre x (pt) and baseline of each published number line
## format: "ne" = "n = N, E <word> events", "n" = "n = N", "e" = "E events"
AT <- list(list("all", "ne", 83.444, 88.60, " DDFS"),
           list("lum", "n", 43.958, 111.33), list("lum", "e", 43.954, 116.13),
           list("her", "n", 84.558, 111.33), list("her", "e", 84.552, 116.13),
           list("tnbc", "n", 123.278, 111.33), list("tnbc", "e", 123.282, 116.13),
           list("early", "ne", 52.192, 141.64, ""),
           list("young", "n", 44.328, 164.65), list("young", "e", 44.322, 169.45),
           list("old", "n", 84.558, 165.32), list("old", "e", 84.552, 170.12))
for (r in AT) {
  v <- A_N[[r[[1]]]]
  lab <- switch(r[[2]], ne = sprintf("n = %d, %d%s events", v[["n"]], v[["e"]], r[[5]]), n = sprintf("n = %d", v[["n"]]),
                e = sprintf("%d events", v[["e"]]))
  txt(lab, r[[3]], r[[4]], 4, hjust = .5)
}
for (k in names(A_N)) addn("a", paste0(k, " n / events"), paste(A_N[[k]], collapse = " / "))

## ------------------------------------------------------------------ (f, g) Kaplan-Meier curves
## Panel frame from the published axes: time -1.5 to 31.5 years (ggsurvplot's limits 0-30 with 5% expansion),
## survival -0.05 to 1.05. Curves 1 pt, confidence bands at opacity 0.298 drawn over the curves, censor marks are
## Arial "+" glyphs (5.207 pt) centred on the curve, legend and log-rank label as published.
GREEN <- "#6AAD8E"; RED <- "#C47A7A"; BAND_A <- 0.298
KMP <- list(
  f = list(x0 = 205.1, leg = list(hi = c(262.00, 321.28, 261.67, 318.36, 263.07, 323.11, 269.34, 322.51),
                                  lo = c(262.03, 330.84, 261.70, 327.92, 263.10, 332.68, 269.42, 332.07)),
           lr = c(210.86, 325.22, 332.27), ytl = c(195.34), ytitle = c(192.17, 326.84), xtitle = 236.59,
           xtl = c(208.17, 234.82, 262.56, 290.29)),
  g = list(x0 = 327.41, leg = list(hi = c(383.51, 319.63, 383.18, 316.72, 384.57, 321.47, 390.84, 320.86),
                                   lo = c(383.54, 329.20, 383.21, 326.28, 384.60, 331.03, 390.92, 330.43)),
           lr = c(333.18, 325.22, 332.27), ytl = c(317.66), ytitle = c(314.12, 325.22), xtitle = 358.34,
           xtl = c(330.49, 357.14, 384.88, 412.61)))
SX <- (296.62 - 205.1) / 33; SY <- (347.98 - 270.45) / 1.1
XT_OFF <- c(0, 3.465, 6.93, 10.01, 13.47, 15.005, 16.545, 19.11, 21.955, 24.8, 26.53, 29.095)
PLUS_DX <- 1.5204; PLUS_DY <- 1.8754                         # origin of the Arial "+" (5.207 pt) from the point it marks:
                                                             # ink centre from the glyph metrics, y checked on the published marks
for (pn in c("f", "g")) {
  P <- KMP[[pn]]; k <- KM[[pn]]; fit <- k$fit
  xm <- function(t) P$x0 + (t + 1.5) * SX; ym <- function(s) 347.98 - (s + 0.05) * SY
  strata <- rep(seq_along(fit$strata), fit$strata)
  cols <- c(GREEN, RED)                                      # XXX = 0 low affinity, XXX = 1 high affinity
  ## curves (steps), then bands, then censor marks
  for (j in 1:2) {
    i <- strata == j
    t <- c(0, fit$time[i]); s <- c(1, fit$surv[i])
    xs <- c(t[1], rep(t[-1], each = 2)); ys <- c(rep(s[-length(s)], each = 2), s[length(s)])
    grid.lines(X(xm(xs)), Y(ym(ys)), gp = gpar(col = cols[j], lwd = LW(1), linejoin = "round", lineend = "butt"))
  }
  for (j in 1:2) {
    i <- strata == j
    t <- c(0, fit$time[i]); lo <- c(1, fit$lower[i]); up <- c(1, fit$upper[i])
    ok <- !(is.na(lo) | is.na(up)); t <- t[ok]; lo <- lo[ok]; up <- up[ok]
    xs <- rep(t, each = 2)[-1]; yl <- rep(lo, each = 2)[-2 * length(lo)]; yu <- rep(up, each = 2)[-2 * length(up)]
    grid.polygon(X(xm(c(xs, rev(xs)))), Y(ym(c(yu, rev(yl)))), gp = gpar(fill = cols[j], col = NA, alpha = BAND_A))
  }
  for (j in 1:2) {
    i <- strata == j & fit$n.censor > 0
    if (any(i)) grid.text(rep("+", sum(i)), X(xm(fit$time[i]) - PLUS_DX), Y(ym(fit$surv[i]) + PLUS_DY), hjust = 0, vjust = 0,
                          gp = gpar(fontsize = 5.207 * S, col = cols[j], fontfamily = "Arial"))
  }
  ## axes
  seg(P$x0, 270.45, P$x0, 347.98, "black", .5, "square")
  seg(P$x0, 347.98, P$x0 + 91.52, 347.98, "black", .5, "square")
  for (s in c(0, .25, .5, .75, 1)) seg(P$x0 - 1.22, ym(s), P$x0, ym(s), "black", .5)
  for (t in c(0, 10, 20, 30)) seg(xm(t), 347.98, xm(t), 349.19, "black", .5)
  for (q in 1:5) txt(c("1.00", "0.75", "0.50", "0.25", "0.00")[q], P$ytl, c(275.37, 292.99, 310.61, 328.23, 345.85)[q], 4)
  for (q in 1:4) txt(c("0", "10", "20", "30")[q], P$xtl[q], 352.95, 4)
  txt("Survival probability", P$ytitle[1], P$ytitle[2], 4, rot = 90)
  ## the published title has irregular letter spacing: each letter at its published offset
  for (q in seq_along(XT_OFF)) txt(substr("DDFS (years)", q, q), P$xtitle + XT_OFF[q], 359.95, 5)
  ## legend: line, band-coloured square, censor glyph, label
  for (j in 2:1) {
    L <- P$leg[[c("lo", "hi")[j]]]
    seg(L[1], L[2], L[1] + 5.171, L[2], cols[j], 1)
    grid.rect(X(L[3]), Y(L[4] + 5.837), width = unit(5.837 * S, "bigpts"), height = unit(5.837 * S, "bigpts"), just = c("left", "bottom"),
              gp = gpar(fill = cols[j], col = NA, alpha = BAND_A))
    txt("+", L[5], L[6], 5.207, cols[j])
    txt(c("low affinity", "high affinity")[j], L[7], L[8], 4)
  }
  txt("Log-rank", P$lr[1], P$lr[2], 4)
  txt(p_lab(k$p), P$lr[1], P$lr[3], 4)
  ## the band must stay clear of the legend block (it does in both versions; stop rather than overlap)
  bl <- min(fit$lower[fit$time >= (P$leg$hi[3] - P$x0) / SX - 1.5], na.rm = TRUE)
  if (ym(bl) > P$leg$hi[4] - 0.5) stop("panel ", pn, ": confidence band reaches the legend")
  addn(pn, "patients (low / high affinity)", sprintf("%d (%d / %d)", k$n, k$n_lo, k$n_hi))
  addn(pn, "optimal cutpoint of delta (0.5 = no preference)", sprintf("%.4f", k$cut))
  addn(pn, "log-rank p", sprintf("%.3g (label: %s)", k$p, p_lab(k$p)))
  addn(pn, "censor marks drawn", sum(fit$n.censor > 0))
}

## ------------------------------------------------------------------ (h) additive Cox model and interaction test
## HR axis of the published panel: x = 487.325 + 42.13 * (HR - 1); rows M1 (y 261.17) and M2 (y 264.25).
## LRT axis: x = 534.81 + 17.758 * (p - 1), lollipop from p = 0 at y 262.52. The circle is drawn in the paper's
## terracotta when the interaction is significant (p < 0.05) and in the published grey-green otherwise, as the
## colour rule of script 06 (significant / additive).
hx <- function(v) 487.325 + 42.13 * (v - 1)
px <- function(p) 534.81 + 17.758 * (p - 1)
rng <- range(hx(c(H$M1_lo, H$M1_hi, H$M2_lo, H$M2_hi)))
if (rng[1] < 461.92 || rng[2] > 513.51) stop("panel h: a confidence interval is outside the published axis")
for (r in list(list(H$M1_HR, H$M1_lo, H$M1_hi, 261.17, 1.02, 2.67, "#7BA3C4"), list(H$M2_HR, H$M2_lo, H$M2_hi, 264.25, 1.03, 3.99, RED))) {
  seg(hx(r[[2]]), r[[4]], hx(r[[3]]), r[[4]], "#4D4D4D", .5)
  for (e in c(r[[2]], r[[3]])) seg(hx(e), r[[4]] - r[[5]] / 2, hx(e), r[[4]] + r[[5]] / 2, "#4D4D4D", .5)
}
circ(hx(H$M1_HR), 261.17, 2.67, "#7BA3C4"); circ(hx(H$M2_HR), 264.25, 3.99, RED)
seg(px(0), 262.52, px(H$LRT_p), 262.52, "#8C8C8C", .5)
circ(px(H$LRT_p), 262.52, 6.14, if (H$LRT_p < .05) RED else "#92A89B")
addn("h", "M1-to-T affinity HR per SD (95% CI), additive model", sprintf("%.2f (%.2f-%.2f) p=%.3g", H$M1_HR, H$M1_lo, H$M1_hi, H$M1_p))
addn("h", "M2-to-T affinity HR per SD (95% CI), additive model", sprintf("%.2f (%.2f-%.2f) p=%.3g", H$M2_HR, H$M2_lo, H$M2_hi, H$M2_p))
addn("h", "LRT p, M1 x M2 interaction", sprintf("%.3g", H$LRT_p))
addn("h", "patients / events", sprintf("%d / %d", H$n, H$events))

invisible(dev.off())
fwrite(rbindlist(nums), file.path(out_dir, "Fig1_numbers.tsv"), sep = "\t")
cat("written:", file.path(out_dir, "Fig1_panels.pdf"), "and Fig1_numbers.tsv\n")
