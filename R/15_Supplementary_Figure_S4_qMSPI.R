########################
## 15 - Supplementary Figure S4: DDFS by qMSPI phenotype in early-stage luminal patients
##
## Kaplan-Meier curves of DDFS for the four qMSPI phenotypes in early-stage luminal patients (stage I-IIb) and in
## young early-stage luminal patients (stage I-IIb, age <40 years), with the log-rank p-value for Q4 vs Q1. The
## quadrants split the M1-to-T and M2-to-T affinity scores at their medians within each subgroup, as in script 05.
## These curves were panels f and g of the submitted Figure 2. In the revision, Figure 2f and g show the curves by
## MSPI quartile (script 08). Data preparation as in scripts 05 and 08.
##
## Inputs : ./data/discovery/ (SAV, affinity scores, T-cell density)
## Outputs: ./output/supplementary_figures/Supplementary_Figure_S4.pdf/.png (152 mm wide, the text width of the
##          supplement, so the font sizes set here are the printed sizes), S4_numbers.tsv
## Usage  : Rscript R/15_Supplementary_Figure_S4_qMSPI.R
########################

suppressMessages({library(foreign); library(data.table); library(survival); library(ggplot2); library(patchwork)})
disc_dir <- "./data/discovery"
OUT <- "./output/supplementary_figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------ data (as script 08)
sav <- list.files(disc_dir, pattern = "workingr_832_updated_211024\\.sav$", full.names = TRUE)[1]
cl <- suppressWarnings(read.spss(sav, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
                                 trim.factor.names = TRUE, trim_values = TRUE, reencode = NA))
if (!"Subtype_BCY21" %in% names(cl)) names(cl)[names(cl) == "Subtype_Discovery21"] <- "Subtype_BCY21"
cl$status_num <- ifelse(as.character(cl$DDFScens2021) == "Distant recurrence" | as.character(cl$BCcens2021) == "Breast cancer death", 1L, 0L)
cl$DDFStime2021 <- as.numeric(cl$DDFStime2021)
cl <- cl[as.character(cl$Subtype_BCY21) != "Unclassified", ]
cl <- as.data.table(cl)[, .(AnnUnit, status_num, DDFStime2021, Subtype = as.character(Subtype_BCY21),
                            Stage = as.character(Stage_TNM_44), Age = as.character(Age_cat_2))]
m1 <- fread(file.path(disc_dir, "discovery_distances_mean_M1_to_CD3_nucleus.txt"))
m2 <- fread(file.path(disc_dir, "discovery_distances_mean_M2_to_CD3_nucleus.txt"))
d <- merge(cl, m1[, .(AnnUnit, M1s = -delta_CD3)], by = "AnnUnit", all.x = TRUE)
d <- merge(d, m2[, .(AnnUnit, M2s = -delta_CD3)], by = "AnnUnit", all.x = TRUE)
s4 <- d[!is.na(M1s) & !is.na(M2s) & !is.na(DDFStime2021) & DDFStime2021 > 0]
early <- s4[Stage %in% c("Stage I", "Stage IIa", "Stage IIb") & Subtype %in% c("Lum A", "Lum B")]
young <- early[Age == "<40"]
cat("patients:", nrow(s4), "| early", nrow(early), sum(early$status_num), "| young", nrow(young), sum(young$status_num), "\n")

## ------------------------------------------------------------------ style (as Figure 2 and Supplementary Figure S3)
QC  <- c(Q1 = "#6AAD8E", Q2 = "#D4B86A", Q3 = "#A0A0A0", Q4 = "#C47A7A")
QL  <- c(Q1 = "Q1: M1-close/M2-far", Q2 = "Q2: both close", Q3 = "Q3: both far", Q4 = "Q4: M2-close/M1-far")
INK <- "#333333"; AX <- "#4D4D4D"
theme_s4 <- function() {
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(axis.text = element_text(colour = AX, size = 7.5), axis.title = element_text(colour = INK, size = 8.5),
          axis.line = element_line(colour = INK, linewidth = 0.3), axis.ticks = element_line(colour = INK, linewidth = 0.3),
          plot.title = element_text(colour = AX, size = 8.5, face = "italic"),
          plot.subtitle = element_text(colour = AX, size = 7.5, margin = margin(b = 4)),
          legend.title = element_text(colour = INK, size = 8, face = "bold"), legend.text = element_text(colour = INK, size = 7.5),
          legend.key.height = unit(3.4, "mm"), legend.key.width = unit(4, "mm"), legend.background = element_blank(),
          plot.tag = element_text(face = "bold", size = 12, colour = "#080504"), plot.margin = margin(4, 8, 4, 4))
}
quad <- function(x, c1, c2) factor(fifelse(x$M1s >= c1 & x$M2s < c2, "Q1", fifelse(x$M1s >= c1 & x$M2s >= c2, "Q2",
                                  fifelse(x$M1s < c1 & x$M2s < c2, "Q3", "Q4"))), levels = names(QC))
p_txt <- function(p) if (p < .001) sprintf("%.1e", p) else sprintf("%#.2g", p)

nums <- list()
km_panel <- function(x, title, sub, key) {
  x <- copy(x); x[, Quad := quad(x, median(x$M1s), median(x$M2s))]
  fit <- survfit(Surv(DDFStime2021, status_num) ~ Quad, data = x)
  p41 <- 1 - pchisq(survdiff(Surv(DDFStime2021, status_num) ~ Quad, data = x[Quad %in% c("Q1", "Q4")])$chisq, 1)
  p4g <- 1 - pchisq(survdiff(Surv(DDFStime2021, status_num) ~ Quad, data = x)$chisq, 3)
  nums[[key]] <<- data.table(panel = key, n = nrow(x), events = sum(x$status_num),
                             n_per_quadrant = paste(table(x$Quad), collapse = "/"),
                             events_per_quadrant = paste(x[, sum(status_num), keyby = Quad]$V1, collapse = "/"),
                             logrank_Q4_vs_Q1 = p41, logrank_four_groups = p4g)
  st <- data.table(time = fit$time, surv = fit$surv, n.censor = fit$n.censor,
                   Quad = factor(rep(sub("Quad=", "", names(fit$strata)), fit$strata), levels = names(QC)))
  steps <- st[, .(time = c(0, time), surv = c(1, surv)), by = Quad]
  ggplot() +
    geom_step(data = steps, aes(time, surv, colour = Quad), linewidth = 0.45, direction = "hv") +
    geom_point(data = st[n.censor > 0], aes(time, surv, colour = Quad), shape = 3, size = 0.9, stroke = 0.35,
               show.legend = FALSE) +
    annotate("text", x = 31, y = 0.42, label = paste0("Q4 vs Q1, log-rank p = ", p_txt(p41)), hjust = 1, vjust = 0,
             size = 7.5 / .pt, fontface = "italic", colour = AX, family = "Arial") +
    scale_colour_manual(values = QC, labels = QL, name = "qMSPI") +
    scale_x_continuous(breaks = seq(0, 30, 10), limits = c(0, 31), expand = expansion(mult = c(.02, .02))) +
    scale_y_continuous(breaks = seq(0, 1, .25), labels = function(v) sprintf("%.2f", v), limits = c(0, 1),
                       expand = expansion(mult = c(.02, .03))) +
    labs(x = "DDFS (years)", y = "Survival probability", title = title, subtitle = sub) +
    theme_s4() + theme(legend.position = c(.02, .02), legend.justification = c(0, 0))
}
pa <- km_panel(early, "Early Luminal (DDFS)", sprintf("Stage I-IIb, n = %d, %d events", nrow(early), sum(early$status_num)), "a")
pb <- km_panel(young, "Young Early Luminal (DDFS)",
               sprintf("Stage I-IIb, age <40 years, n = %d, %d events", nrow(young), sum(young$status_num)), "b")
fig <- (pa | pb) + plot_annotation(tag_levels = "a")   # lowercase panel letters, as in the legend
ggsave(file.path(OUT, "Supplementary_Figure_S4.pdf"), fig, width = 152, height = 72, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Supplementary_Figure_S4.png"), fig, width = 152, height = 72, units = "mm", dpi = 300, bg = "white")
fwrite(rbindlist(nums), file.path(OUT, "S4_numbers.tsv"), sep = "\t")
print(rbindlist(nums))
cat("written:", file.path(OUT, "Supplementary_Figure_S4.pdf"), "\n")
