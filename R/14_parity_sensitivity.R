########################
## 14 - Discovery cohort: parity as a sensitivity analysis of the clinical multivariable models
##
## Parity (parous or nulliparous, SAV variable Parity_1) in the discovery MSPI analysis set (patients with both
## affinity scores and DDFS time > 0, as in script 05). MSPI is compared between parous and nulliparous patients
## (Wilcoxon test), and parity is added to the clinical multivariable Cox models of Figure 2e (luminal) and
## Figure 2h (young early-stage luminal: stage I-IIb, age < 40): MSPI per SD within the subgroup, adjusted for
## age group, TNM stage, systemic treatment and diagnosis period, fitted on the patients with known parity
## without and with parity.
##
## Inputs : ./data/discovery/ SAV, discovery_distances_mean_M1/M2_to_CD3_nucleus.txt
## Outputs: ./output/parity/parity_sensitivity.tsv
## Usage  : Rscript R/14_parity_sensitivity.R
########################

suppressMessages({library(data.table); library(survival); library(foreign)})
disc_dir <- "./data/discovery"
OUT <- "./output/parity"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

sav <- list.files(disc_dir, pattern = "workingr_832_updated_211024\\.sav$", full.names = TRUE)[1]
cl <- suppressWarnings(read.spss(sav, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
                                 trim.factor.names = TRUE, trim_values = TRUE, reencode = NA))
if (!"Subtype_BCY21" %in% names(cl)) names(cl)[names(cl) == "Subtype_Discovery21"] <- "Subtype_BCY21"
cl$status_num <- ifelse(as.character(cl$DDFScens2021) == "Distant recurrence" |
                        as.character(cl$BCcens2021) == "Breast cancer death", 1L, 0L)
cl$DDFStime2021 <- as.numeric(cl$DDFStime2021)
cl <- cl[as.character(cl$Subtype_BCY21) != "Unclassified", ]
cl$AnnUnit <- as.character(cl$AnnUnit)
cl$parous <- factor(ifelse(as.character(cl$Parity_1) == "Missing", NA, as.character(cl$Parity_1)),
                    levels = c("No", "Yes"), labels = c("Nulliparous", "Parous"))
cl <- as.data.table(cl[, c("AnnUnit", "status_num", "DDFStime2021", "Subtype_BCY21", "Stage_TNM_44",
                           "Age_cat_2", "Syst_treat_NE", "diadtm_final", "parous")])

rd <- function(f) fread(file.path(disc_dir, f))[, .(AnnUnit = as.character(AnnUnit), d = delta_CD3)]
sc <- merge(setnames(rd("discovery_distances_mean_M1_to_CD3_nucleus.txt"), "d", "M1d"),
            setnames(rd("discovery_distances_mean_M2_to_CD3_nucleus.txt"), "d", "M2d"), by = "AnnUnit")
d <- merge(cl, sc, by = "AnnUnit")
d <- d[!is.na(M1d) & !is.na(M2d) & !is.na(DDFStime2021) & DDFStime2021 > 0]
d[, MSPI := (-M2d) - (-M1d)]                       # higher score = closer; MSPI = M2 score - M1 score
cat("MSPI analysis set:", nrow(d), "patients,", sum(d$status_num), "events\n")

res <- list(data.table(analysis = "Parity known", subgroup = "All", n = sum(!is.na(d$parous)), total = nrow(d),
                       parous = sum(d$parous == "Parous", na.rm = TRUE),
                       nulliparous = sum(d$parous == "Nulliparous", na.rm = TRUE)))
lum <- c("Lum A", "Lum B"); early <- c("Stage I", "Stage IIa", "Stage IIb")
sg <- list("All" = d, "Luminal" = d[Subtype_BCY21 %in% lum],
           "Young early luminal" = d[Stage_TNM_44 %in% early & Subtype_BCY21 %in% lum & Age_cat_2 == "<40"])
for (g in c("All", "Young early luminal")) {
  x <- sg[[g]][!is.na(parous)]
  res[[length(res) + 1]] <- data.table(analysis = "MSPI, parous vs nulliparous (Wilcoxon p)", subgroup = g,
                                       n = nrow(x), p = wilcox.test(MSPI ~ parous, data = x)$p.value)
}
hr <- function(cx) {
  s <- summary(cx)
  data.table(HR = s$conf.int["MSPI_sd", 1], CI_lo = s$conf.int["MSPI_sd", 3], CI_hi = s$conf.int["MSPI_sd", 4],
             p = s$coefficients["MSPI_sd", "Pr(>|z|)"])
}
for (g in c("Luminal", "Young early luminal")) {
  x <- copy(sg[[g]])
  x[, MSPI_sd := as.numeric(scale(MSPI))]             # per SD within the subgroup, as in script 05
  clin <- c("Age_cat_2", "Stage_TNM_44", "Syst_treat_NE", "diadtm_final")
  clin <- clin[sapply(clin, function(v) uniqueN(na.omit(x[[v]])) > 1)]
  x <- x[complete.cases(x[, c("MSPI_sd", clin), with = FALSE]) & !is.na(parous)]
  f0 <- as.formula(paste("Surv(DDFStime2021, status_num) ~ MSPI_sd +", paste(clin, collapse = " + ")))
  res[[length(res) + 1]] <- cbind(data.table(analysis = "Clinical model, MSPI per SD", subgroup = g, n = nrow(x),
                                             events = sum(x$status_num)), hr(coxph(f0, data = x)))
  res[[length(res) + 1]] <- cbind(data.table(analysis = "Clinical model + parity, MSPI per SD", subgroup = g,
                                             n = nrow(x), events = sum(x$status_num)),
                                  hr(coxph(update(f0, . ~ . + parous), data = x)))
}
res <- rbindlist(res, fill = TRUE)
print(res)
fwrite(res, file.path(OUT, "parity_sensitivity.tsv"), sep = "\t")
cat("written:", file.path(OUT, "parity_sensitivity.tsv"), "\n")
