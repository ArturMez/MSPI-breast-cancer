########################
## 10 - Figure 1a: discovery patient numbers
##
## Counts the patients and events of the MSPI analysis set of script 05 (patients with both affinity scores and DDFS
## time > 0) by merged subtype, early-stage luminal (stage I-IIb) and age < 40 / >= 40 years, the discovery numbers
## shown in Figure 1a. The Kaplan-Meier panels f and g are drawn by script 04, and the additive model of panel h by
## script 06 (discovery_M1M2_interaction_forest.pdf). The figure page was composed from these outputs by separate
## layout scripts, which are not part of this repository.
##
## Inputs : ./data/discovery/ (SAV, affinity scores)
## Outputs: ./output/fig1a/Fig1a_patient_numbers.tsv   group, patients, events
## Usage  : Rscript R/10_Fig1a_patient_numbers.R
########################

suppressMessages({ library(foreign); library(data.table) })
disc_dir <- "./data/discovery"
out_dir  <- "./output/fig1a"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inp <- c(m1 = "discovery_distances_mean_M1_to_CD3_nucleus.txt", m2 = "discovery_distances_mean_M2_to_CD3_nucleus.txt")

## ------------------------------------------------------------------ data
sav <- list.files(disc_dir, pattern = "workingr_832_updated_211024\\.sav$", full.names = TRUE)[1]
cl <- suppressWarnings(read.spss(sav, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
                                 trim.factor.names = TRUE, trim_values = TRUE, reencode = NA))
if (!"Subtype_BCY21" %in% names(cl)) names(cl)[names(cl) == "Subtype_Discovery21"] <- "Subtype_BCY21"
cl <- cl[which(as.character(cl$Subtype_BCY21) != "Unclassified"), ]
m1 <- fread(file.path(disc_dir, inp[["m1"]])); m2 <- fread(file.path(disc_dir, inp[["m2"]]))

## the MSPI analysis set of script 05 (both scores, DDFS time > 0)
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

GN <- c(all = "All patients", lum = "Luminal (A and B)", her = "HER2-positive", tnbc = "Triple-negative",
        early = "Early-stage luminal (stage I-IIb)", young = "Early-stage luminal, age < 40 years",
        old = "Early-stage luminal, age >= 40 years")
out <- data.table(group = GN[names(A_N)], patients = sapply(A_N, `[[`, "n"), events = sapply(A_N, `[[`, "e"))
fwrite(out, file.path(out_dir, "Fig1a_patient_numbers.tsv"), sep = "\t")
cat("written", file.path(out_dir, "Fig1a_patient_numbers.tsv"), "\n")
