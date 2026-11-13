########################
## 00b - Discovery cohort: nucleus classifier, lymphocyte-like or non-lymphocyte-like nucleus
##
## CD4 is also expressed by monocytes, macrophages and dendritic cells, so a CD4+ cell counts as a T cell
## only when its nucleus looks like a lymphocyte nucleus. The classifier is gradient boosting (xgboost) on
## nucleus area, axis ratio, compactness, minor axis, cell area and five nucleus DAPI measures (mean,
## minimum, maximum, standard deviation and total, each standardised per core by median and MAD, plus the
## standard deviation, maximum and minimum relative to the mean, clipped at the training 0.5th and 99.5th
## percentiles).
## Training: unambiguous CD8+ cells (CD8+, CD4-, no macrophage marker) against unambiguous CD68+ or CD163+
## cells (CD4-, CD8-, CD68 > 5 or CD163+), 150,000 cells per class from a random 70% of the cores (seed
## 20261110). The tree depth (4, 6 or 8) and the number of boosting rounds are chosen by early stopping on
## validation cores split off the training cores (15%), and the model is then refit on all training cells.
## Evaluation: the cells of the other 30% of cores (held-out AUC 0.80). p_lymph > 0.5 is read as
## lymphocyte-like.
##
## Inputs:
##   ./output/discovery_nucleus/cells_with_nuclei.tsv (script 00a)
##   inForm export as in script 00a (MSPI_INFORM_DIR), for the minimum, maximum, SD and total DAPI
## Outputs (./output/discovery_nucleus/):
##   cells_for_affinity_nucleus.tsv   one row per cell: Sample.Name, Core.ID, X, Y, CD4, CD8, CD68,
##                                    CD163, p_lymph (input of scripts 01_discovery_affinity_nucleus.py
##                                    and 01c)
##   nucleus_classifier_auc.tsv, nucleus_classifier_tuning.tsv, nucleus_classifier_xgb.ubj (the model),
##   nucleus_classifier_settings.rds (clipping limits, test and validation cores)
########################

suppressMessages({
  library(data.table)
  library(xgboost)
})
setDTthreads(8)
set.seed(20261110)

data_dir   <- "./data/discovery"
inform_dir <- Sys.getenv("MSPI_INFORM_DIR", file.path(data_dir, "inform"))
output_dir <- "./output/discovery_nucleus"

## 1. cells, groups and the train/test split by core
m <- fread(file.path(output_dir, "cells_with_nuclei.tsv"))
m_all_n <- nrow(m)
m[, rid := .I]
m[, ok := !is.na(NucArea) & NucArea > 0]
m[, myel := (CD68 > 5 | CD163 == 1)]
m[ok == TRUE, dapi_z := (NucDAPI - median(NucDAPI)) / mad(NucDAPI), by = Core.ID]
m[, grp := fcase(
  CD8 == 1 & CD4 == 0 & !myel, "CD8+ clean",
  CD4 == 1 & CD8 == 0 & !myel, "CD4+ clean",
  CD4 == 1 & CD8 == 1 & !myel, "CD4+CD8+ clean",
  CD4 == 1 & CD8 == 0 & myel, "CD4+ myeloid+",
  CD8 == 1 & CD4 == 0 & myel, "CD8+ myeloid+",
  CD4 == 1 & CD8 == 1 & myel, "CD4+CD8+ myeloid+",
  CD4 == 0 & CD8 == 0 & CD68 > 5 & CD163 == 0, "M1 clean",
  CD4 == 0 & CD8 == 0 & CD68 > 5 & CD163 == 1, "M2 clean",
  CD4 == 0 & CD8 == 0 & CD163 == 1, "Mx clean",
  default = "other")]
cores <- unique(m$Core.ID)
test_cores <- sample(cores, round(0.3 * length(cores)))
tr_id <- rbind(m[ok == TRUE & grp == "CD8+ clean" & !(Core.ID %in% test_cores)][sample(.N, min(.N, 150000))][, y := 1],
               m[ok == TRUE & grp %in% c("M1 clean", "M2 clean", "Mx clean") & !(Core.ID %in% test_cores)][
                 sample(.N, min(.N, 150000))][, y := 0])[, .(rid, y)]
te_id <- rbind(m[ok == TRUE & grp == "CD8+ clean" & Core.ID %in% test_cores][, y := 1],
               m[ok == TRUE & grp %in% c("M1 clean", "M2 clean", "Mx clean") & Core.ID %in% test_cores][, y := 0])[, .(rid, y)]

## 2. minimum, maximum, SD and total nucleus DAPI from the raw inForm export, joined as in script 00a
cols <- c(2, 8, 9, 20, 22, 23, 24)
nm   <- c("SampleName", "X", "Y", "DAPImin", "DAPImax", "DAPIsd", "DAPItot")
r <- rbind(fread(file.path(inform_dir, "BCY_merged_cell_seg_data.txt"), select = cols, col.names = nm),
           fread(file.path(inform_dir, "BCY_merged_additional_cell_seg_data.txt"), select = cols, col.names = nm))
r[, Core.ID := gsub("]_.*", "]", gsub("_Core...", "_[", SampleName))]
r[, `:=`(Xr = round(X, 1), Yr = round(Y, 1))]
r <- r[!is.na(Xr) & !is.na(Yr)]
dk <- r[, .N, by = .(Core.ID, Xr, Yr)][N > 1]
r <- r[!dk, on = .(Core.ID, Xr, Yr)]
m[, `:=`(Xr = round(Cell.X.Position.um, 1), Yr = round(Cell.Y.Position.um, 1))]
m[r, on = .(Core.ID, Xr, Yr), `:=`(DAPImin = i.DAPImin, DAPImax = i.DAPImax, DAPIsd = i.DAPIsd, DAPItot = i.DAPItot)]
rm(r); invisible(gc())
cat("cells:", nrow(m), "| with nucleus data:", sum(m$ok), "| with the DAPI statistics:", sum(!is.na(m$DAPIsd)), "\n")

## 3. features, clipped at the training 0.5th and 99.5th percentiles
zc <- function(x) (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE)
m[ok == TRUE, `:=`(sd_z = zc(DAPIsd), max_z = zc(DAPImax), min_z = zc(DAPImin), tot_z = zc(DAPItot)), by = Core.ID]
m[ok == TRUE, `:=`(cv = DAPIsd / NucDAPI, maxrel = DAPImax / NucDAPI, minrel = DAPImin / NucDAPI)]
newf <- c("sd_z", "max_z", "min_z", "tot_z", "cv", "maxrel", "minrel")
tr <- m[tr_id, on = "rid"]
lims <- lapply(setNames(newf, newf), function(f) quantile(tr[[f]], c(.005, .995), na.rm = TRUE))
for (f in newf) set(m, j = f, value = pmin(pmax(m[[f]], lims[[f]][1]), lims[[f]][2]))
tr <- m[tr_id, on = "rid"]
te <- m[te_id, on = "rid"]

## 4. model: gradient boosting, depth and rounds chosen on validation cores, then refit; held-out AUC
feats <- c("NucArea", "NucAxisRatio", "NucCompact", "NucMinor", "CellArea", "dapi_z", "sd_z", "cv", "maxrel",
           "minrel", "max_z", "min_z", "tot_z")
auc <- function(p, y) {
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
dmat <- function(d, lab = NULL) xgb.DMatrix(as.matrix(d[, ..feats]), label = lab, nthread = 8)
xparams <- function(depth) list(objective = "binary:logistic", eval_metric = "auc", eta = 0.05, max_depth = depth,
                                subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 50, nthread = 8)
set.seed(20261110)
val_cores <- sample(unique(tr$Core.ID), round(0.15 * uniqueN(tr$Core.ID)))
iv <- tr$Core.ID %in% val_cores
d_fit <- dmat(tr[!iv], tr$y[!iv]); d_val <- dmat(tr[iv], tr$y[iv])
tune <- rbindlist(lapply(c(4, 6, 8), function(dp) {
  set.seed(20261110 + dp)
  f <- xgb.train(params = xparams(dp), data = d_fit, nrounds = 3000, evals = list(val = d_val),
                 early_stopping_rounds = 50, verbose = 0)
  ev <- attributes(f)$evaluation_log
  data.table(depth = dp, rounds = ev$iter[which.max(ev$val_auc)], val_auc = round(max(ev$val_auc), 4))
}))
best <- tune[which.max(val_auc)]
set.seed(20261110)
fit <- xgb.train(params = xparams(best$depth), data = dmat(tr, tr$y), nrounds = best$rounds, verbose = 0)
res <- data.table(model = sprintf("gradient boosting (xgboost), tree depth %d, %d rounds", best$depth, best$rounds),
                  features = "nucleus and cell shape, mean, minimum, maximum, SD and total DAPI",
                  AUC_heldout_cores = round(auc(predict(fit, dmat(te)), te$y), 3))
print(tune); print(res)
fwrite(res, file.path(output_dir, "nucleus_classifier_auc.tsv"), sep = "\t")
fwrite(tune, file.path(output_dir, "nucleus_classifier_tuning.tsv"), sep = "\t")
xgb.save(fit, file.path(output_dir, "nucleus_classifier_xgb.ubj"))
saveRDS(list(clip_limits = lims, test_cores = test_cores, val_cores = val_cores, features = feats),
        file.path(output_dir, "nucleus_classifier_settings.rds"))

## 5. p_lymph for every cell with nucleus data (in chunks, to limit memory)
m[, p_lymph := NA_real_]
idx <- which(m$ok)
for (s in split(idx, ceiling(seq_along(idx) / 1e6)))
  set(m, i = s, j = "p_lymph", value = predict(fit, dmat(m[s])))
fwrite(m[, .(Sample.Name, Core.ID, X = Cell.X.Position.um, Y = Cell.Y.Position.um, CD4, CD8, CD68, CD163,
             p_lymph = round(p_lymph, 4))],
       file.path(output_dir, "cells_for_affinity_nucleus.tsv"), sep = "\t")
cat("cells written:", nrow(m), "of", m_all_n, "\n")
