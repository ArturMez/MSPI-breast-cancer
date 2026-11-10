########################
## 01 - Discovery cohort: M1- and M2-to-T cell affinity scores from cell coordinates
##
## Algorithm (as used for the published discovery scores):
##   * cell classes: T cells = CD4+ or CD8+; macrophage-marker cells = CD68 (cytoplasm mean) > 5
##     or CD163+; M1 = CD68 > 5 and CD163-; M2 = CD68 > 5 and CD163+. The reference pool holds the
##     macrophage-marker cells that are not T cells (M1 + M2 + CD68-CD163+ myeloid cells).
##   * per TMA core with > 1 T cell, >= 1 pool cell and >= 1 M1 (or M2) cell: the pool is randomly
##     subsampled to the M1 (or M2) count, 5 times. For every T cell, the nearest-neighbour
##     distance (NND) to the subsample and to the M1 (or M2) cells is recorded; when the nearest
##     cell lies closer than 1 um, the second nearest is used.
##   * per patient: NNDs of dist_cap um or more are set to missing, separately for pool and subset,
##     and delta = mean(subset NND) / (mean(pool NND) + mean(subset NND)).
##     delta = 0.5 means no preference; below 0.5 the subset lies closer to T cells than the pool.
##     Script 05 negates delta so that higher scores mean higher affinity.
##
## The published scores came from a run without a fixed random seed, so a rerun gives slightly
## different values. Script 01b reproduces the published scores exactly from the stored per-cell
## distance files.
##
## Inputs (./data/discovery/):
##   Merged_BIN_for_spatial.txt  one row per cell: Sample.Name (patient), Core.ID,
##                               Cell.X.Position.um, Cell.Y.Position.um, CD4, CD8 (0/1),
##                               CD68 (cytoplasm mean intensity), CD163 (0/1)
##   Tissue_category_areas.txt   core-to-patient map (Sample.Name, Core.ID)
## Outputs (./output/discovery_affinity/):
##   nnALL_M1_to_CD3.txt, nnALL_M2_to_CD3.txt   per-cell distances (core_ID, MF_to_CD3 = pool,
##                                              Mx_to_CD3 = subset), 5 rows per T cell
##   discovery_distances_mean_M1_to_CD3.txt, discovery_distances_mean_M2_to_CD3.txt
##                                              per-patient scores (AnnUnit, delta_CD3)
########################

suppressMessages({
  library(data.table)
  library(FNN)
})

data_dir   <- "./data/discovery"
output_dir <- "./output/discovery_affinity"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed       <- 12345   # fixed seed for the random pool subsamples
n_resample <- 5       # pool subsamples per core
min_radius <- 1       # um; a nearer neighbour is treated as the same position
dist_cap   <- 25      # um; NNDs at or above this are not used in the per-patient score

cells <- fread(file.path(data_dir, "Merged_BIN_for_spatial.txt"),
               select = c("Sample.Name", "Core.ID", "Cell.X.Position.um", "Cell.Y.Position.um",
                          "CD4", "CD8", "CD68", "CD163"))
cells[, Tcell := CD4 == 1 | CD8 == 1]
cells[, myeloid := CD68 > 5 | CD163 == 1]
cells[, M1 := CD68 > 5 & CD163 == 0]
cells[, M2 := CD68 > 5 & CD163 == 1]
cat("Cells:", nrow(cells), "| cores:", uniqueN(cells$Core.ID), "\n")

nnd <- function(ref, query) {
  nn <- get.knnx(ref, query, k = 2)
  ifelse(nn$nn.dist[, 1] < min_radius, nn$nn.dist[, 2], nn$nn.dist[, 1])
}

by_core <- split(cells, by = "Core.ID", keep.by = TRUE)

compute_nn <- function(subset_col) {
  set.seed(seed)
  out <- vector("list", length(by_core))
  for (i in seq_along(by_core)) {
    d <- by_core[[i]]
    t_xy <- as.matrix(d[Tcell == TRUE, .(Cell.X.Position.um, Cell.Y.Position.um)])
    pool <- d[Tcell == FALSE & myeloid == TRUE]
    s_xy <- as.matrix(pool[get(subset_col) == TRUE, .(Cell.X.Position.um, Cell.Y.Position.um)])
    p_xy <- as.matrix(pool[, .(Cell.X.Position.um, Cell.Y.Position.um)])
    if (nrow(t_xy) <= 1 || nrow(p_xy) == 0 || nrow(s_xy) == 0) next
    pool_d <- unlist(lapply(seq_len(n_resample), function(z)
      nnd(p_xy[sample.int(nrow(p_xy), nrow(s_xy)), , drop = FALSE], t_xy)))
    sub_d <- nnd(s_xy, t_xy)
    out[[i]] <- data.table(core_ID = d$Core.ID[1], MF_to_CD3 = pool_d, Mx_to_CD3 = rep(sub_d, n_resample))
  }
  ## FNN returns NA for a missing second neighbour (one-cell reference); such rows are dropped,
  ## as in the original analysis
  na.omit(rbindlist(out))
}

tissue <- fread(file.path(data_dir, "Tissue_category_areas.txt"))
core_map <- unique(tissue[, .(AnnUnit = get(names(tissue)[1]), core_ID = Core.ID)])[!duplicated(core_ID)]

per_patient <- function(nn) {
  x <- merge(nn, core_map, by = "core_ID", all.x = TRUE)
  x[MF_to_CD3 >= dist_cap, MF_to_CD3 := NA]
  x[Mx_to_CD3 >= dist_cap, Mx_to_CD3 := NA]
  x[, .(delta_CD3 = mean(Mx_to_CD3, na.rm = TRUE) /
          (mean(MF_to_CD3, na.rm = TRUE) + mean(Mx_to_CD3, na.rm = TRUE))), by = AnnUnit]
}

for (ct in c("M1", "M2")) {
  cat("Computing", ct, "to T cell distances...\n")
  nn <- compute_nn(ct)
  fwrite(nn, file.path(output_dir, sprintf("nnALL_%s_to_CD3.txt", ct)), sep = "\t")
  sc <- per_patient(nn)
  fwrite(sc, file.path(output_dir, sprintf("discovery_distances_mean_%s_to_CD3.txt", ct)), sep = "\t")
  cat("  cores used:", uniqueN(nn$core_ID), "| patients with a score:", sum(!is.na(sc$delta_CD3)), "\n")
}
cat("Done. Outputs in", output_dir, "\n")
