########################
## 01b - Discovery cohort: per-patient affinity scores from the stored per-cell distance files
##
## Reproduces the published discovery scores exactly. The per-cell nearest-neighbour distances of
## the published run (format of script 01: core_ID, MF_to_CD3 = pool NND, Mx_to_CD3 = M1 or M2 NND)
## are averaged per patient after dropping NNDs of dist_cap um or more, and
## delta = mean(subset NND) / (mean(pool NND) + mean(subset NND)).
##
## Inputs (./data/discovery/): nnALL_M1_to_CD3.txt, nnALL_M2_to_CD3.txt, Tissue_category_areas.txt,
##   and, if present, the published score files for the check
##   (BCY_distances_mean_M1_to_CD3.txt / BCY_distances_mean_M2_to_CD3.txt or the discovery_* names)
## Outputs (./output/discovery_affinity/): discovery_distances_mean_M1_to_CD3.txt and _M2_
########################

suppressMessages(library(data.table))

data_dir   <- "./data/discovery"
output_dir <- "./output/discovery_affinity"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dist_cap   <- 25   # um

tissue <- fread(file.path(data_dir, "Tissue_category_areas.txt"))
core_map <- unique(tissue[, .(AnnUnit = get(names(tissue)[1]), core_ID = Core.ID)])[!duplicated(core_ID)]

for (ct in c("M1", "M2")) {
  nn <- fread(file.path(data_dir, sprintf("nnALL_%s_to_CD3.txt", ct)))
  x <- merge(nn, core_map, by = "core_ID", all.x = TRUE)
  x[MF_to_CD3 >= dist_cap, MF_to_CD3 := NA]
  x[Mx_to_CD3 >= dist_cap, Mx_to_CD3 := NA]
  sc <- x[, .(delta_CD3 = mean(Mx_to_CD3, na.rm = TRUE) /
                (mean(MF_to_CD3, na.rm = TRUE) + mean(Mx_to_CD3, na.rm = TRUE))), by = AnnUnit]
  fwrite(sc, file.path(output_dir, sprintf("discovery_distances_mean_%s_to_CD3.txt", ct)), sep = "\t")
  pub <- file.path(data_dir, c(sprintf("BCY_distances_mean_%s_to_CD3.txt", ct),
                               sprintf("discovery_distances_mean_%s_to_CD3.txt", ct)))
  pub <- pub[file.exists(pub)]
  if (length(pub)) {
    p <- fread(pub[1])
    m <- merge(p, sc, by = "AnnUnit", suffixes = c(".published", ".recomputed"))
    cat(sprintf("%s: %d patients | max abs difference to %s: %.2g\n", ct, nrow(m), basename(pub[1]),
                max(abs(m$delta_CD3.published - m$delta_CD3.recomputed), na.rm = TRUE)))
  }
}
