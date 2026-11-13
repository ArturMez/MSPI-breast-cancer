########################
## 01c - Discovery cohort: T-cell densities with nucleus-supported T cells
##
## Counts the T-cell density (cells/mm2) of every patient with the T cells of this analysis (CD8+ cells, and
## CD4+ cells with a lymphocyte-like nucleus, p_lymph > 0.5 from script 00b), over the tissue areas of the
## patient's cores (all tissue, tumour and stroma compartments). Only the T-cell columns of the density file
## change: CD3 (all tissue), CD3.1 (tumour) and CD3.2 (stroma). The other densities are kept.
## As an input check, the script first recounts CD4+ or CD8+ cells and stops unless this reproduces the
## T-cell columns of the input density file, which were computed that way.
##
## Inputs:
##   ./output/discovery_nucleus/cells_for_affinity_nucleus.tsv and cells_with_nuclei.tsv (scripts 00a, 00b)
##   ./data/discovery/Tissue_category_areas.txt
##   ./data/discovery/BCY_immune_variables_with_CD3_for spatial.txt (cell densities per patient)
## Outputs (./output/discovery_affinity_nucleus/):
##   discovery_immune_variables_with_CD3_nucleus.txt   the density file with the new T-cell columns
##                                                     (copy it to ./data/discovery/)
########################

suppressMessages(library(data.table))
data_dir   <- "./data/discovery"
nuc_dir    <- "./output/discovery_nucleus"
output_dir <- "./output/discovery_affinity_nucleus"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
find_input <- function(dir, ...) {
  for (f in c(...)) { p <- file.path(dir, f); if (file.exists(p)) return(p) }
  stop("Cannot find any of: ", paste(c(...), collapse = ", "), " in ", dir)
}

cells <- fread(file.path(nuc_dir, "cells_for_affinity_nucleus.tsv"), select = c("Sample.Name", "Core.ID", "CD4", "CD8", "p_lymph"),
               colClasses = c(Sample.Name = "character"))
tc <- fread(file.path(nuc_dir, "cells_with_nuclei.tsv"), select = c("Core.ID", "Tissue.Category"))
stopifnot(nrow(tc) == nrow(cells), all(tc$Core.ID == cells$Core.ID))
cells[, Tissue.Category := tc$Tissue.Category]
rm(tc); invisible(gc())
cells[, `:=`(T_pub = CD4 == 1 | CD8 == 1, T_nuc = CD8 == 1 | (CD4 == 1 & !is.na(p_lymph) & p_lymph > 0.5))]
cnt <- cells[, .(pub_all = sum(T_pub), nuc_all = sum(T_nuc),
                 pub_tum = sum(T_pub & Tissue.Category == "Tumor"), nuc_tum = sum(T_nuc & Tissue.Category == "Tumor"),
                 pub_str = sum(T_pub & Tissue.Category == "Stroma"), nuc_str = sum(T_nuc & Tissue.Category == "Stroma")),
             by = .(AnnUnit = Sample.Name)]

area <- fread(file.path(data_dir, "Tissue_category_areas.txt"))
setnames(area, c("Sample Name", "Tissue Category", "Region.Area..mm2."), c("AnnUnit", "Tissue", "area"))
area <- area[, .(area_all = sum(area), area_tum = sum(area[Tissue == "Tumor"]), area_str = sum(area[Tissue == "Stroma"])),
             by = .(AnnUnit = as.character(AnnUnit))]
cnt <- merge(cnt, area, by = "AnnUnit", all.x = TRUE)

dens_file <- find_input(data_dir, "discovery_immune_variables_with_CD3.txt", "BCY_immune_variables_with_CD3_for spatial.txt")
dens <- fread(dens_file, colClasses = c(AnnUnit = "character"))
x <- merge(dens[, .(AnnUnit, CD3, CD3.1, CD3.2)], cnt, by = "AnnUnit", all.x = TRUE)
stopifnot(!anyNA(x$area_all))
same <- function(a, b) all((is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= 1e-9 * pmax(1, abs(a))) |
                             (is.infinite(a) & is.infinite(b)))
ok <- c(all = same(x$CD3, x$pub_all / x$area_all), tumour = same(x$CD3.1, x$pub_tum / x$area_tum),
        stroma = same(x$CD3.2, x$pub_str / x$area_str))
print(ok)
stopifnot(all(ok))
cat("input T-cell densities reproduced from the cells for", nrow(x), "patients\n")

x[, `:=`(CD3_nuc = nuc_all / area_all, CD3.1_nuc = nuc_tum / area_tum, CD3.2_nuc = nuc_str / area_str)]
out <- copy(dens)
out[x, `:=`(CD3 = i.CD3_nuc, CD3.1 = i.CD3.1_nuc, CD3.2 = i.CD3.2_nuc), on = "AnnUnit"]
fwrite(out, file.path(output_dir, "discovery_immune_variables_with_CD3_nucleus.txt"), sep = "\t")
cat("written:", file.path(output_dir, "discovery_immune_variables_with_CD3_nucleus.txt"), "\n")
