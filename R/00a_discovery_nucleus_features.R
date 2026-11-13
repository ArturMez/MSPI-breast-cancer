########################
## 00a - Discovery cohort: inForm nucleus measures for every analysis cell
##
## The T-cell definition uses the nucleus of each CD4+ cell (script 00b). This script attaches the
## inForm nucleus measures (nucleus area, compactness, minor and major axis, axis ratio, mean DAPI, and
## cell area) from the raw inForm export to the cells of Merged_BIN_for_spatial.txt.
## Join key: Core.ID (the raw sample name turned into the core name used in Merged_BIN_for_spatial.txt)
## and the X/Y position rounded to 0.1 um. Raw keys that occur more than once are left out, so those
## cells have no nucleus data (70 of 9,338,966 cells).
##
## Inputs:
##   ./data/discovery/Merged_BIN_for_spatial.txt
##   inForm export (about 9 GB) in the folder named by the environment variable MSPI_INFORM_DIR
##   (default ./data/discovery/inform/): BCY_merged_cell_seg_data.txt and
##   BCY_merged_additional_cell_seg_data.txt
## Output (./output/discovery_nucleus/):
##   cells_with_nuclei.tsv   the cells of Merged_BIN_for_spatial.txt, same rows and order, with the
##                           nucleus measures
########################

suppressMessages(library(data.table))
setDTthreads(8)

data_dir   <- "./data/discovery"
inform_dir <- Sys.getenv("MSPI_INFORM_DIR", file.path(data_dir, "inform"))
output_dir <- "./output/discovery_nucleus"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## raw inForm columns: sample name, tissue category, X, Y, nucleus area, compactness, minor axis,
## major axis, axis ratio, mean DAPI, cell area
cols <- c(2, 3, 8, 9, 14, 16, 17, 18, 19, 21, 106)
nm   <- c("SampleName", "Tissue", "X", "Y", "NucArea", "NucCompact", "NucMinor", "NucMajor",
          "NucAxisRatio", "NucDAPI", "CellArea")
r <- rbind(fread(file.path(inform_dir, "BCY_merged_cell_seg_data.txt"), select = cols, col.names = nm),
           fread(file.path(inform_dir, "BCY_merged_additional_cell_seg_data.txt"), select = cols, col.names = nm))
cat("raw inForm rows:", nrow(r), "\n")
r[, Core.ID := gsub("]_.*", "]", gsub("_Core...", "_[", SampleName))]
r[, `:=`(Xr = round(X, 1), Yr = round(Y, 1))]
r <- r[!is.na(Xr) & !is.na(Yr)]
dk <- r[, .N, by = .(Core.ID, Xr, Yr)][N > 1]
cat("raw keys (Core.ID, X, Y) that occur more than once:", nrow(dk), "| rows left out:", sum(dk$N), "\n")
r <- r[!dk, on = .(Core.ID, Xr, Yr)]

b <- fread(file.path(data_dir, "Merged_BIN_for_spatial.txt"))
b[, `:=`(Xr = round(Cell.X.Position.um, 1), Yr = round(Cell.Y.Position.um, 1))]
m <- r[, .(Core.ID, Xr, Yr, NucArea, NucCompact, NucMinor, NucMajor, NucAxisRatio, NucDAPI, CellArea)][
  b, on = .(Core.ID, Xr, Yr)]
cat("analysis cells:", nrow(b), "| with nucleus data:", sum(!is.na(m$NucArea)),
    sprintf("(%.2f%%)", 100 * mean(!is.na(m$NucArea))), "\n")
m[, `:=`(Xr = NULL, Yr = NULL)]
fwrite(m, file.path(output_dir, "cells_with_nuclei.tsv"), sep = "\t")
cat("written:", file.path(output_dir, "cells_with_nuclei.tsv"), "\n")
