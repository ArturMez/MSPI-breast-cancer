#!/usr/bin/env python3
"""
01 - Discovery cohort: M1- and M2-to-T cell affinity scores

T cells: CD8+ cells, and CD4+ cells with a lymphocyte-like nucleus (p_lymph > 0.5 from script 00b).
Macrophage-marker cells: CD68 (cytoplasm mean) > 5 or CD163+. Reference pool: macrophage-marker cells that
are not T cells, so CD4+ cells with a macrophage-like nucleus stay in the pool. M1 = CD68 > 5 and CD163-,
M2 = CD68 > 5 and CD163+, both within the pool.

Algorithm (as in version 1 of this repository):
  * per TMA core with > 1 T cell, >= 1 pool cell and >= 1 M1 (or M2) cell: the pool is randomly subsampled
    to the M1 (or M2) count, 5 times; for every T cell the nearest-neighbour distance (NND) to the
    subsample and to the M1 (or M2) cells is recorded; when the nearest cell lies closer than 1 um, the
    second nearest is used
  * per patient: NNDs of 25 um or more are set to missing, separately for pool and subset, and
    delta = mean(subset NND) / (mean(pool NND) + mean(subset NND)); 0.5 means no preference
  * the random generator of each core is seeded from the cell type and core name, so reruns give identical
    scores in any order

Input : ./output/discovery_nucleus/cells_for_affinity_nucleus.tsv (script 00b)
Output: ./output/discovery_affinity_nucleus/
          discovery_distances_mean_M1_to_CD3_nucleus.txt, discovery_distances_mean_M2_to_CD3_nucleus.txt
              per-patient scores (AnnUnit, delta_CD3), read by scripts 04 to 14 (copy them to ./data/discovery/)
          affinity_per_case.tsv     per patient and cell type: delta, mean NNDs, cell and core counts
          cell_counts.tsv           T cells, pool cells and M1/M2 cells in the pool
Usage : python3 python/01_discovery_affinity_nucleus.py
Needs : Python 3.8+, numpy, scipy, pandas
"""
import os
import zlib

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

IN_FILE = "./output/discovery_nucleus/cells_for_affinity_nucleus.tsv"
OUT_DIR = "./output/discovery_affinity_nucleus"
CAP, MIN_RADIUS, NREP = 25.0, 1.0, 5


def nn_dist(points, queries):
    """Nearest-neighbour distance with the 1-um rule (k = 2)."""
    tree = cKDTree(points)
    if len(points) >= 2:
        d, _ = tree.query(queries, k=2)
        return np.where(d[:, 0] < MIN_RADIUS, d[:, 1], d[:, 0])
    d, _ = tree.query(queries, k=1)
    return np.where(d < MIN_RADIUS, np.inf, d)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    cells = pd.read_csv(IN_FILE, sep="\t", dtype={"Sample.Name": str, "Core.ID": str})
    print("cells:", len(cells), flush=True)
    xy = cells[["X", "Y"]].to_numpy(dtype=float)
    cd4 = cells["CD4"].to_numpy() == 1
    cd8 = cells["CD8"].to_numpy() == 1
    cd68 = cells["CD68"].to_numpy() > 5
    cd163 = cells["CD163"].to_numpy() == 1
    myel = cd68 | cd163
    subset_masks = {"M1": cd68 & ~cd163, "M2": cd68 & cd163}
    lymph_like = np.nan_to_num(cells["p_lymph"].to_numpy(dtype=float), nan=0.0) > 0.5
    tmask = cd8 | (cd4 & lymph_like)
    pmask = myel & ~tmask
    counts = pd.DataFrame([{"T_cells": int(tmask.sum()), "pool_cells": int(pmask.sum()),
                            "M1_in_pool": int((pmask & subset_masks["M1"]).sum()),
                            "M2_in_pool": int((pmask & subset_masks["M2"]).sum())}])
    print(counts.to_string(index=False), flush=True)
    counts.to_csv(f"{OUT_DIR}/cell_counts.tsv", sep="\t", index=False)

    core_codes, core_names = pd.factorize(cells["Core.ID"])
    case_of_core = cells.groupby(core_codes)["Sample.Name"].first().to_numpy()
    order = np.argsort(core_codes, kind="stable")
    groups = np.split(order, np.flatnonzero(np.diff(core_codes[order])) + 1)

    rows = []
    for ct, smask in subset_masks.items():
        acc = {}
        for g in groups:
            core = core_names[core_codes[g[0]]]
            t_idx = g[tmask[g]]
            p_idx = g[pmask[g]]
            s_idx = p_idx[smask[p_idx]]
            if len(t_idx) <= 1 or len(p_idx) == 0 or len(s_idx) == 0:
                continue
            txy, pxy, sxy = xy[t_idx], xy[p_idx], xy[s_idx]
            rng = np.random.default_rng(zlib.crc32(f"NUCLEUS|{ct}|{core}".encode()))
            pool_d = np.concatenate([nn_dist(pxy[rng.choice(len(pxy), size=len(sxy), replace=False)], txy)
                                     for _ in range(NREP)])
            sub_d = nn_dist(sxy, txy)
            pool_c = pool_d[pool_d < CAP]
            sub_c = sub_d[sub_d < CAP]
            a = acc.setdefault(case_of_core[core_codes[g[0]]], [0.0, 0, 0.0, 0, 0, 0, 0, 0])
            a[0] += pool_c.sum(); a[1] += len(pool_c)
            a[2] += sub_c.sum() * NREP; a[3] += len(sub_c) * NREP
            a[4] += len(t_idx); a[5] += len(s_idx); a[6] += len(p_idx); a[7] += 1
        for case, a in acc.items():
            mp = a[0] / a[1] if a[1] else np.nan
            ms = a[2] / a[3] if a[3] else np.nan
            delta = ms / (mp + ms) if np.isfinite(mp) and np.isfinite(ms) and (mp + ms) > 0 else np.nan
            rows.append((case, ct, delta, ms, mp, a[4], a[5], a[6], a[7]))
        print(f"{ct}: patients {len(acc)}", flush=True)

    per_case = pd.DataFrame(rows, columns=["Sample.Name", "celltype", "delta", "mean_subset_nnd", "mean_pool_nnd",
                                           "n_T", "n_subset", "n_pool", "n_cores"])
    per_case.to_csv(f"{OUT_DIR}/affinity_per_case.tsv", sep="\t", index=False)
    for ct in subset_masks:
        s = per_case[(per_case.celltype == ct) & per_case.delta.notna()]
        out = pd.DataFrame({"AnnUnit": s["Sample.Name"].astype(int), "delta_CD3": s["delta"]}).sort_values("AnnUnit")
        f = f"{OUT_DIR}/discovery_distances_mean_{ct}_to_CD3_nucleus.txt"
        out.to_csv(f, sep="\t", index=False)
        print("written:", f, len(out), "patients", flush=True)


if __name__ == "__main__":
    main()
