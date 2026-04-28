#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Affinity computation — NK variant (fast, parallel, multi-core per case)
------------------------------------------------------------------------
Speedups:
  • Vectorized 1-NN distance (NumPy, chunked) — no sklearn overhead
  • Parallel per-core with joblib (loky), reproducible RNG per core
  • float64 by default for parity/stability

Math/logic:
  • Per-core 1-NN to POS and to size-matched NEG subsamples (N repeats)
  • Distance cap at >= 50 -> NaN (epsilon-guarded)
  • Region mapping (default): aggregate Tumor+Stroma per core (sum area),
    keep core if region ∈ {T, Peri, Center} AND total_area > threshold
  • Aggregate to per-case delta at the end (rows-level mean by default,
    or core-level weighted averaging)
  • Output wide table: Case × (analyzed_to_centered)

CLI toggles:
  • --region-mode agg|best|first
  • --neg-replace auto|never|always
  • --per-case-agg rows|core and --core-weight equal|area|centered
  • --dtype float32|float64  (default float64)
"""

# ---- limit BLAS threads (works in Spyder) ----
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"          # if you ever switch to MKL
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"   # Apple Accelerate
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from threadpoolctl import threadpool_limits
threadpool_limits(1)   # enforce 1 inner thread for BLAS in this process

import argparse, math, sys, warnings, hashlib
from pathlib import Path
from typing import List, Tuple, Optional

from itertools import product
from tqdm import tqdm
from contextlib import contextmanager

import numpy as np
import pandas as pd
from itertools import product
# Optional parallelism

import joblib
from joblib.parallel import BatchCompletionCallBack

try:
    from joblib import Parallel, delayed
except Exception:
    Parallel = None
    delayed = None


try:
    from scipy.spatial import cKDTree
    HAVE_CKDTREE = True
except Exception:
    HAVE_CKDTREE = False
    


    
# --------------------------
# Config
# --------------------------
#DEFAULT_MARKERS: List[str] = [
#    "CD4_Single", "CD4_Treg", "CD8_Single", "CD8_Treg", "B_cells", "Tumor_cells"
#]


DEFAULT_MARKERS: List[str] = ["CD3e", "M1", "M2", "Myeloid"]

ALLOWED_REGIONS = {"T"}
AREA_COL = "Region.Area..mm2."
REGION_COL = "sample_region"
CORE_COL_CELLS = "Sample Name"   # cells table key
CORE_COL_TISSUE = "Sample.Name"  # tissue table key (R-style export)
CASE_ID_COL = "ID"

X_COL = "Cell X Position"
Y_COL = "Cell Y Position"
CK_COL = "CK"
TISSUE_CAT_COL = "Tissue Category"

DIST_CAP = 50.0
CAP_EPS = 1e-7         # epsilon so 32/64 or tiny rounding land the same side of the cap
NEG_RESAMPLES = 5
MIN_CENTERED = 2
RANDOM_SEED = 42




TARGET_COHORT = "Breast"
TARGET_REGION = "T"



# --------------------------
# Helpers
# --------------------------
# ---- top-level worker wrapper (must be at module scope for spawn) ----
def _worker_pack_fn(args):
    (core_id, cent_xy, pool_xy, pos_mask, neg_mask, seed,
     neg_resamples, dist_cap, cap_eps, chunk, fp_dtype, neg_replace_policy) = args
    dtype_np = np.float64 if fp_dtype == "float64" else np.float32
    return _compute_core_block(
        core_id=core_id,
        cent_core_xy=cent_xy,
        pool_xy=pool_xy,
        pos_mask=pos_mask,
        neg_mask=neg_mask,
        neg_resamples=neg_resamples,
        dist_cap=dist_cap,
        cap_eps=cap_eps,
        seed=seed,
        chunk=chunk,
        dtype=dtype_np,
        neg_replace_policy=neg_replace_policy,
    )
# ---- joblib <-> tqdm bridge ----
from contextlib import contextmanager

@contextmanager
def tqdm_joblib(tqdm_object):
    """
    Patch joblib to report into the given tqdm progress bar.
    """
    from joblib.parallel import BatchCompletionCallBack
    class TqdmBatchCompletionCallback(BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            # update by number of completed tasks
            tqdm_object.update(n=self.parallel.n_completed_tasks - tqdm_object.n)
            return super().__call__(*args, **kwargs)

    # >>> patch on the module, not the class <<<
    old_cb = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_cb
        tqdm_object.close()
        
        
# --------------------------
# Helpers
# --------------------------


def debug_df(df: pd.DataFrame, name: str, n: int = 3) -> None:
    print("\n" + "=" * 80, file=sys.stderr)
    print(f"[DEBUG] {name}: shape={df.shape}", file=sys.stderr)

    # columns + dtypes
    print("[DEBUG] dtypes:", file=sys.stderr)
    print(df.dtypes.sort_index(), file=sys.stderr)

    # non-null counts (like R str shows completeness)
    nn = df.notna().sum().sort_values(ascending=False)
    print("[DEBUG] non-null counts (top 25):", file=sys.stderr)
    print(nn.head(25), file=sys.stderr)

    # show a few rows (transposed so wide tables are readable)
    print(f"[DEBUG] head({n}) (transposed):", file=sys.stderr)
    print(df.head(n).T, file=sys.stderr)

    # small value preview for object columns
    obj_cols = [c for c in df.columns if df[c].dtype == "object"]
    if obj_cols:
        print("[DEBUG] example unique values (object cols, up to 5 each):", file=sys.stderr)
        for c in obj_cols[:20]:
            vals = df[c].dropna().astype(str).unique()[:5]
            print(f"  - {c}: {vals}", file=sys.stderr)

    print("=" * 80 + "\n", file=sys.stderr)
    
    
    
def _ensure_binary(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce").fillna(0.0)
    return (s != 0).astype(int)

def _canon_cols(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generic canonicalization: strip NBSP and whitespace from column names.
    DO NOT rename 'Sample Name' here (cells use that exact header).
    """
    df = df.copy()
    df.columns = [str(c).replace("\u00A0", " ").strip() for c in df.columns]
    return df

def _canon_tissue_cols(df: pd.DataFrame) -> pd.DataFrame:
    """
    Tissue-specific canonicalization: normalize to have CORE_COL_TISSUE ('Sample.Name').
    """
    df = _canon_cols(df)
    if "Sample Name" in df.columns and CORE_COL_TISSUE not in df.columns:
        df = df.rename(columns={"Sample Name": CORE_COL_TISSUE})
    return df



def add_marker_flags(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build binary immunophenotype flags for NK/NKT/M1/M2/Myeloid from base markers.
    Required base columns: CD3, CD56, NKp46, CD68, CD163, CK, Tissue Category

    Also computes:
      - All_NK_groups: union of derived immune populations
      - Tumor_cells: (Tissue Category == "Tumor") OR (CK == 1)
    """
    base_cols = ["CD3", "CD56", "NKp46", "CD68", "CD163", CK_COL, TISSUE_CAT_COL]
    missing = [c for c in base_cols if c not in df.columns]
    if missing:
        raise KeyError(f"Expected columns {missing} in the cells table.")

    # binarize base channels
    for col in ["CD3", "CD56", "NKp46", "CD68", "CD163", CK_COL]:
        df[col] = _ensure_binary(df[col])

    # phenotype logic (from OLD guide; NK uses the effective/overwriting definition)
    df["CD3e"]    = ((df["CD3"] == 1) & (df["CD56"] != 1) & (df["NKp46"] != 1) & (df["CD68"] != 1) & (df["CD163"] != 1) ).astype(int)
    df["NK"]      = ((df["CD3"] == 0) & (df["CD56"] == 1) & (df["NKp46"] == 1)).astype(int)
    df["NKT"]     = ((df["CD3"] == 1) & (df["CD56"] == 1) & (df["NKp46"] == 1)).astype(int)
    df["M1"]      = ((df["CD68"] == 1) & (df["CD163"] == 0) & (df["CD3"] != 1)).astype(int)
    df["M2"]      = ((df["CD68"] == 1) & (df["CD163"] == 1) & (df["CD3"] != 1)).astype(int)
    df["Myeloid"] = ((df["CD68"] == 0) & (df["CD163"] == 1) & (df["CD3"] != 1)).astype(int)

    # convenience groups (kept for downstream compatibility)
    df["All_NK_groups"] = ((df["CD3e"] == 1) | (df["NK"] == 1) | (df["NKT"] == 1) | (df["M1"] == 1) | (df["M2"] == 1) | (df["Myeloid"] == 1)).astype(int)
    df["CK_cells"]       = ((df[CK_COL] == 1) | (df["All_NK_groups"] == 0)).astype(int)

    # Tumor center definition kept as in your pipeline
    df["Tumor_cells"] = ((df[TISSUE_CAT_COL] == "Tumor") | (df[CK_COL] == 1)).astype(int)

    return df



# ---------- Region mappers ----------
def map_regions_best(cells: pd.DataFrame, tissue: pd.DataFrame, area_threshold: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Choose ONE region row per core:
      priority = 1 if region in ALLOWED_REGIONS and area > threshold, else 0;
      pick the row with highest (priority, area).
    Merge that region to cells (m:1). If cells already had sample_region,
    coalesce (prefer the mapped region), then keep only allowed regions.
    """
    def _canon(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df.columns = [str(c).replace("\u00A0", " ").strip() for c in df.columns]
        return df

    cells  = _canon(cells)
    tissue = _canon(tissue)
    if "Sample Name" in tissue.columns and CORE_COL_TISSUE not in tissue.columns:
        tissue = tissue.rename(columns={"Sample Name": CORE_COL_TISSUE})

    for col in [AREA_COL, REGION_COL, CORE_COL_TISSUE, CASE_ID_COL]:
        if col not in tissue.columns:
            raise KeyError(f"Expected column '{col}' in the tissue table. Found: {list(tissue.columns)}")
    if CORE_COL_CELLS not in cells.columns:
        raise KeyError(f"Expected column '{CORE_COL_CELLS}' in the cells table. Found: {list(cells.columns)}")

    tw = tissue.copy()
    tw[REGION_COL] = np.where(tw[AREA_COL] > area_threshold, tw[REGION_COL], "SMALL")
    tw = tw.assign(_prio = tw[REGION_COL].isin(ALLOWED_REGIONS).astype(int))
    tw_sorted = tw.sort_values([CORE_COL_TISSUE, "_prio", AREA_COL], ascending=[True, False, False])
    tissue_core = tw_sorted.drop_duplicates(subset=[CORE_COL_TISSUE], keep="first").drop(columns=["_prio"])

    merged = cells.merge(
        tissue_core[[CORE_COL_TISSUE, REGION_COL]],
        left_on=CORE_COL_CELLS,
        right_on=CORE_COL_TISSUE,
        how="left",
        validate="m:1"
    )

    # Coalesce sample_region_x/y if present
    if f"{REGION_COL}_y" in merged.columns or f"{REGION_COL}_x" in merged.columns:
        right = merged.get(f"{REGION_COL}_y")
        left  = merged.get(f"{REGION_COL}_x")
        if right is not None and left is not None:
            merged[REGION_COL] = right.combine_first(left)
        elif right is not None:
            merged[REGION_COL] = right
        else:
            merged[REGION_COL] = left
        for c in (f"{REGION_COL}_x", f"{REGION_COL}_y"):
            if c in merged.columns:
                merged.drop(columns=c, inplace=True)

    tissue_core_allowed = tissue_core[tissue_core[REGION_COL].isin(ALLOWED_REGIONS)].copy()
    cells_f = merged[merged[REGION_COL].isin(ALLOWED_REGIONS)].copy()
    return cells_f, tissue_core_allowed

def map_regions_r_like(cells: pd.DataFrame, tissue: pd.DataFrame, area_threshold: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    R-like 'first row wins' after tagging SMALL by row (area <= threshold).
    Keeps all cores per case (no dedup by case).
    """
    t = _canon_tissue_cols(tissue)
    c = _canon_cols(cells)

    t[REGION_COL] = np.where(t[AREA_COL] > area_threshold, t[REGION_COL], "SMALL")
    t_first = t.drop_duplicates(subset=[CORE_COL_TISSUE], keep="first")

    region_map = t_first.set_index(CORE_COL_TISSUE)[REGION_COL]
    c[REGION_COL] = c[CORE_COL_CELLS].map(region_map)
    cells_f = c[c[REGION_COL].isin(ALLOWED_REGIONS)].copy()

    tissue_core_allowed = t_first[t_first[REGION_COL].isin(ALLOWED_REGIONS)][[CORE_COL_TISSUE, CASE_ID_COL]].copy()
    return cells_f, tissue_core_allowed

def map_regions_aggregate(cells: pd.DataFrame, tissue: pd.DataFrame, area_threshold: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Per-core aggregation over Tumor+Stroma rows:
      • total_area = sum(Region.Area..mm2.)
      • region_for_core = mode(first non-null) of sample_region
      • If total_area <= area_threshold OR region_for_core ∉ allowed -> SMALL
      • Keep only cores with region ∈ ALLOWED_REGIONS
    Returns:
      cells_f: cells with allowed cores (m:1 mapped region)
      tissue_core_allowed: one row per core with region, Case ID, and total_area
    """
    # Canonicalize (cells keep 'Sample Name'; tissue normalized to 'Sample.Name')
    c = _canon_cols(cells)
    t = _canon_tissue_cols(tissue)

    for col in [AREA_COL, REGION_COL, CORE_COL_TISSUE, CASE_ID_COL]:
        if col not in t.columns:
            raise KeyError(f"Expected column '{col}' in tissue table, found: {list(t.columns)}")

    # Per-core aggregate
    def _mode_non_na(s: pd.Series):
        s = s.dropna()
        return s.mode().iat[0] if not s.empty else np.nan

    gb = (t.groupby(CORE_COL_TISSUE, as_index=False)
            .agg(total_area=(AREA_COL, "sum"),
                 region_mode=(REGION_COL, _mode_non_na),
                 ID=(CASE_ID_COL, lambda s: s.dropna().iat[0] if not s.dropna().empty else np.nan)))

    gb[REGION_COL] = np.where(
        (gb["total_area"] > area_threshold) & (gb["region_mode"].isin(ALLOWED_REGIONS)),
        gb["region_mode"],
        "SMALL"
    )
    tissue_core = gb[[CORE_COL_TISSUE, REGION_COL, "total_area", "ID"]]

    # Map region to cells (m:1)
    merged = c.merge(
        tissue_core[[CORE_COL_TISSUE, REGION_COL]],
        left_on=CORE_COL_CELLS, right_on=CORE_COL_TISSUE,
        how="left", validate="m:1"
    )

    # Coalesce sample_region_y/x -> sample_region (if cells already had sample_region)
    if f"{REGION_COL}_y" in merged.columns or f"{REGION_COL}_x" in merged.columns:
        right = merged.get(f"{REGION_COL}_y")
        left  = merged.get(f"{REGION_COL}_x")
        if right is not None and left is not None:
            merged[REGION_COL] = right.combine_first(left)
        elif right is not None:
            merged[REGION_COL] = right
        else:
            merged[REGION_COL] = left
        for ccol in (f"{REGION_COL}_x", f"{REGION_COL}_y"):
            if ccol in merged.columns:
                merged.drop(columns=ccol, inplace=True)
    elif REGION_COL not in merged.columns:
        raise KeyError(f"After merge, '{REGION_COL}' not found. "
                       f"Cells columns: {list(c.columns)}; Tissue-core: {list(tissue_core.columns)}")

    # Keep only allowed regions
    cells_f = merged[merged[REGION_COL].isin(ALLOWED_REGIONS)].copy()
    tissue_core_allowed = tissue_core[tissue_core[REGION_COL].isin(ALLOWED_REGIONS)].copy()
    return cells_f, tissue_core_allowed

# ---------- Stable per-core seeds ----------
def stable_core_seed(global_seed: int, core_id: str) -> int:
    h = hashlib.blake2b(f"{global_seed}:{core_id}".encode(), digest_size=8).digest()
    return int.from_bytes(h, "little") & 0x7fffffff

# ---------- Fast 1-NN (vectorized, chunked) ----------
def _min_dist_tiled(A: np.ndarray,
                    B: np.ndarray,
                    a_chunk: int = 4096,
                    b_chunk: int = 2048,
                    dtype=np.float64) -> np.ndarray:
    """
    Min distance from each row in A to any row in B, using 2-D tiles.
    Uses ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a·b to avoid 3-D tensors.
    """
    if A.size == 0 or B.size == 0:
        return np.array([], dtype=dtype)
    A = A.astype(dtype, copy=False)
    B = B.astype(dtype, copy=False)

    nA = A.shape[0]
    out = np.full(nA, np.inf, dtype=dtype)

    for ai in range(0, nA, a_chunk):
        a = A[ai:ai + a_chunk]                        # (ta, 2)
        a2 = (a * a).sum(axis=1, keepdims=True)       # (ta, 1)
        best = np.full(a.shape[0], np.inf, dtype=dtype)

        for bj in range(0, B.shape[0], b_chunk):
            b = B[bj:bj + b_chunk]                    # (tb, 2)
            # 2-D GEMM instead of 3-D broadcast
            dot = a @ b.T                             # (ta, tb)
            b2 = (b * b).sum(axis=1)[None, :]         # (1, tb)
            d2 = a2 + b2 - 2.0 * dot                  # (ta, tb)
            # numerical guard: tiny negatives -> 0
            np.maximum(d2, 0, out=d2)
            local_min = d2.min(axis=1)                # (ta,)
            np.maximum(local_min, 0, out=local_min)   # guard
            np.sqrt(local_min, out=local_min)         # in-place sqrt
            best = np.minimum(best, local_min)

        out[ai:ai + a.shape[0]] = best

    return out

# ---------- Core worker (called in parallel) ----------


def _compute_core_block(core_id: str,
                        cent_core_xy: np.ndarray,
                        pool_xy: np.ndarray,
                        pos_mask: np.ndarray,
                        neg_mask: np.ndarray,
                        neg_resamples: int,
                        dist_cap: float,
                        cap_eps: float,
                        seed: int,
                        chunk: int,
                        dtype,
                        neg_replace_policy: str = "auto") -> pd.DataFrame:
    """
    Compute distances for one core. pool_xy/pos_mask/neg_mask refer to the same row order.
    neg_replace_policy: "auto" (default), "never", or "always"
    """
    rng = np.random.default_rng(seed)

    if cent_core_xy.shape[0] < MIN_CENTERED:
        return pd.DataFrame({"core_ID":[core_id], "REF_to_central":[np.nan], "CELL_to_central":[np.nan]})

    # POS = analyzed marker positives
    pos_xy = pool_xy[pos_mask]

    # NEG reference pool = restricted subset (M1|M2|Myeloid)
    neg_xy = pool_xy[neg_mask]

    if pos_xy.shape[0] == 0 or neg_xy.shape[0] == 0:
        return pd.DataFrame({"core_ID":[core_id], "REF_to_central":[np.nan], "CELL_to_central":[np.nan]})

    # POS distances once
    if HAVE_CKDTREE and np.isfinite(dist_cap):
        ub = max(dist_cap - cap_eps, 1e-12)
        tree_pos = cKDTree(pos_xy)
        pos_d, _ = tree_pos.query(cent_core_xy, k=1, distance_upper_bound=ub)
        pos_d = pos_d.astype(dtype, copy=False)
        pos_d[~np.isfinite(pos_d) | (pos_d >= ub)] = np.nan
    else:
        pos_d = _min_dist_tiled(cent_core_xy, pos_xy, a_chunk=chunk, b_chunk=2048, dtype=dtype)

    # NEG distances repeated with size-matched subsamples from neg_xy
    neg_d_list = []
    n_pos = pos_xy.shape[0]
    for _ in range(neg_resamples):
        if neg_replace_policy == "always":
            idx = rng.integers(0, neg_xy.shape[0], size=n_pos)
        elif neg_replace_policy == "never":
            if neg_xy.shape[0] < n_pos:
                raise ValueError(f"NEG < POS in core {core_id}: {neg_xy.shape[0]} < {n_pos}")
            idx = rng.choice(neg_xy.shape[0], size=n_pos, replace=False)
        else:  # auto
            if neg_xy.shape[0] >= n_pos:
                idx = rng.choice(neg_xy.shape[0], size=n_pos, replace=False)
            else:
                idx = rng.integers(0, neg_xy.shape[0], size=n_pos)

        sampled_neg = neg_xy[idx, :]

        if HAVE_CKDTREE and np.isfinite(dist_cap):
            ub = max(dist_cap - cap_eps, 1e-12)
            tree_neg = cKDTree(sampled_neg)
            nd, _ = tree_neg.query(cent_core_xy, k=1, distance_upper_bound=ub)
            nd = nd.astype(dtype, copy=False)
            nd[~np.isfinite(nd) | (nd >= ub)] = np.nan
        else:
            nd = _min_dist_tiled(cent_core_xy, sampled_neg, a_chunk=chunk, b_chunk=2048, dtype=dtype)

        neg_d_list.append(nd)

    neg_d_all = np.concatenate(neg_d_list, axis=0)

    # Repeat POS distances to align with concatenated NEG (R-like recycling)
    reps = math.ceil(len(neg_d_all) / len(pos_d)) if len(pos_d) > 0 else 0
    pos_rep = np.tile(pos_d, reps)[:len(neg_d_all)] if reps > 0 else np.array([], dtype=dtype)

    # Apply distance cap -> NaN (epsilon-guarded)
    if dist_cap is not None and np.isfinite(dist_cap):
        neg_d_all = np.where(neg_d_all >= (dist_cap - cap_eps), np.nan, neg_d_all)
        pos_rep   = np.where(pos_rep   >= (dist_cap - cap_eps), np.nan, pos_rep)

    return pd.DataFrame({"core_ID": core_id,
                         "REF_to_central": neg_d_all,
                         "CELL_to_central": pos_rep})




    if pos_xy.shape[0] == 0 or neg_xy.shape[0] == 0:
        return pd.DataFrame({"core_ID":[core_id], "REF_to_central":[np.nan], "CELL_to_central":[np.nan]})

    # POS distances once
    if HAVE_CKDTREE and np.isfinite(dist_cap):
        ub = max(dist_cap - cap_eps, 1e-12)
        tree_pos = cKDTree(pos_xy)
        pos_d, _ = tree_pos.query(cent_core_xy, k=1, distance_upper_bound=ub)
        pos_d = pos_d.astype(dtype, copy=False)
        pos_d[~np.isfinite(pos_d) | (pos_d >= ub)] = np.nan
    else:
        pos_d = _min_dist_tiled(cent_core_xy, pos_xy, a_chunk=chunk, b_chunk=2048, dtype=dtype)

    # NEG distances repeated with size-matched subsamples
    neg_d_list = []
    n_pos = pos_xy.shape[0]
    for _ in range(neg_resamples):
        if neg_replace_policy == "always":
            idx = rng.integers(0, neg_xy.shape[0], size=n_pos)
        elif neg_replace_policy == "never":
            if neg_xy.shape[0] < n_pos:
                raise ValueError(f"NEG < POS in core {core_id}: {neg_xy.shape[0]} < {n_pos}")
            idx = rng.choice(neg_xy.shape[0], size=n_pos, replace=False)
        else:  # auto
            if neg_xy.shape[0] >= n_pos:
                idx = rng.choice(neg_xy.shape[0], size=n_pos, replace=False)
            else:
                idx = rng.integers(0, neg_xy.shape[0], size=n_pos)

        sampled_neg = neg_xy[idx, :]

        if HAVE_CKDTREE and np.isfinite(dist_cap):
            tree_neg = cKDTree(sampled_neg)
            ub = max(dist_cap - cap_eps, 1e-12)
            nd, _ = tree_neg.query(cent_core_xy, k=1, distance_upper_bound=ub)
            nd = nd.astype(dtype, copy=False)
            nd[~np.isfinite(nd) | (nd >= ub)] = np.nan
        else:
            nd = _min_dist_tiled(cent_core_xy, sampled_neg, a_chunk=chunk, b_chunk=2048, dtype=dtype)

        neg_d_list.append(nd)

    neg_d_all = np.concatenate(neg_d_list, axis=0)

    # Repeat POS distances to align with concatenated NEG (R-like recycling)
    reps = math.ceil(len(neg_d_all) / len(pos_d)) if len(pos_d) > 0 else 0
    pos_rep = np.tile(pos_d, reps)[:len(neg_d_all)] if reps > 0 else np.array([], dtype=dtype)

    # Apply distance cap -> NaN (epsilon-guarded)
    if dist_cap is not None and np.isfinite(dist_cap):
        neg_d_all = np.where(neg_d_all >= (dist_cap - cap_eps), np.nan, neg_d_all)
        pos_rep   = np.where(pos_rep   >= (dist_cap - cap_eps), np.nan, pos_rep)

    return pd.DataFrame({"core_ID": core_id,
                         "REF_to_central": neg_d_all,
                         "CELL_to_central": pos_rep})





def compute_affinity_for_pair(
    cells: pd.DataFrame,
    tissue_core_allowed: pd.DataFrame,
    centered_marker: str,
    analyzed_marker: str,
    neg_resamples: int,
    dist_cap: float,
    cap_eps: float,
    random_seed: int,
    n_jobs: int,
    chunk_size: int,
    fp_dtype: str,
    neg_replace_policy: str,
    per_case_agg: str,
    core_weight: str
) -> pd.DataFrame:
    """Compute per-case delta for one pair (centered -> analyzed)."""
    if centered_marker not in cells.columns or analyzed_marker not in cells.columns:
        raise KeyError(f"Marker not found: {centered_marker} or {analyzed_marker}")

    cols_keep = [CORE_COL_CELLS, X_COL, Y_COL, centered_marker, analyzed_marker, "M1", "M2", "Myeloid"]
    cols_keep = list(dict.fromkeys(cols_keep))

    df = cells[cols_keep].copy()
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]

    cent = df[df[centered_marker] == 1].copy()
    pool = df[df[centered_marker] != 1].copy()
    if cent.empty:
        return pd.DataFrame(columns=["Case", "CELL", "delta_CELL"])

    core_ids = cent[CORE_COL_CELLS].unique().tolist()
    dtype_np = np.float64 if fp_dtype == "float64" else np.float32

    seed_by_core = {cid: stable_core_seed(random_seed, cid) for cid in core_ids}

    # Build compact packs per core (NumPy arrays only)
    packs = []
    for cid in core_ids:
        cent_core = cent.loc[cent[CORE_COL_CELLS] == cid, [X_COL, Y_COL]].to_numpy(dtype=dtype_np, copy=False)
        cols = [X_COL, Y_COL, analyzed_marker, "M1", "M2", "Myeloid"]
        cols = list(dict.fromkeys(cols))  # de-duplicate (critical when analyzed_marker is M1/M2/Myeloid)
        core_pool = pool.loc[pool[CORE_COL_CELLS] == cid, cols]
        
        if cent_core.shape[0] == 0 or core_pool.shape[0] == 0:
            continue
        
        pool_xy = core_pool[[X_COL, Y_COL]].to_numpy(dtype=dtype_np, copy=False)
        
        # force analyzed column to be a 1D vector even if pandas returns a DataFrame
        am = core_pool[analyzed_marker]
        if isinstance(am, pd.DataFrame):
            am = am.iloc[:, 0]
        pos_mask = (am.to_numpy() == 1)
        
        neg_mask = (core_pool[["M1", "M2", "Myeloid"]].to_numpy().any(axis=1))
        if not neg_mask.any():
            continue

        packs.append((
            cid, cent_core, pool_xy, pos_mask, neg_mask, seed_by_core[cid],
            neg_resamples, dist_cap, cap_eps, chunk_size, fp_dtype, neg_replace_policy
        ))

    if not packs:
        return pd.DataFrame(columns=["Case", "CELL", "delta_CELL"])

    # ---- Execute in parallel (prefer processes; robust in Spyder) ----
    parts = None
    if n_jobs != 1:
        try:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor
            workers = mp.cpu_count() if n_jobs == -1 else max(1, int(n_jobs))
            ctx = mp.get_context("spawn")
            # Reasonable granularity so all workers stay busy
            chunk = max(1, len(packs) // (workers * 4))
            with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
                parts = list(ex.map(_worker_pack_fn, packs, chunksize=chunk))
        except Exception:
            # Fallback to joblib if the executor is blocked by the host env
            if Parallel is not None:
                parts = Parallel(n_jobs=n_jobs, backend="loky", max_nbytes=None)(
                    delayed(_worker_pack_fn)(p) for p in packs
                )

    if parts is None:
        # Last-resort sequential
        parts = [_worker_pack_fn(p) for p in packs]

    parts = [p for p in parts if p is not None and not p.empty]
    if not parts:
        return pd.DataFrame(columns=["Case", "CELL", "delta_CELL"])

    nn_all = pd.concat(parts, ignore_index=True)

    # Map to Case / area
    mapping = tissue_core_allowed[[CORE_COL_TISSUE, CASE_ID_COL]].drop_duplicates().rename(
        columns={CORE_COL_TISSUE: "core_ID"}
    )
    nn_all = nn_all.merge(mapping, on="core_ID", how="left")

    if "total_area" in tissue_core_allowed.columns:
        map_area = tissue_core_allowed[[CORE_COL_TISSUE, "total_area"]].drop_duplicates().rename(
            columns={CORE_COL_TISSUE: "core_ID"}
        )
        nn_all = nn_all.merge(map_area, on="core_ID", how="left")

    # Aggregate to per-case delta
    if per_case_agg == "rows":
        g = nn_all.groupby(CASE_ID_COL, dropna=True)[["REF_to_central","CELL_to_central"]].mean()
        denom = g["REF_to_central"] + g["CELL_to_central"]
        delta = g["CELL_to_central"] / denom
        per_case = pd.DataFrame({"Case": g.index, "delta_CELL": delta.values})
    else:
        core_means = (nn_all.groupby(["core_ID", CASE_ID_COL], dropna=True)
                            [["REF_to_central","CELL_to_central"]]
                            .mean().reset_index())
        if core_weight == "area" and "total_area" in nn_all.columns:
            w_map = nn_all[["core_ID","total_area"]].drop_duplicates()
            core_means = core_means.merge(w_map, on="core_ID", how="left")
            w = core_means["total_area"].fillna(0.0).to_numpy()
        elif core_weight == "centered":
            cnt = (nn_all.groupby(["core_ID", CASE_ID_COL])["CELL_to_central"]
                        .apply(lambda s: s.notna().sum())
                        .reset_index(name="n_centered_rows"))
            core_means = core_means.merge(cnt, on=["core_ID", CASE_ID_COL], how="left")
            w = core_means["n_centered_rows"].fillna(0.0).to_numpy()
        else:
            w = np.ones(len(core_means), dtype=float)

        def wmean(df: pd.DataFrame, col: str, wcol: str) -> float:
            v = df[col].to_numpy(); wv = df[wcol].to_numpy()
            m = np.isfinite(v) & np.isfinite(wv) & (wv > 0)
            return np.sum(v[m] * wv[m]) / np.sum(wv[m]) if m.any() else np.nan

        core_means["w"] = w
        rows = []
        for case, grp in core_means.groupby(CASE_ID_COL):
            ref_m  = wmean(grp, "REF_to_central",  "w")
            cell_m = wmean(grp, "CELL_to_central", "w")
            denom = ref_m + cell_m
            delta = (cell_m / denom) if np.isfinite(denom) and denom != 0 else np.nan
            rows.append((case, delta))
        per_case = pd.DataFrame(rows, columns=["Case", "delta_CELL"])

    per_case["CELL"] = f"{analyzed_marker}_to_{centered_marker}"
    return per_case[["Case", "CELL", "delta_CELL"]]




def run_affinity(
    cells_path: Path,
    tissue_path: Path,
    out_path: Path,
    markers: List[str],
    area_threshold: float,
    dist_cap: float,
    cap_eps: float,
    neg_resamples: int,
    random_seed: int,
    n_jobs: int,
    chunk_size: int,
    fp_dtype: str,
    region_mode: str,
    neg_replace_policy: str,
    per_case_agg: str,
    core_weight: str
):
    cells  = pd.read_csv(cells_path, sep="\t", low_memory=False)
    tissue = pd.read_csv(tissue_path, sep="\t", low_memory=False)
    
    
    debug_df(cells, "cells (raw)")
    debug_df(tissue, "tissue (raw)")



    # -------------------------
    # HARD FILTER: cohort + region
    # -------------------------
    # Keep only Breast cohort
    if "cohort" not in cells.columns or "cohort" not in tissue.columns:
        raise KeyError("Expected column 'cohort' in both cells and tissue tables.")
    cells = cells.loc[cells["cohort"] == TARGET_COHORT].copy()
    tissue = tissue.loc[tissue["cohort"] == TARGET_COHORT].copy()
    
    # Keep only sample_region == 'T'
    # Note: sample_region has some missing values, drop them implicitly by the filter
    if REGION_COL not in cells.columns or REGION_COL not in tissue.columns:
        raise KeyError(f"Expected column '{REGION_COL}' in both cells and tissue tables.")
    cells = cells.loc[cells[REGION_COL] == TARGET_REGION].copy()
    tissue = tissue.loc[tissue[REGION_COL] == TARGET_REGION].copy()
    
    # Optional debug
    debug_df(cells, f"cells (after cohort={TARGET_COHORT}, region={TARGET_REGION})")
    debug_df(tissue, f"tissue (after cohort={TARGET_COHORT}, region={TARGET_REGION})")



    for col in [X_COL, Y_COL]:
        if col not in cells.columns:
            raise KeyError(f"Expected XY column '{col}' in the cells table.")

    cells = add_marker_flags(cells)



    # --- NEW: keep only cells that belong to the marker set in use ---
    # 'markers' is the effective list (defaults to DEFAULT_MARKERS or your CLI --markers)
    keep_cols = [m for m in markers if m in cells.columns]
    if not keep_cols:
        raise ValueError(f"No marker columns found among: {markers}")
    keep_mask = cells[keep_cols].any(axis=1)   # row has at least one marker==1
    cells = cells.loc[keep_mask].copy()
    # -----------------------------------------------------------------
    
    
    if region_mode == "first":
        cells_f, tissue_core_allowed = map_regions_r_like(cells, tissue, area_threshold=area_threshold)
    elif region_mode == "best":
        cells_f, tissue_core_allowed = map_regions_best(cells, tissue, area_threshold=area_threshold)
    else:
        cells_f, tissue_core_allowed = map_regions_aggregate(cells, tissue, area_threshold=area_threshold)

    need_cols = set([CORE_COL_CELLS, X_COL, Y_COL, TISSUE_CAT_COL, CK_COL, "All_NK_groups"]) | set(markers)
    cells_f = cells_f[[c for c in need_cols if c in cells_f.columns]].copy()

    all_pairs = []
    # helpful debug


    # ------------------------------------------------------------
    # Parallelize across PAIRS (top level) -> robust in Spyder
    # ------------------------------------------------------------
    pairs = [("CD3e", "M1"), ("CD3e", "M2")]
    print(f"#cores={cells_f[CORE_COL_CELLS].nunique()}, #pairs={len(pairs)}", file=sys.stderr)

    def _run_one_pair(c_center: str, c_an: str) -> pd.DataFrame:
        # IMPORTANT: disable inner parallelism to avoid nested pools
        return compute_affinity_for_pair(
            cells=cells_f,
            tissue_core_allowed=tissue_core_allowed,
            centered_marker=c_center,
            analyzed_marker=c_an,
            neg_resamples=neg_resamples,
            dist_cap=dist_cap,
            cap_eps=cap_eps,
            random_seed=random_seed,
            n_jobs=1,                 # <- no inner pool; outer pool only
            chunk_size=chunk_size,
            fp_dtype=fp_dtype,
            neg_replace_policy=neg_replace_policy,
            per_case_agg=per_case_agg,
            core_weight=core_weight
        )

    if Parallel is not None and n_jobs != 1:
        # Prefer separate processes (loky); also show a tqdm that actually reflects jobs finishing.
        with tqdm_joblib(tqdm(total=len(pairs), desc="Pairs", dynamic_ncols=True)) as _:
            pair_results = Parallel(n_jobs=n_jobs, backend="loky", max_nbytes=None)(
                delayed(_run_one_pair)(c_center, c_an) for (c_center, c_an) in pairs
            )
    else:
        pair_results = []
        for (c_center, c_an) in tqdm(pairs, total=len(pairs), desc="Pairs", dynamic_ncols=True):
            pair_results.append(_run_one_pair(c_center, c_an))
            
    
    
    all_pairs = [p for p in pair_results if p is not None and not p.empty]

    if not all_pairs:
        warnings.warn("No pairwise results were produced. Check markers and input files.")
        pd.DataFrame(columns=["Case"]).to_csv(out_path, sep="\t", index=False)
        return

    long_df = pd.concat(all_pairs, ignore_index=True)
    wide = long_df.pivot_table(index="Case", columns="CELL", values="delta_CELL", aggfunc="mean")

    all_cases = tissue_core_allowed[CASE_ID_COL].dropna().drop_duplicates().rename("Case")
    wide = all_cases.to_frame().merge(wide.reset_index(), on="Case", how="left")

    wide.to_csv(out_path, sep="\t", index=False)
    print(f"Saved: {out_path}")

def main():
    p = argparse.ArgumentParser(description="Affinity computation (NK, fast)")
    p.add_argument("--cells",  type=Path, default=Path("./data/cells.txt"))
    p.add_argument("--tissue", type=Path, default=Path("./data/tissue.txt"))
    p.add_argument("--out",    type=Path, default=Path("./output/affinity_CD3e_M1M2_myeloidpool_50um.txt"))
    p.add_argument("--area-threshold", type=float, default=0.1)
    p.add_argument("--dist-cap", type=float, default=DIST_CAP)
    p.add_argument("--cap-eps", type=float, default=CAP_EPS, help="epsilon for cap comparison")
    p.add_argument("--neg-resamples", type=int, default=NEG_RESAMPLES)
    p.add_argument("--seed", type=int, default=RANDOM_SEED)

    # Performance knobs
    p.add_argument("--n-jobs", type=int, default=-1, help="-1 = use all CPU cores")
    p.add_argument("--chunk-size", type=int, default=8192, help="rows of centers per distance chunk")
    p.add_argument("--dtype", type=str, choices=["float32", "float64"], default="float64",
                   help="float precision for distance math (default float64)")

    # Behavior knobs
    p.add_argument("--region-mode", choices=["agg", "best", "first"], default="agg",
                   help="agg=sum Tumor+Stroma per core (default); best=pick allowed row with largest area; first=R-like first row")
    p.add_argument("--neg-replace", choices=["auto","never","always"], default="auto",
                   help="auto: no-replacement if NEG>=POS else with-replacement; never: strict no-replacement; always: always with-replacement")
    p.add_argument("--per-case-agg", choices=["rows","core"], default="rows",
                   help="rows = mean over all rows (default); core = per-core mean then weighted across cores")
    p.add_argument("--core-weight", choices=["equal","area","centered"], default="equal",
                   help="when --per-case-agg core: equal=unweighted; area=weight by total_area; centered=weight by #centered rows")

    p.add_argument("--markers", type=str, nargs="*", default=DEFAULT_MARKERS)
    args = p.parse_args()

    run_affinity(
        cells_path=args.cells,
        tissue_path=args.tissue,
        out_path=args.out,
        markers=args.markers,
        area_threshold=args.area_threshold,
        dist_cap=args.dist_cap,
        cap_eps=args.cap_eps,
        neg_resamples=args.neg_resamples,
        random_seed=args.seed,
        n_jobs=args.n_jobs,
        chunk_size=args.chunk_size,
        fp_dtype=args.dtype,
        region_mode=args.region_mode,
        neg_replace_policy=args.neg_replace,
        per_case_agg=args.per_case_agg,
        core_weight=args.core_weight
    )

if __name__ == "__main__":
    # make sure loky/multiprocessing can actually spawn workers on macOS/Spyder
    import multiprocessing as mp
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        # it's ok if already set
        pass
    main()
    
    