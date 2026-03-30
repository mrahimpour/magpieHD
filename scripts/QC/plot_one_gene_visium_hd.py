#!/usr/bin/env python3
"""
load Visium HD counts + bin positions, normalize, plot a certain gene on H&E.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


SAMPLE = "1719-MMC-PDAC"
BIN_UM = "002"  
GENE = "P2ry6"

MATRIX = f"input/{SAMPLE}/visium/{BIN_UM}um/filtered_feature_bc_matrix.h5"
POSITIONS = f"input/{SAMPLE}/visium/spatial/tissue_square_{BIN_UM}um_positions.parquet"
SPATIAL_DIR = f"input/{SAMPLE}/visium/spatial"  # for H&E + scalefactors JSON
OUT_PNG = f"output/{SAMPLE}/qc/{GENE}_spatial.png"

# 1) Count matrix: a 10x HDF5 formatted matrix
adata = sc.read_10x_h5(MATRIX)
adata.var_names_make_unique()  

# 2) Find the gene column index in the matrix 
if GENE not in adata.var_names:
    raise SystemExit(f"Gene {GENE!r} not found. Check spelling / species / matrix.")
g_idx = list(adata.var_names).index(GENE)

# 3) Bin positions (full-resolution pixels): Parquet has barcode + x + y columns
pos = pd.read_parquet(POSITIONS)
bc_col = "barcode" if "barcode" in pos.columns else "bin_id"
x_col = "pxl_col_in_fullres" if "pxl_col_in_fullres" in pos.columns else "x"
y_col = "pxl_row_in_fullres" if "pxl_row_in_fullres" in pos.columns else "y"
pos = pos[[bc_col, x_col, y_col]].rename(columns={bc_col: "bc", x_col: "x", y_col: "y"})
pos["bc"] = pos["bc"].astype(str)

# 4) Keep only matrix rows that have a position (same barcode string as in the matrix)
bc_mat = pd.Index(adata.obs_names.astype(str))
adata = adata[bc_mat.isin(pos["bc"])].copy()

# 5) Align coordinates to the same row order as adata
pos = pos.set_index("bc").loc[adata.obs_names.astype(str)]

# 6) Standard normalization: scale each bin to same "library size", then log(1+x)
sc.pp.normalize_total(adata, target_sum=1e4)  # fills adata.X with normalized values
sc.pp.log1p(adata)  
y_plot = adata.X[:, g_idx].toarray().ravel()  # normalized + log1p value for GENE

# 7) Map fullres (x,y) to the hires PNG using Space Ranger scale factor
with open(Path(SPATIAL_DIR) / "scalefactors_json.json", encoding="utf-8") as f:
    sf = json.load(f)
scale = float(sf["tissue_hires_scalef"])
xh = pos["x"].to_numpy() * scale
yh = pos["y"].to_numpy() * scale

# 8) Draw H&E + colored points (only bins with some signal after normalization)
img = plt.imread(Path(SPATIAL_DIR) / "tissue_hires_image.png")
mask = y_plot > 0
Path(OUT_PNG).parent.mkdir(parents=True, exist_ok=True)
fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(img)
sc = ax.scatter(xh[mask], yh[mask], c=y_plot[mask], s=2, cmap="afmhot", linewidths=0)
plt.colorbar(sc, ax=ax, shrink=0.5, label="log1p(norm counts)")
ax.set_title(f"{GENE}  (n bins with signal: {mask.sum()})")
ax.axis("off")
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
plt.close()
print("Saved:", Path(OUT_PNG).resolve())
