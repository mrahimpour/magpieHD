#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collapse transformed MSI to Visium HD bins and write a Space Ranger–style MEX folder.

Note: Each transformed MSI point is assigned to the nearest HD bin center (deterministic
nearest-neighbor in fullres pixels). That is exact under that rule; it does not split
signal across bins or use polygon overlap. For radius-based aggregation or joint gene+MSI
workflows (e.g. hd_bins_per_msi_spot.csv), see scripts/aggregate_to_msi_resolution.py.

Inputs:
  --sample: Sample ID
  --bin_size: HD bin resolution (002, 008, or 016)
  --msi_meta: MSI metadata with x,y and spot_id
  --msi_int: MSI intensities with mz columns
  --transformed_csv: output/<sample>/transformed.csv
  --visium_hd_bins: output/<sample>/visium_hd_square_{bin_size}um_bins.csv

Outputs:
  output/<sample>/{bin_size}um/spaceranger_hd/
    ├── filtered_feature_bc_matrix/
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    └── spatial/
"""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree


# ------------------------------------------------------------
# Geometry: assign each MSI point to nearest HD bin (fullres space)
# ------------------------------------------------------------
def nearest_bins(points_xy, bins_xy):
    """
    Assign each MSI point to nearest bin center.
    points_xy: (N,2), bins_xy: (M,2)
    returns: (distances, indices) - distances (N,) and indices (N,) of nearest bin for each point
    """
    tree = cKDTree(bins_xy)
    distances, idx = tree.query(points_xy, k=1)
    return distances, idx


# ------------------------------------------------------------
# Write Space Ranger–style MEX (matrix + features + barcodes)
# ------------------------------------------------------------
def write_mex(out_dir: Path, X_csr: sparse.csr_matrix, barcodes, features):
    mex = out_dir / "filtered_feature_bc_matrix"
    mex.mkdir(parents=True, exist_ok=True)
    # matrix.mtx.gz
    with gzip.open(mex/"matrix.mtx.gz", "wt") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write("%\n")
        f.write(f"{X_csr.shape[0]} {X_csr.shape[1]} {X_csr.nnz}\n")
        coo = X_csr.tocoo()
        for r, c, v in zip(coo.row, coo.col, coo.data):
            f.write(f"{r+1} {c+1} {float(v)}\n")
    # features.tsv.gz  (gene_id  gene_name  feature_type)
    feats = pd.DataFrame({"gene_id": features,
                          "gene_name": features,
                          "feature_type": "MSI"})
    with gzip.open(mex/"features.tsv.gz", "wt") as f:
        feats.to_csv(f, sep="\t", header=False, index=False)
    # barcodes.tsv.gz
    with gzip.open(mex/"barcodes.tsv.gz", "wt") as f:
        pd.Series(barcodes).to_csv(f, index=False, header=False)

def main():
    ap = argparse.ArgumentParser(description="Collapse MSI to HD bins and write Space Ranger-style MEX")
    ap.add_argument("--sample", required=True, help="Sample ID")
    ap.add_argument("--bin_size", default="008", choices=["002", "008", "016"],
                    help="HD bin resolution: 002 (2µm), 008 (8µm), or 016 (16µm)")
    ap.add_argument("--msi_meta", default=None, help="MSI metadata CSV")
    ap.add_argument("--msi_int", default=None, help="MSI intensities CSV")
    ap.add_argument("--transformed_csv", default=None, help="Transformed MSI coords")
    ap.add_argument("--visium_hd_bins", default=None, help="HD bin positions CSV")
    ap.add_argument("--spaceranger_dir", default=None, help="Original Visium HD Space Ranger dir")
    ap.add_argument("--out_root", default=None, help="Output directory (default: output/<sample>/<bin_size>um)")
    args = ap.parse_args()
    
    sample = args.sample
    bin_size = args.bin_size
    
    # Set default paths based on sample and bin_size
    msi_meta = args.msi_meta or f"input/{sample}/msi/MSI_metadata.csv"
    msi_int = args.msi_int or f"input/{sample}/msi/MSI_intensities.csv"
    transformed_csv = args.transformed_csv or f"output/{sample}/transformed.csv"
    visium_hd_bins = args.visium_hd_bins or f"output/{sample}/visium_hd_square_{bin_size}um_bins.csv"
    spaceranger_dir = args.spaceranger_dir or f"input/{sample}/visium"
    
    # Output directory: organized by resolution
    out_root = Path(args.out_root) if args.out_root else Path(f"output/{sample}/{bin_size}um")
    out_root.mkdir(parents=True, exist_ok=True)
    out_spaceranger = out_root / "spaceranger_hd"
    out_spaceranger.mkdir(exist_ok=True)
    
    print(f"[info] Sample: {sample}")
    print(f"[info] HD bin resolution: {bin_size}µm")
    print(f"[info] Output: {out_spaceranger}")
    
    # Load scale factor for coordinate conversion
    # MSI transformed.csv is in HIRES pixel space (alter_data overlays on tissue_hires_image.png).
    # HD bin CSV from create_visium_hd_index.py is in FULLRES pixels — scale hires → fullres.
    sf_path = Path(spaceranger_dir) / "spatial" / "scalefactors_json.json"
    with open(sf_path) as f:
        sf = json.load(f)
    tissue_hires_scalef = sf.get('tissue_hires_scalef', 0.5)
    print(f"[info] Coordinate scale: HIRES → FULLRES (÷ {tissue_hires_scalef:.4f})")
    
    # Load MSI data
    meta = pd.read_csv(msi_meta)
    ints = pd.read_csv(msi_int)
    assert "spot_id" in meta.columns and "spot_id" in ints.columns, "spot_id not found in MSI tables"
    ints = ints.set_index("spot_id").sort_index()
    meta = meta.set_index("spot_id").loc[ints.index]
    # Load transformed MSI coordinates (in HIRES space)
    T = pd.read_csv(transformed_csv)
    # Accept columns: spot_id, x_new, y_new (MAGPIE outputs transformed.csv)
    # If different names, adapt here:
    cols = {c.lower(): c for c in T.columns}
    xcol = cols.get("x_new") or cols.get("x") or cols.get("pxl_col_in_fullres") or list(T.columns)[1]
    ycol = cols.get("y_new") or cols.get("y") or cols.get("pxl_row_in_fullres") or list(T.columns)[2]
    idcol = cols.get("spot_id") or list(T.columns)[0]
    T = T.set_index(idcol).loc[ints.index]
    points_hires = T[[xcol, ycol]].to_numpy().astype(float)
    
    # Convert MSI coords from HIRES to FULLRES space to match HD bins
    points = points_hires / tissue_hires_scalef
    print(f"[info] MSI coords converted: HIRES [{points_hires[:,0].min():.0f}-{points_hires[:,0].max():.0f}] → FULLRES [{points[:,0].min():.0f}-{points[:,0].max():.0f}]")
    
    # Diagnostic: compare transformed point ranges against H&E image size to catch scaling/origin issues
    try:
        from skimage.io import imread
        he = imread(Path(spaceranger_dir) / "spatial" / "tissue_hires_image.png")
        hh, ww = he.shape[:2]
        # Check in HIRES space (original coords)
        outside = (((points_hires[:, 0] < 0) | (points_hires[:, 0] > ww) | (points_hires[:, 1] < 0) | (points_hires[:, 1] > hh)).mean())
        print(
            f"[diagnostic][collapse] HE size: {ww}x{hh} "
            f"points (hires) x:[{points_hires[:,0].min():.1f},{points_hires[:,0].max():.1f}] "
            f"y:[{points_hires[:,1].min():.1f},{points_hires[:,1].max():.1f}] "
            f"| outside={outside:.3f}"
        )
    except Exception as _e:
        print("[diagnostic][collapse] could not read hires image for bounds check; skipping")
    # Load HD bins
    bins = pd.read_csv(visium_hd_bins)
    bins_xy = bins[["pxl_col_in_fullres","pxl_row_in_fullres"]].to_numpy().astype(float)
    distances, idx = nearest_bins(points, bins_xy)
    bin_ids = bins["bin_id"].astype(str).to_numpy()
    assigned_bins = bin_ids[idx]
    
    # Print mapping statistics
    print(f"[mapping] Stats: n={len(assigned_bins)}, median_dist={np.median(distances):.2f}px, "
          f"max_dist={np.max(distances):.2f}px, mean_dist={np.mean(distances):.2f}px")
    # Build MEX: rows = MSI peaks (features), columns = HD bin barcodes (uniq_bins order).
    # All intensity columns from MSI_intensities.csv (excluding spot_id, already index-only).
    features = list(ints.columns)
    # Aggregate per bin: sum (or mean—here sum is standard for counts)
    # Construct a sparse matrix (features x barcodes)
    # Map bin id -> column index
    uniq_bins, inv = np.unique(assigned_bins, return_inverse=True)
    barcodes = uniq_bins
    # Build rows/cols for COO
    rows = []
    cols = []
    data = []
    for i_peak, peak in enumerate(features):
        vals = ints[peak].to_numpy()
        # sum per bin using bincount
        summed = np.bincount(inv, weights=vals, minlength=len(uniq_bins))
        nz = np.nonzero(summed)[0]
        rows.extend([i_peak]*len(nz))
        cols.extend(nz.tolist())
        data.extend(summed[nz].tolist())
    X = sparse.csr_matrix((data,(rows,cols)), shape=(len(features), len(uniq_bins)))
    # Write MEX under spaceranger_hd/
    write_mex(out_spaceranger, X, barcodes, features)
    # Copy/replicate spatial folder (images + scalefactors) from original Visium HD run
    src_spatial = Path(spaceranger_dir) / "spatial"
    dst_spatial = out_spaceranger / "spatial"
    dst_spatial.mkdir(exist_ok=True)
    # Minimal: copy common files if present
    for fname in ["tissue_hires_image.png","tissue_lowres_image.png","scalefactors_json.json",
                  "aligned_fiducials.jpg","detected_tissue_image.jpg"]:
        p = src_spatial / fname
        if p.exists():
            (dst_spatial / fname).write_bytes(p.read_bytes())
    # tissue_positions_list.csv must align rows with `barcodes` order (Visium-style 6 columns, no header).
    # BUGFIX: do not use bins.loc[isin] — that follows bins file order, not barcodes order.
    bins_by_id = bins.set_index(bins["bin_id"].astype(str))
    if "in_tissue" in bins_by_id.columns:
        in_tissue_vals = bins_by_id.loc[barcodes, "in_tissue"].to_numpy()
    else:
        in_tissue_vals = np.ones(len(barcodes), dtype=int)
    tpos = pd.DataFrame({
        "barcode": barcodes,
        "in_tissue": in_tissue_vals,
        "array_row": 0,  # placeholders (not used by all readers)
        "array_col": 0,
        "pxl_col_in_fullres": bins_by_id.loc[barcodes, "pxl_col_in_fullres"].to_numpy(),
        "pxl_row_in_fullres": bins_by_id.loc[barcodes, "pxl_row_in_fullres"].to_numpy(),
    })
    tpos.to_csv(dst_spatial/"tissue_positions_list.csv", header=False, index=False)
    print(f"Wrote Space Ranger–style MEX to: {out_spaceranger}")

if __name__ == "__main__":
    main()