#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pool Visium HD gene counts onto each MSI observation and build an integrated AnnData.

Pipeline role (vs collapse_to_visiumhd_bins.py)
----------------------------------------------
* ``collapse_to_visiumhd_bins.py`` puts **MSI** onto the **HD barcode grid** (nearest-bin,
  Space Ranger–style MEX) — mainly for HD-aligned visualization / tools.
* **This script** does the opposite direction for **RNA**: each MSI row is a **point** in
  fullres space (coregistered center). All **HD bin centers** within a **search disk**
  (radius ``msi.radius_um`` in µm → fullres px via ``microns_per_pixel``) are selected;
  their Space Ranger **gene counts** are combined with **sum / mean / median**. HD bins are
  **not** modeled as tiles here—only center coordinates. MSI metabolites are attached per row
  in ``obsm['msi']`` (not split across bins).

  *Future extension:* footprint–tile overlap (MSI square/disk ∩ HD square) can replace the
  center-in-disk rule once specified.

Visium HD vs classic Visium (``msi.radius_um``)
-----------------------------------------------
The **same radius-based pooling** idea as Visium+MSI applies, but **units and scale differ**:

* **Classic Visium** uses ~55 µm spots; a large ``radius_um`` matched that geometry.
* **Visium HD** uses a **dense square grid** (2 / 8 / 16 µm); bin centers and your MSI
  points both live in **fullres** pixels after the usual HIRES→FULLRES step.

**Accurate HD use** requires a correct **fullres** ``microns_per_pixel`` (Space Ranger HD
writes ``microns_per_pixel`` in ``scalefactors_json.json``). Radius in µm is converted as
``radius_px = msi_radius_um / microns_per_pixel``, then ``cKDTree.query_ball_point`` selects
all bin **centers** inside that disk.

Set ``msi.radius_um`` to the **search-disk radius** in µm on bin centers (not HD tile overlap).
For ~20 µm ion-image pixels, **10 µm** is a natural default (≈ radius / half-width of the
pixel). Use a larger value if you want to pool more HD bins around each MSI point. Not a
classic Visium spot diameter unless that matches your MSI extent. With 2 µm HD bins, even
~10 µm already includes **many** bin centers;
``aggregation.method`` controls how their gene counts are combined.

Coordinate convention (same as collapse_to_visiumhd_bins.py)
--------------------------------------------------------------
* ``transformed.csv`` columns ``x_new``/``y_new`` (or ``x``/``y``) are in **HIRES** pixels
  (Visium ``tissue_hires_image.png`` space).
* HD bin CSV from ``create_visium_hd_index.py`` uses **FULLRES** ``pxl_col_in_fullres`` /
  ``pxl_row_in_fullres``.
* Conversion: ``fullres = hires / tissue_hires_scalef`` from ``spatial/scalefactors_json.json``.

Row alignment
-------------
* **MSI intensities** define observation order and which ``spot_id`` values must appear in
  ``transformed.csv`` (same rule as ``collapse_to_visiumhd_bins.py``).

Config / Snakemake (see config.yaml + Snakefile)
-------------------------------------------------
* ``msi.radius_um`` → ``--msi_radius_um`` (search-disk radius in µm for ``query_ball_point``
  on **bin centers**; not an MSI tile footprint in the current implementation).
* ``aggregation.method`` → ``--agg_method`` (``sum`` | ``mean`` | ``median``).
* Paths for bins, Space Ranger HD matrix, intensities, and ``spaceranger_dir`` are passed
  explicitly from the workflow so they stay consistent with ``visium.*`` / ``msi.*`` in YAML.

Outputs (under ``output/<sample>/{bin_size}um/`` by default)
------------------------------------------------------------
* ``hd_bins_per_msi_spot.csv`` — long table: each MSI ↔ HD bin within radius (``msi_idx``,
  ``spot_id``, ``bin_id``, ``distance_px``, coordinates).
* ``aggregated_gene_expression.h5ad`` — genes × MSI observations (``obs['n_hd_bins']``);
  ``uns['magpie_hd_aggregate']`` records pooling parameters (radius, ``microns_per_pixel``, etc.).
* ``integrated_multimodal.h5ad`` — same as aggregated plus ``obsm['msi']`` (raw intensity
  matrix) and ``uns['msi_features']``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree


def _hd_bin_pitch_um(bin_size_code: str) -> float | None:
    """Map CLI bin_size '002'|'008'|'016' to nominal square pitch in µm."""
    return {"002": 2.0, "008": 8.0, "016": 16.0}.get(bin_size_code)


def _resolve_microns_per_pixel(sf: dict, cli_mpp: float | None) -> tuple[float, str]:
    """
    Fullres µm per pixel for converting msi_radius_um → fullres pixels.

    Space Ranger **Visium HD** typically provides ``microns_per_pixel`` directly. If absent,
    derive from ``bin_size_um`` / ``spot_diameter_fullres`` (bin extent in fullres px).
    Last resort: classic Visium fiducial scaling (may be wrong for HD — override with CLI).
    """
    if cli_mpp is not None and cli_mpp > 0:
        return float(cli_mpp), "CLI --microns_per_pixel"
    mpp = sf.get("microns_per_pixel")
    if mpp is not None and float(mpp) > 0:
        return float(mpp), "scalefactors_json['microns_per_pixel']"
    bs = sf.get("bin_size_um")
    sd = sf.get("spot_diameter_fullres")
    if bs is not None and sd is not None and float(sd) > 0:
        v = float(bs) / float(sd)
        return v, "bin_size_um / spot_diameter_fullres (HD fallback)"
    fid = sf.get("fiducial_diameter_fullres")
    if fid is not None and float(fid) > 0:
        return (
            120.0 / float(fid),
            "120 µm / fiducial_diameter_fullres (classic Visium fallback; verify for HD)",
        )
    return (
        0.5,
        "hardcoded default 0.5 (unsafe — add microns_per_pixel to JSON or use --microns_per_pixel)",
    )


def _sanity_check_hd_scale(
    sf: dict,
    microns_per_pixel: float,
    bin_size_code: str,
) -> None:
    """Log consistency between JSON and chosen HD resolution when keys exist."""
    pitch = _hd_bin_pitch_um(bin_size_code)
    if pitch is None:
        return
    bs_json = sf.get("bin_size_um")
    sd_json = sf.get("spot_diameter_fullres")
    if bs_json is not None and sd_json is not None:
        try:
            bs_f = float(bs_json)
            sd_f = float(sd_json)
            implied = bs_f / sd_f
            extent_um = sd_f * microns_per_pixel
            print(
                f"[info] HD scale check: JSON bin_size_um={bs_f:g}, spot_diameter_fullres={sd_f:.3f} px "
                f"→ {extent_um:.3f} µm extent along one axis; bin_size_um/spot_diameter_fullres={implied:.6f} µm/px",
                flush=True,
            )
            if abs(implied - microns_per_pixel) / microns_per_pixel > 0.08:
                print(
                    "[warn] microns_per_pixel differs from bin_size_um/spot_diameter_fullres by >8%. "
                    "Prefer explicit microns_per_pixel from JSON or set --microns_per_pixel.",
                    flush=True,
                )
            if abs(bs_f - pitch) > 0.51:
                print(
                    f"[warn] scalefactors bin_size_um ({bs_f:g}) vs --bin_size {bin_size_code} (~{pitch:g} µm) mismatch.",
                    flush=True,
                )
        except (TypeError, ValueError):
            pass


# ------------------------------------------------------------
# Map each MSI point → all HD bins whose centers lie within radius (fullres px)
# ------------------------------------------------------------
def find_hd_bins_within_radius(
    msi_coords: np.ndarray,
    hd_coords: np.ndarray,
    hd_bin_ids: np.ndarray,
    radius_px: float,
) -> pd.DataFrame:
    """
    For each MSI point, find HD bin indices within Euclidean distance ``radius_px`` in
    fullres pixel space (bin table positions are treated as point locations).

    Returns columns: msi_idx, bin_id, distance_px (spot_id / msi_x / msi_y added in main).
    """
    tree = cKDTree(hd_coords)
    records = []
    for i, (x, y) in enumerate(msi_coords):
        indices = tree.query_ball_point(np.array([x, y], dtype=float), r=radius_px)
        for idx in indices:
            dist = float(np.hypot(hd_coords[idx, 0] - x, hd_coords[idx, 1] - y))
            records.append(
                {
                    "msi_idx": i,
                    "bin_id": hd_bin_ids[idx],
                    "distance_px": dist,
                }
            )
    return pd.DataFrame(records)


# ------------------------------------------------------------
# Pool HD gene expression onto each MSI observation
# ------------------------------------------------------------
def _normalize_obs_index(adata: ad.AnnData) -> None:
    """Ensure obs_names are strings for reliable matching to bin_id from CSV."""
    adata.obs_names = adata.obs_names.astype(str)


def aggregate_gene_expression(
    adata_hd: ad.AnnData,
    mapping_df: pd.DataFrame,
    msi_spot_ids: list,
    method: str = "sum",
) -> ad.AnnData:
    """
    For each MSI row, subset HD ``adata_hd`` to mapped bins and reduce along obs (sum/mean/median).

    Duplicate ``bin_id`` rows for the same ``msi_idx`` are deduplicated (first occurrence wins)
    so AnnData indexing stays well-defined.
    """
    n_msi = len(msi_spot_ids)
    n_genes = adata_hd.n_vars

    agg_matrix = np.zeros((n_msi, n_genes), dtype=np.float32)
    n_bins_per_spot = np.zeros(n_msi, dtype=int)

    _normalize_obs_index(adata_hd)
    valid_bins = set(adata_hd.obs_names)

    for msi_idx, group in mapping_df.groupby("msi_idx"):
        bin_ids = group["bin_id"].astype(str).values
        seen: set[str] = set()
        valid_bin_ids: list[str] = []
        for b in bin_ids:
            if b in valid_bins and b not in seen:
                seen.add(b)
                valid_bin_ids.append(b)
        if len(valid_bin_ids) == 0:
            continue

        expr = adata_hd[valid_bin_ids, :].X
        if sparse.issparse(expr):
            expr = expr.toarray()

        if method == "sum":
            vec = expr.sum(axis=0)
        elif method == "mean":
            vec = expr.mean(axis=0)
        elif method == "median":
            vec = np.median(expr, axis=0)
        else:
            raise ValueError(f"Unknown agg method: {method}")

        agg_matrix[msi_idx, :] = np.asarray(vec, dtype=np.float32).ravel()
        n_bins_per_spot[msi_idx] = len(valid_bin_ids)

    return ad.AnnData(
        X=agg_matrix,
        obs=pd.DataFrame(
            {"spot_id": msi_spot_ids, "n_hd_bins": n_bins_per_spot},
            index=msi_spot_ids,
        ),
        var=adata_hd.var.copy(),
    )


# ------------------------------------------------------------
# Load 10x HDF5 or MEX directory (Visium HD Space Ranger)
# ------------------------------------------------------------
def load_visium_hd_expression(path: Path) -> ad.AnnData:
    """Load Visium HD gene expression from ``.h5`` or ``filtered_feature_bc_matrix`` directory."""
    if path.suffix == ".h5":
        return sc.read_10x_h5(path)
    if path.is_dir():
        return sc.read_10x_mtx(path)
    raise ValueError(f"Expected .h5 file or matrix directory, got: {path}")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Pool Visium HD gene counts onto MSI observations: bin centers within msi_radius_um (µm) "
            "of each coregistered MSI point (fullres); sum/mean/median over those bins."
        )
    )
    ap.add_argument("--sample", required=True, help="Sample ID (logging / default paths)")
    ap.add_argument(
        "--bin_size",
        default="008",
        choices=["002", "008", "016"],
        help='HD bin resolution string (matches Space Ranger folder, e.g. "002" for 2 µm bins)',
    )
    ap.add_argument(
        "--transformed_csv",
        default=None,
        help="Coregistered MSI coords (HIRES); default output/<sample>/transformed.csv",
    )
    ap.add_argument(
        "--visium_hd_bins",
        default=None,
        help="HD bin index CSV from create_visium_hd_index.py (FULLRES columns)",
    )
    ap.add_argument(
        "--visium_hd_matrix",
        default=None,
        help="Visium HD gene matrix: filtered_feature_bc_matrix.h5 or MEX directory",
    )
    ap.add_argument(
        "--msi_intensities",
        default=None,
        help="MSI intensities CSV (spot_id + m/z columns); defines obs order",
    )
    ap.add_argument(
        "--spaceranger_dir",
        default=None,
        help="Space Ranger run root (parent of spatial/); default input/<sample>/visium",
    )
    ap.add_argument(
        "--msi_radius_um",
        type=float,
        default=10.0,
        help=(
            "Search radius (µm) on bin centers (default 10 ≈ half of ~20µm MSI pixel width). "
            "config: msi.radius_um"
        ),
    )
    ap.add_argument(
        "--microns_per_pixel",
        type=float,
        default=None,
        help="Fullres µm per pixel; default from scalefactors_json.json key microns_per_pixel",
    )
    ap.add_argument(
        "--agg_method",
        default="sum",
        choices=["sum", "mean", "median"],
        help="How to combine gene counts over HD bins per MSI spot (config: aggregation.method)",
    )
    ap.add_argument(
        "--out_dir",
        default=None,
        help="Output directory; default output/<sample>/<bin_size>um",
    )
    args = ap.parse_args()

    if args.msi_radius_um <= 0:
        raise ValueError("--msi_radius_um must be positive")

    sample = args.sample
    bin_size = args.bin_size

    transformed_csv = (
        Path(args.transformed_csv)
        if args.transformed_csv
        else Path(f"output/{sample}/transformed.csv")
    )
    visium_hd_bins = (
        Path(args.visium_hd_bins)
        if args.visium_hd_bins
        else Path(f"output/{sample}/visium_hd_square_{bin_size}um_bins.csv")
    )
    visium_hd_matrix = (
        Path(args.visium_hd_matrix)
        if args.visium_hd_matrix
        else Path(f"input/{sample}/visium/{bin_size}um/filtered_feature_bc_matrix.h5")
    )
    msi_intensities = (
        Path(args.msi_intensities)
        if args.msi_intensities
        else Path(f"input/{sample}/msi/MSI_intensities.csv")
    )
    spaceranger_dir = Path(
        args.spaceranger_dir if args.spaceranger_dir else f"input/{sample}/visium"
    )

    out_dir = Path(args.out_dir) if args.out_dir else Path(f"output/{sample}/{bin_size}um")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[info] Using HD bin resolution: {bin_size} µm")
    print(f"[info] Output directory: {out_dir}")

    sf_path = spaceranger_dir / "spatial" / "scalefactors_json.json"
    if not sf_path.is_file():
        raise FileNotFoundError(
            f"scalefactors not found: {sf_path}. Set --spaceranger_dir to your Visium HD run root."
        )
    with open(sf_path) as f:
        sf = json.load(f)

    tissue_hires_scalef = float(sf.get("tissue_hires_scalef", 0.5))
    microns_per_pixel, mpp_source = _resolve_microns_per_pixel(
        sf, float(args.microns_per_pixel) if args.microns_per_pixel is not None else None
    )
    if microns_per_pixel <= 0:
        raise ValueError("microns_per_pixel must be positive")

    if "default 0.5" in mpp_source or "unsafe" in mpp_source:
        print(f"[warn] {mpp_source}", flush=True)
    elif "classic Visium fallback" in mpp_source:
        print(f"[warn] {mpp_source}", flush=True)

    _sanity_check_hd_scale(sf, microns_per_pixel, bin_size)

    pitch_um = _hd_bin_pitch_um(bin_size)
    if pitch_um is not None:
        bins_across_diameter = (2.0 * args.msi_radius_um) / pitch_um
        print(
            f"[info] Nominal HD bin pitch (--bin_size {bin_size}): ~{pitch_um:g} µm; "
            f"2×radius_um / pitch ≈ {bins_across_diameter:.1f} bins across disk diameter (informative only)",
            flush=True,
        )

    radius_px = args.msi_radius_um / microns_per_pixel

    print(f"[info] MSI radius: {args.msi_radius_um} µm → {radius_px:.1f} px (fullres)")
    print(f"[info] microns_per_pixel (fullres): {microns_per_pixel}  ({mpp_source})")
    print(f"[info] HIRES → FULLRES: divide by tissue_hires_scalef = {tissue_hires_scalef:.6f}")
    print(f"[info] Aggregation: {args.agg_method}")

    # --- Load MSI intensities first (canonical spot_id order, same as collapse script) ---
    print("\n[step 1] Loading data...")
    msi_int = pd.read_csv(msi_intensities)
    if "spot_id" not in msi_int.columns:
        raise ValueError(f"MSI intensities must contain spot_id; got: {list(msi_int.columns)}")
    msi_int = msi_int.set_index("spot_id")
    msi_int.index = msi_int.index.astype(str)

    msi_df = pd.read_csv(transformed_csv)
    colmap = {c.lower(): c for c in msi_df.columns}
    xcol = colmap.get("x_new") or colmap.get("x") or colmap.get("pxl_col_in_fullres")
    ycol = colmap.get("y_new") or colmap.get("y") or colmap.get("pxl_row_in_fullres")
    idcol = colmap.get("spot_id")
    if xcol is None or ycol is None:
        raise ValueError(
            f"transformed_csv needs x/y (or x_new/y_new); got columns: {list(msi_df.columns)}"
        )
    if idcol is None:
        raise ValueError(f"transformed_csv needs spot_id; got columns: {list(msi_df.columns)}")

    msi_df = msi_df.copy()
    msi_df[idcol] = msi_df[idcol].astype(str)
    msi_df = msi_df.set_index(idcol, drop=True)
    msi_df.index.name = "spot_id"
    missing = msi_int.index.difference(msi_df.index)
    if len(missing) > 0:
        raise ValueError(
            f"{len(missing)} spot_id(s) in MSI intensities have no row in transformed_csv "
            f"(showing up to 5): {list(missing[:5])!r}"
        )
    msi_df = msi_df.loc[msi_int.index]
    msi_spot_ids = msi_df.index.astype(str).tolist()

    msi_coords_hires = msi_df[[xcol, ycol]].values.astype(float)
    msi_coords = msi_coords_hires / tissue_hires_scalef

    print(f"  MSI observations (from intensities order): {len(msi_spot_ids):,}")
    print(
        f"  MSI coords HIRES x:[{msi_coords_hires[:, 0].min():.0f},{msi_coords_hires[:, 0].max():.0f}] "
        f"y:[{msi_coords_hires[:, 1].min():.0f},{msi_coords_hires[:, 1].max():.0f}]"
    )
    print(
        f"  MSI coords FULLRES x:[{msi_coords[:, 0].min():.0f},{msi_coords[:, 0].max():.0f}] "
        f"y:[{msi_coords[:, 1].min():.0f},{msi_coords[:, 1].max():.0f}]"
    )

    hd_bins = pd.read_csv(visium_hd_bins)
    req = {"pxl_col_in_fullres", "pxl_row_in_fullres", "bin_id"}
    if not req.issubset(set(hd_bins.columns)):
        raise ValueError(f"visium_hd_bins CSV must contain columns {req}; got {list(hd_bins.columns)}")
    hd_coords = hd_bins[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    hd_bin_ids = hd_bins["bin_id"].astype(str).values
    print(f"  HD bins: {len(hd_bin_ids):,}")

    print(f"  Visium HD expression: {visium_hd_matrix}")
    adata_hd = load_visium_hd_expression(visium_hd_matrix)
    _normalize_obs_index(adata_hd)
    print(f"  HD matrix: {adata_hd.n_obs:,} barcodes × {adata_hd.n_vars:,} features")

    n_bins_csv = len(set(hd_bin_ids))
    n_intersect = len(set(hd_bin_ids) & set(adata_hd.obs_names))
    print(f"  bin_id overlap (CSV vs expression): {n_intersect:,} / {n_bins_csv:,} unique CSV bin_ids")
    if n_intersect == 0:
        print(
            "[warn] No overlapping bin_id between visium_hd_bins CSV and expression matrix obs_names. "
            "Check barcode string format (e.g. prefixes) and that bin_size matches the .h5 run.",
            flush=True,
        )

    print("\n[step 2] HD bins within radius per MSI observation (cKDTree.query_ball_point)...")
    mapping_df = find_hd_bins_within_radius(msi_coords, hd_coords, hd_bin_ids, radius_px)
    mapping_df["spot_id"] = mapping_df["msi_idx"].map(lambda i: msi_spot_ids[i])
    mapping_df["msi_x"] = mapping_df["msi_idx"].map(lambda i: msi_coords[i, 0])
    mapping_df["msi_y"] = mapping_df["msi_idx"].map(lambda i: msi_coords[i, 1])

    mapping_path = out_dir / "hd_bins_per_msi_spot.csv"
    mapping_df.to_csv(mapping_path, index=False)

    if mapping_df.empty:
        print(
            "[warn] No MSI–HD bin pairs within radius (empty mapping). "
            "Raise msi.radius_um, check registration, scalefactors, and bin_id overlap with the .h5.",
            flush=True,
        )

    bins_per_spot = mapping_df.groupby("msi_idx").size()
    print(
        f"  HD bins per MSI spot: median={bins_per_spot.median():.0f}, mean={bins_per_spot.mean():.1f}"
    )
    print(f"  Total MSI–bin pairs: {len(mapping_df):,}")
    print(f"  Saved: {mapping_path}")

    print(f"\n[step 3] Aggregating genes ({args.agg_method}) onto MSI observations...")
    adata_agg = aggregate_gene_expression(adata_hd, mapping_df, msi_spot_ids, args.agg_method)
    adata_agg.obsm["spatial"] = msi_coords
    # Provenance for methods / debugging (disk-on-centers pooling only)
    adata_agg.uns["magpie_hd_aggregate"] = {
        "pooling": "disk_centers",
        "description": "HD bin centers within radius_px of each MSI point (fullres)",
        "msi_radius_um": float(args.msi_radius_um),
        "radius_px_fullres": float(radius_px),
        "microns_per_pixel": float(microns_per_pixel),
        "microns_per_pixel_source": mpp_source,
        "agg_method": args.agg_method,
        "bin_size": bin_size,
        "tissue_hires_scalef": float(tissue_hires_scalef),
    }

    agg_path = out_dir / "aggregated_gene_expression.h5ad"
    adata_agg.write(agg_path)
    print(f"  Saved: {agg_path}")

    print("\n[step 4] integrated_multimodal.h5ad (obsm['msi'] = raw MSI matrix)...")
    msi_int_aligned = msi_int.loc[msi_spot_ids, :].values
    adata_integrated = adata_agg.copy()
    adata_integrated.obsm["msi"] = msi_int_aligned
    adata_integrated.uns["msi_features"] = list(msi_int.columns)
    adata_integrated.obs["total_msi_intensity"] = msi_int_aligned.sum(axis=1)
    adata_integrated.obs["total_gene_counts"] = np.asarray(adata_integrated.X.sum(axis=1)).ravel()

    integrated_path = out_dir / "integrated_multimodal.h5ad"
    adata_integrated.write(integrated_path)
    print(f"  Saved: {integrated_path}")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"MSI spots with ≥1 HD bin in radius: {(adata_integrated.obs['n_hd_bins'] > 0).sum():,} / {len(msi_spot_ids):,}")
    print(f"Median HD bins per spot (within radius): {adata_integrated.obs['n_hd_bins'].median():.0f}")
    print("Outputs:")
    print(f"  - {mapping_path}")
    print(f"  - {agg_path}")
    print(f"  - {integrated_path}")


if __name__ == "__main__":
    main()
