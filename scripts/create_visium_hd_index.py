#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a single spatial index table (bins or spots) for downstream HD / Visium steps.

Output columns: bin_id, pxl_col_in_fullres, pxl_row_in_fullres, in_tissue, bin_size

Modes
-----
* ``--visium_mode classic`` — use only classic Space Ranger ``tissue_positions*`` (spot-level geometry).
* ``--visium_mode hd`` — prefer HD bin-level position files (glob on ``--bin_size``). If none are found:
  * default: exit with an error (true HD geometry was expected);
  * only if ``--allow_classic_fallback`` is set: write classic spot positions instead (debug / migration only — not equivalent to HD bins).

Snakemake sets ``--out_csv`` (e.g. ``output/{sample}/visium_hd_square_{bin_size}um_bins.csv``)
and ``--visium_mode`` from top-level ``visium_mode``. When ``visium_mode`` is ``hd`` (default
for the HD workflow), classic substitution is off unless the workflow passes
``--allow_classic_fallback`` — e.g. when ``visium.hd.require_hd`` is ``false`` and/or
``visium.hd.allow_classic_fallback`` is ``true`` in ``config.yaml``.
"""
import argparse
from pathlib import Path
import sys
import pandas as pd


# ------------------------------------------------------------
# Discover and load Space Ranger / Visium position tables
# ------------------------------------------------------------
def find_hd_positions(spaceranger_dir: Path, bin_size: str):
    """
    Try common HD output filenames. Space Ranger HD may store per-bin positions
    under spatial/ or the run root. First glob match wins (see main()).
    """
    spatial = spaceranger_dir / "spatial"
    candidates = list(spatial.glob(f"*{bin_size}*pos*.csv")) + \
                 list(spatial.glob(f"*{bin_size}*pos*.tsv")) + \
                 list(spatial.glob(f"*{bin_size}*pos*.parquet")) + \
                 list(spaceranger_dir.glob(f"*{bin_size}*positions*.csv")) + \
                 list(spaceranger_dir.glob(f"*{bin_size}*positions*.tsv")) + \
                 list(spaceranger_dir.glob(f"*{bin_size}*positions*.parquet"))
    return candidates


def load_positions(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported positions file: {path}")


def classic_positions(spatial_dir: Path) -> pd.DataFrame:
    # Classic Visium: no header row (6 columns per Space Ranger convention).
    pos1 = spatial_dir / "tissue_positions_list.csv"
    pos2 = spatial_dir / "tissue_positions.csv"
    for p in (pos1, pos2):
        if p.exists():
            df = pd.read_csv(p, header=None)
            # classic columns: barcode, in_tissue, array_row, array_col, pxl_col, pxl_row
            df.columns = ["barcode","in_tissue","array_row","array_col","pxl_col_in_fullres","pxl_row_in_fullres"]
            df["bin_size"] = "classic"
            df = df.rename(columns={"barcode": "bin_id"})
            return df[["bin_id","pxl_col_in_fullres","pxl_row_in_fullres","in_tissue","bin_size"]]
    raise FileNotFoundError("Could not find classic tissue_positions_{list}.csv")


def write_hd_table(candidates0: Path, bin_size: str, out_csv: str) -> None:
    """Normalize HD positions file to standard columns and write CSV."""
    df = load_positions(candidates0)
    cols = {c.lower(): c for c in df.columns}

    def pick(*keys):
        for k in keys:
            if k in cols:
                return cols[k]
        return None

    col_x = pick("pxl_col_in_fullres", "pxl_col", "x", "col", "pixel_x")
    col_y = pick("pxl_row_in_fullres", "pxl_row", "y", "row", "pixel_y")
    col_id = pick("bin_id", "barcode", "spot_id", "id")
    col_in = pick("in_tissue", "intissue", "in")
    if col_x is None or col_y is None:
        raise ValueError(f"Positions file lacks pixel columns: {candidates0}")
    if col_id is None:
        df["bin_id"] = [f"s_{bin_size}_{i:08d}" for i in range(len(df))]
        col_id = "bin_id"
    if col_in is None:
        df["in_tissue"] = 1
        col_in = "in_tissue"
    out = df[[col_id, col_x, col_y, col_in]].copy()
    out.columns = ["bin_id", "pxl_col_in_fullres", "pxl_row_in_fullres", "in_tissue"]
    out["bin_size"] = bin_size
    print(
        f"[diagnostic][hd_bins] n={len(out)} "
        f"x:[{out['pxl_col_in_fullres'].min():.1f},{out['pxl_col_in_fullres'].max():.1f}] "
        f"y:[{out['pxl_row_in_fullres'].min():.1f},{out['pxl_row_in_fullres'].max():.1f}] "
        f"bin_size={bin_size}"
    )
    out.to_csv(out_csv, index=False)
    print(f"[HD] wrote {out_csv} from {candidates0}")


def main():
    ap = argparse.ArgumentParser(
        description="Build bin/spot index CSV from Visium HD or classic Space Ranger outputs."
    )
    ap.add_argument("--spaceranger_dir", required=True, help="Space Ranger output directory (contains spatial/)")
    ap.add_argument("--sample", required=True, help="Sample id (for logging only)")
    ap.add_argument("--bin_size", default="008", help='HD bin size string for globs: "002"|"008"|"016" (classic mode ignores HD files)')
    ap.add_argument("--out_csv", required=True, help="Output CSV path")
    ap.add_argument(
        "--visium_mode",
        choices=["hd", "classic"],
        default="hd",
        help="hd: require HD bin table unless --allow_classic_fallback. classic: classic tissue_positions only.",
    )
    ap.add_argument(
        "--allow_classic_fallback",
        action="store_true",
        help="When visium_mode=hd: if no HD positions file matches, use classic tissue_positions (not HD geometry).",
    )
    args = ap.parse_args()
    sr = Path(args.spaceranger_dir)
    spatial = sr / "spatial"

    # ----- Classic-only path: never use HD globs -----
    if args.visium_mode == "classic":
        try:
            out = classic_positions(spatial)
            out.to_csv(args.out_csv, index=False)
            print(f"[classic] wrote {args.out_csv} (visium_mode=classic)")
        except Exception as e:
            print(f"ERROR: visium_mode=classic but could not read positions in {sr}: {e}", file=sys.stderr)
            sys.exit(1)
        return

    # ----- HD path -----
    candidates = find_hd_positions(sr, args.bin_size)
    if candidates:
        write_hd_table(candidates[0], args.bin_size, args.out_csv)
        return

    # No HD table found: classic substitute only with explicit opt-in
    if args.allow_classic_fallback:
        try:
            out = classic_positions(spatial)
            out.to_csv(args.out_csv, index=False)
            print(
                f"[warn] No HD positions matched for bin_size={args.bin_size}; "
                f"wrote classic spot geometry to {args.out_csv} (--allow_classic_fallback). "
                "This is not equivalent to true HD bins.",
                file=sys.stderr,
            )
        except Exception as e:
            print(f"ERROR: HD files missing and classic fallback failed: {e}", file=sys.stderr)
            sys.exit(1)
        return

    print(
        f"ERROR: visium_mode=hd but no HD position files matched under {sr} "
        f"(bin_size={args.bin_size}). "
        "Add/install the HD outputs, fix paths, or set visium.hd.require_hd: false "
        "or visium.hd.allow_classic_fallback: true only if you intentionally accept classic spot geometry.",
        file=sys.stderr,
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
