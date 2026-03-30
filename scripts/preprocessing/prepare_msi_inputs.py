#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Normalize vendor MSI tables into MAGPIE's expected CSV layout.


**MAGPIE expects the following files:**

  * ``MSI_metadata.csv`` — at least ``spot_id``, ``x``, ``y`` (one row per MSI pixel/spot).
  * ``MSI_intensities.csv`` — ``spot_id`` plus one column per feature (names like ``mz-<value>``).

  Raw exports often differ: transposed matrices (m/z as rows), alternate separators,
  mixed spot labels, or extra preamble lines. This script maps those inputs to the
  canonical column names and row order, optionally subsets metabolites, and writes
  ``MSI_intensities.csv`` / ``MSI_metadata.csv`` (or paths you choose).

**Alignment rule**
  Rows are **ordered and filtered by metadata**: every ``spot_id`` in metadata is kept;
  intensities are reindexed to that order. Spots present only in intensities are dropped.
  Missing intensity cells after alignment are filled with **0** (not NaN).

**How to run** (from repository root)::

    # Typical: sample folder layout input/<SAMPLE>/msi/MSI_*_org.csv → MSI_*.csv
    python scripts/prepare_msi_inputs.py --sample-dir input/MySample

    # Explicit paths
    python scripts/prepare_msi_inputs.py \\
        --intensities path/to/raw_intensities.csv \\
        --metadata path/to/raw_metadata.csv \\
        --out-intensities MSI_intensities.csv \\
        --out-metadata MSI_metadata.csv

    # Keep only features listed (first column of CSV or one name per line)
    python scripts/prepare_msi_inputs.py --sample-dir input/MySample \\
        --metabolites input/MySample/msi/MSI_selected_peaks.txt

See also: ``docs/inputs.md`` for the MAGPIE MSI table contract.
"""
from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import pandas as pd


def find_header_row(path: Path, prefixes: tuple[str, ...]) -> tuple[int, str | None]:
    """
    Find the first line that looks like a table header (starts with any of ``prefixes``),
    skipping empty lines. If none match, fall back to the **first non-empty line** as the
    header row index (may mis-parse files whose real header is preceded by title text).
    """
    header_row: int | None = None
    header_line: str | None = None
    with open(path, "r", encoding="utf-8-sig") as f:
        for i, line in enumerate(f):
            stripped = line.strip()
            if not stripped:
                continue
            if header_row is None:
                header_line = stripped
                header_row = i
            lower = stripped.lower()
            if any(lower.startswith(p) for p in prefixes):
                return i, stripped
    return header_row or 0, header_line


def detect_sep(header_line: str | None) -> str:
    """Prefer ';' if it appears at least as often as ',' on the header line; else ','."""
    if not header_line:
        return ","
    if header_line.count(";") >= header_line.count(","):
        return ";"
    return ","


def normalize_spot_label(label: object, fallback_index: int) -> str:
    """
    Map a vendor spot label to a stable ``pixel_<n>`` style id when possible.

    Preserves values that already look like ``pixel_*``. Otherwise extracts the first
    integer substring; if none, uses ``fallback_index``.
    """
    value = str(label).strip()
    if value.lower().startswith("pixel_"):
        return value
    match = re.search(r"(\d+)", value)
    if match:
        return f"pixel_{int(match.group(1))}"
    return f"pixel_{fallback_index}"


def normalize_spot_series(series: pd.Series) -> pd.Series:
    """
    Normalize spot labels row-wise; if collisions occur after normalization, fall back to
    ``pixel_1`` … ``pixel_N`` in **current row order** (may desynchronize from another
    file if both needed unique remapping — keep row order consistent across uploads).
    """
    normalized = [
        normalize_spot_label(value, i + 1) for i, value in enumerate(series.tolist())
    ]
    if len(set(normalized)) != len(normalized):
        normalized = [f"pixel_{i}" for i in range(1, len(series) + 1)]
    return pd.Series(normalized, index=series.index)


def normalize_mz(value: object) -> str:
    """Normalize feature names toward ``mz-<...>`` for column matching."""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return text
    lower = text.lower()
    if lower.startswith("mz-"):
        return text
    if lower.startswith("m/z"):
        text = text[3:].strip()
    if lower.startswith("mz"):
        text = text[2:].lstrip("-").strip()
    return f"mz-{text}" if text else text


def load_metadata(path: Path) -> pd.DataFrame:
    """Load metadata CSV; ensure ``spot_id`` and numeric ``x`` / ``y`` when present."""
    header_row, header_line = find_header_row(
        path, ("spot", "spot_id", "x", "y")
    )
    sep = detect_sep(header_line)
    df = pd.read_csv(path, sep=sep, skiprows=header_row, engine="python")
    df.columns = [c.strip() for c in df.columns]
    if "spot_id" not in df.columns:
        df = df.rename(columns={df.columns[0]: "spot_id"})
    df["spot_id"] = normalize_spot_series(df["spot_id"])
    for col in ("x", "y"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    dup = df["spot_id"].duplicated()
    if dup.any():
        n = int(dup.sum())
        warnings.warn(
            f"Metadata contains {n} duplicate spot_id row(s) after normalization; "
            "downstream alignment may be ambiguous.",
            UserWarning,
            stacklevel=2,
        )
    return df


def load_intensities(path: Path) -> pd.DataFrame:
    """
    Load intensity matrix in one of two broad layouts:

    * **Long / standard:** first column is ``spot_id`` (or renamed to it); remaining
      columns are features (m/z).
    * **Wide transposed:** first column is m/z labels; each other column is one spot;
      the table is transposed so rows are spots and columns are features.
    """
    header_row, header_line = find_header_row(
        path, ("m/z", "mz", "spot_id", "spot")
    )
    sep = detect_sep(header_line)
    df = pd.read_csv(path, sep=sep, skiprows=header_row, engine="python")
    df.columns = [c.strip() for c in df.columns]
    first_col = df.columns[0].strip()
    first_lower = first_col.lower()

    # Transposed export: index = m/z, columns = spot ids
    if first_lower in {"m/z", "mz"} or "m/z" in first_lower:
        mz_values = df.iloc[:, 0].apply(normalize_mz)
        data = df.iloc[:, 1:]
        data.columns = [
            normalize_spot_label(name, i + 1)
            for i, name in enumerate(data.columns)
        ]
        data.index = mz_values
        data = data.T
        data.insert(0, "spot_id", data.index)
        return data.reset_index(drop=True)

    if first_col != "spot_id":
        df = df.rename(columns={first_col: "spot_id"})
    df["spot_id"] = normalize_spot_series(df["spot_id"])
    return df


def load_metabolite_list(path: str | None) -> set[str] | None:
    """Load a set of feature names to keep (normalized via ``normalize_mz``)."""
    if not path:
        return None
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Metabolite list not found: {p}")
    try:
        df = pd.read_csv(p, sep=None, engine="python")
        candidates = df.iloc[:, 0].dropna().astype(str).tolist()
    except Exception:
        with open(p, "r", encoding="utf-8-sig") as f:
            candidates = [line.strip() for line in f if line.strip()]
    out = {normalize_mz(v) for v in candidates if normalize_mz(v)}
    return out if out else None


def filter_metabolites(
    intensities_df: pd.DataFrame, metabolites: set[str] | None
) -> pd.DataFrame:
    """Restrict intensity columns to ``metabolites`` ∪ {``spot_id``}."""
    if not metabolites:
        return intensities_df
    keep_cols = ["spot_id"] + [c for c in intensities_df.columns if c in metabolites]
    mz_kept = len(keep_cols) - 1
    if mz_kept == 0:
        warnings.warn(
            "No intensity columns matched the metabolite list; output will contain "
            "only spot_id. Check names in --metabolites vs column headers.",
            UserWarning,
            stacklevel=2,
        )
    missing = metabolites - set(intensities_df.columns)
    if missing and mz_kept > 0:
        sample = ", ".join(sorted(missing)[:5])
        extra = len(missing) - 5 if len(missing) > 5 else 0
        msg = f"{len(missing)} metabolite(s) from list not found in intensities (e.g. {sample}"
        if extra > 0:
            msg += f", ... +{extra} more"
        msg += ")."
        warnings.warn(msg, UserWarning, stacklevel=2)
    return intensities_df.loc[:, keep_cols]


def align_metadata_intensities(
    metadata_df: pd.DataFrame, intensities_df: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reorder intensity rows to match metadata ``spot_id`` order; add NaN for metadata
    spots missing from intensities (filled to 0 in ``main``). Drops intensity-only spots.
    """
    metadata = metadata_df.copy()
    intensities = intensities_df.copy()

    metadata.set_index("spot_id", inplace=True)
    intensities.set_index("spot_id", inplace=True)

    if metadata.index.empty:
        raise ValueError("Metadata has no spot_id values.")

    only_int = intensities.index.difference(metadata.index)
    if len(only_int) > 0:
        warnings.warn(
            f"Dropping {len(only_int)} spot_id(s) present in intensities but not in metadata.",
            UserWarning,
            stacklevel=2,
        )

    intensities = intensities.reindex(metadata.index)

    n_missing = intensities.isna().any(axis=1).sum()
    if n_missing > 0:
        warnings.warn(
            f"{int(n_missing)} metadata spot_id row(s) had no matching intensities "
            "(filled with 0 after alignment).",
            UserWarning,
            stacklevel=2,
        )

    metadata.reset_index(inplace=True)
    intensities.reset_index(inplace=True)

    return metadata, intensities


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    """Resolve input/output paths from ``--sample-dir`` or explicit file arguments."""
    if args.sample_dir:
        sample_dir = Path(args.sample_dir)
        intensities_path = sample_dir / "msi" / "MSI_intensities_org.csv"
        metadata_path = sample_dir / "msi" / "MSI_metadata_org.csv"
        out_intensities = (
            Path(args.out_intensities)
            if args.out_intensities
            else sample_dir / "msi" / "MSI_intensities.csv"
        )
        out_metadata = (
            Path(args.out_metadata)
            if args.out_metadata
            else sample_dir / "msi" / "MSI_metadata.csv"
        )
        return intensities_path, metadata_path, out_intensities, out_metadata

    if not args.intensities or not args.metadata:
        raise ValueError("Provide --sample-dir or both --intensities and --metadata.")

    intensities_path = Path(args.intensities)
    metadata_path = Path(args.metadata)

    out_intensities = (
        Path(args.out_intensities)
        if args.out_intensities
        else intensities_path.with_name("MSI_intensities.csv")
    )
    out_metadata = (
        Path(args.out_metadata)
        if args.out_metadata
        else metadata_path.with_name("MSI_metadata.csv")
    )
    return intensities_path, metadata_path, out_intensities, out_metadata


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize raw MSI metadata and intensity CSVs into MAGPIE-ready "
            "MSI_metadata.csv and MSI_intensities.csv."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/prepare_msi_inputs.py --sample-dir input/MySample
  python scripts/prepare_msi_inputs.py -i raw_meta.csv -I raw_int.csv -o meta_out.csv -O int_out.csv
        """.strip(),
    )
    parser.add_argument(
        "--sample-dir",
        help="Sample directory containing msi/MSI_metadata_org.csv and msi/MSI_intensities_org.csv",
    )
    parser.add_argument("--intensities", "-I", help="Path to raw MSI intensities CSV")
    parser.add_argument("--metadata", "-i", help="Path to raw MSI metadata CSV")
    parser.add_argument(
        "--metabolites",
        "-m",
        help="Optional path: feature names to keep (CSV first column or plain text lines)",
    )
    parser.add_argument(
        "--out-intensities", "-O", help="Output intensities CSV (default: next to input or under msi/)"
    )
    parser.add_argument(
        "--out-metadata", "-o", help="Output metadata CSV (default: next to input or under msi/)"
    )
    args = parser.parse_args()

    intensities_path, metadata_path, out_intensities, out_metadata = resolve_paths(args)

    for label, p in (("Metadata", metadata_path), ("Intensities", intensities_path)):
        if not p.is_file():
            print(f"[error] {label} file not found: {p}", file=sys.stderr)
            sys.exit(1)

    metadata_df = load_metadata(metadata_path)
    intensities_df = load_intensities(intensities_path)

    orig_metadata_shape = metadata_df.shape
    orig_intensities_shape = intensities_df.shape

    metabolites = load_metabolite_list(args.metabolites)
    intensities_df = filter_metabolites(intensities_df, metabolites)

    metadata_df, intensities_df = align_metadata_intensities(
        metadata_df, intensities_df
    )

    intensities_df = intensities_df.fillna(0)

    cleaned_metadata_shape = metadata_df.shape
    cleaned_intensities_shape = intensities_df.shape
    print(
        f"Metadata shape: original={orig_metadata_shape}, cleaned={cleaned_metadata_shape}"
    )
    print(
        f"Intensities shape: original={orig_intensities_shape}, cleaned={cleaned_intensities_shape}"
    )

    out_intensities.parent.mkdir(parents=True, exist_ok=True)
    out_metadata.parent.mkdir(parents=True, exist_ok=True)

    metadata_df.to_csv(out_metadata, index=False)
    intensities_df.to_csv(out_intensities, index=False)

    print(
        f"Saved metadata: {out_metadata} (rows={len(metadata_df)}, cols={metadata_df.shape[1]})"
    )
    print(
        f"Saved intensities: {out_intensities} (rows={len(intensities_df)}, cols={intensities_df.shape[1]})"
    )


if __name__ == "__main__":
    main()
