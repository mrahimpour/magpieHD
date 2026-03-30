#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QC: check whether a gene (e.g. **P2ry6**) is present in the Visium HD feature reference
and whether it carries non‑zero counts in the filtered matrix.

Inputs (same as ``aggregate_to_msi_resolution.load_visium_hd_expression``)
-------------------------------------------------------------------------
* Space Ranger **filtered** output: either

  * ``filtered_feature_bc_matrix.h5``, or
  * directory ``filtered_feature_bc_matrix/`` (``matrix.mtx.gz`` + ``features.tsv.gz`` + ``barcodes.tsv.gz``).

Matching
--------
The query string is compared to ``adata.var_names`` and to every string-like column in
``adata.var`` (e.g. Ensembl id vs symbol). Matching is exact; use ``--case-insensitive``
for symbol-style queries.

“Detected”
----------
* **in_reference**: at least one feature row matches the query (and passes
  ``--gene-expression-only`` when set).
* **expressed**: total UMI/count sum for matched **Gene Expression** rows (or all matched
  rows if type column is missing) is **> 0** after optional ``--min-total-umis`` threshold.

Exit codes (only if ``--strict``): **0** = expressed, **1** = in reference but not expressed,
**2** = not in reference, **3** = usage/runtime error.

JSON file
---------
The full report is **always** written as JSON under **``output/<sample_id>/qc/``** (alongside other
pipeline outputs for that sample). If ``--matrix`` looks like ``.../<sample_id>/visium/...``, the
default file is ``output/<sample_id>/qc/visium_hd_gene_qc_<gene>.json``; otherwise
``output/unscoped/qc/visium_hd_gene_qc_<gene>.json``. Override with ``--out-json``.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def load_visium_hd_matrix(path: Path):
    """Load Visium / Visium HD filtered matrix from ``.h5`` or MEX directory."""
    if path.suffix.lower() == ".h5":
        return sc.read_10x_h5(path)
    if path.is_dir():
        return sc.read_10x_mtx(path)
    raise ValueError(f"Expected .h5 file or filtered_feature_bc_matrix directory, got: {path}")


def _norm(s: str, case_insensitive: bool) -> str:
    s = str(s).strip()
    return s.lower() if case_insensitive else s


def _feature_type_col(var) -> str | None:
    for name in ("feature_types", "feature_type", "Feature Type"):
        if name in var.columns:
            return name
    return None


def _is_gene_expression_row(var, row_i: int, type_col: str | None) -> bool:
    if type_col is None:
        return True
    v = var.iloc[row_i][type_col]
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return True
    return str(v).strip().lower() == "gene expression"


def find_matching_var_indices(
    adata,
    query: str,
    *,
    case_insensitive: bool,
    gene_expression_only: bool,
) -> list[int]:
    """Return ``var`` positions (integer indices) that match ``query``."""
    qn = _norm(query, case_insensitive)
    type_col = _feature_type_col(adata.var)
    hits: list[int] = []

    str_cols = [
        c
        for c in adata.var.columns
        if pd.api.types.is_object_dtype(adata.var[c])
        or pd.api.types.is_string_dtype(adata.var[c])
        or pd.api.types.is_categorical_dtype(adata.var[c])
    ]

    for j in range(adata.n_vars):
        if gene_expression_only and not _is_gene_expression_row(adata.var, j, type_col):
            continue
        candidates: list[str] = [str(adata.var_names[j])]
        for c in str_cols:
            val = adata.var.iloc[j][c]
            if val is None or (isinstance(val, float) and np.isnan(val)):
                continue
            candidates.append(str(val))
        for c in candidates:
            if _norm(c, case_insensitive) == qn:
                hits.append(j)
                break
    return hits


def _summarize_expression(adata, var_indices: list[int]) -> dict:
    """Total counts and bins with >0 for selected ``var`` columns."""
    if not var_indices:
        return {
            "n_feature_rows_matched": 0,
            "total_umis": 0,
            "n_bins_with_signal": 0,
            "frac_bins_with_signal": 0.0,
        }
    # Subset columns; X is obs × var
    sub = adata[:, var_indices].X
    if hasattr(sub, "sum"):
        total = float(sub.sum())
    else:
        total = float(np.asarray(sub).sum())
    # bins with any count on matched features
    if hasattr(sub, "sum"):
        col_sums = np.asarray(sub.sum(axis=1)).ravel()
    else:
        col_sums = np.asarray(sub).sum(axis=1).ravel()
    n_pos = int(np.sum(col_sums > 0))
    n_obs = adata.n_obs
    return {
        "n_feature_rows_matched": len(var_indices),
        "total_umis": total,
        "n_bins_with_signal": n_pos,
        "frac_bins_with_signal": (n_pos / n_obs) if n_obs else 0.0,
    }


def _gene_slug_for_path(gene: str) -> str:
    s = re.sub(r"[^\w.\-]+", "_", str(gene).strip(), flags=re.UNICODE).strip("_")
    return s or "gene"


def _infer_sample_id_from_matrix(matrix: Path) -> str | None:
    """
    If path contains a ``visium`` segment, return the parent directory name (MAGPIE layout:
    ``input/<sample>/visium/...``). Otherwise return None.
    """
    parts = [p for p in matrix.resolve().parts if p]
    for i, p in enumerate(parts):
        if p.lower() == "visium" and i > 0:
            return parts[i - 1]
    return None


def _default_json_out_path(matrix: Path, gene: str) -> Path:
    sample = _infer_sample_id_from_matrix(matrix)
    sub = sample if sample else "unscoped"
    return Path("output") / sub / "qc" / f"visium_hd_gene_qc_{_gene_slug_for_path(gene)}.json"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="QC: gene in Visium HD filtered matrix — reference presence and expression.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example (from repository root, with your conda env activated so scanpy is available):

  python scripts/QC/visium_hd_gene_qc.py \\
    --matrix input/<sample>/visium/002um/filtered_feature_bc_matrix.h5 \\
    --gene P2ry6 \\
    --case-insensitive

Default JSON path: output/<sample>/qc/visium_hd_gene_qc_<gene>.json when --matrix
contains .../<sample>/visium/...; else output/unscoped/qc/.... Override with --out-json.
""",
    )
    ap.add_argument(
        "--matrix",
        required=True,
        type=Path,
        help="Path to filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/ directory",
    )
    ap.add_argument(
        "--gene",
        required=True,
        help="Gene symbol or id to look up (e.g. P2ry6 or Ensembl id, depending on matrix)",
    )
    ap.add_argument(
        "--case-insensitive",
        action="store_true",
        help="Compare query to identifiers case-insensitively (recommended for symbols)",
    )
    ap.add_argument(
        "--gene-expression-only",
        action="store_true",
        help="Ignore non–Gene Expression rows (e.g. Antibody Capture) when matching",
    )
    ap.add_argument(
        "--min-total-umis",
        type=float,
        default=1.0,
        help="Minimum summed UMIs across matched feature rows to count as 'expressed' (default: 1)",
    )
    ap.add_argument(
        "--json",
        action="store_true",
        help="Print one JSON object to stdout (still prints a short human line to stderr unless --json-only)",
    )
    ap.add_argument(
        "--json-only",
        action="store_true",
        help="Print only JSON to stdout",
    )
    ap.add_argument(
        "--out-json",
        type=Path,
        default=None,
        help="QC JSON path (default: output/<sample>/qc/... from matrix path, else output/unscoped/qc/...)",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Non-zero exit: 2 = not in reference, 1 = in reference but not expressed, 0 = expressed",
    )
    args = ap.parse_args(argv)

    if not args.matrix.exists():
        print(f"ERROR: matrix path does not exist: {args.matrix}", file=sys.stderr)
        return 3

    try:
        adata = load_visium_hd_matrix(args.matrix)
    except Exception as e:
        print(f"ERROR: failed to load matrix: {e}", file=sys.stderr)
        return 3

    hits = find_matching_var_indices(
        adata,
        args.gene,
        case_insensitive=args.case_insensitive,
        gene_expression_only=args.gene_expression_only,
    )
    stats = _summarize_expression(adata, hits)
    in_ref = stats["n_feature_rows_matched"] > 0
    expressed = in_ref and stats["total_umis"] >= args.min_total_umis
    inferred_sample = _infer_sample_id_from_matrix(args.matrix)

    report = {
        "matrix": str(args.matrix.resolve()),
        "sample_id_inferred_for_output": inferred_sample,
        "query": args.gene,
        "n_observations": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "in_reference": in_ref,
        "expressed": expressed,
        "case_insensitive": bool(args.case_insensitive),
        "gene_expression_only": bool(args.gene_expression_only),
        "min_total_umis": float(args.min_total_umis),
        **stats,
        "matched_var_names": [str(adata.var_names[j]) for j in hits],
    }

    out_json = args.out_json
    if out_json is None:
        out_json = _default_json_out_path(args.matrix, args.gene)
    out_json = out_json.expanduser().resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    report["json_output_path"] = str(out_json)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    human = (
        f"[visium_hd_gene_qc] gene={args.gene!r} | in_reference={in_ref} | "
        f"expressed={expressed} (total_umis={stats['total_umis']:.0f}, "
        f"bins_with_signal={stats['n_bins_with_signal']}/{adata.n_obs})"
    )

    if args.json_only:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    elif args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
        print(human, file=sys.stderr)
    else:
        print(human)
        if in_ref and hits:
            print(f"  matched var_names: {report['matched_var_names']}")
    print(f"[visium_hd_gene_qc] wrote JSON: {out_json}", file=sys.stderr)

    if args.strict:
        if not in_ref:
            return 2
        if not expressed:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
