#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write manifest.yaml for a run-scoped MAGPIE output directory.

Intended to be invoked from Snakemake (rule write_run_manifest) when
``run.enabled`` is true. Records sample, resolved run subfolder, config
snapshots, and ``selected_msi_peaks`` path (MSI peak list used with Shiny / prep).
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import yaml


def _exists_meta(path: Path) -> dict:
    if not path or not str(path).strip():
        return {"path": "", "exists": False}
    p = Path(path)
    ok = p.is_file()
    out = {"path": str(p).replace("\\", "/"), "exists": ok}
    if ok:
        out["size_bytes"] = p.stat().st_size
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Write run manifest.yaml")
    ap.add_argument("--config", required=True, help="Path to config.yaml used for this run")
    ap.add_argument("--out", required=True, help="Output manifest path")
    ap.add_argument("--sample", required=True, help="Sample id")
    ap.add_argument(
        "--run-subdir-resolved",
        required=True,
        help="Resolved run subfolder (e.g. TPS/002um_run0 or with timestamp suffix)",
    )
    ap.add_argument(
        "--bin-size",
        default="008",
        help="HD bin size string (002|008|016) for key_outputs hints",
    )
    args = ap.parse_args()

    cfg_path = Path(args.config)
    with open(cfg_path, "r", encoding="utf-8") as f:
        full_config = yaml.safe_load(f) or {}

    msi = full_config.get("msi") or {}
    peaks_path = msi.get("selected_peaks") or msi.get("selected_msi_peaks")
    if not peaks_path:
        peaks_path = f"input/{args.sample}/msi/MSI_selected_peaks.txt"

    run_cfg = full_config.get("run") or {}
    out_manifest = Path(args.out)
    output_root = str(out_manifest.parent).replace("\\", "/")
    subdir_raw = run_cfg.get("subdir")
    subdir_explicit = bool(str(subdir_raw or "").strip())

    manifest = {
        "schema": "magpie_run_manifest",
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sample": args.sample,
        "run": {
            "enabled": bool(run_cfg.get("enabled")),
            "subdir_config": subdir_raw,
            "subdir_auto_built": bool(run_cfg.get("enabled")) and not subdir_explicit,
            "subdir_resolved": args.run_subdir_resolved,
            "append_timestamp": bool(run_cfg.get("append_timestamp", False)),
            "run_index": run_cfg.get("run_index", run_cfg.get("index", 0)),
        },
        "output_root": output_root,
        "selected_msi_peaks": _exists_meta(Path(peaks_path)),
        "config_file": str(cfg_path).replace("\\", "/"),
        "visium_mode": full_config.get("visium_mode"),
        "coregistration": full_config.get("coregistration"),
        "visium": full_config.get("visium"),
        "msi": {
            "metadata_csv": msi.get("metadata_csv"),
            "intensities_csv": msi.get("intensities_csv"),
            "radius_um": msi.get("radius_um"),
        },
        "aggregation": full_config.get("aggregation"),
        "spatialdata": full_config.get("spatialdata"),
        "key_outputs_relative_to_output_root": [
            "transformed.csv",
            "transformed.png",
            f"{args.bin_size}um/integrated_multimodal.h5ad",
            f"{args.bin_size}um/spatialdata.zarr",
            f"{args.bin_size}um/spaceranger_hd/",
        ],
    }

    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(out_manifest, "w", encoding="utf-8") as f:
        yaml.safe_dump(
            manifest,
            f,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )
    print(f"[write_run_manifest] wrote {out_manifest}")


if __name__ == "__main__":
    main()
