configfile: "config/config.yaml"

import os
import re
from datetime import datetime

def sample(): return config["visium"]["sample"]
def sr_dir(): return config["visium"]["spaceranger_dir"]
def visium_mode(): return config.get("visium_mode","classic")
def hd_bin_size(): return config.get("visium",{}).get("hd",{}).get("bin_size","008")

# ---- Run-scoped outputs (optional): output/<sample>/runs/<subdir>/... -----------------
def run_enabled() -> bool:
    return bool(config.get("run", {}).get("enabled", False))

def _run_cfg_slug(s: str) -> str:
    t = str(s or "").strip().lower()
    t = re.sub(r"[^a-z0-9]+", "_", t)
    return t.strip("_") or "none"

def run_coreg_slug() -> str:
    """Filesystem-safe tag from the three coregistration transform choices."""
    c = config.get("coregistration") or {}
    parts = [
        _run_cfg_slug(c.get("no_HE_transform", "affine")),
        _run_cfg_slug(c.get("MSI2HE_transform", "affine")),
        _run_cfg_slug(c.get("HE2HE_transform", "affine")),
    ]
    return "_".join(parts)

def run_subdir_auto() -> str:
    """
    Default runs/ path when run.subdir is unset: encode transforms, bin size, and run index.
    Example: affine_affine_tps/002um_run0
    """
    run_cfg = config.get("run") or {}
    idx = run_cfg.get("run_index", run_cfg.get("index", 0))
    try:
        idx = int(idx)
    except (TypeError, ValueError):
        idx = 0
    bin_s = str(hd_bin_size()).strip()
    return f"{run_coreg_slug()}/{bin_s}um_run{idx}"

def run_subdir_resolved() -> str:
    """Effective run subfolder; optional timestamp suffix when run.append_timestamp is true."""
    if not run_enabled():
        return ""
    run_cfg = config.get("run") or {}
    base = str(run_cfg.get("subdir", "")).strip().strip("/")
    if not base:
        base = run_subdir_auto()
    if run_cfg.get("append_timestamp", False):
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"{base}_{ts}"
    return base

RUN_RESOLVED_SUBDIR = run_subdir_resolved()

if run_enabled():
    if visium_mode() != "hd":
        raise ValueError(
            "config run.enabled is True but visium_mode must be 'hd' "
            "(run-scoped paths are implemented for the HD workflow only)."
        )
    if not RUN_RESOLVED_SUBDIR:
        raise ValueError("run.enabled is True but resolved run subdir is empty (internal error).")

def sample_out_base(wildcards) -> str:
    """Per-sample output prefix; under runs/<subdir>/ when run.enabled."""
    smp = wildcards.sample
    if run_enabled() and RUN_RESOLVED_SUBDIR:
        return f"output/{smp}/runs/{RUN_RESOLVED_SUBDIR}"
    return f"output/{smp}"

def sample_out_base_static() -> str:
    """Same as sample_out_base for the configured visium.sample (no wildcards)."""
    smp = sample()
    if run_enabled() and RUN_RESOLVED_SUBDIR:
        return f"output/{smp}/runs/{RUN_RESOLVED_SUBDIR}"
    return f"output/{smp}"

def sample_output_prefix() -> str:
    """Path prefix for rules that use ``{sample}`` wildcards (Snakemake disallows lambdas in ``output``)."""
    if run_enabled() and RUN_RESOLVED_SUBDIR:
        return f"output/{{sample}}/runs/{RUN_RESOLVED_SUBDIR}"
    return "output/{sample}"

HD_BIN_FOR_PATHS = hd_bin_size()
SAMPLE_OUT_PREFIX = sample_output_prefix()

def visium_hd_index_cli_flags(wildcards):
    """Flags for create_visium_hd_index.py: strict HD unless config opts into classic substitute."""
    hm = visium_mode()
    hd = config.get("visium", {}).get("hd", {})
    require_hd = hd.get("require_hd", True)
    allow_fb = hd.get("allow_classic_fallback", False)
    parts = [f"--visium_mode {hm}"]
    if hm == "hd" and (not require_hd or allow_fb):
        parts.append("--allow_classic_fallback")
    return " ".join(parts)

include_barcode_matrix = True

if os.path.isfile('input/selected.txt'):
    file = open("input/selected.txt", "r")
    samples = [line.rstrip() for line in file]
else:
    samples = glob_wildcards("input/{sample}/msi").sample
if os.path.isfile('input/exclude.txt'):
    file = open("input/exclude.txt", "r")
    samples = list(set(samples).difference(set([line.rstrip() for line in file])))

if include_barcode_matrix:
    output_file = "spaceranger_aggregated"
else:
    output_file = "spaceranger"

# Default build targets:
# - HD mode: produce Visium HD MEX, aggregated multimodal data, and spatialdata export
# - Classic mode: produce standard spaceranger(_aggregated) HDF5
if visium_mode() == "hd":
    bin_size = hd_bin_size()
    _base = sample_out_base_static()
    default_targets = [
        f"{_base}/{bin_size}um/spaceranger_hd/filtered_feature_bc_matrix/matrix.mtx.gz",
        f"{_base}/{bin_size}um/integrated_multimodal.h5ad",
        f"{_base}/{bin_size}um/spatialdata.zarr",
    ]
    if run_enabled() and RUN_RESOLVED_SUBDIR:
        default_targets.append(f"{_base}/manifest.yaml")
    ruleorder: visium_hd_index > perform_coreg
else:
    default_targets = expand("output/{sample}/"+output_file+"/filtered_feature_bc_matrix.h5", sample=samples)

rule all:
    input:
        default_targets

rule visium_hd_index:
    input:
        spaceranger_dir = lambda w: sr_dir()
    output:
        bins_csv = f"output/{{sample}}/visium_hd_square_{hd_bin_size()}um_bins.csv"
    params:
        bin_size = lambda w: hd_bin_size(),
        sample = "{sample}",
        index_flags = lambda w: visium_hd_index_cli_flags(w),
    shell:
        """
        python scripts/create_visium_hd_index.py \
        --spaceranger_dir {input.spaceranger_dir} \
        --sample {params.sample} \
        --bin_size {params.bin_size} \
        {params.index_flags} \
        --out_csv {output.bins_csv}
        """
   
rule check_inputs:
    message:
        'Checking inputs.'
    conda: 'magpie'
    input:
        "input/"
    output:
        "output/summary.csv"
    params:
        verbose = True
    script:
        "scripts/check_inputs.py"

rule perform_coreg:
    message:
        "Performing co-registration."
    conda: 'magpie'
    input:
        "input/{sample}/msi/MSI_metadata.csv",
        "output/summary.csv"
    output:
        transformed_csv = f"{SAMPLE_OUT_PREFIX}/transformed.csv",
        transformed_png = f"{SAMPLE_OUT_PREFIX}/transformed.png",
    params:
        no_HE_transform = lambda w: config["coregistration"]["no_HE_transform"],
        MSI2HE_transform = lambda w: config["coregistration"]["MSI2HE_transform"],
        HE2HE_transform = lambda w: config["coregistration"]["HE2HE_transform"],
        sample = "{sample}",
        verbose = True
    script:
        "scripts/alter_data.py"

rule collapse_to_visiumhd_bins:
    input:
        msi_meta = config["msi"]["metadata_csv"],
        msi_int = config["msi"]["intensities_csv"],
        transformed_csv = lambda w: f"{sample_out_base(w)}/transformed.csv",
        visium_hd_bins = f"output/{{sample}}/visium_hd_square_{hd_bin_size()}um_bins.csv"
    output:
        spaceranger_hd = directory(f"{SAMPLE_OUT_PREFIX}/{HD_BIN_FOR_PATHS}um/spaceranger_hd"),
        matrix = f"{SAMPLE_OUT_PREFIX}/{HD_BIN_FOR_PATHS}um/spaceranger_hd/filtered_feature_bc_matrix/matrix.mtx.gz",
        tpos = f"{SAMPLE_OUT_PREFIX}/{HD_BIN_FOR_PATHS}um/spaceranger_hd/spatial/tissue_positions_list.csv"
    params:
        spaceranger_dir = lambda w: sr_dir(),
        bin_size = lambda w: hd_bin_size(),
        out_root = lambda w: f"{sample_out_base(w)}/{hd_bin_size()}um",
        sample = "{sample}"
    shell:
        """
        python scripts/collapse_to_visiumhd_bins.py \
        --sample {params.sample} \
        --bin_size {params.bin_size} \
        --msi_meta {input.msi_meta} \
        --msi_int {input.msi_int} \
        --transformed_csv {input.transformed_csv} \
        --visium_hd_bins {input.visium_hd_bins} \
        --spaceranger_dir {params.spaceranger_dir} \
        --out_root {params.out_root}
        """

rule aggregate_to_msi_resolution:
    input:
        transformed_csv = lambda w: f"{sample_out_base(w)}/transformed.csv",
        visium_hd_bins = f"output/{{sample}}/visium_hd_square_{hd_bin_size()}um_bins.csv",
        visium_hd_matrix = f"input/{{sample}}/visium/{hd_bin_size()}um/filtered_feature_bc_matrix.h5",
        msi_int = config["msi"]["intensities_csv"]
    output:
        mapping = f"{SAMPLE_OUT_PREFIX}/{HD_BIN_FOR_PATHS}um/hd_bins_per_msi_spot.csv",
        integrated = f"{SAMPLE_OUT_PREFIX}/{HD_BIN_FOR_PATHS}um/integrated_multimodal.h5ad"
    params:
        bin_size = lambda w: hd_bin_size(),
        out_dir = lambda w: f"{sample_out_base(w)}/{hd_bin_size()}um",
        spaceranger_dir = lambda w: sr_dir(),
        msi_radius_um = config.get("msi", {}).get("radius_um", 10.0),
        agg_method = config.get("aggregation", {}).get("method", "sum"),
        sample = "{sample}"
    shell:
        """
        python scripts/aggregate_to_msi_resolution.py \
        --sample {params.sample} \
        --bin_size {params.bin_size} \
        --transformed_csv {input.transformed_csv} \
        --visium_hd_bins {input.visium_hd_bins} \
        --visium_hd_matrix {input.visium_hd_matrix} \
        --msi_intensities {input.msi_int} \
        --spaceranger_dir {params.spaceranger_dir} \
        --msi_radius_um {params.msi_radius_um} \
        --agg_method {params.agg_method} \
        --out_dir {params.out_dir}
        """

rule export_spatialdata_zarr:
    conda: 'magpie'
    input:
        spaceranger_hd = f"{sample_out_base_static()}/{hd_bin_size()}um/spaceranger_hd/filtered_feature_bc_matrix/matrix.mtx.gz",
        bins_csv = f"output/{config['visium']['sample']}/visium_hd_square_{hd_bin_size()}um_bins.csv"
    output:
        zarr = directory(f"{sample_out_base_static()}/{hd_bin_size()}um/spatialdata.zarr")
    params:
        sample = config['visium']['sample'],
        bin_size = hd_bin_size(),
        sr_hd = f"{sample_out_base_static()}/{hd_bin_size()}um/spaceranger_hd",
        bins_csv = f"output/{config['visium']['sample']}/visium_hd_square_{hd_bin_size()}um_bins.csv",
        out_zarr = f"{sample_out_base_static()}/{hd_bin_size()}um/spatialdata.zarr",
        features = config.get("spatialdata", {}).get("features", ""),
        he_override = config.get("spatialdata", {}).get("he_override", ""),
        sigma = config.get("spatialdata", {}).get("gaussian_sigma", 0.0)
    shell:
        r"""
        python scripts/export_spatialdata.py \
        --sample {params.sample} \
        --bin_size {params.bin_size} \
        --spaceranger_hd_dir {params.sr_hd} \
        --visium_hd_bins {params.bins_csv} \
        --out_zarr {params.out_zarr} \
        --features "{params.features}" \
        --he_override "{params.he_override}" \
        --gaussian_sigma {params.sigma}
        """

if run_enabled() and RUN_RESOLVED_SUBDIR:
    rule write_run_manifest:
        conda: "magpie"
        input:
            transformed = f"{sample_out_base_static()}/transformed.csv",
            integrated = f"{sample_out_base_static()}/{hd_bin_size()}um/integrated_multimodal.h5ad",
        output:
            manifest = f"{sample_out_base_static()}/manifest.yaml",
        params:
            sample = sample(),
            run_sd = RUN_RESOLVED_SUBDIR,
            bin_sz = hd_bin_size(),
        shell:
            """
            python scripts/write_run_manifest.py \
                --config config/config.yaml \
                --out {output.manifest} \
                --sample {params.sample} \
                --run-subdir-resolved "{params.run_sd}" \
                --bin-size {params.bin_sz}
            """

rule make_spaceranger:
    message:
        "Creating spaceranger formatted output."
    conda: 'magpie'
    input:
        "output/{sample}/transformed.csv",
        "output/{sample}/transformed.png"
    output:
        "output/{sample}/spaceranger/filtered_feature_bc_matrix.h5"
    params:
        sample = "{sample}",
        verbose = True
    script:
        "scripts/create_mock_spaceranger.py"

rule create_barcode_matrix:
    message:
        "Generating aggregated data."
    conda: 'magpie'
    input:
        "output/{sample}/spaceranger/filtered_feature_bc_matrix.h5"
    output:
        "output/{sample}/spaceranger_aggregated/filtered_feature_bc_matrix.h5"
    params:
        sample = "{sample}",
        radius_to_use = 'visium_expanded',
        agg_fn = 'mean',
        verbose = True,
        only_within_tissue = False
    script:
        "scripts/create_perbarcode_matrix.py"
