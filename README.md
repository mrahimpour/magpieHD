# MAGPIE: Multimodal alignment of genes and peaks for integrative exploration
<p align="center">
<img src="figures/magpie_logo.png" width="200">
</p>

Recent developments in spatially resolved -omics have enabled studies linking gene expression and metabolite levels to tissue morphology, offering new insights into biological pathways. By capturing multiple modalities on matched tissue sections, one can better probe how different biological entities interact in a spatially coordinated manner. However, such cross-modality integration presents experimental and computational challenges. 
To align multimodal datasets into a shared coordinate system and facilitate enhanced integration and analysis, we propose *MAGPIE* (Multi-modal Alignment of Genes and Peaks for Integrated Exploration), a framework for co-registering spatially resolved transcriptomics, metabolomics, and tissue morphology from the same or consecutive sections. 
We illustrate the generalisability and scalability of *MAGPIE* on spatial multi-omics data from multiple tissues, combining Visium with both MALDI and DESI mass spectrometry imaging. *MAGPIE* was also applied to newly generated multimodal datasets created using specialised experimental sampling strategy to characterise the metabolic and transcriptomic landscape in an in vivo model of drug-induced pulmonary fibrosis, to showcase the linking of small-molecule co-detection with endogenous responses in lung tissue.
*MAGPIE* highlights the refined resolution and increased interpretability of spatial multimodal analyses in studying tissue injury, particularly in pharmacological contexts, and offers a modular, accessible computational workflow for data integration.

Preprint: https://www.biorxiv.org/content/10.1101/2025.02.26.640381v1

## Installation

The MAGPIE pipeline requires a Python installation and the following package dependencies:
* snakemake
* shiny
* matplotlib
* pandas
* numpy
* scikit-image
* pathlib
* scikit-learn
* scipy
* json
* collections
* shutil
* gzip
* h5py
* scanpy

We recommend to create a conda environment with from which the whole pipeline can be run. You can install all required dependencies using the magpie_environment.yml file within the snakemake folder in this repository using the following command:
```
conda env create -f magpie_environment.yml
```

The pipeline has been previously tested on the following systems:
* macOS: Sequoia (15.3.2)
* Windows: 11 (22H2)

Installation should take up to ~10 minutes on a normal desktop computer.

## Input structure

The MAGPIE pipeline automatically detects the files in your input folder and makes decisions accordingly so you must ensure your files follow the following structure:

    [sample name]
    ├── visium/                              # Spaceranger outputs
    │   ├── spatial/                         # Shared spatial data
    │   │   ├── aligned_fiducials.jpg
    │   │   ├── detected_tissue_image.jpg
    │   │   ├── scalefactors_json.json
    │   │   ├── tissue_hires_image.png
    │   │   ├── tissue_lowres_image.png
    │   │   ├── tissue_square_002um_positions.parquet
    │   │   ├── tissue_square_008um_positions.parquet
    │   │   └── tissue_square_016um_positions.parquet
    │   ├── 002um/                           # 2µm resolution gene expression
    │   │   └── filtered_feature_bc_matrix.h5
    │   ├── 008um/                           # 8µm resolution gene expression
    │   │   └── filtered_feature_bc_matrix.h5
    │   └── 016um/                           # 16µm resolution gene expression
    │       └── filtered_feature_bc_matrix.h5
    ├── msi/                    
    │   ├── MSI_intensities.csv              # Table of intensities with MSI peaks on columns and pixels on rows
    │   ├── MSI_metadata.csv                 # Table of metadata about MSI pixels, including x and y coordinate columns
    │   └── MSI_HE.[jpg,png,tiff]            # (OPTIONAL) intermediate MSI image to assist with coregistration
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Landmarks between MSI image and MSI H&E
    ├── landmarks_HE2HE.csv                  # (OPTIONAL) Landmarks between MSI H&E and Visium H&E
    └── landmarks_noHE.csv                   # (OPTIONAL) Landmarks between MSI image and Visium H&E (direct)
    


## Running the shiny app

To run the pipeline, you need to be in the folder with all files in the _snakemake_ folder in this repository as well as an _input_ folder as described in the previous section.

To start the shiny app for manual landmark selection, run ``` shiny run magpie_shiny_app.py ```

For each sample you will be prompted to select some manual landmarks then download. At the point you download them they will be saved into the file structure described above. If you would prefer to use your own landmarks please save them into that structure instead and you can skip the shiny app step.

## Running the snakemake pipeline

Once landmarks have been selected for each sample, you can switch to the snakemake pipeline to perform the coregistration. Again you must be in the folder with all files in the _snakemake_ folder in this repository as well as an _input_ folder as described in the previous section with your newly selected landmarks. You can then run the pipeline using ``` snakemake --cores [n] ``` where _n_ is the number of cores you would like to use. You can explicitly state which samples you would like to use by listing them in a *selected.txt* file within the *input* folder and equivalently specify some files you would like to exclude using a *exclude.txt* file.

## QC: Visium HD gene check

Optional script: check whether a gene appears in the Space Ranger **filtered** HD matrix and whether it has non-zero UMIs. Run from the repository root with the same Python environment you use for scanpy (e.g. your conda env):

```bash
python scripts/QC/visium_hd_gene_qc.py \
  --matrix input/<sample>/visium/002um/filtered_feature_bc_matrix.h5 \
  --gene P2ry6 \
  --case-insensitive
```

A JSON report is **always** written next to that sample’s outputs: default `output/<sample>/qc/visium_hd_gene_qc_<gene>.json` when `--matrix` looks like `input/<sample>/visium/...`, otherwise `output/unscoped/qc/...`. Use `--out-json` to set the path explicitly. Use `--help` for all options (`--json`, `--strict`, etc.).

Minimal spatial plot (one gene on H&E): edit the CONFIG block in `scripts/QC/plot_one_gene_visium_hd.py`, then `python scripts/QC/plot_one_gene_visium_hd.py` from the repo root.

## Multi-modal Integration (Visium HD)

After registration, use the aggregation script to create an integrated multi-modal dataset:

```bash
python scripts/aggregate_to_msi_resolution.py \
    --sample YOUR_SAMPLE \
    --bin_size 008 \
    --msi_radius_um 10.0 \
    --agg_method sum
```

**Parameters:**
- `--bin_size`: HD resolution — `002` (2µm), `008` (8µm), or `016` (16µm)
- `--msi_radius_um`: Search radius in µm on bin centers (default 10 ≈ half of ~20µm MSI pixel)

See `QuickReview.md` for detailed workflow documentation.

## Output Structure

The pipeline produces two types of outputs:

### Resolution-Independent (Shared)
These files are generated once during registration and apply to all HD resolutions:

| File | Description |
|------|-------------|
| `transformed.csv` | MSI coordinates transformed to HD pixel space |
| `transformed.png` | Registration visualization |
| `transformed_withCoords_VisiumHE.png` | Registration overlay on H&E |
| `visium_hd_square_XXXum_bins.csv` | HD bin positions for each resolution |

### Resolution-Specific
These files are generated per resolution (002um, 008um, 016um):

| File | Description |
|------|-------------|
| `hd_bins_per_msi_spot.csv` | HD bins within each MSI spot's footprint (1:N mapping) |
| `aggregated_gene_expression.h5ad` | Gene expression aggregated to MSI resolution |
| `integrated_multimodal.h5ad` | **Combined dataset with genes + metabolites** |
| `spaceranger_hd/` | MSI in Space Ranger MEX format |
| `spatialdata.zarr/` | SpatialData export for Squidpy |

```
output/{sample}/
├── transformed.csv                      # Shared
├── visium_hd_square_008um_bins.csv      # Shared
├── 002um/                               # 2µm outputs
│   ├── integrated_multimodal.h5ad
│   └── spaceranger_hd/
├── 008um/                               # 8µm outputs
└── 016um/                               # 16µm outputs
```

## Tutorial

We provide extensive [documentation describing the pipeline](https://core-bioinformatics.github.io/magpie/) and a tutorial with example data described [here](https://core-bioinformatics.github.io/magpie/tutorial/SMA_tutorial.html). The tutorial should take around 5-10 minutes to run. 
