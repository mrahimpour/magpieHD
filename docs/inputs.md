(input-structure)=
# MAGPIE directory structure

The directory from which you run the pipeline should have the following general structure:

    ├── Snakefile
    ├── magpie_shiny_app.py
    ├── figures  
    │   ├── magpie_logo.png
    ├── scripts
    │   ├── alter_data.py
    │   ├── create_mock_spaceranger.py
    │   ├── create_perbarcode_matrix.py
    ├── input
    │   ├── (optional) exclude.txt
    │   ├── (optional) selected.txt
    │   ├── ... 

The exclude.txt and selected.txt files are optional and allow the user to list a number of samples to include (rather than using all sub-folders in the input directory by default) or to exclude from analysis. 

The Snakefile, magpie_shiny_app.py and scripts and figures folders should all be copied from the GitHub repository.

(inputstructure)=
## Input folder structure

The MAGPIE pipeline automatically detects the files in your input folder and makes decisions accordingly so you must ensure your files follow the following structure:

    [sample name]
    ├── visium/                              # Spaceranger HD outputs
    │   ├── spatial/                         # Shared spatial data (images, positions)
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
    │   ├── MSI_intensities.csv              # Intensities (peaks on columns, pixels on rows)
    │   ├── MSI_metadata.csv                 # Metadata with spot_id, x, y columns
    │   └── MSI_HE.[jpg,png,tiff]            # (OPTIONAL) Intermediate MSI H&E image
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Landmarks: MSI image → MSI H&E
    ├── landmarks_HE2HE.csv                  # (OPTIONAL) Landmarks: MSI H&E → Visium H&E
    └── landmarks_noHE.csv                   # (OPTIONAL) Landmarks: MSI image → Visium H&E (direct)

**Note:** For Visium HD, gene expression files are organized by resolution in separate folders (`002um/`, `008um/`, `016um/`). 
Set `bin_size` in `config/config.yaml` to select which resolution to use.

The `visium_hd_index` rule runs `scripts/create_visium_hd_index.py` with `--visium_mode` from top-level `visium_mode` (`hd` | `classic`). In `hd` mode the script expects a matching HD bin positions file for `visium.hd.bin_size` and **fails** if none is found, unless you opt in via `visium.hd.allow_classic_fallback: true` or `visium.hd.require_hd: false` (both cause `--allow_classic_fallback` to be passed). Classic spot geometry is not interchangeable with true HD bins—use those flags only for debugging or migration.

`aggregate_to_msi_resolution` reads `visium.spaceranger_dir` for `scalefactors_json.json` (uses **`microns_per_pixel`** for fullres µm→px; HD JSON may also include `bin_size_um` / `spot_diameter_fullres` for cross-checks), the HD bin CSV, Space Ranger `filtered_feature_bc_matrix.h5` at the chosen resolution, and `msi.intensities_csv` + `transformed.csv`. Observation order follows **MSI intensities** (every `spot_id` there must appear in `transformed.csv`). Set **`msi.radius_um`** as the search radius (µm) on bin centers—**10 µm** is a typical default when MSI pixels are **~20 µm** wide (≈ half-width from center). Genes are pooled with `aggregation.method`: `sum` / `mean` / `median`.

## Run-scoped outputs and manifest (Visium HD only)

Set `run.enabled: true` in `config/config.yaml` to write **registration and HD integration** artifacts under:

`output/<sample>/runs/<subdir>/` (e.g. `runs/TPS/002um_run0`).

The HD bin index CSV stays shared at `output/<sample>/visium_hd_square_<bin>um_bins.csv` (Space Ranger geometry only).

- **`run.subdir`**: logical folder under `runs/` (can include slashes).
- **`run.append_timestamp`**: if `true`, Snakemake appends `_YYYYMMDD_HHMMSS` to `subdir` so each invocation gets a **new** directory.
- **`manifest.yaml`**: written in the run root when `run.enabled` is true; includes **`selected_msi_peaks`** (path, `exists`, size) from `msi.selected_peaks` in config (default `input/<sample>/msi/MSI_selected_peaks.txt`), plus registration and pipeline settings.

`run.enabled` requires `visium_mode: hd`. The classic Visium branch is unchanged and does not use run-scoped paths.

## Image Coordinate Spaces

Visium HD uses a multi-resolution image pyramid. Understanding these spaces is important for coordinate alignment:

| Space | Image File | Typical Size | Scale Factor | Usage |
|-------|-----------|--------------|--------------|-------|
| **FULLRES** | (Original scan) | ~11,000 × 12,000 px | 1.0 (reference) | HD bin positions in parquet files |
| **HIRES** | `tissue_hires_image.png` | ~5,000 × 6,000 px | ~0.48 | Registration, visualization |
| **LOWRES** | `tissue_lowres_image.png` | ~500 × 600 px | ~0.05 | Thumbnails |

**Key points:**
- `scalefactors_json.json` contains `tissue_hires_scalef` to convert FULLRES → HIRES
- HD bin positions (`pxl_col_in_fullres`, `pxl_row_in_fullres`) are in **FULLRES** coordinates
- Registration uses `tissue_hires_image.png`, so transformed MSI coordinates are in **HIRES** space
- MAGPIE automatically converts between spaces during mapping

**Important:** `scalefactors_json.json` is **shared** across all HD resolutions (002um, 008um, 016um) because it describes the image scaling, not the bin sizes. The bin-specific information is in each parquet file.

## MSI data

The MSI data should be split into an intensity table and a metadata table. The intensity table should look like this, with peaks on the columns and pixels on the rows. 

|spot_id|mz-101.2345   | mz-102.85754 | mz-303.35855   | mz-344.48575  | mz-321.38583  | mz-112.28485 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---: |
|pixel_0| 5774.812  | 675.361  | 23.555  | 8444.958  | 777.234  | 20.332  |
|pixel_1| 8794.013  | 444.523  | 81.294  | 6775.393  | 899.284  | 10.275  |
|pixel_2| 6777.358  | 857.585  | 14.326  | 9468.367  | 747.385  | 24.521  |
|| ...  | ...  | ...  | ...  | ... | ... |

The metadata table must contains columns spot_id, x and y (and can additionally include others) and the spot_id columns must exactly match the spot_id column in the intensity table.

| spot_id           | x   | y   | Sample   | Treatment | Fibrotic |
| :---:        | :---:       | :---:       | :---:       | :---:   | :---: |
| pixel_0 | 1        | 1        | sample_1        | treatment_1 | normal |
| pixel_1 | 2        | 1        | sample_1        | treatment_1 | normal |
| pixel_2 | 2        | 2        | sample_1       | treatment_1 | fibrotic |
|| ...  | ...  | ...  | 

The first stage of the snakemake pipeline will check all these criteria and save a summary table called summary.txt. 
