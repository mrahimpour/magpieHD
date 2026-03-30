# Pipeline Overview

The pipeline consists of the following steps:

1. **Landmark selection**: Manual landmarks are selected using the MAGPIE shiny app
2. **Coregistration**: Coordinates from modalities are set into the same coordinate system
3. **Creation of spaceranger-style object**: A new spaceranger-style object is created for MSI data using coregistered coordinates
4. **Aligning observations across modalities**: Observations are matched between modalities so that multi-omics methods requiring common samples can be applied.
5. **Gene expression aggregation** (Visium HD): `aggregate_to_msi_resolution.py` finds HD bin **centers** within a fixed radius (`msi.radius_um` in `config.yaml`, µm → fullres pixels via `scalefactors_json.json`) around each coregistered MSI point, then pools Visium gene counts with `aggregation.method` (`sum` / `mean` / `median`). This is **RNA → MSI grid** (contrast with `collapse_to_visiumhd_bins.py`, which places MSI on the HD barcode grid for MEX / SpatialData-style use).
6. **Integrated dataset creation**: Writes `integrated_multimodal.h5ad` (genes in `X`, raw MSI matrix in `obsm['msi']`).

Here is a summary of the MAGPIE pipeline:

![pipeline](figures/pipeline_diagram.jpg)

## Output Files

### Resolution-Independent (Shared)

These files are generated **once** during registration and apply to all HD resolutions:

| File | Description |
|------|-------------|
| `transformed.csv` | MSI coordinates transformed to HD pixel space |
| `transformed.png` | Registration visualization |
| `visium_hd_square_XXXum_bins.csv` | HD bin positions for each resolution |

**Why shared?** Registration aligns MSI to the H&E image, which is the same regardless of HD bin resolution.

### Resolution-Specific (Per Resolution)

These files are generated **separately** for each HD resolution (002um, 008um, 016um):

| File | Description |
|------|-------------|
| `{resolution}/hd_bins_per_msi_spot.csv` | Long table: MSI ↔ each HD bin within radius (`msi.radius_um`) |
| `{resolution}/aggregated_gene_expression.h5ad` | Genes pooled onto MSI observations (`aggregation.method`) |
| `{resolution}/integrated_multimodal.h5ad` | Same + `obsm['msi']` (per-row metabolites; not split across bins) |
| `{resolution}/spaceranger_hd/` | MSI in Space Ranger MEX format |
| `{resolution}/spatialdata.zarr/` | SpatialData export |

**Why separate?** Aggregating gene expression depends on HD bin size—more bins aggregate into each MSI spot at higher resolutions.

## Key Concept: Resolution Matching

Since MSI (often ~10–20 µm pixel pitch) is **coarser** than Visium HD bins (2 / 8 / 16 µm), we **pool HD gene counts onto each MSI observation**: every HD bin **center** within ``msi.radius_um`` (µm, using fullres ``microns_per_pixel`` from ``scalefactors_json.json``) contributes to that MSI row (``aggregation.method``: sum / mean / median). Choose ``radius_um`` from **MSI ground sampling**, not from classic Visium spot size unless they match. This:
- Avoids pretending MSI is natively 2 µm
- Uses many small HD bins as support under each MSI pixel
- Produces one matched observation per MSI row for multimodal ``.h5ad``

## Registration Validation

The pipeline validates registration quality using several metrics:

| Metric | Description | Good Values |
|--------|-------------|-------------|
| **Distance (median)** | Median MSI-to-HD bin distance | < 20 px |
| **Distance (95th pct)** | 95th percentile distance | < 100 px for in-tissue spots |
| **Jacobian (median)** | Local area scaling | ~1.0 (no distortion) |
| **Jacobian (std)** | Variability in local distortion | < 0.3 (uniform transform) |
| **Coverage** | % MSI spots near HD tissue | > 80% |

### Understanding Jacobian Determinant

The Jacobian determinant measures how much local area changes during transformation:

- **|J| = 1.0**: No area change (ideal rigid transform)
- **|J| > 1.0**: Local expansion
- **|J| < 1.0**: Local compression
- **|J| constant**: Pure affine transform (no local warping)
- **High std(J)**: Significant local distortion (TPS transform)

Run validation:
```bash
python scripts/validate_registration.py --sample YOUR_SAMPLE --bin_size 008
```

Output files:
- `registration_metrics.json` - Quantitative metrics
- `registration_validation.png` - Visual summary
