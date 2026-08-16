# Spatial pipeline map

## Contents

- [Samples and entry points](#samples-and-entry-points)
- [Mapping and deconvolution](#mapping-and-deconvolution)
- [Neighborhoods and niches](#neighborhoods-and-niches)
- [Functions and communication](#functions-and-communication)

## Samples and entry points

The notebooks reference 16 samples: `HC01`, `HC02`, `HC03`, `HC05`, `SSc4733`,
`SSc4994`, `SSc5380`, `SSc5722`, `SSc-HL01`, `SSc-HL05`, `SSc-HL06`, `SSc-HL11`,
`SSc-HL13`, `SSc-HL25`, `SSc-HL33`, and `SSc-HL35`. Confirm availability rather than
assuming every downstream notebook uses all 16.

- `01_read_data.ipynb` (R): load each 10x Visium directory with Seurat
  `Load10X_Spatial`, retain spots with at least 100 spatial features, save per-sample RDS.
- `02_viz_marker.ipynb` and `05_viz_marker.ipynb` (R): visualize spatial marker panels.

## Mapping and deconvolution

- `03_run_celltrek.ipynb` (R): map integrated scRNA-seq cells to each Visium sample with
  CellTrek. Upstream settings include `dist_thresh=0.55`, `top_spot=5`, `spot_n=5`,
  `repel_r=20`, and `repel_iter=20`.
- `04_viz_celltrek.ipynb` (R): display mapped cell types with sample-specific point sizes.
- `06_add_cell2location.ipynb` (R): attach earlier cell2location CSV abundance estimates
  to Seurat metadata.
- `09_run_cell2location.ipynb` (Python): train the reference RegressionModel on the
  integrated scRNA-seq atlas, infer reference signatures, fit a Cell2location model per
  spatial sample, and save H5AD/CSV outputs.
- `10_viz_cell2location.ipynb` (R): attach and visualize the newly generated estimates.
- `07_banksy.ipynb` (R): build SpatialExperiment objects and run BANKSY with `k_geom=50`,
  lambda `0.2`, 30 PCs, UMAP, and Leiden clustering.

CellTrek and cell2location are alternative representations: mapped individual reference cells
versus posterior abundance estimates for mixed spots.

## Neighborhoods and niches

- `07_cell_neighborhood.ipynb` / `08_viz_cell_neighborhood.ipynb`: Squidpy Delaunay
  graphs and neighborhood enrichment for the earlier cell2location output.
- `11_cell_neighborhood.ipynb` / `12_viz_cell_neighborhood.ipynb`: repeat neighborhood
  analysis for the notebook-09 cell2location output.
- `13_cell_colocalization.ipynb`, `14_cell_colocalization.ipynb`, and
  `16_cell_colocalization.ipynb`: abundance-correlation/colocalization summaries and
  healthy-versus-SSc heatmaps.
- `15_integrate_with_cell_proportion.ipynb`: concatenate sample-level cell proportions,
  compute a 30-neighbor graph, Leiden niches at resolution 0.5, and niche composition.
- `17_heatmap.ipynb`: visualize niche composition with ComplexHeatmap.
- `18_compare_niches.ipynb`: compare patient-level niche proportions by condition.
- `19_viz_niches_in_spatial.ipynb`: map niche labels back to Seurat spatial objects.
- `22_molecular_niches.ipynb`: integrate molecular spot expression, run Harmony by
  patient, and cluster molecular niches at Leiden resolution 0.3.

## Functions and communication

- `20_function_analysis_of_niches.ipynb`: rank niche markers and score functional gene-set
  collections with decoupler.
- `21_viz_function_of_niches.ipynb`: rank and visualize pathway/source activities across
  niches for multiple collections.
- `23_cellchat.ipynb` / `24_viz_cellchat.ipynb`: run and display spatial CellChat on one
  set of cell2location-enriched Seurat objects.
- `25_cellchat.ipynb` / `26_viz_cellchat.ipynb`: repeat CellChat on the later output set.
- `27_decoupler.ipynb`: calculate per-sample pathway and transcription-factor activities.
