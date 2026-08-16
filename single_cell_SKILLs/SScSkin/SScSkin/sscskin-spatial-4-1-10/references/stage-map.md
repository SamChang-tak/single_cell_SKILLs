# Spatial stages 01–10

The numeric range contains 11 notebooks because stage 07 has two independent branches.

| Stage | Notebook | Language | Purpose and output |
|---|---|---|---|
| 01 | `01_read_data.ipynb` | R | Load each 10x Visium sample with Seurat, retain spots with `nFeature_Spatial >= 100`, and save per-sample RDS files. |
| 02 | `02_viz_marker.ipynb` | R | Display spatial marker expression across healthy and SSc samples. |
| 03 | `03_run_celltrek.ipynb` | R | Transfer the integrated scRNA-seq reference into each Visium section with CellTrek and save mapped-cell RDS files. |
| 04 | `04_viz_celltrek.ipynb` | R | Visualize mapped CellTrek cell types with sample-specific point-size choices. |
| 05 | `05_viz_marker.ipynb` | R | Produce additional spatial marker panels. |
| 06 | `06_add_cell2location.ipynb` | R | Attach precomputed cell2location CSV abundance estimates to Seurat metadata and visualize predicted types/abundances. |
| 07A | `07_banksy.ipynb` | R | Convert counts/coordinates to SpatialExperiment and run BANKSY (`k_geom=50`, lambda `0.2`, 30 PCs, UMAP, Leiden). |
| 07B | `07_cell_neighborhood.ipynb` | Python | Build per-section Squidpy Delaunay graphs and calculate neighborhood enrichment using `pred_cell_type`. |
| 08 | `08_viz_cell_neighborhood.ipynb` | R | Visualize sample-level neighborhood-enrichment CSV matrices with heatmaps. |
| 09 | `09_run_cell2location.ipynb` | Python | Train cell2location reference signatures and per-sample spatial models; write H5AD and CSV outputs. |
| 10 | `10_viz_cell2location.ipynb` | R | Attach and visualize the stage-09 cell2location outputs in Seurat. |

The notebooks reference 16 samples: `HC01`, `HC02`, `HC03`, `HC05`, `SSc4733`,
`SSc4994`, `SSc5380`, `SSc5722`, `SSc-HL01`, `SSc-HL05`, `SSc-HL06`, `SSc-HL11`,
`SSc-HL13`, `SSc-HL25`, `SSc-HL33`, and `SSc-HL35`. Confirm actual availability because
some downstream stages omit selected samples.

CellTrek upstream parameters include `dist_thresh=0.55`, `top_spot=5`, `spot_n=5`,
`repel_r=20`, and `repel_iter=20`. Treat them as study-specific settings and assess
sensitivity where mapping conclusions depend on them.
