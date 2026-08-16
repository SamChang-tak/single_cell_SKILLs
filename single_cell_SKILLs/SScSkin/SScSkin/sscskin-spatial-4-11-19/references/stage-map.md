# Spatial stages 11–19

The inclusive numeric range contains nine notebooks.

| Stage | Notebook | Language | Purpose |
|---|---|---|---|
| 11 | `11_cell_neighborhood.ipynb` | Python | Read stage-09 cell2location H5AD files, build per-section Squidpy Delaunay graphs, and calculate enrichment using `pred_cell_type`. |
| 12 | `12_viz_cell_neighborhood.ipynb` | R | Visualize sample-level neighborhood-enrichment matrices with ComplexHeatmap. |
| 13 | `13_cell_colocalization.ipynb` | Python | Calculate per-sample correlations among inferred cell-type abundances and create clustered heatmaps/condition summaries. |
| 14 | `14_cell_colocalization.ipynb` | R | Display healthy and SSc colocalization summaries with ComplexHeatmap. |
| 15 | `15_integrate_with_cell_proportion.ipynb` | Python | Combine spot-level cell-proportion profiles, build a 30-neighbor graph, define Leiden niches at resolution 0.5, and export niche compositions and spot labels. |
| 16 | `16_cell_colocalization.ipynb` | R | Compare healthy and SSc colocalization matrices. |
| 17 | `17_heatmap.ipynb` | R | Visualize cell-type abundance/proportion across niches using multiple ComplexHeatmap layouts. |
| 18 | `18_compare_niches.ipynb` | Python | Calculate patient-normalized niche proportions and compare healthy versus SSc distributions. |
| 19 | `19_viz_niches_in_spatial.ipynb` | R | Join integrated niche labels back to each sample’s Seurat object and plot spatial niche maps. |

Inputs originate from stage 09 cell2location outputs. The notebooks reference healthy and
SSc samples unevenly; confirm the exact available sample set at each stage and explain every
omission.

Key upstream baseline settings:

- spatial graph: Squidpy Delaunay neighbors;
- neighborhood label: `pred_cell_type`;
- niche feature space: inferred cell-type abundance/proportion;
- integrated niche graph: 30 neighbors;
- niche clustering: Leiden resolution 0.5.
