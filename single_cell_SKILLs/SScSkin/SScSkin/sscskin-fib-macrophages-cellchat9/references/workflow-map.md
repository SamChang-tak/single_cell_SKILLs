# Workflow map

| Stage | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| 01 | Read fibroblast and macrophage H5AD objects, combine count matrices and metadata, create and normalize a Seurat object, run PCA/UMAP | fibroblast stage-01 H5AD; macrophage stage-02 H5AD | integrated `obj.rds` |
| 02 | Create CellChat grouped by `sub_clusters`, use human database, identify overexpressed features/interactions, compute triMean probabilities, filter groups below 10 cells, aggregate pathways | stage-01 Seurat RDS | CellChat `obj.rds` |
| 03 | Visualize interaction counts/weights, sender networks, CXCL pathway layouts, bubbles, and MIF/CXCL expression | stage-02 CellChat RDS | network and expression plots |
| 04 | Export database interactions, classify macrophage↔fibroblast edges, count ligand–receptor pairs, draw subtype heatmap/bubbles | stage-02 CellChat RDS; bundled `all_interaction.csv` | `LR_pair_counts.csv`, `LR_pair_counts_by_subtype.csv`, plots |

## Required upstream objects

- `05_fibroblast/01_clustering/fibroblast.h5ad`, with `layers['counts']` and fibroblast
  `sub_clusters`.
- `07_marchphages/02_annotate/adata.h5ad`, with `layers['counts']` and macrophage `annotation`;
  stage 01 renames this field to `sub_clusters`.
- Both objects must expose compatible `patient_id`, `condition`, `cell_type`, and `data` fields.

## Upstream assumptions to audit

- Stage 01 assigns combined matrix column names from only the fibroblast `var_names`; confirm the
  macrophage genes are identical and identically ordered before combining.
- Stages 03–04 use hard-coded numeric group indices and assume ten Fib states followed by three
  macrophage states. Resolve indices from names in adapted analyses.
- Stage 04 recognizes macrophage labels `Phagocytic`, `Pro-inflammatory`, and
  `Antigen-presenting`, and writes the lineage label as `Macropages`.
