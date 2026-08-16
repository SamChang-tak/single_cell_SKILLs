# Workflow map

| Stage | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| 01 | Read annotated macrophage and consolidated endothelial H5AD objects, combine counts/metadata, normalize Seurat object, run PCA/UMAP | macrophage stage-02 H5AD; endothelial stage-04 H5AD | integrated `obj.rds` |
| 02 | Group by `sub_clusters`, load human CellChatDB, identify overexpressed features/interactions, compute triMean probabilities, filter groups below 10 cells, aggregate pathways | stage-01 Seurat RDS | CellChat `obj.rds` |
| 03 | Visualize interaction counts/weights, per-sender networks, and global ligand–receptor bubbles | stage-02 CellChat RDS | network and bubble plots |

## Input contract

- Macrophage input: `07_marchphages/02_annotate/adata.h5ad`, raw counts in `layers['counts']`,
  subtype in `annotation` renamed to `sub_clusters`.
- Endothelial input: `06_endo/04_merge_clusters/endothelial.h5ad`, raw counts in
  `layers['counts']`, subtype in `sub_clusters_v2` renamed to `sub_clusters`.
- Both inputs must contain compatible `patient_id`, `condition`, `cell_type`, and `data` fields.

## Fragile upstream assumptions

- Stage 01 row-binds the two count matrices and assigns combined columns from macrophage
  `var_names`. Confirm the endothelial gene set and order are identical before combining.
- Stage 03 creates sender-specific plots by matrix row position. Adapt this to subtype names so
  label ordering changes cannot silently alter results.
- The upstream model pools cells across patients and conditions. Use it descriptively only;
  condition claims require replicate-specific models and harmonized network comparison.
