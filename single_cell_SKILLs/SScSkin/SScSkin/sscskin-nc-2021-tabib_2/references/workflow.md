# Upstream workflow and reproducibility notes

## Contents

- [Notebook sequence](#notebook-sequence)
- [Required inputs](#required-inputs)
- [Published preprocessing](#published-preprocessing)
- [Outputs](#outputs)
- [Known issues](#known-issues)

## Notebook sequence

1. `01_loading_data.ipynb`: read GSE138669 SRA metadata; load every 10x HDF5 file;
   filter cells with fewer than 200 genes; add sample and condition; concatenate; write
   `scrna.h5ad`.
2. `02_integration.ipynb`: apply cell and mitochondrial QC; preserve counts; normalize,
   log-transform, select HVGs, calculate PCA/neighbors/UMAP; write
   `scrna.integrated.h5ad`.
3. `03_clustering.ipynb`: run default Leiden; rank markers; inspect curated marker panels;
   map clusters to broad cell types; write `scrna.integrated.clustered.h5ad`.
4. `04_annotation.ipynb`: revise the manual cluster mapping and write
   `scrna.integrated.annotated.h5ad`.
5. `05_annotation.ipynb`: isolate parent cluster 8, re-cluster at resolution 0.1, label
   subcluster 0 as cDC and subcluster 1 as macrophages, merge annotations, and write the
   final `scrna.integrated.annotated.v2.h5ad`.
6. `04_viz.ipynb`: visualize broad marker panels and an upstream 30-cluster marker set.

## Required inputs

- `GSE138669/SraRunTable.txt`, with `Sample Name` and `CONDITION`.
- GSE138669 10x HDF5 count matrices. Sample names are inferred from the filename segment
  before the first underscore.

## Published preprocessing

- Initial load filter: at least 200 detected genes.
- Subsequent filters: at least 300 detected genes and 500 counts.
- Maximum detected genes: strictly below 5,000.
- Maximum mitochondrial percentage: strictly below 5%.
- Normalization: total counts to 10,000 followed by natural-log `log1p`.
- Default Scanpy HVG, ARPACK PCA, neighbor graph, UMAP, and Leiden settings.
- Macrophage/cDC subclustering: Leiden resolution 0.1.

## Outputs

The bundled runner writes raw, processed, clustered, and annotated checkpoints; QC and UMAP
plots; complete and top marker tables; cluster/sample/condition counts; cluster-label
validation summaries; and run metadata under one output directory.

## Known issues

- The “integration” notebook imports `scvi` but never uses it; the clustering notebook
  imports Harmony but never calls it. No batch correction is performed upstream.
- Output paths alternate between `results/nc_2021_Tabib` and
  `results/02_nc_2021_Tabib`.
- `03_clustering.ipynb` and `04_annotation.ipynb` disagree on some keratinocyte/secretory
  assignments; treat the later `04_annotation.ipynb` mapping as the final broad map before
  the macrophage/cDC refinement.
- Default Leiden results depend on Scanpy, igraph/leidenalg versions, graph construction,
  and random state. Published numeric cluster IDs cannot be assumed on a new run.
- Warnings are globally suppressed upstream and package versions/seeds are not recorded.
- Marker visualization includes labels for 30 clusters from an external/upstream clustering
  scheme, while the broad mapping uses Leiden clusters 0–27.
