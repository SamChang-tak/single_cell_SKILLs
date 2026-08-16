# Stage map

| Stage | Kernel | Purpose | Main input | Main output |
|---|---|---|---|---|
| 01 | Python | Subset macrophages, UMAP/Leiden reclustering, Wilcoxon markers | integrated scRNA H5AD | `macrophages.h5ad`, proportions, marker tables |
| 02 | Python | Annotate macrophage states and rerun markers | stage-01 H5AD | `adata.h5ad`, `cell_proportion.csv`, `marker_genes.csv` |
| 03 | R | Perform GO analysis by macrophage cluster | stage-01 top markers | `go_analysis.csv` and plots |
| 04 | R | Score annotated macrophage states spatially with AddModuleScore | spatial RDS, cell2location CSV, stage-02 markers | per-sample RDS and score CSV |
| 05 | Python | Model macrophage–fibroblast spatial associations with LIANA MISTy | filtered spatial H5AD, macrophage/fibroblast scores | H5MU and interaction CSV |
| 06 | R | Aggregate and visualize MISTy outputs | stage-05 CSVs | comparison plots |
| 07 | R | Visualize macrophage-state scores per section | stage-04 RDS | spatial feature plots |
| 08 | R | Summarize macrophage-state scores across samples | stage-04 RDS | aggregate plots |
| 09 | R | Visualize macrophage cell2location estimates | spatial stage-10 RDS | spatial feature plots |
| 10 | R | Visualize macrophage marker expression spatially | spatial RDS | marker plots |
| 11 | Python | Infer TF activities with decoupler ULM | stage-02 H5AD and regulator network | `act.h5ad` and plots |
| 12 | R | Score stage-02 signatures spatially with UCell | spatial RDS, cell2location CSV, stage-02 markers | per-sample RDS and score CSV |
| 13 | R | Derive positive macrophage-state markers with Seurat | stage-02 H5AD | `markers.csv` |
| 14 | R | Score stage-13 signatures spatially with UCell | spatial RDS, cell2location CSV, stage-13 markers | per-sample RDS and score CSV |
| 15 | R | Correlate macrophage and fibroblast UCell states | macrophage stage-14 and fibroblast stage-20 RDS | per-sample and aggregate correlations |

## Dependency branches

- Single-cell backbone: 01 → 02 → 03, 04, 11, 12, or 13.
- AddModuleScore spatial branch: 02 → 04 → 05/07/08; stage 05 also needs fibroblast stage-06
  scores and filtered spatial AnnData.
- UCell branch A: 02 → 12.
- UCell branch B: 02 → 13 → 14 → 15; stage 15 also needs fibroblast stage-20 outputs.
- Stages 04, 09, 10, 12, and 14 depend on products from `04_spatial`.
