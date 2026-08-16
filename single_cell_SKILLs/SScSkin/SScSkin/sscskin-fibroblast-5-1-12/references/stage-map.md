# Stage map

| Stage | Kernel | Purpose | Main input | Main output |
|---|---|---|---|---|
| 01 | Python | Subset fibroblasts; neighbors, UMAP, Leiden, labels, Wilcoxon markers | integrated scRNA H5AD | `fibroblast.h5ad`, proportions, markers |
| 02 | Python | Expand stored rank-genes results | stage-01 H5AD | marker table |
| 03 | R | Plot fibroblast markers spatially | spatial Seurat RDS | spatial feature plots |
| 04 | R | Infer lineages and expression trends | stage-01 H5AD | trajectory RDS and plots |
| 05 | R | Test trajectory genes and GO terms | stage-04 trajectory RDS | GO table and plots |
| 06 | R | Score curated fibroblast states per section | spatial RDS, cell2location CSV, `all_genes.csv` | score CSV and RDS |
| 07 | R | Correlate fibroblast states with cell proportions | stage-06 RDS | per-sample and aggregate correlations |
| 08 | R | Correlate fibroblast and macrophage states | stage-06 plus macrophage RDS | correlation tables |
| 09 | — | No notebook exists upstream | — | — |
| 10 | Python | Fit intra/juxtaview spatial MISTy models | spatial H5AD, cell2location, stage-06 scores | H5MU and interaction CSV |
| 11 | R | Score NABA ECM gene sets | stage-01 H5AD, `NABAgsets.xls` | ECM score plots |
| 12 | Python | Infer TF activity with decoupler ULM | stage-01 H5AD and regulator network | `act.h5ad` and plots |

## Dependency branches

- Single-cell backbone: 01 → 02, 04, 11, or 12.
- Trajectory branch: 01 → 04 → 05.
- Spatial-state branch: 02 → 06 → 07, 08, or 10.
- Stage 08 additionally requires outputs from the separate macrophage workflow.
- Stages 06 and 11 refer to unnumbered source files `all_genes.csv` and `NABAgsets.xls` in the
  upstream `05_fibroblast` directory. They are outside the requested numbered 01–12 subset;
  obtain them from the full skill or original repository and record their checksums.
