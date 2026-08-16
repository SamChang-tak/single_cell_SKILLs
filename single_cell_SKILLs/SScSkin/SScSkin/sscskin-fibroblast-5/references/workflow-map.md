# Workflow map

The upstream folder contains 20 notebooks. Number 09 is absent in the source repository.
Notebook paths assume the repository's relative `results` layout.

| Stage | Kernel | Purpose | Principal inputs | Principal outputs |
|---|---|---|---|---|
| 01 | Python | Subset integrated fibroblasts, neighbors/UMAP, Leiden reclustering, annotation and Wilcoxon markers | `03_integrate_scrna/integrated.h5ad` | `fibroblast.h5ad`, proportions, top markers |
| 02 | Python | Expand stored rank-genes results into a marker table | stage-01 `fibroblast.h5ad` | `all_genes.csv` |
| 03 | R | Display fibroblast markers in a spatial Seurat object | stage-04 spatial RDS | spatial marker plots |
| 04 | R | Infer fibroblast lineages and smooth expression trends | stage-01 `fibroblast.h5ad` | `fibroblast.traj.rds`, trajectory plots |
| 05 | R | Test trajectory-associated genes and GO terms | stage-04 trajectory object | `GO_analysis.csv`, enrichment plots |
| 06 | R | Score curated fibroblast states in each spatial sample | spatial RDS, cell2location CSV, bundled `all_genes.csv` | per-sample RDS and score CSV |
| 07 | R | Correlate stage-06 fibroblast scores with cell proportions | stage-06 RDS files | per-sample/all correlations and mean matrix |
| 08 | R | Correlate fibroblast and macrophage state scores | fibroblast stage-06 and macrophage stage-04 RDS | per-sample and aggregate correlations |
| 10 | Python | Model intra- and juxtaview spatial associations with LIANA MISTy | filtered spatial H5AD, cell2location CSV, stage-06 scores | H5MU models and interaction CSVs |
| 11 | R | Score NABA extracellular-matrix gene sets in fibroblasts | stage-01 H5AD, bundled `NABAgsets.xls` | ECM score plots |
| 12 | Python | Infer transcription-factor activities with decoupler ULM | stage-01 H5AD, regulator network | `act.h5ad`, activity plots |
| 13 | Python | Inspect candidate fibroblast markers | stage-01 H5AD | diagnostic plots |
| 14 | R | Aggregate and visualize MISTy interaction results | stage-10 CSVs | comparison plots |
| 15 | R | Alternate fibroblast-state/cell-proportion colocalization | stage-06 RDS files | per-sample correlation CSVs |
| 16 | R | Test fibroblast marker expression by spatial niche | stage-06 RDS, stage-15 spatial `niches.csv` | niche marker summaries/plots |
| 17 | Python | Differential expression for Fib1 versus Fib2 | stage-01 H5AD | `fib1_vs_fib2.csv` |
| 18 | R | Filter Fib1/Fib2 DE genes and run GO analysis | stage-17 CSV | `fib1_fib2_deg.csv`, `GO_analysis.csv` |
| 19 | R | Find positive markers for fibroblast states in Seurat | stage-01 H5AD | `markers.csv` |
| 20 | R | Score stage-19 signatures spatially with UCell | spatial RDS, cell2location CSV, stage-19 markers | per-sample RDS and score CSV |
| 21 | R | Correlate UCell states with cell proportions | stage-20 RDS files | per-sample/all correlations and mean matrix |

## Branches and dependencies

- The single-cell backbone is 01 → 02, 04, 11, 12, 13, 17, or 19.
- The curated spatial-state branch is 02 → 06 → 07/08/10/15/16.
- The data-derived spatial-state branch is 19 → 20 → 21.
- The Fib1/Fib2 comparison is 17 → 18.
- Stage 08 depends on the separate macrophage workflow; stage 16 depends on spatial niche
  integration. These inputs are not bundled here.

## Bundled source tables

- `NABAgsets.xls`: NABA extracellular-matrix gene sets used by stage 11.
- `all_genes.csv`: source marker table used by the curated scoring branch.
- `fib1_fib2_deg.csv` and `top_200_marker_genes.csv`: source analysis tables retained exactly as
  upstream artifacts; verify their provenance and columns before treating them as fresh output.
