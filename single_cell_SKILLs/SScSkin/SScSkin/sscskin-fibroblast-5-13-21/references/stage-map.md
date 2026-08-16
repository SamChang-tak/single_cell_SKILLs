# Stage map

| Stage | Kernel | Purpose | Main input | Main output |
|---|---|---|---|---|
| 13 | Python | Inspect candidate fibroblast markers | stage-01 `fibroblast.h5ad` | diagnostic plots |
| 14 | R | Aggregate and visualize MISTy interactions | stage-10 MISTy CSVs | comparison plots |
| 15 | R | Correlate fibroblast-state scores with cell proportions | stage-06 spatial RDS | per-sample correlation CSVs |
| 16 | R | Compare fibroblast marker expression across spatial niches | stage-06 RDS, spatial `niches.csv` | niche summaries and plots |
| 17 | Python | Test Fib1 versus Fib2 with Scanpy Wilcoxon | stage-01 H5AD | `fib1_vs_fib2.csv` |
| 18 | R | Filter Fib1/Fib2 DE genes and perform GO analysis | stage-17 CSV | DEG and GO CSVs, plots |
| 19 | R | Derive positive fibroblast-state markers with Seurat | stage-01 H5AD | `markers.csv` |
| 20 | R | Score stage-19 signatures spatially with UCell | spatial RDS, cell2location CSV, stage-19 markers | per-sample RDS and score CSV |
| 21 | R | Correlate UCell state scores with cell proportions | stage-20 RDS | per-sample/all correlations and mean matrix |

## Dependency branches

- Marker diagnostics: stage 01 → 13.
- MISTy visualization: stages 06 and 10 → 14.
- Curated-state spatial analyses: stage 06 → 15 or 16; stage 16 also needs the stage-15 spatial
  niche assignments from the separate `04_spatial` workflow.
- Fib1/Fib2 comparison: stage 01 → 17 → 18.
- Data-derived UCell branch: stage 01 → 19 → 20 → 21; stages 20–21 also require spatial RDS
  objects and cell2location outputs.

Most required checkpoints are produced by stages 01–12 or by `04_spatial`; they are not bundled
in this numbered notebook subset.
