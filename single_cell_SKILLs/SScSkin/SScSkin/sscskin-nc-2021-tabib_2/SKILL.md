---
name: sscskin-nc-2021-tabib
description: Reproduce, adapt, audit, and interpret the GSE138669 systemic-sclerosis skin single-cell RNA-seq workflow from Spatial_SSc_project/02_nc_2021_Tabib using Scanpy and AnnData. Use when Codex needs to load 10x HDF5 matrices, attach Tabib study conditions, perform the published skin-cell QC and preprocessing, cluster with Leiden, apply or review the published cluster-to-cell-type annotations, re-cluster the macrophage/cDC compartment, rank marker genes, generate UMAP or marker plots, or inspect the original nc_2021_Tabib notebooks.
---

# SSc skin nc_2021 Tabib

Reproduce the upstream GSE138669 workflow while preserving counts, recording package
versions and parameters, and separating observed marker evidence from manual cluster labels.

## Route the task

- Read [references/workflow.md](references/workflow.md) before implementing, running,
  or reviewing the workflow.
- Read [references/cluster-annotations.md](references/cluster-annotations.md) before
  transferring the published cluster labels or interpreting marker plots.
- Use `scripts/run_ssc_skin_nc_2021_tabib.py` for a deterministic end-to-end baseline.
  Inspect `--help` and run `--dry-run` before execution.
- Consult the six verbatim notebooks in `assets/upstream/02_nc_2021_Tabib/` when exact
  provenance, exploratory outputs, or marker panels matter.

## Execute safely

1. Require the GSE138669 SRA RunInfo table, directory of 10x `.h5` matrices, and a new
   output directory. Never modify the source files.
2. Validate a unique SRA metadata row for every filename-derived sample name.
3. Preserve raw counts in `layers['counts']` before normalization.
4. Reproduce the published QC thresholds only when requested: initial `min_genes=200`,
   then `min_genes=300`, `min_counts=500`, `n_genes_by_counts < 5000`, and
   `pct_counts_mt < 5`.
5. Normalize to 10,000 counts, apply `log1p`, select HVGs, then compute PCA, neighbors,
   UMAP, and Leiden with an explicit seed.
6. Inspect batch and condition mixing before annotation. The upstream “integration”
   notebook does not actually perform batch correction.
7. Apply the published cluster map only after marker profiles support equivalent cluster
   identities. Cluster numbers are not portable across package versions or reruns.
8. Re-cluster the published macrophage/cDC compartment at Leiden resolution 0.1 only
   when the parent cluster is marker-consistent with upstream cluster 8.
9. Rank markers with Wilcoxon and save complete statistics as well as compact top tables.
10. Write staged `.h5ad` checkpoints, plots, cell-retention summaries, parameters,
    versions, and warnings.

## Validate and interpret

- Report dimensions and relevant AnnData keys after loading, QC, clustering, and annotation.
- Confirm non-negative integer-like counts before normalization and finite embeddings afterward.
- Check cluster sizes and condition/sample composition; flag sparse or sample-specific clusters.
- Validate every manual label with positive and exclusion markers. Do not infer identity from
  cluster number alone.
- Use biological samples, not cells, as replicates for SSc-versus-control inference.
- Treat UMAP separation as visualization, not statistical evidence of batch correction.
- Flag the upstream output-path inconsistency (`results/nc_2021_Tabib` versus
  `results/02_nc_2021_Tabib`) and use one explicit output directory.
- Return checkpoint paths, plots, marker tables, retained-cell counts, package versions,
  parameters, and annotation/batch warnings.

## Provenance

The bundled files come from `lzj1769/Spatial_SSc_project/02_nc_2021_Tabib` and retain
the repository MIT license at `assets/upstream/LICENSE`.
