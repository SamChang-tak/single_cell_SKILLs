---
name: sscskin-cell-2022
description: Reproduce, adapt, audit, and interpret the GSE195452 systemic-sclerosis (SSc) skin single-cell RNA-seq workflow from Spatial_SSc_project/01_cell_2022 using Scanpy and AnnData. Use when Codex needs to load the per-sample count tables, attach SRA and patient metadata, select SSc/control skin CD45+/CD90+ cells, apply the published QC and normalization choices, collapse detailed annotations into broad cell types, identify marker genes, create UMAP or marker plots, or review the original cell_2022 notebooks and patient metadata.
---

# SSc skin cell_2022

Reproduce the upstream `01_cell_2022` workflow while preserving raw counts,
recording deviations, and treating the published notebook parameters as study-specific
choices rather than universal defaults.

## Route the task

- Read [references/workflow.md](references/workflow.md) before implementing, running,
  or reviewing the full workflow.
- Read [references/annotation-map.md](references/annotation-map.md) when collapsing
  detailed annotations into broad cell types or interpreting marker results.
- Use `scripts/run_ssc_skin_cell_2022.py` for a clean end-to-end implementation.
  Inspect `--help` first and run `--dry-run` to validate paths and parameters.
- Consult the verbatim upstream files in `assets/upstream/01_cell_2022/` when exact
  notebook provenance, plotting choices, or exploratory outputs matter.
- Use `assets/upstream/01_cell_2022/patient.csv` only as the upstream cell-metadata
  snapshot; it contains no expression matrix.

## Execute safely

1. Require the annotation table, SRA RunInfo table, per-sample RNA count directory,
   patient workbook, and a new output directory.
2. Preserve the source files and write staged `.h5ad` checkpoints.
3. Confirm cell barcodes and sample names join uniquely before concatenating samples.
4. Report dimensions and relevant `.obs`, `.var`, `.layers`, `.obsm`, `.obsp`, and
   `.uns` keys after each stage.
5. Filter to conditions `SSC` and `Control`, tissues `Skin`/`skin`, selection markers
   `CD45+`/`CD90+`, and remove annotations `UN`, `GBP1`, and `NRXN1` only when
   reproducing the published cohort.
6. Apply the published QC thresholds (`min_genes=300`, `min_counts=500`,
   `n_genes_by_counts < 5000`, `pct_counts_mt < 30`) and report retained cells by
   patient, condition, and selection marker.
7. Save raw counts in `layers['counts']` before total-count normalization and `log1p`.
8. Recompute HVGs, PCA, neighbors, and UMAP; set a random seed and record versions.
9. Collapse annotations with the exact mapping reference, then inspect unmatched labels
   rather than silently accepting partial replacements.
10. Rank broad-cell-type markers with Wilcoxon and save the complete table, not only a
    visually selected subset.

## Validate and interpret

- Check patient and condition mixing on the embedding; flag patient-dominated islands.
- Use patients, not cells, as biological replicates for condition comparisons.
- Do not interpret broad annotation replacement as de novo cell-type discovery.
- Flag the original notebook path inconsistency (`results/cell_2022` versus
  `results/01_cell_2022`) and use one explicit output directory.
- Treat the original mitochondrial cutoff of 30% as permissive; show its distribution
  and justify any stricter sensitivity analysis.
- Record missing metadata joins, duplicated barcodes, sparse groups, and annotations not
  covered by the mapping.
- Return the checkpoint paths, plots, marker table, cell-retention summary, package
  versions, parameters, and all warnings.

## Provenance

The bundled upstream files come from
`lzj1769/Spatial_SSc_project/01_cell_2022` and retain the repository MIT license at
`assets/upstream/LICENSE`.
