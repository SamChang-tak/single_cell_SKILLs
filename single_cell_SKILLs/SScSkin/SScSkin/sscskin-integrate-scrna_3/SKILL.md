---
name: sscskin-integrate-scrna3
description: Integrate, batch-correct, audit, visualize, and export the Gur 2022 GSE195452 and Tabib 2021 GSE138669 systemic-sclerosis skin single-cell RNA-seq datasets following Spatial_SSc_project/03_integrate_scrna. Use when Codex needs to harmonize metadata and cell types across the two AnnData objects, concatenate shared genes and raw-count layers, normalize and embed the joint atlas, run Harmony by study, validate study/condition mixing, rank broad-cell-type markers, visualize canonical markers or NAMPT/NNMT, perform carefully qualified fibroblast comparisons, or convert the integrated AnnData object to Seurat.
---

# SSc skin integrate_scrna3

Integrate the Gur 2022 and Tabib 2021 skin datasets while preserving raw counts and
distinguishing batch alignment from biological evidence.

## Route the task

- Read [references/workflow.md](references/workflow.md) before running or reviewing the
  integration.
- Read [references/markers-and-interpretation.md](references/markers-and-interpretation.md)
  for feature plots, NAMPT/NNMT analyses, and interpretation constraints.
- Use `scripts/run_integrate_ssc_skin.py` for a reproducible Python implementation.
  Inspect `--help` and run `--dry-run` first.
- Use `scripts/anndata_to_seurat.R` only after validating `layers['counts']`, `X_pca`,
  `X_pca_harmony`, and `X_umap` in the integrated checkpoint.
- Consult all seven verbatim notebooks and the upstream marker-heatmap PDF under
  `assets/upstream/03_integrate_scrna/` when exact provenance or exploratory output matters.

## Execute safely

1. Require annotated Gur 2022 and Tabib 2021 `.h5ad` inputs and a new output directory.
2. Preserve both inputs; inspect dimensions, gene identifiers, count layers, observation
   columns, embeddings, and cell-type labels before concatenation.
3. Standardize metadata to `patient_id`, `condition`, `cell_type`, and `data`; map controls
   to `Healthy` and systemic sclerosis to `SSc` with exact value mapping.
4. Report common genes, source-exclusive genes, shared cell types, and source-exclusive cell
   types. Use an inner gene join to reproduce the upstream workflow.
5. Restore `X` from `layers['counts']`; validate non-negative integer-like counts; normalize
   to 10,000, apply `log1p`, and select HVGs with the published mean/dispersion thresholds.
6. Compute PCA, neighbors, and a pre-Harmony UMAP with a fixed seed.
7. Run Harmony on PCA using `data` (study) as the upstream batch variable, store
   `X_pca_harmony`, rebuild neighbors, and compute the final UMAP.
8. Compare pre/post study mixing within each cell type and inspect whether condition is
   confounded with study or patient.
9. Rank complete broad-cell-type markers with Wilcoxon and save canonical marker,
   NAMPT/NNMT, and condition-stratified visualizations when genes exist.
10. Save checkpoints, tables, plots, versions, parameters, and warnings.

## Validate and interpret

- Confirm all transition from counts to normalized expression is explicit and only once.
- Confirm the Harmony matrix has one row per cell, finite values, and the expected components.
- Never claim that UMAP mixing proves successful biological integration.
- Harmony correction by study can remove true biology when study, protocol, condition, or
  patient composition are confounded. Show uncorrected and corrected embeddings.
- Do not transfer a cell type absent from one source without marker review.
- Treat cell-level SSc-versus-Healthy Wilcoxon tests as exploratory pseudoreplication. Use
  patient-aware pseudobulk or mixed models for biological inference.
- Preserve original cell identifiers during Seurat conversion unless renaming is explicitly
  requested, and save an ID mapping if names change.
- Return shared-feature/type summaries, integration QC, checkpoints, plots, marker tables,
  software versions, and all warnings.

## Provenance

The bundled files come from `lzj1769/Spatial_SSc_project/03_integrate_scrna` and retain
the repository MIT license at `assets/upstream/LICENSE`.
