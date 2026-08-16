---
name: integration-lung
description: Integrate lung single-cell RNA-seq datasets with scvi-tools using scVI and label-informed scANVI, following the official scvi-tools 1.3.3 atlas harmonization tutorial. Use when Codex needs to register raw-count AnnData batches, train SCVI/SCANVI models, create latent embeddings and UMAPs, benchmark biological conservation against batch correction with scib-metrics, save reproducible integration outputs, or diagnose overcorrection, confounding, and weak batch mixing.
---

# scvi-tools lung integration

Integrate scRNA-seq batches with explicit data contracts and evaluate whether technical mixing is
achieved without erasing biological structure.

## Route the task

- Read [references/tutorial-workflow.md](references/tutorial-workflow.md) for the versioned
  tutorial steps, parameters, keys, and example dataset.
- Read [references/validation-guardrails.md](references/validation-guardrails.md) before choosing
  covariates or interpreting latent spaces and metrics.
- Use `scripts/run_lung_integration.py` for a reproducible scVI/scANVI run. Inspect `--help` and
  start with `--dry-run`.
- Pin `scvi-tools==1.3.3` when exact API reproducibility is required; record the actual installed
  version because the rendered tutorial reports its last execution with 1.3.2.

## Execute

1. Load AnnData and validate unique cells/genes, nonnegative integer-like raw counts, categorical
   batch labels, and optional cell-type labels.
2. Preserve the source object. Store raw counts in a dedicated layer and never pass scaled or
   log-normalized expression as the scVI count layer.
3. Select the biological sampling unit for `batch_key`. Do not automatically use donor, disease,
   protocol, or study if it is confounded with biology that must be retained.
4. Register data with `SCVI.setup_anndata(adata, layer=count_layer, batch_key=batch_key)`.
5. Train scVI, save `X_scVI`, the model, training history, parameters, versions, and seed.
6. Use scANVI only when labels are sufficiently reliable. Initialize it from an scVI model trained
   on the exact same AnnData, provide `labels_key` and a non-colliding unlabeled category, train,
   and save `X_scANVI` plus its model.
7. Build neighbors and UMAP separately for each latent representation; do not overwrite graphs
   when comparisons are needed.
8. Evaluate PCA, scVI, and scANVI using replicate-aware visual diagnostics and scIB metrics.

## Report

- Report dimensions, genes, cells, batches, donors, labels, count layer, category sizes, missing
  values, software/hardware, seeds, model parameters, epochs, losses, and output paths.
- Compare biological conservation and batch correction rather than optimizing only one side.
- Show batch and cell-type composition, per-label/per-batch coverage, embeddings, quantitative
  metrics, rare groups, warnings, and sensitivity analyses.
- Flag nested/confounded batches, labels unique to one batch, weak label quality, failed training,
  posterior collapse, unstable neighbors, and metric exclusions.

## Provenance

This skill follows the official scvi-tools 1.3.3 tutorial, “Atlas-level integration of lung data”:
https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/harmonization.html
