---
name: ld-vae
description: Fit and interpret linearly decoded variational autoencoders for scRNA-seq with scvi-tools LinearSCVI (LDVAE). Use when Codex needs interpretable latent factors linked directly to gene loadings, count-aware probabilistic factor analysis, LinearSCVI training and convergence checks, latent-space visualization or clustering, top positive and negative loading genes, factor interpretation, or comparison with PCA and standard nonlinear scVI.
---

# Linear decoded VAE

Use `scvi.model.LinearSCVI` to learn count-aware latent factors whose linear decoder provides gene weights for interpretation. Keep the raw counts intact and distinguish exploratory factor interpretation from confirmatory biological evidence.

## Route the task

- Read `references/tutorial-workflow.md` before implementing or reproducing the official workflow.
- Read `references/interpretation-guardrails.md` before naming factors, comparing runs, or making biological claims.
- Use `scripts/run_ldvae.py` for a reproducible H5AD-to-results baseline.

## Workflow

1. Inspect AnnData, identify the raw-count matrix, and verify nonnegative finite values and unique gene/cell identifiers.
2. Preserve raw counts in a layer. Use normalized/log-transformed values only for separate visualization tasks.
3. Select informative genes from raw counts, normally with `seurat_v3`; record the exact gene set. Do not select genes using the condition being tested in a way that leaks outcomes.
4. Register counts with `LinearSCVI.setup_anndata`, including a batch key or other justified covariates when needed.
5. Fit `LinearSCVI` with a documented latent dimension, seed, epochs, learning rate, and accelerator environment.
6. Inspect training and validation ELBO curves. Revisit optimization if they diverge, oscillate severely, or have not stabilized.
7. Export `get_latent_representation()` for cells and `get_loadings()` for genes.
8. Examine both positive and negative high-magnitude loading genes for every factor. Use coherent gene programs, known covariates, and cell-level factor distributions to interpret them.
9. Build neighbors, UMAP, or clusters from the latent representation only when they answer the analysis question.
10. Save the model, aligned AnnData, loadings, latent coordinates, training history, parameters, and software versions.

## Required interpretation rules

- Do not rank factors by index. Unlike PCA, `Z_0` does not necessarily explain more variation than `Z_1`.
- Treat factor sign and ordering as non-identifiable: signs can flip and factors can permute between fits.
- Interpret a loading as a model coefficient, not a differential-expression effect size or causal relationship.
- Check whether factors track batch, donor, library size, mitochondrial fraction, cell cycle, or other technical variables.
- Test robustness across seeds, gene selections, latent dimensions, and held-out biological replicates before emphasizing a factor.
- Prefer a standard nonlinear SCVI model when predictive flexibility or integration quality matters more than decoder interpretability.

## Reproducible baseline

```bash
python scripts/run_ldvae.py \
  --input data.h5ad \
  --output-dir ldvae_output \
  --count-layer counts \
  --n-latent 10 \
  --max-epochs 250 \
  --learning-rate 0.005
```

Add `--batch-key` only when the AnnData column and scientific design justify batch adjustment. Use `--select-hvg 1000` to perform tutorial-style HVG selection; omit it when the input was already deliberately subset. Run `--dry-run` to inspect configuration without importing analysis packages.

## Deliverables

Return latent coordinates, gene loadings, top loading genes in both directions, convergence history, annotated H5AD, saved model, run configuration, and an interpretation that explicitly addresses technical covariates and robustness.
