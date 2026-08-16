# Official tutorial workflow

Source: scvi-tools 1.3.3, **Differential expression on C. elegans data**:
https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/scVI_DE_worm.html

The rendered tutorial reports that its last execution used scvi-tools 1.1.6. Inspect installed signatures and record versions when reproducing it.

## Canonical sequence

1. Start from nonnegative integer count data with cell labels and batch metadata.
2. Use `scvi.data.poisson_gene_selection` or another documented selection method. The tutorial retains selected genes, copies raw counts to a `counts` layer, and registers `batch_key` with `SCVI.setup_anndata`.
3. Fit `SCVI(..., gene_likelihood="nb")` with validation each epoch, early stopping, and an explicit patience. Plot train and validation ELBO after the initial transient.
4. Obtain `X_scVI`, calculate neighbors and UMAP, and inspect biological labels and batch structure.
5. Define two cell sets with boolean indices or query strings.
6. Run `model.differential_expression(idx1=..., idx2=..., mode="change")`.

## DE interpretation

- `lfc_mean > 0`: higher in population 1.
- `proba_de`: posterior probability that the absolute LFC exceeds `delta`.
- `is_de_fdr_0.05`: Bayesian FDR decision at the specified target.
- `weights="importance"` sharpens posterior scores for calibrated gene sets; pair it with `filter_outlier_cells=True`.
- `weights="uniform"` is useful for ranking and sensitivity comparisons.
- `delta` defines the biologically relevant LFC interval around zero; the tutorial example uses the API default of 0.25.
- A small pseudocount can stabilize comparisons dominated by rarely expressed genes.
- `batch_correction=True` is meaningful only when the compared cell sets arise from overlapping batches.

## One-versus-rest

Use `model.differential_expression(groupby="cell_type", mode="change")` for every group versus rest. For targeted reproducibility, prefer an explicit pairwise or one-versus-rest contrast with audited cell counts and direction.

## Important tutorial caveat

The tutorial calls `log10(proba_not_de)` a plotting score. This is not a frequentist p-value. Label transformed axes as posterior non-DE probability scores.
