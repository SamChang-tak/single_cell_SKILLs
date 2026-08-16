---
name: cellassign
description: Annotate scRNA-seq cells probabilistically with scvi-tools CellAssign using a gene-by-cell-type marker matrix. Use for marker-based cell typing without labeled reference cells, preparing counts and size factors, fitting CellAssign, inspecting soft assignment probabilities and convergence, validating annotations across samples or batches, and diagnosing ambiguous, missing, or poorly specified cell types.
---

# CellAssign annotation

Use CellAssign as a probabilistic, marker-guided annotation model. Preserve uncertainty and treat the supplied marker categories as a closed candidate set, not as ground truth.

## Route the task

- Read `references/tutorial-workflow.md` before implementing or reviewing an analysis.
- Read `references/annotation-guardrails.md` when interpreting assignments, choosing markers, or diagnosing low-confidence results.
- Use `scripts/run_cellassign.py` for a reproducible baseline run from an H5AD file and marker CSV.

## Required inputs

Require:

1. An AnnData object containing nonnegative, unlogged counts in `X` or a named layer.
2. A CSV whose rows are gene identifiers, columns are candidate cell types, and values are binary marker indicators.
3. Matching gene identifiers between `adata.var_names` and the marker-matrix index.

Ask which matrix contains raw counts if this is unclear. Do not use scaled, centered, log-normalized, or batch-corrected expression as model input.

## Workflow

1. Inspect AnnData dimensions, layers, observation metadata, gene identifiers, sparsity, and count integrity.
2. Validate marker-matrix uniqueness, binary values, cell-type names, marker counts, and overlap with AnnData genes.
3. Compute each cell's library size from the full count matrix before restricting genes. Store `library_size / mean(library_size)` as the size factor; do not log it.
4. Align and order the marker matrix and AnnData to exactly the same marker genes. Report missing markers and reject any candidate cell type left without markers.
5. Register AnnData with `CellAssign.setup_anndata(..., size_factor_key=...)`, initialize `CellAssign`, and train it.
6. Inspect validation ELBO history for convergence or instability.
7. Export the full probability matrix. Derive the top assignment together with maximum probability, runner-up probability, margin, and entropy.
8. Plot assignments and uncertainty on the embedding. Compare assignments by donor, batch, condition, cluster, and any trusted labels.
9. Save aligned data, probabilities, diagnostics, software versions, parameters, and missing-marker reports.

## Interpretation rules

- Retain soft probabilities; an `idxmax` label alone hides ambiguity.
- Avoid inventing a universal confidence cutoff. Examine probability distributions, marker evidence, and the biological decision being made.
- Flag cells with low maximum probability, a small top-two margin, or high entropy for review.
- Treat unexpected assignments as hypotheses. Check marker specificity, ambient RNA, doublets, low library size, batch effects, and absent candidate types.
- Validate within each biological replicate. A result driven by one donor is not robust evidence of a general population.
- State that CellAssign can only allocate among provided types; it does not discover omitted or novel populations.

## Reproducible baseline

```bash
python scripts/run_cellassign.py \
  --input data.h5ad \
  --markers markers.csv \
  --output-dir cellassign_output \
  --count-layer counts
```

Omit `--count-layer` only when `adata.X` contains raw counts. Run with `--dry-run` to inspect the planned configuration without importing analysis dependencies.

## Deliverables

Provide:

- `cellassign_probabilities.csv`: soft probabilities for every cell and candidate type.
- `cellassign_assignments.csv`: top label and uncertainty summaries.
- `cellassign_annotated.h5ad`: marker-gene AnnData with assignments and probabilities.
- `missing_marker_genes.txt`: markers unavailable in the dataset.
- `run_config.json`: inputs, parameters, dimensions, overlap, and software versions.

Explain limitations and uncertainty alongside biological conclusions.
