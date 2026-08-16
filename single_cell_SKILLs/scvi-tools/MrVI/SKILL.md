---
name: mrvi
description: Analyze comparable multi-sample scRNA-seq datasets with scvi-tools MRVI (Multi-resolution Variational Inference). Use when Codex needs sample-invariant and sample-specific latent representations, local sample-distance matrices, donor-level heterogeneity, sample-covariate-linked differential expression, local differential abundance, multi-donor disease comparisons, nuisance-batch correction, convergence checks, or biologically careful interpretation of cell-state-resolved effects across many samples.
---

# MrVI multi-sample analysis

Use MrVI when cells come from many comparable biological samples from the same tissue or cell line. Preserve the distinction between the target sample identity and nuisance batch variables.

## Route the task

- Read `references/tutorial-workflow.md` before fitting or reproducing MrVI.
- Read `references/inference-guardrails.md` before local distance, differential-expression, or differential-abundance interpretation.
- Use `scripts/run_mrvi.py` for a reproducible baseline that fits a model and exports the `u` representation.

## Workflow

1. Audit sample counts, cells per sample, conditions, nuisance batches, repeated measures, and missing sample metadata.
2. Confirm samples are biologically comparable. Do not combine unrelated tissues or cell lines and expect single-cell-resolution sample comparisons to remain valid.
3. Use raw, nonnegative counts and select HVGs using a method appropriate for multi-sample data. Preserve the complete preprocessing record.
4. Set `sample_key` to biological sample/donor identity. Set `batch_key` only to a nuisance variable that should be corrected.
5. Fit `MRVI`, inspect validation ELBO, and test stability across seeds and reasonable feature/model choices.
6. Extract `u` for broad cell states invariant to sample and nuisance covariates. Remember that the model's `z` representation augments `u` with sample-specific effects while correcting nuisance effects.
7. Compute local sample distances with `keep_cell=False` and a biologically meaningful `groupby` when full cell-specific output would be unnecessarily large.
8. Relate sample distances to sample-level metadata without treating cells as independent sample replicates.
9. For covariate-linked DE or DA, verify sample-level design, contrast/reference coding, minimum group support, and confounding before modeling.
10. Save model, annotated AnnData, latent coordinates, histories, local-distance objects, design tables, parameters, and versions.

## Required interpretation rules

- Treat the biological sample as the replication unit for sample-level covariates.
- Do not infer a condition effect when condition is perfectly confounded with site, batch, donor, or another nuisance variable.
- Distinguish `u`-space cell states, sample-specific effects, local sample distances, DE effect sizes, and DA log-probability contrasts; they answer different questions.
- Interpret proximity or clustering of samples as model-based similarity, not causality or proof of equivalence.
- Validate findings across samples and relevant cell populations; global results can hide sparse or donor-specific effects.
- Correct for multiple testing and report uncertainty for DE/DA rather than selecting only large visual effects.

## Reproducible baseline

```bash
python scripts/run_mrvi.py \
  --input counts.h5ad \
  --output-dir mrvi_output \
  --sample-key patient_id \
  --count-layer counts \
  --select-hvg 10000 \
  --max-epochs 400
```

Add `--batch-key Site` only when site is a genuine nuisance variable. Add `--distance-groupby cell_type` to export group-averaged local sample distances. Run `--dry-run` to inspect configuration without importing analysis dependencies.

## Deliverables

Return `u` coordinates, UMAP, fitted model, convergence histories, annotated H5AD, sample design table, optional grouped sample distances, and an interpretation that respects sample-level replication and contrast coding.
