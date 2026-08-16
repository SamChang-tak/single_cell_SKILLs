---
name: scrna-batch-sysvi
description: Integrate scRNA-seq datasets with substantial cross-system effects using scvi-tools SysVI. Use for cross-species, cross-technology, cell-versus-nucleus, primary-tissue-versus-organoid, or similarly difficult integration; preparing shared normalized log-transformed features; defining system and weaker covariates; tuning VampPrior and latent cycle-consistency; exporting latent embeddings; evaluating batch mixing, biological conservation, overcorrection, convergence, and seed sensitivity.
---

# scRNA-seq integration with sysVI

Use SysVI for integration where the dominant batch effect represents genuinely different systems. Treat integration as a balance between system mixing and preservation of biological structure.

## Route the task

- Read `references/tutorial-workflow.md` before preparing data or fitting SysVI.
- Read `references/evaluation-guardrails.md` before tuning correction strength or interpreting an integrated embedding.
- Use `scripts/run_sysvi.py` for a reproducible baseline on a preprocessed H5AD.

## Workflow

1. Audit the experimental design for confounding among system, condition, donor, tissue, cell type, and technology.
2. Harmonize feature identities across systems. Start from shared genes, select HVGs within each system using within-system batches, and use an intersection or otherwise defensible shared set.
3. Normalize every cell to a fixed library size and log-transform. Fit SysVI on these continuous features; do not pass raw counts merely because count models are common elsewhere in scvi-tools.
4. Use `batch_key` for the strongest system effect. Supply weaker categorical and continuous effects as covariates.
5. Combine multiple strong system dimensions when appropriate, such as species × technology, and use the same grouping during preprocessing.
6. Consider embedded categorical covariates when category cardinality makes one-hot encoding expensive.
7. Train with VampPrior and cycle consistency, inspect reconstruction, KL, and cycle losses, and repeat multiple seeds.
8. Export the latent representation, construct neighbors and UMAP, and inspect system mixing and cell-type conservation jointly.
9. Evaluate within each cell type and system. Compare against unintegrated data and at least one appropriate baseline.
10. Save the model, embedding, losses, preprocessing record, parameters, seed, and software versions.

## Tuning rules

- Start near the tutorial's cycle-consistency weight of 5; the documented usual range is 2–10, not a guarantee.
- Increase cycle weight only when justified by insufficient system correction; excessive correction can erase biology.
- Decrease cycle or KL weight when preservation is poor, then reevaluate both objectives.
- Do not choose a seed or hyperparameter solely for the prettiest UMAP. Use prespecified quantitative and qualitative criteria.
- Confirm loss stabilization; more cells may require fewer epochs, but cell count alone does not prove convergence.

## Reproducible baseline

```bash
python scripts/run_sysvi.py \
  --input preprocessed_shared_hvgs.h5ad \
  --output-dir sysvi_output \
  --system-key system \
  --categorical-covariate batch \
  --cycle-weight 5 \
  --max-epochs 200
```

The runner assumes `adata.X` already contains fixed-depth normalized, log-transformed shared features. Use repeated `--categorical-covariate` and `--continuous-covariate` options as needed. Run `--dry-run` to inspect configuration without importing analysis dependencies.

## Deliverables

Return the integrated coordinates, embedding H5AD with UMAP, model, training histories, run configuration, system-by-covariate tables, and an evaluation that reports both batch correction and biological conservation.
