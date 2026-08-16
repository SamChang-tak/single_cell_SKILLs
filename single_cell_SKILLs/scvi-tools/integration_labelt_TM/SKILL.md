---
name: integration-labelt-tm
description: Integrate labeled reference and unlabeled query scRNA-seq datasets and transfer cell-type labels with scVI/scANVI, following the official scvi-tools 1.3.3 Tabula Muris tutorial. Use when Codex needs to harmonize technologies such as Smart-seq2 and 10x, preserve raw or gene-length-adjusted counts, select batch-aware HVGs, train de novo SCVI/SCANVI models, predict query annotations, visualize latent spaces, or evaluate label transfer with held-out truth and confusion matrices.
---

# Tabula Muris integration and label transfer

Perform de novo reference/query integration and semi-supervised label transfer while keeping
preprocessing, technology effects, label resolution, and evaluation explicit.

## Route the task

- Read [references/tutorial-workflow.md](references/tutorial-workflow.md) for the exact versioned
  tutorial workflow, datasets, parameters, and keys.
- Read [references/label-transfer-guardrails.md](references/label-transfer-guardrails.md) before
  preparing counts, masking labels, or interpreting predictions.
- Use `scripts/run_label_transfer.py` for two preprocessed H5AD inputs. Inspect `--help` and run
  `--dry-run` first.
- Pin `scvi-tools==1.3.3` for API reproducibility, but note that the rendered tutorial reports its
  last execution with scvi-tools 1.2.1.

## Execute

1. Preserve source objects. Validate unique cells/genes, nonnegative integer-like count matrices,
   reference labels, query identity, and compatible species/gene identifiers.
2. Apply technology-specific count correction before concatenation when scientifically required.
   The tutorial gene-length-corrects Smart-seq2 reads, rescales by median gene length, and rounds.
3. Intersect genes explicitly, record losses by dataset, concatenate cells, and store model input
   counts in `layers['counts']` before any log normalization.
4. Select HVGs with `flavor='seurat_v3'`, `layer='counts'`, and the technology/batch key. Preserve
   the full expression space separately when downstream interpretation requires it.
5. Train scVI on all cells; save its model, registration, latent representation, history, seed,
   versions, and parameters.
6. Create a new label column: retain trusted reference labels and set every query cell to the
   declared unlabeled category. Confirm the category is not a real reference label.
7. Initialize scANVI from the scVI model trained on the exact same AnnData, train, predict labels,
   and save the model, `X_scANVI`, predictions, probabilities when available, and output AnnData.
8. Evaluate query predictions separately from reference training cells using appropriate truth,
   per-class metrics, confusion matrices, coverage, uncertainty, and biological marker review.

## Report

- Report reference/query cells, genes before/after intersection and HVG selection, batches,
  donors, label counts, count correction, versions, seeds, epochs, losses, and outputs.
- Distinguish de novo joint integration from online reference mapping; this tutorial demonstrates
  the former.
- Flag absent query cell types, label-resolution mismatches, rare classes, class imbalance,
  technology–biology confounding, low-confidence predictions, and unknown/novel populations.
- Never claim successful transfer from UMAP appearance alone.

## Provenance

Based on the official scvi-tools 1.3.3 tutorial, “Integration and label transfer with Tabula Muris”:
https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/tabula_muris.html
