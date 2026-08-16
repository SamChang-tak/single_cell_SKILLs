# Label-transfer guardrails

- Use the tutorial’s gene-length adjustment only for non-UMI full-length read counts when its
  assumptions and gene-length resource match the organism and annotation release.
- Never feed log-normalized or scaled values into the scVI count layer. Check nonnegativity,
  integer-likeness after justified correction, sparsity, library sizes, and outliers.
- Harmonize gene IDs explicitly and report genes lost from each dataset. Resolve duplicated
  symbols and versioned Ensembl IDs before concatenation.
- Mask query labels without chained assignment. Verify every query cell is unlabeled and every
  retained reference label is trusted. Keep evaluation truth in a separate, inaccessible column
  during training to prevent leakage.
- Reference labels constrain scANVI. Harmonize ontology and resolution; coarse query truth may
  legitimately disagree with finer transferred labels.
- A query cell type absent from the reference cannot be correctly named by closed-set prediction.
  Use uncertainty, neighborhood novelty, reconstruction diagnostics, and marker review to flag
  out-of-reference cells rather than forcing labels.
- Split and evaluate by donor/study, not random cells, when estimating generalization. Report
  macro/per-class precision, recall, F1, confusion, support, rejection/unknown rate, and calibration.
- Inspect class imbalance and `n_samples_per_label`; rare reference types can be unstable or
  oversampled. Run seed and parameter sensitivity analyses.
- Do not infer integration quality from UMAP alone. Check technology mixing within cell types,
  biological separation, marker preservation, latent neighbors, and replicate-aware metrics.
- This is de novo joint training. For frozen-reference or sequential mapping, use the appropriate
  scArches reference-mapping workflow and document query adaptation.
