# Validation and interpretation guardrails

- Use raw count data. Verify nonnegative integer-like values and do not substitute `.X` unless its
  count provenance is established.
- Choose `batch_key` from the unwanted technical variation to model. If disease, donor, tissue,
  protocol, and study are nested or confounded, no integration method can identify their effects
  without assumptions; show the contingency tables and limit claims.
- Preserve rare and batch-specific cell types. Apparent mixing can be overcorrection when genuine
  biology is batch-specific.
- scANVI uses supplied labels as supervision. Audit inconsistent nomenclature, coarse/fine label
  mixtures, low-confidence annotations, and circular evaluation against training labels.
- Do not judge integration from UMAP alone. UMAP is stochastic and distorts global geometry.
  Evaluate latent neighbors, batch/label coverage, patient mixing, marker preservation, and scIB
  metrics.
- Benchmark against unintegrated PCA and, when important, alternative seeds/hyperparameters.
  Report the biology–batch tradeoff rather than a single aggregate score.
- Metric failures are informative. Small or single-batch labels may be skipped; KBET or other
  metrics may be unavailable. List exclusions and reasons.
- Use held-out donors or studies when assessing generalization or label transfer. Cell-level
  random splits overstate performance when donor effects remain.
- Save model weights and AnnData registration metadata together. Loading a model against changed
  genes, cells, category mappings, or layers can invalidate results.
- Inspect training loss, convergence, latent variance, library-size behavior, and sensitivity to
  latent dimension, layers, likelihood, epochs, seed, and highly variable gene selection.
