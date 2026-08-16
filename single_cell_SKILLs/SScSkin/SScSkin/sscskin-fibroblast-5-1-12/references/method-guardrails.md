# Method guardrails

- Stage 01 assumes `X_pca_harmony`; verify its covariates and test clustering stability. The
  upstream removal of Leiden cluster 9 and reassignment of cluster 11 require biological/QC
  justification and sensitivity analysis.
- Cell-level Wilcoxon significance does not provide patient-level replication. Confirm
  condition claims with patient-aware models or pseudobulk analysis.
- Slingshot/tradeSeq results depend on embedding, clusters, roots, lineages, knots, and coverage.
  Cross-sectional trajectories do not establish direction, time, or ancestry.
- GO results require an explicit tested universe, identifier mapping, database version,
  selection threshold, and multiple-testing correction.
- AddModuleScore depends on assay, binning, control genes, coverage, and object composition.
  Report mapped genes and avoid interpreting absolute scores as cell quantities.
- Compute spatial correlations within sections and summarize across patients. Account for
  autocorrelation, abundance, library size, multiple testing, and zero-variance features.
- MISTy importance is conditional predictive association. Record views, spatial radius,
  features, model, seed, performance, and per-patient consistency; do not claim causality.
- NABA ECM scores are gene-set summaries. Decoupler ULM scores depend on target-network
  provenance, weights, coverage, organism mapping, input layer, and minimum target count.
- Report sample/patient counts, software versions, missing inputs, identifier alignment,
  uncertainty, effect sizes, sensitivity checks, and all deviations from upstream notebooks.
