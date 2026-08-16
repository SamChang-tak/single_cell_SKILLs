# Method guardrails

- Verify the integrated representation and test Leiden stability across neighbor settings,
  resolution, seed, and patient holdouts. Require markers and patient distribution to support
  manual macrophage-state labels.
- Cell-level Wilcoxon tests do not establish patient-level replication. Confirm condition claims
  using patient-aware or pseudobulk models.
- For GO analysis, state the tested universe, organism database and version, identifier mapping,
  ontology, gene-selection rule, direction, and multiple-testing correction.
- AddModuleScore depends on assay, expression bins, controls, coverage, and dataset composition.
  UCell is rank-based but still sensitive to marker specificity, coverage, ties, and data quality.
  Do not compare their numerical scales directly.
- Report mapped/total genes for every signature and flag shared or nonspecific markers.
- Cell2location outputs are posterior abundance estimates, not observed cell counts. Confirm
  barcode alignment and record posterior statistic, reference labels, priors, and model version.
- MISTy importance is predictive association conditional on views, features, radius, and model;
  it does not demonstrate intercellular signaling. Require model performance and consistency
  across patients.
- Decoupler ULM activity depends on regulator-network provenance, weights, organism mapping,
  target coverage, minimum targets, and input layer; distinguish it from TF expression.
- Compute macrophage–fibroblast correlations within each section, address spatial
  autocorrelation and abundance, correct multiple tests, and summarize across patients.
- Report sample/patient counts, missing inputs, identifier alignment, versions, seeds, effect
  sizes, uncertainty, sensitivity checks, and all deviations from upstream code.
