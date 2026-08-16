# Method guardrails

- Marker and differential-expression results depend on assay/layer, filtering, contrast,
  detection threshold, and patient composition. Cell-level Wilcoxon tests do not supply
  patient-level replication; confirm biological contrasts with patient-aware or pseudobulk
  models.
- Define Fib1/Fib2 labels before testing, report cells and patients per group, inspect patient
  balance, and distinguish cluster markers from condition effects.
- For GO analysis, report the tested universe, identifier mapping, database version, ontology,
  gene-selection rule, direction of change, and multiple-testing correction.
- MISTy importance is predictive association conditional on views, radius, features, and model.
  Compare only harmonized runs and require replicate-level consistency.
- Niche-marker results inherit uncertainty from the upstream niche definition. Report niche
  membership by patient/condition and avoid spot-level pseudoreplication.
- UCell is rank-based but remains sensitive to marker specificity, signature size and coverage,
  ties, and data quality. Report mapped/total genes and do not interpret scores as abundance.
- Compute spatial correlations per section; address spatial autocorrelation, abundance,
  library-size effects, zero variance, and multiple testing. Summarize across patients.
- Report sample/patient counts, identifier alignment, versions, seeds, missing inputs, effect
  sizes, uncertainty, sensitivity checks, and all deviations from upstream notebooks.
