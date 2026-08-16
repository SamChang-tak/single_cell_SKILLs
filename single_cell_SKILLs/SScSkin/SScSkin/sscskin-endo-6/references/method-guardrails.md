# Method guardrails

- Confirm the integrated representation used for neighbors and whether correction covariates
  could remove condition biology. Test stability across neighbor counts, Leiden resolutions,
  seeds, dataset holdouts, and patient holdouts.
- Document the stage-01 patient-selection list and assess whether exclusions are related to QC,
  dataset, or condition. Do not generalize beyond included patients.
- Validate state labels using coherent panels: lymphatic identity, capillary/venous/arterial
  specialization, endothelial core genes, mesenchymal programs, and exclusion markers.
- EndoMT-like expression can reflect activation, dissociation stress, ambient collagen, doublets,
  or contamination. Avoid claiming lineage transition from cross-sectional scRNA-seq alone.
- Reconcile the stage-02 LEC/VEC variable names with the different stage-04 mapping before
  reusing subcluster results.
- Cell-level Wilcoxon tests can overstate evidence when cells from one patient dominate. Confirm
  condition-specific markers with patient-aware or pseudobulk models and report effect sizes.
- Compare proportions at the patient level with uncertainty; test dataset/batch confounding and
  avoid treating cells as independent replicates.
- Scanpy and Seurat marker outputs may differ because of assay, layer, filtering, transform,
  prevalence threshold, and test implementation. Record these settings before comparison.
- Report patient/cell counts, cluster stability, marker coverage, multiple-testing correction,
  versions, seeds, missing inputs, sensitivity checks, and deviations from upstream notebooks.
