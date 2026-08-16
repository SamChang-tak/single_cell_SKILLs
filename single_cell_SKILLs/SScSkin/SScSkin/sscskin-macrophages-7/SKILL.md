---
name: sscskin-macrophages-7
description: Reproduce, adapt, audit, and interpret the complete systemic-sclerosis skin macrophage workflow from Spatial_SSc_project/07_marchphages. Use for macrophage reclustering and annotation, marker and GO analysis, spatial macrophage-state scoring with Seurat or UCell, macrophage–fibroblast MISTy modeling and colocalization, spatial visualization, cell2location review, or decoupler transcription-factor activity.
---

# SSc skin macrophage workflow

Analyze macrophage states across single-cell and spatial data while separating measured
expression, signature scores, deconvolution estimates, predictive associations, and biological
mechanism.

## Route the task

- Read [references/stage-map.md](references/stage-map.md) to choose a notebook branch and locate
  its upstream checkpoints.
- Read [references/method-guardrails.md](references/method-guardrails.md) before interpreting
  markers, enrichment, scores, MISTy, TF activity, or macrophage–fibroblast colocalization.
- Run `python scripts/audit_upstream.py assets/upstream/07_marchphages` to inventory notebook
  kernels and code cells without execution.
- Consult all 15 verbatim notebooks under `assets/upstream/07_marchphages/` for exact paths,
  parameters, code, and plots.

## Execute safely

1. Preserve source objects and parameterize all `../../results/...` paths into a run-specific
   configuration.
2. Maintain unique cell/spot identifiers plus sample, patient, condition, and batch metadata
   across AnnData, Seurat, cell2location, score tables, and MISTy objects.
3. Record assay/layer, normalization, integration representation, filters, cluster-label rules,
   signature/network provenance, versions, seeds, and sample inclusion.
4. Validate each join by identifiers and report dropped, duplicated, or reordered rows.
5. Save stage checkpoints and verify dimensions, expected columns, finite values, unique IDs,
   nonempty model outputs, and nonblank figures.
6. Fit spatial models within tissue sections and use patients—not cells or spots—as biological
   replicates for condition comparisons.

## Interpret and report

- Treat cluster markers as descriptive and verify macrophage-state labels across patients,
  parameter settings, and external marker evidence.
- Treat module/UCell scores as relative signature summaries, not direct cell abundance or
  absolute pathway activity.
- Treat cell2location values as model estimates and MISTy importance or correlations as
  conditional spatial associations, not causal communication.
- Treat decoupler TF activity as network-based inference distinct from TF expression or protein
  activation.
- Report effect sizes, uncertainty, multiple-testing correction, gene coverage, spatial
  autocorrelation, replicate consistency, confounding, missing inputs, and deviations.

## Provenance

All files from upstream `lzj1769/Spatial_SSc_project/07_marchphages` are bundled verbatim. The
repository uses the misspelling `marchphages`; retain it in source paths. The MIT license is at
`assets/upstream/LICENSE`.
