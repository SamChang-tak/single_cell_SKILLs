---
name: sscskin-fibroblast-5-01-12
description: Reproduce, adapt, audit, and interpret stages 01 through 12 of the systemic-sclerosis skin fibroblast workflow from Spatial_SSc_project/05_fibroblast. Use for fibroblast reclustering and marker discovery, spatial marker visualization, Slingshot/tradeSeq trajectories, GO enrichment, spatial fibroblast-state scoring and colocalization, LIANA MISTy modeling, NABA ECM scoring, or decoupler transcription-factor activity.
---

# SSc skin fibroblast stages 01–12

Run the early fibroblast workflow as explicit single-cell, trajectory, and spatial-analysis
branches. The upstream repository has no stage 09, so this skill contains 11 notebooks.

## Route the task

- Read [references/stage-map.md](references/stage-map.md) to select a stage and its inputs.
- Read [references/method-guardrails.md](references/method-guardrails.md) before interpreting
  trajectories, scores, colocalization, MISTy, ECM, or TF activity.
- Run `python scripts/audit_upstream.py assets/upstream/05_fibroblast_01_12` to inspect kernels
  and code-cell counts without executing notebooks.
- Consult the verbatim notebooks under `assets/upstream/05_fibroblast_01_12/` for exact code,
  parameters, paths, and figures.

## Execute safely

1. Preserve source objects and parameterize every `../../results/...` path into a run-specific
   configuration before execution.
2. Retain unique cell/spot identifiers and explicit sample, patient, condition, and batch
   metadata across AnnData, Seurat, cell2location, and score-table joins.
3. Record package/resource versions, seeds, assay/layer, normalization, filtering, gene IDs,
   gene-set coverage, and every manual cluster-label decision.
4. Save and validate a checkpoint after each stage; test dimensions, identifiers, finite values,
   expected columns, and nonblank plots.
5. Analyze spatial sections separately and use patients—not cells or spots—as biological
   replicates for condition comparisons.

## Interpret and report

- Treat marker lists as descriptive evidence and cross-sectional trajectories as model-dependent
  hypotheses, not demonstrated ancestry or time.
- Treat module scores as relative gene-set summaries; do not compare AddModuleScore values
  across independently processed objects without calibration.
- Treat correlations and MISTy importance as spatial associations, not causal communication.
- Treat decoupler TF activity as network-based inference, distinct from TF expression or protein
  activation.
- Report sensitivity, gene-set/network coverage, multiple-testing correction, replicate-level
  consistency, confounding, missing inputs, and deviations from upstream code.

## Provenance

The bundled notebooks are the complete numbered 01–12 subset from
`lzj1769/Spatial_SSc_project/05_fibroblast`; the repository MIT license is retained at
`assets/upstream/LICENSE`.
