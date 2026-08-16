---
name: sscskin-fibroblast-5
description: Reproduce, adapt, audit, and interpret the complete systemic-sclerosis skin fibroblast workflow from Spatial_SSc_project/05_fibroblast. Use for fibroblast reclustering and markers, Slingshot/tradeSeq trajectories and GO analysis, spatial fibroblast-state scoring, colocalization, MISTy, ECM and transcription-factor programs, molecular-niche markers, or Fib1-versus-Fib2 comparisons in Scanpy, Seurat, and spatial transcriptomics data.
---

# SSc skin fibroblast workflow

Analyze fibroblast states across single-cell and spatial data while keeping measured expression,
module scores, inferred trajectories, spatial associations, and mechanistic claims distinct.

## Route the task

- Read [references/workflow-map.md](references/workflow-map.md) to choose the relevant notebook
  branch, input checkpoint, and output.
- Read [references/method-guardrails.md](references/method-guardrails.md) before interpreting
  trajectories, state scores, correlations, MISTy, enrichment, or regulator activity.
- Run `python scripts/audit_upstream.py assets/upstream/05_fibroblast` to inventory the bundled
  notebooks, kernels, code cells, and local data files before adapting them.
- Consult the 20 verbatim notebooks and four source tables in
  `assets/upstream/05_fibroblast/` for exact code, parameters, paths, and figures.

## Execute safely

1. Work from a copy of each input and write new checkpoints under a run-specific output root.
2. Preserve unique cell or spot identifiers across AnnData, Seurat, cell2location, niche, and
   module-score tables; report dropped and reordered identifiers.
3. Recreate the upstream relative directory layout or parameterize every `../../results/...`
   path before execution. Do not silently substitute similarly named files.
4. Record package and resource versions, random seeds, gene identifiers, normalization layer,
   assay, covariates, filtering, and sample inclusion.
5. Treat patients or tissue sections—not cells or spots—as biological replicates for condition
   comparisons. Model repeated samples explicitly.
6. Validate every saved table for dimensions, unique identifiers, finite values, and expected
   columns; validate every figure as nonblank.
7. Keep exploratory upstream label edits and thresholds visible in the report, especially the
   removal of Leiden cluster 9 and reassignment of cluster 11 in notebook 01.

## Interpret and report

- Separate descriptive marker evidence from stable fibroblast-state annotation.
- Report trajectory sensitivity to root, clustering, embedding, filtering, and lineage count.
- Describe module scores as relative gene-set summaries and MISTy or TF outputs as
  model-dependent associations, not causal regulation.
- For spatial correlations, show per-sample estimates and patient-level summaries; do not pool
  spots across sections as independent observations.
- Report gene-set coverage, zero-variance features, multiple-testing method, low-abundance
  states, missing samples, and condition–patient–batch confounding.
- Return checkpoint paths, machine-readable tables, plots, parameters, versions, provenance,
  warnings, and a concise biological interpretation with limitations.

## Provenance

All files under the requested upstream `05_fibroblast` directory are bundled verbatim from
`lzj1769/Spatial_SSc_project`; the repository MIT license is retained at
`assets/upstream/LICENSE`.
