---
name: sscskin-fibroblast-5-13-21
description: Reproduce, adapt, audit, and interpret stages 13 through 21 of the systemic-sclerosis skin fibroblast workflow from Spatial_SSc_project/05_fibroblast. Use for checking fibroblast markers, visualizing MISTy associations, state–cell-type colocalization, testing markers across spatial niches, Fib1-versus-Fib2 differential expression and GO enrichment, deriving Seurat markers, spatial UCell scoring, or correlating UCell states with cell proportions.
---

# SSc skin fibroblast stages 13–21

Run the late fibroblast workflow as explicit marker, niche, differential-expression, and
spatial-scoring branches.

## Route the task

- Read [references/stage-map.md](references/stage-map.md) to select the relevant notebook and
  upstream checkpoint.
- Read [references/method-guardrails.md](references/method-guardrails.md) before interpreting
  markers, MISTy, niches, enrichment, module scores, or spatial correlations.
- Run `python scripts/audit_upstream.py assets/upstream/05_fibroblast_13_21` to inventory the
  bundled notebooks without executing them.
- Consult the nine verbatim notebooks under `assets/upstream/05_fibroblast_13_21/` for exact
  code, paths, parameters, and plots.

## Execute safely

1. Preserve all inputs and parameterize each `../../results/...` path into a run-specific
   configuration before execution.
2. Validate identifier-based joins among AnnData, Seurat, spatial niches, cell2location
   proportions, state scores, and MISTy tables; report missing and reordered identifiers.
3. Retain sample, patient, condition, and batch metadata, and use patients rather than cells or
   spots as biological replicates.
4. Record assay/layer, normalization, filters, gene IDs, thresholds, correction method,
   package/resource versions, seeds, and signature coverage.
5. Save machine-readable outputs and verify dimensions, unique identifiers, finite values,
   expected columns, and nonblank figures.

## Interpret and report

- Treat marker results as dataset- and contrast-dependent; confirm Fib1/Fib2 findings across
  patients and with patient-aware models.
- Treat MISTy importance and spatial correlations as conditional associations, not causal
  communication or regulation.
- Treat niche differences as conditional on the upstream niche model and verify that niches are
  represented across patients rather than driven by one section.
- Treat UCell scores as relative rank-based signature summaries, not cell abundance or absolute
  activity.
- Report effect sizes, uncertainty, multiple-testing correction, signature coverage, spatial
  autocorrelation, replicate consistency, confounding, missing inputs, and deviations.

## Provenance

The bundled notebooks are the complete numbered 13–21 subset from
`lzj1769/Spatial_SSc_project/05_fibroblast`; the repository MIT license is retained at
`assets/upstream/LICENSE`.
