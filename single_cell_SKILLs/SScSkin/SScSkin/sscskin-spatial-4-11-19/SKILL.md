---
name: sscskin-spatial-4-11-19
description: Reproduce, adapt, audit, and interpret stages 11 through 19 of the systemic-sclerosis skin spatial transcriptomics workflow from Spatial_SSc_project/04_spatial. Use when Codex needs to compute or visualize Squidpy neighborhood enrichment from cell2location outputs, quantify inferred-cell-type colocalization, integrate spot-level cell-abundance profiles, define and characterize composition-based spatial niches, compare niche prevalence between healthy and SSc patients, create niche heatmaps, or map niche labels back onto Visium tissue sections.
---

# SSc skin spatial stages 11–19

Analyze cell2location-derived spatial organization from per-section neighborhoods through
cross-sample niche definition and patient-level niche comparison.

## Route the task

- Read [references/stage-map.md](references/stage-map.md) to identify stage inputs and
  outputs.
- Read [references/niche-analysis.md](references/niche-analysis.md) before computing
  colocalization, clustering niches, or testing condition differences.
- Use `scripts/run_squidpy_neighborhood.py` for deterministic stage-11 neighborhood
  enrichment.
- Use `scripts/summarize_niche_prevalence.py` for patient-level niche proportions and
  qualified healthy-versus-SSc comparisons.
- Inspect `--help` and run `--dry-run` before executing either script.
- Consult the nine verbatim notebooks under `assets/upstream/04_spatial_11_19/` for exact
  upstream code and plotting choices.

## Execute safely

1. Require stage-09 cell2location H5AD outputs with spatial coordinates, posterior abundance
   estimates, hard predicted labels when used, and explicit sample/patient/condition metadata.
2. Preserve all per-sample inputs and create new outputs for every stage.
3. Build spatial graphs independently within each tissue section. Report disconnected
   components, edge counts, rare groups, and permutation settings.
4. Save neighborhood adjacency counts and enrichment z-scores, not only heatmaps.
5. Define colocalization explicitly: abundance correlation, spatial adjacency, or another
   statistic. Do not mix these quantities under one label.
6. When concatenating cell-abundance profiles, align cell-type columns by name, fill only
   structurally absent values with justified zeros, and retain original spot and sample IDs.
7. Record feature scaling, neighbor count, representation, Leiden resolution, random seed,
   and cluster sizes for niche definition. The upstream baseline uses 30 neighbors and
   Leiden resolution 0.5.
8. Characterize each niche by normalized cell-type abundance and show patient/condition
   composition before assigning biological names.
9. Calculate niche prevalence within each patient before condition comparison.
10. Join niche labels back to Seurat spatial objects by barcode and sample, never row order.

## Validate and interpret

- Treat cell2location values as posterior abundance estimates, not observed cell counts.
- Recognize compositional closure: increased abundance/proportion of one type can induce
  negative associations with others.
- Distinguish correlation of abundances from physical cell–cell adjacency.
- Check niche stability across seeds, neighbor counts, resolutions, abundance summaries,
  and sample leave-one-out analyses.
- Align niche labels by composition before comparing reruns; numeric Leiden labels are not
  inherently equivalent.
- Use patients, not spots, as replicates. Report effect sizes and patient-level points even
  when formal testing is underpowered.
- Flag sample-specific niches, sparse conditions, missing cell types, tissue edge effects,
  disconnected graphs, and condition–batch confounding.
- Return checkpoints, matrices, patient-level tables, plots, parameters, versions, and all
  warnings.

## Provenance

The bundled notebooks come from `lzj1769/Spatial_SSc_project/04_spatial`, stages 11–19,
and retain the repository MIT license at `assets/upstream/LICENSE`.
