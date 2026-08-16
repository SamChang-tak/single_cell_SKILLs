---
name: sscskin-spatial-4-1-10
description: Reproduce, adapt, audit, and interpret stages 01 through 10 of the systemic-sclerosis skin spatial transcriptomics workflow from Spatial_SSc_project/04_spatial. Use when Codex needs to load or QC 10x Visium samples, visualize spatial markers, map integrated single cells with CellTrek, attach or run cell2location deconvolution, cluster spatial domains with BANKSY, compute or visualize initial Squidpy neighborhood enrichment, or manage the stage-01–10 Seurat/AnnData handoffs.
---

# SSc skin spatial stages 01–10

Run the early spatial workflow as explicit per-sample modules. This skill contains every
notebook numbered 01 through 10, including both notebooks numbered 07.

## Route the task

- Read [references/stage-map.md](references/stage-map.md) to select the notebook and its
  required inputs.
- Read [references/data-contracts.md](references/data-contracts.md) before transferring
  data among Seurat, AnnData, CellTrek, cell2location, BANKSY, and Squidpy.
- Use `scripts/run_squidpy_neighborhood.py` for deterministic stage-07 neighborhood
  enrichment. Inspect `--help` and run `--dry-run` first.
- Consult the 11 verbatim notebooks under `assets/upstream/04_spatial_01_10/` for exact
  parameters and original R/Python code.

## Execute safely

1. Inventory 10x Visium matrices, spatial positions, scale-factor JSON files, tissue images,
   sample metadata, and integrated scRNA-seq reference checkpoints.
2. Preserve all source data and save a new checkpoint after every stage.
3. Validate barcode alignment, count integrity, gene identifiers, coordinate orientation,
   coordinate units, and image/spot registration before mapping.
4. Run Visium loading and QC per tissue section; retain explicit `sample`, `patient_id`, and
   `condition` metadata.
5. Choose CellTrek for mapped reference cells or cell2location for posterior spot-level cell
   abundance. Do not describe one method's output as the other.
6. For CellTrek, record anchors and mapping settings, including `dist_thresh`, `top_spot`,
   `spot_n`, repulsion radius, and iterations.
7. For cell2location, record reference labels/counts, training history, posterior statistic,
   cells-per-location prior, detection alpha, seeds, and output keys.
8. For BANKSY, record coordinate scale, `k_geom`, lambda, component count, clustering method,
   resolution, and domain stability.
9. Construct Squidpy graphs within one section at a time; report graph components, cluster
   sizes, permutation count, adjacency counts, and enrichment z-scores.
10. Verify identifier-based joins when adding CSV deconvolution results to Seurat objects.

## Validate and interpret

- Confirm nonblank tissue images and spatial plots visually.
- Flag missing samples, low-feature spots, disconnected graphs, rare predicted cell types,
  edge effects, reference mismatch, and failed or unstable model convergence.
- Treat deconvolution abundances as model estimates, not direct cell counts.
- Treat neighborhood enrichment as adjacency relative to label permutations, not causal
  cell–cell interaction.
- Keep samples separate for graph inference and use patients as replicates for condition
  comparisons.
- Return checkpoints, plots, tables, parameters, package versions, sample coverage, and all
  warnings.

## Provenance

The bundled notebooks come from `lzj1769/Spatial_SSc_project/04_spatial`, stages 01–10,
and retain the repository MIT license at `assets/upstream/LICENSE`.
