---
name: sscskin-spatial
description: Reproduce, adapt, audit, and interpret the systemic-sclerosis skin spatial transcriptomics workflow from Spatial_SSc_project/04_spatial. Use when Codex needs to load or QC 10x Visium samples, visualize spatial markers, map single cells with CellTrek, deconvolve spots with cell2location, cluster spatial domains with BANKSY, compute Squidpy neighborhood enrichment or cell colocalization, derive cell-composition or molecular niches, compare niches between healthy and SSc samples, score pathways or transcription factors with decoupler, run spatial CellChat, or trace data handoffs across the 28 upstream R and Python notebooks.
---

# SSc skin spatial

Run the spatial workflow as explicit, checkpointed modules. Preserve raw Visium data,
histology coordinates, sample identity, and deconvolution uncertainty throughout.

## Route the task

- Read [references/pipeline-map.md](references/pipeline-map.md) to select the correct
  notebook stage and its upstream dependencies.
- Read [references/data-contracts.md](references/data-contracts.md) before moving data
  between Seurat, AnnData, CellTrek, cell2location, BANKSY, Squidpy, or CellChat.
- Read [references/interpretation.md](references/interpretation.md) before comparing
  neighborhoods, niches, pathways, or communication networks.
- Use `scripts/run_squidpy_neighborhood.py` for reproducible spatial-neighbor and
  neighborhood-enrichment analysis. Inspect `--help` and run `--dry-run` first.
- Consult the 28 verbatim upstream notebooks in `assets/upstream/04_spatial/` for exact
  method parameters and the original mixed R/Python workflow.

## Execute safely

1. Inventory sample directories, filtered matrices, tissue images, scale factors, spatial
   coordinates, reference scRNA-seq checkpoints, package versions, and available compute.
2. Preserve source data and write a checkpoint after every module. Never overwrite a
   Seurat or AnnData input in place.
3. Validate coordinate orientation, units, image scale factors, spot/cell identifiers,
   gene identifiers, count integrity, and sample identity before modeling.
4. Process samples independently through loading/QC and deconvolution; combine samples only
   after adding explicit `patient_id`, `condition`, and `sample` fields.
5. Match each requested result to one primary mapping method: CellTrek for mapped single
   cells or cell2location for spot-level abundance. Do not conflate their outputs.
6. Record model versions, seeds, training epochs, priors, detection-alpha values, reference
   signatures, convergence diagnostics, and output keys for cell2location.
7. Construct spatial graphs within samples only. Never create edges across tissue sections.
8. Define niches from stated features—cell proportions or molecular expression—and record
   scaling, graph, batch correction, resolution, and cluster stability.
9. Use patients as replicates for SSc-versus-healthy comparisons. Retain spot-level analyses
   only as descriptive or exploratory.
10. Save tables, nonblank plots, checkpoints, parameters, versions, and warnings with an
    explicit output manifest.

## Validate each method

- Visium: confirm in-tissue spots, count distributions, mitochondrial fraction, image/spot
  registration, and minimum-feature filtering.
- CellTrek: inspect anchor transfer, prediction distances, mapped-cell density, collisions,
  and sensitivity to `dist_thresh`, `top_spot`, `spot_n`, and repulsion settings.
- cell2location: inspect reference signatures, training histories, posterior abundance
  quantiles, total abundance, reconstruction, and stability across seeds/prior choices.
- BANKSY: inspect coordinate scale, `k_geom`, lambda, resolution, and domain stability.
- Squidpy: validate graph connectivity and cluster counts; report enrichment z-scores and
  permutation settings, not only heatmaps.
- Niche analysis: check patient/batch composition and sensitivity to feature scaling,
  neighbor count, Harmony correction, and Leiden resolution.
- CellChat: verify spatial-distance units and scale factors, minimum group sizes, database
  subset, and per-sample network comparability.
- decoupler: report network source/version, statistic, minimum targets, and whether scores
  are inferred activities rather than measured pathway flux.

## Guardrails

- Do not interpret deconvolution abundance as direct cell counts without calibration.
- Do not treat spatial adjacency, correlation, or ligand–receptor scores as causal contact.
- Do not pool spots across patients as independent biological replicates.
- Do not compare raw niche numbers across reruns unless labels are aligned by composition.
- Flag missing samples, disconnected graphs, rare cell types, edge effects, low-resolution
  Visium mixtures, source/reference mismatch, and condition–batch confounding.

## Provenance

The bundled files come from `lzj1769/Spatial_SSc_project/04_spatial` and retain the
repository MIT license at `assets/upstream/LICENSE`.
