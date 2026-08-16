---
name: sscskin-spatial-4-20-27
description: Reproduce, adapt, audit, and interpret stages 20 through 27 of the systemic-sclerosis skin spatial transcriptomics workflow from Spatial_SSc_project/04_spatial. Use when Codex needs to rank composition-niche markers, infer or visualize pathway and regulator activities with decoupler, integrate expression-based molecular niches with Harmony and Leiden, run or compare spatial CellChat networks, summarize ligand–receptor pathways, or inspect per-sample spatial pathway and transcription-factor activity outputs.
---

# SSc skin spatial stages 20–27

Analyze niche molecular programs and spatial communication while preserving the distinction
between measured expression, inferred activity, spatial association, and causal biology.

## Route the task

- Read [references/stage-map.md](references/stage-map.md) to select the correct stage and
  required upstream checkpoint.
- Read [references/method-guardrails.md](references/method-guardrails.md) before interpreting
  decoupler, molecular-niche, or CellChat results.
- Use `scripts/summarize_activity_by_group.py` to validate and summarize an activity AnnData
  matrix by niche, condition, or another observation group.
- Inspect `--help` and run `--dry-run` before using the script.
- Consult the eight verbatim notebooks under `assets/upstream/04_spatial_20_27/` for exact
  upstream code, resources, parameters, and visualization choices.

## Execute safely

1. Require stage-15 niche labels or filtered per-sample spatial AnnData inputs, with unique
   spot IDs and explicit `sample`, `patient_id`, and `condition` metadata.
2. Preserve all inputs and save new checkpoints for expression, activity, niche, and
   communication results.
3. For decoupler, record network/resource name, version, organism, source/target columns,
   weights, scoring method, minimum targets, gene-ID mapping, transform, and output layer.
4. Validate inferred activities against contributing target genes and spatial localization.
5. For molecular niches, record feature selection, scaling, PCA, Harmony variables, neighbor
   graph, Leiden resolution, seed, and cluster stability. The upstream baseline corrects by
   `patient_id` and uses resolution 0.3.
6. Align molecular-niche labels by expression/marker profiles before comparing reruns.
7. For CellChat, validate Seurat Spatial assay, hard group labels, tissue coordinates,
   scale-factor JSON, physical distance conversion, database subset, trimming, minimum group
   size, and output network slots.
8. Run CellChat per sample; never create communication edges across tissue sections.
9. Compare networks only after harmonizing cell-type labels, database version, filtering,
   scale factors, and normalization.
10. Use patients as biological replicates for healthy-versus-SSc claims.

## Validate and interpret

- Treat decoupler scores as model-inferred regulator/pathway activities, not direct pathway
  flux or protein activity measurements.
- Treat CellChat probabilities as database- and model-conditional communication hypotheses,
  not observed ligand–receptor binding or causal signaling.
- Treat molecular niches as unsupervised model states conditional on preprocessing and
  clustering choices.
- Report missing gene-set targets, low-coverage activities, rare groups, disconnected graphs,
  sample-specific niches, and unstable communication pathways.
- Avoid spot-level pseudoreplication; show patient-level summaries and effect sizes.
- Flag condition–patient–batch confounding and cases where Harmony may remove biology.
- Return checkpoints, activity matrices, network objects, summary tables, nonblank plots,
  versions, parameters, database/resource provenance, and warnings.

## Provenance

The bundled notebooks come from `lzj1769/Spatial_SSc_project/04_spatial`, stages 20–27,
and retain the repository MIT license at `assets/upstream/LICENSE`.
