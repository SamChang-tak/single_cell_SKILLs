---
name: sscskin-fib-macrophages-9
description: Reproduce, adapt, audit, and interpret the complete systemic-sclerosis skin fibroblast–macrophage workflow from Spatial_SSc_project/09_fib_macrophages. Use for combining fibroblast and macrophage AnnData counts in Seurat, constructing a CellChat network, visualizing pathway and sender–receiver interactions, or summarizing ligand–receptor counts between fibroblast and macrophage subtypes.
---

# SSc skin fibroblast–macrophage workflow

Analyze inferred communication between fibroblast and macrophage states while keeping expression,
database-supported ligand–receptor compatibility, model probability, and biological mechanism
distinct.

## Route the task

- Read [references/workflow-map.md](references/workflow-map.md) for inputs, outputs, and stage
  dependencies.
- Read [references/cellchat-guardrails.md](references/cellchat-guardrails.md) before running or
  interpreting CellChat.
- Run `python scripts/audit_upstream.py assets/upstream/09_fib_macrophages` to inventory all
  bundled source files without executing notebooks.
- Consult the four verbatim notebooks and `all_interaction.csv` under
  `assets/upstream/09_fib_macrophages/` for exact upstream code and parameters.

## Execute safely

1. Preserve the fibroblast and macrophage input objects; parameterize all relative paths and
   write to a new run directory.
2. Confirm both AnnData objects use identical gene identifiers, unique cell IDs, raw count
   layers, compatible count orientation, and harmonized metadata before row-binding counts.
3. Preserve `patient_id`, `condition`, original lineage, and subtype. Validate subtype names and
   ordering rather than relying on numeric CellChat indices.
4. Record CellChat and database versions, database subset, expression assay, normalization,
   overexpression settings, probability estimator, population-size option, minimum cells, and
   random seed.
5. Save the integrated Seurat object, CellChat object, full interaction table, filtered tables,
   plots, parameters, versions, and validation logs.

## Interpret and report

- Treat CellChat edges as database- and model-conditional communication hypotheses, not observed
  binding, direction of causality, or signaling flux.
- Report expression support for ligand and receptor subunits and distinguish interaction counts
  from interaction strengths or probabilities.
- Avoid pooled-cell condition claims from the upstream combined object. Run models per patient or
  tissue section, harmonize labels/settings, then summarize replicate-level effects.
- Report group sizes, rare groups removed by `min.cells = 10`, missing genes, database coverage,
  uncertain complexes, multiple comparisons, and sensitivity to model settings.
- Flag the upstream `Macropages` spelling in derived source/target labels and normalize it only
  with an explicit mapping.

## Provenance

All files from `lzj1769/Spatial_SSc_project/09_fib_macrophages` are bundled verbatim, with the
repository MIT license retained at `assets/upstream/LICENSE`.
