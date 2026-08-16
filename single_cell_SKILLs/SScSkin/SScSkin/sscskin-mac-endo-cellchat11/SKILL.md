---
name: sscskin-mac-endo-cellchat11
description: Reproduce, adapt, audit, and interpret the complete systemic-sclerosis skin macrophage–endothelial CellChat workflow from Spatial_SSc_project/11_mac_endo. Use for combining annotated macrophage and endothelial AnnData counts in Seurat, constructing a human CellChat network, or visualizing macrophage–endothelial interaction counts, strengths, senders, receivers, and ligand–receptor bubbles.
---

# SSc macrophage–endothelial CellChat workflow

Analyze inferred communication between macrophage and endothelial states while distinguishing
expression, database compatibility, model probability, and biological mechanism.

## Route the task

- Read [references/workflow-map.md](references/workflow-map.md) for stages, input contracts, and
  upstream dependencies.
- Read [references/cellchat-guardrails.md](references/cellchat-guardrails.md) before running or
  interpreting CellChat.
- Run `python scripts/audit_upstream.py assets/upstream/11_mac_endo` to inventory the bundled
  notebooks without executing them.
- Consult the three verbatim R notebooks under `assets/upstream/11_mac_endo/` for exact code,
  parameters, paths, and plots.

## Execute safely

1. Preserve macrophage and endothelial inputs; parameterize relative paths and write a new run.
2. Prove identical gene identifiers and ordering, unique cell IDs, compatible raw-count layers,
   and matching matrix orientation before combining AnnData matrices.
3. Harmonize subtype labels while retaining lineage, patient, condition, dataset, and original
   annotations; detect name collisions across lineages.
4. Record CellChat and database versions, database subset, assay, normalization, overexpression
   settings, probability estimator, population-size option, minimum cells, and seed.
5. Use subtype names rather than numeric positions for sender/receiver selection and plotting.
6. Save and validate the integrated Seurat object, CellChat object, interaction tables, plots,
   parameters, group sizes, database provenance, and warnings.

## Interpret and report

- Treat CellChat edges as database- and model-conditional hypotheses, not observed binding,
  causal direction, or signaling flux.
- Distinguish edge counts, communication weights/probabilities, pathway aggregates, and gene
  expression; do not compare them as equivalent quantities.
- Run inference per patient or tissue section for biological comparisons. Do not use pooled cells
  as independent replicates.
- Verify ligand and all receptor/cofactor subunits in the corresponding source and target groups.
- Report rare groups removed by `min.cells = 10`, missing genes, group sizes, database coverage,
  multiple testing, replicate consistency, sensitivity, and condition–dataset confounding.

## Provenance

All files from `lzj1769/Spatial_SSc_project/11_mac_endo` are bundled verbatim, with the repository
MIT license at `assets/upstream/LICENSE`. Macrophage upstream paths retain the repository spelling
`07_marchphages`.
