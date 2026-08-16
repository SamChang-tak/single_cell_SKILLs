---
name: sscskin-endo-6
description: Reproduce, adapt, audit, and interpret the complete systemic-sclerosis skin endothelial workflow from Spatial_SSc_project/06_endo. Use for endothelial extraction and Scanpy reclustering, LEC/VEC subclustering, endothelial functional marker analysis, consolidation into LEC/CapEC/VenEC/ArtEC/EndoMT states, patient-level composition comparisons, or marker export with Scanpy and Seurat.
---

# SSc skin endothelial workflow

Analyze endothelial heterogeneity while separating exploratory clusters, manually assigned cell
states, marker evidence, and patient-level disease associations.

## Route the task

- Read [references/workflow-map.md](references/workflow-map.md) to choose the correct stage and
  checkpoint.
- Read [references/method-guardrails.md](references/method-guardrails.md) before interpreting
  endothelial identities, EndoMT, markers, or composition differences.
- Run `python scripts/audit_upstream.py assets/upstream/06_endo` to inventory notebook kernels and
  code cells without executing them.
- Consult all six verbatim notebooks under `assets/upstream/06_endo/` for exact paths,
  parameters, manual mappings, and plots.

## Execute safely

1. Preserve the integrated AnnData source and parameterize all `../../results/...` paths into a
   run-specific configuration.
2. Require unique cell IDs and explicit patient, condition, dataset, batch, and original cell-type
   metadata; report exclusions and cells per patient.
3. Record input layer, normalization, integration representation, neighbor/UMAP settings, Leiden
   resolution and seed, marker method, filtering, and every manual label mapping.
4. Save separate checkpoints for exploratory Endo1–Endo8 clusters and consolidated biological
   labels; never overwrite their provenance.
5. Validate outputs for dimensions, expected metadata, finite values, unique identifiers,
   patient coverage, and nonblank plots.

## Interpret and report

- Treat cluster labels as hypotheses supported by marker combinations, not single-marker proof.
- Distinguish lymphatic, capillary, venous, arterial, and EndoMT-like programs; test for doublets,
  ambient RNA, stress, and pericyte/fibroblast contamination.
- Treat EndoMT labels as transcriptional resemblance unless orthogonal evidence supports a
  transition or lineage claim.
- Use patients—not cells—as replicates for condition comparisons and report patient-level
  proportions, effect sizes, uncertainty, and dataset confounding.
- Report cluster stability, marker consistency across patients, multiple-testing correction,
  missing genes, sensitivity analyses, and all deviations from upstream code.

## Provenance

All files from `lzj1769/Spatial_SSc_project/06_endo` are bundled verbatim, with the repository
MIT license retained at `assets/upstream/LICENSE`.
