# Integration evaluation guardrails

## Design audit

- Tabulate system against donor, condition, tissue, cell type, chemistry, and sample. Do not claim correction can resolve perfect confounding.
- Require meaningful shared biology or anchors across systems. A model cannot validate correspondences absent from the design.
- Match gene identifiers and feature characteristics before modeling; cross-species mappings need documented orthology choices.

## Evaluate two competing objectives

Assess system removal with metrics such as batch ASW, graph connectivity, iLISI, or kBET where appropriate. Assess biological conservation with cell-type ASW, cLISI, label transfer, trajectory conservation, marker coherence, and donor-held-out checks. Use metrics suitable for the design rather than a single aggregate score.

Inspect metrics by cell type, system, and replicate. Global averages can hide rare-population loss, inappropriate merging, or a failed system.

## Diagnose failure

- Poor mixing: review feature harmonization, system definition, optimization, cycle weight, and whether populations truly overlap.
- Biological collapse: lower correction strength, inspect confounding, compare unintegrated structure, and verify markers within systems.
- Seed instability: run several prespecified seeds and report variability rather than selecting only a favorable visualization.
- Unexpected alignment: test whether it is supported by markers and independent biological knowledge; do not infer equivalence from proximity alone.

## Reporting

Report shared-feature construction, normalization target, transform, system and covariate definitions, category counts, weights, epochs, seeds, versions, losses, baselines, metrics, and system-stratified visual checks. Do not use an integrated embedding as direct evidence of differential abundance, differential expression, lineage, or causality.
