# Spatial interpretation and statistics

## Neighborhood enrichment

Squidpy neighborhood enrichment compares observed cluster adjacencies with label
permutations on a fixed spatial graph. Report graph construction, permutations, cluster
sizes, z-scores, and sample-level results. Analyze sections independently and combine
patient-level summaries afterward.

## Colocalization

Correlations among inferred cell-type abundances can reflect true co-occurrence, shared spot
mixtures, compositional constraints, reference cross-talk, or tissue architecture. Avoid
calling correlation direct cellular interaction.

## Niches

Cell-composition niches and molecular-expression niches answer different questions. Record
the feature matrix, scaling, batch correction, neighbor graph, resolution, and label-matching
procedure. Compare niche prevalence using patient-level proportions and show all patients.

## Functional activity

decoupler scores are model-derived gene-set or regulator activities. Report the resource,
version, scoring method, minimum targets, filtering, and direction. Validate key results with
underlying genes and spatial localization.

## CellChat

Ligand–receptor probability is conditional on expression, database, group labels, spatial
distance model, and scale factors. Validate coordinates and physical units, avoid comparing
raw network weights across samples without normalization, and summarize at the patient level.

## Condition comparisons

Use healthy and SSc patients as replicates. Do not use individual spots as independent
replicates. When the sample set is small or unbalanced, report effect sizes, patient-level
distributions, and uncertainty rather than relying on nominal spot-level p-values.
