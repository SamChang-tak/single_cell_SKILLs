# Colocalization and niche analysis

## Required metadata

Every spot must retain a unique identifier plus `sample`, `patient_id`, and `condition`.
Every abundance column must name one reference cell type and identify the cell2location
posterior statistic used. Spatial coordinates must remain in `obsm['spatial']`.

## Neighborhood enrichment

Compute spatial graphs per section. Squidpy neighborhood enrichment compares observed label
adjacencies with permutations on the fixed graph. Report cluster sizes, permutations,
adjacency counts, z-scores, and disconnected components. Do not merge sections into one graph.

## Colocalization

Abundance correlation within Visium spots can reflect biological co-occurrence, shared spot
mixtures, posterior cross-talk, tissue architecture, and compositional constraints. Report
the correlation method, transformation, missing-value handling, and sample-level matrices.
Aggregate at patient level before condition comparison.

## Composition-based niches

State whether features are posterior abundance, normalized proportions, or transformed
values. Record scaling and graph construction. Check cluster size, sample composition, and
stability before naming niches. Compare alternative neighbor counts and Leiden resolutions,
and perform sample leave-one-out sensitivity checks when feasible.

## Niche prevalence

For each patient, divide the number of spots in a niche by that patient’s total analyzed
spots. Plot every patient. Statistical tests with very few patients are low-powered and
should accompany effect sizes and confidence intervals rather than replace them.

## Spatial visualization

Join niche labels to Seurat by sample plus barcode. Confirm that the number and order of
matched spots are exact, then visually inspect tissue registration and missing labels.
