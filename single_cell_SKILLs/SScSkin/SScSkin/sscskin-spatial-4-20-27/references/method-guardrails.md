# Functional, molecular-niche, and CellChat guardrails

## decoupler

Record the exact gene-set or regulator resource, retrieval date/version, organism, scoring
method, minimum targets, weight handling, gene-ID mapping, and expression transform. Report
target coverage per activity. Validate key activities by showing contributing target genes.
Do not equate an inferred score with measured signaling flux or regulator protein activity.

## Molecular niches

Expression-based niches and cell-composition niches are different constructs. Record the
input feature space and never reuse labels as if they were equivalent. Harmony correction by
patient can remove patient-linked biology; show pre/post embeddings and perform sensitivity
analysis. Check stability across seeds, neighbor counts, resolutions, and sample leave-one-out
runs. Align labels by markers before comparison.

## CellChat

Require correct spatial coordinates and physical scale factors. Record CellChat version,
database version/subset, group labels, minimum cells, `computeCommunProb` settings, distance
parameters, and filtering. Small groups and differing tissue area or spot counts affect network
density. Compare per-patient normalized summaries rather than raw edge totals alone.

CellChat identifies expression-supported ligand–receptor hypotheses conditional on its model
and database. Spatial proximity and expression do not establish binding, signal direction,
functional response, or causality.

## Condition comparisons

Use patients as replicates. Present patient-level distributions, effect sizes, and uncertainty.
If condition is confounded with batch, tissue region, or cell-type composition, report that the
effect cannot be uniquely attributed to disease.
