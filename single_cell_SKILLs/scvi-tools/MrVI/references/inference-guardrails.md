# MrVI inference guardrails

## Experimental design

- Use biological samples, not cells, as independent replicates for sample-level exposure inference.
- Tabulate sample covariates before fitting. Identify perfect or near-perfect confounding, repeated measures, unequal group sizes, and samples with few cells.
- Confirm comparable tissue and sampling context. Strong composition differences can eliminate local overlap needed for meaningful comparisons.
- Keep target sample identity separate from nuisance batch identity.

## Local sample distances

- Report the cell-state or grouping over which distances were summarized.
- Inspect cell counts and representation across every sample within each group.
- Avoid hierarchical-clustering narratives unsupported by uncertainty or stability analysis.
- Compare results across seeds, groupings, and feature selections; do not interpret one dendrogram as a phylogeny.

## Differential expression

- Define reference levels before fitting the contrast and state the coefficient direction.
- Use sample covariates measured at the sample level; do not attach cell-level outcomes and call them donor covariates.
- Correct gene-level tests for multiplicity and summarize cell-resolved effects within prespecified populations.
- Separate model-attributed expression changes from causal effects.

## Differential abundance

- State the order of the log-probability contrast and what positive and negative values mean.
- Verify that local cell states are represented across comparison groups.
- Avoid interpreting smooth local DA patterns as discrete cell-type frequency tests without aggregation or validation.
- Compare with sample-level compositional methods when claims concern conventional cell-type proportions.

## Reporting

Report sample counts, cells per sample, sample and nuisance keys, HVG method, count layer, category order, missing metadata, epochs, seed, versions, convergence, effect uncertainty, multiplicity correction, and sensitivity analyses.
