# Annotation guardrails

## Before fitting

- Confirm the expression matrix is raw, nonnegative counts. Integer-like values are expected for ordinary UMI data.
- Make gene identifiers consistent and unique; resolve symbols versus Ensembl IDs and duplicated symbols explicitly.
- Remove ubiquitous, contradictory, or poorly supported markers only with a documented biological rationale.
- Review marker counts per type and overlaps between types. Highly overlapping marker sets reduce identifiability.
- Include all plausible major populations. CellAssign is closed-set and may force an omitted type into the nearest supplied category.
- Compute size factors from all measured genes, before selecting marker genes.

## After fitting

- Check validation ELBO for stabilization, divergence, or abrupt behavior.
- Inspect maximum probability, top-two margin, and entropy by cell, cluster, donor, batch, and library size.
- Confirm assigned cells express several coherent markers, not only one ambient or stress-associated transcript.
- Investigate mixed profiles as possible doublets, transitional states, low-quality cells, or inadequate marker definitions.
- Compare with unsupervised clusters and trusted orthogonal evidence. Do not optimize markers solely to reproduce an existing annotation.
- Refit after materially changing the marker matrix and preserve each version for provenance.

## Reporting

Report candidate cell types, marker source and version, overlap statistics, missing markers, count layer, size-factor definition, training parameters, convergence evidence, and assignment uncertainty. Separate model assignments from manual biological interpretation.

Do not claim that CellAssign establishes lineage, causality, state transitions, or the absence of a cell type. It provides probabilities conditional on the data, model, and supplied marker categories.
