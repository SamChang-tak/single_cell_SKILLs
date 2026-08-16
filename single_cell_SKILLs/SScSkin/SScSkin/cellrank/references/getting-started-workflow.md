# CellRank getting-started workflow

Source: <https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/general/100_getting_started.html>

Read this reference when implementing or reviewing a CellRank analysis based on the official introductory workflow.

## Data and preprocessing

Use an `AnnData` object containing the biological signal required by the selected kernel. The tutorial filters genes detected in at least five cells, normalizes each cell to 10,000 counts, log-transforms, selects 3,000 highly variable genes, computes PCA, and builds a nearest-neighbor graph. Treat those settings as an example, not universal defaults.

## Kernel sequence

1. Construct a directional kernel, such as `PseudotimeKernel(adata, time_key=...)`.
2. Call `compute_transition_matrix()`.
3. Inspect random walks when they provide a useful directionality diagnostic.
4. Optionally compute a `ConnectivityKernel` and combine it with a weighted sum. The tutorial illustrates `0.8 * pseudotime + 0.2 * connectivity`.
5. Persist a kernel with `write_to_adata()`. Reload it with the matching `from_adata()` constructor and verify the transition key, commonly `T_fwd`.

Select `VelocityKernel`, `CytoTRACEKernel`, `RealTimeKernel`, or `PrecomputedKernel` when the corresponding modality is the defensible source of transitions.

## Estimator sequence

1. Initialize an estimator from the computed kernel; the tutorial demonstrates `GPCCA`.
2. Compute the Schur decomposition and inspect the spectrum.
3. Compute macrostates with a justified state count.
4. predict or set initial and terminal states, then verify them biologically.
5. Compute fate probabilities.
6. Compute lineage drivers for named lineages.
7. Fit lineage-weighted expression trends with `GAM`, or `GAMR` when its R dependencies are available.

## Interpretation constraints

- Treat transition direction, kernel weights, macrostate count, and terminal-state selection as model assumptions.
- Treat fate probabilities as conditional on the fitted transition model.
- Treat lineage-driver correlations as associations, not causal evidence.
- Check current API documentation because signatures and recommendations can change across CellRank versions.
