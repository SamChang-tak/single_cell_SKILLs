---
name: cellrank
description: Analyze single-cell state transitions and differentiation trajectories with CellRank and AnnData. Use when Codex needs to select or combine CellRank kernels, compute transition matrices, identify initial or terminal states, calculate fate probabilities or lineage-driver genes, visualize lineage-specific expression trends, or persist and reload CellRank results.
---

# CellRank

Build reproducible CellRank analyses around the separation between kernels, which estimate cell-cell transitions, and estimators, which interpret those transitions.

## Use bundled resources

- Read [references/getting-started-workflow.md](references/getting-started-workflow.md) when implementing or reviewing the official introductory workflow.
- Run `scripts/run_pseudotime_workflow.py` for a deterministic pseudotime-plus-connectivity analysis when its assumptions match the data. Inspect `--help` first and do not use `--preprocess` on already normalized data.
- Copy [assets/analysis-report-template.md](assets/analysis-report-template.md) into the output directory and replace every `{{PLACEHOLDER}}` when producing a formal report.

## Establish inputs

- Inspect the `AnnData` object before analysis: report dimensions and relevant keys in `.obs`, `.obsm`, `.obsp`, `.layers`, and `.uns`.
- Preserve the original object or save a checkpoint before mutating it.
- Confirm that the chosen kernel's required signal exists. Examples include pseudotime, RNA velocity, experimental time, developmental potential, or a precomputed transition matrix.
- Ensure that preprocessing matches the kernel. For a basic expression workflow, filter genes, normalize counts, log-transform, select highly variable genes, compute PCA, and build a neighbor graph. Do not repeat normalization or overwrite a biologically appropriate existing graph without justification.
- Set random seeds for reproducibility and record package versions, parameters, input keys, and output paths.
- Check the installed CellRank API when versions may differ. Do not assume tutorial defaults or estimator recommendations remain current.

## Choose and compute a kernel

Choose the kernel from the biological evidence available:

- Use `PseudotimeKernel` when a defensible pseudotime ordering exists.
- Use `VelocityKernel` for suitably processed spliced/unspliced RNA velocity data.
- Use `CytoTRACEKernel` for developmental-potential information.
- Use `RealTimeKernel` for experimental or spatial time-course data.
- Use `PrecomputedKernel` when transitions come from another method.
- Use `ConnectivityKernel` as an expression-similarity component, commonly combined with a directional kernel rather than treated as directional evidence by itself.

Compute the transition matrix explicitly and inspect the kernel representation and matrix properties. Plot random walks or a projection-based diagnostic when useful. Combine evidence with a justified weighted sum, for example:

```python
import cellrank as cr

pk = cr.kernels.PseudotimeKernel(
    adata, time_key="pseudotime"
).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
kernel = 0.8 * pk + 0.2 * ck
```

Explain kernel weights as modeling assumptions and assess sensitivity to plausible alternatives when conclusions depend on them.

## Estimate macrostates and fates

Use an estimator supported by the installed CellRank version. GPCCA provides the tutorial's general workflow:

```python
g = cr.estimators.GPCCA(kernel)
g.compute_schur()
g.plot_spectrum(real_only=True)
g.compute_macrostates(n_states=10, cluster_key="clusters")
g.plot_macrostates(which="all", basis="umap")
```

- Choose `n_states` from the spectrum, stability, biological resolution, and sensitivity analysis; do not accept a single arbitrary value silently.
- Predict or manually set initial and terminal states using biological context. Plot them and verify that they are coherent with known markers, sampling design, and directionality.
- Compute fate probabilities only after terminal states are established:

```python
g.compute_fate_probabilities()
g.plot_fate_probabilities(basis="umap", legend_loc="right")
```

- Treat fate probabilities as model-based probabilities conditional on the transition model, not direct lineage-tracing observations.
- Compute lineage drivers for explicit lineages and report correlation, effect direction, multiple-testing-adjusted significance, and biological caveats:

```python
drivers = g.compute_lineage_drivers(lineages="target_lineage")
```

## Visualize lineage programs

Fit lineage-weighted gene trends using a model available in the installed environment. Prefer the Python `GAM` when an R dependency is unnecessary; use `GAMR` only when its R and `rpy2` requirements are available.

```python
model = cr.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    genes=["GENE1", "GENE2"],
    time_key="pseudotime",
    same_plot=True,
)
```

Use heatmaps to compare multiple genes or lineages. Select genes from a stated hypothesis or a reproducible driver-ranking rule rather than only visual preference.

## Persist and verify

- Store the kernel in `AnnData` with `write_to_adata()` before writing an `.h5ad` checkpoint.
- Reload through the matching `from_adata()` constructor and verify the expected transition key; the tutorial commonly uses `T_fwd`, but the actual key can vary.
- Validate transition-matrix dimensions, finite values, non-negativity, and row normalization.
- Confirm that fate probabilities are finite, bounded, and sum approximately to one across terminal lineages for each cell.
- Save tables and plots with explicit filenames and return a concise methods summary containing preprocessing, kernel type and parameters, weights, estimator settings, state definitions, lineage names, software versions, and warnings.

## Guardrails

- Never infer directionality from cluster order alone.
- Do not interpret correlation-based lineage drivers as causal regulators.
- Flag disconnected neighbor graphs, sparse terminal populations, weak or conflicting direction signals, unstable macrostate choices, and batch-confounded trajectories.
- Prefer current official CellRank API documentation for exact signatures. Use the official getting-started tutorial as the conceptual baseline: <https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/general/100_getting_started.html>.
