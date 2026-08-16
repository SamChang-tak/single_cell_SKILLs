# LDVAE interpretation guardrails

## Data and fitting

- Fit on unlogged counts, using a count layer registered through `setup_anndata`.
- Preserve gene identifiers and the exact HVG set; loading interpretation is conditional on included genes.
- Model known batch structure only when supported by the experimental design. Batch correction can remove biology when batch and condition are confounded.
- Compare training and validation ELBO and record stopping behavior, seed, versions, and hardware.

## Factor identity

- Factor order is arbitrary and factors may permute across fits.
- Factor orientation is arbitrary; an equivalent solution may reverse the sign of a factor and all its loadings.
- Align factors across runs by loading or score correlation, allowing sign flips, rather than by column name.
- Correlated or weakly separated programs can rotate or split when `n_latent`, genes, or regularization changes.

## Biological interpretation

- Examine both loading tails, expression coherence, and factor scores across cells.
- Test associations with donor, batch, count depth, detected genes, mitochondrial percentage, cell cycle, stress, and doublet scores.
- Use pathway enrichment as supporting context, with an explicit background gene universe and multiple-testing correction.
- Validate putative programs across donors and, where possible, independent datasets or orthogonal assays.
- Do not call loading coefficients differential expression, causal regulators, trajectories, or pathway activity without an appropriate downstream analysis.

## Model choice

Choose LinearSCVI when interpretable gene-factor links are central. Compare with standard SCVI when nonlinear reconstruction, integration, or predictive performance is important. A simpler decoder is an interpretability tradeoff, not automatically a better biological model.
