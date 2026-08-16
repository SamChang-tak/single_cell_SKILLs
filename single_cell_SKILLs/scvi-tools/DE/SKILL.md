---
name: den
description: Run single-cell differential expression with donor-aware edgeR or limma-voom pseudobulk for replicated condition comparisons and scvi-tools SCVI as a Bayesian cell-level sensitivity workflow. Use when Codex needs to compare disease, treatment, genotype, or other conditions within cell types; aggregate raw AnnData counts by donor; fit edgeR quasi-likelihood or limma-voom models; export FDR-controlled or explicitly exploratory nominal-P results; train or load scVI; estimate posterior LFC and DE probabilities; or create per-cell-type CSV, volcano, combined-plot, audit, and HTML outputs.
---

# Donor-aware pseudobulk and scVI differential expression

Read `references/interpretation-guardrails.md` before choosing a design or
interpreting results. For scVI, also read `references/tutorial-workflow.md`.
For replicated condition comparisons, read `references/pseudobulk-edger.md`.
For limma-voom, also read `references/pseudobulk-limma.md`.

Use `scripts/run_scvi_de.py` for a targeted scVI comparison. Use
`scripts/run_scvi_de_all_celltypes.py` for one scVI condition contrast within
every annotated cell type. Use `scripts/run_edger_pseudobulk.R` for a prepared
gene-by-sample pseudobulk matrix.

Use `scripts/run_limma_voom_pseudobulk.R` for a prepared pseudobulk matrix and
`scripts/run_limma_all_celltypes.py` to generate FDR and nominal limma tables,
individual volcano plots, combined plots, summaries, and HTML-ready outputs.

Use `scripts/refine_scvi_and_run_pseudobulk.py` to aggregate raw AnnData counts
by donor for every audited cell type, run edgeR, regenerate finite-sampling-
capped scVI plots, and export combined results. Use
`scripts/generate_nominal_pseudobulk.py` only when explicitly requested to add
an exploratory nominal-P result set. Finish combined reports with
`scripts/finalize_scvi_de_report.py`.

## Workflow selection

- Use donor-aware pseudobulk for disease, treatment, genotype, exposure, or
  other sample-level condition effects when donor/sample identifiers exist.
- Use scVI DE for cell-population contrasts, ranking, or sensitivity analysis;
  do not substitute cell-level evidence for biological replication.
- If biological replicates do not exist, report that population-level condition
  DE is not identifiable. Do not treat cells as replicates.
- Keep FDR-controlled and nominal-P outputs separate. Label nominal results
  exploratory in plots, tables, HTML, and conclusions.

## Workflow

1. Inspect AnnData shape, gene identifiers, layers, group, cell-type, donor, and batch counts.
2. Require unique cells and genes and nonnegative integer-like raw counts.
3. Define the contrast and sign explicitly. Positive LFC means group 1 minus group 2.
4. Check biological replication; cells are not independent biological replicates.
5. For sample-level conditions, sum raw counts by cell type and donor and retain a sample audit with donor, condition, cell count, and library size.
6. Filter with edgeR `filterByExpr`, calculate TMM factors, estimate robust negative-binomial dispersions, and test the condition coefficient with `glmQLFit`/`glmQLFTest`.
7. Define primary pseudobulk calls with FDR and a prespecified effect threshold. Default to FDR < 0.05 and |log2FC| >= 1 unless specified otherwise.
8. If requested, export nominal calls separately using raw P < 0.05 and the same effect threshold. Preserve FDR in every nominal table and state that calls are not multiple-testing controlled.
9. For limma-voom, reuse the filtered/TMM-normalized pseudobulks, estimate precision weights with `voom`, fit with `lmFit`, and moderate with robust `eBayes`.
10. For secondary scVI DE, select genes, register raw counts, train with validation ELBO and early stopping, and review convergence.
11. Run scVI change-mode DE. Prefer `weights="importance"` with outlier filtering for FDR-calibrated calls; use `uniform` for ranking or sensitivity comparison.
12. Enable scVI batch correction only when both groups contain overlapping batch levels. Report confounding when no overlap exists.
13. Export full and filtered tables, sample/group audits, model diagnostics, individual and combined plots, configuration, and HTML methodology.

## scVI baseline

```bash
python scripts/run_scvi_de.py \
  --input data.h5ad --output-dir results/scvi_de \
  --groupby cell_type --group1 Type_A --group2 Type_B \
  --count-layer counts --batch-key donor \
  --weights importance --batch-correction auto
```

Omit `--group2` for group 1 versus all remaining cells. Use `--dry-run`
before training when keys or group definitions are uncertain.

## Quality rules

- Prefer at least three biological replicates per condition. Flag any group with fewer than three.
- Use donor/sample pseudobulks as independent units. Do not add a donor blocking term when each donor belongs to only one condition.
- Plot FDR-controlled edgeR results with `-log10(FDR)` and an FDR threshold line.
- Plot nominal results with `-log10(raw P)` and a raw-P threshold line. Never mix axis and coloring criteria.
- Reject empty groups and warn for small groups; do not silently downsample.
- Treat `proba_de`, Bayes factors, LFC, expression prevalence, and FDR calls as complementary evidence.
- Do not call `proba_not_de` a frequentist p-value; label it as a posterior score.
- Cap scVI posterior volcano scores at finite Monte Carlo resolution. With 5,000 samples, use a visualization floor of `proba_not_de=1e-4`, cap at 4, and mark saturated genes. Do not alter saved posterior probabilities.
- Preserve software versions, random seed, trained genes, model parameters, contrast direction, and thresholds.

## Outputs

- Primary pseudobulk: per-cell-type `edger_results.csv`, `edger_significant.csv`, `pseudobulk_counts.csv.gz`, `sample_metadata.csv`, and `volcano.png`.
- Combined pseudobulk: summary CSV, compressed full results, and combined FDR volcano plot.
- Optional nominal pseudobulk: separately named annotated/full and filtered CSVs, individual volcano plots, summary CSV, compressed combined tables, and combined nominal plot.
- limma-voom pseudobulk: full, FDR, and nominal CSVs; per-cell-type FDR and nominal volcano plots; summaries; and combined plots/tables.
- scVI: complete and Bayesian-FDR CSVs, group/batch audit, training history, convergence, latent UMAP, model, processed AnnData, and capped volcano plots.

In HTML, put donor-aware pseudobulk first, optional nominal results second, and
scVI sensitivity results last. State that scVI results cover only genes included
during model training.
