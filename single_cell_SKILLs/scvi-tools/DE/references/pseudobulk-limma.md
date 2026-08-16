# Donor-aware limma-voom pseudobulk

Use the same cell type × donor raw-count pseudobulks, sample audit, condition
factor levels, and contrast direction as the edgeR analysis.

## Workflow

1. Construct an edgeR `DGEList` and design matrix with control as reference.
2. Run `filterByExpr(y, design=design)`.
3. Run `calcNormFactors(y, method="TMM")`.
4. Run `voom(y, design, plot=FALSE)` to estimate precision weights.
5. Fit `lmFit(v, design)` and moderate with `eBayes(fit, robust=TRUE)`.
6. Test the disease/treatment coefficient and export gene, logFC, AveExpr, t,
   raw P, adjusted P, B statistic, and mean CPM per group.

Use adjusted P < 0.05 and |log2FC| >= 1 for primary FDR-controlled calls. Plot
`-log10(adjusted P)`. If nominal results are explicitly requested, separately
export raw P < 0.05 and |log2FC| >= 1 and plot `-log10(raw P)`. Retain both raw
and adjusted P-values in all full tables, and label nominal outputs exploratory.
