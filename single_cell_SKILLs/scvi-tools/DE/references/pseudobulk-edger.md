# Donor-aware edgeR pseudobulk

## Design

For sample-level condition effects, aggregate raw integer counts separately for
each cell type and biological donor. Use one column per donor pseudobulk and one
metadata row containing sample ID, donor, condition, cell count, and library
size. Omit cell types lacking usable samples in both groups and flag any group
with fewer than three donors.

Do not add donor as a blocking factor when donors are nested within condition.
The condition coefficient compares independent donor pseudobulks. Add paired or
repeated-measures terms only when the same biological unit occurs in multiple
conditions.

## edgeR quasi-likelihood workflow

1. Construct `DGEList` from the gene-by-pseudobulk count matrix.
2. Set control as the factor reference and build `~ condition`.
3. Run `filterByExpr(y, design=design)` and drop filtered genes.
4. Run `calcNormFactors(y, method="TMM")`.
5. Run `estimateDisp(y, design, robust=TRUE)`.
6. Fit `glmQLFit(y, design, robust=TRUE)`.
7. Test the disease/treatment coefficient with `glmQLFTest`.
8. Export gene, logFC, logCPM, F statistic, raw P, FDR, and mean CPM per group.

Positive logFC must mean group 1 minus group 2 and must be stated in plots,
tables, and HTML.

## Calls and plots

Use FDR < 0.05 and |log2FC| >= 1 as the default primary call unless another
threshold is prespecified. Plot `-log10(FDR)` so the horizontal line matches the
coloring rule.

If explicitly requested, create a separate exploratory result tree using raw P
< 0.05 and the same effect threshold. Plot `-log10(raw P)`, retain FDR in every
table, and label results unadjusted and not multiple-testing controlled. Never
overwrite or relabel primary FDR results.

## Interpretation

- Treat donors, not cells, as biological replicates.
- Inspect library sizes, cell contributions, and sample-level outliers.
- Treat low-replication cell types as unstable even if nominal P-values are low.
- Report condition/donor or batch confounding the design cannot resolve.
- Use nominal calls for ranking and hypothesis generation only.
