# Interpretation guardrails

## Study design

- Cell-level scVI DE borrows information across cells but does not create biological replication.
- For disease, treatment, genotype, or other condition claims, inspect donor balance and corroborate with donor-aware pseudobulk when possible.
- Avoid a global condition contrast when condition and batch/donor are confounded. Restrict, stratify, or redesign the comparison.
- For within-cell-type condition DE, subset or define indices so cell identity does not dominate the contrast.

## Batch correction

- Enable only if both groups contain the same batch categories or meaningful overlap.
- If batches are disjoint, correction requires unsupported counterfactual extrapolation; report the confounding instead.
- Record the exact overlap and the effective batch-correction choice.

## Gene universe and effect size

- scVI tests only genes used to train the model. An HVG-only model cannot support claims about untested genes.
- Rank with posterior probability and effect size together. Large LFC for extremely sparse genes can be unstable.
- Inspect raw means and nonzero proportions in both groups.
- Define the sign in every plot and table: positive LFC is group 1 over group 2.

## Statistical language

- `proba_de` and Bayes factors are posterior quantities, not p-values.
- `is_de_fdr_*` is a Bayesian FDR decision conditional on model assumptions.
- Do not convert cell counts into replicate counts or report cell-level evidence as population-level causality.
- Treat pathway enrichment and marker interpretation as downstream hypotheses requiring independent validation.

## Diagnostics

- Require stable validation ELBO and no divergence.
- Check latent UMAP by group and batch for separation, residual batch structure, and confounding.
- Compare importance and uniform weights if conclusions depend on borderline genes.
- Review top hits for mitochondrial, ribosomal, dissociation, ambient-RNA, and low-detection artifacts.
