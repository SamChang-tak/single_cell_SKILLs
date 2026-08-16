# Method and interpretation guardrails

## Reclustering and markers

- Confirm that `X_pca_harmony` exists and that its integration variables do not erase the
  condition signal of interest. Recompute neighbors only after documenting preprocessing.
- Test label stability across neighbor counts, Leiden resolutions, seeds, and patient holdouts.
- The upstream manual removal/relabeling of clusters is an analytical decision, not a neutral
  preprocessing step. Justify it with markers, QC, and patient composition.
- Wilcoxon tests across cells can overstate significance when patient replication is ignored.
  Prefer patient-aware or pseudobulk confirmation for condition-level claims.

## Trajectory and enrichment

- Slingshot/tradeSeq trajectories are geometry- and root-dependent hypotheses; they do not
  establish time, direction, or lineage ancestry in cross-sectional tissue.
- Inspect lineage weights, coverage by patient/condition, convergence, and sensitivity to the
  selected embedding and knots. Do not interpret sparsely populated trajectory tails strongly.
- For GO analysis, state the tested universe, identifier mapping, organism database version,
  ontology, correction method, and gene-selection rule. Collapse redundant terms cautiously.

## Scores and spatial association

- Seurat `AddModuleScore` depends on assay, expression bins, control sampling, gene-set size,
  and dataset composition. UCell is rank-based but still sensitive to signature coverage and
  data quality. Scores from the two methods are not numerically interchangeable.
- Report mapped/total genes for every signature and flag overlapping or nonspecific markers.
- Cell2location proportions and deconvolved state scores are estimated quantities. Their
  correlation may reflect shared abundance, tissue structure, library size, or smoothing.
- Compute spatial correlations within each section, summarize across patients, correct for
  multiple testing, and compare against spatially aware nulls when making localization claims.

## MISTy and regulator activity

- MISTy importance describes predictive contribution conditional on views, features, spatial
  geometry, and model settings. It is not proof of cell–cell communication or causality.
- Fit per tissue section, document view construction and radius, inspect model performance,
  and require consistency across biological replicates.
- Decoupler ULM activity is inferred from the supplied regulator–target network. Record network
  provenance, weights, minimum target count, species mapping, coverage, and input layer.
- Distinguish inferred TF activity from TF transcript abundance and protein activation.

## Minimum reporting

Report sample and patient counts, cells/spots per group, missing inputs, software/resource
versions, seeds, assay/layer, filters, identifier alignment, gene-set coverage, effect sizes,
uncertainty, multiple-testing correction, sensitivity analyses, and all deviations from the
upstream notebooks.
