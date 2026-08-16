# CellChat guardrails

- Verify raw counts, gene symbols, complex subunits, metadata, cell identities, lineage labels,
  and group sizes before inference. Never merge matrices by position without proving identical
  gene order.
- Record CellChat/CellChatDB versions and database checksum. Document `subsetDB`, overexpression
  functions, `computeCommunProb(type = "triMean")`, population-size setting, trimming, and
  `filterCommunication(min.cells = 10)`.
- Fit networks separately by patient or tissue section. Compare only after harmonizing subtype
  labels, preprocessing, database, filters, and model parameters.
- Confirm ligand expression in senders and every required receptor/cofactor subunit in receivers.
  Database presence alone is not expression or signaling evidence.
- Use subtype names instead of numeric indices and confirm sender/receiver direction after all
  table reshaping and visualization selections.
- Treat count, probability/weight, pathway aggregation, and gene expression as distinct outputs.
- Assess sensitivity to group size, estimator, downsampling, expression thresholds, database
  subset, rare states, and patient exclusion.
- Report replicate-level effects and uncertainty, missing genes, multiple comparisons, group
  coverage, convergence/warnings, and external validation. Avoid causal wording and cell-level
  pseudoreplication.
