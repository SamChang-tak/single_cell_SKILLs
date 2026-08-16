# CellChat guardrails

- Confirm raw counts, gene symbols, complex subunits, metadata, cell identities, and group sizes
  before inference. Never combine matrices by position without proving identical gene order.
- Record CellChat/CellChatDB versions and checksum the database export. Results can change with
  database content and package implementation.
- Document `subsetDB`, `identifyOverExpressedGenes`, `identifyOverExpressedInteractions`,
  `computeCommunProb(type = "triMean")`, `population.size`, distance settings, trimming, and
  `filterCommunication(min.cells = 10)`.
- Do not infer across unrelated patients or tissue sections. Fit separate networks for each
  biological replicate; compare only after harmonizing labels, filtering, normalization, and
  database settings.
- Check ligand and every required receptor/cofactor subunit in the appropriate source/target
  groups. A database match without expression support is not biological evidence.
- Distinguish edge count, communication probability/weight, pathway aggregate, and gene
  expression. They answer different questions and are not interchangeable.
- Use subtype names rather than numeric indices in visualization and export code. Validate that
  sender/receiver direction is preserved after reshaping tables.
- Assess sensitivity to minimum group size, probability estimator, downsampling, rare states,
  expression thresholds, and patient exclusion.
- Report replicate-level effect sizes and uncertainty. Avoid cell-level pseudoreplication and
  causal wording; validate prioritized interactions with independent spatial/protein evidence.
