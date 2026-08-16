# Early spatial workflow data contracts

## Visium input

Require a 10x-compatible sample directory containing a filtered feature-barcode matrix,
spatial positions, `spatial/scalefactors_json.json`, and tissue image files. Spot barcodes
must match across the matrix, coordinates, images, Seurat object, AnnData, and exported CSVs.

## Single-cell reference

CellTrek uses the integrated Seurat object from the preceding scRNA-seq integration skill.
cell2location uses the integrated AnnData reference with raw counts in `layers['counts']`
and labels in `obs['cell_type']`. Confirm compatible gene identifiers and retain patient,
condition, and study metadata.

## Spatial AnnData

Require unique spot IDs, non-negative counts, and finite two-dimensional coordinates in
`obsm['spatial']`. Add `sample`, `patient_id`, and `condition` before cross-sample summaries.
Document the exact cell2location posterior key/quantile and the rule used to create any hard
`pred_cell_type` label.

## Seurat and CSV handoffs

Require a `Spatial` assay, image object, tissue coordinates, and sample metadata. Join
cell2location CSVs to Seurat metadata by barcode, never by row position. Save a schema with
the posterior statistic, units, cell-type naming, and source checkpoint.

## Interpretation boundaries

- CellTrek outputs mapped reference-cell coordinates.
- cell2location outputs posterior cell-abundance estimates for mixed spatial locations.
- BANKSY outputs spatially informed expression domains.
- Squidpy neighborhood enrichment outputs permutation-based adjacency z-scores.

These objects are related but not interchangeable.
