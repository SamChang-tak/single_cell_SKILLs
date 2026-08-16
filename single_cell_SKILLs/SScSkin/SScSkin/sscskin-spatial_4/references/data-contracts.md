# Spatial data contracts

## Visium source directories

Require a 10x-compatible sample directory containing the filtered feature-barcode matrix,
`spatial/tissue_positions*`, `spatial/scalefactors_json.json`, and tissue image files. Verify
that barcodes align across the matrix, spatial positions, image, Seurat object, and AnnData.

## Integrated scRNA-seq reference

CellTrek expects the integrated Seurat RDS produced by the preceding integration workflow.
cell2location expects the integrated AnnData object with raw counts in `layers['counts']` and
cell labels in `obs['cell_type']`. Retain `patient_id`, `condition`, and study/source fields.

## Spatial AnnData

Require:

- non-negative count matrix or a documented counts layer;
- `obsm['spatial']` with finite two-dimensional coordinates;
- unique spot identifiers;
- `obs['patient_id']`, `obs['condition']`, and `obs['sample']` before integration;
- cell2location abundance columns plus a documented `pred_cell_type` rule when neighborhood
  enrichment is performed on a hard label.

cell2location typically writes posterior abundance summaries into `obsm` and may also export
selected values to `obs`/CSV. Record the exact posterior key and quantile used.

## Seurat objects

Require a `Spatial` assay, image object, tissue coordinates, and sample metadata. After
attaching deconvolution CSVs, validate row order by barcode; never assign columns by position
without an explicit identifier join.

## Cross-language handoffs

- Preserve original barcodes when exporting CSV or moving between AnnData and Seurat.
- Write a schema/manifest beside every CSV: row identifier, units, transform, posterior
  statistic, and cell-type naming convention.
- Validate dimensions and identifiers after every handoff.
- Avoid `as.matrix()` on large sparse counts unless memory requirements are acceptable.
