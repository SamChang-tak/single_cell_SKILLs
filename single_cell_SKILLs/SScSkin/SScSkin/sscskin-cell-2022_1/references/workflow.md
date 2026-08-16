# Upstream workflow and reproducibility notes

## Contents

- [Notebook sequence](#notebook-sequence)
- [Required inputs](#required-inputs)
- [Published filters](#published-filters)
- [Outputs](#outputs)
- [Known issues](#known-issues)

## Notebook sequence

1. `01_loading_data.ipynb`: read GSE195452 annotation and SRA tables, convert each
   transposed count table to AnnData, attach sample metadata, concatenate, discard `_`
   annotations, and write `scrna.h5ad`.
2. `02_filter_by_tissue.ipynb`: attach disease labels from `mmc1.xlsx`; retain SSc and
   controls from skin with CD45+/CD90+ selection; run cell/gene QC; save counts; normalize,
   log-transform, calculate HVGs/PCA/neighbors/UMAP; write `scrna.skin.h5ad`.
3. `03_merge_annotation.ipynb`: collapse detailed immune and stromal annotations into
   broad `cell_type` categories; write `scrna.skin.annotated.h5ad`.
4. `05_viz_genes.ipynb`: rank broad-cell-type markers using Wilcoxon; write a marker CSV;
   display a marker matrix plot and a curated dot plot.

## Required inputs

- `GSE195452_Cell_metadata_v26_anno.txt`: cell annotation table keyed by cell barcode.
- `SraRunTable.txt`: sample fields including `Sample Name`, `selection_marker`,
  `source_name`, `Tissue`, `PATIENT_ID`, and `Organism`.
- `rna/`: tab-separated per-sample count matrices; notebooks transpose each file so cells
  are observations and genes are variables.
- `mmc1.xlsx`: patient disease metadata; the notebook reads `header=1`, promotes its first
  data row to headers, and uses `PID` and `Disease`.

The upstream `patient.csv` bundled as an asset is a cell-level metadata snapshot and cannot
replace the count matrices.

## Published filters

- Conditions: `SSC`, `Control`
- Tissue: `Skin`, `skin`
- Selection: `CD45+`, `CD90+`
- Excluded annotations: `UN`, `GBP1`, `NRXN1`
- Minimum detected genes: 300
- Minimum counts: 500
- Maximum detected genes: strictly below 5,000
- Maximum mitochondrial fraction: strictly below 30%
- Gene filter: detected in at least one retained cell
- Normalization target: 10,000 counts per cell followed by natural-log `log1p`

## Outputs

The bundled runner writes staged checkpoints, UMAP/QC plots, cell-retention tables,
unmatched annotation labels, a complete Wilcoxon marker table, a compact top-marker table,
and run metadata under one explicit output directory.

## Known issues

- Upstream notebooks alternate between `../../results/cell_2022/` and
  `../../results/01_cell_2022/`.
- The annotation-table indexing is implicit in the first notebook. Validate the barcode
  index explicitly before `.loc` joins.
- Sequential regex replacement can be order-sensitive for labels such as `NK` and
  `NK_XCL1`; use exact mapping instead.
- The first notebook derives sample names by splitting filenames at `_`; confirm this rule
  against the actual count filenames.
- The notebooks suppress warnings globally. The reusable runner does not.
- No batch correction or sample-aware differential expression is performed upstream.
- The notebooks use unpinned package versions and do not record random seeds.
