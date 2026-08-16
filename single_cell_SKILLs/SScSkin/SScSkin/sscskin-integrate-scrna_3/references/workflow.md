# Upstream integration workflow and reproducibility notes

## Contents

- [Notebook sequence](#notebook-sequence)
- [Input contract](#input-contract)
- [Integration procedure](#integration-procedure)
- [Required outputs](#required-outputs)
- [Known issues](#known-issues)

## Notebook sequence

1. `01_integrate.ipynb`: load annotated Gur 2022 and Tabib 2021 AnnData objects;
   standardize metadata and conditions; concatenate common features; restore counts;
   normalize, calculate HVGs/PCA/neighbors/UMAP; run Harmony by dataset; rebuild the graph
   and UMAP; write `integrated.h5ad`.
2. `02_viz_CD31_CD34_ACTA2.ipynb`: compare `PECAM1`, `CD34`, `ACTA2`, and cell types in
   each source and the integrated atlas.
3. `03_markers.ipynb`: rank broad-cell-type markers with Wilcoxon and display the top five
   positive-log-fold-change genes per type.
4. `04_feature_plot.ipynb`: visualize canonical marker panels for broad cell types.
5. `05_anndata_to_seurat.ipynb`: create a Seurat RNA assay from counts and attach PCA,
   Harmony, and UMAP reductions before saving `integrated.rds`.
6. `06_NAMPT_NNMT.ipynb`: visualize `NAMPT` and `NNMT`, stratify dot plots by condition
   and cell type, and run an exploratory cell-level fibroblast Wilcoxon comparison.
7. `07.ipynb`: create a standardized marker dot plot; the bundled PDF is its saved output.

## Input contract

Gur 2022 input:

- `.obs`: `patient_id`, `condition`, `cell_type`
- `.layers['counts']`: raw counts

Tabib 2021 input:

- `.obs`: `sample`, `condition`, `cell_type`; rename `sample` to `patient_id`
- `.layers['counts']`: raw counts

Both inputs must use compatible gene identifiers. The upstream concatenation uses the common
feature set. Standardized source labels are `Gur2022` and `Tabib2021`; standardized condition
labels are `Healthy` and `SSc`.

## Integration procedure

- Restore joint `X` from the concatenated count layer.
- Normalize each cell to 10,000 counts and apply natural-log `log1p`.
- Select HVGs with `min_mean=0.0125`, `max_mean=3`, and `min_disp=0.5`.
- Calculate PCA and an uncorrected graph/UMAP.
- Run Harmony on the PCA matrix using only the study label `data`.
- Store corrected coordinates in `obsm['X_pca_harmony']`.
- Recompute neighbors with `use_rep='X_pca_harmony'` and compute the final UMAP.

## Required outputs

- Integrated `.h5ad` with counts, normalized `X`, PCA, Harmony, UMAP, neighbors, source,
  patient, condition, and cell type.
- Pre- and post-Harmony embeddings.
- Source/condition/patient counts by cell type.
- Shared and source-exclusive gene/cell-type tables.
- Complete and top marker tables plus requested feature plots.
- Run metadata with package versions, random seed, parameters, and warnings.

## Known issues

- Upstream file paths alternate between numbered and unnumbered result directories.
- Source notebooks suppress warnings and do not record package versions or random seeds.
- `ad.concat` behavior is not stated explicitly; use an explicit inner feature join and unique
  source-aware observation names.
- Harmony corrects only the two-study label. Patient is not used as a batch covariate.
- The custom `kmeans2` callback in the notebook does not set a seed. A deterministic runner
  should set the available Harmony random-state option.
- Several notebooks load older or differently named intermediate checkpoints.
- The Seurat notebook renames all cells to `Cell1...CellN`, discarding original IDs.
- The fibroblast condition test treats cells as independent replicates.
