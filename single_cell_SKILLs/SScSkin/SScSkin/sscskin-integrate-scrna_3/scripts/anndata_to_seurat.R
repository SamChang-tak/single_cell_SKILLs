#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(anndataR)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: anndata_to_seurat.R INPUT.h5ad OUTPUT.rds")
}

adata <- read_h5ad(args[[1]])
required_obsm <- c("X_pca", "X_pca_harmony", "X_umap")
missing_obsm <- setdiff(required_obsm, names(adata$obsm))
if (!("counts" %in% names(adata$layers)) || length(missing_obsm) > 0) {
  stop(paste("Missing counts layer or reductions:", paste(missing_obsm, collapse = ", ")))
}

counts <- t(adata$layers[["counts"]])
colnames(counts) <- adata$obs_names
rownames(counts) <- adata$var_names
obj <- CreateSeuratObject(counts = counts, meta.data = adata$obs, assay = "RNA")

add_reduction <- function(obj, matrix, name, key) {
  matrix <- as.matrix(matrix)
  rownames(matrix) <- colnames(obj)
  colnames(matrix) <- paste0(key, seq_len(ncol(matrix)))
  obj[[name]] <- CreateDimReducObject(embeddings = matrix, key = key, assay = DefaultAssay(obj))
  obj
}

obj <- add_reduction(obj, adata$obsm[["X_pca"]], "pca", "PCA_")
obj <- add_reduction(obj, adata$obsm[["X_pca_harmony"]], "harmony", "Harmony_")
obj <- add_reduction(obj, adata$obsm[["X_umap"]], "umap", "UMAP_")
saveRDS(obj, args[[2]])
