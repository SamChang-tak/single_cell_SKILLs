# Official tutorial workflow

Source: scvi-tools 1.3.3, **Linearly decoded VAE**:
https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/linear_decoder.html

The rendered 1.3.3 page reports that the notebook was last run with scvi-tools 1.1.6. Confirm APIs when reproducing it under another release.

## Model idea

Standard scVI uses nonlinear neural-network decoders. LDVAE replaces that decoder with a linear function, directly connecting cell latent coordinates with gene weights. It resembles probabilistic PCA or factor analysis while modeling counts and retaining scalable variational inference.

## Tutorial sequence

The example uses PBMC 10K data (6,855 cells × 16,727 genes) from:
https://github.com/YosefLab/scVI-data/raw/master/pbmc_10k_protein_v3.h5ad?raw=true

```python
adata.layers["counts"] = adata.X.copy()

# Normalization below is for other Scanpy uses; model fitting uses counts.
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", layer="counts", n_top_genes=1000, subset=True
)

scvi.model.LinearSCVI.setup_anndata(adata, layer="counts")
model = scvi.model.LinearSCVI(adata, n_latent=10)
model.train(
    max_epochs=250,
    plan_kwargs={"lr": 5e-3},
    check_val_every_n_epoch=10,
)
```

Inspect convergence using `model.history["elbo_train"]` and `model.history["elbo_validation"]`.

## Extract results

```python
latent = model.get_latent_representation()
loadings = model.get_loadings()
adata.obsm["X_scVI"] = latent
```

For each loading column, sort genes and inspect both tails. Large positive and negative coefficients identify genes whose modeled expression varies in opposite directions along that factor.

The tutorial builds a Scanpy neighbor graph from `X_scVI`, calculates UMAP and Leiden clusters, and colors the UMAP by each latent coordinate.

## Critical tutorial caveat

The factors are not ordered by explained variance as PCA components are. A lower factor number carries no automatic priority.
