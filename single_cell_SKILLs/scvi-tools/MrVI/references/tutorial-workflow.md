# Official MrVI tutorial workflow

Source: scvi-tools 1.3.3, **MrVI Quick Start Tutorial**:
https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/MrVI_tutorial.html

The rendered tutorial reports that it was last run with scvi-tools 1.3.0.

## Intended design

MrVI analyzes multi-sample scRNA-seq with comparable observations across many samples, such as the same tissue or cell line. The tutorial uses 30,000 PBMCs from 16 donors in one cohort of a COVID-19 study:
https://figshare.com/ndownloader/files/46017615

## Preprocessing and fitting

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=10000,
    inplace=True,
    subset=True,
    flavor="seurat_v3",
)

MRVI.setup_anndata(adata, sample_key="patient_id")
model = MRVI(adata)
model.train(max_epochs=400)
```

Use `sample_key` for the target sample/donor identity. Use `batch_key` for a nuisance variable such as site only when correction is appropriate. Inspect `model.history["elbo_validation"]`, ignoring a few initial epochs only for visualization rather than to conceal instability.

## Representations and sample distances

MrVI learns two representations:

- `u`: broad cell state, invariant to sample and nuisance covariates.
- `z`: augments `u` with sample-specific effects while remaining corrected for nuisance covariates.

```python
u = model.get_latent_representation()
adata.obsm["u"] = u
sc.pp.neighbors(adata, use_rep="u")
sc.tl.umap(adata, min_dist=0.3)
```

Compute local sample relationships with:

```python
dists = model.get_local_sample_distances(
    keep_cell=False,
    groupby="initial_clustering",
    batch_size=32,
)
```

`keep_cell=False` reduces memory use. Grouped matrices support cell-population-specific sample clustering and metadata comparison.

## Covariate-linked DE and DA

For a sample-level covariate such as disease status:

```python
de_res = model.differential_expression(
    sample_cov_keys=["Status"],
    store_lfc=True,
)
da_res = model.differential_abundance(sample_cov_keys=["Status"])
```

Category order determines coefficient interpretation for categorical covariates. DE returns cell-resolved effect estimates and gene-level LFC information. DA returns covariate-specific cell-state log probabilities; subtracting two category log probabilities gives a local log-probability ratio whose sign depends on contrast order.
