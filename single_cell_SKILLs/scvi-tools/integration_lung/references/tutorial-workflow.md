# Official tutorial workflow

Source: scvi-tools 1.3.3 documentation, “Atlas-level integration of lung data.” The rendered page
states that its last run used scvi-tools 1.3.2, so record the runtime version.

## Example data contract

- Download: `https://figshare.com/ndownloader/files/52859312`
- Filename: `lung_atlas_preprocessed.h5ad`
- Shape shown by the tutorial: 32,472 cells × 2,000 genes.
- Raw counts: `layers['counts']`.
- Batch key: `obs['batch']`.
- Label key: `obs['cell_type']`.
- Other metadata include dataset, location, patient group, protocol, donor, sampling method, UMI
  and gene counts, and mitochondrial percentage.

## scVI baseline

```python
scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()
adata.obsm["X_scVI"] = model.get_latent_representation()
```

The tutorial uses non-default `n_layers=2`, `n_latent=30`, and negative-binomial likelihood.
Build a neighbors graph from `X_scVI`, run Leiden, then UMAP with `min_dist=0.3`.

## scANVI refinement

```python
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="cell_type",
    unlabeled_category="Unknown",
)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)
```

Initialize scANVI only from an scVI model trained on the exact same AnnData. The unlabeled category
must not collide with a real label. For fully labeled data, the tutorial uses an otherwise unused
name; for partially labeled data, it must identify the actual unlabeled cells.

## Benchmarking

The tutorial benchmarks `X_pca`, `X_scVI`, and `X_scANVI` with:

```python
Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"],
    n_jobs=-1,
)
```

Interpret label conservation metrics (isolated labels, NMI/ARI, silhouette label, cLISI) jointly
with batch metrics (silhouette batch, iLISI, KBET, graph connectivity, PCR). `NaN` or skipped
metrics must be reported, not silently dropped.
