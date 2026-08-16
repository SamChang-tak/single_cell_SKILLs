# Official tutorial workflow

Source: scvi-tools 1.3.3 documentation, “Integration and label transfer with Tabula Muris.” The
rendered page reports its last run with scvi-tools 1.2.1.

## Tutorial datasets

- 10x droplet H5AD: `https://figshare.com/ndownloader/files/23938934`
- Smart-seq2 FACS H5AD: `https://figshare.com/ndownloader/files/23939711`
- Gene-length table: `https://raw.githubusercontent.com/chenlingantelope/HarmonizationSCANVI/master/data/gene_len.txt`
- Subset both datasets to female, labeled marrow cells.
- Tutorial subset shapes before gene harmonization: 11,707 × 20,138 (10x) and 4,200 × 22,966
  (Smart-seq2).
- Technology key: `tech`, with values `10x` and `SS2`.
- Biological label: `cell_ontology_class`.

## Technology-aware preprocessing

Smart-seq2 read counts are divided by average gene length, multiplied by the median gene length,
and rounded. Do not apply this correction blindly to UMI data. Concatenate datasets, copy model
input into `layers['counts']`, normalize/log only for conventional expression views, preserve full
dimension in `.raw`, then select 2,000 HVGs using Seurat v3 on counts with `batch_key='tech'`.

## scVI and scANVI

```python
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="tech")
scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
scvi_model.train()
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
```

The tutorial treats Smart-seq2 as labeled reference and 10x as unlabeled query:

```python
adata.obs["celltype_scanvi"] = "Unknown"
mask = adata.obs["tech"] == "SS2"
adata.obs.loc[mask, "celltype_scanvi"] = adata.obs.loc[mask, "cell_ontology_class"]
scanvi = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    labels_key="celltype_scanvi",
    unlabeled_category="Unknown",
)
scanvi.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = scanvi.get_latent_representation(adata)
adata.obs["C_scANVI"] = scanvi.predict(adata)
```

Build neighbors/UMAP from each latent representation with `min_dist=0.3`. The tutorial evaluates
predictions against known 10x annotations with a row-normalized confusion matrix and notes that
reference/query label resolutions differ.
