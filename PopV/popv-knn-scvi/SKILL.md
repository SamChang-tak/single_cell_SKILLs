---
name: popv-knn-scvi
description: Configure and run popV KNN_SCVI for scVI latent integration followed by KNN cell-type transfer, saved model query mapping, probabilities, and UMAP. Use for popv_knn_on_scvi_prediction, X_scvi_popv, scVI batch correction, or KNN-on-scVI troubleshooting.
---

# Run KNN_SCVI

Require raw integer counts in `adata.layers["scvi_counts"]`, categorical labels and batches, and aligned genes across reference and query.

```python
from popv.algorithms import KNN_SCVI

method = KNN_SCVI(
    model_kwargs={"n_latent": 20, "n_layers": 3},
    classifier_dict={"n_neighbors": 15, "weights": "uniform"},
    train_kwargs={"max_epochs": 20, "batch_size": 512},
    embedding_kwargs={"min_dist": 0.3},
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Expect latent values in `adata.obsm["X_scvi_popv"]`, predictions in `adata.obs["popv_knn_on_scvi_prediction"]`, and UMAP in `adata.obsm["X_umap_scvi_popv"]`. Retraining can save the `scvi` model; inference maps query data with `SCVI.load_query_data`. Fast mode uses one epoch. Check count integrity, gene order, batch covariates, latent mixing, biological conservation, neighbor count, and accelerator/device settings before tuning the network.
