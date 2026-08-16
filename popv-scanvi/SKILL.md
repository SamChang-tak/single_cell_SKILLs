---
name: popv-scanvi
description: Configure and run popV SCANVI_POPV for semi-supervised scANVI integration and cell-type label transfer, including retraining, query mapping, latent embeddings, probabilities, and UMAP. Use for popv_scanvi_prediction, X_scanvi_popv, scvi-tools model settings, or scANVI troubleshooting in popV.
---

# Run SCANVI_POPV

Require raw integer counts in `adata.layers["scvi_counts"]`, categorical batch and label fields, and the unknown label category configured by popV.

```python
from popv.algorithms import SCANVI_POPV

method = SCANVI_POPV(
    model_kwargs={"n_latent": 20, "n_layers": 3},
    classifier_kwargs={"n_layers": 3, "dropout_rate": 0.1},
    train_kwargs={"max_epochs": 20, "batch_size": 512},
    embedding_kwargs={"min_dist": 0.3},
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Expect latent values in `adata.obsm["X_scanvi_popv"]`, predictions in `adata.obs["popv_scanvi_prediction"]`, and UMAP in `adata.obsm["X_umap_scanvi_popv"]`. Retraining initializes from a saved scVI model when available, then saves `scanvi`; inference uses `SCANVI.load_query_data` and freezes the classifier. Fast mode limits training to one epoch. Verify gene order, count integrity, unknown-label spelling, category order, device selection, and model compatibility before tuning architecture.
