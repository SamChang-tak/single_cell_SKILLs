---
name: popv-knn-scanorama
description: Configure and run popV KNN_SCANORAMA for Scanorama integration followed by KNN cell-type transfer and optional UMAP. Use for Scanorama-based popV annotation, X_pca_scanorama_popv, popv_knn_scanorama_prediction, or legacy integration comparisons.
---

# Run KNN_SCANORAMA

This method is listed by current popV as outdated; use deliberately for compatibility or benchmarking.

```python
from popv.algorithms import KNN_SCANORAMA

method = KNN_SCANORAMA(
    method_kwargs={},
    classifier_kwargs={"n_neighbors": 15, "weights": "uniform"},
    embedding_kwargs={"min_dist": 0.1},
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Expect the integrated representation in `adata.obsm["X_pca_scanorama_popv"]`, predictions in `adata.obs["popv_knn_scanorama_prediction"]`, and UMAP in `adata.obsm["X_umap_scanorama_popv"]`. The implementation sorts cells by batch for integration and restores original order; verify indices remain unique. It builds a joint embedding and is not a pretrained fast-mode method. Assess batch mixing and biological neighborhood preservation before trusting transferred labels.
