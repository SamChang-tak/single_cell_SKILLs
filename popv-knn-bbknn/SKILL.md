---
name: popv-knn-bbknn
description: Configure and run popV KNN_BBKNN for cell-type transfer using BBKNN batch-balanced neighbor integration followed by KNN classification and optional UMAP. Use for BBKNN-specific popV annotation, neighbor tuning, batch-mixing diagnostics, or popv_knn_bbknn_prediction outputs.
---

# Run KNN_BBKNN

Require a popV-processed AnnData with PCA in `obsm["X_pca"]`, categorical labels, and a valid batch field.

```python
from popv.algorithms import KNN_BBKNN

method = KNN_BBKNN(
    batch_key="_batch_annotation",
    labels_key="_labels_annotation",
    method_kwargs={"n_pcs": 50, "neighbors_within_batch": 3},
    classifier_kwargs={"n_neighbors": 15, "weights": "uniform"},
    embedding_kwargs={"min_dist": 0.1},
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Expect predictions in `adata.obs["popv_knn_bbknn_prediction"]` and UMAP in `adata.obsm["X_umap_bbknn_popv"]`. BBKNN writes the neighbor graph to `adata.obsp`; the classifier uses precomputed graph distances. popV reduces `neighbors_within_batch` for more than 100 batches and may disable RAPIDS above 200 batches. Check connectivity and per-batch cell counts before interpreting poor predictions. BBKNN is not a pretrained fast-mode method because it builds a joint graph.
