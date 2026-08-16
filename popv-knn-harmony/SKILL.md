---
name: popv-knn-harmony
description: Configure and run popV KNN_HARMONY for Harmony batch integration, latent-space KNN cell-type transfer, optional UMAP, and saved FAISS classifier inference. Use for Harmony-based popV annotation, X_pca_harmony_popv, popv_knn_harmony_prediction, or batch-correction diagnostics.
---

# Run KNN_HARMONY

Require a popV-processed AnnData with `X_pca`, batch annotations, and categorical reference labels.

```python
from popv.algorithms import KNN_HARMONY

method = KNN_HARMONY(
    method_kwargs={"dimred": 50},
    classifier_dict={"n_neighbors": 15, "weights": "uniform"},
    embedding_kwargs={"min_dist": 0.1},
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Expect the corrected representation in `adata.obsm["X_pca_harmony_popv"]`, predictions in `adata.obs["popv_knn_harmony_prediction"]`, and UMAP in `adata.obsm["X_umap_harmony_popv"]`. Retraining can save `harmony_knn_classifier`; inference may position query cells from nearest reference PCA neighbors when embeddings are reused. Fast mode is unsupported. Evaluate biological conservation alongside batch mixing; overcorrection can erase real cell-state structure.
