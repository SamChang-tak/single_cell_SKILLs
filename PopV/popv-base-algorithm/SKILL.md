---
name: popv-base-algorithm
description: Implement, extend, or review custom popV algorithm classes that inherit BaseAlgorithm and operate on processed AnnData. Use when adding a predictor to popv.algorithms, defining integration/prediction/UMAP hooks, managing result keys, or debugging the shared algorithm contract.
---

# Extend BaseAlgorithm

Subclass `popv.algorithms.BaseAlgorithm`. Set stable keys through `super().__init__`, and implement all three hooks: `compute_integration(adata)`, `predict(adata)`, and `compute_umap(adata)`.

```python
from popv.algorithms import BaseAlgorithm

class MyAlgorithm(BaseAlgorithm):
    def __init__(self, result_key="popv_my_prediction"):
        super().__init__(result_key=result_key)

    def compute_integration(self, adata):
        pass

    def predict(self, adata):
        # Write predictions to adata.obs[self.result_key].
        ...

    def compute_umap(self, adata):
        pass
```

Honor global settings copied by the base class: `enable_cuml`, `return_probabilities`, and `compute_umap_embedding`. Write embeddings to `adata.obsm[self.embedding_key]`, UMAP to `adata.obsm[self.umap_key]`, and probabilities to both a maximum-probability `obs` column and a class-aligned `obsm` matrix when supported.

Export the class from `popv/algorithms/__init__.py`; popV discovers exported classes for `AlgorithmsNT`. Ensure it works with the private fields created by `Process_Query`, preserves cell order, handles `retrain`/`inference`/`fast` explicitly, and persists any model required for later inference.
