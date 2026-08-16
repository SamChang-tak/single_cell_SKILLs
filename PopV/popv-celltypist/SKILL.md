---
name: popv-celltypist
description: Configure and run popV CELLTYPIST for supervised single-cell cell-type annotation, model retraining or inference, majority voting, and probability outputs. Use for CellTypist within popV, popv_celltypist_prediction, over-clustering behavior, or CellTypist model troubleshooting.
---

# Run CELLTYPIST

Use a popV-processed AnnData with the expected normalized expression representation and categorical reference labels.

```python
from popv.algorithms import CELLTYPIST

method = CELLTYPIST(
    method_kwargs={"n_jobs": 10, "max_iter": 500},
    classifier_dict={"mode": "best match", "majority_voting": True},
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Expect `adata.obs["popv_celltypist_prediction"]`; probability-enabled runs also populate maximum and class-wise probabilities. In `retrain`, save/load `celltypist.pkl` through popV's trained-model path. In `fast`, majority voting is disabled. Other modes may compute or reuse `over_clustering`; verify Leiden dependencies and cluster alignment. For more than 100,000 training cells without GPU support, popV switches CellTypist to SGD mini-batches. Diagnose expression-scale or feature mismatches before tuning iterations.
