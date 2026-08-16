---
name: popv-xgboost
description: Configure and run popV XGboost for multiclass single-cell annotation, CPU or GPU training, saved-model inference, and probability outputs. Use for popv_xgboost_prediction, XGBoost parameter tuning, device selection, or boosted-tree troubleshooting in popV.
---

# Run XGboost

Use popV-processed expression with identical genes and order across training and inference.

```python
from popv.algorithms import XGboost

method = XGboost(
    layer_key=None,
    classifier_dict={
        "tree_method": "hist",
        "device": "cpu",
        "objective": "multi:softprob",
    },
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Use `adata.X` by default or a compatible layer. Retraining constructs `DMatrix` objects, sets `num_class`, trains 300 boosting rounds, and saves `xgboost_classifier.model`; inference loads it on CPU. Expect predictions in `adata.obs["popv_xgboost_prediction"]` and optional class-wise probabilities. Verify labels exclude the unknown category from `num_class`, inspect class imbalance and overfitting, and confirm XGBoost/device compatibility before changing tree parameters.
