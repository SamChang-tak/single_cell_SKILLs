---
name: popv-random-forest
description: Configure and run popV Random_Forest for supervised cell-type annotation from expression features, including retraining, model persistence, and probabilities. Use for popv_rf_prediction, random-forest feature or class-weight tuning, or legacy popV classifier comparisons.
---

# Run Random_Forest

This method is listed by current popV as outdated; include it deliberately for compatibility or benchmarking.

```python
from popv.algorithms import Random_Forest

method = Random_Forest(
    layer_key=None,
    classifier_dict={
        "class_weight": "balanced_subsample",
        "max_features": 200,
        "n_jobs": 8,
    },
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Use `adata.X` by default or a compatible expression layer. In `retrain`, popV fits scikit-learn on `_ref_subsample` and saves `rf_classifier.joblib`; inference loads it. Expect predictions in `adata.obs["popv_rf_prediction"]` and optional probabilities. Keep gene order identical between training and inference. Diagnose sparse/dense memory use, class imbalance, overly restrictive `max_features`, and model-path mismatches before changing the estimator.
