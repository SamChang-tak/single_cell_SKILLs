---
name: popv-support-vector
description: Configure and run popV Support_Vector for linear SVM cell-type annotation, calibrated probabilities, retraining, and saved-model inference. Use for popv_svm_prediction, LinearSVC hyperparameters, class balancing, or SVM performance and convergence troubleshooting in popV.
---

# Run Support_Vector

Use popV-processed expression with identical feature order between training and inference.

```python
from popv.algorithms import Support_Vector

method = Support_Vector(
    layer_key=None,
    classifier_dict={"C": 1, "max_iter": 5000, "class_weight": "balanced"},
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Use `adata.X` by default or select a compatible layer. Retraining fits on `_ref_subsample`, calibrates LinearSVC for probabilities, and saves `svm_classifier.joblib`; inference loads it. Expect predictions in `adata.obs["popv_svm_prediction"]` and optional probabilities. The CPU path densifies training data, so estimate memory first. Address feature scaling, class imbalance, convergence warnings, absent classes, and gene-order mismatches before increasing `max_iter` or changing `C`.
