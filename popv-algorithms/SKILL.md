---
name: popv-algorithms
description: Coordinate popV single-cell RNA-seq cell-type annotation, choose among its predictors, configure prediction modes, run annotate_data, and interpret per-method and consensus outputs. Use for end-to-end popV algorithm selection, methods lists, methods_kwargs, consensus voting, or comparison of popV predictors.
---

# Coordinate popV algorithms

Use `popv.preprocessing.Process_Query` before annotation; the algorithm classes expect popV's private `adata.obs` and `adata.uns` fields rather than an arbitrary AnnData object.

## Select methods

- Use `popv.annotation.algorithms_nt.CURRENT_ALGORITHMS` by default.
- Use `FAST_ALGORITHMS` for pretrained CPU-friendly inference: `KNN_SCVI`, `SCANVI_POPV`, `Support_Vector`, `XGboost`, `ONCLASS`, and `CELLTYPIST`.
- Pass `methods="all"` only when explicitly including outdated `Random_Forest` and `KNN_SCANORAMA` is acceptable.
- Omit `ONCLASS` when `adata.uns["_cl_obo_file"] is False`.

```python
from popv.annotation import annotate_data

annotate_data(
    adata,
    methods=["KNN_SCVI", "SCANVI_POPV", "CELLTYPIST"],
    methods_kwargs={"KNN_SCVI": {"classifier_dict": {"n_neighbors": 25}}},
    save_path="popv_results",
)
```

## Check outputs

Inspect `adata.uns["prediction_keys"]`, `adata.uns["methods"]`, per-method columns in `adata.obs`, `popv_majority_vote_prediction`, `popv_majority_vote_score`, `popv_prediction`, and `popv_prediction_score`. Treat agreement as confidence evidence, not calibrated correctness. Preserve the processed AnnData and trained-model directory for reproducible inference.

## Troubleshoot

Verify `_prediction_mode` (`retrain`, `inference`, or `fast`), model paths, categorical label fields, raw counts in `layers["scvi_counts"]` for scVI methods, and compatible reference prediction keys before changing model hyperparameters.
