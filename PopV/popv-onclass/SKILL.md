---
name: popv-onclass
description: Configure and run popV ONCLASS for ontology-aware cell-type prediction of seen and unseen labels. Use for OnClass integration with popV, Cell Ontology files, popv_onclass_prediction, popv_onclass_seen, ontology aggregation, or unseen-cell-type troubleshooting.
---

# Run ONCLASS

Require popV preprocessing with a valid Cell Ontology OBO file and label-to-ontology mapping. Skip this method when `adata.uns["_cl_obo_file"] is False`.

```python
from popv.algorithms import ONCLASS

method = ONCLASS(
    layer_key=None,
    max_iter=30,
    cell_ontology_obs_key="cell_ontology_class",
)
method.compute_integration(adata)
method.predict(adata)
method.compute_umap(adata)
```

Expect ontology-aware predictions in `adata.obs["popv_onclass_prediction"]` and predictions restricted to seen labels in `adata.obs["popv_onclass_seen"]`. Use `layer_key` only when the selected layer has the expected gene order and transform. Retraining persists OnClass artifacts in popV's trained-model path; inference requires matching ontology and feature metadata. Review unmapped labels, obsolete ontology terms, gene overlap, and ontology depth before treating unseen-class predictions as discoveries.
