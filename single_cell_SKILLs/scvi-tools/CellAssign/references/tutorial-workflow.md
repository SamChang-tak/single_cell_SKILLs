# Official tutorial workflow

Source: scvi-tools 1.3.3, **Annotation with CellAssign**:
https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/cellassign_tutorial.html

The rendered 1.3.3 documentation notes that its notebook was last executed with scvi-tools 1.2.0. Check installed API behavior if reproducing it under another version.

## Model contract

CellAssign annotates cells from a gene-by-cell-type marker matrix without labeled training cells. Rows are genes, columns are candidate cell types, and binary values indicate marker membership.

The tutorial demonstrates two datasets:

- Follicular lymphoma AnnData (9,156 cells × 33,694 genes): https://figshare.com/ndownloader/files/27458798
- High-grade serous carcinoma AnnData (4,848 cells × 33,694 genes): https://figshare.com/ndownloader/files/27458822
- HGSC marker matrix: https://figshare.com/ndownloader/files/27458828
- Follicular lymphoma marker matrix: https://figshare.com/ndownloader/files/27458831
- Original data record: https://zenodo.org/records/3372746

## Canonical sequence

```python
import numpy as np
from scvi.external import CellAssign

# Use the full raw-count matrix here, before selecting marker genes.
library_size = np.asarray(adata.X.sum(axis=1)).ravel()
adata.obs["size_factor"] = library_size / library_size.mean()

# Ensure identical genes and ordering in both objects.
bdata = adata[:, marker_matrix.index].copy()

CellAssign.setup_anndata(bdata, size_factor_key="size_factor")
model = CellAssign(bdata, marker_matrix)
model.train()
probabilities = model.predict()
bdata.obs["cellassign_label"] = probabilities.idxmax(axis=1)
```

The size factor stays on the original scale because the model performs its own logarithmic transformation. Computing it after marker-gene subsetting changes its meaning and is incorrect.

## Tutorial checks

- Plot `model.history["elbo_validation"]` to assess optimization.
- Inspect a heatmap of assignment probabilities rather than only hard labels.
- Plot assignments on UMAP.
- When prior annotations exist, compare them with a confusion matrix while recognizing that neither label set is automatically ground truth.

## Version-aware reproduction

Record `scvi-tools`, `anndata`, `scanpy`, Python, and accelerator versions. Fix random seeds and retain the exact marker matrix. Avoid assuming that numerical results will be identical across scvi-tools, PyTorch, CUDA, or hardware versions.
