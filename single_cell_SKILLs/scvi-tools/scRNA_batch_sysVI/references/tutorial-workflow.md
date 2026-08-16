# Official sysVI tutorial workflow

Source: scvi-tools 1.3.3, **Integration of scRNA-seq data with substantial batch effects using sysVI**:
https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/sysVI.html

The rendered page reports that the notebook was last run with scvi-tools 1.3.0.

## Intended use

SysVI augments a conditional VAE with a VampPrior and latent cycle-consistency loss. It targets substantial system differences such as cross-species, cross-technology, cell/nucleus, and primary/in-vitro integration.

## Preprocessing contract

SysVI assumes Gaussian feature noise. For scRNA-seq, use fixed-library-size normalized and log-transformed expression, subset to shared informative genes. The paper workflow selected HVGs per system with within-system batches as `batch_key`, then intersected the system-specific HVG sets to obtain about 2,000 shared genes.

The tutorial data are 10,000 mouse/human pancreas cells and 1,768 genes:
https://github.com/theislab/cross_system_integration/raw/main/tutorials/data/mouse-human_pancreas_subset10000.h5ad

## Setup and fitting

```python
from scvi.external import SysVI

SysVI.setup_anndata(
    adata=adata,
    batch_key="system",
    categorical_covariate_keys=["batch"],
)

model = SysVI(adata=adata)
model.train(
    max_epochs=200,
    check_val_every_n_epoch=1,
    plan_kwargs={"z_distance_cycle_weight": 5},
)
```

Use `embed_categorical_covariates=True` at initialization when embedding categorical covariates is preferable to large one-hot encodings.

The tutorial recommends VampPrior plus cycle consistency. Setting `z_distance_cycle_weight=0` disables cycle consistency and yields a vanilla conditional VAE with the chosen prior. Increasing cycle weight strengthens correction; decreasing cycle or KL weight may improve biological preservation. The documented preferred cycle-weight range was usually 2–10, with rare cases up to 50.

## Diagnostics and embedding

Inspect training and validation versions of reconstruction, local KL, and cycle losses. Random train/validation cells may produce similar curves and do not constitute external biological validation.

```python
latent = model.get_latent_representation(adata=adata)
embed = sc.AnnData(latent, obs=adata.obs.copy())
sc.pp.neighbors(embed, use_rep="X")
sc.tl.umap(embed)
```

Plot the embedding by system and biological labels. Plot each system separately with consistent cell-type colors to expose missing populations and apparent overcorrection.
