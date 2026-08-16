#!/usr/bin/env python3
"""Run de novo scVI/scANVI integration for labeled reference and unlabeled query H5AD files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference", type=Path)
    p.add_argument("query", type=Path)
    p.add_argument("output_dir", type=Path)
    p.add_argument("--reference-label-key", required=True)
    p.add_argument("--count-layer", help="Count layer in both inputs; default uses X")
    p.add_argument("--batch-key", default="tech")
    p.add_argument("--reference-batch", default="reference")
    p.add_argument("--query-batch", default="query")
    p.add_argument("--unlabeled-category", default="Unknown")
    p.add_argument("--n-top-genes", type=int, default=2000)
    p.add_argument("--n-layers", type=int, default=2)
    p.add_argument("--n-latent", type=int, default=30)
    p.add_argument("--scanvi-epochs", type=int, default=20)
    p.add_argument("--samples-per-label", type=int, default=100)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> int:
    a = args()
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    if a.dry_run:
        print(json.dumps(config, indent=2))
        return 0
    if a.output_dir.exists():
        raise SystemExit(f"output directory already exists: {a.output_dir}")
    if not a.reference.is_file() or not a.query.is_file():
        raise SystemExit("reference and query H5AD files must exist")

    import anndata as ad
    import scanpy as sc
    import scvi

    ref, qry = sc.read_h5ad(a.reference), sc.read_h5ad(a.query)
    for name, obj in (("reference", ref), ("query", qry)):
        if not obj.obs_names.is_unique or not obj.var_names.is_unique:
            raise SystemExit(f"{name} cell and gene identifiers must be unique")
        if a.count_layer and a.count_layer not in obj.layers:
            raise SystemExit(f"{name} lacks count layer {a.count_layer!r}")
    if a.reference_label_key not in ref.obs:
        raise SystemExit(f"reference lacks label key {a.reference_label_key!r}")
    if ref.obs[a.reference_label_key].isna().any():
        raise SystemExit("reference labels contain missing values")
    if a.unlabeled_category in set(ref.obs[a.reference_label_key].astype(str)):
        raise SystemExit("unlabeled category collides with a reference label")

    ref.obs[a.batch_key], qry.obs[a.batch_key] = a.reference_batch, a.query_batch
    ref.obs["_scanvi_label"] = ref.obs[a.reference_label_key].astype(str)
    qry.obs["_scanvi_label"] = a.unlabeled_category
    if a.count_layer:
        ref.X, qry.X = ref.layers[a.count_layer].copy(), qry.layers[a.count_layer].copy()
    data = ad.concat([ref, qry], join="inner", label="_source", keys=["reference", "query"])
    data.layers["counts"] = data.X.copy()
    sc.pp.highly_variable_genes(
        data,
        flavor="seurat_v3",
        n_top_genes=a.n_top_genes,
        layer="counts",
        batch_key=a.batch_key,
        subset=True,
    )

    a.output_dir.mkdir(parents=True)
    scvi.settings.seed = a.seed
    scvi.model.SCVI.setup_anndata(data, layer="counts", batch_key=a.batch_key)
    vae = scvi.model.SCVI(data, n_layers=a.n_layers, n_latent=a.n_latent)
    vae.train()
    data.obsm["X_scVI"] = vae.get_latent_representation()
    vae.save(a.output_dir / "scvi_model")
    scanvi = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=data,
        labels_key="_scanvi_label",
        unlabeled_category=a.unlabeled_category,
    )
    scanvi.train(max_epochs=a.scanvi_epochs, n_samples_per_label=a.samples_per_label)
    data.obsm["X_scANVI"] = scanvi.get_latent_representation(data)
    data.obs["C_scANVI"] = scanvi.predict(data)
    scanvi.save(a.output_dir / "scanvi_model")
    data.write_h5ad(a.output_dir / "integrated_label_transfer.h5ad", compression="gzip")
    config["scvi_tools_version"] = scvi.__version__
    (a.output_dir / "run_config.json").write_text(json.dumps(config, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
