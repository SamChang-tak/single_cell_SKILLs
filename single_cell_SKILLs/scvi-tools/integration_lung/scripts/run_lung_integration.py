#!/usr/bin/env python3
"""Train scVI and optional scANVI on an AnnData count layer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path, help="Input .h5ad")
    p.add_argument("output_dir", type=Path, help="New output directory")
    p.add_argument("--count-layer", default="counts")
    p.add_argument("--batch-key", default="batch")
    p.add_argument("--labels-key", help="Enable scANVI using this obs column")
    p.add_argument("--unlabeled-category", default="Unknown")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-layers", type=int, default=2)
    p.add_argument("--n-latent", type=int, default=30)
    p.add_argument("--scanvi-epochs", type=int, default=20)
    p.add_argument("--samples-per-label", type=int, default=100)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> int:
    a = parse_args()
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    if a.dry_run:
        print(json.dumps(config, indent=2))
        return 0
    if not a.input.is_file():
        raise SystemExit(f"input does not exist: {a.input}")
    if a.output_dir.exists():
        raise SystemExit(f"output directory already exists: {a.output_dir}")

    import scanpy as sc
    import scvi

    adata = sc.read_h5ad(a.input)
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise SystemExit("cell and gene identifiers must be unique")
    if a.count_layer not in adata.layers:
        raise SystemExit(f"missing count layer: {a.count_layer}")
    if a.batch_key not in adata.obs:
        raise SystemExit(f"missing batch key: {a.batch_key}")
    if adata.obs[a.batch_key].isna().any():
        raise SystemExit("batch key contains missing values")
    if a.labels_key and a.labels_key not in adata.obs:
        raise SystemExit(f"missing labels key: {a.labels_key}")
    if a.labels_key and a.unlabeled_category in set(adata.obs[a.labels_key].dropna()):
        print(f"note: {a.unlabeled_category!r} occurs in labels and will be treated as unlabeled")

    a.output_dir.mkdir(parents=True)
    scvi.settings.seed = a.seed
    scvi.model.SCVI.setup_anndata(adata, layer=a.count_layer, batch_key=a.batch_key)
    model = scvi.model.SCVI(
        adata, n_layers=a.n_layers, n_latent=a.n_latent, gene_likelihood="nb"
    )
    model.train()
    adata.obsm["X_scVI"] = model.get_latent_representation()
    model.save(a.output_dir / "scvi_model", overwrite=False)

    if a.labels_key:
        scanvi = scvi.model.SCANVI.from_scvi_model(
            model,
            adata=adata,
            labels_key=a.labels_key,
            unlabeled_category=a.unlabeled_category,
        )
        scanvi.train(max_epochs=a.scanvi_epochs, n_samples_per_label=a.samples_per_label)
        adata.obsm["X_scANVI"] = scanvi.get_latent_representation(adata)
        scanvi.save(a.output_dir / "scanvi_model", overwrite=False)

    adata.write_h5ad(a.output_dir / "integrated.h5ad", compression="gzip")
    config["scvi_tools_version"] = scvi.__version__
    (a.output_dir / "run_config.json").write_text(json.dumps(config, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
