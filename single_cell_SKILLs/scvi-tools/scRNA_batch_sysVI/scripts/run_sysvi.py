#!/usr/bin/env python3
"""Fit SysVI to preprocessed continuous features and export an integrated embedding."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Preprocessed H5AD")
    parser.add_argument("--output-dir", required=True, type=Path, help="New output directory")
    parser.add_argument("--system-key", required=True, help="Strongest system covariate")
    parser.add_argument("--categorical-covariate", action="append", default=[])
    parser.add_argument("--continuous-covariate", action="append", default=[])
    parser.add_argument("--embed-categorical-covariates", action="store_true")
    parser.add_argument("--cycle-weight", type=float, default=5.0)
    parser.add_argument("--kl-weight", type=float, default=1.0)
    parser.add_argument("--max-epochs", type=int, default=200)
    parser.add_argument("--validation-frequency", type=int, default=1)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plan = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scvi
    from scipy import sparse
    from scvi.external import SysVI

    if args.output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {args.output_dir}")
    if args.max_epochs <= 0 or args.validation_frequency <= 0 or args.n_neighbors <= 0:
        raise ValueError("Epochs, validation frequency, and neighbor count must be positive")
    if args.cycle_weight < 0 or args.kl_weight < 0:
        raise ValueError("Loss weights must be nonnegative")

    adata = ad.read_h5ad(args.input)
    keys = [args.system_key, *args.categorical_covariate, *args.continuous_covariate]
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise KeyError(f"Missing adata.obs columns: {missing}")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell and feature identifiers must be unique")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if values.size and not np.isfinite(values).all():
        raise ValueError("adata.X must contain finite preprocessed features")
    if adata.obs[args.system_key].nunique() < 2:
        raise ValueError("The system key must contain at least two systems")

    SysVI.setup_anndata(
        adata=adata,
        batch_key=args.system_key,
        categorical_covariate_keys=args.categorical_covariate or None,
        continuous_covariate_keys=args.continuous_covariate or None,
    )
    scvi.settings.seed = args.seed
    model = SysVI(
        adata=adata,
        embed_categorical_covariates=args.embed_categorical_covariates,
    )
    model.train(
        max_epochs=args.max_epochs,
        check_val_every_n_epoch=args.validation_frequency,
        plan_kwargs={
            "z_distance_cycle_weight": args.cycle_weight,
            "kl_weight": args.kl_weight,
        },
    )

    latent = model.get_latent_representation(adata=adata)
    columns = [f"sysVI_{i}" for i in range(latent.shape[1])]
    latent_df = pd.DataFrame(latent, index=adata.obs_names, columns=columns)
    embed = ad.AnnData(X=latent, obs=adata.obs.copy())
    sc.pp.neighbors(embed, use_rep="X", n_neighbors=args.n_neighbors)
    sc.tl.umap(embed)
    adata.obsm["X_sysVI"] = latent
    adata.obsm["X_sysVI_umap"] = embed.obsm["X_umap"]

    args.output_dir.mkdir(parents=True)
    latent_df.to_csv(args.output_dir / "sysvi_latent.csv")
    embed.write_h5ad(args.output_dir / "sysvi_embedding.h5ad")
    adata.write_h5ad(args.output_dir / "sysvi_annotated.h5ad")
    model.save(args.output_dir / "model", overwrite=False)
    history = model.trainer.logger.history
    for key, frame in history.items():
        frame.to_csv(args.output_dir / f"history_{key}.csv")

    crosstabs = {}
    for key in args.categorical_covariate:
        table = pd.crosstab(adata.obs[args.system_key], adata.obs[key], dropna=False)
        table.to_csv(args.output_dir / f"system_by_{key}.csv")
        crosstabs[key] = list(table.shape)
    config = {
        **plan,
        "n_cells": adata.n_obs,
        "n_features": adata.n_vars,
        "n_systems": int(adata.obs[args.system_key].nunique()),
        "crosstab_shapes": crosstabs,
        "versions": {
            "anndata": ad.__version__,
            "scanpy": sc.__version__,
            "scvi_tools": scvi.__version__,
        },
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
