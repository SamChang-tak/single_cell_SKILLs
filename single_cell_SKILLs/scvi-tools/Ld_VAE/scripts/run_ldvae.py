#!/usr/bin/env python3
"""Fit LinearSCVI and export latent coordinates, gene loadings, and diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input H5AD")
    parser.add_argument("--output-dir", required=True, type=Path, help="New output directory")
    parser.add_argument("--count-layer", help="Raw-count layer; default uses adata.X")
    parser.add_argument("--batch-key", help="Optional categorical obs column")
    parser.add_argument("--select-hvg", type=int, help="Select this many seurat_v3 HVGs")
    parser.add_argument("--n-latent", type=int, default=10)
    parser.add_argument("--max-epochs", type=int, default=250)
    parser.add_argument("--learning-rate", type=float, default=5e-3)
    parser.add_argument("--validation-frequency", type=int, default=10)
    parser.add_argument("--top-genes", type=int, default=20)
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

    if args.output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {args.output_dir}")
    for name in ("n_latent", "max_epochs", "validation_frequency", "top_genes"):
        if getattr(args, name) <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if args.learning_rate <= 0 or (args.select_hvg is not None and args.select_hvg <= 0):
        raise ValueError("Learning rate and requested HVG count must be positive")

    adata = ad.read_h5ad(args.input)
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell and gene identifiers must be unique")
    if args.batch_key and args.batch_key not in adata.obs:
        raise KeyError(f"Batch key not found in adata.obs: {args.batch_key}")
    if args.count_layer:
        if args.count_layer not in adata.layers:
            raise KeyError(f"Count layer not found: {args.count_layer}")
        counts = adata.layers[args.count_layer]
    else:
        counts = adata.X
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if values.size and (not np.isfinite(values).all() or values.min() < 0):
        raise ValueError("Counts must be finite and nonnegative")

    work = adata.copy()
    work.layers["_ldvae_counts"] = counts.copy()
    if args.select_hvg and args.select_hvg < work.n_vars:
        sc.pp.highly_variable_genes(
            work,
            flavor="seurat_v3",
            layer="_ldvae_counts",
            n_top_genes=args.select_hvg,
            subset=True,
        )

    scvi.settings.seed = args.seed
    scvi.model.LinearSCVI.setup_anndata(
        work, layer="_ldvae_counts", batch_key=args.batch_key
    )
    model = scvi.model.LinearSCVI(work, n_latent=args.n_latent)
    model.train(
        max_epochs=args.max_epochs,
        plan_kwargs={"lr": args.learning_rate},
        check_val_every_n_epoch=args.validation_frequency,
    )

    latent = model.get_latent_representation()
    loadings = model.get_loadings()
    latent_columns = [f"Z_{i}" for i in range(latent.shape[1])]
    latent_df = pd.DataFrame(latent, index=work.obs_names, columns=latent_columns)
    loadings.columns = latent_columns

    rows = []
    for factor in loadings.columns:
        ranked = loadings[factor].sort_values()
        n = min(args.top_genes, len(ranked))
        rows.extend(
            {"factor": factor, "direction": "negative", "gene": gene, "loading": ranked[gene]}
            for gene in ranked.head(n).index
        )
        rows.extend(
            {"factor": factor, "direction": "positive", "gene": gene, "loading": ranked[gene]}
            for gene in ranked.tail(n).sort_values(ascending=False).index
        )
    top_loadings = pd.DataFrame(rows)

    args.output_dir.mkdir(parents=True)
    latent_df.to_csv(args.output_dir / "latent_coordinates.csv")
    loadings.to_csv(args.output_dir / "gene_loadings.csv")
    top_loadings.to_csv(args.output_dir / "top_gene_loadings.csv", index=False)
    work.obsm["X_ldvae"] = latent
    for column in latent_columns:
        work.obs[column] = latent_df[column]
    work.write_h5ad(args.output_dir / "ldvae_annotated.h5ad")
    model.save(args.output_dir / "model", overwrite=False)
    for key, history in model.history.items():
        history.to_csv(args.output_dir / f"history_{key}.csv")

    config = {
        **plan,
        "n_cells": work.n_obs,
        "n_genes_modeled": work.n_vars,
        "versions": {
            "anndata": ad.__version__,
            "scanpy": sc.__version__,
            "scvi_tools": scvi.__version__,
        },
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
