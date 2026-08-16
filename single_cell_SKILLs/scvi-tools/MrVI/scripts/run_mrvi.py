#!/usr/bin/env python3
"""Fit MRVI and export its sample-invariant representation and diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input H5AD")
    parser.add_argument("--output-dir", required=True, type=Path, help="New output directory")
    parser.add_argument("--sample-key", required=True, help="Biological sample/donor obs key")
    parser.add_argument("--batch-key", help="Optional nuisance batch obs key")
    parser.add_argument("--count-layer", help="Raw-count layer; default uses adata.X")
    parser.add_argument("--select-hvg", type=int, help="Select this many seurat_v3 HVGs")
    parser.add_argument("--distance-groupby", help="Optional obs key for grouped sample distances")
    parser.add_argument("--distance-batch-size", type=int, default=32)
    parser.add_argument("--max-epochs", type=int, default=400)
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
    from scvi.external import MRVI

    if args.output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {args.output_dir}")
    if args.max_epochs <= 0 or args.n_neighbors <= 0 or args.distance_batch_size <= 0:
        raise ValueError("Epoch, neighbor, and distance batch sizes must be positive")
    if args.select_hvg is not None and args.select_hvg <= 0:
        raise ValueError("--select-hvg must be positive")

    adata = ad.read_h5ad(args.input)
    required = [args.sample_key]
    if args.batch_key:
        required.append(args.batch_key)
    if args.distance_groupby:
        required.append(args.distance_groupby)
    missing = [key for key in required if key not in adata.obs]
    if missing:
        raise KeyError(f"Missing adata.obs columns: {missing}")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell and gene identifiers must be unique")
    if adata.obs[args.sample_key].isna().any() or adata.obs[args.sample_key].nunique() < 2:
        raise ValueError("Sample key must be complete and contain at least two samples")

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
    work.layers["_mrvi_counts"] = counts.copy()
    if args.select_hvg and args.select_hvg < work.n_vars:
        sc.pp.highly_variable_genes(
            work,
            n_top_genes=args.select_hvg,
            subset=True,
            flavor="seurat_v3",
            layer="_mrvi_counts",
        )

    scvi.settings.seed = args.seed
    MRVI.setup_anndata(
        work,
        layer="_mrvi_counts",
        sample_key=args.sample_key,
        batch_key=args.batch_key,
    )
    model = MRVI(work)
    model.train(max_epochs=args.max_epochs)

    u = model.get_latent_representation()
    columns = [f"u_{i}" for i in range(u.shape[1])]
    u_df = pd.DataFrame(u, index=work.obs_names, columns=columns)
    work.obsm["u"] = u
    sc.pp.neighbors(work, use_rep="u", n_neighbors=args.n_neighbors)
    sc.tl.umap(work, min_dist=0.3)

    args.output_dir.mkdir(parents=True)
    u_df.to_csv(args.output_dir / "mrvi_u_coordinates.csv")
    work.write_h5ad(args.output_dir / "mrvi_annotated.h5ad")
    model.save(args.output_dir / "model", overwrite=False)
    for key, frame in model.history.items():
        frame.to_csv(args.output_dir / f"history_{key}.csv")

    sample_fields = [args.sample_key] + ([args.batch_key] if args.batch_key else [])
    sample_design = work.obs[sample_fields].drop_duplicates().sort_values(args.sample_key)
    sample_design.to_csv(args.output_dir / "sample_design.csv", index=False)
    cells_per_sample = work.obs[args.sample_key].value_counts().rename("n_cells")
    cells_per_sample.to_csv(args.output_dir / "cells_per_sample.csv")

    if args.distance_groupby:
        distances = model.get_local_sample_distances(
            keep_cell=False,
            groupby=args.distance_groupby,
            batch_size=args.distance_batch_size,
        )
        distances.to_netcdf(args.output_dir / "local_sample_distances.nc")

    config = {
        **plan,
        "n_cells": work.n_obs,
        "n_genes_modeled": work.n_vars,
        "n_samples": int(work.obs[args.sample_key].nunique()),
        "versions": {
            "anndata": ad.__version__,
            "scanpy": sc.__version__,
            "scvi_tools": scvi.__version__,
        },
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
