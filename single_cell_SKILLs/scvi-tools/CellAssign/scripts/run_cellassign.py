#!/usr/bin/env python3
"""Run a guarded CellAssign annotation workflow from H5AD and marker CSV."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input AnnData H5AD")
    parser.add_argument("--markers", required=True, type=Path, help="Binary gene-by-cell-type CSV")
    parser.add_argument("--output-dir", required=True, type=Path, help="New output directory")
    parser.add_argument("--count-layer", help="Raw-count layer; default uses adata.X")
    parser.add_argument("--size-factor-key", default="size_factor")
    parser.add_argument("--max-epochs", type=int, help="Optional CellAssign training limit")
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
    import scvi
    from scipy import sparse
    from scvi.external import CellAssign

    result_names = {
        "cellassign_probabilities.csv", "cellassign_assignments.csv",
        "cellassign_annotated.h5ad", "run_config.json",
    }
    if args.output_dir.exists() and any((args.output_dir / name).exists() for name in result_names):
        raise FileExistsError(f"CellAssign outputs already exist: {args.output_dir}")
    if args.max_epochs is not None and args.max_epochs <= 0:
        raise ValueError("--max-epochs must be positive")

    adata = ad.read_h5ad(args.input, backed="r")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("AnnData observation and gene identifiers must be unique")
    if args.count_layer:
        if args.count_layer not in adata.layers:
            raise KeyError(f"Count layer not found: {args.count_layer}")
        counts = adata.layers[args.count_layer]
    else:
        counts = adata.X

    # Stream the full matrix so size factors use every gene without loading a
    # large H5AD into memory. This also validates every stored count value.
    library_size = np.zeros(adata.n_obs, dtype=np.float64)
    chunk_size = 4096
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        chunk = counts[start:stop]
        count_values = chunk.data if sparse.issparse(chunk) else np.asarray(chunk).ravel()
        if count_values.size and (not np.isfinite(count_values).all() or count_values.min() < 0):
            raise ValueError("Counts must be finite and nonnegative")
        library_size[start:stop] = np.asarray(chunk.sum(axis=1)).ravel()

    marker_matrix = pd.read_csv(args.markers, index_col=0)
    if marker_matrix.empty or not marker_matrix.index.is_unique or not marker_matrix.columns.is_unique:
        raise ValueError("Marker matrix must be nonempty with unique genes and cell types")
    marker_matrix = marker_matrix.apply(pd.to_numeric, errors="raise")
    values = marker_matrix.to_numpy()
    if not np.isin(values, [0, 1]).all():
        raise ValueError("Marker matrix values must be binary (0 or 1)")

    if not np.isfinite(library_size).all() or library_size.mean() <= 0:
        raise ValueError("Cannot compute size factors from empty or invalid libraries")

    missing = marker_matrix.index.difference(adata.var_names)
    present = marker_matrix.index.intersection(adata.var_names, sort=False)
    if len(present) == 0:
        raise ValueError("No marker genes match adata.var_names")
    aligned_markers = marker_matrix.loc[present]
    empty_types = aligned_markers.columns[aligned_markers.sum(axis=0) == 0].tolist()
    if empty_types:
        raise ValueError(f"Cell types without available markers: {empty_types}")

    # Subset before copying so large source objects are not duplicated in full.
    # If a layer supplies counts, replace the marker-subset X explicitly.
    bdata = adata[:, present].to_memory()
    if args.count_layer:
        marker_positions = adata.var_names.get_indexer(present)
        bdata.X = counts[:, marker_positions].copy()
    bdata.obs[args.size_factor_key] = library_size / library_size.mean()

    scvi.settings.seed = args.seed
    CellAssign.setup_anndata(bdata, size_factor_key=args.size_factor_key)
    model = CellAssign(bdata, aligned_markers)
    train_kwargs = {} if args.max_epochs is None else {"max_epochs": args.max_epochs}
    model.train(**train_kwargs)
    probabilities = model.predict()
    probabilities.index = bdata.obs_names

    ordered = np.sort(probabilities.to_numpy(), axis=1)
    top = probabilities.idxmax(axis=1)
    maximum = ordered[:, -1]
    second = ordered[:, -2] if probabilities.shape[1] > 1 else np.zeros(bdata.n_obs)
    clipped = np.clip(probabilities.to_numpy(), 1e-12, 1.0)
    entropy = -(clipped * np.log(clipped)).sum(axis=1)
    assignments = pd.DataFrame(
        {
            "cellassign_label": top,
            "max_probability": maximum,
            "second_probability": second,
            "probability_margin": maximum - second,
            "entropy": entropy,
        },
        index=bdata.obs_names,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    probabilities.to_csv(args.output_dir / "cellassign_probabilities.csv")
    assignments.to_csv(args.output_dir / "cellassign_assignments.csv")
    for column in assignments.columns:
        bdata.obs[column] = assignments[column]
    bdata.obsm["cellassign_probabilities"] = probabilities.to_numpy()
    bdata.uns["cellassign_probability_columns"] = probabilities.columns.to_list()
    bdata.write_h5ad(args.output_dir / "cellassign_annotated.h5ad")
    (args.output_dir / "missing_marker_genes.txt").write_text(
        "\n".join(map(str, missing)) + ("\n" if len(missing) else ""), encoding="utf-8"
    )

    config = {
        **plan,
        "n_cells": bdata.n_obs,
        "n_input_genes": adata.n_vars,
        "n_marker_genes_used": bdata.n_vars,
        "n_marker_genes_missing": len(missing),
        "candidate_cell_types": probabilities.columns.to_list(),
        "versions": {"anndata": ad.__version__, "scvi_tools": scvi.__version__},
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
