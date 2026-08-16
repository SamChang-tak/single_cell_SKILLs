#!/usr/bin/env python3
"""Run a reproducible CellRank pseudotime workflow on an AnnData file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input .h5ad file")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--time-key", required=True, help="Pseudotime column in adata.obs")
    parser.add_argument("--cluster-key", required=True, help="Cluster column in adata.obs")
    parser.add_argument("--n-states", type=int, default=10)
    parser.add_argument("--n-terminal-states", type=int, default=6)
    parser.add_argument("--pseudotime-weight", type=float, default=0.8)
    parser.add_argument("--lineage", help="Optional lineage for driver-gene calculation")
    parser.add_argument(
        "--preprocess",
        action="store_true",
        help="Filter, normalize, log-transform, compute PCA and neighbors",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Validate arguments without importing CellRank"
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.n_states < 2:
        raise ValueError("--n-states must be at least 2")
    if not 1 <= args.n_terminal_states <= args.n_states:
        raise ValueError("--n-terminal-states must be between 1 and --n-states")
    if not 0 <= args.pseudotime_weight <= 1:
        raise ValueError("--pseudotime-weight must be between 0 and 1")
    if not args.dry_run and not args.input.is_file():
        raise FileNotFoundError(args.input)


def main() -> None:
    args = parse_args()
    validate_args(args)
    config = {
        "input": str(args.input),
        "time_key": args.time_key,
        "cluster_key": args.cluster_key,
        "n_states": args.n_states,
        "n_terminal_states": args.n_terminal_states,
        "pseudotime_weight": args.pseudotime_weight,
        "connectivity_weight": round(1 - args.pseudotime_weight, 10),
        "lineage": args.lineage,
        "preprocess": args.preprocess,
    }
    if args.dry_run:
        print(json.dumps(config, indent=2, sort_keys=True))
        return

    import cellrank as cr
    import scanpy as sc

    adata = sc.read_h5ad(args.input)
    missing = [key for key in (args.time_key, args.cluster_key) if key not in adata.obs]
    if missing:
        raise KeyError(f"Missing adata.obs columns: {missing}")

    if args.preprocess:
        sc.pp.filter_genes(adata, min_cells=5)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=3000)
        sc.tl.pca(adata, random_state=0)
        sc.pp.neighbors(adata, random_state=0)
    elif "connectivities" not in adata.obsp:
        raise KeyError("adata.obsp['connectivities'] is required unless --preprocess is used")

    pk = cr.kernels.PseudotimeKernel(
        adata, time_key=args.time_key
    ).compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    kernel = args.pseudotime_weight * pk + (1 - args.pseudotime_weight) * ck

    estimator = cr.estimators.GPCCA(kernel)
    estimator.compute_schur()
    estimator.compute_macrostates(n_states=args.n_states, cluster_key=args.cluster_key)
    estimator.predict_terminal_states(method="top_n", n_states=args.n_terminal_states)
    estimator.compute_fate_probabilities()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.lineage:
        drivers = estimator.compute_lineage_drivers(lineages=args.lineage)
        drivers.to_csv(args.output_dir / "lineage_drivers.tsv", sep="\t")

    kernel.write_to_adata()
    adata.write_h5ad(args.output_dir / "cellrank_results.h5ad")
    config["cellrank_version"] = cr.__version__
    config["scanpy_version"] = sc.__version__
    (args.output_dir / "run_metadata.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
