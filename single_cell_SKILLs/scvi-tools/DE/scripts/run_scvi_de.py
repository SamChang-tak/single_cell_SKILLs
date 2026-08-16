#!/usr/bin/env python3
"""Train/load SCVI and run a guarded posterior differential-expression contrast."""
from __future__ import annotations

import argparse
import importlib.metadata
import inspect
import json
import platform
from html import escape
from pathlib import Path


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--groupby", required=True)
    p.add_argument("--group1", required=True)
    p.add_argument("--group2", help="Omit for group1 versus all remaining cells")
    p.add_argument("--batch-key")
    p.add_argument("--count-layer")
    p.add_argument("--gene-selection", choices=("auto", "poisson", "existing", "all"), default="auto")
    p.add_argument("--n-top-genes", type=int, default=4000)
    p.add_argument("--weights", choices=("importance", "uniform"), default="importance")
    p.add_argument("--batch-correction", choices=("auto", "true", "false"), default="auto")
    p.add_argument("--delta", type=float, default=0.25)
    p.add_argument("--fdr-target", type=float, default=0.05)
    p.add_argument("--pseudocounts", type=float, default=0.0)
    p.add_argument("--max-epochs", type=int, default=400)
    p.add_argument("--early-stopping-patience", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--model-dir", type=Path, help="Load a compatible saved SCVI model instead of training")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main():
    a = arguments()
    plan = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    if a.dry_run:
        print(json.dumps(plan, indent=2)); return
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scvi
    from scipy import sparse

    if a.output_dir.exists() and any(a.output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {a.output_dir}")
    a.output_dir.mkdir(parents=True, exist_ok=True)
    if not 0 < a.fdr_target < 1: raise ValueError("--fdr-target must be between 0 and 1")
    if a.delta < 0 or a.pseudocounts < 0: raise ValueError("delta and pseudocounts must be nonnegative")
    data = ad.read_h5ad(a.input)
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Cell and gene identifiers must be unique")
    if a.groupby not in data.obs: raise KeyError(f"Missing group key: {a.groupby}")
    if a.batch_key and a.batch_key not in data.obs: raise KeyError(f"Missing batch key: {a.batch_key}")
    if a.count_layer:
        if a.count_layer not in data.layers: raise KeyError(f"Missing count layer: {a.count_layer}")
        counts = data.layers[a.count_layer]
    else: counts = data.X
    vals = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
    probe = vals[:min(len(vals), 500_000)]
    if len(probe) and (np.any(~np.isfinite(probe)) or np.any(probe < 0) or not np.allclose(probe, np.round(probe))):
        raise ValueError("Selected matrix does not resemble nonnegative integer counts")
    groups = data.obs[a.groupby].astype(str)
    idx1 = (groups == a.group1).to_numpy()
    idx2 = (groups != a.group1).to_numpy() if a.group2 is None else (groups == a.group2).to_numpy()
    if not idx1.any() or not idx2.any(): raise ValueError("Both comparison groups must contain cells")
    if np.any(idx1 & idx2): raise ValueError("Comparison groups overlap")
    label2 = "Rest" if a.group2 is None else a.group2

    batches1 = batches2 = set()
    if a.batch_key:
        batch = data.obs[a.batch_key].astype(str)
        batches1, batches2 = set(batch[idx1]), set(batch[idx2])
    overlap = sorted(batches1 & batches2)
    requested = a.batch_correction
    if requested == "auto": batch_correction = bool(a.batch_key and overlap)
    else: batch_correction = requested == "true"
    if batch_correction and not a.batch_key: raise ValueError("Batch correction requested without --batch-key")
    if batch_correction and not overlap: raise ValueError("Batch correction requested but groups share no batch levels")
    audit = pd.DataFrame([
        {"group": a.group1, "cells": int(idx1.sum()), "batches": "; ".join(sorted(batches1)), "shared_batches": "; ".join(overlap)},
        {"group": label2, "cells": int(idx2.sum()), "batches": "; ".join(sorted(batches2)), "shared_batches": "; ".join(overlap)},
    ])
    audit.to_csv(a.output_dir / "group_batch_audit.csv", index=False)

    data.layers["_scvi_de_counts"] = counts.copy()
    selection = a.gene_selection
    if selection == "auto": selection = "existing" if "highly_variable" in data.var and data.var["highly_variable"].sum() else "poisson"
    if selection == "poisson":
        scvi.data.poisson_gene_selection(data, n_top_genes=min(a.n_top_genes, data.n_vars), batch_key=a.batch_key, inplace=True)
        data = data[:, data.var["highly_variable"]].copy()
    elif selection == "existing":
        if "highly_variable" not in data.var or not data.var["highly_variable"].sum(): raise ValueError("No existing highly_variable mask")
        data = data[:, data.var["highly_variable"].astype(bool)].copy()
    else: data = data.copy()
    if not sparse.isspmatrix_csr(data.layers["_scvi_de_counts"]):
        data.layers["_scvi_de_counts"] = sparse.csr_matrix(data.layers["_scvi_de_counts"])
    data.var_names.to_series().rename("gene").to_csv(a.output_dir / "trained_genes.csv", index=False)

    scvi.settings.seed = a.seed
    scvi.model.SCVI.setup_anndata(data, layer="_scvi_de_counts", batch_key=a.batch_key)
    if a.model_dir:
        model = scvi.model.SCVI.load(a.model_dir, adata=data)
        trained = False
    else:
        model = scvi.model.SCVI(data, gene_likelihood="nb")
        model.train(max_epochs=a.max_epochs, check_val_every_n_epoch=1, early_stopping=True,
                    early_stopping_patience=a.early_stopping_patience,
                    early_stopping_monitor="elbo_validation")
        model.save(a.output_dir / "model", overwrite=False, save_anndata=False)
        trained = True
    history_parts = {}
    for key, value in model.history.items():
        if isinstance(value, pd.DataFrame):
            history_parts[key] = value.iloc[:, 0].reset_index(drop=True)
        elif isinstance(value, pd.Series):
            history_parts[key] = value.reset_index(drop=True)
    history = pd.DataFrame(history_parts)
    history.to_csv(a.output_dir / "training_history.csv", index_label="epoch")
    fig, ax = plt.subplots(figsize=(8,5))
    for key in ("elbo_train", "elbo_validation"):
        if key in history: ax.plot(history.index, history[key], label=key)
    ax.set_yscale("log"); ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO loss (log scale)"); ax.set_title("SCVI convergence"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(a.output_dir / "convergence.png", dpi=180); plt.close(fig)

    data.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(data, use_rep="X_scVI"); sc.tl.umap(data, random_state=a.seed)
    fig, axes = plt.subplots(1, 2 if a.batch_key else 1, figsize=(12 if a.batch_key else 7, 6))
    axes = np.atleast_1d(axes)
    sc.pl.umap(data, color=a.groupby, ax=axes[0], show=False, title=f"Latent UMAP: {a.groupby}")
    if a.batch_key: sc.pl.umap(data, color=a.batch_key, ax=axes[1], show=False, title=f"Latent UMAP: {a.batch_key}")
    fig.tight_layout(); fig.savefig(a.output_dir / "latent_umap.png", dpi=180); plt.close(fig)

    de = model.differential_expression(
        idx1=idx1, idx2=idx2, mode="change", weights=a.weights,
        filter_outlier_cells=(a.weights == "importance"), batch_correction=batch_correction,
        delta=a.delta, fdr_target=a.fdr_target, pseudocounts=a.pseudocounts,
    )
    de.index.name = "gene"; de.to_csv(a.output_dir / "de_results.csv")
    fdr_cols = [c for c in de if c.startswith("is_de_fdr_")]
    if fdr_cols: de.loc[de[fdr_cols[0]].astype(bool)].to_csv(a.output_dir / "de_genes_fdr.csv")
    score = -np.log10(np.clip(pd.to_numeric(de["proba_not_de"], errors="coerce"), 1e-300, 1))
    is_de = de[fdr_cols[0]].astype(bool) if fdr_cols else de["proba_de"] >= .95
    fig, ax = plt.subplots(figsize=(9,7))
    ax.scatter(de.loc[~is_de,"lfc_mean"], score[~is_de], s=7, c="#bdbdbd", alpha=.5, linewidths=0)
    ax.scatter(de.loc[is_de,"lfc_mean"], score[is_de], s=9, c="#d62728", alpha=.7, linewidths=0)
    ax.axvline(a.delta, ls="--", c="black", lw=1); ax.axvline(-a.delta, ls="--", c="black", lw=1)
    ax.set_xlabel(f"Posterior mean LFC: {a.group1} minus {label2}")
    ax.set_ylabel("-log10 posterior non-DE probability"); ax.set_title("scVI change-mode differential expression")
    fig.tight_layout(); fig.savefig(a.output_dir / "volcano.png", dpi=180); plt.close(fig)
    top = de.assign(abs_lfc=de["lfc_mean"].abs()).sort_values([fdr_cols[0],"proba_de","abs_lfc"] if fdr_cols else ["proba_de","abs_lfc"], ascending=False).head(30)
    heat = top[["scale1","scale2"]].copy(); heat.columns=[a.group1,label2]
    fig, ax = plt.subplots(figsize=(7, max(6, len(heat)*.25)))
    im=ax.imshow(np.log1p(heat), cmap="viridis", aspect="auto")
    ax.set_xticks([0,1],heat.columns); ax.set_yticks(range(len(heat)),heat.index); ax.set_title("Top DE posterior expression scales (log1p)")
    fig.colorbar(im,ax=ax,label="log1p posterior scale"); fig.tight_layout(); fig.savefig(a.output_dir / "top_de_heatmap.png",dpi=180); plt.close(fig)
    data.write_h5ad(a.output_dir / "processed.h5ad", compression="gzip")

    sig = int(is_de.sum())
    config = {**plan, "group2_effective": label2, "group1_cells": int(idx1.sum()), "group2_cells": int(idx2.sum()),
              "shared_batches": overlap, "batch_correction_effective": batch_correction, "gene_selection_effective": selection,
              "trained_genes": data.n_vars, "trained_new_model": trained, "de_genes": sig,
              "versions": {"python": platform.python_version(), "scvi_tools": importlib.metadata.version("scvi-tools"),
                           "anndata": importlib.metadata.version("anndata"), "scanpy": importlib.metadata.version("scanpy")},
              "de_signature": str(inspect.signature(model.differential_expression))}
    (a.output_dir / "run_config.json").write_text(json.dumps(config, indent=2, default=str)+"\n")
    html=f"""<!doctype html><html><head><meta charset='utf-8'><title>scVI DE</title><style>body{{font:15px/1.5 sans-serif;max-width:1100px;margin:30px auto;padding:0 20px}}img{{max-width:100%}}table{{border-collapse:collapse}}th,td{{border:1px solid #ddd;padding:6px}}</style></head><body><h1>scVI differential expression</h1><p><b>Contrast:</b> {escape(a.group1)} minus {escape(label2)}. Positive LFC is higher in {escape(a.group1)}.</p><p>{idx1.sum():,} vs {idx2.sum():,} cells; {data.n_vars:,} trained/tested genes; {sig:,} DE calls. Batch correction: {batch_correction}.</p><h2>Method</h2><p>Negative-binomial SCVI; change-mode Bayesian DE; weights={a.weights}; delta={a.delta}; FDR target={a.fdr_target}. Posterior probabilities are not frequentist p-values. Cells are not biological replicates; condition-level claims require donor-aware validation.</p><img src='convergence.png'><img src='latent_umap.png'><img src='volcano.png'><img src='top_de_heatmap.png'><h2>Batch audit</h2>{audit.to_html(index=False)}<h2>Limitations</h2><p>Results apply only to genes used to train the model. Batch correction is defensible only for overlapping batches. Review donor balance, raw expression prevalence, and sparse or technical top hits.</p></body></html>"""
    (a.output_dir / "index.html").write_text(html)


if __name__ == "__main__": main()
