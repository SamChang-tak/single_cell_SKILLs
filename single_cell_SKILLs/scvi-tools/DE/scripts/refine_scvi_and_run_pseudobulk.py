#!/usr/bin/env python3
"""Correct finite-sampling scVI plots and run donor-aware edgeR pseudobulk DE."""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")


def save_rgb(fig, path: Path, dpi: int) -> None:
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)


def corrected_scvi_plot(table: pd.DataFrame, cell_type: str, path: Path, cap: float = 4.0) -> dict:
    fdr_col = next(c for c in table.columns if c.startswith("is_de_fdr_"))
    sig = table[fdr_col].astype(bool)
    up = sig & table.lfc_mean.gt(0)
    down = sig & table.lfc_mean.lt(0)
    saturated = table.proba_not_de.le(0)
    score = -np.log10(np.maximum(table.proba_not_de.astype(float), 10 ** -cap))
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    ax.scatter(table.loc[~sig, "lfc_mean"], score[~sig], s=7, c="#bdbdbd", alpha=.5, linewidths=0)
    ax.scatter(table.loc[up, "lfc_mean"], score[up], s=9, c="#d62728", alpha=.7, linewidths=0,
               label=f"Bayesian FDR up ({up.sum():,})")
    ax.scatter(table.loc[down, "lfc_mean"], score[down], s=9, c="#1f77b4", alpha=.7, linewidths=0,
               label=f"Bayesian FDR down ({down.sum():,})")
    if saturated.any():
        ax.scatter(table.loc[saturated, "lfc_mean"], np.full(saturated.sum(), cap), marker="v", s=13,
                   facecolors="none", edgecolors="#111111", linewidths=.45,
                   label=f"Saturated at cap ({saturated.sum():,})")
    ax.axvline(.25, ls="--", c="black", lw=1)
    ax.axvline(-.25, ls="--", c="black", lw=1)
    ax.axhline(cap, ls=":", c="#555555", lw=.8)
    ax.set_ylim(-.08, cap + .25)
    ax.set_xlabel("Posterior mean LFC: Fibrotic liver minus Healthy")
    ax.set_ylabel(f"Posterior evidence −log10(P non-DE), capped at {cap:g}")
    ax.set_title(f"{cell_type} — scVI cell-level sensitivity analysis")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    save_rgb(fig, path, 170)
    return {"cell_type": cell_type, "up": int(up.sum()), "down": int(down.sum()),
            "saturated": int(saturated.sum())}


def pseudobulk_plot(table: pd.DataFrame, cell_type: str, path: Path) -> dict:
    sig = table.FDR.lt(.05) & table.logFC.abs().ge(1)
    up = sig & table.logFC.gt(0)
    down = sig & table.logFC.lt(0)
    score = -np.log10(np.maximum(table.FDR.astype(float), np.nextafter(0, 1)))
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    ax.scatter(table.loc[~sig, "logFC"], score[~sig], s=7, c="#bdbdbd", alpha=.5, linewidths=0)
    ax.scatter(table.loc[up, "logFC"], score[up], s=10, c="#d62728", alpha=.75, linewidths=0,
               label=f"FDR<0.05, up ({up.sum():,})")
    ax.scatter(table.loc[down, "logFC"], score[down], s=10, c="#1f77b4", alpha=.75, linewidths=0,
               label=f"FDR<0.05, down ({down.sum():,})")
    ax.axvline(1, ls="--", c="black", lw=1)
    ax.axvline(-1, ls="--", c="black", lw=1)
    ax.axhline(-np.log10(.05), ls=":", c="#555555", lw=.8, label="FDR = 0.05")
    ax.set_xlabel("edgeR log2FC: Fibrotic liver minus Healthy")
    ax.set_ylabel("−log10(edgeR FDR)")
    ax.set_title(f"{cell_type} — donor-level pseudobulk")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    save_rgb(fig, path, 170)
    return {"cell_type": cell_type, "genes_tested": len(table), "up": int(up.sum()),
            "down": int(down.sum()), "total": int(sig.sum())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--r-script", required=True, type=Path)
    args = parser.parse_args()
    out = args.results_dir
    cfg = json.loads((out / "run_config.json").read_text())
    audit = pd.read_csv(out / "celltype_audit.csv")
    run_types = audit.loc[audit.status.eq("run"), "cell_type"].tolist()

    # Preserve the original misleading plots and replace them with capped versions.
    scvi_rows = []
    scvi_tables = []
    for cell_type in run_types:
        folder = out / "celltypes" / slug(cell_type)
        original = folder / "volcano_uncapped_original.png"
        if not original.exists():
            shutil.copy2(folder / "volcano.png", original)
        table = pd.read_csv(folder / "de_results.csv")
        scvi_rows.append(corrected_scvi_plot(table, cell_type, folder / "volcano.png"))
        scvi_tables.append((cell_type, table))
    pd.DataFrame(scvi_rows).to_csv(out / "scvi_plot_saturation_audit.csv", index=False)
    original_combined = out / "combined_volcano_uncapped_original.png"
    if not original_combined.exists():
        shutil.copy2(out / "combined_volcano.png", original_combined)
    cols = 3
    rows = math.ceil(len(scvi_tables) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(15, 4.5 * rows), squeeze=False)
    for ax, (cell_type, table) in zip(axes.flat, scvi_tables):
        fdr_col = next(c for c in table.columns if c.startswith("is_de_fdr_"))
        sig = table[fdr_col].astype(bool)
        up = sig & table.lfc_mean.gt(0)
        down = sig & table.lfc_mean.lt(0)
        saturated = table.proba_not_de.le(0)
        score = -np.log10(np.maximum(table.proba_not_de.astype(float), 1e-4))
        ax.scatter(table.loc[~sig, "lfc_mean"], score[~sig], s=3, c="#c8c8c8", alpha=.45, linewidths=0)
        ax.scatter(table.loc[up, "lfc_mean"], score[up], s=4, c="#d62728", alpha=.65, linewidths=0)
        ax.scatter(table.loc[down, "lfc_mean"], score[down], s=4, c="#1f77b4", alpha=.65, linewidths=0)
        if saturated.any():
            ax.scatter(table.loc[saturated, "lfc_mean"], np.full(saturated.sum(), 4), marker="v", s=6,
                       facecolors="none", edgecolors="#111111", linewidths=.3)
        ax.axvline(.25, ls="--", c="black", lw=.7)
        ax.axvline(-.25, ls="--", c="black", lw=.7)
        ax.set_ylim(-.08, 4.2)
        ax.set_title(f"{cell_type}\nup {up.sum():,} / down {down.sum():,} / capped {saturated.sum():,}")
        ax.set_xlabel("Posterior mean LFC")
        ax.set_ylabel("Posterior score (cap 4)")
    for ax in axes.flat[len(scvi_tables):]:
        ax.axis("off")
    fig.suptitle("Secondary scVI cell-level sensitivity analysis", fontsize=16)
    fig.tight_layout()
    save_rgb(fig, out / "combined_volcano.png", 180)

    # Read raw counts and retain only the intended liver disease/control comparison.
    data = ad.read_h5ad(cfg["input"])
    obs = data.obs
    keep = obs[cfg["tissue_key"]].astype(str).eq(cfg["tissue"]) & obs[cfg["condition_key"]].astype(str).isin([cfg["disease"], cfg["control"]])
    data = data[keep].copy()
    data.X = sparse.csr_matrix(data.X)
    pbulk = out / "pseudobulk"
    (pbulk / "celltypes").mkdir(parents=True, exist_ok=True)
    summaries = []
    plot_tables = []

    for cell_type in run_types:
        cell_mask = data.obs[cfg["celltype_key"]].astype(str).eq(cell_type).to_numpy()
        subset = data[cell_mask]
        meta = subset.obs[[cfg["donor_key"], cfg["condition_key"]]].copy()
        meta.columns = ["donor", "condition_label"]
        meta["sample_id"] = meta.donor.astype(str).map(slug)
        sample_rows = meta.drop_duplicates("sample_id").sort_values("sample_id")
        sample_rows["condition"] = np.where(sample_rows.condition_label.eq(cfg["disease"]), "disease", "control")
        sample_ids = sample_rows.sample_id.tolist()
        vectors = []
        cells = []
        for sample_id in sample_ids:
            mask = meta.sample_id.eq(sample_id).to_numpy()
            vectors.append(np.asarray(subset.X[mask].sum(axis=0)).ravel())
            cells.append(int(mask.sum()))
        counts = np.column_stack(vectors).astype(np.int64, copy=False)
        nonzero = counts.sum(axis=1) > 0
        counts_df = pd.DataFrame(counts[nonzero], index=subset.var_names[nonzero], columns=sample_ids)
        folder = pbulk / "celltypes" / slug(cell_type)
        folder.mkdir(parents=True, exist_ok=True)
        counts_file = folder / "pseudobulk_counts.csv.gz"
        samples_file = folder / "sample_metadata.csv"
        results_file = folder / "edger_results.csv"
        counts_df.to_csv(counts_file, compression="gzip", index_label="gene")
        sample_rows = sample_rows[["sample_id", "donor", "condition"]].copy()
        sample_rows["cells"] = cells
        sample_rows.to_csv(samples_file, index=False)
        subprocess.run(["Rscript", str(args.r_script), str(counts_file), str(samples_file), str(results_file)], check=True)
        result = pd.read_csv(results_file)
        call = result.FDR.lt(.05) & result.logFC.abs().ge(1)
        result["significant_fdr05_abs_log2fc1"] = call
        result.to_csv(results_file, index=False)
        result.loc[call].to_csv(folder / "edger_significant.csv", index=False)
        stats = pseudobulk_plot(result, cell_type, folder / "volcano.png")
        stats.update(disease_donors=int((sample_rows.condition == "disease").sum()),
                     control_donors=int((sample_rows.condition == "control").sum()),
                     disease_cells=int(sample_rows.loc[sample_rows.condition == "disease", "cells"].sum()),
                     control_cells=int(sample_rows.loc[sample_rows.condition == "control", "cells"].sum()),
                     low_replication=bool((sample_rows.groupby("condition").size() < 3).any()))
        summaries.append(stats)
        plot_tables.append((cell_type, result))

    summary = pd.DataFrame(summaries).sort_values("total", ascending=False)
    summary.to_csv(pbulk / "pseudobulk_summary_by_celltype.csv", index=False)
    pd.concat([x.assign(cell_type=ct) for ct, x in plot_tables], ignore_index=True).to_csv(
        pbulk / "all_celltypes_edger_results.csv.gz", index=False, compression="gzip")

    cols = 3
    rows = math.ceil(len(plot_tables) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(15, 4.5 * rows), squeeze=False)
    for ax, (cell_type, table) in zip(axes.flat, plot_tables):
        sig = table.FDR.lt(.05) & table.logFC.abs().ge(1)
        up = sig & table.logFC.gt(0)
        down = sig & table.logFC.lt(0)
        score = -np.log10(np.maximum(table.FDR.astype(float), np.nextafter(0, 1)))
        ax.scatter(table.loc[~sig, "logFC"], score[~sig], s=3, c="#c8c8c8", alpha=.45, linewidths=0)
        ax.scatter(table.loc[up, "logFC"], score[up], s=4, c="#d62728", alpha=.65, linewidths=0)
        ax.scatter(table.loc[down, "logFC"], score[down], s=4, c="#1f77b4", alpha=.65, linewidths=0)
        ax.axvline(1, ls="--", c="black", lw=.7)
        ax.axvline(-1, ls="--", c="black", lw=.7)
        ax.axhline(-np.log10(.05), ls=":", c="#555555", lw=.7)
        ax.set_title(f"{cell_type}\nup {up.sum():,} / down {down.sum():,}")
        ax.set_xlabel("edgeR log2FC")
        ax.set_ylabel("−log10 edgeR FDR")
    for ax in axes.flat[len(plot_tables):]:
        ax.axis("off")
    fig.suptitle("Primary donor-level pseudobulk: Fibrotic liver vs Healthy", fontsize=16)
    fig.tight_layout()
    save_rgb(fig, pbulk / "combined_pseudobulk_volcano.png", 180)


if __name__ == "__main__":
    main()
