#!/usr/bin/env python3
"""Run limma-voom on existing donor pseudobulks and create FDR/nominal outputs."""
from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")


def plot_panel(ax, table: pd.DataFrame, cell_type: str, mode: str, annotate: bool) -> tuple[int, int]:
    if mode == "fdr":
        probability = table["adj.P.Val"]
        call = probability.lt(.05) & table.logFC.abs().ge(1)
        ylabel = "−log10(limma adjusted P)"
        label = "FDR"
    else:
        probability = table["P.Value"]
        call = probability.lt(.05) & table.logFC.abs().ge(1)
        ylabel = "−log10(limma raw P)"
        label = "Nominal P"
    up = call & table.logFC.gt(0)
    down = call & table.logFC.lt(0)
    score = -np.log10(np.maximum(probability.astype(float), np.nextafter(0, 1)))
    ax.scatter(table.loc[~call, "logFC"], score[~call], s=7 if annotate else 3,
               c="#c8c8c8", alpha=.45, linewidths=0)
    ax.scatter(table.loc[up, "logFC"], score[up], s=11 if annotate else 5,
               c="#d62728", alpha=.75, linewidths=0, label=f"{label} up ({up.sum():,})")
    ax.scatter(table.loc[down, "logFC"], score[down], s=11 if annotate else 5,
               c="#1f77b4", alpha=.75, linewidths=0, label=f"{label} down ({down.sum():,})")
    ax.axvline(1, ls="--", c="black", lw=1 if annotate else .7)
    ax.axvline(-1, ls="--", c="black", lw=1 if annotate else .7)
    ax.axhline(-np.log10(.05), ls=":", c="#555555", lw=.9 if annotate else .7,
               label=f"{label} = 0.05" if annotate else None)
    ax.set_xlabel("limma log2FC: Fibrotic liver minus Healthy" if annotate else "limma log2FC")
    ax.set_ylabel(ylabel)
    title_mode = "FDR-controlled" if mode == "fdr" else "exploratory nominal P"
    ax.set_title(f"{cell_type} — limma-voom {title_mode}" if annotate else
                 f"{cell_type}\nup {up.sum():,} / down {down.sum():,}")
    if annotate:
        labels = table.loc[call].nsmallest(8, "adj.P.Val" if mode == "fdr" else "P.Value")
        value_col = "adj.P.Val" if mode == "fdr" else "P.Value"
        for _, row in labels.iterrows():
            y = -np.log10(max(float(row[value_col]), np.nextafter(0, 1)))
            ax.annotate(row["gene"], (row["logFC"], y), xytext=(3, 3),
                        textcoords="offset points", fontsize=7)
        ax.legend(frameon=False, fontsize=9)
    return int(up.sum()), int(down.sum())


def combined_plot(plot_tables: list[tuple[str, pd.DataFrame]], mode: str, path: Path) -> None:
    cols = 3
    rows = math.ceil(len(plot_tables) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(15, 4.5 * rows), squeeze=False)
    for ax, (cell_type, table) in zip(axes.flat, plot_tables):
        plot_panel(ax, table, cell_type, mode, annotate=False)
    for ax in axes.flat[len(plot_tables):]:
        ax.axis("off")
    descriptor = "FDR<0.05" if mode == "fdr" else "exploratory raw P<0.05"
    fig.suptitle(f"limma-voom donor pseudobulk ({descriptor}): Fibrotic liver vs Healthy", fontsize=16)
    fig.tight_layout()
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--r-script", required=True, type=Path)
    args = parser.parse_args()
    root = args.results_dir
    source = root / "pseudobulk"
    out = root / "pseudobulk_limma"
    (out / "celltypes").mkdir(parents=True, exist_ok=True)
    source_summary = pd.read_csv(source / "pseudobulk_summary_by_celltype.csv")
    summaries = []
    combined = []
    combined_fdr = []
    combined_nominal = []
    plot_tables = []

    for meta in source_summary.itertuples(index=False):
        cell_type = meta.cell_type
        source_folder = source / "celltypes" / slug(cell_type)
        folder = out / "celltypes" / slug(cell_type)
        folder.mkdir(parents=True, exist_ok=True)
        samples_file = folder / "sample_metadata.csv"
        shutil.copy2(source_folder / "sample_metadata.csv", samples_file)
        results_file = folder / "limma_results.csv"
        subprocess.run(["Rscript", str(args.r_script), str(source_folder / "pseudobulk_counts.csv.gz"),
                        str(samples_file), str(results_file)], check=True)
        table = pd.read_csv(results_file)
        fdr = table["adj.P.Val"].lt(.05) & table.logFC.abs().ge(1)
        nominal = table["P.Value"].lt(.05) & table.logFC.abs().ge(1)
        table["significant_fdr05_abs_log2fc1"] = fdr
        table["nominal_p05_abs_log2fc1"] = nominal
        table.to_csv(results_file, index=False)
        table.loc[fdr].sort_values("adj.P.Val").to_csv(folder / "limma_fdr_significant.csv", index=False)
        table.loc[nominal].sort_values("P.Value").to_csv(folder / "limma_nominal_significant.csv", index=False)
        fig, ax = plt.subplots(figsize=(8.5, 6.2))
        fdr_up, fdr_down = plot_panel(ax, table, cell_type, "fdr", annotate=True)
        fig.tight_layout(); fig.savefig(folder / "volcano_fdr.png", dpi=170, facecolor="white"); plt.close(fig)
        fig, ax = plt.subplots(figsize=(8.5, 6.2))
        nominal_up, nominal_down = plot_panel(ax, table, cell_type, "nominal", annotate=True)
        fig.tight_layout(); fig.savefig(folder / "volcano_nominal.png", dpi=170, facecolor="white"); plt.close(fig)
        summaries.append({"cell_type": cell_type, "genes_tested": len(table),
                          "fdr_up": fdr_up, "fdr_down": fdr_down, "fdr_total": int(fdr.sum()),
                          "nominal_up": nominal_up, "nominal_down": nominal_down,
                          "nominal_total": int(nominal.sum()), "disease_donors": meta.disease_donors,
                          "control_donors": meta.control_donors, "low_replication": meta.low_replication})
        combined.append(table.assign(cell_type=cell_type))
        combined_fdr.append(table.loc[fdr].assign(cell_type=cell_type))
        combined_nominal.append(table.loc[nominal].assign(cell_type=cell_type))
        plot_tables.append((cell_type, table))

    summary = pd.DataFrame(summaries).sort_values("fdr_total", ascending=False)
    summary.to_csv(out / "limma_summary_by_celltype.csv", index=False)
    pd.concat(combined, ignore_index=True).to_csv(out / "all_celltypes_limma_results.csv.gz", index=False, compression="gzip")
    pd.concat(combined_fdr, ignore_index=True).to_csv(out / "all_celltypes_limma_fdr_significant.csv.gz", index=False, compression="gzip")
    pd.concat(combined_nominal, ignore_index=True).to_csv(out / "all_celltypes_limma_nominal_significant.csv.gz", index=False, compression="gzip")
    combined_plot(plot_tables, "fdr", out / "combined_limma_fdr_volcano.png")
    combined_plot(plot_tables, "nominal", out / "combined_limma_nominal_volcano.png")


if __name__ == "__main__":
    main()
