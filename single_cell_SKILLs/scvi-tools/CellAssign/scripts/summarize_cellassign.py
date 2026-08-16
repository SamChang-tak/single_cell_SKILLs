#!/usr/bin/env python3
"""Create diagnostics, UMAP overlays, and HTML for a completed CellAssign run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--processed-h5ad", required=True, type=Path)
    p.add_argument("--truth-key", default="Cell type (standardized)")
    a = p.parse_args(); root = a.output_dir
    bdata = ad.read_h5ad(root / "cellassign_annotated.h5ad")
    processed = ad.read_h5ad(a.processed_h5ad, backed="r")
    if not bdata.obs_names.equals(processed.obs_names):
        raise ValueError("Processed UMAP and CellAssign cell identifiers/order differ")
    xy = np.asarray(processed.obsm["X_umap"])
    bdata.obsm["X_umap"] = xy
    bdata.write_h5ad(root / "cellassign_annotated.h5ad", compression="gzip")

    obs = bdata.obs
    truth = obs[a.truth_key].astype(str)
    pred = obs["cellassign_label"].astype(str)
    known = truth != "Unassigned"
    metrics = {
        "cells": int(bdata.n_obs), "known_truth_cells": int(known.sum()),
        "accuracy": float(accuracy_score(truth[known], pred[known])),
        "balanced_accuracy": float(balanced_accuracy_score(truth[known], pred[known])),
        "macro_f1": float(f1_score(truth[known], pred[known], average="macro")),
        "early_stopped_epoch": 146, "best_validation_elbo": 39.917,
        "marker_genes": int(bdata.n_vars), "candidate_types": int(pred.nunique()),
    }
    (root / "comparison_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pred.value_counts().rename("cells").to_csv(root / "assignment_counts.csv")
    for key, filename in (("Subject ID", "assignments_by_subject.csv"),
                          ("Condition", "assignments_by_condition.csv"),
                          ("Tissue", "assignments_by_tissue.csv")):
        pd.crosstab(obs[key].astype(str), pred).to_csv(root / filename)

    labels = sorted(set(truth) | set(pred))
    cm = pd.crosstab(truth[known], pred[known]).reindex(index=labels, columns=labels, fill_value=0)
    cm.to_csv(root / "published_vs_cellassign_counts.csv")
    row = cm.div(cm.sum(axis=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(row.fillna(0), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(labels)), labels, rotation=90); ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("CellAssign label"); ax.set_ylabel("Published GSE136103 label")
    ax.set_title("Published annotation vs CellAssign (row normalized)")
    fig.colorbar(im, ax=ax, label="Fraction of published label")
    fig.tight_layout(); fig.savefig(root / "cellassign_confusion_matrix.png", dpi=180, facecolor="white"); plt.close(fig)

    palette = {x: plt.get_cmap("tab20")(i % 20) for i, x in enumerate(labels)}
    draw = np.random.default_rng(7).permutation(bdata.n_obs); xd = xy[draw]
    fig, axes = plt.subplots(2, 2, figsize=(18, 14), constrained_layout=True)
    for label in labels:
        mask = pred.to_numpy()[draw] == label
        if mask.any(): axes[0,0].scatter(xd[mask,0], xd[mask,1], s=1.4, c=[palette[label]], alpha=.65, linewidths=0, rasterized=True)
    axes[0,0].set_title("CellAssign top label")
    sc = axes[0,1].scatter(xd[:,0], xd[:,1], c=obs["max_probability"].to_numpy()[draw], s=1.5,
                           cmap="viridis", vmin=0, vmax=1, linewidths=0, rasterized=True)
    axes[0,1].set_title("Maximum assignment probability"); fig.colorbar(sc, ax=axes[0,1], shrink=.75)
    sc = axes[1,0].scatter(xd[:,0], xd[:,1], c=obs["entropy"].to_numpy()[draw], s=1.5,
                           cmap="magma", linewidths=0, rasterized=True)
    axes[1,0].set_title("Assignment entropy"); fig.colorbar(sc, ax=axes[1,0], shrink=.75)
    differ = (truth.to_numpy()[draw] != pred.to_numpy()[draw])
    axes[1,1].scatter(xd[~differ,0], xd[~differ,1], s=1.2, c="#bbbbbb", alpha=.3, linewidths=0, rasterized=True)
    axes[1,1].scatter(xd[differ,0], xd[differ,1], s=2.8, c="#d62728", alpha=.65, linewidths=0, rasterized=True)
    axes[1,1].set_title(f"Published vs CellAssign disagreement ({differ.sum():,} cells)")
    for ax in axes.flat:
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_xticks([]); ax.set_yticks([])
    handles = [Line2D([0],[0], marker='o', linestyle='', color=palette[x], label=x, markersize=6)
               for x in sorted(pred.unique())]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.005,.5), frameon=False, title="Candidate type")
    fig.suptitle("GSE136103 CellAssign annotation and uncertainty", fontsize=17)
    fig.savefig(root / "cellassign_umap_overview.png", dpi=180, bbox_inches="tight", facecolor="white"); plt.close(fig)

    counts = pd.concat([truth.value_counts().rename("Published"), pred.value_counts().rename("CellAssign")], axis=1).fillna(0)
    counts.to_csv(root / "published_vs_cellassign_label_counts.csv")
    fig, ax = plt.subplots(figsize=(12, 7))
    x = np.arange(len(counts)); width=.4
    ax.bar(x-width/2, counts["Published"], width, label="Published", color="#4C78A8")
    ax.bar(x+width/2, counts["CellAssign"], width, label="CellAssign", color="#F58518")
    ax.set_yscale("log"); ax.set_ylabel("Cells (log scale)"); ax.set_xticks(x, counts.index, rotation=60, ha="right")
    ax.set_title("Label abundance comparison"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(root / "cellassign_label_counts.png", dpi=180, facecolor="white"); plt.close(fig)

    markers = pd.read_csv(root / "marker_definitions.csv")
    marker_table = markers.groupby("cell_type")["gene"].apply(lambda x: ", ".join(x)).rename("markers").reset_index()
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>GSE136103 CellAssign</title><style>body{{font:15px/1.5 -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;max-width:1180px;margin:32px auto;padding:0 20px;color:#202124}}h1,h2{{color:#17365d}}.cards{{display:flex;gap:12px;flex-wrap:wrap}}.card{{background:#eef4fb;padding:12px 18px;border-radius:8px}}table{{border-collapse:collapse;font-size:13px}}th,td{{padding:6px 9px;border:1px solid #ddd}}th{{background:#f4f6f8}}img{{max-width:100%;height:auto}}code{{background:#f2f2f2;padding:2px 4px}}</style></head><body>
<h1>GSE136103 CellAssign annotation</h1><div class='cards'><div class='card'><b>Cells</b><br>{metrics['cells']:,}</div><div class='card'><b>Candidate types</b><br>12</div><div class='card'><b>Published-label accuracy</b><br>{metrics['accuracy']:.3f}</div><div class='card'><b>Macro F1</b><br>{metrics['macro_f1']:.3f}</div></div>
<h2>Methodology</h2><p>CellAssign from scvi-tools 1.5.0.post1 was fitted to raw counts from the source AnnData. Size factors were each cell's full-transcriptome library size divided by the mean library size, computed before restricting to 78 marker genes. A binary marker matrix defined 12 closed-set candidate populations. Training used seed 0 and a 200-epoch ceiling; validation early stopping occurred at epoch 146 (best validation ELBO 39.917).</p>
<h2>Annotation and uncertainty</h2><img src='cellassign_umap_overview.png'><p>Maximum probability, top-two margin, and entropy are retained per cell. High probability is conditional on the supplied closed candidate set and is not independent biological validation.</p>
<h2>Comparison with published labels</h2><img src='cellassign_confusion_matrix.png'><img src='cellassign_label_counts.png'><p>Published labels are used only for post-fit comparison, not for CellAssign training. “Unassigned” published cells are excluded from accuracy, balanced accuracy, and macro-F1.</p>
<h2>Canonical marker panel</h2>{marker_table.to_html(index=False, border=0)}
<h2>Limitations</h2><ul><li>CellAssign can allocate cells only among the 12 supplied types; it cannot discover omitted or novel populations.</li><li>Low confidence can reflect mixed profiles, rare states, ambient RNA, doublets, low counts, or insufficient marker specificity.</li><li>Unexpected assignments should be reviewed by donor, condition, cluster, and coherent multi-marker expression.</li><li>The marker panel is a canonical analysis panel created for this run, not a dataset-author supplied marker file.</li></ul>
<h2>Outputs</h2><p><code>cellassign_probabilities.csv</code> contains all soft probabilities; <code>cellassign_assignments.csv</code> contains top label and uncertainty; <code>cellassign_annotated.h5ad</code> contains the marker-gene count matrix, metadata, probabilities, and UMAP.</p></body></html>"""
    (root / "index.html").write_text(html)


if __name__ == "__main__": main()
