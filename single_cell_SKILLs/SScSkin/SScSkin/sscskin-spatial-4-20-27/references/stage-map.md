# Spatial stages 20–27

The inclusive numeric range contains eight notebooks.

| Stage | Notebook | Language | Purpose |
|---|---|---|---|
| 20 | `20_function_analysis_of_niches.ipynb` | Python | Add stage-15 niche labels to spatial expression, rank niche markers, score multiple functional gene-set collections with decoupler, and save activity H5AD files. |
| 21 | `21_viz_function_of_niches.ipynb` | Python | Visualize pathway/source activities by niche, rank sources, and inspect collections including Hallmark, Reactome, KEGG, GO, WikiPathways, TF targets, and perturbation signatures. |
| 22 | `22_molecular_niches.ipynb` | Python | Integrate spot-level molecular expression, correct PCA with Harmony by `patient_id`, cluster molecular niches at Leiden resolution 0.3, compare patient niche proportions, and rank niche markers. |
| 23 | `23_cellchat.ipynb` | R | Run spatial CellChat per sample using cell2location-derived Seurat objects, tissue coordinates, Visium scale factors, the human database, truncated-mean probabilities, and minimum group size 10. |
| 24 | `24_viz_cellchat.ipynb` | R | Visualize stage-23 per-sample CellChat networks and pathway summaries. |
| 25 | `25_cellchat.ipynb` | R | Run a second CellChat branch on the later cell2location/visualization outputs. |
| 26 | `26_viz_cellchat.ipynb` | R | Visualize stage-25 CellChat objects across samples. |
| 27 | `27_decoupler.ipynb` | Python | Calculate and save per-sample pathway and transcription-factor activity AnnData objects. |

Stage 20 depends on the composition-niche labels produced by stage 15. Stage 22 instead
defines molecular niches from expression. Stages 23–26 depend on Seurat objects containing
spatial expression, coordinates, and predicted cell-type labels from earlier cell2location
stages. Stage 27 operates on filtered per-sample spatial AnnData.
