# Workflow map

| Stage | Kernel | Purpose | Main input | Main output |
|---|---|---|---|---|
| 01 | Python | Select endothelial cells, restrict patients, build neighbors/UMAP, Leiden at 0.5, assign Endo1–Endo8, rank markers, compare proportions | integrated scRNA `integrated.h5ad` | `endothelial.h5ad`, top marker CSV |
| 02 | Python | Subcluster selected endothelial compartments labeled Endo6 and Endo5 | stage-01 H5AD | exploratory LEC/VEC subcluster plots and markers |
| 03 | Python | Visualize functional endothelial gene programs by exploratory clusters | stage-01 H5AD | dot plots/functional summaries |
| 04 | Python | Map Endo1–Endo8 into biological states, rerun markers and patient composition | stage-01 H5AD | consolidated `endothelial.h5ad`, `endo_prop_in_patient.csv` |
| 05 | Python | Export Scanpy Wilcoxon markers for consolidated states | stage-04 H5AD | `markers.csv` |
| 06 | R | Convert/read endothelial H5AD in Seurat and export positive markers | stage-04 H5AD | `markers_endo.csv` |

## Manual mapping in stage 04

| Exploratory cluster | Consolidated state |
|---|---|
| Endo1 | EndoMT-2 |
| Endo2 | VenEC |
| Endo3 | VenEC |
| Endo4 | ArtEC |
| Endo5 | VenEC |
| Endo6 | CapEC |
| Endo7 | LEC |
| Endo8 | EndoMT-1 |

Stage 02 names its Endo6 subset `adata_LEC` and Endo5 subset `adata_VEC`, whereas stage 04 later
maps Endo6 to CapEC and Endo5 to VenEC. Treat this mismatch as an upstream exploratory naming
artifact and resolve identities from marker evidence rather than variable names.
