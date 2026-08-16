# Marker panels and interpretation

Use only genes present in the integrated feature set and report missing markers.

| Group | Marker genes |
|---|---|
| B cells | `CD79A` |
| Endothelial | `VWF`, `PECAM1`, `CDH5`, `ERG`, `KDR`, `TEK` |
| Fibroblasts | `DCN`, `COL1A1`, `PDGFRA`, `PCOLCE`, `THY1`, `ACKR3` |
| Keratinocytes | `KRT1`, `KRT5` |
| Langerhans cells | `CD1A` |
| Lymphatic | `CCL21`, `PDPN`, `PROX1`, `FLT4` |
| Macrophages | `CD68`, `CD163` |
| Mast cells | `TPSAB1` |
| Melanocytes | `PMEL` |
| Monocytes | `FCGR3A`, `C5AR1`, `FCN1` |
| NK cells | `NKG7` |
| Neural | `PLP1` |
| Pericytes | `RGS5`, `PDGFRB` |
| Plasma cells | `SDC1` |
| Secretory cells | `DCD`, `SCGB2A2`, `MUCL1`, `CALML5` |
| Smooth muscle | `ACTA2`, `DES`, `SMTN`, `TAGLN`, `TPM2` |
| T cells | `CD2`, `CD3D` |
| cDC | `CD1C`, `CD1B`, `CD1E` |
| pDC | `NRP1`, `IRF8` |

The source-specific comparison highlights `PECAM1` (CD31), `CD34`, and `ACTA2`.
The metabolic analysis highlights `NAMPT` and `NNMT`, including cell-type and
condition-stratified dot plots.

For SSc-versus-Healthy comparisons, aggregate raw counts by patient and cell type before
differential expression. Require adequate patients and cells per patient, include study as a
design factor when estimable, and report when study and condition cannot be separated.
