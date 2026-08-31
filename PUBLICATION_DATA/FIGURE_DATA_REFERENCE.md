# FIGURE_DATA_REFERENCE

Reference table linking manuscript figures to processed data files and source datasets.

## Main figure mapping

| Figure | Panel | Dataset | File mapping |
| :--- | :--- | :--- | :--- |
| Figure 2 | a | Wheat (D02266B1) | `WHEAT_MILLSTEED_2025/PLOTS/QC/PCA_ELBOW.png` |
| Figure 2 | b | Wheat (D02266B1) | `WHEAT_MILLSTEED_2025/PLOTS/QC/PRE_NORM_COUNT.png` and `WHEAT_MILLSTEED_2025/PLOTS/QC/POS_NORM_COUNT.png` |
| Figure 3 | a | Wheat (D02266B1) | `WHEAT_MILLSTEED_2025/NETWORK/tan_EDGE.txt` + `WHEAT_MILLSTEED_2025/NETWORK/tan_NODE.txt` |
| Figure 3 | b | Wheat (D02266B1) | `WHEAT_MILLSTEED_2025/NETWORK/green_EDGE.txt` + `WHEAT_MILLSTEED_2025/NETWORK/green_NODE.txt` |
| Figure 3 | c | Wheat (D02266B1) | `WHEAT_MILLSTEED_2025/NETWORK/red_EDGE.txt` + `WHEAT_MILLSTEED_2025/NETWORK/red_NODE.txt` |
| Figure 3 | d | Wheat (D02266B1) | `WHEAT_MILLSTEED_2025/NETWORK/salmon_EDGE.txt` + `WHEAT_MILLSTEED_2025/NETWORK/salmon_NODE.txt` |
| Figure 4 | a | Wheat | `WHEAT_MILLSTEED_2025/CLUSTERING/LEIDEN_CLUSTERS.png` and `WHEAT_MILLSTEED_2025/CLUSTERING/SPATIAL_LEIDEN_CLUSTERS.png` |
| Figure 4 | b | Mouse benchmark | data not included |
| Figure 4 | c | Arabidopsis benchmark | data not included |
| Figure 4 | d | Rice benchmark | data not included |

## Supplementary figure mapping

- Cluster-level supplementary spatial panels:
  - `WHEAT_MILLSTEED_2025/CLUSTERING/LEIDEN/`
  - `WHEAT_MILLSTEED_2025/CLUSTERING/LOUVAIN/`
  - `WHEAT_MILLSTEED_2025/CLUSTERING/SPATIAL_LEIDEN/`
- Volcano supplementary panels:
  - `WHEAT_MILLSTEED_2025/PLOTS/LEIDEN_VOLCANO_PLOTS_COMPLETE/`
  - `WHEAT_MILLSTEED_2025/PLOTS/LOUVAIN_VOLCANO_PLOTS_COMPLETE/`
  - `WHEAT_MILLSTEED_2025/PLOTS/SPATIAL_LEIDEN_VOLCANO_PLOTS_COMPLETE/`
- Additional supplementary summary outputs:
  - `WHEAT_MILLSTEED_2025/STATISTICAL_ANALYSIS/`
  - `WHEAT_MILLSTEED_2025/LOGS/`
  - `WHEAT_MILLSTEED_2025/REPORTS/ANALYSIS_REPORT.txt`

## Raw data policy

Raw `.gef` files are not included in this repository and are intentionally excluded from versioned publication packages.

## External dataset references

- Wheat: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298021
- Rice: https://ftp.cngb.org/pub/stomics/STT0000026/Analysis/STSA0000251/STTS0000395/
- Arabidopsis: https://db.cngb.org/stomics/datasets/STDS0000104/
- Mouse: https://db.cngb.org/stomics/datasets/STDS0000058/
