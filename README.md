# **PDB Structure Analysis Program**

## Overview

This R program performs analysis on **PDB (Protein Data Bank)** structures using entropy and energy data. It leverages various statistical and machine learning techniques to cluster and analyze structural data. The program can handle single and multi-chain PDB files, providing detailed results like silhouette scores, convex hull areas, RMSE, and SER values. Two visualizations are also generated to assist in understanding the analysis: an **Enhanced Convex Hill Plot** and a **UMAP-KMeans Clustering Plot**.

## Features

- Automatic chain selection for PDB structures with multiple chains.
- Computes important metrics like **Silhouette Score**, **Hull Area**, **PEA-Variance**, **RMSE**, **SER**, and **SER-scaled** values.
- Generates two plots:
  - `PDBid-enhanced-convex-hill-plot.png`
  - `PDBid-umap-kmeans-plot.png`
- Supports further biological and cancer-related research using PDB analysis.

## Prerequisites

### Required R Libraries

The following R packages are required to run the program. You can install them using the command below:

```r
install.packages(c("bio3d", "cluster", "ggplot2", "umap", "readxl", "dplyr", 
                   "geometry", "MASS", "FactoMineR", "dbscan", "sf", 
                   "ggrepel", "neuralnet", "gridExtra", "reshape2"))
```
## Package Descriptions:
- `bio3d`: For PDB structure analysis.
- `cluster`: For clustering algorithms.
- `ggplot2`: For data visualization.
- `umap`: For dimension reduction and clustering.
- `readxl`: For reading Excel files.
- `dplyr`: For data manipulation.
- `geometry`: For convex hull and geometry-based calculations.
- `MASS`: For Mahalanobis distance calculation.
- `FactoMineR`: For PCA analysis.
- `dbscan`: For DBSCAN clustering.
- `sf`: For minimum bounding rectangle (MBR).
- `ggrepel`: For better label placement in plots.
- `neuralnet`: For neural network computations.
- `gridExtra`: For arranging multiple ggplots.
- `reshape2`: For data reshaping.




 
