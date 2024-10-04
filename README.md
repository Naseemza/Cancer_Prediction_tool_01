# **PDB Structure Analysis Program**

## Overview

This R program performs analysis on **PDB (Protein Data Bank)** structures using entropy and energy data. It leverages various statistical and machine learning techniques to cluster and analyze structural data. The program can handle single and multi-chain PDB files, providing detailed results like silhouette scores, convex hull areas, RMSE, and SER values. Two visualizations are also generated to assist in understanding the analysis: an **Enhanced Convex Hill Plot** and a **UMAP-KMeans Clustering Plot**.

## Features

- Automatic chain selection for PDB structures with multiple chains.
- Computes important metrics like **Silhouette Score**, **Hull Area**, **PEA-Variance**, **RMSE**, **SER**, and **SER-scaled** values.
- Generates two plots:
  - `PDBid-enhanced-convex-hull-plot.png`
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

# Generating `entropy.xls` File
Follow these steps to generate the `entropy.xls` file:

- ### Step 1: Access the Packing Entropy Webserver
  - Go to the Packing Entropy Webserver ([link to the actual webserver](https://packing-entropy.bb.iastate.edu/)).
  - Input your PDB ID in the designated dialog box.
  - Complete the CAPTCHA and press the Submit button.
- ### Step 2: Download the output.txt File
  - After some processing time, the server will generate an output.txt file.
  - Download this file to your local system.
- ### Step 3: Import output.txt to Excel
  - Open Excel and import the output.txt file.
  - Choose delimited when importing and set the delimiter based on the structure of the file.
- ###BStep 4: Create a Residue Column
  - In the Excel sheet, create a new column titled residue.
- ### Step 5: Apply Formula for Bio3D-Compatible Residue Names
  - In the first cell under the residue column (e.g., `F2`), enter the following formula:
          
            =C2&"."&B2&"."&D2
  
  - This will combine the chain ID, residue position, and residue name similar to Bio3D residue format.

- ### Step 6: Save the File as `entropy_pdbid.xlsx` 
  - Once the formula is applied across all relevant rows, save the file as entropy-pdbid.xlsx, replacing pdbid with your actual PDB ID (e.g., entropy_4o33.xlsx).

# Generating `energy.xlsx` File
Follow these steps to generate the `energy.xlsx` file:

- ### Step 1: Download the PDB File
  - Go to the Protein Data Bank and download the PDB file for your desired PDB ID.
- ### Step 2: Install and Open `Swiss-PDB_Viewer`
  - Download and install `Swiss-PDB-Viewer` (SPDBV) from here ([download link](https://spdbv.unil.ch/disclaim.html#)).
  - Open the `PDB file` using `Swiss-PDB-Viewer`.
- ### Step 3: Select All Residues
  - In Swiss-PDB-Viewer, go to the Select menu and choose All residues.
- ### Step 4: Compute Energy
  - Go to the Tools menu and select `Compute Energy (Force Field)`.
- ### Step 5: Confirm Energy Calculation
  - A dialog box will appear. Click OK to proceed with the energy calculation.
- ### Step 6: Locate the Output File
  - After calculation, navigate to the installation folder of Swiss-PDB-Viewer.
  - Inside the temp folder, you will find an `energy.E1` file.
  - Open this file using `Notepad`.
- ### Step 7: Transfer Data to Excel
  - Copy all the table data from the `energy.E1` file and paste it into an `Excel sheet`.
- ### Step 8: Modify the Excel Sheet
  - Remove the columns titled:
    - `Bonds, Angles , Torsion, Improper, Electrostatic, Constraint`. 
  - Insert a column to the right of the residue column.
    - Use the formula to combine relevant columns for Bio3D compatibility in B2 cell and scroll it down till end.
      ```excel
      =IF(LEN(A2) = 3, A2, LEFT(A2, LEN(A2)-1))
      ```
  - Insert another column to the left of the Total column and name it Variable.
    - Use the following formula in cell E2:
      ```
      =D2&"."&C2&"."&B2
      ```
- Step 9: Save the Excel File as `energy_pdbid.xlsx`
  - Once the modifications are complete, save the file as `energy_pdbid.xlsx`, replacing pdbid with your actual PDB ID (e.g., `energy_4o33.xlsx`).

# Usage
- ### Step 1: Load the Program
To start the program, run the R script containing the main program code using the source() function:

```r
source("run_program.R")
```
- ### Step 2: Execute the Program
To execute the analysis, call the `run_program` function with the following parameters:

- PDBID – Protein Data Bank identifier (e.g., `"4o33"`)
- Path to the entropy file (`entropy_4o33.xls`)
- Path to the energy file (`energy_4o33.xls`)
```r
run_program(PDBID, "path/to/entropy.xls", "path/to/energy.xls")
```

# Program Behavior
- Single Chain PDB: The program will automatically execute for single-chain PDB files.
- Multiple Chains: If multiple chains are present, the user will be prompted to select a chain ID (e.g., "A" or "B").
## Output
The program produces the following outputs:
- PDB ID
- Final_result: A data frame that contains the following analysis details:
- Silhouette score
- Convex hull area
- MBR area
- Ratio of MNR_area/Hull_area
- RMSE
- SER and SER-scaled values
## Plots:
- `PDBID-enhanced-convex-hull-plot.png`: Convex hull and MBR area plot visualizing the structure.
- `PDBID-umap-kmeans-plot.png`: UMAP and k-means clustering plot.
Both the data frame and plots will be stored in the working directory for further analysis.

# Example
To analyze a PDB file, you can use the following example code:

```r
run_program("4o33", "/path/to/entropy_4o33.xls", "/path/to/energy_4o33.xls")
```
If 4o33 contains multiple chains, you will be prompted to select a chain:

```r
Select chain ID (e.g., A, B): A
```
# Analysis Results
The output Final_result can be used to check if the analyzed PDB location corresponds to known cancer types. Statistical metrics such as `silhouette scores`, `ratio` and `scale_SER` score results provide insight into protein misfolding behavior.


 
