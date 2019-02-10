# Primate CellType DGE
Scripts for RNA seq analysis.

### Usage
# 1) Run DGE analysis: 
**R CMD BATCH --vanilla Dge_NeuN_LM.R**

**R CMD BATCH --vanilla Dge_OLIG2_LM.R**

### Details
The script involves 4 steps
- **Data transformation** using *log2(CPM)*
- **QC** of tranformed counts
- **Modeling** of transformed counts based on *Linear Modeling*
- **Filtering** of the results detecting *Species specific DGE*

# 2) After DGE analysis use:
**R CMD BATCH --vanilla Make_DGE_Data.R**

### Details
The script simplify the outputs from DGE analysis

# 3) For basic visualization use: 
**R CMD BATCH --vanilla Viz_NeuN.R**

**R CMD BATCH --vanilla Viz_OLIG2.R**

### Details
The scripts contains codes for visualization based on *ggplot2*
