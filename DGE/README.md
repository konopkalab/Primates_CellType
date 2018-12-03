# 1) Run DGE analysis: 
R CMD BATCH --vanilla Dge_NeuN_LM.R

R CMD BATCH --vanilla Dge_OLIG2_LM.R

# 2) After DGE analysis use:
R CMD BATCH --vanilla Make_DGE_Data.R

# 3) For basic visualization use: 
R CMD BATCH --vanilla Viz_NeuN.R

R CMD BATCH --vanilla Viz_OLIG2.R

