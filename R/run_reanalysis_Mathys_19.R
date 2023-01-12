#Code to rerun Mathys 2019 scRNA-Seq analysis differential expression on the
#reprocessed data using a pseudobulk approach

#Ran on R 4.2
#Get package dependencies
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("RcppParallel")) install.packages("RcppParallel")
if (!require("pacman")) install.packages("pacman")
if (!require("qs")) install.packages("qs")
if (!require("SingleCellExperiment")) 
  BiocManager::install("SingleCellExperiment")
if (!require("biomaRt")) BiocManager::install("biomaRt")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("cowplot")) install.packages("cowplot")
if (!require("viridis")) install.packages("viridis")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("Hmisc")) install.packages("Hmisc")
if (!require("edgeR")) BiocManager::install("edgeR")
if (!require("stats")) install.packages("stats")
if (!require("data.table")) install.packages("data.table")
if (!require("reshape2")) install.packages("reshape2")
if (!require("wesanderson")) install.packages("wesanderson")
if (!require("EnsDb.Hsapiens.v79")){
  #This package uses dependencies which are large so increase time for download
  options(timeout=2000)
  BiocManager::install("EnsDb.Hsapiens.v79")
}  

#source the necessary function to run pb DE 
source("sc_cell_type_de.R")

#library necessary packages
library(ggplot2)
library(data.table)

#load dataset and exclude endothelial cells to match Tsai group analysis
Mathys_SCE <- 
  qs::qread("../data/sce.qs")
#remove endothelial cells to match Mathys 2019's approach
Mathys_SCE <- Mathys_SCE[,!Mathys_SCE$allan_celltype=="Endo"]

#run DE analysis
sc_cell_type_de_return <- 
  sc_cell_type_de(Mathys_SCE, design = ~ pathological_diagnosis_original, 
                  pseudobulk_ID="sample_id", celltype_ID="allan_celltype",
                  folder = "../results/sc_cell_type_de_graphs/",
                  verbose=TRUE)
#TO DO - add gene name to DEGs in results.......
#save results
#Cell Counts
data.table::fwrite(data.table(cell=names(
  sc_cell_type_de_return$celltype_counts),
  count=sc_cell_type_de_return$celltype_counts),
  "../results/cell_counts.csv")
#DEGs
data.table::fwrite(
  data.table::rbindlist(sc_cell_type_de_return$celltype_DEGs,idcol='Cell'),
  "../results/DEGs.csv")
#All Genes
data.table::fwrite(
  data.table::rbindlist(sc_cell_type_de_return$celltype_all_genes,idcol='Cell'),
  "../results/all_genes.csv")
#R object
save(sc_cell_type_de_return,file="../results/Mathys_pb_res.RData")
