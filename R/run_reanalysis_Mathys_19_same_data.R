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
library(SingleCellExperiment)

#load dataset and exclude endothelial cells to match Tsai group analysis
#load dataset - filtered, processed data from author's original work
Mathys_org <- 
  Matrix::readMM("../data/filtered_count_matrix.mtx")
#load metadata
meta <- fread("../data/filtered_column_metadata.txt")
clinicalMeta <- fread("../data/ROSMAP_Clinical_2019-05_v2.csv")
idmap <- fread("../data/id_mapping.csv")
idmap[,fastq:=NULL]
idmap <- unique(idmap)
#we know ROS 1-24 no AD, >24 AD from Manuscript
idmap[ , Subject2 := as.numeric(gsub(".*ROS?(\\d+).*", "\\1", Subject))]
idmap[Subject2>24,AD_pathology:='AD_pathology']
idmap[Subject2<=24,AD_pathology:='non-AD_pathology']
#join to clinical features
clinicalMeta[idmap,AD_pathology:=i.AD_pathology, on="projid"]
#filter to 48 patients
clinicalMeta <- clinicalMeta[!is.na(AD_pathology)]
#I checked values against supp table 1 and they match
setorder(clinicalMeta,AD_pathology)
#add AD_pathology to meta too
meta[idmap,AD_pathology:=i.AD_pathology, on="projid"]
#add data to count matrix
colnames(Mathys_org)<-meta$TAG
rownames(Mathys_org)<- readLines("../data/filtered_gene_row_names.txt")
#create SCE object to hold meta data and expression data
sce_mathys <- 
  SingleCellExperiment::SingleCellExperiment(list(counts=Mathys_org),
                                             colData=meta)
#remove endothelial cells to match Mathys 2019's approach
sce_mathys <- sce_mathys[,!sce_mathys$broad.cell.type %in% c("End","Per")]

sc_cell_type_de_return <- 
  sc_cell_type_de(sce_mathys, design = ~ AD_pathology, 
                  pseudobulk_ID="projid", 
                  celltype_ID="broad.cell.type",
                  folder = "../results/sc_cell_type_de_graphs_same_data/",
                  control = "non-AD_pathology",
                  verbose=TRUE)

#save results
#Cell Counts
data.table::fwrite(data.table(cell=names(
  sc_cell_type_de_return$celltype_counts),
  count=sc_cell_type_de_return$celltype_counts),
  "../results/cell_counts_same_data.csv")
#DEGs
data.table::fwrite(
  data.table::rbindlist(sc_cell_type_de_return$celltype_DEGs,idcol='Cell'),
  "../results/DEGs_same_data.csv")
#All Genes
data.table::fwrite(
  data.table::rbindlist(sc_cell_type_de_return$celltype_all_genes,idcol='Cell'),
  "../results/all_genes_same_data.csv")
#R object
save(sc_cell_type_de_return,file="../results/Mathys_pb_res_same_data.RData")
