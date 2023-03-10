library(ggplot2)
library(data.table)
library(Seurat)


#function to run pseudoreplication analysis --------------------

#' Run differential expression pseudoreplication analysis using Wilcoxon
#' rank-sum test and multiple test corrections with FDR approach
#'
#' @param count_mat scRNA-Seq count matrix
#' @param meta metadat for the scRNA-Seq count matrix
#' @param AD_discrim_col column in the metadata to discriminate between cases
#' and controls
#' @return A datatable of all genes and their DE values for each cell type 
#' analysed.
run_pseudoreplication <- function(count_mat,meta,AD_discrim_col){
  #Wilcoxon rank-sum test and FDR multiple-testing correction
  #implement using Seurat findmarkers function
  rownames(meta) <- colnames(count_mat)
  sc <- Seurat::CreateSeuratObject(count_mat, meta.data = meta)
  Idents(sc) <- sc$broad.cell.type
  #log normalise data
  sc <- NormalizeData(sc)
  #loop through cell types to find DEGs
  DEG = list()
  for (cell_type_i in seq_along(unique(meta$broad.cell.type))) {
    cell <- unique(meta$broad.cell.type)[cell_type_i]
    #ignore Endothelial/pericytes cells to match authors' approach
    if (!cell %in% c("End","Per")){
      print(cell)
      # subset to the right cell type
      sc_cell_i <- subset(sc,idents = cell)
      # NOTE: CURRENTLY NOT THRESHOLDING BY NO. EXPRESSING GENES OR NO. CELLS
      #keep <- rowSums(sc_cell_i) >= min_exp_genes
      #sc_cell_i <- sc_cell_i[keep,]
      # run DE analysis
      res <- FindMarkers(sc_cell_i, ident.1 = labels[2], #AD pathology
                         ident.2 = labels[1], #non-AD pathology
                         assay = 'RNA', min.pct = -Inf, 
                         min.cells.feature = 0,
                         #in.cells.group = min_cells, #number cells minimum
                         logfc.threshold = -Inf,
                         group.by = AD_discrim_col, 
                         subset.ident = cell,
                         test.use = 'wilcox')
      #adj p-value base on bonferroni, need FDR
      res$p_val_adj <- p.adjust(res$p_val_adj, method = 'BH')
      res$gene <- rownames(res)
      DEG[[cell]] <- res
    }
  }
  
  DEG <- rbindlist(DEG,idcol='celltype')
  DEG <- as.data.table(DEG)
  return(DEG)
}

#' Shuffle id's at a patient level
#'
#' @param meta metadat for the scRNA-Seq count matrix
#' @param seed set the seed with this value to get reproducible results
#' and controls
#' @return vector of patient id's to be added to meta
randomly_shuffle <- function(meta,seed){
  #we want to randomly shuffle the patient id's 
  #set seed for reporducibility
  set.seed(seed)
  pat_ids <- unique(meta$projid)
  shuff_pat_ids <- sample(pat_ids)
  #swap out old id with rand picked new one
  dict <- as.list(shuff_pat_ids)
  #add char so no reassigning
  pat_ids <- paste0(pat_ids,'*')
  names(dict) <- pat_ids
  vals <- meta$projid
  #add char so no reassigning
  pat_ids <- paste0(pat_ids,'*')
  vals <- paste0(vals,'*')
  #update ids
  for (i in 1:length(dict)){
    vals <- replace(vals, vals == names(dict[i]), dict[[i]])
  }
  return(vals)
}

# -------------------------------------------------------------

#load dataset 
Mathys_org <- 
  Matrix::readMM("data/filtered_count_matrix.mtx")#"../data/filtered_count_matrix.mtx")
#load metadata
meta <- fread("data/filtered_column_metadata.txt")#"../data/filtered_column_metadata.txt")
clinicalMeta <- fread("data/ROSMAP_Clinical_2019-05_v2.csv")#"../data/ROSMAP_Clinical_2019-05_v2.csv")
idmap <- fread("data/id_mapping.csv")#"../data/id_mapping.csv")
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
rownames(Mathys_org)<- readLines("data/filtered_gene_row_names.txt")#../data/filtered_gene_row_names.txt)
#we want to randomly shuffle the patient id's 
#set seed for reporducibility
#do this three times
vals <- randomly_shuffle(meta,seed=101)
vals2 <- randomly_shuffle(meta,seed=1)
vals3 <- randomly_shuffle(meta,seed=9999)
#add to meta
meta[,rand_projid:=vals]
meta[,rand_projid2:=vals2]
meta[,rand_projid3:=vals3]
#join on info from idmap
idmap[,rand_projid:=as.character(projid)]
idmap[,rand_projid2:=as.character(projid)]
idmap[,rand_projid3:=as.character(projid)]
#add AD_pathology to meta too
meta[idmap,AD_pathology_rand:=i.AD_pathology, on="rand_projid"]
meta[idmap,AD_pathology_rand2:=i.AD_pathology, on="rand_projid2"]
meta[idmap,AD_pathology_rand3:=i.AD_pathology, on="rand_projid3"]
labels <- unique(meta$AD_pathology)

#now we can run both actual ids and randomly shuffled patient ids with a 
#pseudoreplication approach the same as the one in the paper -
#Wilcoxon rank-sum test and FDR multiple-testing correction
DEG <- run_pseudoreplication(Mathys_org,meta,
                             AD_discrim_col='AD_pathology')
DEG_rand_perm <- run_pseudoreplication(Mathys_org,meta,
                                       AD_discrim_col='AD_pathology_rand')
DEG_rand_perm2 <- run_pseudoreplication(Mathys_org,meta,
                                        AD_discrim_col='AD_pathology_rand2')
DEG_rand_perm3 <- run_pseudoreplication(Mathys_org,meta,
                                        AD_discrim_col='AD_pathology_rand3')

#compare to what authors reported
mathys_degs <- data.table::fread("data/Mathys_DEGs.csv") #../data/Mathys_DEGs.csv")
mathys_degs
#convert cell names to match
DEG[celltype=="Ex",celltype:="Exc"]
DEG[celltype=="Oli",celltype:="Oligo"]
DEG[celltype=="In",celltype:="Inh"]
DEG[celltype=="Mic",celltype:="Micro"]
DEG[celltype=="Opc",celltype:="OPC"]
DEG[celltype=="Ast",celltype:="Astro"]
#join to ours
DEG[mathys_degs,mathys_p_val_adj:=i.adj_p_val,on=.(celltype,gene)]
cor(DEG[!is.na(mathys_p_val_adj),]$p_val_adj,
    DEG[!is.na(mathys_p_val_adj),]$mathys_p_val_adj)
DEG[mathys_degs,mathys_avg_log2FC:=i.lfc,on=.(celltype,gene)]
cor(DEG[!is.na(mathys_avg_log2FC),]$avg_log2FC,
    DEG[!is.na(mathys_avg_log2FC),]$mathys_avg_log2FC)
mathys_degs[DEG,our_adj_p_val:=i.p_val_adj,on=.(celltype,gene)]
nrow(mathys_degs[adj_p_val<0.05])
nrow(mathys_degs[adj_p_val<0.05 & our_adj_p_val<0.05])



