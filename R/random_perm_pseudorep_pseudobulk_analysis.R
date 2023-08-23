library(data.table)
library(Seurat)
library(cowplot)
library(wesanderson)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)


#source functions ----------------------------------------------
source("./R/sc_cell_type_de.R")

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


#' Run differential expression pseudobulk
#' uses other functions for running DE pb normally but adapted to handle 
#' the random permutations
#'
#' @param count_mat scRNA-Seq count matrix
#' @param meta metadat for the scRNA-Seq count matrix
#' @param AD_discrim_col column in the metadata to discriminate between cases
#' and controls
#' @return A datatable of all genes and their DE values for each cell type 
#' analysed.
run_pseudobulk <- function(count_mat,meta,AD_discrim_col){
  #using 
  #make SCE object to pass to the function
  #filter metadata to the rand perm of interest
  meta[,y:=get(AD_discrim_col)]
  sce_mathys <- 
    SingleCellExperiment::SingleCellExperiment(list(counts=count_mat),
                                               colData=meta)
  #remove endothelial cells to match Mathys 2019's approach
  sce_mathys <- sce_mathys[,!sce_mathys$broad.cell.type %in% c("End","Per")]
  
  sc_cell_type_de_return <- 
    sc_cell_type_de(sce_mathys, design = ~ y, 
                    pseudobulk_ID="projid", 
                    celltype_ID="broad.cell.type",
                    folder = FALSE,#no plots
                    verbose=TRUE,get_gene_names=FALSE)
  
  DEG <- rbindlist(sc_cell_type_de_return$celltype_DEGs,idcol='celltype')
  #just so naming matches pseudoreplication for comparison
  DEG[,p_val_adj:=adj_pval]
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
run_act_pr <- FALSE
#load dataset - filtered, processed data from author's original work
Mathys_org <- 
  Matrix::readMM("data/filtered_count_matrix.mtx")
#load metadata
meta <- fread("data/filtered_column_metadata.txt")
clinicalMeta <- fread("data/ROSMAP_Clinical_2019-05_v2.csv")
idmap <- fread("data/id_mapping.csv")
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
rownames(Mathys_org)<- readLines("data/filtered_gene_row_names.txt")

#First let's check running it normally against author's results ---------------
if (run_act_pr){
  #now we can run actual ids with a pseudoreplication approach the same as the 
  #one in the paper - Wilcoxon rank-sum test and FDR multiple-testing correction
  DEG <- run_pseudoreplication(Mathys_org,meta,
                               AD_discrim_col='AD_pathology')
  #load to compare to what authors reported
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
  print(cor(DEG[!is.na(mathys_p_val_adj),]$p_val_adj,
            DEG[!is.na(mathys_p_val_adj),]$mathys_p_val_adj))
  DEG[mathys_degs,mathys_avg_log2FC:=i.lfc,on=.(celltype,gene)]
  print(cor(DEG[!is.na(mathys_avg_log2FC),]$avg_log2FC,
            DEG[!is.na(mathys_avg_log2FC),]$mathys_avg_log2FC))
  mathys_degs[DEG,our_adj_p_val:=i.p_val_adj,on=.(celltype,gene)]
  print(nrow(mathys_degs[adj_p_val<0.05]))
  print(nrow(mathys_degs[adj_p_val<0.05 & our_adj_p_val<0.05]))
  
  #save res - overlap to what author's found isn't identical
  #methods section was very sparse for this though
  fwrite(DEG,"./data/DEGs_pseudorep.csv")
}  
# -----------------------------------------------------------------------------

#Now run random permutations
rand_perm_tests <- 100
#we want to randomly shuffle the patient id's 
#set seed for reproducibility
set.seed(101) # so we can get the same random seeds again
rand_seeds <- sample(1:100000, rand_perm_tests, replace=FALSE)
rand_projids <- lapply(rand_seeds,randomly_shuffle,meta=meta)
#add as column to meta
for(i in seq_len(rand_perm_tests)){
  rand_projid_i <- paste0('rand_projid',i)
  meta[,(rand_projid_i):=rand_projids[[i]]]
  #join on info from idmap to get AD pathology
  idmap[,(rand_projid_i):=as.character(projid)]
  #add AD_pathology to meta too
  AD_rand_projid_i <- paste0('AD_pathology_rand',i)
  meta[idmap,(AD_rand_projid_i):=i.AD_pathology, on=rand_projid_i]
  #save space delete ids
  meta[,(rand_projid_i):=NULL]
  idmap[,(rand_projid_i):=NULL]
}
labels <- unique(meta$AD_pathology)

#now we can run randomly shuffled patient ids with a 
#pseudoreplication approach the same as the one in the paper -
#Wilcoxon rank-sum test and FDR multiple-testing correction
all_rand_perm_degs <- vector(mode='list',length=rand_perm_tests)
#set name as index
names(all_rand_perm_degs)<-paste0('AD_pathology_rand',seq_len(rand_perm_tests))
for(i in seq_len(rand_perm_tests)){
  print(paste0("Random Perm: Pseudoreplication ",i))
  AD_rand_projid_i <- paste0('AD_pathology_rand',i)
  DEG_i <- run_pseudoreplication(Mathys_org,meta,
                          AD_discrim_col=AD_rand_projid_i)
  all_rand_perm_degs[[i]] <- DEG_i[p_val_adj<0.05,]
}

#combine results
all_rand_perm_degs_dt <- rbindlist(all_rand_perm_degs,idcol='run')

#save results
fwrite(all_rand_perm_degs_dt,"./data/DEGs_pseudorep_rand_perm_seed101.csv")


#now run the same for pseudobulk

#source the necessary function to run pb DE 
source("./R/sc_cell_type_de.R")

all_rand_perm_degs_pb <- vector(mode='list',length=rand_perm_tests)
#set name as index
names(all_rand_perm_degs_pb)<-paste0('AD_pathology_rand',
                                     seq_len(rand_perm_tests))
for(i in seq_len(rand_perm_tests)){
  print(paste0("Random Perm: Pseudobulk ",i))
  AD_rand_projid_i <- paste0('AD_pathology_rand',i)
  all_rand_perm_degs_pb[[i]] <-
    run_pseudobulk(Mathys_org,meta,
                   AD_discrim_col=AD_rand_projid_i)
}

#combine results
all_rand_perm_degs_pb_dt <- rbindlist(all_rand_perm_degs_pb,idcol='run')

#save results
fwrite(all_rand_perm_degs_pb_dt,"./data/DEGs_pseudobulk_rand_perm_seed101.csv")




#let's plot ---------------
#plot pseudoreplication results
#load results
all_rand_perm_degs_pr_dt <- 
  fread("./data/DEGs_pseudoreplication_rand_perm_seed101.csv")
#get N count 
all_rand_perm_degs_pr_dt <- all_rand_perm_degs_pr_dt[,.N,by=.(run,celltype)]
runs <- sort(unique(all_rand_perm_degs_pr_dt$run))
#plot correlation between cell counts and rand perm DEGs
#get counts
mathys_cell_counts <- data.table::fread("data/Mathys_cell_counts.csv")
setnames(mathys_cell_counts,"cell","celltype")
#add to dt
all_rand_perm_degs_pr_dt[mathys_cell_counts,cell_counts:=i.count,on="celltype"]

#get corr for each run
rand_perm_cors_pr <- vector(mode='list',length=length(runs))
names(rand_perm_cors_pr) <- runs
for(perm_i in runs){
  #can only test if more than two values
  if (nrow(all_rand_perm_degs_pr_dt[run==perm_i,])>2){
    cor_i <- cor.test(all_rand_perm_degs_pr_dt[run==perm_i,]$cell_counts,
                      all_rand_perm_degs_pr_dt[run==perm_i,]$N,
                      method = "pearson")
    rand_perm_cors_pr[[perm_i]] <- cor_i$estimate
  }
}
#get rid of any null entries
rand_perm_cors_pr[sapply(rand_perm_cors_pr, is.null)] <- NULL
comp_mathys_text_rand_perm_pr <- paste0('mean r = ', 
                                        round(mean(unlist(rand_perm_cors_pr)),2))

#plot correlation between count cells and number DEGs
pal=c(wesanderson::wes_palette("Royal2"),
      wesanderson::wes_palette("Moonrise3")[1])

rand_perm_fig_pr <-
  ggplot(data=all_rand_perm_degs_pr_dt,
         aes(x = cell_counts, 
             y = N,colour=celltype))+ 
  geom_smooth(data=all_rand_perm_degs_pr_dt,
              method = "lm", se = T,formula='y ~ x',
              color = "grey20", alpha = 0.1,fill="lightgrey") +  
  geom_boxplot(outlier.colour = NULL,outlier.shape = NA)+
  #geom_point(size=3,color=pal) +
  annotate('text', x = -Inf, y = Inf, hjust = -4.35, vjust = 1.1, 
           label=comp_mathys_text_rand_perm_pr) +
  scale_x_continuous(labels = scales::label_number(suffix = "K", 
                                                   scale = 1e-3))+
  ylim(0,60)+
  labs(y= "Number of DEGs", x = "Cell Counts",colour="") +
  ggtitle("Pseudobulk - 100 random permutations")+
  theme_cowplot()+
  theme(axis.text = element_text(size=9),
        plot.title = element_text(size = 13, face = "plain"),
        axis.title = element_text(size=11))+
  scale_colour_manual(values=pal)

#--------------------------
#now plot pseudobulk results
#load results
all_rand_perm_degs_pb_dt <- 
  fread("./data/DEGs_pseudobulk_rand_perm_seed101.csv")
#get N count 
all_rand_perm_degs_pb_dt <- all_rand_perm_degs_pb_dt[,.N,by=.(run,celltype)]
runs <- sort(unique(all_rand_perm_degs_pb_dt$run))
#plot correlation between cell counts and rand perm DEGs
#get counts
mathys_cell_counts <- data.table::fread("data/Mathys_cell_counts.csv")
setnames(mathys_cell_counts,"cell2","celltype")
#add to dt
all_rand_perm_degs_pb_dt[mathys_cell_counts,cell_counts:=i.count,on="celltype"]

#get corr for each run
rand_perm_cors_pb <- vector(mode='list',length=length(runs))
names(rand_perm_cors_pb) <- runs
for(perm_i in runs){
  #can only test if more than two values
  if (nrow(all_rand_perm_degs_pb_dt[run==perm_i,])>2){
    cor_i <- cor.test(all_rand_perm_degs_pb_dt[run==perm_i,]$cell_counts,
                      all_rand_perm_degs_pb_dt[run==perm_i,]$N,
                      method = "pearson")
    rand_perm_cors_pb[[perm_i]] <- cor_i$estimate
  }
}
#get rid of any null entries
rand_perm_cors_pb[sapply(rand_perm_cors_pb, is.null)] <- NULL
comp_mathys_text_rand_perm_pb <- paste0('mean r = ', 
                                     round(mean(unlist(rand_perm_cors_pb)),2))

#plot correlation between count cells and number DEGs
pal=c(wesanderson::wes_palette("Royal2"),
      wesanderson::wes_palette("Moonrise3")[1])

rand_perm_fig_pb <-
  ggplot(data=all_rand_perm_degs_pb_dt,
         aes(x = cell_counts, 
             y = N,colour=celltype))+ 
  geom_smooth(data=all_rand_perm_degs_pb_dt,
              method = "lm", se = T,formula='y ~ x',
              color = "grey20", alpha = 0.1,fill="lightgrey") +  
  geom_boxplot(outlier.colour = NULL,outlier.shape = NA)+
  #geom_point(size=3,color=pal) +
  annotate('text', x = -Inf, y = Inf, hjust = -4.35, vjust = 1.1, 
           label=comp_mathys_text_rand_perm_pb) +
  scale_x_continuous(labels = scales::label_number(suffix = "K", 
                                                   scale = 1e-3))+
  ylim(0,60)+
  labs(y= "Number of DEGs", x = "Cell Counts",colour="") +
  ggtitle("Pseudobulk - 100 random permutations")+
  theme_cowplot()+
  theme(axis.text = element_text(size=9),
        plot.title = element_text(size = 13, face = "plain"),
        axis.title = element_text(size=11))+
  scale_colour_manual(values=pal)

