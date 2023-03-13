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

#First let's check running it normally against author's results ---------------

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
fwrite(DEG,"./data/DEGs_pseudorep.csv")#"../data/DEGs_pseudorep.csv")
# -----------------------------------------------------------------------------

#Now run random permutations
rand_perm_tests <- 100
#we want to randomly shuffle the patient id's 
#set seed for reporducibility
#do this three times
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
  AD_rand_projid_i <- paste0('AD_pathology_rand',i)
  all_rand_perm_degs[[i]] <-
    run_pseudoreplication(Mathys_org,meta,
                          AD_discrim_col=AD_rand_projid_i)
}

#combine results
all_rand_perm_degs_dt <- rbindlist(all_rand_perm_degs,idcol='run')

#save results
fwrite(DEG_rand_perm,"./data/DEGs_pseudorep_rand_perm_seed101.csv")#"../data/DEGs_pseudorep_rand_perm_seed101.csv")





#number of DEGs
print(paste0("Number of DEGs based on AD pathology: ",
             nrow(DEG[p_val_adj<0.05])))
print(paste0("Number of DEGs based on random permutation AD pathology: ",
             nrow(DEG_rand_perm[p_val_adj<0.05])))
print(paste0("Number of DEGs based on random permutation AD pathology: ",
             nrow(DEG_rand_perm2[p_val_adj<0.05])))
print(paste0("Number of DEGs based on random permutation AD pathology: ",
             nrow(DEG_rand_perm3[p_val_adj<0.05])))

#plot correlation between cell counts and DEGs
#get DEG counts
dist_deg <- DEG[p_val_adj<0.05,.N,by=celltype]$N
names(dist_deg) <- DEG[p_val_adj<0.05,.N,by=celltype]$celltype
dist_deg_rand_perm <- DEG_rand_perm[p_val_adj<0.05,.N,by=celltype]$N
names(dist_deg_rand_perm) <- 
  DEG_rand_perm[p_val_adj<0.05,.N,by=celltype]$celltype
dist_deg_rand_perm2 <- DEG_rand_perm2[p_val_adj<0.05,.N,by=celltype]$N
names(dist_deg_rand_perm2) <- 
  DEG_rand_perm2[p_val_adj<0.05,.N,by=celltype]$celltype
#missing Mic, add 0
dist_deg_rand_perm3 <- c(DEG_rand_perm3[p_val_adj<0.05,.N,by=celltype]$N,0)
names(dist_deg_rand_perm3) <- 
  c(DEG_rand_perm3[p_val_adj<0.05,.N,by=celltype]$celltype,'Mic')
#get counts
mathys_cell_counts <- data.table::fread("./data/Mathys_cell_counts.csv")#../data/Mathys_cell_counts.csv")
#make a list
mathys_all_gene_counts <- mathys_cell_counts$count
names(mathys_all_gene_counts) <- mathys_cell_counts$cell
#add DEG counts for perms and norm
rand_perm_mathys <- 
  data.table("celltype"=sort(names(mathys_all_gene_counts)),
             "cell_counts"=mathys_all_gene_counts[
               order(names(mathys_all_gene_counts))],
             "degs"=dist_deg[order(names(dist_deg))],
             "degs_rand_perm"=dist_deg_rand_perm[
               order(names(dist_deg_rand_perm))],
             "degs_rand_perm2"=dist_deg_rand_perm2[
               order(names(dist_deg_rand_perm2))],
             "degs_rand_perm3"=dist_deg_rand_perm3[
               order(names(dist_deg_rand_perm3))]
             )
rand_perm_corr_txts <- vector(mode='list',
                              length=4)
cor_mathys <- 
  cor.test(rand_perm_mathys$cell_counts,
           rand_perm_mathys$degs,method = "pearson")
comp_mathys_text <- paste0('r = ', round(cor_mathys$estimate,2))
rand_perm_corr_txts[[1]]<-comp_mathys_text
cor_mathys_rand_perm <- 
  cor.test(rand_perm_mathys$cell_counts,
           rand_perm_mathys$degs_rand_perm,method = "pearson")
comp_mathys_text_rand_perm <- paste0('r = ', 
                                     round(cor_mathys_rand_perm$estimate,2))
rand_perm_corr_txts[[2]]<-comp_mathys_text_rand_perm
cor_mathys_rand_perm2 <- 
  cor.test(rand_perm_mathys$cell_counts,
           rand_perm_mathys$degs_rand_perm2,method = "pearson")
comp_mathys_text_rand_perm2 <- paste0('r = ', 
                                     round(cor_mathys_rand_perm2$estimate,2))
rand_perm_corr_txts[[3]]<-comp_mathys_text_rand_perm2
cor_mathys_rand_perm3 <- 
  cor.test(rand_perm_mathys$cell_counts,
           rand_perm_mathys$degs_rand_perm3,method = "pearson")
comp_mathys_text_rand_perm3 <- paste0('r = ', 
                                     round(cor_mathys_rand_perm3$estimate,2))
rand_perm_corr_txts[[4]]<-comp_mathys_text_rand_perm3
#wide to long
setnames(rand_perm_mathys,c('celltype','cell_counts',
                            'No random permutation','Random permutation 1',
                            'Random permutation 2','Random permutation 3'))
rand_perm_mathys_long <-
  data.table::melt(rand_perm_mathys,
                   id.vars = c('celltype', 'cell_counts')
                   )

#plot correlation between count cells and number DEGs
library(data.table)
library(cowplot)
library(wesanderson)
library(ggplot2)
library(patchwork)

pal=c(wesanderson::wes_palette("Royal2"),
      wesanderson::wes_palette("Moonrise3")[1])

rand_perms_var <- unique(rand_perm_mathys_long$variable)
rand_perm_plts <- vector(mode='list',
                         length=length(rand_perms_var))
for(i in seq_along(rand_perms_var)){
  rand_perm_i <- rand_perms_var[[i]]
  #exclude pal 2 and 4 - used for main two colours
  pal_i <- pal[-4][-2][[i]]
  rand_perm_plts[[i]]<-
    ggplot(data=rand_perm_mathys_long[variable==rand_perm_i,],
         aes(x = cell_counts, 
             y = value,colour=variable))+ 
    geom_smooth(data=rand_perm_mathys_long[variable==rand_perm_i,],
                method = "lm", se = T,formula='y ~ x',
                color = pal_i, alpha = 0.4,fill="white") +  
    geom_point(size=3,color=pal_i) +
    annotate('text', x = -Inf, y = Inf, hjust = -1.7, vjust = 1, 
             label=rand_perm_corr_txts[[i]],color = pal_i) +
    scale_y_continuous(labels = scales::label_number(suffix = "K", 
                                                     scale = 1e-3))+
    scale_x_continuous(labels = scales::label_number(suffix = "K", 
                                                     scale = 1e-3))+
    labs(y= "Number of DEGs", x = "Cell Counts",colour="") +
    ggtitle(paste0("pseudoreplication - ", rand_perm_i))+
    theme_cowplot()+
    theme(axis.text = element_text(size=9),
          plot.title = element_text(size = 13, face = "plain"),
          axis.title = element_text(size=11),
          legend.position="none")
}

combn_rnd_perm_fig <- 
  rand_perm_plts[[1]]+rand_perm_plts[[2]]+
  rand_perm_plts[[3]]+rand_perm_plts[[4]]
