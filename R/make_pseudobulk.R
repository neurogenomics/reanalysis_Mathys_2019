# make_pseudobulk
# Calculate the summed pseudobulk values for an SCE object based on one 
# single cell type only. Ensure to filter SCE to pass one cell type's data.
#' @param SCE SingleCellExperiment object, a specialised S4 class for storing 
#' data from single-cell experiments. 
#' @param design the design formula of class type `formula`. Equation used to 
#' fit the model- data for the generalised linear model.
#' @param pseudobulk_ID the column name in the SCE object to perform pseudobulk 
#' on, usually the patient identifier. This column is used for grouping in the 
#' pseudobulk approach
#' @param pb_columns vector, list of annotation column names in the SCE object 
#' to be returned in annot_pb rolled up to pseudobulk level. Default is NULL 
#' which won't return any information in annot_pb.
#' @param region The column name in the SCE object for the study region. If 
#' there are multiple regions in the study (for example two brain regions). 
#' Pseudobulk values can be derived separately. Default is "single_region" 
#' which will not split by region.
#' @param rmv_zero_count_genes whether genes with no count values in any cell 
#' should be removed. Default is TRUE
#' @return A list containing
#' \itemize{
#'   \item \code{sumDat}: matrix of the summed pseudobulk count values
#'   \item \code{annot_pb}: dataframe of the annotation data from the SCE 
#'   rolled up based on the pseudobulk aggregation.
#' }
#' @export
make_pseudobulk <- function(data,pseudobulk_ID, pb_columns=NULL,
                            region="single_region",rmv_zero_count_genes=TRUE){
    allAnnot <- colData(data)
    if(region=="single_region") #constant value for all samples
        allAnnot[[region]] <- "one_region"
    indvs <- as.character(unique(allAnnot[[pseudobulk_ID]]))
    regions <- as.character(unique(allAnnot[[region]]))
    sumDat <- 
        matrix(0,nrow=dim(data)[1],
               ncol=length(indvs)*length(regions))
    colnames(sumDat) <- rep(indvs,length(regions))
    rownames(sumDat) <- rownames(data)
    count=0
    #for a single cell type, look at a single individual in single brain region
    #create list to hold annot values for each region, individual combination 
    annot_list <- vector(mode="list",length=length(regions)*length(indvs))
    names(annot_list) <- paste0(indvs,"_",regions)
    for(region_i in regions){
        for(indv in indvs){
            whichCells = 
                allAnnot[[pseudobulk_ID]]==indv & allAnnot[[region]]==region_i
            theData <- data[,whichCells]
            count <- count+1
            sumDat[,count] <- Matrix::rowSums(counts(theData))
            colnames(sumDat)[count] = sprintf("%s_%s",indv,region_i)
            #Get annotation data
            print_warning <- FALSE
            if(!is.null(pb_columns)){
                annot_i <- colData(theData)
                #there should only be a single value for each variable since 
                #this is single cell type, brain region and person
                #throw warning if not the case for numeric & aggregate 
                #accordingly. Throw error if categorical
                #first restrict data to just those from the design formula
                annot_i <- annot_i[,names(annot_i) %in% pb_columns]
                #change all factor variables to characters
                facts <- sapply(annot_i, is.factor)
                annot_i[facts] <- lapply(annot_i[facts], as.character)
                #check if values are unique
                cols_sum <- 
                    names(annot_i)[sapply(annot_i, 
                                          function(x) length(unique(x))>1)]
                if(length(cols_sum)>=1){
                    print_warning <- TRUE
                    cols_sum_warn <- cols_sum
                    #get numeric columns and categorical separately to deal with 
                    num_cols <- 
                        unlist(lapply(annot_i[,cols_sum,drop = FALSE], 
                                        is.numeric))  
                    #if any categorical throw error, shouldn't be lower level
                    cat_cols <- 
                        names(annot_i[,cols_sum,drop=FALSE][,!num_cols,drop=FALSE])
                    if(length(cat_cols)>0)
                        stop(paste0(c("Your design for the DE analysis contain",
                                      "s ",length(cat_cols)," categorical",
                                      " variable(s) which don't have a unique",
                                      " value at the celltype-region-individua",
                                      "l level. The non-unique column(s) are: ",
                                      cat_cols)))
                    annot_i_pb <- c(
                        as.list(unique(annot_i[,!(names(annot_i) %in% 
                                                        cols_sum)])),
                        sapply(annot_i[,cols_sum][,num_cols], 
                                function(x) mean(x))
                    )
                    #rearrange order to input order
                    annot_i_pb<-annot_i_pb[names(annot_i)]
                }
                else{
                    #all values unique, just get unique values
                    annot_i_pb <- as.list(unique(annot_i))
                }
                annot_list[[paste0(indv,"_",region_i)]] <- unlist(annot_i_pb)
            }    
        }
    }
    if(print_warning){
        print(paste0("Warning: The following cell level data were aggregated ", 
                     "to brain region, cell type, individual level where ",
                     "the mean for these numeric values will be ",
                     "taken:"))
        print(cols_sum_warn)
    }
    #create annotation dataframe from list
    annot_df <- as.data.frame(do.call(rbind, annot_list))
    annot_df$group_sample <- rownames(annot_df)
    rownames(annot_df) <- NULL
    #remove genes with 0 counts
    if(rmv_zero_count_genes)
        sumDat <- sumDat[rowSums(sumDat)!=0,]
    
    return(list("sumDat"=sumDat,"annot_pb"=annot_df))
}
