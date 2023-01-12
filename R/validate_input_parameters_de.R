#' Validate input parameters differential expression analysis
#' Internal function called by sc_cell_type_de()
#'
#' Validate that there are no user input errors for the differential expression 
#' analysis - sc.cell.type.de
#' @param SCE SingleCellExperiment object, a specialised S4 class for storing 
#' data from single-cell experiments
#' @param design the design formula of class type `formula`. Equation used to 
#' fit the model- data for the generalised linear model.
#' @param pseudobulk_ID the column name in the SCE object to perform pseudobulk 
#' on, usually the patient identifier. This column is used for grouping in the 
#' pseudobulk approach
#' @param celltype_ID the column name in the SCE object for the cell type 
#' variable 
#' @param y the column name in the SCE object for the return variable e.g. 
#' "diagnosis" - Case or disease
#' @param region The column name in the SCE object for the study region. If 
#' there are multiple regions in the study (for example two brain regions). 
#' Pseudobulk values can be derived separately. Default is "single_region" 
#' which will not split by region.
#' @param coef character specifying which level to investigate for the 
#' differential expression analysis e.g. in a case/control study use "case" if 
#' case is the identifier in the y column to get positive fold changes to 
#' relate to case samples. leave as default value for continuous y. 
#' @param control character specifying which control level for the differential 
#' expression analysis e.g. in a case/control/other study use "control" in the 
#' y column to compare against. NOTE only need to specify if more than two 
#' groups in y, leave as default value for two groups or continuous y.
#' @param pval_adjust_method the adjustment method for the p-value in the 
#' differential expression analysis. Default is benjamini hochberg "BH". See  
#' stats::p.adjust for available options
#' @param adj_pval the adjusted p-value cut-off for the differential expression 
#' analysis, 0-1 range
#' @param folder the folder where the graphs from the differential expression 
#' analysis are saved. Default will create a folder in the current working 
#' directory "sc_cell_type_de_graphs". False will skip plotting.
#' @param rmv_zero_count_genes whether genes with no count values in any cell 
#' should be removed. Default is TRUE
#' @param verbose logical indicating if extra information about the 
#' differential expression analysis should be printed
#' @return NULL
validate_input_parameters_de<-function(SCE, design, pseudobulk_ID, celltype_ID, 
                                       y, region, coef, control, 
                                       pval_adjust_method, adj_pval, folder, 
                                       rmv_zero_count_genes, verbose){
    if(class(SCE)[1]!="SingleCellExperiment")
        stop("Please input the data as a SingleCellExperiment object")
    if(!is.character(celltype_ID))
        stop(paste0("Please input a character for celltype_ID indicating the ",
                    "column holding the cell type information"))
    if(is.null(SCE[[celltype_ID]]))
        stop(paste0("The inputted celltype_ID: ",celltype_ID,
                    " is not present in the SCE object, perhaps check for ",
                    "spelling is correct"))
    #check design is a formula, formatting is correct and variables exist
    if(!inherits(design,"formula"))
        stop(paste0("Please input a formula for the design variable specifying",
                    " the comparison. See examples."))
    design_txt <- paste0(deparse(design,width.cutoff = 500),collapse=',')
    # check if formula contains `~` to ensure user inputting correct format
    if(!grepl( "~", design_txt, fixed = TRUE))
        stop(paste0("Please input a correctly formated formula for the design",
                    " variable containing a `~`. See examples."))
    design_txt <- gsub(".*~","",design_txt)
    design_txt <- gsub("^\\s+|\\s+$","",strsplit(design_txt, "[+]")[[1]])
    #check for duplicate entries
    if(length(design_txt)!=length(unique(design_txt)))
        stop("There are duplicate entries in the design formula")
    #check each variable to see if they are in SCE object
    for(i in design_txt){
        if(is.null(SCE[[i]]))
            stop(paste0("The inputted value: ",i,
                        " in the design formula is not present in the SCE ",
                        "object, perhaps check for spelling is correct"))
    }
    #only check y if user is setting it themselves
    if(!is.null(y)){
        if(!is.character(y))
            stop(paste0("Please input a character for y indicating the column ",
                        "holding the response variable information"))
        if(is.null(SCE[[y]]))
            stop(paste0("The inputted y value: ",y,
                        " is not present in the SCE object, perhaps check for ",
                        "spelling is correct"))
    }
    else{
        #if y not specified take last value in design matrix
        y <- design_txt[[length(design_txt)]]
        if(is.null(SCE[[y]]))
            stop(paste0("The inputted y value taken from your formula (the ",
                        "last variable): ",y,
                        " \nis not present in the SCE object, perhaps check ",
                        "the spelling is correct"))
    }
    #only check region if user is setting it themselves
    if(region!="single_region"){
        if(!is.character(region))
            stop(paste0("Please input a character for region indicating the ",
                        "column holding the variable information"))
        if(is.null(SCE[[region]]))
            stop(paste0("The inputted region value: ",region,
                        " is not present in the SCE object, perhaps check for ",
                        "spelling is correct"))
    }
    if(!is.null(coef)){
        if(!coef %in% unique(SCE[[y]]))
            stop(paste0("The inputted coef value: ",coef,
                        " is not present in y ", y,
                        "\nin the SCE object, perhaps check for spelling ",
                        "is correct"))
    }
    if(!is.null(control)){
        if(!control %in% unique(SCE[[y]]))
            stop(paste0("The inputted control value: ",control,
                        " is not present in y ", y,
                        "\nin the SCE object, perhaps check for spelling is",
                        " correct"))
    }
    if(!is.character(pseudobulk_ID))
        stop(paste0("Please input a character for pseudobulk_ID indicating the",
                    " column holding the patient identifier"))
    if(is.null(SCE[[pseudobulk_ID]]))
        stop(paste0("The inputted patient_ID value: ",pseudobulk_ID,
                    " is not present in the SCE object, ",
                    "\nperhaps check for spelling is correct"))
    if(!is.character(pval_adjust_method))
        stop(paste0("Please input a character for the pval_adjust_method ",
                    "indicating the method to be used"))
    if(class(adj_pval)[1]!="numeric")
        stop(paste0("Please input a 0-1 range number for the adj_pval ",
                    "indicating the cut off of significance for the ",
                    "differential expression analysis"))
    if(adj_pval<0|adj_pval>1)
        stop(paste0("Please input a 0-1 range number for the adj_pval ",
                    "indicating the cut off of significance for the ",
                    "differential expression analysis"))
    #if the user doesn't want to save plots from DE analysis they can pass 
    #in a false for the folder input parameter
    if(!isFALSE(folder))
        if(!is.character(folder))
            stop(paste0("For the folder input variable, please pass a director",
                        "y for the plots to be saved in.",
                        "\nGo with the default or pass in FALSE stop plotting"))
    if(!is.logical(verbose))
        stop("Please input TRUE/FALSE for verbose")
    if(!is.logical(rmv_zero_count_genes))
        stop("Please input TRUE/FALSE for rmv_zero_count_genes")
}
