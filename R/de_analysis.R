#' Differential Expression Analysis using edgeR LRT on pseudobulk data
#'
#' @param pb_dat A list containing
#' \itemize{
#'   \item \code{sumDat}: matrix of the summed pseudobulk count values
#'   \item \code{annot_pb}: dataframe of the annotation data from the SCE 
#'   rolled up based on the pseudobulk aggregation.
#' }
#' @param formula the validated design formula of class type `formula`. 
#' Equation used to fit the model- data for the generalised linear model 
#' e.g. ~ sex + pmi + disease.
#' @param y_name the column name in the SCE object for the return variable 
#' e.g. "diagnosis" - Case or disease. y can be discrete (logisitc regression) 
#' or continuous (linear regression)
#' @param y_contin is the variable being modelled continuous e.g. if 
#' case/control then TRUE if level of Tau in AD study then FALSE
#' @param coef character specifying which level to investigate for the 
#' differential expression analysis e.g. in a case/control study use "case" if 
#' case is the identifier in the y column to get positive fold changes to 
#' relate to case samples. Leave as default value for continuous y.
#' @param control character specifying which control level for the differential 
#' expression analysis e.g. in a case/control/other study use "control" in the 
#' y column to compare against. NOTE only need to specify if more than two 
#' groups in y, leave as default value for two groups or continuous y. 
#' @param pval_adjust_method the adjustment method for the p-value in the 
#' differential expression analysis. Default is benjamini hochberg "BH". See  
#' stats::p.adjust for available options
#' @param adj_pval the adjusted p-value cut-off for the differential expression 
#' analysis, 0-1 range
#' @param verbose logical indicating if extra information about the 
#' differential expression analysis should be printed
#' @return A list containing differential expression data (a dataframe) for 
#' each cell type. The dataframe contains log fold change (logFC), log counts 
#' per million (logCPM), log ratio (LR), p-value (PValue), adjusted p-value 
#' (adj_pval)
de_analysis <- function(pb_dat,formula,y_name,y_contin,coef, control,
                            pval_adjust_method, adj_pval,verbose){
    #Use Likelihood ratio test from edgeR, found to have best perf in a recent
    #benchmark: https://www.biorxiv.org/content/10.1101/2021.03.12.435024v1.full
    #Other option for edgeR is quasi-likelihood F-tests
    #Run each cell type through
    DE_ct_rtrn <- vector(mode="list",length=length(names(pb_dat)))
    names(DE_ct_rtrn) <- names(pb_dat)
    for(ct_i in names(pb_dat)){
        if(verbose)
            message("Analysing ",ct_i)
        targets <- pb_dat[[ct_i]]$annot_pb
        if(!y_contin){#adjust levels based on user input if not contin
            #chec if factor inputted, if not convert
            if(!is.factor(targets[[y_name]])){
                targets[[y_name]] <- as.factor(targets[[y_name]])
            }
            #reorder levels so coef is last so won't be intercept and chosen
            old_levels <- levels(targets[[y_name]])
            if(length(old_levels)>2){
                if(!is.null(coef)){#if user inputted value for case samples
                    if(is.null(control)){
                        message(paste0("WARNING: No control passed so coef may",
                                       " not be compared against desired ",
                                       "variable"))
                        new_levels <- 
                            c(old_levels[which(coef!=old_levels)],coef)
                    }
                    else{
                        new_levels <- 
                            c(control,
                              old_levels[which(coef!=old_levels & 
                                                   control!=old_levels)],coef)
                    }
                }
                else{
                    if(is.null(control)){
                        message(paste0("WARNING: No control passed so may not ",
                                       "be compared against desired variable"))
                        new_levels <- old_levels #do nothing no coef or control
                    }
                    else{
                        new_levels <- 
                            c(control,old_levels[which(control!=old_levels)])
                    }  
                }
            }
            else{
                if(!is.null(coef)){#if user inputted value for case samples
                    new_levels <- c(old_levels[which(coef!=old_levels)],coef)
                }
                else{
                    new_levels <- old_levels #do nothing no coef
                }
            }
            #assign the new levels
            targets[[y_name]] <- factor(targets[[y_name]],levels=new_levels)
        }
        # create design
        design <- model.matrix(formula, data = targets)   
        #slight diff approach for contin or cat
        if(y_contin){
            y <- edgeR::DGEList(counts = pb_dat[[ct_i]]$sumDat)
        }
        else{#categorical i.e. case/control
            y <- edgeR::DGEList(counts = pb_dat[[ct_i]]$sumDat,
                                    group = targets[[y_name]])
        }
        #Calculate normalization factors to scale the raw library sizes.
        y <- edgeR::calcNormFactors(y,method = 'TMM')
        #Maximizes the negative binomial likelihood to give the estimate of the 
        #common,trended and tagwise dispersions across all tags.
        y <- edgeR::estimateDisp(y, design)
        #LRT
        fit <- edgeR::glmFit(y, design = design)
        #slight diff approach for contin or cat
        if(y_contin){#need to specify coefficient target as not done in design
            test <- edgeR::glmLRT(fit,coef=y_name)
        }
        else{#categorical i.e. case/control, no need
            if(!is.null(coef)){#if user inputted value for case samples
                test <- edgeR::glmLRT(fit,coef=length(new_levels))
            }
            else{#they didn't direction may not be as they expect
                message(paste0("WARNING: No coefficient passed so direction ma",
                               "y not relate to desired variable"))
                test <- edgeR::glmLRT(fit)
            }
        }
        pvals <- test$table
        #add adj p-values
        pvals$adj_pval <- stats::p.adjust(pvals$PValue, 
                                          method = pval_adjust_method)
        pvals$name <- rownames(pvals)
        rownames(pvals) <- NULL
        if(verbose){
            numDEGs <- nrow(pvals[pvals$adj_pval<adj_pval,])
            message(numDEGs," DEGs found")
        }
        DE_ct_rtrn[[ct_i]] <- pvals
    }
    return(DE_ct_rtrn)
}
