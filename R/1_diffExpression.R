# https://www.rpubs.com/jlubieni/13450

# purpose: a function to conduct differential expression analysis
    # for a given dataset - beginning with ExpressionSet class

# this function works with an ExpressionSet with a lockedEnvironment storage mode

# pre-processing completed prior to using function
    # (normalization, removing control sequences, log2 transformation, collapse to genes tc.)


# write function to do DEA ------------------------------------------------

#' @title Conduct differential gene expression analysis for a given gene expression dataset.
#' @description Calculates log2 fold-changes and associated p-values for a
#' given matrix of gene expression data.
#' @importFrom magrittr %>%
#' @param data An expression matrix of intensity values for gene probes,
#' where the row names are the probe identifiers. Pre-processing should be completed,
#' including normalization, removing control sequences, log2 transformation, collapse to genes.
#' @param gene_ids A vector of strings containing gene IDs corresponding to each row of the
#' expression matrix.
#' @param design_mat A matrix specifying the experimental design,
#' generated with model.matrix().
#' #' @param contrast_mat A matrix specifying the preferred contrasts, generated
#' with makeContrasts().
#' @export
diffExpression <- function(data,
                           design_mat,
                           contrast_mat) {

    # fit linear model using provided design and contrast matrices
    initial_fit <- limma::lmFit(data, design_mat)
    temp_fit <- limma::contrasts.fit(initial_fit, contrast_mat)
    fit <- limma::eBayes(temp_fit, robust = TRUE, trend = TRUE) # robust, trend

    contrast_names <- colnames(contrast_mat)

    results <-c()

    for(i in 1:length(contrast_names)){
        cont <- contrast_names[i]
        result <- limma::topTable(fit, coef = i, number = Inf)
        result$canon_entrez <- rownames(result) # added
        results[[i]] <- result
        # original: results[[i]] <- topTable(fit, coef = i, number = Inf)
    }

    names(results) <- contrast_names

    return(results)

}

