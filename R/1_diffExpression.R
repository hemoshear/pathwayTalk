
library(tidyverse)
library(magrittr)
library(Biobase)
library(limma)
library(WGCNA)
library(gsap)
library(caret)


# https://www.rpubs.com/jlubieni/13450

# purpose: a function to conduct differential expression analysis
    # for a given dataset - beginning with ExpressionSet class

# this function works with an ExpressionSet with a lockedEnvironment storage mode

# pre-processing completed prior to using function (normalization, removing control sequences, etc.)


# write function to do DEA ------------------------------------------------

#' @title Conduct differential gene expression analysis for a given gene expression dataset.
#' @description Collapses probe-level expression to genes using collapseRows()
#' from WGCNA, and calculates log2 fold-changes and associated p-values for a
#' given matrix of gene expression data.
#' @importFrom magrittr %>%
#' @param data An expression matrix of intensity values for gene probes,
#' where the row names are the probe identifiers.
#' @param collapse_method Specifies the method by which the probes will be collapsed
#' to gene expression. The default value is 'MaxMean'. See the WGCNA collapseRows()
#' documentation for more options.
#' @param gene_ids A vector of strings containing gene IDs corresponding to each row of the
#' expression matrix.
#' @param design_mat A matrix specifying the experimental design,
#' generated with model.matrix().
#' #' @param contrast_mat A matrix specifying the preferred contrasts, generated
#' with makeContrasts().
#' @export
diffExpression <- function(data,
                           collapse_method = 'MaxMean',
                           gene_ids,
                           design_mat,
                           contrast_mat) {

    # begin with log2 transformation of intensity values
    data <- log2(data)

    # collapse probe rows to genes
    probes <- data

    genes <- WGCNA::collapseRows(datET = probes,
                                 rowGroup = gene_ids,
                                 rowID = rownames(probes),
                                 method = collapse_method)

    data <- genes$datETcollapsed

    # fit linear model using provided design and contrast matrices
    initial_fit <- lmFit(data, design_mat)
    temp_fit <- contrasts.fit(initial_fit, contrast_mat)
    fit <- eBayes(temp_fit)

    contrast_names <- colnames(contrast_mat)

    results <-c()

    for(i in 1:length(contrast_names)){
        cont <- contrast_names[i]
        result <- topTable(fit, coef = i, number = Inf)
        result$canon_entrez <- rownames(result) # added
        results[[i]] <- result
        # original: results[[i]] <- topTable(fit, coef = i, number = Inf)
    }

    names(results) <- contrast_names

    return(results)

}

