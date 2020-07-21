# https://www.rpubs.com/jlubieni/13450

# purpose: a function to conduct differential expression analysis
    # for a given dataset - beginning with ExpressionSet class

# this function works with an ExpressionSet with a lockedEnvironment storage mode

# pre-processing completed prior to using function (normalization, removing control sequences, etc.)


# write function to do DEA ------------------------------------------------

#' @export
diffExpression <- function(data,
                           collapse_method,
                           gene_ids,
                           design_mat,
                           contrast_mat) {

    # begin with log2 transformation of intensity values
    data <- log2(data)

    # collapse probe rows to genes
    probes <- data

    genes <<- WGCNA::collapseRows(datET = probes,
                                 rowGroup = gene_ids,
                                 rowID = rownames(probes),
                                 method = collapse_method)

    data <- genes$datETcollapsed

    # fit linear model using provided design and contrast matrices
    initial_fit <- limma::lmFit(data, design_mat)
    temp_fit <- limma::contrasts.fit(initial_fit, contrast_mat)
    fit <- limma::eBayes(temp_fit)

    contrast_names <- colnames(contrast_mat)

    # output
    # results <- decideTests(fit)

    results <-c()

    for(i in 1:length(contrast_names)){

        cont <- contrast_names[i]

        results[[i]] <- limma::topTable(fit, coef = i, number = Inf)

    }

    names(results) <- contrast_names

    results

}

