#' @title networkWrapper
#' @description Wrapper function for the pathwayTalk workflow.
#' @importFrom magrittr %>%
#' @param expression_matrix An matrix of gene expression data where the
#' row names are the gene probe identifiers and column names are the sample identifiers.
#' For RNAseq data, a matrix of un-normalized integer counts.
#' For microarray data, a matrix of intensity values for gene probes with
#' pre-processing completed, including normalization, removing control sequences,
#' log2 transformation, and collapse to genes.
#' @param groups A dataframe containing the sample identifiers ('sample_id', a factor with
#' the reference condition as the first level) and associated treatment conditions ('group').
#' The sample identifiers must be in the same order as the columns of the count_matrix.
#' @param DEPs A dataframe of pathways with associated p-values in a column 'p'.
#' The output of pathwayEnrichment().
#' @param pathways A list in which each element is a named list of genes
#' corresponding to a particular pathway.
#' @param pathway_alpha Significance level for enriched pathways. Pathways with a Fisher's exact
#'  test p-value less than `pathway_alpha` will be labeled as enriched for the purposes
#'  of downstream analyses.
#' @param ... Additional arguments to crosstalkNetwork() function.
#' @return Returns a named list.
#'      \itemize{
#'      \item full_network - An igraph object representing the full
#'       network of top-performing pathway pairs.
#'      \item full_network_results - A dataframe describing the full network.
#'      \item pruned_network - An igraph object representing the pruned
#'       network of top-performing pathway pairs.
#'      \item pruned_network_results - A dataframe describing the pruned network.
#'      }
#' @export
#'

networkWrapper <- function(expression_matrix, groups, DEPs, pathways, pathway_alpha, lambda, ...){

    # remove any duplicate gene ids from the pathways list
    pathways <- purrr::map(pathways, ~ unique(.))

    # step 3 - generate matrix of discriminating scores
    crosstalk_matrix <- pathwayCrosstalk(expression_matrix = expression_matrix,
                                         DEPs = DEPs,
                                         groups = groups,
                                         pathways = pathways,
                                         pathway_alpha = pathway_alpha)

    # step 4 - perform lasso classification and generate pathway-pair network
        # from non-zero coefficients
    network_results <- crosstalkNetwork(crosstalk_matrix = crosstalk_matrix,
                                        groups = groups,
                                        lambda = lambda,
                                        output_network = TRUE,
                                        output_model = FALSE,
                                        ...)
    network <- network_results[['network']]


    # step 6 - perform simulated crosstalk inhibition to identify
        # significant crosstalks between pairs of pathways
    crosstalk_inhibition_results <- crosstalkInhibition(network = network)

    # save results
    # output_list <- c()
    # output_list[['crosstalk_matrix']] <- crosstalk_matrix
    # output_list[['crosstalk_network']] <- network_results
    # output_list[['crosstalk_inhibition']] <- crosstalk_inhibition_results


    # return output
    return(crosstalk_inhibition_results)

}


