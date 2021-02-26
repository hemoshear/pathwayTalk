#' @title pathwayCrosstalk
#' @param expression_matrix An matrix of gene expression data where the
#' row names are gene probe identifiers and column names are sample identifiers.
#' For RNAseq data, a matrix of un-normalized integer counts. For microarray data,
#' a matrix of intensity values for gene probes with pre-processing completed,
#' including log2 transformation, normalization, removal of control sequences.
#' @param DEPs A data frame of results for the Fisher's exact tests with a column 'p'.
#' The output of fisherPathwayEnrichment().
#' @param groups A dataframe containing the mappings between sample identifiers
#' ('sample_id', a factor with the reference condition as the first level) and
#' associated treatment conditions ('group'). The sample identifiers must be in
#' the same order as the columns of the count_matrix.
#' @param pathways A named list in which each element is a pathway containing a
#' character vector of the corresponding gene IDs.
#' @param pathway_alpha Significance level for pathways. Pathways where the Fisher's exact
#'  test p-value is less than `pathway_alpha` will be labeled as enriched for the purposes
#'  of downstream analyses.
#'  @return A matrix of discriminating scores in which each column is a pair of
#'  enriched pathways and each row is a sample.
#' @export

pathwayCrosstalk <- function(expression_matrix, DEPs,
                             groups, pathways, pathway_alpha) {

    # Remove any duplicate gene ids from the pathways list
    unique_pathways <- purrr::map(pathways, ~ unique(.))

    # Filter DEPs to only the significantly enriched pathways
    DEPs <- filter(DEPs, p < pathway_alpha)

    # Filter pathways list to relevant pathways and convert to list
    e_pathways <- unique_pathways[names(unique_pathways) %in% DEPs$pathway]
    e_pathways <- purrr::map(e_pathways, ~ .[. %in% rownames(expression_matrix)])
    pathways_list <- purrr::map(names(e_pathways), ~ expression_matrix[e_pathways[[.]], ] )
    names(pathways_list) <- names(e_pathways)

    # Calculate sample means and standard deviations for genes in each pathway
    pathway_means <- purrr::map(pathways_list, colMeans)
    names(pathway_means) <- names(pathways_list)
    pathway_stdevs <- purrr::map(pathways_list, ~ apply(., 2, sd))
    names(pathway_stdevs) <- names(pathways_list)
    combs <- gtools::combinations(n = length(pathways_list), 2, names(pathways_list))

    # Calculate discriminating scores
    ds <- apply(combs, 1, function(x) {
        (pathway_means[[x[1]]] - pathway_means[[x[2]]]) /
            (pathway_stdevs[[x[1]]] + pathway_stdevs[[x[2]]])
        })

    # Create paste of in format PATHWAY1___PATHWAY2 and set row names
    comb_ids <- purrr::map2_chr(combs[,1], combs[,2],  ~ paste0(.x, '___', .y))
    colnames(ds) <- comb_ids

    # Return dataframe of scores with pathway pairs as features and sample ids as observations
    return(data.frame(ds))

}
