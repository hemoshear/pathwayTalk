#' @title fisherPathwayEnrichment
#' @param importFrom magrittr %<>%
#' @param DEGs A dataframe of differential gene expression results in with
#' associated p-values in a column 'pvalue'. Extracted from the output of diffExpression().
#' @param pathways A named list in which each element is a pathway containing a
#' character vector of the corresponding gene IDs.
#' @param gene_alpha Significance level for differentially expressed genes. Genes with a
#'  p-value less than `gene_alpha` will be labeled as differentially expressed for the purpose
#'  doing a Fisher's exact test for pathway enrichment.
#' @return A data frame of results for the Fisher's exact tests.
#'     \itemize{
#'     \item p - Unadjusted p-value.
#'     \item adj_p - FDR adjusted p-value.
#'     \item estimate - Odds ratio point estimate.
#'     \item pathway - Pathway ID.
#'     }
#' @export
fisherPathwayEnrichment <- function(DEGs, gene_alpha, pathways) {

    # if DEGs is a dataframe, convert to a list:
    if (class(DEGs) != 'list'){

        DEGs <- list(DEGs)
    }

    # Remove any duplicate gene ids from the pathways list
    pathways <- purrr::map(pathways, ~ unique(.))

    # define function to do fisher test
    .fisherTest <- function(...){
        result <- fisher.test(...)
        output <- data.frame(p = result$p.value,
                             estimate = result$estimate)
        return(output)
    }

    # nested lapply
    fisher_enrichment <- function(deg, pathways){

        tests <- lapply(deg, function(x)

            lapply(pathways, function(y)

                .fisherTest(factor(x$pvalue < gene_alpha,
                               levels = c('FALSE', 'TRUE')),
                        factor(x$canon_entrez %in% y,
                               levels = c('FALSE', 'TRUE')),
                        alternative = 'greater')))
        tests <- lapply(tests, function(x) dplyr::bind_rows(x, .id = 'pathway'))
        tests <- lapply(tests, function(x) dplyr::mutate(x, adj_p = p.adjust(p, method='fdr')))

        return(tests)

    }

    output <- fisher_enrichment(DEGs, pathways)

    return(output)

}


