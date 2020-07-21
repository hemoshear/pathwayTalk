# pathwayEnrichment.R

#' Pathway enrichment using Fisher's exact test
#'
#' @importFrom magrittr %>%
#' @importFrom stats fisher.test p.adjust sd
#' @param deg Output of diffExpression function.
#' @param alpha numeric vector indicating significance level for the differential expression analysis. Genes in the differential
#' expression analysis with an adjusted p-value less than alpha are considered to be differentially expressed for the purpose
#' of doing Fisher's exact test for every pathway.
#' @param pathways Named list where names are pathway identifiers and elements are sets of Entrez gene identifiers (character)
#' @return data frame of results for the Fisher's exact tests.
#'     `p` Unadjusted p-value.
#'     `adj_p` FDR adjusted p-value.
#'     `estimate` odds ratio point estimate.
#'     `pathway` pathway identifier.
#' @export

fisherPathwayEnrichment <- function(deg, alpha, pathways) {

    sig <- deg %>% dplyr::filter(adj.P.Val < alpha)
    nonsig <- deg %>% dplyr::filter(adj.P.Val >= alpha)

    # Enumerate contingency table, this can be cleaned up a lot.
    sig_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(sig, canon_entrez %in% .) %>% nrow)

    non_sig_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(nonsig, canon_entrez %in% .) %>% nrow)

    sig_not_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(sig, !(canon_entrez %in% .)) %>% nrow)

    non_sig_not_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(nonsig, !(canon_entrez %in% .)) %>% nrow)

    # Do Fisher's exact test for each  pathway.
    contingency <- fisher <- list()
    for (i in 1:length(pathways)){
        contingency[[i]] <- matrix(c(sig_in_pathway[[i]],
                                     non_sig_in_pathway[[i]],
                                     sig_not_in_pathway[[i]],
                                     non_sig_not_in_pathway[[i]]),
                                   byrow = TRUE, nrow = 2, ncol=2)
        res <- fisher.test(contingency[[i]], alternative = 'greater')
        # Change structure of fisher.test output to data frame
        fisher[[i]] <- data.frame(p = res$p.value, adj_p = p.adjust(res$p.value, method='fdr'),
                         estimate = res$estimate, pathway=names(pathways)[i])
    }
    dplyr::bind_rows(fisher) %>%
        dplyr::arrange(p)
}


