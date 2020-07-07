# pathwayEnrichment.R

#' Do pathway enrichment with Fisher's exact test
#' 
#' @param deg diffExpression --> entrez_to_hgnc
#' @param alpha numeric vector indicating significance level
#' @return data frame of results for the Fisher's exact tests.
#'     `p` Unadjusted p-value
#'     `conf_int_lb` Lower bound of odds ratio estimate
#'     `conf_int_ub` Upper bound of odds ratio estimate
#'     `estimate` odds ratio point estimate
#'     `pathway` Paste of Reactome pathway identifier and description.
#' @export

fisherPathwayEnrichment <- function(deg, alpha) {
    pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    pathways <- pathways[grep('HSA', names(pathways))]
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
    
    # Do Fisher's exact test for each reactome pathway.
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


