# pathwayEnrichment.R

#' Do pathway enrichment with Fisher's exact test
#'
#' @param deg diffExpression --> entrez_to_hgnc
#' @param gene_alpha numeric vector indicating significance level for genes.
#' @return data frame of results for the Fisher's exact tests.
#'     `p` Unadjusted p-value
#'     `adj_p` FDR adjusted p-value
#'     `estimate` odds ratio point estimate
#'     `pathway` Paste of Reactome pathway identifier and description.
#' @export

.fisherPathwayEnrichment <- function(deg, gene_alpha) {
    pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    pathways <- pathways[grep('HSA', names(pathways))]
    sig <- deg %>% dplyr::filter(adj.P.Val < gene_alpha)
    nonsig <- deg %>% dplyr::filter(adj.P.Val >= gene_alpha)

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


#' @param importFrom magrittr %<>%
#' @param gene_alpha Significance level for genes. Genes with a differential expression
#'  p-value less than `gene_alpha` will be labeled as differentially expressed for the purpose
#'  doing a Fisher's exact test for pathway enrichment.
#' @param pathway_alpha Significance level for pathways. Pathways where the Fisher's exact
#'  test p-value is less than `pathway_alpha` will be labeled as enriched for the purposes
#'  of downstream analyses. Pathways with a p-value greater than or equal to `pathway_alpha`
#'  will be returned.
#' @return data frame of results for the Fisher's exact tests.
#'     `p` Unadjusted p-value
#'     `adj_p` FDR adjusted p-value
#'     `estimate` odds ratio point estimate
#'     `pathway` Paste of Reactome pathway identifier and description.
#'     `contrast` Character vector representing contrasts.
#' @export

fisherPathwayEnrichment <- function(deg, gene_alpha, pathway_alpha) {
    tests <- purrr::map(deg, ~ .fisherPathwayEnrichment(., gene_alpha=gene_alpha))
    names(tests) <- names(DEG)
    # keep enriched pathways for each cancer subtype.
    sig_pathways <- purrr::map(tests, ~ dplyr::filter(., adj_p < pathway_alpha))
    # create one data frame with enriched pathways across cancer subtypes
    for (i in 1:length(sig_pathways)) {
        sig_pathways[[i]]$contrast <- names(sig_pathways)[i]
    }
    sig_pathways %<>% dplyr::bind_rows()
    sig_pathways

}
