# pathwayEnrichment.R

#' Get Reactome pathways and corresponding Hugo Gene Symbols
#' for Homo sapiens
#' 
#' @param importFrom magrittr %<>%
#' @param species Character vector to indicate which species should be used. Can be more than 1.
#' @return list where each element is a vector of Hugo gene symbols that correspond to a Reactome pathway.
#' @export 

getReactomePathways <- function(species = 'Homo sapiens') {
    # reactome_infile <- system.file('extdata', 'ReactomePathways.gmt', package = 'pathwayTalk')
    # This file contains the Hugo gene symbols that correspond to each Reactome Pathway
    reactome <- GSA::GSA.read.gmt('inst/extdata/ReactomePathways.txt')
    # This file contains the reactome pathway identifiers, descriptions, and species.
    pathways <- read.delim('https://reactome.org/download/current/ReactomePathways.txt',
                           header=FALSE)
    colnames(pathways) <- c('reaction_pathway', 'reactome_description', 'species')
    # Filter pathways by species
    pathways %<>% dplyr::filter(species %in% species)
    # Build up a list of pathways
    index <- purrr::map(pathways$reaction_pathway, 
                        ~ which(reactome$geneset.descriptions == .))
    names(index) <- paste0(pathways$reaction_pathway, '___', pathways$reactome_description)
    index %<>% Filter(function(x) (length(x) > 0), .)
    for (i in 1:length(index)) {
        index[[i]] <- reactome$genesets[[index[[i]]]]
    }
    index
}


#' Link Entrez identifiers to Hugo Gene Symbols for the object
#' produced by diffExpression
#' 
#' @importFrom magrittr %<>%
#' @param deg Output of diffExpression, which is a list of topTables
#'   for each contrast. 
#' @return Original input is returned with two new columns: `entrezgene_id` and `hgnc_symbol`.
#'   
entrez_to_hgnc <- function(deg, mart='hsapiens_gene_ensembl') {
    ensembl <- biomaRt::useMart('ensembl')
    datasets <- biomaRt::listDatasets(ensembl)
    hensembl <- biomaRt::useDataset(mart=ensembl, dataset=mart)
    mapping <- biomaRt::getBM(attributes=c('hgnc_symbol', 'entrezgene_id'),
                              mart=hensembl)
    mapping %<>% dplyr::filter(!is.na(entrezgene_id), !is.na(hgnc_symbol),
                               hgnc_symbol != "")
    mapping %<>% dplyr::arrange(entrezgene_id)
    ids <- purrr::map(deg, ~ rownames(.) %>%
                          strsplit(split='///') %>%
                          lapply(., function(x) (x[1])) %>%
                          unlist)
    mapping$entrezgene_id %<>% as.character()
    for (i in 1:length(ids)) {
        deg[[i]]$entrezgene_id <- ids[[i]]
        deg[[i]] %<>% dplyr::inner_join(mapping)
    }
    deg
}


#' Do pathway enrichment with Fisher's exact test
#' 
#' @param deg diffExpression --> entrez_to_hgnc
#' @param fdr Numeric vector to indicate FDR for calling significance
#' @return data frame of results for the Fisher's exact tests.
#'     `p` Unadjusted p-value
#'     `conf_int_lb` Lower bound of odds ratio estimate
#'     `conf_int_ub` Upper bound of odds ratio estimate
#'     `estimate` odds ratio point estimate
#'     `pathway` Paste of Reactome pathway identifier and description.
#' @export

fisherPathwayEnrichment <- function(deg, fdr) {
    pathways <- getReactomePathways()
    sig <- deg %>% dplyr::filter(adj.P.Val < fdr)
    nonsig <- deg %>% dplyr::filter(adj.P.Val >= fdr)
    
    # Enumerate contingency table, this can be cleaned up a lot. 
    sig_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(sig, hgnc_symbol %in% .) %>% nrow)
    
    non_sig_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(nonsig, hgnc_symbol %in% .) %>% nrow)
    
    sig_not_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(sig, !(hgnc_symbol %in% .)) %>% nrow)
    
    non_sig_not_in_pathway <- purrr::map(
        pathways, ~ dplyr::filter(nonsig, !(hgnc_symbol %in% .)) %>% nrow)
    
    # Do Fisher's exact test for each reactome pathway.
    contingency <- fisher <- list()
    for (i in 1:length(pathways)){
        contingency[[i]] <- matrix(c(sig_in_pathway[[i]],
                                     non_sig_in_pathway[[i]],
                                     sig_not_in_pathway[[i]],
                                     non_sig_not_in_pathway[[i]]),
                                   byrow = TRUE, nrow = 2, ncol=2)
        res <- fisher.test(contingency[[i]])
        # Change structure of fisher.test output to data frame
        fisher[[i]] <- data.frame(p = res$p.value, conf_int_lb = res$conf.int[1],
                         conf_int_ub = res$conf.int[2],
                         estimate = res$estimate, pathway=names(pathways)[i])
    }
    dplyr::bind_rows(fisher) %>%
        dplyr::arrange(p)
}


# Run diffExpression.R first
# a <- entrez_to_hgnc(b)
# tests <- purrr::map(a, ~ fisherEnrichment(., fdr=0.01))
# names(tests) <- names(a)
