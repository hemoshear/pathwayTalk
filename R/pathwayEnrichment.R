# pathwayEnrichment.R

#' Get Reactome pathways and corresponding Hugo Gene Symbols
#' for Homo sapiens
#' 
#' @param importFrom magrittr %<>%
#' @param species Character vector to indicate which species should be used
#' @return list T
#' @export 

getReactomePathways <- function(species = 'Homo sapiens') {
    reactome <- GSA::GSA.read.gmt('ReactomePathways.gmt')
    pathways <- read.delim('https://reactome.org/download/current/ReactomePathways.txt',
                           header=FALSE)
    colnames(pathways) <- c('reaction_pathway', 'reactome_description', 'species')
    pathways %<>% dplyr::filter(species %in% species)
    index <- purrr::map(pathways$reaction_pathway, 
                        ~ which(reactome$geneset.descriptions == .))
    names(index) <- pathways$reaction_pathway
    index %<>% Filter(function(x) (length(x) > 0), .)
    for (i in 1:length(index)) {
        index[[i]] <- reactome$genesets[[index[[i]]]]
    }
    index
}
