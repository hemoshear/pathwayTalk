#' Calculate matrix of discriminating scores for pathway crosstalk
#' 
#' @param enriched_pathways
#' @param eset 
#' 

pathwayCrosstalk <- function(enriched_pathways, eset, processes=1) {
    pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    pathways <- pathways[grep('HSA', names(pathways))] # get pathways for homo sapiens
    # Get data frame of expression values
    d <- as.data.frame(genes$datETcollapsed)
    d$entrez <- rownames(d)
    # I take the first Entrez ID if there are multiple to make this simpler
    d$canon_entrez <- strsplit(d$entrez, split = '///') %>%
        purrr::map_chr(~ .[1])
    
    pdata <- phenoData(eset)
    pdata@data$sample_id <- rownames(pdata@data)

    do_pathway_crosstalk <- function(d, enriched_pathways) {
        s <- list()
        for (pathway in unique(enriched_pathways$pathway)) { 
            # Filter data frame to genes in pathway
            to_keep <- which(d$canon_entrez %in% pathways[[pathway]])
            s[[pathway]] <- c('mean' = mean(as.numeric(d[to_keep,1])), 
                              'sd' = sd(as.numeric(d[to_keep,1])))
        }
        pathway_pairs <- gtools::combinations(n=length(unique(enriched_pathways$pathway)),
                                              r=2, v=unique(enriched_pathways$pathway))
        scores <- purrr::map2(as.character(pathway_pairs[,1]), 
                              as.character(pathway_pairs[,2]),
                              ~ (s[[.x]]['mean'] - s[[.y]]['mean']) /
                                  (s[[.x]]['sd'] + s[[.y]]['sd']))
        names(scores) <- apply(pathway_pairs, 1, function(x) (paste0(x[1], '___', x[2], collapse='___')))
        return(scores)
    }
    
    sample_data <- purrr::map(pdata@data$sample_id, ~ d[,c(., 'entrez', 'canon_entrez')])
    res <- parallel::mclapply(sample_data, function(x) (do_pathway_crosstalk(x, enriched_pathways)),
                              mc.cores=processes)
    pathway_pairs <- names(res[[1]]) # get pathway pairs for first sample, which will be same as the rest
    res <- purrr::map(res, ~ unlist(.))
    m <- do.call('cbind', res)
    colnames(m) <- pdata@data$sample_id
    rownames(m) <- pathway_pairs
    return(m)
}
