#' Calculate matrix of discriminating scores for pathway crosstalk
#'
#' @param enriched_pathways Character vector of enriched pathways
#' @param genes Output of WGCNA::collapseRows
#' @param pathways List where names are pathway names and elements are sets of Entrez
#'     gene identifiers.
#' @param processes Number of processes (integer) that should be used
#' @return Matrix of dimensions $m x n$ where $m$ is the number of pathway pairs and $n$ is the number
#'     of samples in the data. Each element of this matrix is the discriminating score for a particular sample,
#'     which is defined as:
#'     $\frac{mean(g_x) - mean(g_y)}{sd(g_x) - sd(g_y)}$
#'
#' @export

pathwayCrosstalk <- function(enriched_pathways, genes, pathways, processes=1) {

    # Get data frame of expression values
    d <- as.data.frame(genes$datETcollapsed)
    d$entrez <- rownames(d)
    # I take the first Entrez ID if there are multiple to make this simpler
    d$canon_entrez <- strsplit(d$entrez, split = '///') %>%
        purrr::map_chr(~ .[1])

    do_pathway_crosstalk <- function(d, enriched_pathways) {
        s <- list()
        for (pathway in enriched_pathways) {
            # Filter data frame to genes in pathway
            to_keep <- which(d$canon_entrez %in% pathways[[pathway]])
            s[[pathway]] <- c('mean' = mean(as.numeric(d[to_keep,1])),
                              'sd' = sd(as.numeric(d[to_keep,1])))
        }
        # Enumerate pairs of pathways
        pathway_pairs <- gtools::combinations(n=length(enriched_pathways),
                                              r=2, v=unique(enriched_pathways))
        # Calculate "discriminating scores" for all pairs
        scores <- purrr::map2(as.character(pathway_pairs[,1]),
                              as.character(pathway_pairs[,2]),
                              ~ (s[[.x]]['mean'] - s[[.y]]['mean']) /
                                  (s[[.x]]['sd'] + s[[.y]]['sd']))
        names(scores) <- apply(pathway_pairs, 1, function(x) (paste0(x[1], '___', x[2], collapse='___')))
        return(scores)
    }

    sample_data <- purrr::map(colnames(d), ~ d[,c(., 'entrez', 'canon_entrez')])
    # Calculate discriminating scores in parallel
    res <- parallel::mclapply(sample_data, function(x) (do_pathway_crosstalk(x, enriched_pathways)),
                              mc.cores=processes)
    pathway_pairs <- names(res[[1]]) # get pathway pairs for first sample, which will be same as the rest
    # Change structure of discriminating scores so we have
    # matrix with rows as pathway pairs and columns as df
    res <- purrr::map(res, ~ unlist(.))
    m <- do.call('cbind', res)
    colnames(m) <- colnames(d)
    rownames(m) <- pathway_pairs
    return(m)
}
