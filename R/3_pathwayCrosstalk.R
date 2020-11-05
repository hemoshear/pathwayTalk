#' Calculate matrix of discriminating scores for pathway crosstalk
#'
#' @param enriched_pathways Output of fisherPathwayEnrichment
#' @param treatment_map
#' @export

pathwayCrosstalk <- function(enriched_pathways, treatment_map, processes=4) {
    pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    pathways <- pathways[grep('HSA', names(pathways))] # get pathways for homo sapiens
    # Get data frame of expression values
    do_pathway_crosstalk <- function(expr, enriched_pathways) {
        score <- list()
        for (pathway in unique(enriched_pathways$pathway)) {
            # Filter data frame to genes in pathway
            to_keep <- which(expr$entrez %in% pathways[[pathway]])
            this_d <- expr[to_keep,]
            # Get gene expression values in the pathway
            this_d <- this_d[[1]] %>% as.numeric
            score[[pathway]] <- c('mean' = mean(this_d),
                              'sd' = sd(this_d))
        }
        # Enumerate pathway pairs
        combs <- gtools::combinations(n=length(score),
                                      r=2, v=names(score))
        combs <- lapply(1:nrow(combs), function(x) (c(combs[x,1], combs[x,2])))
        # Calculate "discriminating score" for each pair of enriched pathways
        crosstalk <- sapply(combs, function(x) (
            (score[[x[1]]]['mean'] - score[[x[2]]]['mean']) /
                (score[[x[1]]]['sd'] + score[[x[2]]]['sd'])
            )
        )
        # Create paste of in format PATHWAY1___PATHWAY2
        comb_id <- purrr::map_chr(combs, ~ paste0(.[1], '___', .[2]))
        names(crosstalk) <- comb_id
        # Return named vector
        return(crosstalk)

    }
    # Split sample data
    sample_data <- purrr::map(names(treatment_map), ~ expr[,c(., 'entrez')])
    # rm(expr)
    # Calculate discriminating scores in parallel
    res <- parallel::mclapply(sample_data, function(x) (do_pathway_crosstalk(x, enriched_pathways)),
                              mc.cores=processes)
    res %<>% do.call('cbind', .)
    colnames(res) <- names(treatment_map)
    # names(res) <- names(treatment_map)
    return(data.frame(t(res)))
}

