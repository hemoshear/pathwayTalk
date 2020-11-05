#' Calculate matrix of discriminating scores for pathway crosstalk
#'
#' @param enriched_pathways Output of fisherPathwayEnrichment
#' @param treatment_map
#' @export pathwayCrosstalkParallel

pathwayCrosstalk <- function(enriched_pathways, treatment_map, processes=4) {
    pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    pathways <- pathways[grep('HSA', names(pathways))] # get pathways for homo sapiens
    # Get data frame of expression values
    d <- as.data.frame(genes$datETcollapsed)
    d$entrez <- rownames(d)
    # I take the first Entrez ID if there are multiple to make this simpler
    d$canon_entrez <- strsplit(d$entrez, split = '///') %>%
        purrr::map_chr(~ .[1])

    do_pathway_crosstalk <- function(d, enriched_pathways) {
        score <- list()
        for (pathway in unique(enriched_pathways$pathway)) {
            # Filter data frame to genes in pathway
            to_keep <- which(d$canon_entrez %in% pathways[[pathway]])
            this_d <- d[to_keep,]
            # Get gene expression values in the pathway
            this_d <- this_d[[1]] %>% as.numeric
            score[[pathway]] <- c('mean' = mean(this_d),
                              'sd' = sd(this_d))
        }
        combs <- gtools::combinations(n=length(score),
                                      r=2, v=names(score))
        combs <- lapply(1:nrow(combs), function(x) (c(combs[x,1], combs[x,2])))
        crosstalk <- sapply(combs, function(x) (
            (score[[x[1]]]['mean'] - score[[x[2]]]['mean']) -
                (score[[x[1]]]['sd'] + score[[x[2]]]['sd'])
            )
        )
        comb_id <- purrr::map_chr(combs, ~ paste0(.[1], '___', .[2]))
        names(crosstalk) <- comb_id
        return(crosstalk)
#        crosstalk <- data.frame(pathway_pair = comb_id,
#                                score = crosstalk)
#        return(crosstalk)
    }

    sample_data <- purrr::map(names(treatment_map), ~ d[,c(., 'entrez', 'canon_entrez')])
    res <- parallel::mclapply(sample_data, function(x) (do_pathway_crosstalk(x, enriched_pathways)),
                              mc.cores=processes)
    res %<>% do.call('cbind', .)
    colnames(res) <- names(treatment_map)


    # names(res) <- names(treatment_map)
    return(data.frame(t(res)))
}

