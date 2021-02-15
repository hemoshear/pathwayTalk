#' Calculate matrix of discriminating scores for pathway crosstalk
#'
#' @param DEPs The output of fisherPathwayEnrichment.
#' @param groups A dataframe containing the sample identifiers ('sample_id')
#' and associated treatment conditions ('group').
#' @param pathways A list in which each element is a named list of entrez gene IDs
#' corresponding to a particular pathway.
#' @export

pathwayCrosstalk <- function(DEPs, expression_matrix,
                             groups, pathways) {

    # ****
    # filter to remove duplicate genes within elements of pathway list
    # add a warning to flag duplicates
    #


    # pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    # pathways <- pathways[grep('HSA', names(pathways))] # get pathways for homo sapiens

    sample_ids <- groups$sample_id

    # shouldn't need the entrez column in expression_matrix

    e_pathways <- pathways[names(pathways) %in% DEPs$pathway]
    e_pathways <- purrr::map(e_pathways, ~ .[. %in% rownames(expression_matrix)])
    pathways_list <- purrr::map(names(e_pathways),~ expression_matrix[e_pathways[[.]], ] )
    names(pathways_list) <- names(e_pathways)

    pathway_means <- purrr::map(pathways_list, colMeans)
    names(pathway_means) <- names(pathways_list)
    pathway_stdevs <- purrr::map(pathways_list, ~ apply(., 2, sd))
    names(pathway_stdevs) <- names(pathways_list)
    combs <- gtools::combinations(n = length(pathways_list), 2, names(pathways_list))

    # combs <- lapply(1:nrow(combs), function(x) (c(combs[x,1], combs[x,2])))
    # (pathway_means[[combs[1,1]]] - pathway_means[[combs[1,2]]]) /
    #     (pathway_stdevs[[combs[1,1]]] + pathway_stdevs[[combs[1,2]]])

    ds <- apply(combs, 1, function(x) {
        (pathway_means[[x[1]]] - pathway_means[[x[2]]]) /
            (pathway_stdevs[[x[1]]] + pathway_stdevs[[x[2]]])
        })

    # combs  %<>% cbind(., t(ds)) %>% as.data.frame
    # combs  %<>% cbind(., as.data.frame(t(ds)))


    # # Get data frame of expression values
    # do_pathway_crosstalk <- function(expr, enriched_pathways) {
    #     score <- list()
    #     for (pathway in unique(enriched_pathways$pathway)) {
    #         # Filter data frame to genes in pathway
    #         to_keep <- which(expr$entrez %in% pathways[[pathway]])
    #         this_d <- expr[to_keep,]
    #         # Get gene expression values in the pathway
    #         this_d <- this_d[[1]] %>% as.numeric
    #         score[[pathway]] <- c('mean' = mean(this_d),
    #                           'sd' = sd(this_d))
    #     }
    #     # Enumerate pathway pairs
    #     combs <- gtools::combinations(n=length(score),
    #                                   r=2, v=names(score))
    #     combs <- lapply(1:nrow(combs), function(x) (c(combs[x,1], combs[x,2])))
    #     # Calculate "discriminating score" for each pair of enriched pathways
    #     crosstalk <- sapply(combs, function(x) (
    #         (score[[x[1]]]['mean'] - score[[x[2]]]['mean']) /
    #             (score[[x[1]]]['sd'] + score[[x[2]]]['sd'])
    #         )
    #     )

    # Format correction to prevent breaking:

    # Create paste of in format PATHWAY1___PATHWAY2 and set row names
    comb_ids <- purrr::map2_chr(combs[,1], combs[,2],  ~ paste0(.x, '___', .y))
    colnames(ds) <- comb_ids

    # Return dataframe of scores with pathway pairs as features and sample ids as observations
    return(data.frame(ds))

    # }

    # # Split sample data
    # treatment_map <- groups$group
    # names(treatment_map) <- groups$sample_id
    # sample_data <- purrr::map(names(treatment_map), ~ expr[,c(., 'entrez')])
    # # rm(expr)
    # # Calculate discriminating scores in parallel
    # res <- parallel::mclapply(sample_data, function(x) (do_pathway_crosstalk(x, enriched_pathways)),
    #                           mc.cores=processes)
    # res %<>% do.call('cbind', .)
    # colnames(res) <- names(treatment_map)
    # # names(res) <- names(treatment_map)
    # return(data.frame(t(res)))
}

