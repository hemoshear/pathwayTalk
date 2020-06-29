#' Calculate matrix of discriminatory scores for pathway crosstalk
#' 
#' @param enriched_pathways
#' @param eset 
#' 

pathwayCrosstalk <- function(enriched_pathways, eset) {
    pathways <- getReactomePathways()
    # Convert Entrez to Hugo in data frame where we have probes collapsed to gene-level
    d <- entrez_to_hgnc(list(as.data.frame(genes$datETcollapsed)))
    d <- d[[1]] 
    # Get phenoData slot of the original ExpressionSet so we can map
    # sample identifiers to treatments
    pdata <- phenoData(eset)
    pdata@data$sample_id <- rownames(pdata@data)# opinionated dplyr...
    pdata@data %<>% dplyr::filter(title != 'GFP') # ** remove controls for now, come back after dicussion
    treatment_map <- pdata@data$title
    names(treatment_map) <- pdata@data$sample_id

    ds <- list()
    print(names(enriched_pathways))
    for (i in 1:3) { # hack to only do subset for now
        # Enriched pathways for the subtype that 
        print(treatment_map[i])
        # length(treatment_map)
        
        this_pathway_set <- enriched_pathways[[treatment_map[i]]]
        s <- list()
        for (pathway in this_pathway_set$pathway) {
            this_d <- dplyr::filter(d, hgnc_symbol %in% pathways[[pathway]])
            this_d <- this_d[c(names(treatment_map)[i])][[1]] %>% as.numeric
            s[[pathway]] <- c('mean' = mean(this_d), 
                        'sd' = sd(this_d))
            
        }
        this_matrix <- matrix(NA, nrow=length(names(s)),
                              ncol=length(names(s)))
        rownames(this_matrix) <- colnames(this_matrix) <- names(s)
        # Create matrix with discriminatory scores for the pathway pairs
        for (n in 1:nrow(this_matrix)) {
            for (m in 1:ncol(this_matrix)) {
                pathway_n <- rownames(this_matrix)[n]
                pathway_m <- colnames(this_matrix)[m]
                this_matrix[n,m] <- (s[[pathway_n]]['mean'] - s[[pathway_m]]['mean']) / 
                    (s[[pathway_n]]['sd'] + s[[pathway_m]]['sd'])
            }
        }
        ds[[names(treatment_map)[i]]] <- this_matrix
    }
    ds
}
