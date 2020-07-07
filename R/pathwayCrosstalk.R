#' Calculate matrix of discriminatory scores for pathway crosstalk
#' 
#' @param enriched_pathways
#' @param eset 
#' 

pathwayCrosstalk <- function(enriched_pathways, eset) {
    pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    pathways <- pathways[grep('HSA', names(pathways))]    # Convert Entrez to Hugo in data frame where we have probes collapsed to gene-level
    d <- as.data.frame(genes$datETcollapsed)
    d$entrez <- rownames(d)
    d$canon_entrez <- strsplit(d$entrez, split = '///') %>%
        purrr::map_chr(~ .[1])
    # Get phenoData slot of the original ExpressionSet so we can map
    # sample identifiers to treatments
    pdata <- phenoData(eset)
    pdata@data$sample_id <- rownames(pdata@data)# opinionated dplyr...
    treatment_map <- pdata@data$title
    names(treatment_map) <- pdata@data$sample_id
    ds <- list()
    for (i in 1:length(treatment_map)) { 
        s <- list()
        for (pathway in unique(enriched_pathways$pathway)) {
            # Filter data frame with hugo gene symbols that co
            this_d <- dplyr::filter(d, canon_entrez %in% pathways[[pathway]])
            # Get gene expression values in the pathway
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
                if (m > n) {
                    this_matrix[n,m] <- NA
                }
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


#' Calculate matrix of discriminatory scores for pathway crosstalk
#' 
#' @param enriched_pathways
#' @param eset 
#' 

pathwayCrosstalkParallel <- function(enriched_pathways, eset) {
    pathways <- as.list(reactome.db::reactomePATHID2EXTID)
    pathways <- pathways[grep('HSA', names(pathways))]    # Convert Entrez to Hugo in data frame where we have probes collapsed to gene-level
    
    d <- as.data.frame(genes$datETcollapsed)
    d$entrez <- rownames(d)
    d$canon_entrez <- strsplit(d$entrez, split = '///') %>%
        purrr::map_chr(~ .[1])
    # Get phenoData slot of the original ExpressionSet so we can map
    # sample identifiers to treatments
    pdata <- phenoData(eset)
    pdata@data$sample_id <- rownames(pdata@data)# opinionated dplyr...
    treatment_map <- pdata@data$title
    names(treatment_map) <- pdata@data$sample_id

    do_pathway_crosstalk <- function(d, enriched_pathways) {
        s <- list()
        for (pathway in unique(enriched_pathways$pathway)[1:300]) { # 1:400 is a hack for debugging right now
            # Filter data frame 
            to_keep <- which(d$canon_entrez %in% pathways[[pathway]])
            this_d <- d[to_keep,]
            # Get gene expression values in the pathway
            this_d <- this_d[[1]] %>% as.numeric
            s[[pathway]] <- c('mean' = mean(this_d), 
                              'sd' = sd(this_d))
        }
        this_matrix <- matrix(NA, nrow=length(names(s)),
                              ncol=length(names(s)))
        rownames(this_matrix) <- colnames(this_matrix) <- names(s)
        # Create matrix with discrimination scores for the pathway pairs
        for (n in 1:nrow(this_matrix)) {
            for (m in 1:ncol(this_matrix)) {
                if (m > n) {
                    this_matrix[n,m] <- NA
                }
                pathway_n <- rownames(this_matrix)[n]
                pathway_m <- colnames(this_matrix)[m]
                this_matrix[n,m] <- (s[[pathway_n]]['mean'] - s[[pathway_m]]['mean']) / 
                    (s[[pathway_n]]['sd'] + s[[pathway_m]]['sd'])
            }
        }
        return(this_matrix)
    }
    
    sample_data <- purrr::map(names(treatment_map), ~ d[,c(., 'entrez', 'canon_entrez')])
    # do_pathway_crosstalk(sample_data[[1]], enriched_pathways)
    res <- parallel::mclapply(sample_data, function(x) (do_pathway_crosstalk(x, enriched_pathways)),
                              mc.cores=40)
    names(res) <- names(treatment_map)
    
    return(res)
}
