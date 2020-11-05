
# step 6: network efficiency and cross-talk inhibition  --------------------------------------------------


networkEfficiency <- function(network){

    # shortest path lengths between all nodes
    paths <- igraph::shortest.paths(network)
    paths_df <- setNames(melt(paths), c('V1', 'V2', 'path_length'))

    # select only unique combinations, remove self-references
    paths_df %<>% filter(V1 != V2)

    pairs <- c()

    paths_df$V1 %<>% as.character
    paths_df$V2 %<>% as.character

    for (i in 1:nrow(paths_df)){
        sorted <- mixedsort(c(paths_df[i, 'V1'], paths_df[i, 'V2']))
        pasted <- paste(sorted, sep = '|')
        pairs[i] <- pasted
    }

    paths_df$pairs <- pairs

    # remove infinite values (unconnected nodes)
    paths_df %<>% filter(!(is.infinite(path_length)))

    # number of nodes in the network
    N <- as.numeric(length(V(network)))

    # calculate network efficiency
    NE <- (sum(1/paths_df$path_length)) / (N*(N-1))
    return(NE)

}



crosstalkInhibition <- function(networks_list){

    results <- c()

    for(i in 1:length(networks_list)){

        network_name <- names(networks_list)[i]
        network <- networks_list[[network_name]]

        print(network_name)

        # initial network efficiency of the network
        NE <- networkEfficiency(network)

        edge_list <- igraph::as_edgelist(network)

        # list to store subtype-specific results
        significant_crosstalks <- c()

        # iterate over edges of the network
        for (i in 1:nrow(edge_list)){

            edge <- paste0(edge_list[i,1], '|', edge_list[i,2])

            # delete one edge
            new_network <- delete.edges(network, i)

            # calculate new network efficiency
            nNE <- networkEfficiency(new_network)

            # store results if network efficiency is decreased
            # by the removal of the edge
            if (nNE < NE) {

                percent_change <- ((nNE - NE) / NE) * 100
                significant_crosstalks[edge] <- percent_change
            }
        }

        results[[network_name]] <- significant_crosstalks
    }

        return(results)
}


# testing and characterization -----------------------------------------------------------------


characterizeResults <- function(sig_crosstalks_results, output_dir = NULL){

    results <- c()

    # using reactome.db:
    rpathways <- as.list(reactome.db::reactomePATHNAME2ID)
    rpathways_df <- tibble::enframe(rpathways) %>% unnest
    colnames(rpathways_df) <- c('pathway_name', 'name')
    rpathways_df <- rpathways_df[grepl('HSA', rpathways_df$name),]
    rpathways_df$name %<>% gsub('\\-', '\\.', .)

    for (i in 1:length(sig_crosstalks_results)){

        subtype <- names(sig_crosstalks_results)[i]
        pathways <- sig_crosstalks_results[[i]]

        PCI_results <- stack(sig_crosstalks_results[[i]])
        colnames(PCI_results) <- c('PCI', 'pathway_pair')

        subtype_results <- c()

        for (j in 1:length(pathways)){

            pathway_pair <- names(pathways)[j]
            pathway_pair_split <-  unlist(str_split(pathway_pair, '\\|'))

            pathway1 <- pathway_pair_split[1]
            pathway2 <- pathway_pair_split[2]

            pathway1_name <- rpathways_df$pathway_name[rpathways_df$name == pathway1 ]
            pathway2_name <- rpathways_df$pathway_name[rpathways_df$name == pathway2 ]

            pathway_names_pasted <- paste0(pathway1_name, ' | ', pathway2_name)

            subtype_results[[pathway_pair]] <- pathway_names_pasted

        }

        # openxlsx::write.xlsx(subtype_results, glue::glue('{`output_dir`}/{`subtype`}_summary.xlsx'))

        subtype_results_df <- data.frame('pathway_pair' = names(subtype_results),
                                         'pathway_names' = subtype_results)
        subtype_results_df %<>% left_join(PCI_results, by = 'pathway_pair')

        results[[subtype]] <- subtype_results_df
    }

    return(results)
}




