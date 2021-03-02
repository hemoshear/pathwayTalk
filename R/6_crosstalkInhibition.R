#' @title Conduct simulated network crosstalk inhibition and prune network.
#' @description Given a network, calculate the network efficiency. Iterate over
#' network edges and conduct simulated network crosstalk inhibition, pruning the network
#' of edges that do not contribute to network efficiency.
#' @importFrom magrittr %<>%
#' @param network An igraph object in which the nodes represent enriched pathways and
#' the edges represent significant pathway crosstalks. The output of crosstalkNetwork().
#' @return Returns a named list.
#'      \itemize{
#'      \item full_network - An igraph object representing the full
#'       network of top-performing pathway pairs.
#'      \item full_network_results - A dataframe describing the full network.
#'      \item pruned_network - An igraph object representing the pruned
#'       network of top-performing pathway pairs.
#'      \item pruned_network_results - A dataframe describing the pruned network.
#'      }
#' @export

networkEfficiency <- function(network){

    # shortest path lengths between all nodes
    paths <- data.frame(igraph::shortest.paths(network))
    paths$pathway <- rownames(paths)
    paths_df <- setNames(tidyr::pivot_longer(paths, -pathway),
                         c('V1', 'V2', 'path_length'))

    # select only unique combinations, remove self-references
    paths_df <- paths_df[paths_df$V1 != paths_df$V2,]

    pairs <- c()

    paths_df$V1 %<>% as.character
    paths_df$V2 %<>% as.character

    for (i in 1:nrow(paths_df)){
        these_paths <- as.character(c(paths_df[i, 'V1'], paths_df[i, 'V2']))
        sorted <- gtools::mixedsort(these_paths)
        pasted <- paste(sorted[1], sorted[2], sep = '|')
        pairs[i] <- pasted
    }

    paths_df$pairs <- pairs

    # remove duplicates
    paths_df <- paths_df[!(duplicated(paths_df$pairs)),]

    # remove infinite values (unconnected nodes)
    paths_df <- paths_df[!(is.infinite(paths_df$path_length)),]

    # number of nodes in the network
    N <- as.numeric(length(igraph::V(network)))

    # calculate network efficiency
    NE <- (sum(1/paths_df$path_length)) / (N*(N-1))
    return(NE)

}



crosstalkInhibition <- function(network){

    results <- c()

    # initial network efficiency of the network
    NE <- networkEfficiency(network)

    full_edge_df <- data.frame(igraph::as_edgelist(network)) %>%
        setNames(c('pathway1', 'pathway2'))
    full_edge_df$edge <- paste0(full_edge_df$pathway1, '|', full_edge_df$pathway2)
    full_edge_df$NE <- NE

    # NE and PCI
    nNE_results <- c()
    PCI_results <- c()
    significant_crosstalks <- c()

    # iterate over edges of the network
    for (i in 1:nrow(full_edge_df)){

        # delete one edge
        new_network <- igraph::delete.edges(network, i)

        # calculate new network efficiency and percent change and store
        nNE <- networkEfficiency(new_network)
        nNE_results[i] <- nNE

        percent_change <- ((nNE - NE) / NE) * 100
        PCI_results[i] <- percent_change

        # # store results if network efficiency is decreased
        # # by the removal of the edge
        # if (nNE < NE) {
        #     significant_crosstalks[i] <- full_edge_df[i, 'edge']
        # }
    }

    # add network efficiency and PCI information to edge dataframe
    full_edge_df$nNE <- nNE_results
    full_edge_df$PCI <- PCI_results

    # select significant crosstalks (nNE < NE)
    pruned_edge_df <- full_edge_df[full_edge_df$PCI < 0, ]

    # generate new network from significant crosstalks and associated edge df
    pruned_network <- igraph::graph_from_edgelist(as.matrix(pruned_edge_df[,c(1,2)]),
                                                  directed = FALSE)


    results[['full_network']] <- network
    results[['full_network_results']] <- full_edge_df
    results[['pruned_network']] <- pruned_network
    results[['pruned_network_results']] <- pruned_edge_df

    return(results)
}





