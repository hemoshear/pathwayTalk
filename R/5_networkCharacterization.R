# extension of step 4: classification and network generation  --------------------------------------------------
# Taylor Derby Pourtaheri

#' @title Characterize reactome pathway networks.
#' @description Saves xlsx files containing the common names for Reactome pathways
#' # representing nodes in a given pathway network.
#' @importFrom magrittr %<>% %>%
#' @param networks_list A named list of igraph networks, where the names are the
#' subtypes represented in each network.
#' @param output_dir The directory to which the xlsx files will be saved.
#' @export
characterizeNetworks <- function(networks_list, output_dir = NULL){

    # using reactome.db:
    rpathways <- as.list(reactome.db::reactomePATHNAME2ID)
    rpathways_df <- tibble::enframe(rpathways) # %>% unnest
    colnames(rpathways_df) <- c('pathway_name', 'name')
    rpathways_df <- rpathways_df[grepl('HSA', rpathways_df$name),]
    rpathways_df$name %<>% gsub('\\-', '\\.', .)

    results <- c()

    for (i in 1:length(networks_list)){

        subtype <- names(networks_list[i])
        network <- networks_list[[i]]

        v_list <- igraph::as_data_frame(network, 'vertices')
        v_list %<>% dplyr::left_join(rpathways_df)

        # openxlsx::write.xlsx(v_list, glue::glue('{`output_dir`}/{`subtype`}_summary.xlsx'))

        results[[subtype]] <- v_list

    }

    return(results)

}









