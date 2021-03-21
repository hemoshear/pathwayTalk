library(caret)
library(future)
library(furrr)

# a series of test/training partitions are created using createDataPartition
# using the groups as the outcome (y) ensures even splitting

cvWrapper <- function(expression_matrix,
                      groups,
                      platform,
                      gene_alpha,
                      pathways,
                      pathway_alpha,
                      lambda,
                      sampling_method = 'partition',
                      times = 10,
                      p = 0.6){

    # output list
    output_list <- c()

    expression_matrix_t <- t(expression_matrix)

    # split data into resamples
    if (sampling_method == 'partition'){

        # generate indices to subset data
        indices <- caret::createDataPartition(y  = groups$group,
                                              times = times,
                                              p = p)
        # subset data
        exprs_list <- purrr::map(indices, ~ t(expression_matrix_t[., ])) # weird format to preserve column names
        groups_list <- purrr::map(indices, ~ groups[., ])
        resample <- purrr::map(groups_list, ~.$sample_id)
        resample_summary <- purrr::map(groups_list, ~table(.$group))
    }


    # else {
    #     return(paste0("Sampling method should be a string ",
    #                   "with the value 'partition', 'down', or 'up'"))
    # }

    # save sampling results
    output_list[['resample']] <- resample
    output_list[['resample_summary']] <- resample_summary

    # # # #### TEMP ####
    # resample_exprs <- exprs_list[[1]]
    # resample_groups <- groups_list[[1]]

    internalWrapper <- function(resample_exprs, resample_groups, platform,
                                gene_alpha, pathways, pathway_alpha){

        internal_return_list <- c()

        # step 1: do differential expression analysis
        DEG <- diffExpression(expression_matrix = resample_exprs,
                                groups = resample_groups,
                                platform = platform)
        internal_return_list[['DEG']] <- DEG
        # extract DEGs from diffExpression output
        DEG_dfs <- purrr::map(DEG, ~.$DEGs)

        # step 2: do Fisher test for pathway enrichment analysis
        enriched <- fisherPathwayEnrichment(DEGs = DEG_dfs,
                                            gene_alpha=gene_alpha,
                                            pathways = reactome_pathways)
        # remove infinite values, if present
        enriched_short <- purrr::map(enriched, ~filter(., !is.infinite(estimate)))
        internal_return_list[['enriched']] <- enriched_short

        # step 3 - generate matrix of discriminating scores
        # remove any duplicate gene ids from the pathways list
        pathways <- purrr::map(pathways, ~ unique(.))
        crosstalk_matrices <- purrr::map2(DEG, enriched_short,
                                         ~ pathwayCrosstalk(expression_matrix = .x$data,
                                                            groups = .x$group,
                                                            DEPs = .y,
                                                            pathways = pathways,
                                                            pathway_alpha = pathway_alpha))
        internal_return_list[['crosstalk_matrices']] <- crosstalk_matrices


        # # ### TEMP
        # sub_groups <- DEG$LUMB$groups
        # crosstalk_matrix <- crosstalk_matrices$LUMB

        # step 4 - classification
        model_results <- purrr::map2(DEG, crosstalk_matrices,
                                     ~ crosstalkNetwork(crosstalk_matrix = .y,
                                                        groups = .x$groups,
                                                        lambda = lambda,
                                                        alpha = 1,
                                                        output_network = TRUE,
                                                        output_model = FALSE))
        # internal_return_list[['model_results']] <- purrr::map(model_results, ~.$model)
        internal_return_list[['network_results']] <- purrr::map(model_results, ~.$dataframe)


        return(internal_return_list)


    }

    # future::plan(multisession, workers = 4)
    internal_wrapper_results <- purrr::map2(exprs_list, groups_list, ~internalWrapper(.x, .y,
                                            platform, gene_alpha, pathways, pathway_alpha))

    # get vertex list for each network
    extractVertices <- function(resample_results){

        net_results <- resample_results$network_results

        net_df <- purrr::map(net_results, ~as.data.frame(.)) %>%
            bind_rows(.id = 'contrast')
    }
    network_results <- purrr::map(internal_wrapper_results,
                                           ~extractVertices(.)) %>%
        bind_rows(.id = 'resample')
    network_results_list <- split(network_results, network_results$contrast)
    keeps <- purrr::map(network_results_list,
                                       ~ duplicated(.[,c('V1', 'V2')]))
    network_results_list <- purrr::map2(network_results_list, keeps,
                                        ~ .x[.y,c('V1', 'V2')])
    network_results_list <- purrr::map(network_results_list,
                                        ~ unique(.))


    # generate full network
    generateFullNetwork <- function(network_df){
        unique_vertices <- unique(network_df)
        network <- igraph::graph_from_data_frame(unique_vertices,
                                                 directed = FALSE)
        return(network)

    }

    full_networks <- purrr::map(network_results_list, ~generateFullNetwork(.))

    # step 6: crosstalk inhibition
    crosstalk_inhibition_results <- purrr::map(full_networks,
                                               ~ crosstalkInhibition(network = .))

    # save components
    # output_list[['internal_wrapper_results']] <- internal_wrapper_results
    # output_list[['full_networks']] <- full_networks
    output_list[['crosstalk_inhibition_results']] <- crosstalk_inhibition_results

    return(output_list)

}


# sampling methods:

# else if(sampling_method == 'down'){
#     full_list <- purrr::map(1:times,
#                             ~ caret::downSample(x = expression_matrix_t,
#                                                 y = groups$group,
#                                                 list = TRUE))
#     exprs_list <- purrr::map(full_list, ~t(.$x))
#     groups_list <- purrr::map(full_list, ~.$y)
#     sample_summary <- purrr::map(groups_list, ~table(.))
# }
#
# else if(sampling_method == 'up'){
#     full_list <- purrr::map(1:times,
#                             ~ caret::upSample(x = expression_matrix_t,
#                                                 y = groups$group,
#                                                 list = TRUE))
#     exprs_list <- purrr::map(full_list, ~t(.$x))
#     groups_list <- purrr::map(full_list, ~.$y)
#     sample_summary <- purrr::map(groups_list, ~table(.))
#
# }

# further test/train splitting:
#
# # train/test split
# generateTestTrain <- function(sub_groups, crosstalk_matrix){
#
#     testTrainOutput <- c()
#
#     if (sampling_method == 'partition'){
#         train_indices <- caret::createDataPartition(y  = sub_groups$group,
#                                                     times = 1,
#                                                     p = 0.6)
#         # training set
#         xt_train <- purrr::map(train_indices, ~ crosstalk_matrix[., ])
#         groups_train <- purrr::map(train_indices, ~ sub_groups[., ])
#         sample_summary_train <- purrr::map(groups_train, ~table(.$group))
#
#         # testing set
#         xt_test <- purrr::map(train_indices, ~ crosstalk_matrix[-., ])
#         groups_test <- purrr::map(train_indices, ~ sub_groups[-., ])
#         sample_summary_test <- purrr::map(groups_test, ~table(.$group))
#     }
#
#     testTrainOutput[['xt_train']] <- xt_train[[1]]
#     testTrainOutput[['groups_train']] <- groups_train[[1]]
#     testTrainOutput[['sample_summary_train']] <- sample_summary_train[[1]]
#     testTrainOutput[['xt_test']] <- xt_test[[1]]
#     testTrainOutput[['groups_test']] <- groups_test[[1]]
#     testTrainOutput[['sample_summary_test']] <- sample_summary_test[[1]]
#
#     return(testTrainOutput)
# }
#
# # # test/train splits for each contrast
# # a <- generateTestTrain(sub_groups, crosstalk_matrix)
# # iterate over contrasts
# splits <- purrr::map2(DEG, crosstalk_matrices, ~generateTestTrain(.x$groups, .y))


# model evaluation

# # #### TEMP
# # model <- model_results$LUMB$model
# # test_features <- splits$LUMB$xt_test
# # test_groups <- splits$LUMB$groups_test
#
# # testing
# modelPredict <- function(model, test_features, test_groups){
#
#     model_evaluation <- c()
#
#     preds <- predict(object = model,
#                      family = 'binomial',
#                      newx = as.matrix(test_features),
#                      type = 'response')
#     table <- table(preds, test_groups$group)
#     conf_mat <- table(preds>0.5, test_groups$group)
#
#     eval <- glmnet::assess.glmnet(model, newx = as.matrix(test_features),
#                                   newy = test_groups$group)
#
#     model_evaluation[['predictions']] <- table
#     model_evaluation[['confusion_matrix']] <- conf_mat
#     model_evaluation[['evaluation']] <- eval
#
#
#     return(model_evaluation)
#
# }
#
# evaluation_results <- purrr::map2(model_results, splits,
#                                   ~modelPredict(model = .x$model,
#                                                 test_features = .y$xt_test,
#                                                 test_groups = .y$groups_test))
# internal_return_list[['evaluation_results']] <- evaluation_results


# selecting best network by mse:
# # #### TEMP
# # resample_groups <- groups_list[[1]]
# # resample_crosstalk <- internal_wrapper_results$crosstalk_matrices[[1]]
# # resample_results <- internal_wrapper_results[[1]]
#
# extractMSE <- function(resample_results){
#
#     mod_results <- resample_results$evaluation_results
#
#     eval_results <- purrr::map(mod_results, ~.$evaluation)
#     eval_df <- purrr::map(eval_results, ~as.data.frame(.)) %>%
#         bind_rows(.id = 'contrast')
# }
#
# model_evaluation_results <- purrr::map(internal_wrapper_results,
#                                        ~extractMSE(.)) %>%
#     bind_rows(.id = 'resample')
#
# resample_df <- model_evaluation_results %>% dplyr::group_by(contrast) %>%
#     dplyr::slice(which.min(mse))
#
#
# # select best model for each contrast
# best_resample <- resample_df$resample
# names(best_resample) <- resample_df$contrast
#
# selectBestNetwork <- function(condition){
#
#     resample <- best_resample[[condition]]
#     best_resample <- internal_wrapper_results[[resample]]
#     best_network <- best_resample$network_results[[condition]]
#
#     return(best_network)
#
# }
#
# conditions <- names(internal_wrapper_results[[1]]$DEG)
# networks <- purrr::map(conditions, ~selectBestNetwork(.))
# names(networks) <- conditions
#
# best_resample <- resample_df$resample
# names(best_resample) <- resample_df$contrast
#
# selectBestNetwork <- function(condition){
#
#     resample <- best_resample[[condition]]
#     best_resample <- internal_wrapper_results[[resample]]
#     best_network <- best_resample$network_results[[condition]]
#
#     return(best_network)
#
# }
#
#
