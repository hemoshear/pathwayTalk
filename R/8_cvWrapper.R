library(caret)
library(future)
library(furrr)

# a series of test/training partitions are created using createDataPartition
# using the groups as the outcome (y) ensures even splitting

cvWrapper <- function(sampling_method = 'partition',
                      times = 10,
                      p = 0.6,
                      expression_matrix,
                      groups,
                      platform,
                      gene_alpha,
                      pathways,
                      pathway_alpha){

    # output list
    output_list <- c()

    expression_matrix_t <- t(expression_matrix)

    if (sampling_method == 'partition'){
        train_indices <- caret::createDataPartition(y  = groups$group,
                                              times = times,
                                              p = p)
        # training set
        exprs_list_train <- purrr::map(train_indices, ~ t(expression_matrix_t[., ]))
        groups_list_train <- purrr::map(train_indices, ~ groups[., ])
        sample_summary_train <- purrr::map(groups_list_train, ~table(.$group))

        # testing set
        exprs_list_test <- purrr::map(train_indices, ~ t(expression_matrix_t[-., ]))
        groups_list_test <- purrr::map(train_indices, ~ groups[-., ])
        sample_summary_test <- purrr::map(groups_list_test, ~table(.$group))

    }


    else{
        return(paste0("Sampling method should be a string ",
                      "with the value 'partition', 'down', or 'up'"))
    }

    # save sampling results
    output_list[['training_sample']] <- sample_summary_train
    output_list[['testing_sample']] <- sample_summary_test

    internalWrapper <- function(exprs_sub, groups_sub, platform,
                                gene_alpha, pathways, pathway_alpha, model=TRUE){

        return_list <- c()

        # exprs_sub <- exprs_list_train[[1]]
        # groups_sub <- groups_list_train[[1]]

        # step 1: do differential expression analysis
        DEG <- diffExpression(expression_matrix = exprs_sub,
                                groups = groups_sub,
                                platform = platform)
        return_list[['DEG']] <- DEG

        # extract DEGs from diffExpression output
        DEG_dfs <- purrr::map(DEG, ~.$DEGs)

        # step 2: do Fisher test for pathway enrichment analysis
        enriched <- fisherPathwayEnrichment(DEGs = DEG_dfs,
                                            gene_alpha=gene_alpha,
                                            pathways = reactome_pathways)

        # remove infinite values, if present
        enriched_short <- purrr::map(enriched, ~filter(., !is.infinite(estimate)))

        # step 3 - generate matrix of discriminating scores

        # remove any duplicate gene ids from the pathways list
        pathways <- purrr::map(pathways, ~ unique(.))

        crosstalk_matrices <- purrr::map2(DEG, enriched_short,
                                         ~ pathwayCrosstalk(expression_matrix = .x$data,
                                                            groups = .x$groups,
                                                            DEPs = .y,
                                                            pathways = pathways,
                                                            pathway_alpha = pathway_alpha))
        return_list[['crosstalk_matrices']] <- crosstalk_matrices

        # step 4 - classification
        if (model == TRUE){
        model_results <- purrr::map2(crosstalk_matrices, DEG,
                                     ~ crosstalkNetwork(crosstalk_matrix = .x,
                                                        groups = .y$groups,
                                                        lambda = lambda,
                                                        alpha = 1,
                                                        output_network = TRUE,
                                                        output_model = TRUE))
        return_list[['model_results']] <- purrr::map(model_results, ~.$model)

        }

        return(return_list)
    }

    # generate training and testing sets
    train_models <- purrr::map2(exprs_list_train, groups_list_train,
                                ~internalWrapper(exprs_sub = .x,
                                                 groups_sub = .y,
                                                 platform = 'microarray',
                                                 gene_alpha = gene_alpha,
                                                 pathways = pathways,
                                                 pathway_alpha = pathway_alpha,
                                                 model = TRUE))


    # testing
    modelPredict <- function(model, test_data, test_groups){

        model_evaluation <- c()

        preds <- predict(model, as.matrix(test_data),
                         type = 'class') # finish this function

        conf_mat <- table(preds, test_groups$group)

        eval <- glmnet::assess.glmnet(model, newx = as.matrix(test_data),
                                  newy = test_groups$group)

        model_evaluation[['confusion_matrix']] <- conf_mat
        model_evaluation[['evaluation']] <- eval

    }

    # define conditions
    conditions <- levels(groups$group) %>% .[2:length(.)]

    a <- runEval(conditions, train_models){


    }




    }

    # select best model


}















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
