#' @title Generate a network containing the top-performing pathway pairs.
#' @description Generates a glmnet binomial classification model for a given matrix of
#' pathway-pair discriminating scores and outputs a network object containing the
#' top-performing pathway pairs in the model.
#' @importFrom magrittr %>%
#' @param crosstalk_matrix  A matrix of discriminating scores in which each column
#' is a pair of enriched pathways and each row is a sample.
#' The output of pathwayCrosstalk().
#' @param groups A dataframe containing the mappings between sample identifiers
#' ('sample_id', a factor with the reference condition as the first level) and
#' associated treatment conditions ('group'). The sample identifiers must be in
#' the same order as the columns of the count_matrix.
#' @param alpha The alpha value applied to the glmnet classification model.
#' Use alpha = 1 for lasso regression and alpha = 0 for ridge regression.
#' @param lambda The lambda value applied to the glmnet classification model.
#' @param output_network Logical. If TRUE, the function outputs the network object.
#' @param output_model Logical. If TRUE, the function outputs the model.
#' @return Depending on the values of the `output_network` and `model_evaluation`
#' arguments, either a list containing the network object or a list containing
#' the model evaluation metrics.
#' @export
crosstalkNetwork <- function(crosstalk_matrix, groups, alpha = 1, lambda = 'best',
                            output_network = TRUE, output_model = FALSE){

    matrix_sample_ids <- rownames(crosstalk_matrix)
    ordered_group_vector <- groups[match(matrix_sample_ids, groups$sample_id), 'group'] %>%
        as.character

    crosstalk_matrix <- as.matrix(crosstalk_matrix)
    sample_phenotype <- as.matrix(ordered_group_vector)

    output_list <- list()

    if (lambda == 'best'){

        lasso_model <- glmnet::cv.glmnet(x = crosstalk_matrix, y = sample_phenotype, family = 'binomial',
                                alpha = alpha)

        lambda <- lasso_model$lambda.min

    }


    # fit glmnet model
    model <- glmnet::glmnet(x = crosstalk_matrix, y = sample_phenotype, family = 'binomial',
                            alpha = alpha, lambda = lambda)
    # sort(abs(model$beta), decreasing = TRUE) %>% View

    # generate network:
    if (output_network){

        # identify top pathway pairs (non-zero coefficients)
        coef_matrix <- as.matrix(abs(model$beta))
        keeps <- which(coef_matrix != 0) + 1 # offset by 1 due to intercept

        # ridge regression: only return top 15 pathway pairs
        if (alpha != 1){

            keeps <- keeps[1:15]
        }

        values <- coef(model)[keeps]
        pathways_pairs <- as.character(rownames(coef(model))[keeps])


        # generate dataframe of top-performing pathway pairs for 'edges' in network
        coefs <- data.frame(pathways_pairs, values)

        coefs$pathway1 <- purrr::map_chr(coefs$pathways_pairs,
                                     ~unlist(strsplit(x = ., split = '___'))[1])
        coefs$pathway2 <- purrr::map_chr(coefs$pathways_pairs,
                                         ~unlist(strsplit(x = ., split = '___'))[2])

        edges <- as.matrix(coefs[c('pathway1', 'pathway2')])

        # generate network
        network <- igraph::graph_from_edgelist(edges, directed = FALSE)

        # generate dataframe
        df <- igraph::as_data_frame(network) %>% setNames(c('V1', 'V2'))

        output_list[['network']] <- network
        output_list[['dataframe']] <- df

    }

    if (output_model){

        output_list[['model']] <- model
        output_list[['alpha']] <- alpha
        output_list[['lambda']] <- lambda

    }

    return(output_list)

}




