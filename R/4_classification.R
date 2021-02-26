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
#' @param model_evaluation Logical. If TRUE, the function outputs the model evaluation metrics.
#' @return Depending on the values of the `output_network` and `model_evaluation`
#' arguments, either a list containing the network object or a list containing
#' the model evaluation metrics.
#' @export
crosstalkNetwork <- function(crosstalk_matrix, groups, alpha = 1, lambda = 'best',
                            output_network = TRUE, model_evaluation = FALSE){

    matrix_sample_ids <- rownames(crosstalk_matrix)
    ordered_group_vector <- groups[match(matrix_sample_ids, groups$sample_id), 'group'] %>%
        as.character

    crosstalk_matrix <- as.matrix(crosstalk_matrix)

    sample_phenotype <- ordered_group_vector
    sample_phenotype <- as.matrix(ordered_group_vector)

    output_list <- list()

    if (lambda == 'best'){

        full_data <- as.data.frame(crosstalk_matrix)

        full_data$phenotype <- factor(sample_phenotype)

        train_control <- caret::trainControl(method = "cv",
                                      number = 10,
                                      savePredictions = TRUE,
                                      classProbs = TRUE)

        lasso_model <- caret::train(phenotype ~ .,
                             data = full_data,
                             method = 'glmnet',
                             family = 'binomial',
                             tuneGrid = expand.grid(.alpha = alpha,
                                                    .lambda=seq(0.05, 1, 0.05)),
                             trControl = train_control)

        lambda <- lasso_model$bestTune$lambda

    }

    # generate network:
    if (output_network){

        # fit glmnet model
        model <- glmnet::glmnet(x = crosstalk_matrix, y = sample_phenotype, family = 'binomial',
                        alpha = alpha, lambda = lambda)
        # sort(abs(model$beta), decreasing = TRUE)

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

        output_list[['alpha']] <- alpha
        output_list[['lambda']] <- lambda
        output_list[['network']] <- network
        output_list[['dataframe']] <- df

        return(output_list)

    }

    if (model_evaluation){

        # split into test and train subsets
        set.seed(8)
        index <- sample.int(nrow(crosstalk_matrix), 0.8*nrow(crosstalk_matrix))
        train <- crosstalk_matrix[index,]
        test <- crosstalk_matrix[-index,]

        train_y <- as.matrix(sample_phenotype[index])
        train_x <- as.matrix(train)

        test_y <- as.matrix(sample_phenotype[-index])
        test_x <- as.matrix(test)

        # fit glmnet model
        model <- glmnet(x = train_x, y = train_y, family = 'binomial',
                        alpha = alpha, lambda = lambda)

        test_y %<>% factor

        # evaluating with assess.glmnet:
        a <- glmnet::assess.glmnet(model, newx = test_x, newy = test_y, family = 'binomial')
        c <- glmnet::confusion.glmnet(model, newx = test_x, newy = test_y, family = 'binomial')
        r <- glmnet::roc.glmnet(model, newx = test_x, newy = test_y, family = 'binomial')

        output_list[['model_evaluation_metrics']] <- a
        output_list[['modelconfusion_matrix']] <- c
        output_list[['model_ROC']] <- r
        output_list[['alpha']] <- alpha
        output_list[['lambda']] <- lambda

        return(output_list)

    }

}




