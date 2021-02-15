
# step 4: classification and network generation  --------------------------------------------------

# phenotype splitting function --------------------------------------------

#' @title Outputs list of matrices, each containing samples from a single disease phenotype
#' and the reference condition.
#' @description Splits a given matrix of discriminating scores by phenotype.
#' The output can be used to generate a binary classification model.
#' @importFrom magrittr %>%
#' @param input_matrix A matrix of discriminating scores for pathway pairs,
#' where each row represents a sample and each columns is a pair of pathway.
#' @param phenotypes A vector of strings containing the phenotype associated with each sample.
#' Must be the same length as the number of rows in the feature_matrix input.
#' @param reference_condition String. The phenotype identifier for the reference condition in *phenotypes*.
#' For example, 'GFP'.
#' @export
byPhenotype <- function(input_matrix, phenotypes, reference_condition){

    subtypes <- groups %>% setdiff(reference_condition)

    phenotypes %<>% factor %>% relevel(ref = reference_condition)

    # pairs_index: list of indices for each phenotype + reference cond.
    pairs_index <- c()

    for (subtype in subtypes){
        this_subtype <- which(phenotypes == subtype)
        ref_con <- which(phenotypes == reference_condition)
        pairs_index[[subtype]] <- c(ref_con, this_subtype)
    }

    # pairs: list of sub-matrices containing samples for each phenotype + reference cond.
    pairs <- c()

    for (subtype in names(pairs_index)){
        # print(subtype)
        index <- pairs_index[[subtype]]
        pairs[[subtype]] <- input_matrix[index,]
    }

    # phenotypes_list: list of phenotypes corresponding to samples in each sub-matrix

    phenotypes_list <- c()

    for (subtype in subtypes){
        subtype_indices <- pairs_index[[subtype]]
        phenotypes_list[[subtype]] <- phenotypes[subtype_indices]
    }


    output_list <- c()

    output_list[['matrices']] <- pairs
    output_list[['phenotypes']] <- phenotypes_list

    return(output_list)
}


# select features with greatest variance -------------------------------------------------------

#' @title Allows for selection of features with greatest variation.
#' @description Selects the features of a given input matrix with the greatest variance.
#' The n_top_features argument allows for the specification of the number of features returned.
#' Returns a matrix containing only the top features of the original matrix.
#' @importFrom magrittr %>%
#' @param input_matrix A matrix of discriminating scores for pathway pairs,
#' where each row represents a sample and each columns is a pair of pathways.
#' @param n_top_features Integer. The number of top features returned.
#' @export
selectTopFeatures <- function(input_matrix, n_top_features = 20000){

    variances <- input_matrix %>% summarise_if(is.numeric, var)
    ordered_matrix <- input_matrix[, order(variances, decreasing = TRUE)]
    short_matrix <- ordered_matrix[, 1:n_top_features]

    return(short_matrix)
}


# binary classification and graph generation ------------------------------

#' @title Outputs a graph object of the top-performing pathway pairs for a given subtype.
#' @description Generates a glmnet binomial classification model for a given matrix of
#' pathway-pair discriminating scores.
#' Outputs a graph object containing the top-performing pathway pairs in the model.
#' @importFrom magrittr %>%
#' @param feature_matrix A matrix of discriminating scores for pathway pairs,
#' where each row represents a sample and each columns is a pair of pathways.
#' #' The matrix should include samples from the reference condition and one disease phenotype.
#' @param groups A dataframe containing the sample identifiers and associated treatment conditions.
#' Must be the same length as the number of rows in the feature_matrix input.
#' @param alpha The alpha value applied to the glmnet classification model.
#' Use alpha = 1 for lasso regression and alpha = 0 for ridge regression.
#' @param lambda The lambda value applied to the glmnet classification model.
#' @param output_graph Logical. If TRUE, the function outputs the graph object.
#' @param model_evaluation Logical. If TRUE, the function outputs the model evaluation metrics.
#' @export
subtypeNetwork <- function(feature_matrix, groups, alpha = 1, lambda = 'best',
                            output_graph = TRUE, model_evaluation = FALSE){

    matrix_sample_ids <- rownames(feature_matrix)
    ordered_group_vector <- groups[match(matrix_sample_ids, groups$sample_id), 'group'] %>%
        as.character

    feature_matrix <- as.matrix(feature_matrix)

    sample_phenotype <- ordered_group_vector
    sample_phenotype <- as.matrix(ordered_group_vector)

    output_list <- list()

    if (lambda == 'best'){

        full_data <- as.data.frame(feature_matrix)

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

    # generate graph:
    if (output_graph){

        # fit glmnet model
        model <- glmnet::glmnet(x = feature_matrix, y = sample_phenotype, family = 'binomial',
                        alpha = alpha, lambda = lambda)
        sort(abs(model$beta), decreasing = TRUE)

        # identify top pathway pairs (non-zero coefficients)
        coef_matrix <- as.matrix(abs(model$beta))
        keeps <- which(coef_matrix != 0) + 1 # these indices offset by 1

        # ridge regression: only return top 15 pathway pairs
        if (alpha != 1){

            keeps <- keeps[1:15]
        }

        values <- coef(model)[keeps]
        pathways_pairs <- as.character(rownames(coef(model))[keeps])


        # generate dataframe of top-performing pathway pairs for 'edges' in network
        coefs <- data.frame(pathways_pairs, values)

        coefs$pathway1 <- stringr::str_extract(string = coefs$pathways_pairs, pattern = '^R.HSA.\\d+')
        coefs$pathway2 <- stringr::str_extract(string = coefs$pathways_pairs, pattern = 'R.HSA.\\d+$')

        edges <- as.matrix(coefs[c('pathway1', 'pathway2')])

        # generate graph
        graph <- igraph::graph_from_edgelist(edges, directed = FALSE)

        # generate dataframe
        df <- igraph::as_data_frame(graph)

        output_list[['alpha']] <- alpha
        output_list[['lambda']] <- lambda
        output_list[['graph']] <- graph
        output_list[['dataframe']] <- df

        return(output_list)

    }

    if (model_evaluation){

        # split into test and train subsets
        set.seed(8)
        index <- sample.int(nrow(feature_matrix), 0.8*nrow(feature_matrix))
        train <- feature_matrix[index,]
        test <- feature_matrix[-index,]

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




