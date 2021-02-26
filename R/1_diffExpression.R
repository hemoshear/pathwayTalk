#' @title Conduct differential gene expression analysis for a given gene expression dataset.
#' @description Calculates log2 fold-changes and associated p-values for a
#' given matrix of gene expression data. diffExpression() will split the data
#' according to the levels of the 'group' column in the supplied `groups` object,
#' and generate appropriately cross-sectioned expression matrices, groups dataframes,
#' and design matrices, which are returned with the differential expression results.
#' These objects can be supplied to subsequent functions in the pathwayTalk pipeline.
#' @importFrom magrittr %>%
#' @param expression_matrix An matrix of gene expression data where the
#' row names are gene probe identifiers and column names are sample identifiers.
#' For RNAseq data, a matrix of un-normalized integer counts. For microarray data,
#' a matrix of intensity values for gene probes with pre-processing completed,
#' including log2 transformation, normalization, removal of control sequences.
#' @param groups A dataframe containing the mappings between sample identifiers
#' ('sample_id', a factor with the reference condition as the first level) and
#' associated treatment conditions ('group'). The sample identifiers must be in
#' the same order as the columns of the count_matrix.
#' @param platform A string specifying the data type. Either 'rnaseq' or 'microarray'.
#' @return A named list of contrasts, with each element containing the following objects:
#'      `data` The expression matrix cross-section relevant to the contrast.
#'      `groups` The groups dataframe cross-section relevant to the contrast.
#'      `design` The design matrix relevant to the contrast.
#'      `DEG` A dataframe containing the results of the differential expression analysis.
#' @export
diffExpression <- function(expression_matrix,
                           groups,
                           platform,
                           processes = 4) {

    # split data into list of dataframes by contrast
    ref_condition <- levels(groups$group)[1]
    conditions <- unique(groups$group) %>% setdiff(ref_condition)


    split_data <- function(condition, ref_condition, groups_df, counts){

        results <- c()

        groups_sub <- groups_df[groups_df$group %in% c(ref_condition, condition),]
        groups_sub$group %<>% factor(levels = c(ref_condition, condition))
        these_counts <- counts[,groups_sub$sample_id]

        # generate design matrix:
        this_design <- model.matrix(~ groups_sub$group)
        colnames(this_design) %<>% gsub('groups_sub\\$group', '', .)

        results <- list(these_counts, groups_sub, this_design)

    }

    data_list <- c()
    data_list <- purrr::map(conditions,
                            ~ split_data(., ref_condition, groups, expression_matrix))
    names(data_list) <- conditions

    do_rnaseq_DEA <- function(condition){

        # construct DESeqDataSet object

        dea_data <- data_list[[condition]][[1]]
        dea_groups <- data_list[[condition]][[2]]
        dea_design <- data_list[[condition]][[3]]

        dds <- DESeq2::DESeqDataSetFromMatrix(countData = dea_data,
                                              colData = dea_groups,
                                              design = dea_design)

        # pre-filter low-count genes
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep,]

        # generate DEA results
        dds <- DESeq(dds)
        res <- results(dds)
        res_df <- as.data.frame(res)
        res_df$canon_entrez <- rownames(res)

        # filter out NA values in padj
        res_df %<>% filter(!is.na(padj))

        DEA_results_list <- c()
        DEA_results_list[['data']] <- dea_data
        DEA_results_list[['groups']] <- dea_groups
        DEA_results_list[['design']] <- dea_design
        DEA_results_list[['DEGs']] <- res_df

        return(DEA_results_list)

    }

    do_microarray_DEA <- function(condition){

        dea_data <- data_list[[condition]][[1]]
        dea_groups <- data_list[[condition]][[2]]
        dea_design <- data_list[[condition]][[3]]

        # fit linear model using provided design and contrast matrices
        fit <- limma::lmFit(dea_data, dea_design)
        fit2 <- limma::eBayes(fit, robust = TRUE, trend = TRUE) # robust, trend

        results <- limma::decideTests(fit2)
        summary <- summary(results)

        final <- limma::topTable(fit2, number = Inf)
        final$canon_entrez <- rownames(final)
        rownames(final) <- NULL
        final %<>% dplyr::rename(., pvalue = P.Value)

        # save results
        DEA_results_list <- c()
        DEA_results_list[['data']] <- dea_data
        DEA_results_list[['groups']] <- dea_groups
        DEA_results_list[['design']] <- dea_design
        DEA_results_list[['DEGs']] <- final
        DEA_results_list[['summary']] <- summary

        return(DEA_results_list)
    }

    if (platform == 'rnaseq'){

        future::plan(multisession, workers = processes)
        output_list <- furrr::future_map(conditions, ~ do_rnaseq_DEA(.))
        names(output_list) <- conditions

    }

    if (platform == 'microarray'){

        future::plan(multisession, workers = processes)
        output_list <- furrr::future_map(conditions, ~ do_microarray_DEA(.))
        names(output_list) <- conditions

    }

    return(output_list)
}

