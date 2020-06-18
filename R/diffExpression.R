
library(tidyverse)
library(magrittr)
library(limma)
library(WGCNA)

# https://www.rpubs.com/jlubieni/13450

# purpose: a function to conduct differential expression analysis
    # for a given dataset - beginning with ExpressionSet class

# this function works with an ExpressionSet with a lockedEnvironment storage mode

# pre-processing completed prior to using function (normalization, removing control sequences, etc.)


# write function to do DEA ------------------------------------------------


diffExpression <- function(data,
                           collapse_method,
                           gene_ids,
                           design_mat,
                           contrast_mat) {

    # begin with log2 transformation of intensity values
    data <- log2(data)

    # collapse probe rows to genes
    probes <- data

    genes <- WGCNA::collapseRows(datET = probes,
                                 rowGroup = gene_ids,
                                 rowID = rownames(probes),
                                 method = collapse_method)

    data <- genes$datETcollapsed

    # fit linear model using provided design and contrast matrices
    initial_fit <- lmFit(data, design_mat)
    temp_fit <- contrasts.fit(initial_fit, contrast_mat)
    fit <- eBayes(temp_fit)

    contrast_names <- colnames(contrast_mat)

    # output
    # results <- decideTests(fit)

    results <-c()

    for(i in 1:length(contrast_names)){

        cont <- contrast_names[i]

        results[[i]] <- topTable(fit, coef = i, number = Inf)

    }

    names(results) <- contrast_names

    results

}


# testing -----------------------------------------------------------------

data_test <- readRDS(file = 'onco_data.RDS')

# trim sample id numbers from treatment groups:
data_test$title <- gsub('\\-\\d+', '', data_test$title)
data_test$title %<>% toupper()
table(data_test$title)

data_test$title %<>% factor() %>% relevel(ref = 'GFP')

# select a normalization method
norm_method_test <- 'quantile'
# single channel: none, scale, quantile, cyclicloess
# two-channel: above, plus Aquantile, Gquantile, Rquantile, Tquantile

# normalize using provided normalization method
exprs(data_test) <- normalizeBetweenArrays(exprs(data_test), method=norm_method_test)

# select collapse method
collapse_method_test <- 'MaxMean'

# generate vector of gene IDs:
gene_ids_test <- data_test@featureData@data$ENTREZ_GENE_ID

# generate design and contrast matrices:

design_mat_test <- model.matrix(~ 0 + data_test$title)
colnames(design_mat_test) %<>% gsub('data_test\\$title', '', .)

contrast_mat_test <- makeContrasts(BCAT-GFP,
                              E2F3-GFP,
                              MYC-GFP,
                              RAS-GFP,
                              SRC-GFP,
                              levels = design_mat_test)

# select the expression matrix from the ExpressionSet data
data_test <- exprs(data_test)


b <- diffExpression(data = data_test,
                    collapse_method = collapse_method_test,
                    gene_ids = gene_ids_test,
                    design_mat = design_mat_test,
                    contrast_mat = contrast_mat_test)


data <- data_test
collapse_method <- collapse_method_test
gene_ids <- gene_ids_test
design_mat <- design_mat_test
contrast_mat <- contrast_mat_test
