# testing -----------------------------------------------------------------
library(Biobase)
library(reactome.db)
source('R/diffExpression.R')
source('R/pathwayEnrichment.R')
source('R/pathwayCrosstalk.R')

data_test <- readRDS(file = 'onco_data.RDS')

# trim sample id numbers from treatment groups:
data_test$title <- gsub('\\-\\d+', '', data_test$title)
data_test$title %<>% toupper()
# table(data_test$title)

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
b <- diffExpression(data = exprs(data_test),
                    collapse_method = collapse_method_test,
                    gene_ids = gene_ids_test,
                    design_mat = design_mat_test,
                    contrast_mat = contrast_mat_test)

for (i in 1:length(b)){
    entrez_id <- rownames(b[[i]])
    b[[i]] %<>% as.data.frame()
    b[[i]]$entrez <- entrez_id
    b[[i]]$canon_entrez <- strsplit(b[[i]]$entrez, split = '///') %>%
        purrr::map_chr(~ .[1])
}

# do pathway enrichment with Fisher's exact test
tests <- purrr::map(b, ~ fisherPathwayEnrichment(., alpha=0.001))
names(tests) <- names(b)

# keep enriched pathways for each cancer subtype.
sig_pathways <- purrr::map(tests, ~ dplyr::filter(., adj_p < 0.001))
names(sig_pathways) %<>% gsub(' \\- GFP', '', .)

# create one data frame with enriched pathways across cancer subtypes
for (i in 1:length(sig_pathways)) {
    sig_pathways[[i]]$subtype <- names(sig_pathways)[i]
}
sig_pathways %<>% dplyr::bind_rows()

# Each element is a matrx of discriminating scores.
ct <- pathwayCrosstalkParallel(sig_pathways, data_test, processes=40)

# Convert each matrix into a single vector and combine these to create a matrix
# where rows are pathway pairs and columns are samples.
ct_feature <- list()
matrix2vec <- function(m) {
    vec <- c()
    pathway_pair <- c()
    for (i in 1:nrow(m)) {
        for (j in 1:ncol(m)) {
            if (i > j) {
                vec <- c(vec, m[i,j])
                pathway_pair <- c(pathway_pair, paste0(rownames(m)[i], ',', colnames(m)[j]))
            }
        }
    }

    names(vec) <- pathway_pair
    vec
}
for (rna_sample in 1:length(ct)) {
    ct_feature[[rna_sample]] <- matrix2vec(ct[[rna_sample]])
} 
ct_feature_matrix <- do.call('cbind', ct_feature)
rownames(ct_feature_matrix) <- names(ct_feature[[1]])
colnames(ct_feature_matrix) <- names(ct)
# ct_feature_matrix will be input to classifier in next step


