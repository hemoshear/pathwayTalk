# testing -----------------------------------------------------------------
library(Biobase)
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

b %<>% entrez_to_hgnc()
# do pathway enrichment with Fisher's exact test
tests <- purrr::map(b, ~ fisherPathwayEnrichment(., fdr=0.01))
names(tests) <- names(b)

# keep enriched pathways for each cancer subtype.
sig_pathways <- purrr::map(tests, ~ dplyr::filter(., p < 0.05))
names(sig_pathways) %<>% gsub(' \\- GFP', '', .)
sig_pathways[[4]] %>% head

ct <- pathwayCrosstalk(sig_pathways, data_test)

# ct[[1]] %>% dim
# ctalk <- purrr::map(sig_pathways, ~ pathwayCrosstalk(., data_test))
