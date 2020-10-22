
source('R/diffExpression.R') # diffExpression()
source('R/pathwayEnrichment.R') # fisherPathwayEnrichment()
source('R/pathwayCrosstalk.R') # pathwayCrosstalk(), pathwayCrosstalkParallel()
source('R/classification.R') # byPhenotype(), selectTopFeatures(), subtypeNetwork()
source('R/networkCharacterization.R') # characterizeNetworks()
source('R/crosstalkInhibition.R') # networkEfficiency(), crosstalkInhibition, characterizeResults()

options(stringsAsFactors = FALSE)

# step 0: data pre-processing ---------------------------------------------

data <- readRDS(file = 'onco_data.RDS')

# trim sample id numbers from treatment groups:
data$title <- gsub('\\-\\d+', '', data$title)
data$title %<>% toupper()
table(data$title)

data$title %<>% factor() %>% relevel(ref = 'GFP')

# select a normalization method
norm_method_test <- 'quantile'
# single channel: none, scale, quantile, cyclicloess
# two-channel: above, plus Aquantile, Gquantile, Rquantile, Tquantile

# normalize using provided normalization method
exprs(data) <- normalizeBetweenArrays(exprs(data), method=norm_method_test)


# step 1: differential expression analysis  -------------------------------

# select collapse method
collapse_method <- 'MaxMean'

# generate vector of gene IDs:
gene_ids <- data@featureData@data$ENTREZ_GENE_ID

# generate design and contrast matrices:
design_mat <- model.matrix(~ 0 + data$title)
colnames(design_mat) %<>% gsub('data\\$title', '', .)

contrast_mat <- makeContrasts(BCAT-GFP,
                                   E2F3-GFP,
                                   MYC-GFP,
                                   RAS-GFP,
                                   SRC-GFP,
                                   levels = design_mat)

# select the expression matrix from the ExpressionSet data
expression_data <- exprs(data)

DEG <- diffExpression(data = expression_data,
                    collapse_method = collapse_method,
                    gene_ids = gene_ids,
                    design_mat = design_mat,
                    contrast_mat = contrast_mat)
# ^ runtime: about 2 minutes

# saveRDS(DEG, 'onco_DEA.RDS')


# step 2: pathway enrichment analysis -------------------------------------

# do pathway enrichment with Fisher's exact test
tests <- purrr::map(DEG, ~ fisherPathwayEnrichment(., alpha=0.001))
# ^ runtime: about 2 minutes
names(tests) <- names(DEG)

# keep enriched pathways for each cancer subtype.
sig_pathways <- purrr::map(tests, ~ dplyr::filter(., adj_p < 0.001))
names(sig_pathways) %<>% gsub(' \\- GFP', '', .)

purrr::map(sig_pathways, ~nrow(.))

# create one data frame with enriched pathways across cancer subtypes
for (i in 1:length(sig_pathways)) {
    sig_pathways[[i]]$subtype <- names(sig_pathways)[i]
}
sig_pathways %<>% dplyr::bind_rows()

count(sig_pathways, subtype)

# step 3: pathway crosstalk -----------------------------------------------

# # ----
#
# # need to add 'genes' for pathwayCrosstalkParallel?
#     # currently only defined within the diffExpression function
#     # for now, define here:
#
# # # begin with log2 transformation of intensity values
# probes <- log2(exprs(data))
#
# # # collapse probe rows to genes
# genes <- WGCNA::collapseRows(datET = probes,
#                              rowGroup = gene_ids,
#                              rowID = rownames(probes),
#                              method = collapse_method)
#
# # ----
#
# # Each element is a matrx of discriminating scores.
# ct <- pathwayCrosstalkParallel(sig_pathways, data)
# # ^ runtime: > 10 minutes
#
# # Convert each matrix into a single vector and combine these to create a matrix
# # where rows are pathway pairs and columns are samples.
# ct_feature <- list()
# matrix2vec <- function(m) {
#     vec <- c()
#     pathway_pair <- c()
#     for (i in 1:nrow(m)) {
#         for (j in 1:ncol(m)) {
#             if (i > j) {
#                 vec <- c(vec, m[i,j])
#                 pathway_pair <- c(pathway_pair, paste0(rownames(m)[i], ',', colnames(m)[j]))
#             }
#         }
#     }
#
#     names(vec) <- pathway_pair
#     vec
# }
# for (rna_sample in 1:length(ct)) {
#     ct_feature[[rna_sample]] <- matrix2vec(ct[[rna_sample]])
# }
#
# ct_feature_matrix <- do.call('cbind', ct_feature)
# rownames(ct_feature_matrix) <- names(ct_feature[[1]])
# colnames(ct_feature_matrix) <- names(ct)

# ct_feature_matrix will be input to classifier in next step

# steps 4 and 5: classification and network construction -------------------------------

# purpose: select the best discriminating pairs of crosstalking pathways

# beginning with a matrix of discrimining scores for each pathway pair, for each sample
# discriminating score: quantification of pathway cross-talk between pathways
# enriched with subtype-derived DEGs

# input: ct_feature matrix and transpose so that features are columns:

ct_feature_matrix <- readRDS('pathway_crosstalk_feature_matrix.RDS')
ct_feature_matrix_t <- ct_feature_matrix %>% t() %>% data.frame()


groups <- data$title # phenotypes


# split feature matrix by subtype
matrices_by_phenotype <- byPhenotype(input_matrix = ct_feature_matrix_t,
                                     phenotypes = groups,
                                     reference_condition = 'GFP')


# generate networks for all subtypes - use subtypeNetwork()
subtype_list <- names(matrices_by_phenotype$phenotypes)

networks_by_subtype <- function(subtypes){

    networks_list <- c()

    for (subtype in subtypes) {

        subtype_matrix <- matrices_by_phenotype$matrices[[subtype]]
        subtype_phenotypes <- matrices_by_phenotype$phenotypes[[subtype]]

        subtype_network <- subtypeNetwork(feature_matrix = subtype_matrix,
                                          sample_phenotype = subtype_phenotypes,
                                          alpha = 1,
                                          lambda = 0.01,
                                          output_graph = TRUE,
                                          model_evaluation = FALSE)

        subtype_graph <- subtype_network$graph
        networks_list[[subtype]] <- subtype_graph

    }

    return(networks_list)
}

networks <- networks_by_subtype(subtype_list)

plot(networks[['MYC']])
plot(networks[['SRC']])
plot(networks[['BCAT']])
plot(networks[['E2F3']])
plot(networks[['RAS']])



# step 5b: characterize networks ------------------------------------------

# purpose: characterize the vertices present in the pathways and export as xslx file

network_results <- characterizeNetworks(networks_list = networks)

# step 6: pathway crosstalk inhibition ------------------------------------

# purpose: select the best discriminating pairs of crosstalking pathways

significant_crosstalks <- crosstalkInhibition(networks)
final_results <- characterizeResults(significant_crosstalks)










