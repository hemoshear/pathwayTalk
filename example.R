library(magrittr)
# library(pathwayTalk)

source('R/1_diffExpression.R')
source('R/2_pathwayEnrichment.R')
source('R/3_pathwayCrosstalk.R')
source('R/4_classification.R')
source('R/5_networkCharacterization.R')
source('R/6_crosstalkInhibition.R')

options(stringsAsFactors = FALSE)

# step 0: data pre-processing ---------------------------------------------

data <- readRDS(file = 'data/onco_data.RDS')

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
exprs(data) <- limma::normalizeBetweenArrays(exprs(data), method=norm_method_test)

# log2 transformation of intensity values
exprs(data) <- log2(exprs(data))

# select the expression matrix from the ExpressionSet data (probe-level)
probes <- exprs(data)

# generate vector of gene IDs:
gene_ids <- data@featureData@data$ENTREZ_GENE_ID

# select collapse method
collapse_method <- 'MaxMean'

# collapse probe rows to genes
genes <- WGCNA::collapseRows(datET = probes,
                             rowGroup = gene_ids,
                             rowID = rownames(probes),
                             method = collapse_method)

processed_expression_data <- genes$datETcollapsed

# step 1: differential expression analysis  -------------------------------

# generate design and contrast matrices:
design_mat <- model.matrix(~ 0 + data$title)
colnames(design_mat) %<>% gsub('data\\$title', '', .)

contrast_mat <- limma::makeContrasts(BCAT-GFP,
                              E2F3-GFP,
                              MYC-GFP,
                              RAS-GFP,
                              SRC-GFP,
                              levels = design_mat)


DEG <- diffExpression(data = processed_expression_data,
                    gene_ids = gene_ids,
                    design_mat = design_mat,
                    contrast_mat = contrast_mat)
# ^ runtime: about 2 minutes

# saveRDS(DEG, 'onco_DEA.RDS')


# step 2: pathway enrichment analysis -------------------------------------
enriched <- fisherPathwayEnrichment(DEG, gene_alpha=0.001, pathway_alpha=0.001)

# step 3: pathway crosstalk -----------------------------------------------

# # ----
#
# # need to add 'genes' for pathwayCrosstalkParallel?
#     # currently only defined within the diffExpression function
#     # for now, define here:
#
# # # begin with log2 transformation of intensity values

# Need a named vector that maps sample identifiers to treatment groups
pdata <- phenoData(data)
pdata@data$sample_id <- rownames(pdata@data)
treatment_map <- pdata@data$title
names(treatment_map) <- pdata@data$sample_id

# Get data frame of expression values
processed_expression_data %<>% as.data.frame()
processed_expression_data$entrez_full <- rownames(processed_expression_data)
# I take the first Entrez ID if there are multiple to make this simpler -- entrez is
# a required column in this Data Frame
processed_expression_data$entrez <- strsplit(processed_expression_data$entrez_full, split = '///') %>%
    purrr::map_chr(~ .[1])

# expr is not an argument to pathwayCrosstalk, but needs to be defined in the global environment.
# This function runs orders of magnitude faster this way
expr <- processed_expression_data

# # Each element is a matrx of discriminating scores.
ct_feature_matrix <- pathwayCrosstalk(enriched, treatment_map,
                                      processes = 4)


# steps 4 and 5: classification and network construction -------------------------------

# purpose: select the best discriminating pairs of crosstalking pathways

# beginning with a matrix of discrimining scores for each pathway pair, for each sample
# discriminating score: quantification of pathway cross-talk between pathways
# enriched with subtype-derived DEGs

# input: ct_feature matrix and transpose so that features are columns:

groups <- data$title # phenotypes


# split feature matrix by subtype
matrices_by_phenotype <- byPhenotype(input_matrix = ct_feature_matrix,
                                     phenotypes = groups,
                                     reference_condition = 'GFP')


# generate networks for all subtypes - use subtypeNetwork()
subtype_list <- names(matrices_by_phenotype$phenotypes)
matrices_by_phenotype[[1]][1:5]
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

# MYC example
final_results[['MYC']]

