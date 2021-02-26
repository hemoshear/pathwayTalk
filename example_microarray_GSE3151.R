# Taylor Derby Pourtaheri

# purpose: import public microarray dataset from NCBI via GEO query

library(GEOquery)
library(Biobase)
library(biomaRt)
library(dplyr)
library(magrittr)
library(future)

options(stringsAsFactors = FALSE)

# source('R/1_diffExpression.R')
# source('R/2_pathwayEnrichment.R')
# source('R/3_pathwayCrosstalk.R')
# source('R/4_classification.R')
# source('R/5_networkCharacterization.R')
# source('R/6_crosstalkInhibition.R')
# source('R/7_networkWrapper.R')

library(devtools)
devtools::install_github('hemoshear/pathwayTalk', ref = 'wrapper_dev')
library(pathwayTalk)


# import data object ------------------------------------------------------

# # dataset on NCBI: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3151
#
# # define accension ID for the dataset
# GEO_acc_ID <- 'GSE3151'
#
# # import data as ExpressionSet object - see ?ExpressionSet
# gse <- getGEO(GEO = GEO_acc_ID)
# gse <- gse$GSE3151_series_matrix.txt.gz
#
# # save RDS for future use
# saveRDS(object = gse, file = 'data/GSE3151_expression_set.RDS')
#
# # reading and subsetting --------------------------------------------------

# read in the data
E <- readRDS(file = 'data/GSE3151_expression_set.RDS')

# create a cleaned metadata object
metadata <- pData(E)

# we only need two columns from the metadata
metadata <- dplyr::select(metadata, title, geo_accession)

# create a new treatment column in meta
metadata <- dplyr::mutate(metadata,
                          treatment = gsub("\\-.*", "", title),
                          treatment = toupper(treatment))

#
# # log2 transformation of the expression matrix
# Biobase::exprs(E) <- log2(Biobase::exprs(E))
#
# # quantile normalization
# E_qt <- limma::normalizeBetweenArrays(exprs(E), method = "quantile")
#
# # collapse probe rows to genes
# E_genes <- WGCNA::collapseRows(datET = E_qt,
#                              rowGroup = feature_data$ENTREZ_GENE_ID,
#                              rowID = rownames(E_qt),
#                              method = 'MaxMean')
#
# saveRDS(E_genes, 'data/GSE3151_processed_e_mat.RDS')

# step 1: differential expression analysis  -------------------------------

# read in the processed data
processed_data <- readRDS('data/GSE3151_processed_e_mat.RDS')

processed_exprs <- processed_data$datETcollapsed

# define groups
groups <- metadata[,c('geo_accession', 'treatment')] %>%
    setNames(c('sample_id', 'group'))
groups$group %<>% factor %>% relevel(ref = 'GFP')

# do differential expression analysis
DEG <- diffExpression(expression_matrix = processed_exprs,
                      groups = groups,
                      platform = 'microarray')

# step 2: pathway enrichment analysis -------------------------------------

# extract DEGs from diffExpression output
DEG_dfs <- purrr::map(DEG, ~.$DEGs)

# import pathways list
reactome_pathways <- readRDS('data/reactome_pathways.RDS')

enriched <- fisherPathwayEnrichment(DEG_dfs, gene_alpha=0.05,
                                    pathways = reactome_pathways)

purrr::map(enriched, ~nrow(.))
purrr::map(enriched, ~sum(is.infinite(.$estimate)))
enriched_short <- purrr::map(enriched, ~filter(., !is.infinite(estimate)))
purrr::map(enriched_short, ~nrow(.))

# try wrapper --------------------------------------------------------

all(names(DEG) == names(enriched_short))

results <- purrr::map2(DEG, enriched_short,
                       ~ networkWrapper(expression_matrix = .x$data,
                                        groups = .x$groups,
                                        DEPs = .y,
                                        pathways = reactome_pathways,
                                        pathway_alpha = 0.001,
                                        lambda = 0.01))

# extract results ---------------------------------------------------------

full_networks <- purrr::map(results, ~ .$crosstalk_inhibition$full_network)
pruned_networks <- purrr::map(results, ~ .$crosstalk_inhibition$pruned_network)

full_edges <- purrr::map(results, ~ .$crosstalk_inhibition$full_network_results)
pruned_edges <- purrr::map(results, ~ .$crosstalk_inhibition$pruned_network_results)

plot(full_networks[[1]])
plot(pruned_networks[[1]])

purrr::map(full_edges, ~nrow(.)) %>% unlist
purrr::map(pruned_edges, ~nrow(.)) %>% unlist

full_edges_all_contrasts <- bind_rows(full_edges, .id = 'contrast')
pruned_edges_all_contrasts <- bind_rows(pruned_edges, .id = 'contrast')

# annotate results ------------------------------------

# using reactome.db:
r_pathway_names <- as.list(reactome.db::reactomePATHNAME2ID)
r_pathway_names_df <- tibble::enframe(r_pathway_names) %>%
    setNames(c('pathway_name', 'name'))
r_pathway_names_df <- r_pathway_names_df[grepl('HSA', r_pathway_names_df$name),]
r_pathway_names_df$name %<>% gsub('\\-', '\\.', .)

pruned_edges_all_contrasts$pathway1 %in% r_pathway_names_df$name
pruned_edges_all_contrasts$pathway2 %in% r_pathway_names_df$name

names1 <- r_pathway_names_df[match(pruned_edges_all_contrasts$pathway1,
                                   r_pathway_names_df$name), 'pathway_name']
names2 <- r_pathway_names_df[match(pruned_edges_all_contrasts$pathway2,
                                   r_pathway_names_df$name), 'pathway_name']

pruned_edges_all_contrasts$pathway1_name <- names1$pathway_name
pruned_edges_all_contrasts$pathway2_name <- names2$pathway_name


# add associated genes ----------------------------------------------------

pathways <- reactome_pathways
names(pathways) %<>% gsub('\\-', '\\.', .)

genes <- purrr::map2_chr(pruned_edges_all_contrasts$pathway1,
                        pruned_edges_all_contrasts$pathway2,
                        ~ ifelse(length(intersect(pathways[[.x]], pathways[[.y]])) == 0,
                                 'None',
                                 intersect(pathways[[.x]], pathways[[.y]])))


genes %<>% unique %>% setdiff('None')


mart = useMart('ensembl', 'hsapiens_gene_ensembl')

key = getBM(mart = mart, attributes = c('entrezgene_id', 'external_gene_name',
                                        'gene_biotype', 'description',
                                        'arrayexpress'))
key$entrezgene_id %<>% as.character

genes_df <- data.frame('entrezgene_id' = genes)
genes_df %<>% dplyr::left_join(key)




