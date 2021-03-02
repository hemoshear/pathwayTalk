# Taylor Derby Pourtaheri

# purpose: import public microarray dataset from NCBI via GEO query

library(GEOquery)
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

full_networks <- purrr::map(results, ~ .$full_network)
pruned_networks <- purrr::map(results, ~ .$pruned_network)

full_edges <- purrr::map(results, ~ .$full_network_results)
pruned_edges <- purrr::map(results, ~ .$pruned_network_results)

plot(full_networks[[1]])
plot(pruned_networks[[1]])

purrr::map(full_edges, ~nrow(.)) %>% unlist
purrr::map(pruned_edges, ~nrow(.)) %>% unlist

# full_edges_all_contrasts <- bind_rows(full_edges, .id = 'contrast')
# pruned_edges_all_contrasts <- bind_rows(pruned_edges, .id = 'contrast')

# annotate results ------------------------------------

# using reactome.db:
r_pathway_names <- as.list(reactome.db::reactomePATHNAME2ID)
r_pathway_names_df <- tibble::enframe(r_pathway_names) %>%
    setNames(c('pathway_name', 'name'))
r_pathway_names_df <- r_pathway_names_df[grepl('HSA', r_pathway_names_df$name),]
r_pathway_names_df$name %<>% gsub('\\-', '\\.', .)

annotateNetwork <- function(network_df, pathway_df){

    res <- network_df

    names1 <- pathway_df[match(res$pathway1,
                               pathway_df$name), 'pathway_name']
    names2 <- pathway_df[match(res$pathway2,
                               pathway_df$name), 'pathway_name']

    res$pathway1_name <- names1$pathway_name
    res$pathway2_name <- names2$pathway_name

    return(res)
}


pathway_results <- purrr::map(pruned_edges,
                              ~annotateNetwork(network_df = .,
                                               pathway_df = r_pathway_names_df))


# pruned_edges_all_contrasts$pathway1 %in% r_pathway_names_df$name
# pruned_edges_all_contrasts$pathway2 %in% r_pathway_names_df$name

# add associated genes ----------------------------------------------------

mart = useMart('ensembl', 'hsapiens_gene_ensembl')
key = getBM(mart = mart, attributes = c('entrezgene_id', 'external_gene_name',
                                        'gene_biotype', 'description',
                                        'arrayexpress'))
key$entrezgene_id %<>% as.character

pathways <- reactome_pathways
names(pathways) %<>% gsub('\\-', '\\.', .)


extractNetworkGenes <- function(network_df, pathway_list, gene_key){

    results <- c()

    for (i in 1:nrow(network_df)){

        pathway1 <- network_df$pathway1[i]
        pathway2 <- network_df$pathway1[i]

        intersection <- intersect(pathway_list[[pathway1]],
                                  pathway_list[[pathway2]])

        if(length(intersection) != 0){
            results[[i]] <- intersection
        }
        else{
            results[[i]] <- 'None'
        }

    }

    genes <- unlist(results)
    genes %<>% unique %>% setdiff('None')
    genes_df <- data.frame('entrezgene_id' = genes)
    genes_df %<>% left_join(key)

    return(genes_df)


}


gene_results <- purrr::map(pruned_edges,
                           ~extractNetworkGenes(network_df = .,
                                                pathway_list = pathways,
                                                gene_key = key))


# save results

openxlsx::write.xlsx(full_edges, 'results/microarray/final/example_GS3151_full_edges.xlsx')
openxlsx::write.xlsx(pruned_edges, 'results/microarray/final/example_GS3151_pruned_edges.xlsx')
openxlsx::write.xlsx(pathway_results, 'results/microarray/final/example_GS3151_pathway_results.xlsx')
openxlsx::write.xlsx(gene_results, 'results/microarray/final/example_GS3151_gene_results.xlsx')












