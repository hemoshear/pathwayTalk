library(magrittr)
library(dplyr)
library(tidyr)
library(gtools)
library(GEOquery)
library(Rsubread)
library(igraph)
library(limma)
library(edgeR)
library(DESeq2)
library(pathwayTalk)
library(furrr)
library(future)
library(biomaRt)

options(stringsAsFactors = FALSE)


source('R/1_diffExpression.R')
source('R/2_pathwayEnrichment.R')
source('R/3_pathwayCrosstalk.R')
source('R/4_classification.R')
source('R/5_networkCharacterization.R')
source('R/6_crosstalkInhibition.R')
source('R/7_networkWrapper.R')

# data import ----------------------------------------------------

# import expression data file (BDG FBX experiment, courtesy of ML)
exp_list <- readRDS('data/BDG0202_1.experiment_data.RDS')

meta <- exp_list[['metadata']]
counts <- exp_list[['gene_counts']]
contrast_list <- exp_list[['contrast_list']]

# select certain cell type (SMCs) and waveform (ICS)
meta %<>% filter(cell_type == 'SMC'  & waveform == 'ICS')
ids <- unique(meta$ea_code) %>% as.character
cols <- colnames(counts) %in% ids
counts <- counts[,cols]

# step 0: data pre-processing ---------------------------------------------

# select sample phenotypes and create phenotype xref dataframe
groups <- meta[,c('ea_code', 'treatment_label')] %>% setNames(c('sample_id', 'group'))
groups$sample_id %<>% as.character

# set reference level for reatment condition (group)
groups$group %<>% factor %>% relevel(ref = 'VEH_DMSO_ICS')

# DESeq2 expects *un-normalized* counts
# counts need to be positive integers - according to ML, we can round the values
counts_rounded <- round(counts)

# step 1: data splitting and differential expression  -------------------------------

DEG <- diffExpression(expression_matrix = counts_rounded,
                      groups = groups,
                      platform = 'rnaseq',
                      processes = 4)


# step 2: pathway enrichment analysis -------------------------------------

# select DEG results from output list
DEG_dfs <- purrr::map(DEG, ~.$DEGs)

# import pathways source and define gene significance level
reactome_pathways <- readRDS('data/reactome_pathways.RDS')


# run pathway enrichment analysis with Fisher's test
enriched <- fisherPathwayEnrichment(DEGs = DEG_dfs,
                                    gene_alpha = 0.05,
                                    pathways = reactome_pathways)


# replace infinite values?
# fix <- which(is.infinite(enriched$estimate))
# enriched$estimate[fix] <- 100

purrr::map(enriched, ~nrow(.))
purrr::map(enriched, ~sum(is.infinite(.$estimate)))
enriched_short <- purrr::map(enriched, ~filter(., !is.infinite(estimate)))
purrr::map(enriched_short, ~nrow(.))


# network wrapper --------------------------------------------------------

all(names(DEG) == names(enriched))

results <- purrr::map2(DEG, enriched_short,
                              ~ networkWrapper(expression_matrix = .x$data,
                                               groups = .x$groups,
                                               DEPs = .y,
                                               pathways = reactome_pathways,
                                               pathway_alpha = 0.05,
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
pruned_edges_all_contrasts$pathway1 %in% r_pathway_names_df$name

names1 <- r_pathway_names_df[match(pruned_edges_all_contrasts$pathway1,
                                          r_pathway_names_df$name), 'pathway_name']
names2 <- r_pathway_names_df[match(pruned_edges_all_contrasts$pathway2,
                                   r_pathway_names_df$name), 'pathway_name']

pruned_edges_all_contrasts$pathway1_name <- names1$pathway_name
pruned_edges_all_contrasts$pathway2_name <- names2$pathway_name


# add associated genes ----------------------------------------------------

pathways <- reactome_pathways
names(pathways) %<>% gsub('\\-', '\\.', .)

a <- pruned_edges_all_contrasts[1, 'pathway1']
b <- pruned_edges_all_contrasts[1, 'pathway2']

intersect(pathways[[a]], pathways[[b]])

genes_list <- c()

for (i in 1:nrow(pruned_edges_all_contrasts)){
    a <- pruned_edges_all_contrasts[i, 'pathway1']
    b <- pruned_edges_all_contrasts[i, 'pathway2']

    res <- intersect(pathways[[a]], pathways[[b]])

}

temp <- purrr::map2(pruned_edges_all_contrasts$pathway1,
                        pruned_edges_all_contrasts$pathway2,
                        ~ ifelse(length(intersect(pathways[[.x]], pathways[[.y]])) == 0,
                                 'None',
                                 paste0(intersect(pathways[[.x]], pathways[[.y]]))))

names(temp) <- pruned_edges_all_contrasts$contrast

genes <- unique(temp) %>% setdiff('None')


mart = useMart('ensembl', 'hsapiens_gene_ensembl')

key = getBM(mart = mart, attributes = c('entrezgene_id', 'external_gene_name',
                                        'gene_biotype', 'description',
                                        'arrayexpress'))
key$entrezgene_id %<>% as.character

genes_df <- data.frame('entrezgene_id' = genes)
genes_df %<>% dplyr::left_join(key)
