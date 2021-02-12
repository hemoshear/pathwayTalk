library(magrittr)
library(dplyr)
library(tidyr)
library(gtools)
library(GEOquery)
library(Rsubread)
library(igraph)
library(limma)
library(edgeR)
library(pathwayTalk)

options(stringsAsFactors = FALSE)

# RNAseq workflow for pathwayTalk

# data import ----------------------------------------------------

# import expression data file (BDG FBX experiment, courtesy of ML)
exp_list <- readRDS('data/BDG0202_1.experiment_data.RDS')

meta <- exp_list[['metadata']]
counts <- exp_list[['gene_counts']]
contrast_list <- exp_list[['contrast_list']]

# inspect counts data
# View(counts[1:10, 1:10])
# rows are genes, columns are samples


# select certain cell type (SMCs) and waveform (ICS)
meta %<>% filter(cell_type == 'SMC'  & waveform == 'ICS')
ids <- unique(meta$ea_code) %>% as.character
cols <- colnames(counts) %in% ids
counts <- counts[,cols]

unique(meta$treatment_label)

contrast_list <- contrast_list[c(2,3,4,5,6,7,8,9,11)]

# step 00: data pre-processing ---------------------------------------------

# select sample phenotypes and create phenotype xref dataframe
xref <- meta[,c('ea_code', 'treatment_label')]

# create DGEList from matrix of count values using edgeR package
colnames(counts) == xref$ea_code

# # save matrix if needed later
# counts_mat <- as.matrix(counts)
# class(counts_mat) <- 'numeric'

dge <- DGEList(counts = counts, samples = xref, group = xref$treatment_label)

# generate design matrix:
levels(dge$samples$group)

controls <- c('VEH_NaOH', 'VEH_DMSO_ICS')
tx_levels <- levels(dge$samples$group) %>% setdiff(controls) %>% c(controls,.)
dge$samples$treatment_label %<>% factor(levels = tx_levels)

design_mat <- model.matrix(~ 0 + dge$samples$treatment_label)
colnames(design_mat) %<>% gsub('dge\\$samples\\$treatment_label', '', .)
rownames(design_mat) <- dge$samples$ea_code

# # remove rows with consistently zero or low counts
nrow(dge$counts) # 19,756
keep <- filterByExpr(dge$counts, design_mat)
dge <- dge[keep,,keep.lib.sizes=FALSE]
nrow(dge$counts) # 14,386

# scale normalization with the TMM method
# "TMM normalization is applied to this dataset to account for
# compositional difference between the libraries."
dge <- calcNormFactors(dge)

# save RDS files for use in development down the line
# saveRDS(dge$counts, 'data/wrapper_dev_data/processed_BDG_expression_mat.RDS')

# save treatment condition metadata
# temp <- dge$samples
# temp$ea_code <- rownames(temp)
# saveRDS(temp[,c('ea_code', 'group')], 'data/wrapper_dev_data/processed_BDG_groups.RDS')

# step 0: data exploration --------------------------------------------------------

# # plotMDS:
#     # "examine the samples for outliers and for other relationships"
#     # "The function plotMDS produces a plot in which distances between samples
#     # correspond to leading biological coefficient of variation (BCV)
#     # between those samples"
#
# pdf("results/RNAseq/MDS.pdf", width = 15, height = 20)
# plotMDS(dge, col = c(rep(1, 12), rep(2, 14))) # color corr to phenotype
# dev.off()
#
# # Distances on an MDS plot of a DGEList object correspond to leading
# # log-fold-change between each pair of samples.
#
# # plotBCV:
#     # The square root of the common dispersion gives the
#     # coefficient of variation of biological variation
#
# # estimate NB dispersion:
#
# dge <- estimateDisp(dge, design_mat, robust = TRUE)
# dge$common.dispersion # 0.156
#
# pdf("results/RNAseq/BCV.pdf", width = 15, height = 20)
# plotBCV(dge)
# dev.off()


# step 1: DEA with voom  -------------------------------

# "When the library sizes are quite variable between samples,
# the voom approach is theoretically more powerful than limma-trend"

# voom transformation
v <- voom(dge, design_mat, plot = TRUE, span = 1)


# "If the data are very noisy, one can apply the same between-array
# normalization methods as would be used for microarrays, for example:"
# v_norm <- voom(dge, design_mat, plot=FALSE, normalize="quantile")


contrast_mat <- limma::makeContrasts(contrasts = contrast_list,
                                     levels = design_mat)

# after using voom, the dge object is compatible with limma

# use local version of diffExpression() - new arguments to ebayes
source('R/1_diffExpression.R')

DEG <- diffExpression(data = v,
                      design_mat = design_mat,
                      contrast_mat = contrast_mat)



# step 2: pathway enrichment analysis -------------------------------------

# some contrasts are showing no enriched pathways - need to address this in the function

names(DEG)

hist(DEG[[1]]$P.Value)
hist(DEG[[2]]$P.Value)
hist(DEG[[3]]$P.Value)
hist(DEG[[4]]$P.Value)
hist(DEG[[5]]$P.Value)
hist(DEG[[6]]$P.Value)
hist(DEG[[7]]$P.Value)
hist(DEG[[8]]$P.Value)
hist(DEG[[9]]$P.Value)

# DEG[3:4] <- NULL

# save a contrast for wrapper development
# saveRDS(DEG, 'data/wrapper_dev_data/processed_BDG_DEGs.RDS')

gene_alpha <- 0.01
pathway_alpha <- 0.01


# local function definition
source('R/2_pathwayEnrichment.R')

# define functions in local pathway enrichment script - use p val instead of adj p val
enriched <- fisherPathwayEnrichment(DEG, gene_alpha=gene_alpha,
                                    pathway_alpha=pathway_alpha)

# replace infinite values
# fix <- which(is.infinite(enriched$estimate))
# enriched$estimate[fix] <- 100

sum(is.infinite(enriched$estimate))
enriched %<>% filter(!is.infinite(estimate))

# save a contrast for wrapper development
# saveRDS(enriched, 'data/wrapper_dev_data/processed_BDG_EPs.RDS')

# step 3: pathway crosstalk -----------------------------------------------

# Need a named vector that maps sample identifiers to treatment groups
treatment_map <- xref$treatment_label
names(treatment_map) <- xref$ea_code

# I take the first Entrez ID if there are multiple to make this simpler -- entrez is
# a required column in this Data Frame

counts$entrez <- rownames(counts)

# expr is not an argument to pathwayCrosstalk, but needs to be defined in the global environment.
# This function runs orders of magnitude faster this way
expr <- counts

# # Each element is a matrx of discriminating scores.
ct_feature_matrix <- pathwayCrosstalk(enriched, treatment_map,
                                      processes = 6)

# check for NAs
sum(is.na(ct_feature_matrix))


# steps 4 and 5: classification and network construction -------------------------------

# purpose: select the best discriminating pairs of crosstalking pathways

# beginning with a matrix of discrimining scores for each pathway pair, for each sample
# discriminating score: quantification of pathway cross-talk between pathways
# enriched with subtype-derived DEGs

# input: ct_feature matrix and transpose so that features are columns:

groups <- dge$samples$group # phenotypes

# split feature matrix by subtype
matrices_by_phenotype <- byPhenotype(input_matrix = ct_feature_matrix,
                                     phenotypes = groups,
                                     reference_condition = 'VEH_DMSO_ICS')


# generate networks for all subtypes - use subtypeNetwork()
subtype_list <- as.character(purrr::map(strsplit(names(DEG), '-'), ~.[1]))

# source local classification function - fixed non-zero coef selection (which)
source('R/4_classification.R')

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

plot(networks[[1]])
plot(networks[[2]])
plot(networks[[3]])
plot(networks[[4]])
plot(networks[[5]])
plot(networks[[6]])
plot(networks[[7]])


# step 5b: characterize networks ------------------------------------------

# purpose: characterize the vertices present in the pathways and export as xslx file

source('R/5_networkCharacterization.R')
# unnest not found > removed

network_results <- characterizeNetworks(networks_list = networks)

# step 6: pathway crosstalk inhibition ------------------------------------

source('R/6_crosstalkInhibition.R')
# melt not found > gather
# issues with mixedsort as well


# purpose: select the best discriminating pairs of crosstalking pathways

significant_crosstalks <- crosstalkInhibition(networks)
final_results <- characterizeResults(significant_crosstalks)


names(final_results)

final_results[['FBX_10']] %>% View

final_results[['FBX_10']]$pathway_names
final_results[['FBX_100']]$pathway_names
final_results[['TMX_10']]$pathway_names
final_results[['TMX_100']]$pathway_names
final_results[['TPS_10']]$pathway_names
final_results[['TPS_100']]$pathway_names
final_results[['VEH_NaOH']]$pathway_names



