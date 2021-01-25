library(magrittr)
library(dplyr)
library(GEOquery)
library(Rsubread)
library(limma)
library(edgeR)
library(pathwayTalk)

options(stringsAsFactors = FALSE)

# RNAseq workflow for pathwayTalk

# data import ----------------------------------------------------

# import expression data file (BDG FBX experiment, courtesy of ML)
exp_list <- readRDS('BDG0202_1.experiment_data.RDS')

meta <- exp_list[['metadata']]
counts <- exp_list[['gene_counts']]
contrast_list <- exp_list[['contrast_list']]

# inspect counts data
# View(counts[1:10, 1:10])
# rows are genes, columns are samples

# subset contrasts to comparisons to VEH-DMSO-ICS
contrast_list <- contrast_list[1:11]

# # select only samples from main waveform treatment (ICS)
# meta %<>% filter(waveform == 'ICS')
# ids <- unique(meta$ea_code) %>% as.character
# counts <- counts[,ids]

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

controls <- c('VEH_NaOH', 'VEH_DMSO_ICS', 'VEH_DMSO_CCA')
tx_levels <- levels(dge$samples$group) %>% setdiff(controls) %>% c(controls,.)

design_mat <- model.matrix(~ 0 + dge$samples$treatment_label)
colnames(design_mat) %<>% gsub('dge\\$samples\\$treatment_label', '', .)
rownames(design_mat) <- dge$samples$ea_code

# # remove rows with consistently zero or low counts
nrow(dge$counts) # 19,756
keep <- filterByExpr(dge$counts, design_mat)
dge <- dge[keep,,keep.lib.sizes=FALSE]
nrow(dge$counts) # 15,075

# scale normalization with the TMM method
# "TMM normalization is applied to this dataset to account for
# compositional difference between the libraries."
dge <- calcNormFactors(dge)


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

# use local version of diffExpression()
source('R/1_diffExpression.R')

DEG <- diffExpression(data = v,
                      gene_ids = gene_ids,
                      design_mat = design_mat,
                      contrast_mat = contrast_mat)


# alternative: DEA with edgeR (negative binomial)  ----------------------

# # calculate dispersions - trended is the default
# # "edgeR uses the Cox-Reid profile-adjusted likelihood (CR) method in estimating dispersions"
# dge <- estimateGLMTrendedDisp(dge, design_mat)
#
# fit <- glmQLFit(dge, design_mat)
#
# qlf.2vs1 <- glmQLFTest(fit, coef=2)
# topTags(qlf.2vs1)
#
# # QL (quasi-likelihood) dispersions visualization:
# pdf("results/RNAseq/QLdisp.pdf", width = 15, height = 20)
# plotQLDisp(fit)
# dev.off()
#
# topTags(qlf.2vs1)
#
# # enriched <- fisherPathwayEnrichment(DEG, gene_alpha=0.01, pathway_alpha=0.001)

# step 2: pathway enrichment analysis -------------------------------------

# some contrasts are showing no enriched pathways - need to address this in the function

hist(DEG[[1]]$P.Value)
hist(DEG[[2]]$P.Value)
hist(DEG[[3]]$P.Value)
hist(DEG[[4]]$P.Value)
hist(DEG[[5]]$P.Value)
hist(DEG[[6]]$P.Value)
hist(DEG[[7]]$P.Value)
hist(DEG[[8]]$P.Value)
hist(DEG[[9]]$P.Value)
hist(DEG[[10]]$P.Value)
hist(DEG[[11]]$P.Value)

# DEG[['ALP_50-VEH_DMSO_ICS']] <- NULL
# DEG[['FBX_10-VEH_DMSO_ICS']] <- NULL
# DEG[['OXP_50-VEH_DMSO_ICS']] <- NULL
# DEG[['OXP_500-VEH_DMSO_ICS']] <- NULL
# DEG[['TMX_10-VEH_DMSO_ICS']] <- NULL
# DEG[['TPS_10-VEH_DMSO_ICS']] <- NULL

gene_alpha <- 0.05
pathway_alpha <- 0.05

# local function definition
source('R/2_pathwayEnrichment.R')

# define functions in local pathway enrichment script - use p val instead of adj p val
enriched <- fisherPathwayEnrichment(DEG, gene_alpha=gene_alpha,
                                    pathway_alpha=pathway_alpha)

# replace infinite values
fix <- which(is.infinite(enriched$estimate))
enriched$estimate[fix] <- 100

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
                                      processes = 4)


# getting lots of NAs here - may be due to infinite values in enriched estimates?
# does not appear to be the case


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

plot(networks[['TPS_100']])
plot(networks[['FBX_100']])
plot(networks[['TMX_100']])
plot(networks[['OXP_500']])
plot(networks[['OXP_500']])


# step 5b: characterize networks ------------------------------------------

# purpose: characterize the vertices present in the pathways and export as xslx file

network_results <- characterizeNetworks(networks_list = networks)

# step 6: pathway crosstalk inhibition ------------------------------------

# purpose: select the best discriminating pairs of crosstalking pathways

significant_crosstalks <- crosstalkInhibition(networks)
final_results <- characterizeResults(significant_crosstalks)

# MYC example
final_results[['MYC']]



