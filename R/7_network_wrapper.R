
# pathways <- as.list(reactome.db::reactomePATHID2EXTID)
# pathways <- pathways[grep('HSA', names(pathways))]
# saveRDS(pathways, 'data/wrapper_dev_data/reactome_pathways.RDS')

# import dev materials
expression_matrix <- readRDS('data/wrapper_dev_data/processed_BDG_expression_mat.RDS')
DEGs <- readRDS('data/wrapper_dev_data/processed_BDG_DEGs.RDS')
DEPs <- readRDS('data/wrapper_dev_data/processed_BDG_EPs.RDS')
pathways <- readRDS('data/wrapper_dev_data/reactome_pathways.RDS')
groups <- readRDS('data/wrapper_dev_data/processed_BDG_groups.RDS') %>%
    setNames(c('sample_id', 'group'))


# subset to one contrast
DEPs %<>% filter(contrast == 'FBX_100-VEH_DMSO_ICS')
DEGs <- DEGs[['FBX_100-VEH_DMSO_ICS']]
groups %<>% filter(group %in% c('VEH_DMSO_ICS', 'FBX_100'))
expression_matrix <- expression_matrix[,c(groups$sample_id)]

# pathwayCrosstalk function has been heavily edited here
    # need to compare results between the two workflows

# first, generate original DS matrix using original workflow
# library(pathwayTalk)
# counts <- data.frame(expression_matrix)
# colnames(counts) <- groups$sample_id
# counts$entrez <- rownames(counts)
# expr <- counts
# treatment_map <- groups$group
# names(treatment_map) <- groups$sample_id
# original_xtalk_mat <- pathwayCrosstalk(enriched_pathways = DEPs, treatment_map = treatment_map)
# saveRDS(original_xtalk_mat, 'data/BDG_FBX_100_contrast_DS_matrix_original_workflow.RDS')


# load local function files, as they have been updated
source('R/3_pathwayCrosstalk.R')
source('R/4_classification.R')
source('R/6_crosstalkInhibition.R')

#' @title
#' @description
#' @importFrom magrittr %>%
#' @param expression_matrix An expression matrix of intensity values for gene probes,
#' where the row names are the probe identifiers and the column names are the sample identifiers. Pre-processing should be completed,
#' including normalization, removing control sequences, log2 transformation, collapse to genes.
#' @param DEGs A vector of strings containing differentially expressed genes.
#' @param DEPs A vector of strings containing enriched pathways.
#' @param groups A dataframe containing the sample identifiers and associated treatment conditions.
#' @param pathways A list in which each element is a named list of genes corresponding to a particular pathway.
#' @export
network_wrapper <- function(expression_matrix, DEGs, DEPs, pathways, groups, ...){
        # ^ network characterization output can be added optionally?
        # characterization = FALSE,
        # this is annotation of the pathways with common names
        # may be better outside the wrapper, with the final results

    output_list <- c()

    # filter pathways to only include elements in DEPS
    # e_pathways <- pathways[names(pathways) %in% DEPs$pathway]
    # ^ moved to pathwayCrosstalk() function

    # step 3
    crosstalk_matrix <- pathwayCrosstalk(DEPs, expression_matrix,
                                         groups, pathways)

    # step 4
    network_results <- subtypeNetwork(crosstalk_matrix, groups) # documentation
    output_list[['full_network_results']] <- network_results

    # # step 5
    network <- network_results['graph']
    # if(characterization){
    #
    # }

    # step 6
    significant_crosstalks <- crosstalkInhibition(network)[[1]]
    # fix the format here^



    output_list[['significant_crosstalks']] <- significant_crosstalks



    # output
    # df of step 4 network and associated nNE values (PCI values)
    # df of step 6 network + ^

    return(output_list)

}



# save new matrix
# saveRDS(crosstalk_matrix, 'data/BDG_FBX_100_contrast_DS_matrix_wrapper_workflow.RDS')


# test new workflow

test <- network_wrapper(expression_matrix, DEGs, DEPs, pathways, groups)



# compare step 3 between new wrapper workflow and original version
old <- readRDS('data/BDG_FBX_100_contrast_DS_matrix_original_workflow.RDS')
new <- readRDS('data/BDG_FBX_100_contrast_DS_matrix_wrapper_workflow.RDS')

x <- (old[1,] + new[1,])/2
y <- (old[1,] - new[1,])

qplot(unlist(x),unlist(y))

dim(old) == dim(new)
all(colnames(old) == colnames(new))

any(old == new)
sum(old == new)

rowMeans(old) - rowMeans(new)
mean(rowMeans(old) - rowMeans(new))
mean(colMeans(old) - colMeans(new))
