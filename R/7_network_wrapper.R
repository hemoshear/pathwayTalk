
# pathways <- as.list(reactome.db::reactomePATHID2EXTID)
# pathways <- pathways[grep('HSA', names(pathways))]
# saveRDS(pathways, 'data/wrapper_dev_data/reactome_pathways.RDS')

# import dev materials
expression_matrix <- readRDS('data/wrapper_dev_data/processed_BDG_expression_mat.RDS')
groups <- readRDS('data/wrapper_dev_data/processed_BDG_groups.RDS')
DEGs <- readRDS('data/wrapper_dev_data/processed_BDG_DEGs.RDS')
DEPs <- readRDS('data/wrapper_dev_data/processed_BDG_EPs.RDS')
pathways <- readRDS('data/wrapper_dev_data/reactome_pathways.RDS')

# subset to one contrast
DEPs %<>% filter(contrast == 'FBX_100-VEH_DMSO_ICS')
DEGs <- DEGs[['FBX_100-VEH_DMSO_ICS']]
groups %<>% filter(group %in% c('VEH_DMSO_ICS', 'FBX_100'))
expression_matrix <- expression_matrix[,c(groups$ea_code)]

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
network_wrapper <- function(expression_matrix, DEGs, DEPs, pathways, groups,
                            processes = 4, characterization = FALSE, ...){
    # ^ network characterization output can be added optionally?

    # filter pathways to only include elements in DEPS
    e_pathways <- pathways[names(pathways) %in% DEPs$pathway] # correct length - yes

    # step 3
    crosstalk_matrix <- pathwayCrosstalk(DEPs, pathways, groups, processes = processes)
    # pathwayCrosstalk - add pathways as argument ^

    # step 4
    network <- subtypeNetwork(crosstalk_matrix, groups, ...) # documentation

    # step 6
    ne_network <- crosstalkInhibition(list(network))

    # output
    # df of step 4 network and associated nNE values (PCI values)
    # df of step 6 network + ^



}
