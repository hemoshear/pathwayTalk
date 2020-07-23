# 
# library(limma)
# library(GEOquery)
# library(dplyr)
# library(magrittr)
# library(affy)
# library(GiANT)
# library(WGCNA)
# 
# # browseVignettes('limma')
# # limmaUsersGuide()
# 
# # Read in microarray expression profiling data:
# 
# # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3151
# 
# # oncogene signature dataset (pathway deregulation markers)
# # RNA from human mammary epithelial cells expressing oncogenes
# # platform: Affymetrix Human Genome U133 Plus 2.0 Array
# 
# 
# # read in GEO ExpressionSet -----------------------------------------------
# 
# # GEO <- 'GSE3151'
# # # data <- getGEO(GEO)
# #
# # data <- data[[1]]
# # show(data)
# # class(data) # ExpressionSet
# #
# # # pData - retrieve information on experimental phenotypes in set
# # dim(pData(data))
# # head(pData(data))
# 
# # access GSE data tables if applicable
# # df <- getGSEDataTables(GEO) # none
# 
# # save dataset
# # saveRDS(data, file = 'onco_data.RDS')
# 
# data <- readRDS(file = 'onco_data.RDS')
# 
# 
# # inspect data ------------------------------------------------------------
# 
# # inspect data
# 
# experimentData(data)
# 
# colnames(data) # sample IDs
# table(data$organism_ch1) # homo sapiens
# table(data$type) # RNA
# table(data$label_ch1) # biotin
# 
# data$title # sample subtypes
# # GFP expression - control
# 
# data$channel_count # this is one-channel array data
# 
# head(data@featureData@data$`Gene Title`) # gene names
# 
# # trim sample id numbers from treatment groups:
# 
# data$title <- gsub('\\-\\d+', '', data$title)
# 
# data$title %<>% toupper()
# 
# table(data$title)
# 
# # qc / pre-processing  ----------------------------------------------------------------------
# 
# # first, log-transform intensity values:
# ?`AssayData-class`
# storageMode(data) # lockedEnvironment
# # because it's a locked environment, copy exprs and make changes, then re-assign
# 
# dim(data@assayData$exprs) # 54675    55
# 
# # create log2 version of intensity matrix:
# log2exp <- log2(data@assayData$exprs)
# 
# # reassign using exprs function:
# exprs(data) <- log2exp
# 
# 
# # array density plots - colored by treatment:
# 
# # unnormalized:
# # plotMD(data, main = "Unnormalized")
# # plotDensities(data,  main = "Unnormalized")
# 
# 
# # quantile normalization with normalizeBetweenArrays
# 
# # data.cor <- backgroundCorrect(data,method="minimum")
# # no background - correction already done for this data?
# 
# data.quant <- normalizeBetweenArrays(exprs(data), method='quantile')
# data.cl <- normalizeBetweenArrays(exprs(data), method='cyclicloess')
# 
# # plotMD(data.quant, main = "Quantile Normalized")
# # plotDensities(data.quant,  main = "Quantile Normalized")
# #
# # plotMD(data.cl, main = "Cyclic Loess Normalized")
# # plotDensities(data.cl,  main = "Cyclic Loess Normalized")
# 
# # data.quant and data.cl are equivalent format to exprs(data)
# 
# 
# 
# # filter out probes that don't appear to be expressed --------------------------
# 
# # how many arrays do we have here? does this need to be factored in?
# 
# exp_quant <- rowSums(data.quant > 0) > 7 # what threshold should we use here?
# exp_cl <- rowSums(data.cl > 0) > 7
# 
# # select probes we want to keep:
# 
# data.quant.exp <- data.quant[exp_quant, ]
# data.cl.exp <- data.cl[exp_cl, ]
# 
# dim(exprs(data))
# 
# dim(data.quant.exp) # same
# dim(data.cl.exp) # same
# 
# # do we need to filter out control sequences? how do we find these in this set?
# table(data@featureData@data$`Sequence Type`) # 62
# 
# controls <- data@featureData@data$`Sequence Type` == 'Control sequence'
# data[!controls,] # doesn't work
# 
# # gene annotation ---------------------------------------------------------
# 
# # move from probe level to gene level
# # bioconductor message boards - subscribe
# # wgcna
# 
# # AnnotationDb - virtual base class for all annotation packages (for Bioconductor)
# 
# # this array is 'Platform GPL570'
# # [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# 
# # annotation file: hgu133plus2.db
# library(hgu133plus2.db)
# 
# ?mapIds
# 
# data@featureData@data$ENTREZ_GENE_ID
# length(data@featureData@data$ENTREZ_GENE_ID) # 54675
# length(unique(data@featureData@data$ENTREZ_GENE_ID)) # 21880
# 
# data@featureData@data$`Gene Symbol`
# length(data@featureData@data$`Gene Symbol`) # 54675
# length(unique(data@featureData@data$`Gene Symbol`)) # 23521
# 
# # Entrez gene IDs and gene symbols already present
# 
# any(is.na(data@featureData@data$ENTREZ_GENE_ID)) # FALSE
# any(is.na(data@featureData@data$`Gene Symbol`)) # FALSE
# 
# # fit model with limma ----------------------------------------------------
# 
# # https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
# 
# data$title %<>% factor() %>% relevel(ref = 'GFP')
# 
# # generate design matrix
# design <- model.matrix(~ 0 + data$title)
# colnames(design) %<>% gsub('data\\$title', '', .)
# 
# # fit linear model
# fit <- lmFit(data, design)
# 
# colnames(fit$coefficients)
# 
# # fit with contrasts
# 
# contrast_mat <- makeContrasts(BCAT-GFP,
#                               E2F3-GFP,
#                               MYC-GFP,
#                               RAS-GFP,
#                               SRC-GFP,
#                               levels = design)
# 
# fit2 <- contrasts.fit(fit, contrast_mat)
# 
# fit <- eBayes(fit)
# fit2 <- eBayes(fit2)
# 
# # View(topTable(fit))
# # View(topTable(fit2))
# 
# # in fit2, coefs are lfcs between groups
# # p-values represent sig. of contrasts
# 
# table <- topTable(fit2, coef = 1, number = Inf)
# 
# 
# # decide fit
# results <- decideTests(fit)
# 
# # look at cyclic loess norm also
# # compare one row (probe set) for each type of normalization
# 
# # fits with quantile normalized data -----------------------------------------------
# 
# data_q <- data
# 
# exprs(data_q) <- data.quant
# 
# # generate design matrix
# design_q <- model.matrix(~ 0 + data_q$title)
# colnames(design_q) %<>% gsub('data_q\\$title', '', .)
# 
# # fit linear model
# fit_q <- lmFit(data_q, design_q)
# 
# colnames(fit_q$coefficients)
# 
# # fit with contrasts
# 
# contrast_mat_q <- makeContrasts(BCAT-GFP,
#                               E2F3-GFP,
#                               MYC-GFP,
#                               RAS-GFP,
#                               SRC-GFP,
#                               levels = design_q)
# 
# fit2_q <- contrasts.fit(fit_q, contrast_mat_q)
# 
# fit_q <- eBayes(fit_q)
# fit2_q <- eBayes(fit2_q)
# 
# # View(topTable(fit_q))
# # View(topTable(fit2_q))
# 
# table_q <- topTable(fit2_q, coef = 1, number = Inf)
# 
# # in fit2, coefs are lfcs between groups
# # p-values represent sig. of contrasts
# 
# # decide fit
# results <- decideTests(fit_q)
# 
# # gene IDs
# 
# any(is.na(topTable(fit_q)$ENTREZ_GENE_ID)) # but we do see missing values here
# any(is.na(topTable(fit2_q)$ENTREZ_GENE_ID))
# 
# # fits with cyclic loess normalized data ----------------------------------
# 
# data_cl <- data
# 
# exprs(data_cl) <- data.cl
# 
# # generate design matrix
# design_cl <- model.matrix(~ 0 + data_cl$title)
# colnames(design_cl) %<>% gsub('data_cl\\$title', '', .)
# 
# # fit linear model
# fit_cl <- lmFit(data_cl, design_cl)
# 
# colnames(fit_cl$coefficients)
# 
# # fit with contrasts
# 
# contrast_mat_cl <- makeContrasts(BCAT-GFP,
#                                 E2F3-GFP,
#                                 MYC-GFP,
#                                 RAS-GFP,
#                                 SRC-GFP,
#                                 levels = design_cl)
# 
# fit2_cl <- contrasts.fit(fit_cl, contrast_mat_cl)
# 
# fit_cl <- eBayes(fit_cl)
# fit2_cl <- eBayes(fit2_cl)
# 
# # View(topTable(fit_cl))
# # View(topTable(fit2_cl))
# 
# table_cl <- topTable(fit2_cl, coef = 1, number = Inf)
# 
# # in fit2, coefs are lfcs between groups
# # p-values represent sig. of contrasts
# 
# # gene IDs
# 
# any(is.na(topTable(fit_cl)$ENTREZ_GENE_ID))
# any(is.na(topTable(fit2_cl)$ENTREZ_GENE_ID))
# 
# # comparison between normalization methods --------------------------------
# 
# # top-ranked genes:
# 
# # topTable(fit)[,17:26]
# # topTable(fit_q)[,17:26]
# # topTable(fit_cl)[,17:26]
# #
# # topTable(fit2)[,17:25]
# # topTable(fit2_q)[,17:25]
# # topTable(fit2_cl)[,17:25]
# 
# # one probe-set:
# 
# fit$coefficients['200606_at',]
# fit_q$coefficients['200606_at',]
# fit_cl$coefficients['200606_at',]
# 
# fit2$coefficients['200606_at',]
# fit2_q$coefficients['200606_at',]
# fit2_cl$coefficients['200606_at',]
# 
# 
# # decideTest --------------------------------------------------------------
# 
# # unnormalized
# results <- decideTests(fit2)
# summary(results)
# vennDiagram(results)
# 
# # quantile normalized
# results_q <- decideTests(fit2_q)
# summary(results_q)
# vennDiagram(results_q)
# 
# # cyclic loess normalized
# results_cl <- decideTests(fit2_cl)
# summary(results_cl)
# vennDiagram(results_cl)
# 
# 
# 
# # strip charts  -----------------------------------------------------------
# 
# probe_id <- '1007_s_at'
# 
# # write function to generate strip charts:
# generateStripCharts <- function(probe_id){
# 
#      # select associated data from each normalized dataset
#     un_data <- data.frame(exprs(data)[probe_id,])
#     q_data <- data.frame(exprs(data)[probe_id,])
#     cl_data <- data.frame(exprs(data)[probe_id,])
# 
#     expression_sets <- c(un_data, q_data, cl_data)
#     methods <- c('unnormalized', 'quantile', 'cyclic_loess')
# 
#     # annotate with sample condition
# 
#     for (i in 1:length(expression_sets)){
# 
#         set <- expression_sets[[i]]
# 
#         set$sample <- rownames(set)
# 
#         set %<>% bind_rows
# 
#     }
# 
#     combined_data <- bind_cols(expression_sets)
#     combined_data$condition <- data@phenoData@data$title
# }
# 
# 
# test_probe <-  '1007_s_at'
# 
# generateStripCharts(test_probe)
# 
# 
# 
# 
# 
# 
# 
# 
# # combine probe sets --> gene level ---------------------------------------
# 
# # name rows of expression matrix with gene names, then collapse rows
# 
# # options for collapsing rows:
# 
# # mergeProbesForGenes from GiANT
#     # options - mean, max, min, median
#     # * no option for a weighted average *
#     # works on expression matrix where row names are gene identifiers
# 
# expression <- exprs(data)
# rownames(expression) <- data@featureData@data$ENTREZ_GENE_ID
# 
# length(unique(names(expression))) # 21881
# length(unique(data@featureData@data$ENTREZ_GENE_ID)) # 21880
# 
# expression1 <- mergeProbesForGenes(expression, method = c('median'))
# nrow(expression1) # 21880
# 
# # collapseRows from WGCNA
#     # options - Average, MaxMean, MinMean, maxRowVariance, absMaxMean, absMinMean, ME (eigenrow)
#     # paper from developers at UCLA - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3166942/
# 
# # issues with WGCNA install - in terminal, run ('xcode-select --install')
# 
# expression <- exprs(data)
# expression2 <- WGCNA::collapseRows()
# 
# 
# # comparing DE across methods ---------------------------------------------
# 
# data_genes <- data
# exprs(data_genes) <- expression1
# 
# # how to incorporate new expression matrix back into the ExpressionSet object?
# 
# # create new eset object?
# 
# exprs <- expression1
# 
# data_genes <- ExpressionSet(assayData = exprs)
# data_genes@phenoData <- data@phenoData
# 
# 
# # repeat analysis at gene level --------------------------------------------
# 
# data_genes_q <- normalizeBetweenArrays(exprs(data_genes), method='quantile')
# data_genes_cl <- normalizeBetweenArrays(exprs(data_genes), method='cyclicloess')
# 
# # fit model with limma ----------------------------------------------------
# 
# # https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
# 
# data <- data_genes
# 
# data$title %<>% factor() %>% relevel(ref = 'GFP')
# 
# # generate design matrix
# design <- model.matrix(~ 0 + data$title)
# colnames(design) %<>% gsub('data\\$title', '', .)
# 
# # fit linear model
# fit <- lmFit(data, design)
# 
# colnames(fit$coefficients)
# 
# # fit with contrasts
# 
# contrast_mat <- makeContrasts(BCAT-GFP,
#                               E2F3-GFP,
#                               MYC-GFP,
#                               RAS-GFP,
#                               SRC-GFP,
#                               levels = design)
# 
# fit2 <- contrasts.fit(fit, contrast_mat)
# 
# fit <- eBayes(fit)
# fit2 <- eBayes(fit2)
# 
# # View(topTable(fit))
# # View(topTable(fit2))
# 
# # in fit2, coefs are lfcs between groups
# # p-values represent sig. of contrasts
# 
# table <- topTable(fit2, coef = 1, number = Inf)
# 
# 
# # decide fit
# results <- decideTests(fit)
# 
# # look at cyclic loess norm also
# # compare one row (probe set) for each type of normalization
# 
# # fits with quantile normalized data -----------------------------------------------
# 
# data_q <- data_genes_q
# 
# # generate design matrix
# design_q <- model.matrix(~ 0 + data_q$title)
# colnames(design_q) %<>% gsub('data_q\\$title', '', .)
# 
# # fit linear model
# fit_q <- lmFit(data_q, design_q)
# 
# colnames(fit_q$coefficients)
# 
# # fit with contrasts
# 
# contrast_mat_q <- makeContrasts(BCAT-GFP,
#                                 E2F3-GFP,
#                                 MYC-GFP,
#                                 RAS-GFP,
#                                 SRC-GFP,
#                                 levels = design_q)
# 
# fit2_q <- contrasts.fit(fit_q, contrast_mat_q)
# 
# fit_q <- eBayes(fit_q)
# fit2_q <- eBayes(fit2_q)
# 
# # View(topTable(fit_q))
# # View(topTable(fit2_q))
# 
# table_q <- topTable(fit2_q, coef = 1, number = Inf)
# 
# # in fit2, coefs are lfcs between groups
# # p-values represent sig. of contrasts
# 
# # decide fit
# results <- decideTests(fit_q)
# 
# # gene IDs
# 
# any(is.na(topTable(fit_q)$ENTREZ_GENE_ID)) # but we do see missing values here
# any(is.na(topTable(fit2_q)$ENTREZ_GENE_ID))
# 
# # fits with cyclic loess normalized data ----------------------------------
# 
# data_cl <- data_genes_cl
# 
# # generate design matrix
# design_cl <- model.matrix(~ 0 + data_cl$title)
# colnames(design_cl) %<>% gsub('data_cl\\$title', '', .)
# 
# # fit linear model
# fit_cl <- lmFit(data_cl, design_cl)
# 
# colnames(fit_cl$coefficients)
# 
# # fit with contrasts
# 
# contrast_mat_cl <- makeContrasts(BCAT-GFP,
#                                  E2F3-GFP,
#                                  MYC-GFP,
#                                  RAS-GFP,
#                                  SRC-GFP,
#                                  levels = design_cl)
# 
# fit2_cl <- contrasts.fit(fit_cl, contrast_mat_cl)
# 
# fit_cl <- eBayes(fit_cl)
# fit2_cl <- eBayes(fit2_cl)
# 
# # View(topTable(fit_cl))
# # View(topTable(fit2_cl))
# 
# table_cl <- topTable(fit2_cl, coef = 1, number = Inf)
# 
# # in fit2, coefs are lfcs between groups
# # p-values represent sig. of contrasts
# 
# # gene IDs
# 
# any(is.na(topTable(fit_cl)$ENTREZ_GENE_ID))
# any(is.na(topTable(fit2_cl)$ENTREZ_GENE_ID))
# 
# # comparison between normalization methods --------------------------------
# 
# # top-ranked genes:
# 
# # topTable(fit)[,17:26]
# # topTable(fit_q)[,17:26]
# # topTable(fit_cl)[,17:26]
# #
# # topTable(fit2)[,17:25]
# # topTable(fit2_q)[,17:25]
# # topTable(fit2_cl)[,17:25]
# 
# # one probe-set:
# 
# fit$coefficients['200606_at',]
# fit_q$coefficients['200606_at',]
# fit_cl$coefficients['200606_at',]
# 
# fit2$coefficients['200606_at',]
# fit2_q$coefficients['200606_at',]
# fit2_cl$coefficients['200606_at',]
# 
# 
# # decideTest --------------------------------------------------------------
# 
# # unnormalized
# results <- decideTests(fit2)
# summary(results)
# vennDiagram(results)
# 
# # quantile normalized
# results_q <- decideTests(fit2_q)
# summary(results_q)
# vennDiagram(results_q)
# 
# # cyclic loess normalized
# results_cl <- decideTests(fit2_cl)
# summary(results_cl)
# vennDiagram(results_cl)
# 

# venn diagrams / tables
# comparing the normalization methods vs unnormalized for each contrast

