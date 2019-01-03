options(stringsAsFactors = FALSE)
library(dplyr)
library(corrplot)
library(pdist)
library(tibble)
library(openxlsx)
library(qtl2convert)
source('function_lm_MUGA_miQTL.R')




### Load in data
load('final_dataset_mESC_v2.RData')                 # QTL Viewer formatted data
load("~/Desktop/DO_mESC/Proteomics/Original Data/Choi_probs_interp69k.RData")                          # Load in gigaMuga array
pbid    <- read.xlsx('~/Desktop/DO_mESC/Proteomics/Original Data/PBID_SampleID_2017_05_01.xlsx')       # Load in sample annotation for gigaMuga array
dataset <- get('dataset.esc.proteins')                      





### Get required data
lod.peaks <- dataset.esc.proteins$lod.peaks$additive
lod.peaks$qtl.chr <- as.character(lod.peaks$qtl.chr)
expr      <- dataset[['rankz']][dimnames(genoprobs[[1]])[[1]],]
covar     <- dataset$covar[dimnames(genoprobs[[1]])[[1]],]
probs     <- probs_doqtl_to_qtl2(probs, map = markers)







### Just making sure name match
#
#   For genoprobs 195 x 8 x 69005
#   For probs (gigaMuga) 969 x 8 x 690005
pbid$SampleID <- gsub('-','.',pbid$SampleID)


for(i in 2:length(genoprobs)){
    stopifnot(dimnames(genoprobs[[i]])[[1]] == dimnames(genoprobs[[1]])[[1]])
}
stopifnot(rownames(expr) == rownames(covar), rownames(covar) == rownames(genoprobs[[1]][[1]]), rownames(expr) == rownames(genoprobs[[1]][[1]]))


for(i in 2:length(probs)){
  stopifnot(dimnames(probs[[i]])[[1]] == dimnames(probs[[1]])[[1]])
  
}









### Get predicted and observed matrix for threshold
lm_MUGA_miQTL(c(10, 15, 20, 25))











### Get correlations between predicted vs observed matrix
cor_mat_10  <- cor(t(observed_matrix_10), t(predicted_matrix_10), use = 'pairwise.complete.obs')
cor_mat_15  <- cor(t(observed_matrix_15), t(predicted_matrix_15), use = 'pairwise.complete.obs')
cor_mat_20  <- cor(t(observed_matrix_20), t(predicted_matrix_20), use = 'pairwise.complete.obs')
cor_mat_25  <- cor(t(observed_matrix_25), t(predicted_matrix_25), use = 'pairwise.complete.obs')

cor_matrix_column_names <- gsub('X', 'PB', colnames(cor_mat_10))
which.sampleID <- which(cor_matrix_column_names %in% pbid$SampleID)
cor_matrix_column_names[which.sampleID] <- pbid$PBID[match(cor_matrix_column_names[which.sampleID], pbid$SampleID)]
colnames(cor_mat_10) <- cor_matrix_column_names
colnames(cor_mat_15) <- cor_matrix_column_names
colnames(cor_mat_20) <- cor_matrix_column_names









#   Get sample name that is highly correlated to sample0
cor_max_10 <- colnames(cor_mat_10)[apply(cor_mat_10, 1, which.max)]
cor_max_15 <- colnames(cor_mat_15)[apply(cor_mat_15, 1, which.max)]
cor_max_20 <- colnames(cor_mat_20)[apply(cor_mat_20, 1, which.max)]
cor_max_25 <- colnames(cor_mat_25)[apply(cor_mat_25, 1, which.max)]













### Combine the best correlated samples across 3 thresholds
stopifnot(rownames(cor_mat_10) == rownames(observed_matrix_10), rownames(cor_mat_10) == rownames(observed_matrix_15),  rownames(cor_mat_10) == rownames(observed_matrix_20))
stopifnot(rownames(cor_mat_15) == rownames(observed_matrix_10), rownames(cor_mat_15) == rownames(observed_matrix_15),  rownames(cor_mat_15) == rownames(observed_matrix_20))
stopifnot(rownames(cor_mat_20) == rownames(observed_matrix_10), rownames(cor_mat_20) == rownames(observed_matrix_15),  rownames(cor_mat_20) == rownames(observed_matrix_20))
best_correlation_df <- data.frame(thres_10 = cor_max_10, 
                                  thres_15 = cor_max_15, 
                                  thres_20 = cor_max_20)
rownames(best_correlation_df) <- rownames(cor_mat_10)







### See how many samples are highly correlated with one sample across threshold
num_of_correlated_samples <- apply(best_correlation_df, 1, FUN = function(x) length(unique(x)))













### Find samples that are not mix-ups, mix-ups, and weird mixups
no_mixup    <- best_correlation_df %>%
                                   rownames_to_column('sample') %>%
                                   filter(num_of_correlated_samples == 1 & sample == thres_10) %>%     # Only check column thres_10 since the rows should be unique if num of correlated sample is 1
                                   column_to_rownames('sample') 
mixup       <- best_correlation_df %>%
                                   rownames_to_column('sample') %>%
                                   filter(num_of_correlated_samples == 1 & sample != thres_10) %>%     # Only check column thres_10 since the rows should be unique if num of correlated sample is 1
                                   column_to_rownames('sample')
weird_mixup <- best_correlation_df %>%
                                   rownames_to_column('sample') %>%
                                   filter(num_of_correlated_samples !=1) %>%                                               # Get rows where highly correlated samples is not unique across threshold
                                   column_to_rownames('sample')



mixup_samples_not_in_genoprobs <- mixup %>%
                                        rownames_to_column('sample') %>%
                                        filter(!(paste0(thres_10, '_rep1') %in% dimnames(genoprobs[[1]])[[1]])) %>%
                                        column_to_rownames('sample') 
mixup_samples_in_genoprobs     <- mixup %>%
                                        rownames_to_column('sample') %>%
                                        filter(paste0(thres_10, '_rep1') %in% dimnames(genoprobs[[1]])[[1]]) %>%
                                        column_to_rownames('sample') 






### Swapping Genotypes
#   Only these are replicates in original genotype: 'PB358.02_rep2', 'PB358.05_rep2', 'PB360.07_rep2'

genoprobs2 <- probs_qtl2_to_doqtl(genoprobs)
genoprobs3 <- probs_qtl2_to_doqtl(genoprobs)
probs2     <- probs_qtl2_to_doqtl(probs)


swap_samples_in_genoprobs <- c('PB361.76_rep1','PB361.77_rep1','PB361.81_rep1','PB361.88_rep1','PB358.04_rep1','PB357.22_rep1')
swap_to_this_in_genoprobs <- c('PB361.88_rep1','PB361.76_rep1','PB361.77_rep1','PB361.81_rep1','PB357.22_rep1','PB358.04_rep1')
dimnames(genoprobs2)[[1]][match(swap_samples_in_genoprobs,dimnames(genoprobs2)[[1]])] <- swap_to_this_in_genoprobs


genoprobs2['PB359.22_rep1',,] <- genoprobs2['PB359.33_rep1',,]
genoprobs2['PB358.05_rep1',,] <- probs2['PB358.50',,]
genoprobs2['PB359.06_rep1',,] <- probs2['A903DI1B.A09',,]
genoprobs2['PB359.07_rep1',,] <- probs2['A903DI1B.F03',,]
genoprobs2['PB359.09_rep1',,] <- probs2['PB359.90',,]
genoprobs2['PB359.33_rep1',,] <- probs2['A903DI1B.B05',,]
genoprobs2['PB360.05_rep1',,] <- probs2['PB360.50',,]
genoprobs2['PB360.07_rep2',,] <- probs2['A903DI1C.E09',,]
genoprobs2['PB361.07_rep1',,] <- probs2['A903DI12.E07',,]
genoprobs2['PB361.09_rep1',,] <- probs2['A903DI12.G03',,]
genoprobs2['PB362.04_rep1',,] <- probs2['PB362.40',,]
genoprobs2['PB366.72_rep1',,] <- probs2['A903DI12.E09',,]
genoprobs2['PB359.25_rep1',,] <- probs2['PB359.25',,]     ### A case where it's more correlated to a duplicated one







### Save data
genoprobs <- probs_doqtl_to_qtl2(genoprobs2, map = markers, pos_column = 'pos')
K <- calc_kinship(genoprobs, type = 'loco, cores = 0)
save(dataset.esc.proteins, K, map, markers, genoprobs, file = 'final_dataset_mESC_v2.RData')
