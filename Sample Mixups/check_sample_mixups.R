options(stringsAsFactors = FALSE)
library(dplyr)
library(corrplot)
library(pdist)
library(tibble)
library(miqtl)       #devtools::install_github('gkeele/miqtl')
#source('function_lmQTL.R')





### Load in data
#load("DO_mESC/do_mESC_proteomics_qtl_viewer_v2.2.RData")                 # QTL Viewer formatted data
dataset <- get('dataset.islet.proteins')                      







### Get required data
expr      <- dataset[['rankz']]
covar     <- dataset$covar
lod.peaks <- dataset$lod.peaks$additive
lod.peaks$qtl.chr <- as.character(lod.peaks$qtl.chr)





### Making sure names match
for(i in 2:length(genoprobs)){
    stopifnot(dimnames(genoprobs[[i]])[[1]] == dimnames(genoprobs[[1]])[[1]])
    stopifnot(dimnames(K[[i]])[[1]] == dimnames(K[[1]])[[1]])
}




### Get predicted and observed matrix for threshold
lm_miQTL(threshold = c(10, 15, 20), 
         genoprobs = genoprobs, 
         expr = expr, 
         covar = covar, 
         K = K,
         lod.peaks = lod.peaks,
         covar_in_pred = FALSE)







### Get correlations between predicted vs observed matrix
cor_mat_10  <- cor(t(observed_matrix_10), t(predicted_matrix_10), use = "pairwise.complete.obs")
cor_mat_15  <- cor(t(observed_matrix_15), t(predicted_matrix_15), use = "pairwise.complete.obs")
cor_mat_20  <- cor(t(observed_matrix_20), t(predicted_matrix_20), use = "pairwise.complete.obs")





#   Get sample name that is highly correlated to sample0
cor_max_10 <- colnames(cor_mat_10)[apply(cor_mat_10, 1, which.max)]
cor_max_15 <- colnames(cor_mat_15)[apply(cor_mat_15, 1, which.max)]
cor_max_20 <- colnames(cor_mat_20)[apply(cor_mat_20, 1, which.max)]







### Combine the best correlated samples across 3 thresholds
stopifnot(rownames(cor_mat_10) == rownames(observed_matrix_10) & rownames(cor_mat_10) == rownames(observed_matrix_15) &  rownames(cor_mat_10) == rownames(observed_matrix_20))
stopifnot(rownames(cor_mat_15) == rownames(observed_matrix_10) & rownames(cor_mat_15) == rownames(observed_matrix_15) &  rownames(cor_mat_15) == rownames(observed_matrix_20))
stopifnot(rownames(cor_mat_20) == rownames(observed_matrix_10) & rownames(cor_mat_20) == rownames(observed_matrix_15) &  rownames(cor_mat_20) == rownames(observed_matrix_20))

best_correlation_df <- data.frame(thres_10 = cor_max_10, thres_15 = cor_max_15, thres_20 = cor_max_20)
rownames(best_correlation_df) <- rownames(cor_mat_10)






#   See how many samples are highly correlated with one sample across threshold
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
multi_mixup <- best_correlation_df %>%
                                   rownames_to_column('sample') %>%
                                   filter(num_of_correlated_samples !=1) %>%                           # Get rows where highly correlated samples is not unique across threshold
                                   column_to_rownames('sample')

