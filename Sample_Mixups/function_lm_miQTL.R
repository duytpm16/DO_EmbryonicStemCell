### Helper lm function. Runs lm on all pQTL that pass threshold
#    1.) threshold to subset the lod.peaks data
lm_miQTL <- 
  function(threshold){
    
    
    ### Make threshold a vector to loop through
    threshold <- as.vector(threshold)
    
    
    
    ### For each threshold...
    for(i in threshold){
      
        high_peaks <- lod.peaks %>%
                                filter(lod >= i & cis == TRUE)
      
      
      
        # Creating matrix to store predicted values
        resid_matrix <- predicted_matrix <- as.data.frame(matrix(0, nrow = length(dimnames(genoprobs[[1]])[[1]]), ncol = nrow(high_peaks)))
        rownames(predicted_matrix) <- dimnames(genoprobs[[1]])[[1]]
      
      
      
      
      
      
        for(j in 1:nrow(high_peaks)) {
            # Extract pQTL info
            expr_column <- high_peaks$annot.id[j]
            qtl_chr     <- high_peaks$qtl.chr[j]
            qtl_marker  <- high_peaks$marker.id[j]  
        
        
            this_y <- expr[,expr_column]
            this_K <- K[[qtl_chr]]
            this_X_Q <- genoprobs[[qtl_chr]][,,qtl_marker]
            this_X_Q_2 <- cbind(rep(1, nrow(this_X_Q)), this_X_Q[,-which.max(colSums(this_X_Q))])
        
            colnames(this_X_Q_2)[1] <- "Intercept"
            
            this_X_covar <- dataset$covar
        
            if(j == 1) {
              null_fit <- lmmbygls(y = this_y, X = cbind(rep(1, nrow(this_X_covar)), this_X_covar), K = this_K)
            }
        
            alt_fit <- lmmbygls(y = this_y, X = cbind(this_X_Q_2, this_X_covar), K = this_K[names(this_y), names(this_y)])
        
            predicted_matrix[,j]  <- this_X_Q_2 %*% alt_fit$coefficients[colnames(this_X_Q_2)]
            resid_matrix[,j] <- this_y - this_X_covar %*% null_fit$coefficients[colnames(this_X_covar)]
        }
      
      
      
        # Get correlation matrices
        raw_cor <- zscale_cor <- rep(NA, nrow(expr))
        cor_mat <- matrix(0, nrow = nrow(expr), ncol = nrow(expr))
        for(k in 1:nrow(expr)) {
          these_cor     <- cor(t(rbind(resid_matrix[k,], predicted_matrix)))[-1]    # Correlation between same individual
          these_cor2    <- cor(t(rbind(resid_matrix[k,], predicted_matrix)))[,-1]   # Correlation between sample vs all
          raw_cor[k]    <- these_cor[k]
          cor_mat[k,]   <- these_cor2[1,]
          zscale_cor[k] <- scale(these_cor)[k]
          
        }
      
      
      
        assign(paste0('cor_mat_',i), cor_mat, envir = .GlobalEnv)
        assign(paste0('raw_cor_',i), raw_cor, envir = .GlobalEnv)
        assign(paste0('zscale_cor_',i), zscale_cor, envir = .GlobalEnv)
        
        
        
        
        ### Rename columns of predicted matrix and get obsevered expression values in expr matrix
        colnames(predicted_matrix) <- high_peaks$annot.id
        observed_matrix <- as.data.frame(expr[,high_peaks$annot.id])

        assign(paste0('predicted_matrix_',i), predicted_matrix, envir = .GlobalEnv)
        assign(paste0('observed_matrix_',i), observed_matrix, envir = .GlobalEnv)
      
    } # for i in threshold
}
