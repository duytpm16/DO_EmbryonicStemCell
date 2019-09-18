### Helper lm function. Runs lm on all pQTL that pass threshold
#    1.) threshold to subset the lod.peaks data
lm_miQTL <- 
  function(threshold, lod.peaks, genoprobs, expr, K, covar, covar_in_pred = FALSE){
    
    ### Make threshold parameter a vector to loop over
    threshold <- as.vector(threshold)
    
    
    
    ### For each threshold...
    for(i in threshold){
      
      
      
        # Subset lod.peaks to contain QTLs above threshold i and make sure they are cis
        high_peaks <- lod.peaks %>%
                                filter(lod >= i & cis == TRUE)
      
      
      
      
      
        # Creating matrix to store predicted values
        predicted_matrix           <- as.data.frame(matrix(0, nrow = length(dimnames(genoprobs[[1]])[[1]]), ncol = nrow(high_peaks),
                                                           dimnames = list(dimnames(genoprobs[[1]])[[1]], high_peaks$annot.id)))
      
      
      
      
      
        ### For each QTL...
        for(j in 1:nrow(high_peaks)) {
          
            # Extract QTL info
            expr_column <- high_peaks$annot.id[j]
            qtl_chr     <- high_peaks$qtl.chr[j]
            qtl_marker  <- high_peaks$marker.id[j]  
            
            
            
            # Get required data for model fitting and prediction
            this_y <- expr[,expr_column]
            this_y <- this_y[!is.na(this_y)]                # Observed expressions
            this_K <- K[[qtl_chr]]
            this_K <- this_K[names(this_y),names(this_y)]   # Kinship matrix at QTL chromosome
            this_X_Q   <- genoprobs[[qtl_chr]][,,qtl_marker]                                            # Genoprob at QTL chromosome and marker 
            this_X_Q_2 <- cbind(rep(1, nrow(this_X_Q)), this_X_Q[,-which.max(colSums(this_X_Q))])     # Combine intercept to genoprob, used for model fitting
            this_X_Q_2 <- this_X_Q_2[names(this_y),]
            this_predict <- genoprobs[[qtl_chr]][,,qtl_marker]                                        # Get same probabilities in full array at QTL chromosome
            if(covar_in_pred){
               this_predict_2 <-  cbind(rep(1, nrow(this_predict)), this_predict[,-which.max(colSums(this_X_Q))], covar)  # Add intercept and covariates to probabilities used for prediction
            }else{
               this_predict_2 <-  cbind(rep(1, nrow(this_predict)), this_predict[,-which.max(colSums(this_X_Q))])         # Add intercept to probabilities used for prediction
            }
            
            
            
            # Rename first column of probabilities to Intercept
            colnames(this_X_Q_2)[1] <- "Intercept"
            colnames(this_predict_2)[1] <- "Intercept"
            
            
            # Get covariates
            this_X_covar <- covar[names(this_y),]
            
            stopifnot(rownames(this_y) == rownames(this_K))
            stopifnot(rownames(this_y) == rownames(cbind(this_X_Q_2, this_X_covar)))
            
            
            
            # Fit Alt model
            alt_fit <- lmmbygls(y = this_y, X = cbind(this_X_Q_2, this_X_covar), K = this_K, na.action = na.omit)
            
            
            # Get predicted value 
            predicted_matrix[,j]  <- this_predict_2 %*% alt_fit$coefficients[colnames(this_predict_2)]
          
        } # for j
      
      
      
      
      
      
      ### Rename columns of predicted matrix and get obsevered expression values in expr matrix
      observed_matrix <- as.data.frame(expr[,high_peaks$annot.id])
      
      
      
      ### Save predicted and observed matrix to global environment
      assign(paste0('predicted_matrix_',i), predicted_matrix, envir = .GlobalEnv)
      assign(paste0('observed_matrix_',i), observed_matrix, envir = .GlobalEnv)
      
      
      
    } # for i in threshold
}