### Helper lm function. Runs lm on all pQTL that pass threshold
#    1.) threshold to subset the lod.peaks data
lm_QTL <- 
  function(threshold){
   
    ### Make threshold parameter a vector to loop over
    threshold <- as.vector(threshold)
    
    
    
    ### For each threshold...
    for(i in threshold){
      
      
      
        # Subset lod.peaks to contain QTLs above threshold i and make sure they are cis
        high_peaks <- lod.peaks %>%
                                filter(lod >= i & cis == TRUE)
      
      
        
      
      
        # Creating matrix to store predicted values
        predicted_matrix           <- as.data.frame(matrix(0, nrow = length(dimnames(probs[[1]])[[1]]), ncol = nrow(high_peaks)))
        rownames(predicted_matrix) <- dimnames(probs[[1]])[[1]]

      
      
      
      
        ### For each QTL...
        for(j in 1:nrow(high_peaks)) {
            
            # Extract QTL info
            expr_column <- high_peaks$annot.id[j]
            qtl_chr     <- high_peaks$qtl.chr[j]
            qtl_marker  <- high_peaks$marker.id[j]  
            
            
          
            # Get required data for model fitting and prediction
            this_y <- expr[,expr_column]                                                                        # Observed expressions
            this_K <- K[[qtl_chr]]                                                                              # Kinship matrix as QTL chromosome
            this_X_Q <- genoprobs[[qtl_chr]][,,qtl_marker]                                                      # Genoprob at QTL chromosome and marker 195 x 8
            this_X_Q_2 <- cbind(rep(1, nrow(this_X_Q)), this_X_Q[,-which.max(colSums(this_X_Q))])               # Combine intercept to genoprob, used for model fitting
            this_predict <- probs[[qtl_chr]][,,qtl_marker]                                                      # Get same probabilities in full MUGA array at QTL chromosome and marker 969 x 8
            this_predict_2 <-  cbind(rep(1, nrow(this_predict)), this_predict[,-which.max(colSums(this_X_Q))])  # Add intercept to probabilities used for prediction
          
          
            # Rename first column of probabilities to Intercept
            colnames(this_X_Q_2)[1] <- "Intercept"
            colnames(this_predict_2)[1] <- "Intercept"
          
          
            # Get covariates
            this_X_covar <- dataset$covar

          
          
            # Fit Null model
            if(j == 1) {
               null_fit <- lmmbygls(y = this_y, X = cbind(rep(1, nrow(this_X_covar)), this_X_covar), K = this_K)
            }
          
            # Fit Alt model
            alt_fit <- lmmbygls(y = this_y, X = cbind(this_X_Q_2, this_X_covar), K = this_K[names(this_y), names(this_y)])
          
          
            # Get predicted value 
            predicted_matrix[,j]  <- this_predict_2 %*% alt_fit$coefficients[colnames(this_predict_2)]
     
        }
        
      
      
      
      
      
      ### Rename columns of predicted matrix and get obsevered expression values in expr matrix
      colnames(predicted_matrix) <- high_peaks$annot.id
      observed_matrix <- as.data.frame(expr[,high_peaks$annot.id])
      
       
      
      ### Save predicted and observed matrix to global environment
      assign(paste0('predicted_matrix_',i), predicted_matrix, envir = .GlobalEnv)
      assign(paste0('observed_matrix_',i), observed_matrix, envir = .GlobalEnv)
      
    } # for i in threshold
}
