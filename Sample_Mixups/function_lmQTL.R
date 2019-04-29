### Helper lm function. Runs lm on all pQTL that pass threshold
#    1.) threshold to subset the lod.peaks data
lm_QTL <- 
  function(threshold, id, genoprobs, predict_genoprobs){
    
    threshold <- as.vector(threshold)
    
    
    ### For each threshold...
    for(i in threshold){
      
      
        # Subset QTLs above threshold i
        high_peaks <- lod.peaks %>%
          filter(lod >= i & cis == TRUE)
      
      
      
        # Creating matrix to store predicted values
        predicted_matrix           <- as.data.frame(matrix(0, nrow = length(dimnames(predict_genoprobs[[1]])[[1]]), 
                                                           ncol = nrow(high_peaks), 
                                                           dimnames = list(dimnames(predict_genoprobs[[1]])[[1]], high_peaks[,id][[1]])))
      
      
      
      
      
        ### For each QTL
        for(j in 1:nrow(high_peaks)){
          
            # Extract QTL info
            expr_column <- high_peaks[,id][[1]][j]
            qtl_chr     <- high_peaks$qtl.chr[j]
            qtl_marker  <- high_peaks$marker.id[j]
        
        
        
        
        
            # Get expression of one of the pQTL
            y <- expr[ ,expr_column]
            y <- y[!is.na(y)]
            geno <- genoprobs[[qtl_chr]][names(y),,qtl_marker]         # Get the genotype probability at the pQTL in matrix form
            stopifnot(names(y) == rownames(geno))
            full_matrix  <- data.frame(cbind(pheno = y, geno))  # Combine expression, genotype probability at pQTL, and covariates together for lm function
            
            
            predict_geno <- predict_genoprobs[[qtl_chr]][,,qtl_marker]
            predict_this <- data.frame(Intercept = 1, predict_geno)
            colnames(predict_this)[1] <- '(Intercept)'
            predict_this <- predict_this[,-9]
           
        
        
          # lm begins
          fit <- lm(pheno ~ ., data = full_matrix)                   # Fit lm
  
          predicted_matrix[,j] <- as.matrix(predict_this) %*% coefficients(fit)[colnames(predict_this)] # Get predicted value
          
      
    } # for j
    
    
    
    
    
    
    ### Rename columns of predicted matrix and get obsevered expression values in expr matrix
    observed_matrix <- as.data.frame(expr[,high_peaks[,id][[1]]])
    


    
    ### Save predicted and observed matrix to global environment
    assign(paste0('predicted_matrix_',i), predicted_matrix, envir = .GlobalEnv)
    assign(paste0('observed_matrix_',i), observed_matrix, envir = .GlobalEnv)
    
  } # for i in threshold
}
