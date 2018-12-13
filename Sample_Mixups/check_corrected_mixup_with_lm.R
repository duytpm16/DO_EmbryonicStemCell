### Helper lm function. Runs lm on all pQTL that pass threshold
#    1.) threshold to subset the lod.peaks data
lm_QTL <- 
  function(threshold){
    
    threshold <- as.vector(threshold)
    
    
    ### Perform for each threshold
    for(i in threshold){
      
      high_peaks <- lod.peaks %>%
        filter(lod >= i & cis == TRUE)
      
      
      
      # Creating matrix to store predicted values
      predicted_matrix           <- as.data.frame(matrix(0, nrow = length(dimnames(genoprobs[[1]])[[1]]), ncol = nrow(high_peaks)))
      rownames(predicted_matrix) <- dimnames(genoprobs[[1]])[[1]]
      column_names               <- character()
      
      
      
      
      
      ### Run regression on each pQTL
      for(j in 1:nrow(high_peaks)){
        
        # Extract pQTL info
        expr_column <- high_peaks$annot.id[j]
        qtl_chr     <- high_peaks$qtl.chr[j]
        qtl_marker  <- high_peaks$marker.id[j]
        
        
        
        
        
        # Get expression of one of the pQTL
        y <- expr[ ,expr_column]
        
        
        
        
        
        
        # Run lm on pQTL if no missing values are present
        if(all(!is.na(y))){
          column_names <- c(column_names, expr_column)               # Storing name of expression if it passes the if statement
          geno         <- genoprobs[[qtl_chr]][,,qtl_marker]         # Get the genotype probability at the pQTL in matrix form
          
          
          
          full_matrix  <- data.frame(cbind(pheno = y, geno, covar))          # Combine expression, genotype probability at pQTL, and covariates together for lm function
          
          
          # lm begins
          fit <- lm(pheno ~ ., data = full_matrix)                           # Fit lm          
          predicted_matrix[,j] <- predict(fit)
          
        } # if
        
      } # for j
      
      
      
      
      
      
      ### Rename columns of predicted matrix and get obsevered expression values in expr matrix
      colnames(predicted_matrix) <- column_names
      observed_matrix <- as.data.frame(expr[,column_names])
      
      
      stopifnot(rownames(predicted_matrix) == rownames(observed_matrix))
      stopifnot(colnames(predicted_matrix) == colnames(observed_matrix))
      
      
      
      ### Save predicted and observed matrix to global environment
      assign(paste0('predicted_matrix_',i), predicted_matrix, envir = .GlobalEnv)
      assign(paste0('observed_matrix_',i), observed_matrix, envir = .GlobalEnv)
      
    } # for i in threshold
  }
