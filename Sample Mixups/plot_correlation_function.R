### Helper plot function
#   1.) name1        = sample name to plot
#   2.) name2        = second sample name to plot. Used to see the mix-up between name1 and name2 across matrices in (3)
#   3.) cor_matrices = vector of characters corresponding to the name of the correlation matrices
#   4.) threshold    = character of threshold used to subset pQTL
plot_mix_up_cor <- 
  function(name1, name2 = NULL, cor_matrices = c('cor_mat_10', 'cor_mat_15', 'cor_mat_20'), threshold = c('10','15','20')){
    
    cor_matrices <- as.vector(cor_matrices)
    
    # Plot display
    if(!is.null(name2)){
      par(mfrow = c(2,length(cor_matrices)))
    }else{
      par(mfrow = c(1,length(cor_matrices)))
    }
    
    
    
    # Assign matrices to temporary variable p*, where * = 1..length cor matrices
    for(i in 1:length(cor_matrices)){
      assign(paste0('p',i), get(cor_matrices[i]))
    }
    
    
    
    #### Plot results for name1 across each correlation matrix
    for(i in 1:length(cor_matrices)){
      
      df <- get(paste0('p',i))
      
      
      
      ### Plot all correlations for name1
      plot(df[name1,], main = paste(name1, 'LOD >=', threshold[i]), cex.main = 2.5, cex.axis = 2, cex.lab = 1.7, ylab = 'Correlation Value', xlab = 'Samples', col = 'grey57', ylim = c(-.8,.95))
      abline(h = 0.5)
      
      index_max_cor_1 <- which.max(df[name1,])
      points(index_max_cor_1, df[name1, index_max_cor_1], pch=21, cex = 2, bg='darkred')   
      
      ### Plotting maximum correlation point
      if(index_max_cor_1 >= 100 | index_max_cor_1 <= 750){
        if(!is.null(name2) | length(cor_matrices == 1)){
          text(index_max_cor_1, df[name1,index_max_cor_1] -.1 , labels = names(which.max(df[name1,])), cex = 2, col = 'darkred') 
        }else{
          text(index_max_cor_1, df[name1,index_max_cor_1] -.05 , labels = names(which.max(df[name1,])), cex = 2, col = 'darkred') 
        }
      }else if(index_max_cor_1 < 100){
        text(index_max_cor_1 + 155, df[name1, index_max_cor_1], labels = names(which.max(df[name1,])), cex = 2, col = 'darkred')  
      }else{                                           
        text(index_max_cor_1 - 155, df[name1, index_max_cor_1], labels = names(which.max(df[name1,])), cex = 2, col = 'darkred')  
      }
      
      
      index_name1 <- which(colnames(df) == gsub('_rep[12]','',name1))
      points(index_name1, df[name1, gsub('_rep[12]', '', name1)], pch = 21, cex = 2, bg = 'darkblue')        
      if(index_name1 >= 100 & index_name1 <= 750){       
        if(!is.null(name2) | length(cor_matrices == 1)){
          text(index_name1,   df[name1, gsub('_rep[12]', '', name1)] -.1, labels = gsub('_rep[12]','',name1), cex = 2, col = 'darkblue', font = 2)
        }else{
          text(index_name1,   df[name1, gsub('_rep[12]', '', name1)] -.05, labels = gsub('_rep[12]','',name1), cex = 2, col = 'darkblue', font = 2)
        }
      }else if(index_name1 < 100){                          
        text(index_name1 + 190, df[name1, gsub('_rep[12]', '', name1)], labels = gsub('_rep[12]','',name1), cex = 2, col = 'darkblue', font = 2)
      }else{                            
        text(index_name1 - 190, df[name1, gsub('_rep[12]', '', name1)], labels = gsub('_rep[12]','',name1), cex = 2, col = 'darkblue', font = 2)
      }
      
      
      
    }
    
    # Plot results for name2 across each correlation matrix if it exists. Good for checking sample mix-ups
    if(!is.null(name2)){
      
      for(i in 1:length(cor_matrices)){
        
        df <- get(paste0('p',i))
        plot(df[name2,], main = paste(name2, 'LOD >=', threshold[i]), cex.main = 2.5, cex.axis = 2, cex.lab = 1.7, ylab = 'Correlation Value', xlab = 'Samples', col = 'grey57', ylim = c(-.8,.95))
        abline(h = 0.5)
        index_max_cor_2 <- which.max(df[name2,])
        points(index_max_cor_2, df[name2, index_max_cor_2], pch=21, cex = 2, bg='darkred')   
        
        ### Plotting maximum correlation point
        if(index_max_cor_2 >= 100 | index_max_cor_2 <= 750){
          text(index_max_cor_2, df[name2,index_max_cor_2] - 0.075 , labels = names(which.max(df[name2,])), cex = 2, col = 'darkred')                       # Label point to sample that is best correlated to name1
        }else if(index_max_cor_2 < 100){
          text(index_max_cor_2 + 20, df[name2, index_max_cor_2], labels = names(which.max(df[name2,])), cex = 2, col = 'darkred')  
        }else{                                           
          text(index_max_cor_2 - 20, df[name2, index_max_cor_2], labels = names(which.max(df[name2,])), cex = 2, col = 'darkred')  
        }
        
        
        index_name2 <- which(rownames(df) == name2)
        points(index_name2, df[name2, gsub('_rep[12]', '', name2)], pch = 21, cex = 2, bg = 'darkblue')        
        if(index_name2 >= 100 & index_name2 <= 750){          
          if(!is.null(name2)){
            text(index_name2, df[name2, gsub('_rep[12]', '', name2)] -.075, labels = gsub('_rep[12]','',name2), cex = 2, col = 'darkblue', font = 2)
          }else{
            text(index_name2, df[name2, gsub('_rep[12]', '', name2)] -.05, labels = gsub('_rep[12]','',name2), cex = 2, col = 'darkblue', font = 2)
          }
        }else if(index_name2 < 100){                          
          text(index_name2 + 20, df[name2, gsub('_rep[12]', '', name2)] , labels = gsub('_rep[12]','',name2), cex = 2, col = 'darkblue', font = 2)
          
        }else{                            
          text(index_name2 - 20, df[name2, gsub('_rep[12]', '', name2)] , labels = gsub('_rep[12]','',name2), cex = 2, col = 'darkblue', font = 2)
        }
        
      }
    }
    
    
    par(mfrow = c(1,1))
  }
