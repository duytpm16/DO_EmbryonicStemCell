### Options and Libraries
options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(intermediate)
library(pcaMethods)
library(qtl2)








### Load and extract data
#load("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Version 1/munger_esc_proteomics_qtl_viewer_v1.RData")
transband   <- readRDS("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Version 1/munger_esc_proteomics_transband.rds")
target_ds   <- 'dataset.esc.proteins.174'
mediator_ds <- 'dataset.esc.mrna.174'
lodpeaks_ds <- 'dataset.esc.proteins'
target_id   <- 'protein.id'
mediator_id <- 'protein.id'




lod.peaks <- get(lodpeaks_ds)$lod.peaks$additive
lod.peaks <- lod.peaks[!lod.peaks$cis,]
expr      <- get(target_ds)$data$rz
covar     <- get(target_ds)$covar.matrix



annot <- switch(get(mediator_ds)$datatype, "mrna" = 'annot.mrna', "protein" = 'annot.protein') 
med_expr  <- get(mediator_ds)$data$rz
med_annot <- get(mediator_ds)[[annot]] %>% dplyr::rename(pos = start)
med_expr  <- med_expr[rownames(expr), med_annot[,mediator_id, drop = TRUE]]
med_type  <- get(mediator_ds)$datatype












### Get transband with counts above 60
transband <- transband[transband$distant_esc_prot > 60,]











### Plot PCA Mediation
interval <- 5
color    <- 'brown1'
for(i in 5:5){
  
  
  
  # Get QTLs in hotspot interval
  sub_lodpeaks <- lod.peaks %>% subset(qtl.chr == transband$chr[i] & qtl.pos >= transband$start[i] & qtl.pos <= transband$end[i])
  sub_lodpeaks <- sub_lodpeaks[,c(target_id,'gene.symbol', LETTERS[1:8]), drop = FALSE]
  sub_lodpeaks$gene.symbol[duplicated(sub_lodpeaks$gene.symbol)] <- paste0(sub_lodpeaks$gene.symbol[duplicated(sub_lodpeaks$gene.symbol)], '-1')
  
  
  
  
  print(paste0(length(intersect(colnames(expr), sub_lodpeaks$protein.id)), nrow(sub_lodpeaks)))
  # Get Expressions of proteins in hotspot interval and compute the PCA
  sub_expr <- expr[, intersect(colnames(expr), sub_lodpeaks$protein.id)]
  expr_pca <- pca(sub_expr, method = 'svdImpute')
  expr_pca <- scores(expr_pca)[,1,drop = FALSE]
  
  
  
  
  
  ### LOD scan of PCA and mediate PCA against all mediators
  pca_scan1 <- scan1(genoprobs = genoprobs, pheno = expr_pca, kinship = K, addcovar = covar)
  pca_max   <- max_scan1(scan1_output = pca_scan1, map = map, chr = transband$chr[i])
  
  
  med   <- mediation.scan(target     = expr_pca,
                          mediator   = med_expr,
                          annotation = med_annot,
                          qtl.geno   = genoprobs[[pca_max$chr]][,,rownames(pca_max)],
                          covar      = covar,
                          verbose    = FALSE)    
  med$z <- scale(med$LOD)    
  med   <- med %>% filter(chr %in% c(1:19,'X'))
  xpos  <- xpos_scan1(map = map , thechr = med$chr, thepos = med$middle)
  med   <- cbind(med,xpos)
  
  
  
  
  
  
  
  
  
  ### Formatting plot background  
  rect_border <- numeric(length = 21)
  middle_pos  <- numeric(length = 20)
  chromosomes <- c(1:19,'X')
  for( j in 1:length(chromosomes)){
    
    if(j == 20){
      left  <- med %>% subset(chr == chromosomes[j])
      right <- med %>% subset(chr == chromosomes[20])
      rect_border[21] <- max(right$xpos)
      middle_pos[j]   <- (min(left$xpos) + max(left$xpos)) / 2
    }
    else{
      left  <- med %>% subset(chr == chromosomes[j])
      right <- med %>% subset(chr == chromosomes[j+1])
      rect_border[j+1] <- (max(left$xpos) + min(right$xpos)) / 2
      middle_pos[j] <- (min(left$xpos) + max(left$xpos)) / 2
    }
  }
  
  
  
  
  
  
  
  
  
  
  ### Plot result
  plot(x = med$xpos, y = med$z, 
       ylim = c(min(med$z -2), max(med$z)), 
       xaxs = 'i', yaxs='i', 
       xaxt = 'n', las = 1,
       xlab = 'Chromosomes', ylab = 'Mediation LOD (Scaled)', 
       main = paste0('Chr', transband$chr[i], ' HotSpot PC1 Global Mediation by ', med_type))
  axis(1, at = middle_pos, labels = c(1:19,'X'))
  for(j in 2:21){
    if(j %% 2 != 0){
      rect(xleft = rect_border[j-1], ybottom = min(med$z) - 2, xright = rect_border[j], ytop = max(med$z), border = 'white', col = 'white')
    }else{
      rect(xleft = rect_border[j-1], ybottom = min(med$z) - 2, xright = rect_border[j], ytop = max(med$z), border = 'white', col = 'gray90')
    }
  }
  
  
  
  
  
  
  
  ### Color and label points
  cl <- character(length = nrow(med))
  cl[med$chr == transband$chr[i]] <- color
  cl[med$chr != transband$chr[i]] <- 'grey50'
  points(xpos, med$z, col = 'black', pch = 21, bg = cl)
  abline(h = -4)
  sig_genes <- which(med$z < -4 & med$chr == transband$chr[i] & med$middle >= transband$start[i] - interval & med$middle <= transband$end[i] + interval)
  
  if(length(sig_genes) != 0){
    text(med$xpos[sig_genes] + c(70,110,110) , med$z[sig_genes], labels = med$symbol[sig_genes])    
  }
  
}
