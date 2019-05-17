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
load("~/Desktop/Munger Embryonic Stem Cells/Overlap/Version 2/munger_esc_174_mrna_protein_qtl_viewer_v2.RData")
load("~/Desktop/Munger Embryonic Stem Cells/Viewer/Version 2/munger_esc_paired_mrna_proteins_full_qtl_viewer.RData")
transband   <- readRDS("~/Desktop/Munger Embryonic Stem Cells/Viewer/Version 2/munger_esc_proteomics_transband_v2.rds")
target_ds   <- 'dataset.esc.proteins.174'
mediator_ds <- 'dataset.esc.mrna.174'
lodpeaks_ds <- 'dataset.esc.proteins'
target_id   <- 'protein.id'
mediator_id <- 'gene.id'




lod.peaks <- get(lodpeaks_ds)$lod.peaks$additive
lod.peaks <- lod.peaks[!lod.peaks$cis,]
expr      <- get(target_ds)$data$rz
covar     <- get(target_ds)$covar.matrix



annot <- switch(get(mediator_ds)$datatype, "mrna" = 'annot.mrna', "protein" = 'annot.protein') 
med_expr  <- get(mediator_ds)$data$rz
med_annot <- get(mediator_ds)[[annot]] %>% dplyr::rename(pos = start)
med_expr  <- med_expr[rownames(expr), med_annot[,mediator_id, drop = TRUE]]
med_type  <- get(mediator_ds)$datatype














### Plot PCA Mediation
interval <- 20
color    <- 'brown1'
for(i in 9:9){
  
  
  
  # Get QTLs in hotspot interval
  sub_lodpeaks <- lod.peaks %>% subset(qtl.chr == transband$chr[i] & qtl.pos >= transband$start[i] & qtl.pos <= transband$end[i])
  sub_lodpeaks <- sub_lodpeaks[,c(target_id,'gene.symbol'), drop = FALSE]
  sub_lodpeaks$gene.symbol[duplicated(sub_lodpeaks$gene.symbol)] <- paste0(sub_lodpeaks$gene.symbol[duplicated(sub_lodpeaks$gene.symbol)], '-1')
  
  
  
  
  print(paste0(length(intersect(colnames(expr), sub_lodpeaks$protein.id)),',', nrow(sub_lodpeaks)))
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
                          qtl.geno   = genoprobs[[pca_max$chr]][rownames(expr_pca),,rownames(pca_max)],
                          covar      = covar,
                          verbose    = FALSE)    
  med$z <- scale(med$LOD)    
  med   <- med %>% filter(chr %in% transband$chr[i])
  
  
 
  
  
  
  
  
  
  ### Color and label points
  par(xpd = FALSE) 
  plot(med$middle, med$z, col = 'black', pch = 21, bg = color,
       ylim = c(min(med$z) - 2, max(med$z)),
       xlab = paste('Chromosome', transband$chr[i], '(Mbp)'),
       ylab = 'Mediation LOD (Scaled)',
       main = paste0('Chr', transband$chr[i], ' Hotspot PC1 Mediated by Proteins'))
  abline(h = -4)
  abline(v = (transband$start[i] + transband$end[i]) / 2, lwd = 2, lty = 2)
  sig_genes <- which(med$z < -4 & med$chr == transband$chr[i] & med$middle >= transband$start[i] - interval & med$middle <= transband$end[i] + interval)
  
  par(xpd = TRUE) 
  text((transband$start[i] + transband$end[i]) / 2, max(med$z) + .7, 'Hotspot')
  if(length(sig_genes) != 0){
    text(med$middle[sig_genes] + c(-3,-3,3.5), 
         med$z[sig_genes] , labels = med$symbol[sig_genes])    
  }
  
  
}

