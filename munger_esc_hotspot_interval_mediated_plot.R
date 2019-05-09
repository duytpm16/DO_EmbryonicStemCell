options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(intermediate)





load("~/Desktop/Munger Embryonic Stem Cells/Viewer/Version 1/munger_esc_paired_mrna_proteins_full_qtl_viewer_v1.RData")
transband <- readRDS("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Version 1/munger_esc_proteomics_transband.rds")




target_ds   <- 'dataset.esc.proteins'
mediator_ds <- 'dataset.esc.proteins'
target_id   <- 'protein.id'
mediator_id <- 'protein.id'




lod.peaks <- get(target_ds)$lod.peaks$additive
lod.peaks <- lod.peaks[!lod.peaks$cis,]
expr      <- get(target_ds)$data$rz
covar     <- get(target_ds)$covar.matrix



annot <- switch(get(mediator_ds)$datatype, "mrna" = 'annot.mrna', "protein" = 'annot.protein') 
med_expr  <- get(mediator_ds)$data$rz
med_annot <- get(mediator_ds)[[annot]] %>% dplyr::rename(pos = start)
med_expr  <- med_expr[rownames(expr), med_annot[,mediator_id, drop = TRUE]]
med_type  <- get(mediator_ds)$datatype

















interval <- 5
for(i in 13:13){
  
  # Get transband info
  hs_chr <- as.character(transband$chr[i])
  hs_start <- transband$start[i]
  hs_end   <- transband$end[i]
  
  
  
  # Subset lodpeaks to proteins in hotspot interval
  sub_lodpeaks <- lod.peaks %>% subset(qtl.chr == hs_chr & qtl.pos >= hs_start & qtl.pos <= hs_end)
  sub_expr     <- expr[,sub_lodpeaks[,target_id, drop = TRUE]]
  stopifnot(sub_lodpeaks[,target_id, drop = TRUE] == colnames(sub_expr))
  
  
  
  # Get mediators in the hotspot interval
  med_in_interval <- med_annot %>% subset(chr == hs_chr & middle >= hs_start - interval & middle <= hs_end + interval) %>% dplyr::arrange(pos)
  med_in_interval_expr <- med_expr[,med_in_interval[,mediator_id, drop = TRUE]]
  stopifnot(colnames(med_in_interval_expr) == med_in_interval[,mediator_id, drop = TRUE])
  
  
  
  
  # Arranging proteins in hotspot interval based on correlation
  corr  <- cor(sub_expr, use = 'pairwise.complete.obs')
  cl    <- hclust(as.dist(1.0 - corr), method = "average") 
  order <- sub_lodpeaks[cl$order,]
  
  
  
  
  # Create mediation result matrix
  med_result  <- as.data.frame(matrix(0, nrow = nrow(order), ncol = ncol(med_in_interval_expr),dimnames = list(order[,target_id, drop = TRUE], colnames(med_in_interval_expr))))
  stopifnot(rownames(sub_expr) == rownames(med_in_interval_expr))
  stopifnot(colnames(med_in_interval_expr) == colnames(med_result))
  stopifnot(rownames(sub_expr) == rownames(covar))
  
  
  
  
  
  
  for(j in 1:nrow(order)){
    med <- mediation.scan(target     = sub_expr[, order[,target_id, drop = TRUE][j], drop = FALSE],
                          mediator   = med_in_interval_expr,
                          qtl.geno   = genoprobs[[order$qtl.chr[j]]][rownames(sub_expr),,order$marker.id[j]],
                          annotation = med_in_interval,
                          covar      = covar,
                          method     = 'double-lod-diff',
                          verbose    = FALSE)
    
    med$dp <- (order$lod[j] - med$LOD) / order$lod[j]
    med$dp[med$dp < 0] <- 0
    med_result[order[,target_id, drop = TRUE][j], med[,mediator_id, drop = TRUE]] <- med$dp
    
  }
  stopifnot(rownames(med_result) == order[,target_id, drop = TRUE])
  stopifnot(colnames(med_result) == med_in_interval[,mediator_id, drop = TRUE])
  
  
  
  
  
  ### Converting Ensembl IDs to symbol. Add '-1' to duplicate symbols
  med_in_interval$symbol[duplicated(med_in_interval$symbol)] <- paste0(med_in_interval$symbol[duplicated(med_in_interval$symbol)],'-1')
  order$gene.symbol[duplicated(order$gene.symbol)]           <- paste0(order$gene.symbol[duplicated(order$gene.symbol)],'-1')
  rownames(med_result) <- order$gene.symbol
  colnames(med_result) <- med_in_interval$symbol
  
  
  
  
  
  
  ### Melt for ggplot2
  med_melt <- cbind(target = rownames(med_result), med_result)
  med_melt <- melt(med_melt)
  med_melt$target <- factor(med_melt$target, levels = rownames(med_result))
  med_melt$variable <- factor(med_melt$variable, levels = colnames(med_result))
  med_melt <- med_melt %>% plyr::rename(c("value" = "Proportion Drop"))
  
  
  
  
  
  
  
  print(ggplot(med_melt, aes(x = variable, y = target)) + geom_tile(aes(fill = med_melt$`Proportion Drop`), color = 'black')  + 
          xlab(paste0('Chr',hs_chr, '  ~ ', round(hs_start - interval), '-', round(hs_end + interval), ' Mbp')) + 
          ylab(paste0('Proteins in Chr',hs_chr,' HotSpot')) +
          theme(plot.title = element_text(hjust = .5, size = 18),
                axis.title = element_text(size= 15),
                axis.text.y = element_text(face = 'bold', color =  'black', size = 6),
                axis.text.x = element_text(angle = 90, vjust = .5, face = 'bold', color = 'black'),
                panel.background = element_blank(),
                legend.title = element_blank()) +
          scale_fill_gradientn(colors = c('dodgerblue','orange','yellow')) + ggtitle(paste0('Chr',hs_chr,' HotSpot Interval Mediation')))
  
}
