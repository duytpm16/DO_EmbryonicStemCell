### Options and Libraries
options(stringsAsFactors = FALSE)
library(plyr)
library(dplyr)
library(intermediate)
library(data.table)
library(gplots)
library(tidyverse)
setwd('~/Desktop')










### Load and get data
#load("~/Desktop/DO_mESC/Proteomics/Version 3.1/prelim_proteomics_correctedIDs_v3.1.RData")
load("~/Desktop/munger_esc_viewer_174_v2.RData")
targ_id <- 'gene_id'
med_id <- 'gene_id'
targ_data <- 'dataset.esc.mrna'
med_data <- 'dataset.esc.mrna'
lod.peaks   <- get(targ_data)$lod.peaks$additive
targ.expr   <- get(targ_data)$rankz
med.expr    <- get(med_data)$rankz
targ.covar  <- get(targ_data)$covar
med.covar   <- get(med_data)$covar
targ.annots <- get(targ_data)$annots
med.annots  <- get(med_data)$annots















### Counting number of LOD above 6 within a 4MB window across each chromosome
lod.peaks <- lod.peaks[lod.peaks$cis == FALSE,]
lod_df <- list()
slide  <- 1
window <- 4
for(i in unique(markers$chr)){
  
  # Finding floor of minimum marker position and ceiling of maximum marker position
  min <- round_any(min(map[[i]]), 1, f = floor)
  max <- round_any(max(map[[i]]), 4, f = ceiling)
  
  # Creating x-axis scale. min to max with slide (or 'by' in seq function)
  x_axis_scale <- seq(min, max, slide)
  chr <- rep(i, length(x_axis_scale))
  
  # Getting LOD peaks from chromosome i
  sub <- subset(lod.peaks, qtl.chr == i)
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}

lod_df <- rbindlist(lod_df)
lod_df <- lod_df[lod_df$count > 50,]







hotspots <- list(chr1 = c('1', '189'),
                 chr2 = c('2', '173'),
                 chr4 = c('4', '150'),
                 chr5 = c('5', '148'),
                 chr7 = c('7', '110'),
                 chr8 = c('8', '123'),
                 chr9 = c('9', '85'),
                 chr10 = c('10','119'),
                 chr11 = c('11','112'),
                 chr12 = c('12','35'),
                 chr14 = c('14','33'),
                 chr15 = c('15','8'))








hs_med_result <- list(chr1 = data.frame(),
                      chr2 = data.frame(),
                      chr4 = data.frame(),
                      chr5 = data.frame(),
                      chr7 = data.frame(),
                      chr8 = data.frame(),
                      chr9 = data.frame(),
                      chr10 = data.frame(),
                      chr11 = data.frame(),
                      chr12 = data.frame(),
                      chr14 = data.frame(),
                      chr15 = data.frame())







### Do mediation for each hostspot
for(i in 1:length(hotspots)){
  
  
  # Get hotspot info
  hs_chr <- hotspots[[i]][1]  
  hs_pos <- as.numeric(hotspots[[i]][2])
  
  
  
  
  # Extract targets and mediators in the hotspot
  sub_targ_peaks  <- lod.peaks %>% subset(qtl.chr == hs_chr & abs(qtl.pos - hs_pos) <= 4/2)
  sub_targ_annots <- targ.annots[targ.annots[,targ_id] %in% sub_targ_peaks$annot.id,]
  sub_targ_expr   <- targ.expr[,colnames(targ.expr) %in% sub_targ_annots[,targ_id], drop = FALSE]
  
  sub_med_annots <- med.annots %>% subset(chr == hs_chr & abs(start - hs_pos) <= 5) %>% dplyr::rename(pos = start)
  sub_med_expr   <- med.expr[,colnames(med.expr) %in% sub_med_annots[,med_id]]
  stopifnot(colnames(sub_med_expr) == sub_med_annots[,med_id])
  
  
  
  
  
  
  # Create data.frame to store results
  med_results <- as.data.frame(matrix(0, 
                                      nrow = nrow(sub_targ_peaks), 
                                      ncol = nrow(sub_med_annots),
                                      dimnames = list(sub_targ_peaks$annot.id, sub_med_annots[,med_id])))
  stopifnot(rownames(med_results) == sub_targ_peaks$annot.id)
  stopifnot(colnames(med_results) == sub_med_annots[,med_id])
  stopifnot(rownames(sub_med_expr) == rownames(sub_targ_expr) & rownames(sub_med_expr) == rownames(targ.covar))
  
  
  
  
  
  # Get QTL info and run mediation
  qtl.annots <- sub_targ_peaks$annot.id
  qtl.chr    <- sub_targ_peaks$qtl.chr
  qtl.marker <- sub_targ_peaks$marker.id
  qtl.lod    <- sub_targ_peaks$lod
  
  
  
  
  
  for(j in 1:nrow(sub_targ_peaks)){
    
    
    # Mediation: Y ~ Q  + M + covar
    med <- mediation.scan(target     = sub_targ_expr[, qtl.annots[j], drop = FALSE],
                          mediator   = sub_med_expr,
                          annotation = sub_med_annots,
                          qtl.geno   = genoprobs[[qtl.chr[j]]][,,qtl.marker[j]],
                          covar      = targ.covar,
                          method     = "double-lod-diff",
                          verbose    = FALSE)
    med  <- med[colnames(med_results),]
    
    med_results[qtl.annots[j],] <- (qtl.lod[j] - med$LOD) / qtl.lod[j]
    stopifnot(rownames(med) == colnames(med_results))
  }
  
  
  
  
  # Arrange targets by how correlated they are
  expr.cor <- cor(sub_targ_expr, use = 'pairwise.complete.obs')
  expr.cor <- expr.cor[rownames(med_results),]
  cl       <- hclust(as.dist(1.0 - expr.cor), method = "average")
  
  stopifnot(rownames(expr.cor) == rownames(med_results))
  med_results <- med_results[cl$order,]
  
  
  # Arrange mediators by position
  med_results <- med_results[,sub_med_annots[,med_id]]
  genomic.pos <- order(sub_med_annots$pos)
  stopifnot(colnames(med_results) == sub_med_annots[,med_id])   
  
  
  med_results <- med_results[,genomic.pos]
  
  
  prot_name <- sub_targ_annots$symbol[match(rownames(med_results), sub_targ_annots[,targ_id])]
  prot_name[duplicated(prot_name)] <- paste0(prot_name[duplicated(prot_name)],'-1')
  
  med_name  <- sub_med_annots$symbol[match(colnames(med_results), sub_med_annots[,med_id])]
  
  rownames(med_results) <- prot_name
  colnames(med_results) <- med_name
  
  
  hs_med_result[[i]] <- med_results
}







hmcols<-colorRampPalette(c("white","red"))(256)
for(i in 1:length(hs_med_result)){
  hs_chr <- hotspots[[i]][1]  
  hs_pos <- as.numeric(hotspots[[i]][2])
  data <- hs_med_result[[i]]
  data[data < 0] <- 0
  data <- data %>%
               rownames_to_column()
  m_data <- melt(data)
  m_data$rowname <- factor(m_data$rowname, levels = unique(m_data$rowname))
  
  if(med_id == 'gene_id'){
     dot_color <- 'darkred'
  }else{
     dot_color = 'blue'}
  
  print(ggplot(m_data, aes(x = variable, y = rowname, col = value)) + geom_point(size = m_data$value * 10)  + 
          xlab(paste0('Chr',hs_chr,' ~',hs_pos-5,'-',hs_pos+5,' Mb')) + 
          ylab(paste0('Proteins in Chr',hs_chr,' HotSpot')) +
          scale_color_gradient(low = 'green', high = 'red') + ggtitle(paste0('Protein Chr',hs_chr,' HotSpot Mediation')))
}

