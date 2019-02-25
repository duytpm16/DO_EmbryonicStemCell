### Options and Libraries
options(stringsAsFactors = FALSE)
library(clusterProfiler)
library(ensimplR)
library(plyr)
library(dplyr)
setwd('~/Desktop')




### Load and get data
load("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Version 3.1/prelim_proteomics_correctedIDs_v3.1.RData")
dataset <- 'dataset.esc.proteins'
expr.annots <- get(dataset)$annots
lod.peaks <- get(dataset)$lod.peaks$additive







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
rm(min,max,x_axis_scale,chr,sub,count,pos)
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




### Get entrez IDs for each protein in a hotspot
hs_entrez <- list()
for(i in 1:length(hotspots)){
  
    hs_chr <- hotspots[[i]][1]
    hs_pos <- as.numeric(hotspots[[i]][2])
    
    
    
    sub_peaks <- lod.peaks %>% subset(qtl.chr == hs_chr & abs(hs_pos - qtl.pos) <= 2)
    sub_peaks$gene_id <- expr.annots$gene_id[match(sub_peaks$annot.id, expr.annots$protein_id)]
    
    
    entrez <- batchGenes(as.list(unique(sub_peaks$gene_id)), version = 90)
  
    
    hs_entrez[[paste0('chr',hs_chr,'_',hs_pos)]] <- entrez$entrez_id    
  
}




### Cluster Profiler
bg   <- batchGenes(as.list(unique(expr.annots$gene_id)), version = 90)
bp   <- compareCluster(hs_entrez, fun = 'enrichGO', OrgDb = 'org.Mm.eg.db',ont = 'BP', universe = bg$entrez_id, pAdjustMethod = 'fdr', minGSSize = 1, readable = TRUE)
cc   <- compareCluster(hs_entrez, fun = 'enrichGO', OrgDb = 'org.Mm.eg.db',ont = 'CC', universe = bg$entrez_id, pAdjustMethod = 'fdr', minGSSize = 1,readable = TRUE)
kegg <- compareCluster(hs_entrez, fun = 'enrichKEGG', organism = 'mmu',keyType = 'kegg', universe = bg$entrez_id, minGSSize = 1, pAdjustMethod = 'fdr')






### Plot enrichment results
dotplot(bp, title = 'GO - Biological Processes') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))

dotplot(cc, title = 'GO - Cellular Component') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))

dotplot(kegg, title = 'KEGG Pathway') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))
