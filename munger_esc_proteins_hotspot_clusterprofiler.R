options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(plyr)
library(data.table)
library(ensimplR)


load("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_proteomics_qtl_viewer_v2.RData")
transband <- readRDS("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_proteomics_transband.rds")
lod.peaks <- dataset.esc.proteins$lod.peaks$additive
lod.peaks <- lod.peaks[!lod.peaks$cis,]
annots    <- dataset.esc.proteins$annot.protein







transband <- transband[transband$distant_esc_prot > 60,]










bg <- batchGenes(ids = unique(as.list(annots$gene.id)), version = 91)


for(i in 1:nrow(transband)){
  
  sub_lodpeaks <- lod.peaks %>% subset(qtl.chr == transband$chr[i] & qtl.pos >= transband$start[i] & qtl.pos <= transband$end[i])
  sub_lodpeaks$gene.id <- annots$gene.id[match(sub_lodpeaks$protein.id, annots$protein.id)]
  sub_lodpeaks$entrez.id <- bg$entrez_id[match(sub_lodpeaks$gene.id, bg$gene_id)]
  
  assign(paste0('chr_', transband$chr[i],'_',round(transband$start[i]),'_BP'),
         enrichGO(sub_lodpeaks$gene.id, universe = annots$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL', ont = 'BP', pAdjustMethod = 'fdr', readable = TRUE))
  
  assign(paste0('chr_', transband$chr[i],'_',round(transband$start[i]),'_CC'),
         enrichGO(sub_lodpeaks$gene.id, universe = annots$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL', ont = 'CC', pAdjustMethod = 'fdr', readable = TRUE))
  
  assign(paste0('chr_', transband$chr[i],'_',round(transband$start[i]),'_MF'),
         enrichGO(sub_lodpeaks$gene.id, universe = annots$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL', ont = 'MF', pAdjustMethod = 'fdr', readable = TRUE))
  
  
  kegg <- enrichKEGG(sub_lodpeaks$entrez.id, universe = bg$entrez_id, organism = 'mmu', keyType = 'kegg', pAdjustMethod = 'fdr')
  for(j in 1:nrow(kegg@result)){
    kegg@result$geneID[j] <- paste0(bg$symbol[match(strsplit(kegg@result$geneID[j], split = '/')[[1]], bg$entrez_id)], collapse = '/')
  }
  assign(paste0('chr_', transband$chr[i],'_',round(transband$start[i]),'_KEGG'), kegg)
}







rm(list = ls()[!ls() %in% grep('chr_',ls(),value = TRUE)])
save.image("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_protein_hotspot_clusterprofiler.RData")





save.image("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/Version 2/attie_islet_protein_hotspot_clusterprofiler.RData")