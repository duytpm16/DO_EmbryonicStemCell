### Options and Libraries
options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(plyr)
library(data.table)
library(plyr)






### Load and get required data
load("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_proteomics_qtl_viewer_v2.RData")
transband <- readRDS("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_proteomics_transband.rds")
lod.peaks <- dataset.esc.proteins$lod.peaks$additive
lod.peaks <- lod.peaks[!lod.peaks$cis,]
annots    <- dataset.esc.proteins$annot.protein
expr      <- dataset.esc.proteins$data$rz






### Find hotspots
transband <- transband[transband$distant_esc_prot > 60,]











### Plot begins
for(i in 1:nrow(transband)){
  
    # Get QTLs in hotspot interval
    sub_lodpeaks <- lod.peaks %>% subset(qtl.chr == transband$chr[i] & qtl.pos >= transband$start[i] & qtl.pos <= transband$end[i])
    sub_lodpeaks <- sub_lodpeaks[,c('protein.id','gene.symbol', LETTERS[1:8]), drop = FALSE]
    sub_lodpeaks$gene.symbol[duplicated(sub_lodpeaks$gene.symbol)] <- paste0(sub_lodpeaks$gene.symbol[duplicated(sub_lodpeaks$gene.symbol)], '-1')
    
    
    # Get Expressions of proteins in hotspot interval and compute the PCA
    sub_expr <- expr[, sub_lodpeaks$protein.id]
    expr_pca <- pca(sub_expr, method = 'svdImpute')
    expr_pca <- scores(expr_pca)[,1,drop = FALSE]
    
    
    # Scan1 of PCA
    pca_scan1 <- scan1(genoprobs = genoprobs, pheno = expr_pca, kinship = K, addcovar = dataset.esc.proteins$covar.matrix)
    max_pca   <- max_scan1(pca_scan1, map = map, chr = transband$chr[i])
    stopifnot((max_pca$pos >= transband$start[i] & max_pca$pos <= transband$end[i]))
    
    # BLUP scan of the PCA
    gp <- genoprobs[,max_pca$chr]
    gp[[1]]  <- gp[[1]][,,rownames(max_pca), drop = FALSE]
    pca_blup <- scan1blup(genoprobs = gp[,transband$chr[i]], pheno = expr_pca, kinship = K[[transband$chr[i]]], addcovar = dataset.esc.proteins$covar.matrix)
    
    
    # Create a df of blup results for ggplot
    blup_df  <- data.frame(founders = c('   AJ', '   B6', '  129','  NOD','  NZO',' CAST',' PWK',' WSB'), blup = c(pca_blup[,LETTERS[1:8]]))
    blup_df$founders <- factor(blup_df$founders, levels = rev(blup_df$founders))
    
    
    
    
    
    
  
    
    
    # Find PCA allele effect range
    r <- c(-round_any(max(abs(blup_df$blup)), 1, f = ceiling), round_any(max(abs(blup_df$blup)), 1, f = ceiling))
    
    # Plot allele effect Greg Keele style
    l <- ggplot(data = blup_df, aes(x = founders, y = blup, group = 1)) + 
            geom_point(size = 4, col = 'grey40') + 
            geom_line(size = 1.2, col = 'cyan') +
            ylab('Allele Effect') +
            ggtitle('PC1') + 
            scale_y_continuous(limits = c(r[1],r[2]), breaks = seq(r[1],r[2], round(r[2]/3))) +
            scale_x_discrete(position = 'top') +
            theme(panel.background = element_rect(fill= 'white'), 
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_line(color = 'white', linetype = 2),
                  axis.line.x = element_line(size = .5),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_text(size = 15, color = 'black'),
                  axis.text.x = element_text(size = 12, color = 'black'),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(face = 2, color =  'black'),
                  plot.title = element_text(hjust = .5, face = 2)) + 
            geom_hline(yintercept = 0, lty = 2, col = 'grey') +
            coord_flip()
    

    
    
    
    
    
    
    
    
    
    
    # Correlate expression of proteins in hotspot interval and reorder
    sub_expr <- expr[,sub_lodpeaks$protein.id]
    expr_cor <- cor(sub_expr, use = 'pairwise.complete.obs')
    cl <- hclust(as.dist(1 - expr_cor), method = 'average')
    sub_lodpeaks <- sub_lodpeaks[cl$order,]

    
    
    # Melt for ggplot
    m <- melt(sub_lodpeaks)
    m$gene.symbol <- factor(m$gene.symbol, levels = sub_lodpeaks$gene.symbol)
    m$variable    <- factor(m$variable, levels = rev(LETTERS[1:8]))
    
    
    # Plot heatmap of proteins allele effect
    hm <- ggplot(m, aes(x = gene.symbol, y = variable, fill = value)) + geom_tile(color = 'grey60') + ylab('Genes in Hotspot') +
              theme(axis.text.y = element_blank(),
                    axis.text.x = element_text(angle = 90, color = 'black'),
                    axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    plot.title = element_text(color = 'black', hjust = .5, face = 2)) + 
              ggtitle(paste(paste0('Chr',transband$chr[i]), '~', paste0(round(transband$start[i]),'-',round(transband$end[i])), 'Allele Effects')) +
              scale_fill_gradientn(colors = c('cadetblue2','azure1','brown2')) + theme(panel.background = element_blank(),legend.position = 'none')
    
    
    
    
    
    
    
    
    
    
    # Arrange both plots side by side
    cowplot::plot_grid(l, hm, align = 'h', rel_widths = c(.25, 1))

}

