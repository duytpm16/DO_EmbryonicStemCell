### Options and Libraries
options(stringsAsFactors = FALSE)
library(clusterProfiler)
library(ggplot2)







load("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_protein_hotspot_clusterprofiler.RData")






cp_data <- 'chr_1_187_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 1 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))









cp_data <- 'chr_4_41_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 41 Hotspot GO - Biological Processes') + theme(plot.title = element_text(hjust = 0.5))


cp_data <- 'chr_7_104_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 7 104 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))






cp_data <- 'chr_9_81_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$x[which(data$name == 'translational initiation')] <- -1
data$x[which(data$name == 'Eif4g3')] <- -0.5119739
data$y[which(data$name == 'Eif4g3')] <- -3.289395e-02
data$x[which(data$name %in% c('Iars','Rps24','Rps9','Gapdh'))] <- 0.3939123
data$y[which(data$name %in% c('Iars','Rps24','Rps9','Gapdh'))] <- -1.317800e-01
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 9 Hotspot GO - Biological Processes') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_9_81_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 9 81 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))


cp_data <- 'chr_9_81_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 9 81 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5))



cp_data <- 'chr_9_93_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 9 93 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5))









cp_data <- 'chr_13_63_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 13 63 Hotspot GO - Biological Processes') + theme(plot.title = element_text(hjust = 0.5))



cp_data <- 'chr_13_63_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 13 63 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5))



cp_data <- 'chr_13_63_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 13 63 Hotspot KEGG Pathway') + theme(plot.title = element_text(hjust = 0.5))













cp_data <- 'chr_15_3_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 15 3 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5))



cp_data <- 'chr_15_3_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 15 3 Hotspot KEGG Pathway') + theme(plot.title = element_text(hjust = 0.5))










cp_data <- 'chr_17_33_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 17 33 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_17_33_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 17 33 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5))



cp_data <- 'chr_17_33_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('dodgerblue', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 17 33 Hotspot KEGG Pathway') + theme(plot.title = element_text(hjust = 0.5))
