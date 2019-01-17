### Options and Libraries
options(stringsAsFactors = FALSE)
library(pcaMethods)










### Load and get data
load("~/Desktop/prelim_proteomics_correctedIDs_v2.RData")
raw     <- dataset.esc.proteins$raw
norm    <- dataset.esc.proteins$norm
samples <- dataset.esc.proteins$samples
stopifnot(rownames(raw)  == samples$mouse.id)
stopifnot(rownames(norm) == samples$mouse.id)




### Get PCA of raw data
pca.raw   <- pca(raw, method = 'svdImpute', nPcs = 5)
pca.raw.scores <- scores(pca.raw) 


pca.log.raw    <- pca(log(raw), method = 'bpca', nPcs = 5)
pca.log.scores <- scores(pca.log.raw) 


pca.norm        <- pca(norm, method = 'bpca', nPcs = 5)
pca.norm.scores <- scores(pca.norm) 








### Set up colors by sex
colors <- as.numeric(factor(samples$sex))
colors[colors == 1] <- 'pink'
colors[colors == 2] <- 'blue'





### Identify Males that look like they should be Females. I saw plots beforehand
raw.index  <- which(pca.raw.scores[,1] > 0 & colors == 'blue')
log.index  <- which(pca.log.scores[,1] > 0 & colors == 'blue')
norm.index <- which(pca.norm.scores[,1] > 0 & colors == 'blue')


raw.pc1.index <- pca.raw.scores[raw.index,'PC1']
raw.pc2.index <- pca.raw.scores[raw.index,'PC2']
log.pc1.index <- pca.log.scores[log.index,'PC1']
log.pc2.index <- pca.log.scores[log.index,'PC2']
norm.pc1.index <- pca.norm.scores[norm.index,'PC1']
norm.pc2.index <- pca.norm.scores[norm.index,'PC2']


raw.label  <- rownames(pca.raw.scores)[raw.index]
log.label  <- rownames(pca.log.scores)[log.index]
norm.label <- rownames(pca.norm.scores)[norm.index]







plot(pca.raw.scores[,1],  pca.raw.scores[,2],  pch = 16, col = colors, xlab = 'PC1', ylab = 'PC2', main = 'PCA of Raw by Sex')
text(x = raw.pc1.index + 2, y = raw.pc2.index, labels = raw.label)
plot(pca.log.scores[,1],  pca.log.scores[,2],  pch = 16, col = colors, xlab = 'PC1', ylab = 'PC2', main = 'PCA of Log Raw by Sex')
text(x = log.pc1.index + 2, y = log.pc2.index, labels = log.label)
plot(pca.norm.scores[,1], pca.norm.scores[,2], pch = 16, col = colors, xlab = 'PC1', ylab = 'PC2', main = 'PCA of Norm by Sex')
text(x = norm.pc1.index + 1.5, y = norm.pc2.index, labels = norm.label)








### Looking at the position of the 4 weird samples
raw.index  <- which(rownames(pca.raw.scores) %in% c('PB360.45_B_repB','PB360.76_B_repB','PB366.18_B_repB','PB360.93_B_repB'))
log.index  <- which(rownames(pca.log.scores) %in% c('PB360.45_B_repB','PB360.76_B_repB','PB366.18_B_repB','PB360.93_B_repB'))
norm.index <- which(rownames(pca.norm.scores) %in% c('PB360.45_B_repB','PB360.76_B_repB','PB366.18_B_repB','PB360.93_B_repB'))


raw.pc1.index <- pca.raw.scores[raw.index,'PC1']
raw.pc2.index <- pca.raw.scores[raw.index,'PC2']
log.pc1.index <- pca.log.scores[log.index,'PC1']
log.pc2.index <- pca.log.scores[log.index,'PC2']
norm.pc1.index <- pca.norm.scores[norm.index,'PC1']
norm.pc2.index <- pca.norm.scores[norm.index,'PC2']


raw.label  <- rownames(pca.raw.scores)[raw.index]
log.label  <- rownames(pca.log.scores)[log.index]
norm.label <- rownames(pca.norm.scores)[norm.index]







plot(pca.raw.scores[,1],  pca.raw.scores[,2],  pch = 16, col = colors, xlab = 'PC1', ylab = 'PC2', main = 'PCA of Raw by Sex')
text(x = raw.pc1.index + 2, y = raw.pc2.index, labels = raw.label)
plot(pca.log.scores[,1],  pca.log.scores[,2],  pch = 16, col = colors, xlab = 'PC1', ylab = 'PC2', main = 'PCA of Log Raw by Sex')
text(x = log.pc1.index + 2, y = log.pc2.index, labels = log.label)
plot(pca.norm.scores[,1], pca.norm.scores[,2], pch = 16, col = colors, xlab = 'PC1', ylab = 'PC2', main = 'PCA of Norm by Sex')
text(x = norm.pc1.index + 1.5, y = norm.pc2.index, labels = norm.label)
