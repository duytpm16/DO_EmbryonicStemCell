## Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(ensimplR)
library(qtl2)
library(dplyr)
library(sva)








### From Selcan Aydin in Steve Munger lab
load("~/Downloads/DO_mESC_paired_eQTL_forMapping.RData")
mrna_samples <- covarTidy %>% filter(!grepl('repB', sampleid))
mrna_raw     <- esc.raw.expr[,mrna_samples$sampleid]













### Annotate gene IDs
annots    <- batchGenes(ids = as.list(unique(rownames(mrna_raw))), version = 91, species = 'Mm')
na_annots <- annots[is.na(annots$symbol),]
na_annots <- batchGenes(ids = as.list(unique(rownames(na_annots))), version = 82, species = 'Mm')


annots <- annots[!is.na(annots$symbol),]
annots <- rbind(annots, na_annots)
annots <- annots %>%
            mutate(start  = as.numeric(start),
                   end    = as.numeric(end),
                   middle = (start + end) / 2,
                   start  = start / 1e6,
                   end    = end / 1e6,
                   middle = middle / 1e6) %>%
            dplyr::rename(chr = chromosome, description = name, gene.id = gene_id, entrez.id = entrez_id) %>%
            select(gene.id, symbol, chr, start, middle, end, strand, synonyms, entrez.id) %>% 
            arrange(gene.id)













### Get raw data
mrna_raw <- mrna_raw[annots$gene.id,]









### Normalized and Combat
mrna_norm <- mrna_raw
q75       <- apply(mrna_norm, 2, quantile, 0.75)
ratio     <- mean(q75)/q75
mrna_norm <- sweep(mrna_norm, 2, ratio, FUN = "*")
mrna_norm <- log(mrna_norm + 1)

mod <- model.matrix(~sex, data = mrna_samples)
exprComBat <- ComBat(dat = mrna_norm, batch = mrna_samples$libraryprep, mod = mod, par.prior = TRUE, prior.plots = FALSE)
mrna_norm  <- exprComBat









### RankZ
# rankZ transformation of mrna data
rankZ <- function (x) {
  x <- rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
  qnorm(x)
}

mrna_rankz <- apply(mrna_norm, 1, rankZ)













### Samples and Covariate 
mrna_samples <- mrna_samples %>% dplyr::rename(mouse.id = sampleid)

covar.info <- data.frame(sample.column = c('sex'),
                         covar.column  = c('sexM'),
                         display.name  = c('Sex'),
                         interactive   = c(TRUE),
                         primary       = c(TRUE),
                         lod.peaks     = c('sex_int'))
covar.matrix <- model.matrix(~sex, data = mrna_samples)[,-1, drop = FALSE]
rownames(covar.matrix) <- mrna_samples$mouse.id






### QTL2 formatted data
genoprobs <- esc.probs
K         <- calc_kinship(probs = genoprobs, type = 'loco', cores = 10)
markers   <- map_dat2 %>% 
                mutate(pos = as.numeric(pos),
                       pos_bp = as.numeric(pos_bp)) %>%
                dplyr::rename(cM = pos, bp = pos_bp, marker.id = marker) %>%
                mutate(pos = bp / 1e6) %>% 
                select(marker.id, chr, pos, bp, cM) %>%
                as_tibble()
map <- map_df_to_list(map = as.data.frame(markers), 
                      chr_column    = 'chr', 
                      pos_column    = 'pos', 
                      marker_column = 'marker.id')










### Get nearest marker to gene
nearest.marker.id <- character(length = nrow(annots))
for(i in 1:nrow(annots)){
  
  if(annots$chr[i] %in% c(1:19,'X')){
    sub <- markers %>% filter(chr == annots$chr[i])
    nearest.marker.id[i] <- sub$marker.id[which.min(abs(annots$middle[i] - sub$pos))]
  }
}










### QTL2 viewer format
dataset.esc.mrna <- list(annot.mrna    = as_tibble(annots),
                         annot.samples = as_tibble(mrna_samples),
                         covar.matrix  = covar.matrix,
                         covar.info    = as_tibble(covar.info),
                         data          = list(norm = mrna_norm,
                                              raw  = mrna_raw,
                                              rz   = mrna_rankz),
                         datatype      = 'mrna',
                         display.name  = 'Embryonic Stem Cell Transcripts',
                         lod.peaks     = list())







### Save
rm(list = ls()[!ls() %in% c(grep('dataset[.]', ls(), value = TRUE),'genoprobs','K','map','markers')])
save.image(file = '~/Desktop/Munger Embryonic Stem Cells/RNA/Version 2 - Without Duplicates/munger_esc_mrna_qtl_viewer_v2.RData')


