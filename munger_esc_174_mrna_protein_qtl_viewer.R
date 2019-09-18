### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(ensimplR)
library(abind)
library(dplyr)
library(qtl2)
library(sva)





### rankZ transformation helper function
rankZ <- function (x) {
  x <- rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
  qnorm(x)
}







### Load proteomics formatted data
load('~/Desktop/Munger Embryonic Stem Cells/Proteomics/Version 2 - No TMT/munger_esc_proteomics_qtl_viewer_v2.RData')
prot_genoprobs <- genoprobs
prot_samples   <- dataset.esc.proteins$annot.samples
prot_raw       <- dataset.esc.proteins$data$raw
prot_ref       <- readRDS("~/Desktop/Ensembl/Ensembl_90/ref_protein_info.rds")


### Load mrna formatted data
load('~/Desktop/Munger Embryonic Stem Cells/RNA/Version 2/munger_esc_mrna_qtl_viewer_v2.RData')
mrna_genoprobs <- genoprobs
mrna_samples   <- dataset.esc.mrna$annot.samples
mrna_raw       <- dataset.esc.mrna$data$raw








### Get overlapping sample information
overlap_prot_samples <- prot_samples %>% filter(mouse.id %in% mrna_samples$mouse.id)
overlap_mrna_samples <- mrna_samples %>% filter(mouse.id %in% prot_samples$mouse.id)
stopifnot(overlap_prot_samples$sex == overlap_mrna_samples$sex)
stopifnot(overlap_prot_samples$mouse.id == overlap_mrna_samples$mouse.id)

overlap_samples <- overlap_prot_samples
overlap_samples <- merge(overlap_samples, mrna_samples[,c('mouse.id','libraryprep','flowcell')], by = 'mouse.id', all.x = TRUE)




### Get overlapping genoprobs
prot_genoprobs <- probs_qtl2_to_doqtl(probs = prot_genoprobs)
mrna_genoprobs <- probs_qtl2_to_doqtl(probs = mrna_genoprobs)
stopifnot(prot_genoprobs[overlap_samples$mouse.id,,] == mrna_genoprobs[overlap_samples$mouse.id,,])


overlap_genoprobs <- prot_genoprobs[overlap_samples$mouse.id,,]
overlap_genoprobs <- probs_doqtl_to_qtl2(probs = overlap_genoprobs, 
                                         map   = as.data.frame(markers), 
                                         chr_column    = 'chr', 
                                         pos_column    = 'pos', 
                                         marker_column = 'marker.id')
genoprobs <- overlap_genoprobs
K         <- calc_kinship(probs = genoprobs, type = 'loco')










### Get protein annotations
prot_annots <- data.frame(protein_id = colnames(dataset.esc.proteins$data$raw))
prot_annots <- prot_annots %>% filter(grepl('^ENSMUSP', protein_id))
prot_annots <- merge(prot_annots, prot_ref[,c('protein_id','gene_id','transcript_id','gene_biotype','description')], by = 'protein_id', all.x = TRUE)
prot_annots$gene_id <- gsub('\\..*', '', prot_annots$gene_id)

gene_annots <- batchGenes(ids = as.list(unique(prot_annots$gene_id)), species = 'Mm', version = 90)
prot_annots <- merge(prot_annots, gene_annots[,c('gene_id','symbol','chromosome','start','end','strand','synonyms','entrez_id')], by = 'gene_id', all.x = TRUE)
prot_annots <- prot_annots %>% 
                  mutate(middle = (as.numeric(start) + as.numeric(end)) / 2, 
                         middle = middle / 1e6, 
                         start = as.numeric(start) / 1e6, 
                         end = as.numeric(end) / 1e6) %>%
                  dplyr::select(protein_id, gene_id, transcript_id, symbol, chromosome, start, middle, end, strand, entrez_id, gene_biotype, synonyms) %>%
                  dplyr::rename(chr = chromosome, protein.id = protein_id, gene.id = gene_id, transcript.id = transcript_id, entrez.id = entrez_id, gene.biotype = gene_biotype) %>%
                  filter(chr %in% c(1:19,'X','Y','MT'))





### Get protein raw, normalized, and rankz data
prot_raw <- dataset.esc.proteins$data$raw


# Gygi lab already normalized so we will just filter out proteins with more than half the sample missing
prot_norm <- prot_raw[,colnames(prot_raw) %in% prot_annots$protein.id]
prot_norm <- prot_norm[rownames(prot_norm) %in% prot_samples$mouse.id,]
keep <- apply(prot_norm, 2, FUN = function(x) mean(!is.na(x)) > 0.5)


# Get filtered protein set, rankZ, and filter annot dataframe to have same proteins
prot_norm   <- prot_norm[,keep]
prot_rankz  <- apply(prot_norm, 2, rankZ)
prot_annots <- prot_annots[prot_annots$protein.id %in% colnames(prot_norm),]




### Covariate Info and covariate matrix
prot_covarinfo <- data.frame(sample.column = c('sex'),
                             covar.column  = c('sexM'),
                             display.name  = c('Sex'),
                             interactive   = c(TRUE),
                             primary       = c(TRUE),
                             lod.peaks     = c('sex_int'))
prot_covarmatrix <- model.matrix(~ sex, data = overlap_samples)[, -1, drop = FALSE]
rownames(prot_covarmatrix) <- overlap_samples$mouse.id










### Get mrna annots and raw mrna data. 
mrna_annots <- dataset.esc.mrna$annot.mrna
mrna_raw    <- dataset.esc.mrna$data$raw %>% as.data.frame()
mrna_raw    <- mrna_raw[,colnames(mrna_raw) %in% overlap_samples$mouse.id]
mrna_raw    <- as.matrix(mrna_raw)



### Normalizing using Dan Skelly and Selcan Aydin method
q75   <- apply(mrna_raw, 2, quantile, 0.75)
ratio <- mean(q75)/q75
mrna_norm <- sweep(mrna_raw, 2, ratio, FUN = "*")
mrna_norm <- log(mrna_norm + 1)


# ComBat
mod <- model.matrix(~sex, data = overlap_samples)
exprComBat <- ComBat(dat = mrna_norm, batch = overlap_samples$libraryprep, mod = mod, par.prior = TRUE, prior.plots = FALSE)


mrna_norm  <- exprComBat
mrna_rankz <- apply(mrna_norm, 1, rankZ)




### Making the mrna covar matrix like proteins since we are just adjusting for sex
mrna_covarmatrix <- prot_covarmatrix
mrna_covarinfo   <- prot_covarinfo









### Ordering the annots and samples dataframe
prot_annots <- prot_annots %>% arrange(protein.id)
mrna_annots <- mrna_annots %>% arrange(gene.id)
overlap_samples <- overlap_samples %>% arrange(mouse.id)










### QTL Viewer format
dataset.esc.proteins.174 <- list(annot.protein  = as_tibble(prot_annots),
                                 annot.samples  = as_tibble(overlap_samples),
                                 covar.info     = as_tibble(prot_covarinfo),
                                 covar.matrix   = as.matrix(prot_covarmatrix),
                                 data           = list(norm = as.matrix(prot_norm)[overlap_samples$mouse.id, prot_annots$protein.id],
                                                       raw  = as.matrix(prot_raw),
                                                       rz   = as.matrix(prot_rankz)[overlap_samples$mouse.id, prot_annots$protein.id]),
                                 datatype       = 'protein',
                                 display.name   = 'DO 174 Embryonic Stem Cell Proteome',
                                 lod.peaks      = list())




dataset.esc.mrna.174 <- list(annot.mrna     = as_tibble(mrna_annots),
                             annot.samples  = as_tibble(overlap_samples),
                             covar.info     = as_tibble(mrna_covarinfo),
                             covar.matrix   = as.matrix(mrna_covarmatrix),
                             data           = list(norm = as.matrix(t(mrna_norm)[overlap_samples$mouse.id, mrna_annots$gene.id]),
                                                   raw  = as.matrix(t(mrna_raw)),
                                                   rz   = as.matrix(mrna_rankz)[overlap_samples$mouse.id, mrna_annots$gene.id]),
                             datatype       = 'mrna',
                             display.name   = 'DO 174 Embryonic Stem Cell Transcriptome',
                             lod.peaks      = list())






### Save
rm(list = ls()[!ls() %in% c(grep('dataset.esc.*.174', ls(), value = TRUE), 'genoprobs', 'K', 'map', 'markers')])
save.image(file = '~/Desktop/Munger Embryonic Stem Cells/Overlap/Version 2 - No TMT - Duplicates/munger_174_esc_qtl_viewer_v2.RData')
