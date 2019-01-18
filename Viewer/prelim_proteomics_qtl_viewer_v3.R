### Options and Libraries 
options(stringsAsFactors = FALSE)
library(ensimplR)                 # devtools::install_github('https://github.com/churchill-lab/ensimplR')
library(dplyr)
library(preprocessCore)           # For quantile normalization
library(qtl2)










### Read in data with corrected sampled IDs
load("~/Desktop/prelim_proteomics_correctedIDs.RData")                # Preliminary pQTL mapping .RData file, by Selcan, /projects/munger-lab/projects/DO_eQTL_comparison/data
load("~/Desktop/do_mesc_corrected_genoprobs_v1.1.RData")              # Created by Duy, swapped genotypes but not sample IDs, /projects/munger-lab/projects/DO_eQTL_comparison/data
load("~/Desktop/Corrected_PBID_SampleMatch_Key.RData")                # Correct matching IDs file by Selcan, /projects/munger-lab/projects/DO_eQTL_comparison/data
ensembl_90_prot  <- readRDS("reference_protein_info_ensembl_90.rds")  # Created from https://github.com/duytpm16/Proteomics-Utils/blob/master/extract_protein_fasta_info.R using Ensembl 94
orig_samples_key <- readRDS('do_esc_sample_ids.rds')










### Get correct ID labels for samples.
#   Lines 35-55 is taken from Selcan (Munger Lab).
matches <- matches %>%
                   filter(!is.na(sampleid.esc_prot)) %>%
                   select(-ends_with(".npc_rna"), -ends_with(".esc_rna")) %>%
                   distinct() %>%
                   arrange(PBID)
colnames(matches) <- gsub("\\.esc_prot","",colnames(matches))
duplicates        <- matches[duplicated(matches$top_muga),]

mixed   <- matches %>%
                   filter(mixup ==T & !endsWith(matches$sampleid,"_rep2") ) %>%
                   mutate(correct_id = ifelse(PBID %in% duplicates$PBID, paste0(PBID,"_B_repB"), paste0(PBID,"_repA"))) %>%
                   mutate(correct_id = ifelse(cor <0.5, sampleid, correct_id))

matches <- matches %>%
                   filter(.,!sampleid %in% mixed$sampleid) %>%
                   mutate(correct_id= sampleid) %>%
                   rbind(mixed) %>%
                   mutate(correct_id=gsub("_rep1","_repA", correct_id)) %>%
                   mutate(correct_id=gsub("_rep2","_B_repB", correct_id)) %>%
                   mutate(correct_id=ifelse(! PBID %in% duplicates$PBID, gsub("_B_repB","_repA", correct_id),correct_id)) %>%
                   arrange(PBID)

#   Correcting a typo error in matches data frame
matches$correct_id[90] <- 'PB360.70_repA'











### Creating covar.factos
#
#   covar.factors: 2 x 5
covar.factors <- data.frame(column.name = c('sex','tmt_label'),
                            display.name = c('Sex','TMT Label'),
                            int.covar = c('factor','NA'),
                            lod.peaks = c('sex_int','NA'),
                            covar.name = c('sex','NA'))














### Changing genotype sample IDs to correct labels
#
#   Making sure genoprobs sample id are consistent for every chromosome
for(i in 2:length(genoprobs)){
    stopifnot(dimnames(genoprobs[[i]])[[1]] == dimnames(genoprobs[[1]])[[1]])
}



#   Changin genoprobs sampled id
gn <- dimnames(genoprobs[[1]])[[1]]
cc <- match(dimnames(genoprobs[[1]])[[1]], matches$sampleid)

for(i in 1:length(genoprobs)){
    dimnames(genoprobs[[i]])[[1]] <- matches$correct_id[cc]
}



#  Should be 1 due to typo
sum(colnames(probs) != dimnames(genoprobs[[1]])[[1]])













### Computing Kinship
K <- calc_kinship(genoprobs, type = 'loco', cores = 10)













### Creating new raw protein expression
#   raw.expr: 7,945 x 195
#   new_raw:  7,853 x 195
new_raw   <- read.delim("~/Desktop/DO_mESC/Proteomics/Original Data/DO_mESC_protein_quant_new.tsv")
rownames(new_raw) <- new_raw$Protein.Id
col_name  <- gsub('.' ,'~', colnames(new_raw), fixed = TRUE)
col_name  <- gsub('X' ,'', col_name, fixed = TRUE)
colnames(new_raw) <- col_name
new_raw <- new_raw[,colnames(new_raw) %in% orig_samples_key$colname]



#   Fixing column names
oc <- match(colnames(new_raw), orig_samples_key$colname)
colnames(new_raw) <- orig_samples_key$full_id[oc]
cc <- match(colnames(new_raw), matches$sampleid)
colnames(new_raw) <-  matches$correct_id[cc]



#   Make 0 NAs, remove samples, and remove proteins with more than 50% NAs
#   Raw: 7,949 x 191
new_raw[new_raw == 0] <- NA
new_raw <- new_raw[order(rownames(new_raw)), order(colnames(new_raw))]
new_raw <- new_raw %>%
                   as.data.frame() %>%
                   select(-PB360.45_B_repB,
                          -PB360.76_B_repB,
                          -PB366.18_B_repB,
                          -PB360.93_B_repB) %>%
                   as.matrix()
rm_prots <- which(rowSums(is.na(new_raw)) > ncol(new_raw) * .5 )
new_raw  <- new_raw[-rm_prots,]













### Create new annots
# 
#   annot    : 9,788 x 3
new_annots <- annot %>%
                    select(protein_id) %>%                                             # Get only protein ID column
                    filter(protein_id %in% rownames(new_raw)) %>%                      # Get proteins in rownames of raw.expr because there is less than 50% NAs for these protein expressions
                    filter(!(grepl('##', protein_id) | grepl('sp',protein_id)))  %>%   # Filter out proteins that start with '##' or 'sp'
                    merge(., ensembl_90_prot, by.x = 'protein_id', all.x = TRUE) %>%   # Merge protein ID with ensembl 94 annotation
                    filter(prot.chr %in% c(1:19,'X')) %>%                              # Get all proteins in 1-19, X chromosomes
                    mutate(gene_id = gsub("\\..*","", gene_id))                        # Remove .* in gene ids
table(new_annots$prot.chr)



#  Get gene annotations
gene_annots <- batchGenes(as.list(unique(new_annots$gene_id)), version = 90, species = 'Mm')
new_annots  <- merge(new_annots, gene_annots[,c('gene_id','chromosome','start','end','strand', 'name','synonyms','symbol')], by = 'gene_id', all.x = TRUE)
stopifnot(new_annots$prot.chr == new_annots$chromosome | new_annots$gene_symbol == new_annots$synonyms)



#  new_annots: 7,858 x 11
new_annots <- new_annots %>%
                         mutate(start  = as.numeric(start) / 1e6,
                                end    = as.numeric(end) / 1e6,
                                middle = (start + end) * 0.5) %>%
                         dplyr::rename(chr=chromosome) %>%
                         select(protein_id, gene_id, transcript_id, symbol, chr, 
                                start, end, middle, strand, gene_biotype, description) %>%
                         arrange(protein_id)



#  Find nearest marker
#
#  new_annots: 7,858 x 12
nearest_marker <- character(length = nrow(new_annots))
for(i in 1:nrow(new_annots)){
  
    sub <- subset(markers, chr == new_annots$chr[i])
    nearest_marker[i] <- sub$marker[which.min(abs(sub$pos - new_annots$start[i]))]
  
}
new_annots$nearest.marker.id <- nearest_marker 
rownames(new_annots) <- new_annots$protein_id















### Filter raw to have  proteins from 1-19, X chromosomes
#   new_raw:  7,949 x 195
#   new_norm: 7,858 x 191
new_raw <- new_raw[rownames(new_raw) %in% new_annots$protein_id,]









### Create new norm dataset by log transformation
#   expr.qnorm: 7,945 x 195
#   new_norm:   7,858 x 191
new_norm <- normalize.quantiles(new_raw)
dimnames(new_norm) <- dimnames(new_raw)








### Creating new rankz proteins
#   expr:      195 x 7,945
#   new_rankz: 191 x 7,858
rankZ <- function (x) {
  x <- rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
  qnorm(x)
}
new_rankz <- apply(new_norm, 1, rankZ)












### Creating new samples
#   covarTidy: 195 x 7
#   new_samples: 191 x 7
id.index <- match(orig_samples_key$full_id, matches$sampleid)
orig_samples_key$full_id  <- matches$correct_id[id.index]
orig_samples_key$array_id <- matches$top_muga[id.index]
orig_samples_key$id       <- matches$PBID[id.index]


new_samples <- orig_samples_key %>%
                                dplyr::rename(mouse.id=full_id,
                                              PBID=id,
                                              MUGA.id=array_id,
                                              orig.name=colname) %>%
                                select(mouse.id,PBID,MUGA.id, sex, tmt_label,tmt_num,orig.name) %>%
                                arrange(mouse.id) %>%
                                filter(!mouse.id %in% c('PB360.45_B_repB','PB360.76_B_repB','PB366.18_B_repB','PB360.93_B_repB'))
rownames(new_samples) <- new_samples$mouse.id


#  Wrong sex labels
new_samples['PB359.70_repA', 'sex']   <- 'F'
new_samples['PB361.72_repA', 'sex']   <- 'F'
new_samples['PB360.52_repA', 'sex']   <- 'F'












### Creating new covariates
#   covar: 195 x 10
#   new_covar: 191 x 10
new_covar <- model.matrix(~sex + tmt_label, data = new_samples)[,-1]
colnames(new_covar)[1] <- 'sex'
rownames(new_covar)    <- new_samples$mouse.id















### QTL Viewer dataset
dataset.esc.proteins <- list(annots        = new_annots,
                             covar         = new_covar,
                             covar.factors = covar.factors,
                             datatype      = 'protein',
                             display.name  = 'DO ESC Proteins',
                             lod.peaks     = list(),
                             norm          = t(new_norm),
                             rankz         = new_rankz,
                             raw           = t(new_raw),
                             samples       = new_samples)

ensembl.version = 90
rm(list = ls()[!ls() %in% c('dataset.esc.proteins','K','genoprobs','map','markers','ensembl.version')])
save.image("~/Desktop/prelim_proteomics_correctedIDs_v3.RData")
