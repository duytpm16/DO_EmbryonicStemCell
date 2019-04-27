### Code taken from Selcan (Steven Munger lab)




### Options and Libraries
options(stringsAsFactors = FALSE)
library(GenomicRanges)
library(dplyr)
library(assertthat)
library(tidyr)
library(tibble)
library(purrr)








### Calling dplyr functions
rename    <- dplyr::rename
select    <- dplyr::select
summarize <- dplyr::summarize













### Load viewer data
load("~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_proteomics_qtl_viewer_v2.RData")
lod.peaks <- dataset.esc.proteins$lod.peaks$additive












# Prep the objects
markers2 <- markers %>% 
                mutate(chrom = chr, 
                       chromF = factor(chr, levels=c(as.character(1:19), "X")), 
                       n = 1:nrow(.),
                       pos_bp = round(1e6 * pos), 
                       pos_cM = cM, 
                       pos_bp=as.numeric(pos_bp),
                       pos_cM=as.numeric(pos_cM))


peaks <- lod.peaks %>% 
                   filter(!cis) %>%
                   rename(chrom = qtl.chr) %>% 
                   mutate(start = round(qtl.pos * 1e6), end = start) %>%
                   select(chrom, start, end) %>% 
                   GRanges()













### Get counts within windows
x  <- nearest(peaks, markers2 %>% select(chrom, pos_bp) %>% rename(start=pos_bp) %>% mutate(end=start) %>% GRanges())
chrom_markers <- markers2 %>% select(chromF, n) %>% rename(chrom=chromF) %>% group_by(chrom) %>% summarize(start=min(n), end=max(n)) %>% GRanges()
markers_bynum <- markers2 %>% select(chrom, n) %>% rename(start=n) %>% mutate(end=start) %>% GRanges()
windows       <- unlist(slidingWindows(chrom_markers, width=50, step=10))
windows$distant_esc_prot <- countOverlaps(windows, markers_bynum[x])










### Creating tibble of counts in windows
window_counts <- tibble(chrom = as.character(seqnames(windows)),
	                      start = start(windows), 
                        end   = end(windows),
                        distant_esc_prot = windows$distant_esc_prot)
mm <- match(window_counts$start, markers2$n)
m2 <- match(window_counts$end,   markers2$n)
window_counts$pos_cM_start <- markers2$pos_cM[mm]
window_counts$pos_bp_start <- markers2$pos_bp[mm]
window_counts$pos_cM_end   <- markers2$pos_cM[m2]
window_counts$pos_bp_end   <- markers2$pos_bp[m2]
window_counts <- window_counts %>% mutate(midpoint=(pos_cM_end + pos_cM_start)/2, 4)









### Collapse windows
bands.esc.prot <- window_counts %>% select(chrom, starts_with("pos_bp"), starts_with("distant")) %>% filter(distant_esc_prot > 20 )
bands.esc.prot <- bands.esc.prot %>% dplyr::rename(start=pos_bp_start, end=pos_bp_end) %>% GRanges() %>% GenomicRanges::reduce()
bands.esc.prot$distant_esc_prot <- countOverlaps(bands.esc.prot, peaks)










### Final transband count
final.bands.prots <- bands.esc.prot %>% as.data.frame() %>% rename(chr=seqnames) %>% mutate(start = start / 1e6, end = end / 1e6)









### Save
rm(list = ls()[!ls() %in% 'final.bands.prots'])
saveRDS(object = final.bands.prots, file = '~/Desktop/Munger Embryonic Stem Cells/Proteomics/Viewer/munger_esc_proteomics_transband.rds')
