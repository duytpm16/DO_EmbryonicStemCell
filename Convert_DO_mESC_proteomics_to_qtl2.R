### Options and Libraries
options(stringsAsFactors = FALSE)
library(dplyr)
library(qtl2)
library(qtl2convert)




### Read in .RData
load("~/Desktop/prelim_proteomics.RData")
load("~/Desktop/map_pQTL_rankZ.RData")










### Changing map_data for qtl2
# Getting positions in mB
pos_mb <- unlist(strsplit(map_dat$marker,'_'))
pos_mb <- pos_mb[seq(2,length(pos_mb), by = 2)]
pos_mb <- as.numeric(pos_mb) / 1e6


# Changing the map_data dataframe
new_map     <- map_dat %>%
                       rename(cM  = pos) %>%          # Changing original pos column to cM
                       mutate(pos = pos_mb,           # Making pos the pos_mb vector
                              cM  = as.numeric(cM))   # Making cM column numeric


# Convert the new map dataframe to list
new_map_list <- map_df_to_list(new_map, chr_column = 'chr', 
                               pos_column          = 'pos', 
                               marker_column       = 'marker')









### Converting genoprobs to qtl2 format
new_genoprobs <- probs_doqtl_to_qtl2(probs         = genoprobs, 
                                     map           = new_map, 
                                     chr_column    = 'chr',
                                     pos_column    = 'pos',
                                     marker_column = 'marker')








### Calculate kinship
K <- calc_kinship(probs = new_genoprobs, 
                  type  = 'loco',
                  cores = 0)




### Formatting for qtl viewer
dataset.esc.proteins <- list(annots       = annot,
                             covar        = covar,
                             datatype     = 'protein',
                             display.name = 'DO mESC Proteome',
                             lod.peaks    = list(),
                             rankz        = expr,
                             samples      = as.data.frame(covarTidy))

genoprobs      <- new_genoprobs
map            <- new_map_list
markers        <- new_map
original.peaks <- peaks

rm(gmap, annot, covarTidy, covar, expr, new_genoprobs, new_map, new_map_list,peaks, map_dat, map_pQTL_rankZ.cmd, prelim_proteomics.cmd, pos_mb)






### Save
save.image(file = 'do_mESC_proteomics_qtl_viewer.RData')



