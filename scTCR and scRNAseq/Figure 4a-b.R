library(dplyr)
library(readr)
library(RColorBrewer)
library(circlize)

setwd('scTCR and scRNAseq')
data <- read_csv(file = 'data/scTCR_data_merge.csv.gz')

sub_name <- 'CD8T' # 'CD8T' for Figure 4a, 'CD4T' for Figure 4b 
data_vb <- data %>% 
  mutate(V = sub('\\*.*$', '', TCR_Beta_Delta_V_gene_Dominant))
data_vb_sub <- data_vb %>% filter(cell_type == sub_name)
clone_id_map <- data %>% select(CDR3_concat, clone_id) %>% unique()
exp_id <- data_vb_sub %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat') %>%
  filter(clone_count > 1) %>% select(clone_id) %>% unlist() #id of expanded clones
data_vb_sub <- data_vb_sub %>% mutate(
  clone_id = if_else(clone_id %in% exp_id, paste0('TCR', clone_id), 'unique')) %>%
  mutate(group = case_when(
    grepl('LC', Sample_Name) ~ 'LongCOVID',
    grepl('CP', Sample_Name) ~ 'Convalescent',
    TRUE ~ 'others'
  ))

# Figure 4a/4b ---------
pdf(file = 'figures/4a.pdf', # change the name accordingly
    width = 20, height =10)
par(mfrow = c(1, 2))
df <- data_vb_sub %>% 
  filter(
    clone_id != 'unique',
    group == 'LongCOVID') %>%
  mutate(V = if_else(V=='TRBV11-2', V, 'Others')) %>%
  group_by(V, clone_id) %>%
  tally()
chord_order <- c(df %>% group_by(V) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(V) %>% unlist(use.names = F), 
                 df %>% group_by(clone_id) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(clone_id) %>% unlist(use.names = F)
)
circos.clear()
grid.col = colorRampPalette(brewer.pal(10, "Blues"))(length(chord_order))
names(grid.col) = chord_order
grid.col['Others'] = 'darkgrey'
grid.col['TRBV11-2'] = 'red'
chordDiagram(df, 
             grid.col = grid.col,
             annotationTrack = "grid", 
             order = chord_order,
             small.gap = 0,
             preAllocateTracks = list(track.height = max(strwidth(unique(c(df[[1]],df[[2]]))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
title('LongCOVID')

df <- data_vb_sub %>% 
  filter(
    clone_id != 'unique',
    group == 'Convalescent') %>%
  mutate(V = if_else(V=='TRBV11-2', V, 'Others')) %>%
  group_by(V, clone_id) %>%
  tally()

chord_order <- c(df %>% group_by(V) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(V) %>% unlist(use.names = F), 
                 df %>% group_by(clone_id) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(clone_id) %>% unlist(use.names = F)
)
circos.clear()
grid.col = colorRampPalette(brewer.pal(10, "Blues"))(length(chord_order))
names(grid.col) = chord_order
grid.col['Others'] = 'darkgrey'
grid.col['TRBV11-2'] = 'red'
chordDiagram(df, 
             grid.col = grid.col,
             annotationTrack = "grid", 
             order = chord_order,
             small.gap = 0,
             preAllocateTracks = list(track.height = max(strwidth(unique(c(df[[1]],df[[2]]))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
title('Convalescent')
dev.off()
