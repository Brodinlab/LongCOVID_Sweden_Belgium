# FlowSOM clustering

# 0. Load packages ----------

library(FlowSOM)
library(flowCore)

# 0. Load data ----------
setwd('~/Documents/LongCovid/Analysis/FlowSOM') 

ff = read.FCS('cytof_49marker_outlier99_arcsinh5_combat_data_wCC.fcs', 
              transformation = FALSE, truncate_max_range = FALSE)
label = read.csv('cytof_49marker_outlier99_arcsinh5_combat_label_wCC.csv', 
                 check.names = F, stringsAsFactors = F) 

# 1. Run FlowSOM ----------
common_marker = colnames(ff@exprs)

# Step 1: generate 30 clusters, separate Neutrophil clusters
fSOM_30 = FlowSOM(ff,
                  compensate = F,
                  transform = F,
                  scale = F,
                  colsToUse = common_marker,
                  nClus = 10,
                  xdim = 5, ydim = 6,
                  seed = 824)

# FlowSOM results summary
FlowSOMmary(fSOM_30, plotFile = 'flowsom_30.pdf')

# FlowSOM clusters median marker expression
MFI_30 = GetClusterMFIs(fSOM_30, colsUsed = T)
clusters_30 = GetClusters(fSOM_30)
freqClusters_30 = data.frame(clusters = clusters_30) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_30 = cbind(MFI_30, freqClusters_30)

# Heatmap of FlowSOM clusters
q99 = quantile(as.matrix(MFI_30), 0.99)
pheatmap::pheatmap(as.matrix(MFI_30),
                   scale = 'none',
                   labels_row = paste(rownames(MFI_30),' (', round(freqClusters_30$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "30 clusters: Median marker expression per cluster")

# Separate Neutrophil clusters
nrow_neutrophils = which(res_30$clusters %in% c(7,17,18,11,1,3,15,6,12,13,8,19,10,2,9,14))
length(nrow_neutrophils) 
label_neutrophils = label[nrow_neutrophils,]  
clusters_neutrophils = res_30$clusters[nrow_neutrophils]
label_neutrophils$lineage = 'Neutrophils'
label_neutrophils$subtype = ''
label_neutrophils$clusters = clusters_neutrophils
data_neutrophils = ff@exprs[nrow_neutrophils,]
res_neutrophils = res_30 %>% dplyr::filter(clusters %in% c(7,17,18,11,1,3,15,6,12,13,8,19,10,2,9,14))

# Remaining cells after removing Neutrophils
data_remaining = ff@exprs[-nrow_neutrophils,]
label_remaining = label[-nrow_neutrophils,]

# Step 2: generating 100 clusters on remaining cells
ff_remove_neutrophils = ff
ff_remove_neutrophils@exprs = data_remaining
fSOM_100 = FlowSOM(ff_remove_neutrophils,
                   compensate = F,
                   transform = F,
                   scale = F,
                   colsToUse = common_marker,
                   nClus = 10,
                   xdim = 10, ydim = 10,
                   seed = 824)

# FlowSOM results summary
FlowSOMmary(fSOM_100, plotFile = 'flowsom_100.pdf')

# FlowSOM clusters median marker expression
MFI_100 = GetClusterMFIs(fSOM_100, colsUsed = T)
clusters_100 = GetClusters(fSOM_100)
freqClusters_100 = data.frame(clusters = clusters_100) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_100 = cbind(MFI_100, freqClusters_100)

#____________________________________________________________________________________________________________

# Heatmap of FlowSOM clusters
q99 = quantile(as.matrix(MFI_100), 0.99)

col_order = c("CD24", "CD16", "CD15","CD55","CD39","Siglec-8","CD9","CD11c","CD64","CD14","CD57","CD123","HLA-DR","CD85j",
                                               "TCRgd","CD95","CD22","CD38","CD25","CD10","CD20","IgD","CD1c","CD161","CD137","CD26","CD28",
                                               "CD29","CD33","CD147","CD141","CD71","CX3CR1","CD99","CD3e","CD27","CD127","CD5","CD4","CD8a",
                                               "CD45RA","CD49d","CD81","CD52","CD7","CD45","CD43","KIR2DL2","CD56")
# Manual column reordering
MFI_100_reordered <- MFI_100[, col_order]

# Plotting
pheatmap::pheatmap(as.matrix(MFI_100_reordered),
                   scale = 'none',
                   labels_row = paste(rownames(MFI_100_reordered), ' (', round(freqClusters_100$percentage, 1), '%', ')', sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   cluster_cols = FALSE,
                   breaks = seq(0, q99, q99/90),
                   main = "100 clusters removing Neutrophils: Median marker expression per cluster")


# Manual annotation
annotation = read.csv2('LC_cluster_annotation.csv', check.names = F, stringsAsFactors = F, sep=';') 
annotationUni <- unique(annotation)
res_100 %<>% left_join(annotationUni)
label_remaining$clusters = clusters_100
label_remaining %<>% left_join(annotationUni)

#Neutrophils
annotationNeut = read.csv2('LC_cluster_annotationNeut.csv', check.names = F, stringsAsFactors = F, sep=';')
res_neutrophils %<>% left_join(annotationNeut)

# 4. Bind remaining cells with neutrophils ----------
sample_info = read.csv('~/Documents/LongCovid/Analysis/FlowSOM/sample_info_metav2.csv', sep=';', stringsAsFactors = FALSE) 

data_bind = rbind(data_neutrophils, data_remaining)
label_bind = rbind(label_neutrophils, label_remaining)
label_bind %<>% mutate(cell_cluster = paste(lineage,'_', subtype, '_', clusters, sep = '')) %>% left_join(sample_info)
all_bind = cbind(data_bind, label_bind) 
res_bind = rbind(res_neutrophils, res_100)
res_bind %<>% mutate(frequency = (n/sum(res_bind$n))*100)
res_bind %<>% mutate(cell_cluster = paste(lineage,'_', subtype, '_', clusters, sep = ''))

# Z-score
res_bind_z <- res_bind %>% mutate_at(c(1:49), funs(c(scale(.))))

# Heatmap of all cell clusters
q99 = quantile(as.matrix(res_bind_z[,1:length(common_marker)]), 0.99)
pheatmap::pheatmap(as.matrix(res_bind_z[,1:length(common_marker)]),
                   scale = 'none',
                   labels_row = paste(res_bind$cell_cluster,' (', round(res_bind$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "Median marker expression per cluster")

# Calculate relative frequency of each celltype in each sample
freq_lineage = label_bind %>%
  group_by(sample_id) %>%
  count(.data$lineage) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  tidyr::spread(lineage, percentage)
freq_subtype = label_bind %>%
  group_by(sample_id) %>%
  count(.data$subtype) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  tidyr::spread(subtype, percentage)
freq_cellcluster = label_bind %>%
  group_by(sample_id) %>%
  count(.data$cell_cluster) %>%
  mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  dplyr::select(-n) %>%
  tidyr::spread(cell_cluster, percentage)
ct_subtype = label_bind %>%
  group_by(sample_id) %>%
  count(.data$subtype) %>%
  mutate(count = .data$n) %>%
  dplyr::select(-n) %>%
  tidyr::spread(subtype, count)
head(ct_subtype)

# PCA of all cell clusters
pca = prcomp(res_bind[,1:length(common_marker)], scale. = T)
pca_df = as.data.frame(pca$x)
pca_df = cbind(pca_df, res_bind[,-(1:length(common_marker))])
ggplot(pca_df, aes(x=PC2, y=PC1)) + 
  geom_point(aes(color=lineage, size = percentage), alpha = 0.8) + 
  theme_bw() 

# 5. Save results ----------
save(all_bind, file = 'flowsom_results.RData')
write.csv(ct_subtype, file = 'ct_subtype.csv', row.names = F)
write.csv(freq_lineage, file = 'freq_lineage.csv', row.names = F)
write.csv(freq_subtype, file = 'freq_subtype.csv', row.names = F)
write.csv(freq_cellcluster, file = 'freq_cellcluster.csv', row.names = F)

# Save median marker expression data for Network plot

write.table(res_bind, file = 'flowsom_clustered_LC.txt', 
            sep="\t", row.names=F, col.names = T, quote = F) #_zscore






