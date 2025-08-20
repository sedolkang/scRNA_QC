sample_meta <- seurat.data@meta.data$orig.ident
cluster_meta <- seurat.data@meta.data$celltype_gene

sample <- unique(sample_meta)
cluster <- unique(cluster_meta)

ratio_df <- data.frame(
 cluster = character(),
 sample = character(),
 ratio = numeric()
)

for (C in cluster) {
  total_counts <- sum(cluster_meta ==C)

  for (S in sample){
  counts <- sum(cluster_meta == C & sample_meta == S)
  ratio <- counts / total_counts

  ratio_df <- rbind(ratio_df, data.frame(cluster = C, sample =S, ratio = ratio))

  }
}

write.table(ratio_df, file = '/lustre/daystar/Postech/Results/Final_ratio.csv')

library(ggplot2)
library(scales)

ggplot(ratio_df, aes(x = ratio, y = cluster, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +   # fill로 클러스터별 정규화
  theme_minimal() +
  labs(title = "Sample Ratios per Cluster",
       x = "Proportion",
       y = "Cluster") +
  scale_x_continuous(labels = percent_format()) +    # x축 비율을 %로 표시
  scale_fill_brewer(palette = "Set3")





library(dplyr)

for (S in sample) {
  total_counts <- sum(sample_meta == S)
  
  for (C in cluster) {
    cell_counts <- sum(sample_meta == S & cluster_meta == C)
    ratio <- cell_counts / total_counts
    
    ratio_df <- rbind(ratio_df, data.frame(cluster = C, sample = S, ratio = ratio))
  }
}

ratio_df$cluster <- as.numeric(as.character(ratio_df$cluster))
ratio_df$cluster <- factor(ratio_df$cluster, levels = sort(unique(ratio_df$cluster)))

write.table(ratio_df, file = '/lustre/daystar/Postech/Results/total_ratio.csv')

library(ggplot2)
library(scales)

ggplot(ratio_df, aes(x = ratio, y = cluster, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +   # fill로 클러스터별 정규화
  theme_minimal() +
  labs(title = "Normalized Sample Ratios per Cluster",
       x = "Proportion",
       y = "Cluster") +
  scale_x_continuous(labels = percent_format()) +    # x축 비율을 %로 표시
  scale_fill_brewer(palette = "Set3")

