
library(dplyr)
library(pheatmap)

FAM_30 <- read_xlsx("/lustre/tell123/Projects/Postech/Seurat/3_Differential_expression_analysis/Differential_expression.xlsx")
FAM_20 <- read.csv("/lustre/daystar/Postech/Annotation/Seurat.cluster.markers.16.0.4.csv")

# 각 클러스터별 상위 N개 마커만 선택
topN <- 100  # 클러스터당 상위 100개 마커 사용
glist30 <- 30_FAM %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = topN, with_ties = FALSE) %>%
  summarise(genes = list(unique(gene)))

glist20 <- 20_FAM %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = topN, with_ties = FALSE) %>%
  summarise(genes = list(unique(gene)))


  # Jaccard 유사도 계산
jaccard_mat <- outer(seq_len(nrow(glist30)), seq_len(nrow(glist20)), Vectorize(function(i, j){
  g30 <- glist30$genes[[i]]
  g20 <- glist20$genes[[j]]
  length(intersect(g30, g20)) / length(union(g30, g20))
}))

rownames(jaccard_mat) <- paste0("res30_", glist30$cluster)
colnames(jaccard_mat) <- paste0("res20_", glist20$cluster)

# 히트맵 시각화
pheatmap(jaccard_mat, main = "Jaccard similarity (Marker overlap)")
