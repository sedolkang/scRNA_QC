

clusters <- unique(data.16@meta.data$seurat_clusters)
out.dir <- "/lustre/daystar/Postech/Results/DEG_clusters/"
data.16 <- JoinLayers(data.16)

for (cluster in clusters) {
  cells_in_cluster <- rownames(data.16@meta.data)[data.16@meta.data$seurat_clusters == cluster]
  data_subset <- subset(data.16, cells = cells_in_cluster)

  markers <- FindMarkers(data_subset, 
                         ident.1 = "WT", ident.2 = "MCE",
                         group.by = "orig.ident",
                         method = "MAST")
  

 write.csv(markers, file = paste0(out.dir, cluster, ".csv"), row.names = TRUE)
}
