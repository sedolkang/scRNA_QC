library(Seurat)
library(dplyr)
library(plyr)


data.harmony <-readRDS("/lustre/daystar/Postech/seurat.harmony.16.0.4.RDS")
x.outdir <- x.dirs
  
x.celltype.file <- "/lustre/daystar/Postech/Annotation/cluster_gene_celltype.csv"

##### add celltypes
seurat.data <- data.harmony

x.celltype <- read.table(file = x.celltype.file, header= F ,sep = "\t")
colnames(x.celltype) <- c("seurat_clusters", "celltype_gene")


x.meta <- seurat.data@meta.data
x.meta <- join(x = x.meta, y = x.celltype, by = "seurat_clusters")

seurat.data@meta.data$celltype_gene <- x.meta$celltype_gene


##### Save
saveRDS(object = seurat.data, file = paste0(x.outdir, "/Annotation.Final.16.0.4.RDS"))

#Dimplot
dev.off()
DimPlot(seurat.data, reduction = "umap", group.by = "cluster_with_types", label = TRUE)
