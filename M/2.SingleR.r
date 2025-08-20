x.outdir <- paste0(x.dir, "/Annotation/SingleR")
print(x.outdir)
if (!dir.exists(x.outdir)) { dir.create(x.outdir) }

x.outfile <- paste0(x.outdir, "/SingleR.pred.celltypes")

x.temp.dir <- paste0(x.outdir, "/temp")
if (!dir.exists(x.temp.dir)) { dir.create(x.temp.dir) }

#####	Read data
x.obj <- data.harmony
DefaultAssay(x.obj) = "RNA"

#####	Read mouse reference cells
x.ref.list <- list()

# 마우스 전용 reference dataset만 사용
x.ref.list[["MouseRNAseqData"]] <- MouseRNAseqData()
x.ref.list[["ImmGenData"]] <- ImmGenData()

#####	Run per cluster
x.cluster <- as.character(unique(x.obj@meta.data$seurat_clusters))
x.obj.list <- list()

for (i in 1:length(x.cluster)) {
  x.obj.list[[i]] <- as.SingleCellExperiment(
    subset(x.obj, subset = seurat_clusters == x.cluster[i])
  )
}

x.pred.df <- foreach (i = 1:length(x.cluster), .combine = rbind) %dopar% {
  temp.cluster <- x.cluster[i]
  temp.sce <- x.obj.list[[i]]
  
  temp.pred.df <- data.frame()
  
  for (temp.db in names(x.ref.list)) {
    temp.ref <- x.ref.list[[temp.db]]
    
    temp.genes <- intersect(rownames(temp.sce), rownames(temp.ref))
    
    temp.sce.subset <- temp.sce[temp.genes,]
    temp.ref.subset <- temp.ref[temp.genes,]
    
    temp.sce.subset <- logNormCounts(temp.sce.subset)
    
    #####	Predict main cell types
    temp.pred.main <- SingleR(test = temp.sce.subset,
                              ref = temp.ref.subset,
                              labels = temp.ref.subset$label.main)
    
    temp.pred.main.cell.df <- data.frame(
      Barcodes = rownames(temp.pred.main),
      Ref = temp.db,
      labels = "main",
      CellTypes = temp.pred.main$labels
    )
    write.table(
      x = temp.pred.main.cell.df,
      file = paste0(x.temp.dir, "/", temp.db, ".main.", temp.cluster),
      quote = F, sep = "\t", row.names = F
    )
    
    temp.pred.main.count <- data.frame(
      Clusters = temp.cluster,
      Ref = temp.db,
      labels = "main",
      table(temp.pred.main$labels)
    )
    colnames(temp.pred.main.count) <- c("Clusters", "Ref", "Lables", "CellTypes", "Counts")
    temp.pred.main.count$Ratio <- temp.pred.main.count$Counts /
                                  sum(temp.pred.main.count$Counts)
    temp.pred.main.count <- temp.pred.main.count[
      order(temp.pred.main.count$Ratio, decreasing = T), ]
    
    temp.pred.df <- rbind(temp.pred.df, temp.pred.main.count)
  }
  
  return(temp.pred.df)
}

x.pred.df <- na.omit(x.pred.df)

#####	Write result
write.table(x = x.pred.df, file = x.outfile, quote = F, sep = "\t", row.names = F)













