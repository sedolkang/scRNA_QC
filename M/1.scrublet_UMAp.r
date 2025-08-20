library(dplyr)
library(Seurat)
library(patchwork)

data_dir <- "/lustre/daystar/CellRanger/agg_MCE_WT/outs/count/filtered_feature_bc_matrix.h5"
WT_dir <- "/lustre/daystar/Postech/Resources/CellRanger/WT/filtered_feature_bc_matrix/"
MCE_dir <- "/lustre/daystar/Postech/Resources/CellRanger/MCE/filtered_feature_bc_matrix/"

WT <- Read10X(data.dir = WT_dir)
MCE <- Read10X(data.dir = MCE_dir)
Aggre_data <- Read10X_h5(data_dir)

wt.seurat <- CreateSeuratObject(counts = WT, project = "WT", min.cells = 3, min.features = 200)
mce.seurat <- CreateSeuratObject(counts = MCE, project = "MCE", min.cells = 3, min.features = 200)

#QC_doublet
scrublet.dir <- "/lustre/daystar/Postech/scrublet/"
scrublet.files <- Sys.glob(paste0(scrublet.dir, "scrublet_results*"))

MCE_scr <- scrublet.files[[1]]
MCE <- fread(MCE_scr)
WT_scr <- scrublet.files[[2]]
WT <- fread(WT_scr)

WT <- subset(WT, predicted_doublet == "FALSE")
wt.seurat@meta.data$barcode <- rownames(wt.seurat)
wt.seurat <- subset(x = wt.seurat, barcode %in% WT$barcode)

MCE <- subset(MCE, predicted_doublet == "FALSE")
mce.seurat@meta.data$barcode <- rownames(mce.seurat@meta.data)
mce.seurat <- subset(x = mce.seurat, barcode %in% MCE$barcode)


colnames(wt.seurat) <- paste0("WT_", colnames(wt.seurat))
colnames(mce.seurat) <- paste0("MCE_", colnames(mce.seurat))

#merge
merged <- merge(mce.seurat, wt.seurat)


#QC_MT
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data<- subset(x = data, subset = percent.mt < 20)

#QC
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data<-subset(data, subset = nFeature_RNA > 200)
data<-subset(data, subset = nCount_RNA > 1000)

VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalize data
data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)

#FindVariableFeatures/Scaling
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("percent.mt"))

#Run PCA
data <- RunPCA(data, features = VariableFeatures(object = data))

saveRDS(data, file='~/bladder/Resource/QC_raw.RDS')


###############
#Clustering
data <- FindNeighbors(data, dims = 1:13)
data <- FindClusters(data, resolution = 0.4) #Higher resolution, smaller cluster

data <- RunUMAP(data, dims = 1:13)

DimPlot(data, reduction = "umap", group.by = "orig.ident")

#Correcting Batch effect
library(harmony)

data.harmony <- RunHarmony(data, group.by.vars = "orig.ident")
data.harmony <- FindNeighbors(data.harmony, dims = 1:30, reduction = 'harmony',graph.name = "harmony_snn")
data.harmony <- FindClusters(data.harmony, resolution = 0.3, graph.name = "harmony_snn")

data.harmony <- RunUMAP(data.harmony, reduction = 'harmony',dims = 1:30)

#plotting
dev.new()
DimPlot(data.harmony, reduction = "umap", group.by = 'cluster_type', label = TRUE)

x.dir <- '/lustre/daystar/Postech/'
saveRDS(data.harmony, file = paste0(x.dir, '/seurat.harmony.13.0.4.RDS'))
