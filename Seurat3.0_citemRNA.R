#!/usr/bin/env Rscript
#Seurat 3.0   Human
###Example: Rscript count.file stat.file barcode.file outdir/pre project_Name

args<-commandArgs(T)
countfile <- args[1]
statfile <- args[2]
barcodef <- args[3]
outdir <- args[4]
pro <- args[5]

## 1.Setup the Seurat Object
library(dplyr)
library(Seurat)
library(ggplot2)

## 2.Load the snRNA dataset
snRNA.data <- read.table(countfile,sep = "\t",row.names = 1,check.names=F,stringsAsFactors=F,header=T)
snRNA.data[is.na(snRNA.data)] <- 0
dim(snRNA.data)

cell_stat <- read.table(statfile,sep = "\t",header=T,check.names=F)
filtered_cell_stat <- cell_stat[which(cell_stat$nUMI>500),]
filtered_cell_stat$CELL_BARCODE <- gsub("-",".",filtered_cell_stat$CELL_BARCODE)
colnames(snRNA.data) <- gsub("-",".",colnames(snRNA.data))
snRNA.data <- snRNA.data[,intersect(colnames(snRNA.data),filtered_cell_stat$CELL_BARCODE)]
Humanbar <- read.table(barcodef,header=F)
Humanbar[,1] <- gsub("-",".",Humanbar[,1])
snRNA.data <- snRNA.data[,intersect(as.vector(Humanbar[,1]),colnames(snRNA.data))]

# Initialize the Seurat object with the raw (non-normalized data).
snRNA <- CreateSeuratObject(counts = snRNA.data, project = "snRNA", min.cells=0, min.features = 200)

## 3.Normalizing the data
# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in snRNA[["RNA"]]@data.

snRNA[["percent.mt"]] <- PercentageFeatureSet(snRNA,pattern = "^MT-")
pdf(paste(outdir,pro,".FeaturePlot.pdf",sep=''))
VlnPlot(snRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(snRNA, features = c("nFeature_RNA", "nCount_RNA"),ncol=2)
snRNA <- subset(snRNA, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA <- NormalizeData(snRNA, normalization.method = "LogNormalize", scale.factor = 10000)
dev.off()

## 4.Identification of highly variable features (feature selection)
snRNA <- FindVariableFeatures(snRNA, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(snRNA), 10)

# 5.Scaling the data
all.genes <- rownames(snRNA)
snRNA <- ScaleData(snRNA, features = all.genes)

## 6.Perform linear dimensional reduction
snRNA <- RunPCA(snRNA, features = VariableFeatures(object = snRNA))

## 7.Determine the ‘dimensionality’ of the dataset
snRNA <- JackStraw(snRNA, num.replicate = 100)
snRNA <- ScoreJackStraw(snRNA, dims = 1:20)
pdf(paste(outdir,pro,".JackStrawPlot.pdf",sep=''),width = 10,height = 8)
JackStrawPlot(snRNA, dims = 1:15)
ElbowPlot(snRNA)
dev.off()

## 8.Cluster the cells
snRNA <- FindNeighbors(snRNA, dims = 1:16)
# We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
# Optimal resolution often increases for larger datasets. 
snRNA <- FindClusters(snRNA, resolution = 0.5)
# The clusters can be found using the Idents function.
# Look at cluster IDs of the first 5 cells
head(Idents(snRNA), 5)

## 9.Run non-linear dimensional reduction (UMAP/tSNE)
snRNA <- RunTSNE(snRNA, dims = 1:16)
snRNA <- RunUMAP(snRNA, dims = 1:16)
#snRNA <- RunTSNE(snRNA, dims = 1:16, method = "FIt-SNE")

pdf(paste(outdir,pro,".tSNE_UMAP.pdf",sep=''))
p1 <- DimPlot(object = snRNA, reduction = "tsne", no.legend = FALSE, pt.size = 0.5) + ggtitle(label = paste0(pro,"_","tSNE"))
p2 <- DimPlot(object = snRNA, reduction = "umap", no.legend = FALSE, pt.size = 0.5) + ggtitle(label = paste0(pro,"_","UMAP"))
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
#CombinePlots(plots = list(p1, p2), legend = "right")
CombinePlots(plots = list(p1, p2), legend = "none")
DimPlot(object = snRNA, reduction = "tsne", no.legend = FALSE, pt.size = 0.5) + ggtitle(label = paste0(pro,"_","tSNE"))
DimPlot(object = snRNA, reduction = "umap", no.legend = FALSE, pt.size = 0.5) + ggtitle(label = paste0(pro,"_","UMAP"))
dev.off()

## ADT-gene plot
genes <- c("CD44","EPCAM","PROM1","CD24","THY1","KRT19","ICAM1")
pdf(paste(outdir,pro,".RNAsignal_UMAP.pdf",sep=''))
DimPlot(object = snRNA, reduction = "umap", no.legend = FALSE, pt.size = 0.5) + ggtitle(label = paste0(pro,"_","UMAP"))
FeaturePlot(snRNA, features = genes,reduction = "umap",cols = c("grey90","yellow","red3"),pt.size=0.005,label.size=0.02)+ggtitle(label = paste0(pro,"_","UMAP"))
dev.off()

## 10.Find markers
snRNA.markers <- FindAllMarkers(snRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
snRNA.markers <- snRNA.markers[snRNA.markers$p_val_adj < 0.05,]
top10 <- snRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top50 <- snRNA.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top100 <- snRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
pdf(paste(outdir,pro,".markers_heatmap.pdf",sep=''))
#FeaturePlot(snRNA, features = top10$gene, min.cutoff = "q05", max.cutoff = "q95", ncol = 4)+ggtitle(label = pro)
DoHeatmap(snRNA, features = top10$gene) + NoLegend()
DoHeatmap(snRNA, features = top50$gene) + NoLegend()+theme(text = element_text(size = 3))
DoHeatmap(snRNA, features = top100$gene) + NoLegend()+theme(text = element_text(size = 1.5))
DoHeatmap(snRNA, features = c("S100A6","VIM","CTSE","KRT20","CD24","AFP","ANPEP","CD44","DLK1","KRT19","EPCAM","PROM1","ICAM1","CD47","LGR5","SO9","POU5F1","NANOG","THY1")) + NoLegend()
dev.off()
write.table(top50,paste0(outdir,pro,"_","top50markers.txt"),sep="\t",col.names=T,row.names=F,quote=F)
write.table(top100,paste0(outdir,pro,"_","top100markers.txt"),sep="\t",col.names=T,row.names=F,quote=F)
write.table(snRNA.markers,paste0(outdir,pro,"_","allmarkers.txt"),sep="\t",col.names=T,row.names=F,quote=F)
saveRDS(snRNA, file = paste(outdir,pro,".snRNA_clustering_CITE.rds",sep=''))
