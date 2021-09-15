
####################  scRNA  ################################
library(Seurat,lib.loc="/hwfssz1-tmp/ST_MCHRI/STEMCELL/USER/wuliang2/BackupData/wuliang2/biosoftware/anaconda3/lib/R/library")

library(ggplot2,lib.loc="/hwfssz1-tmp/ST_MCHRI/STEMCELL/USER/wuliang2/BackupData/wuliang2/biosoftware/anaconda3/lib/R/library")
library(RColorBrewer)
library(dplyr,lib.loc="/hwfssz1-tmp/ST_MCHRI/STEMCELL/USER/wuliang2/BackupData/wuliang2/biosoftware/anaconda3/lib/R/library")
library(matrixStats,lib.loc="/hwfssz1-tmp/ST_MCHRI/STEMCELL/USER/wuliang2/BackupData/wuliang2/biosoftware/anaconda3/lib/R/library")



setwd("/hwfssz5/ST_PRECISION/TOMCAT/zouxuanxuan/3.cite-seq/03.Seurat/05.fiveCeLines/")

data <- readRDS("./data/cellline5.snRNA_CITE.rds")
Idents(data) <- data@meta.data$orig.ident
data <- RenameIdents(data,`Huh-7`="Huh7",`HepG2`="HepG2",`MHCC97L`="MHCC97L",`MHCC97H`="MHCC97H",`SK-Hep`="SK-HEP-1")
Idents(data) <- factor(Idents(data),levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK-HEP-1"))


data <- RunTSNE(data, dims = 1:15)

######## mRNA
####  tsne
pdf("figures/5celllines_tsne_Set1.pdf")
DimPlot(data, reduction="tsne" ,pt.size = 0.4,cols=brewer.pal(5,"Set1"))+ ggtitle(label = "5celllines_tsne(7311 cells)")+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 14),axis.title = element_text(size = 16))
dev.off()

pdf("figures/5celllines_tsne.pdf")
DimPlot(data, reduction="tsne" ,pt.size = 0.4)+ ggtitle(label = "5celllines_tsne(7311 cells)")+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 14),axis.title = element_text(size = 16))
dev.off()
pdf("figures/5celllines_QC.pdf",width=10,height=8)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, pt.size=0)
dev.off()


############# DEG
snRNA.markers <- FindAllMarkers(data, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25)
snRNA.markers <- snRNA.markers[!grepl("MT-",snRNA.markers$gene),]
top10 <- snRNA.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)

pdf("figures/5celllines_top10DEG.pdf")
DoHeatmap(data, angle=30, features = top10$gene)
dev.off()
saveRDS(data,"cellline5_20210727.rds")



########################################   expression density  ############################################ 

object <- NormalizeData(data, normalization.method="LogNormalize", scale.factor=1000000)

############ 1. data(Mean)
Huh7 <- rowMeans(as.matrix(object@assays$RNA@data[,grepl("Huh7",colnames(object))]))
HepG2 <- rowMeans(as.matrix(object@assays$RNA@data[,grepl("HepG2",colnames(object))]))
C97L <- rowMeans(as.matrix(object@assays$RNA@data[,grepl("MHCC97L",colnames(object))]))
C97H <- rowMeans(as.matrix(object@assays$RNA@data[,grepl("MHCC97H",colnames(object))]))
SK <- rowMeans(as.matrix(object@assays$RNA@data[,grepl("SK-Hep",colnames(object))]))
MeanDa <- rbind(cbind(Gene=names(Huh7),Exp=as.numeric(Huh7),Class="Huh7"),
				cbind(Gene=names(HepG2),Exp=as.numeric(HepG2),Class="HepG2"),
				cbind(Gene=names(C97L),Exp=as.numeric(C97L),Class="MHCC97L"),
				cbind(Gene=names(C97H),Exp=as.numeric(C97H),Class="MHCC97H"),
				cbind(Gene=names(SK),Exp=as.numeric(SK),Class="SK-HEP-1"))
MeanDa <- as.data.frame(MeanDa, stringsAsFactors=F)
MeanDa$Exp <- as.numeric(MeanDa$Exp)
MeanDa$Class <- factor(MeanDa$Class,levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK-HEP-1"))

############ 2. scale.data(Mean)
Huh7 <- rowMeans(as.matrix(object@assays$RNA@scale.data[,grepl("Huh7",colnames(object))]))
HepG2 <- rowMeans(as.matrix(object@assays$RNA@scale.data[,grepl("HepG2",colnames(object))]))
C97L <- rowMeans(as.matrix(object@assays$RNA@scale.data[,grepl("MHCC97L",colnames(object))]))
C97H <- rowMeans(as.matrix(object@assays$RNA@scale.data[,grepl("MHCC97H",colnames(object))]))
SK <- rowMeans(as.matrix(object@assays$RNA@scale.data[,grepl("SK-Hep",colnames(object))]))
MeanScaleDa <- rbind(cbind(Gene=names(Huh7),Exp=as.numeric(Huh7),Class="Huh7"),
				cbind(Gene=names(HepG2),Exp=as.numeric(HepG2),Class="HepG2"),
				cbind(Gene=names(C97L),Exp=as.numeric(C97L),Class="MHCC97L"),
				cbind(Gene=names(C97H),Exp=as.numeric(C97H),Class="MHCC97H"),
				cbind(Gene=names(SK),Exp=as.numeric(SK),Class="SK-HEP-1"))
MeanScaleDa <- as.data.frame(MeanScaleDa, stringsAsFactors=F)
MeanScaleDa$Exp <- as.numeric(MeanScaleDa$Exp)
MeanScaleDa$Class <- factor(MeanScaleDa$Class,levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK-HEP-1"))

pdf("figures/5celllines_density_logCPM.Mean.pdf",width=10, height=8)
ggplot(MeanDa,aes(x=Exp,colour=Class))+geom_density()+labs(x="log(CPM+1)",title="data") + theme_bw()+ theme(plot.title=element_text(hjust=0.5,size=22),axis.text.x=element_text(angle=0,hjust = 0.8,colour="black",family="Times",size=18),axis.text.y=element_text(family="Times",size=18,face="plain"),axis.title=element_text(size=17))
ggplot(MeanScaleDa,aes(x=Exp,colour=Class))+geom_density()+labs(x="Standardized log(CPM+1)",title="Scale.data") + theme_bw()+ theme(plot.title=element_text(hjust=0.5,size=22),axis.text.x=element_text(angle=0,hjust = 0.8,colour="black",family="Times",size=18),axis.text.y=element_text(family="Times",size=18,face="plain"),axis.title=element_text(size=17))
dev.off()



library(matrixStats,lib.loc="/hwfssz1-tmp/ST_MCHRI/STEMCELL/USER/wuliang2/BackupData/wuliang2/biosoftware/anaconda3/lib/R/library")

############ 1. data(Median)
Huh7 <- rowMedians(as.matrix(object@assays$RNA@data[,grepl("Huh7",colnames(object))]))
HepG2 <- rowMedians(as.matrix(object@assays$RNA@data[,grepl("HepG2",colnames(object))]))
C97L <- rowMedians(as.matrix(object@assays$RNA@data[,grepl("MHCC97L",colnames(object))]))
C97H <- rowMedians(as.matrix(object@assays$RNA@data[,grepl("MHCC97H",colnames(object))]))
SK <- rowMedians(as.matrix(object@assays$RNA@data[,grepl("SK-Hep",colnames(object))]))
MedianDa <- rbind(cbind(Gene=names(Huh7),Exp=as.numeric(Huh7),Class="Huh7"),
				cbind(Gene=names(HepG2),Exp=as.numeric(HepG2),Class="HepG2"),
				cbind(Gene=names(C97L),Exp=as.numeric(C97L),Class="MHCC97L"),
				cbind(Gene=names(C97H),Exp=as.numeric(C97H),Class="MHCC97H"),
				cbind(Gene=names(SK),Exp=as.numeric(SK),Class="SK-HEP-1"))
MedianDa <- as.data.frame(MedianDa, stringsAsFactors=F)
MedianDa$Exp <- as.numeric(MedianDa$Exp)
MedianDa$Class <- factor(MedianDa$Class,levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK-HEP-1"))

############ 2. scale.data(Median)
Huh7 <- rowMedians(as.matrix(object@assays$RNA@scale.data[,grepl("Huh7",colnames(object))]))
HepG2 <- rowMedians(as.matrix(object@assays$RNA@scale.data[,grepl("HepG2",colnames(object))]))
C97L <- rowMedians(as.matrix(object@assays$RNA@scale.data[,grepl("MHCC97L",colnames(object))]))
C97H <- rowMedians(as.matrix(object@assays$RNA@scale.data[,grepl("MHCC97H",colnames(object))]))
SK <- rowMedians(as.matrix(object@assays$RNA@scale.data[,grepl("SK-Hep",colnames(object))]))
MedianScaleDa <- rbind(cbind(Gene=names(Huh7),Exp=as.numeric(Huh7),Class="Huh7"),
				cbind(Gene=names(HepG2),Exp=as.numeric(HepG2),Class="HepG2"),
				cbind(Gene=names(C97L),Exp=as.numeric(C97L),Class="MHCC97L"),
				cbind(Gene=names(C97H),Exp=as.numeric(C97H),Class="MHCC97H"),
				cbind(Gene=names(SK),Exp=as.numeric(SK),Class="SK-HEP-1"))
MedianScaleDa <- as.data.frame(MedianScaleDa, stringsAsFactors=F)
MedianScaleDa$Exp <- as.numeric(MedianScaleDa$Exp)
MedianScaleDa$Class <- factor(MedianScaleDa$Class,levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK-HEP-1"))

pdf("figures/5celllines_density_logCPM.Median.pdf",width=10, height=8)
ggplot(MedianDa,aes(x=Exp,colour=Class))+geom_density()+labs(x="log(CPM+1)",title="data") + theme_bw()+ theme(plot.title=element_text(hjust=0.5,size=22),axis.text.x=element_text(angle=0,hjust = 0.8,colour="black",family="Times",size=18),axis.text.y=element_text(family="Times",size=18,face="plain"),axis.title=element_text(size=17))
ggplot(MedianScaleDa,aes(x=Exp,colour=Class))+geom_density()+labs(x="Standardized log(CPM+1)",title="Scale.data") + theme_bw()+ theme(plot.title=element_text(hjust=0.5,size=22),axis.text.x=element_text(angle=0,hjust = 0.8,colour="black",family="Times",size=18),axis.text.y=element_text(family="Times",size=18,face="plain"),axis.title=element_text(size=17))
dev.off()


################# sample cell


p <- list()
for(i in c(1,3,5,7)){
set.seed(i)
Huh7 <- as.matrix(object@assays$RNA@data[,sample(grep("Huh7",colnames(object)),1)])
HepG2 <- as.matrix(object@assays$RNA@data[,sample(grep("HepG2",colnames(object)),1)])
C97L <- as.matrix(object@assays$RNA@data[,sample(grep("MHCC97L",colnames(object)),1)])
C97H <- as.matrix(object@assays$RNA@data[,sample(grep("MHCC97H",colnames(object)),1)])
SK <- as.matrix(object@assays$RNA@data[,sample(grep("SK-Hep",colnames(object)),1)])
MedianDa <- rbind(cbind(Gene=names(Huh7),Exp=as.numeric(Huh7),Class="Huh7"),
				cbind(Gene=names(HepG2),Exp=as.numeric(HepG2),Class="HepG2"),
				cbind(Gene=names(C97L),Exp=as.numeric(C97L),Class="MHCC97L"),
				cbind(Gene=names(C97H),Exp=as.numeric(C97H),Class="MHCC97H"),
				cbind(Gene=names(SK),Exp=as.numeric(SK),Class="SK-HEP-1"))
MedianDa <- as.data.frame(MedianDa, stringsAsFactors=F)
MedianDa$Exp <- as.numeric(MedianDa$Exp)
MedianDa$Class <- factor(MedianDa$Class,levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK-HEP-1"))

############ 2. scale.data(Median)
Huh7 <- as.matrix(object@assays$RNA@scale.data[,sample(grep("Huh7",colnames(object)),1)])
HepG2 <- as.matrix(object@assays$RNA@scale.data[,sample(grep("HepG2",colnames(object)),1)])
C97L <- as.matrix(object@assays$RNA@scale.data[,sample(grep("MHCC97L",colnames(object)),1)])
C97H <- as.matrix(object@assays$RNA@scale.data[,sample(grep("MHCC97H",colnames(object)),1)])
SK <- as.matrix(object@assays$RNA@scale.data[,sample(grep("SK-Hep",colnames(object)),1)])
MedianScaleDa <- rbind(cbind(Gene=names(Huh7),Exp=as.numeric(Huh7),Class="Huh7"),
				cbind(Gene=names(HepG2),Exp=as.numeric(HepG2),Class="HepG2"),
				cbind(Gene=names(C97L),Exp=as.numeric(C97L),Class="MHCC97L"),
				cbind(Gene=names(C97H),Exp=as.numeric(C97H),Class="MHCC97H"),
				cbind(Gene=names(SK),Exp=as.numeric(SK),Class="SK-HEP-1"))
MedianScaleDa <- as.data.frame(MedianScaleDa, stringsAsFactors=F)
MedianScaleDa$Exp <- as.numeric(MedianScaleDa$Exp)
MedianScaleDa$Class <- factor(MedianScaleDa$Class,levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK-HEP-1"))

p[[i]] <- ggplot(MedianDa,aes(x=Exp,colour=Class))+geom_density()+labs(x="log(CPM+1)",title=paste0("data_seed",i)) + theme_bw()+ theme(plot.title=element_text(hjust=0.5,size=22),axis.text.x=element_text(angle=0,hjust = 0.8,colour="black",family="Times",size=18),axis.text.y=element_text(family="Times",size=18,face="plain"),axis.title=element_text(size=17))
p[[(i+1)]] <- ggplot(MedianScaleDa,aes(x=Exp,colour=Class))+geom_density()+labs(x="Standardized log(CPM+1)",title=paste0("Scale.data_seed",i)) + theme_bw()+ theme(plot.title=element_text(hjust=0.5,size=22),axis.text.x=element_text(angle=0,hjust = 0.8,colour="black",family="Times",size=18),axis.text.y=element_text(family="Times",size=18,face="plain"),axis.title=element_text(size=17))
}

pdf("figures/5celllines_density_logCPM.sample.pdf",width=10, height=8)
for(i in 1:length(p)){
print(p[[i]])
}
dev.off()

