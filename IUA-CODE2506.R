

library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(reshape2)
#library(SeuratData)
library(patchwork)

#=======================================================================================
IUA1 <- read.table("WBB-20210421_matrix.tsv")
IUA1 <- CreateSeuratObject(counts= IUA1, project = "IUA1", min.cells = 5)
IUA1$stim <- "IUA1"
IUA1[["percent.mt"]] <- PercentageFeatureSet(IUA1, pattern = "^mt-")
summary(IUA1$nCount_RNA)
summary(IUA1$nFeature_RNA)
VlnPlot(IUA1, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA1<- subset(IUA1, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA1, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA1<- NormalizeData(IUA1, verbose = FALSE)
IUA1 <- FindVariableFeatures(IUA1, selection.method = "vst", nfeatures = 2000)
#=======================================================================================

IUA2 <- read.table("WBB-20210422_matrix.tsv")
# Set up control object
IUA2<- CreateSeuratObject(counts= IUA2, project = "IUA2", min.cells = 5)
IUA2$stim <- "IUA2"
IUA2[["percent.mt"]] <- PercentageFeatureSet(IUA2, pattern = "^MT-")
summary(IUA2$nCount_RNA)
summary(IUA2$nFeature_RNA)
VlnPlot(IUA2, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA2<- subset(IUA2, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA2, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA2<- NormalizeData(IUA2, verbose = FALSE)
IUA2<- FindVariableFeatures(IUA2, selection.method = "vst", nfeatures = 2000)

#############################################################################
IUA3 <- read.table("WBB-20210423_matrix.tsv")
IUA3 <- CreateSeuratObject(counts= IUA3, project = "IUA3", min.cells = 5)
IUA3$stim <- "IUA3"
IUA3[["percent.mt"]] <- PercentageFeatureSet(IUA3, pattern = "^MT-")
summary(IUA3$nCount_RNA)
summary(IUA3$nFeature_RNA)
VlnPlot(IUA3, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA3<- subset(IUA3, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA3, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA3<- NormalizeData(IUA3, verbose = FALSE)
IUA3<- FindVariableFeatures(IUA3, selection.method = "vst", nfeatures = 2000)


#IUA4=======================================================================================
IUA4 <- CreateSeuratObject(Read10X('E-20210819'),"IUA4")
IUA4$stim <- "IUA4"
IUA4[["percent.mt"]] <- PercentageFeatureSet(IUA4, pattern = "^mt-")
summary(IUA4$nCount_RNA)
summary(IUA4$nFeature_RNA)
VlnPlot(IUA4, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA4<- subset(IUA4, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA4, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA4<- NormalizeData(IUA4, verbose = FALSE)
IUA4<- FindVariableFeatures(IUA4, selection.method = "vst", nfeatures = 2000)


#IUA5=======================================================================================
IUA5<- CreateSeuratObject(Read10X('E-20210831'),"IUA5")
IUA5$stim <- "IUA5"
IUA5[["percent.mt"]] <- PercentageFeatureSet(IUA5, pattern = "^mt-")
summary(IUA5$nCount_RNA)
summary(IUA5$nFeature_RNA)
VlnPlot(IUA5, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA5<- subset(IUA5, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA5, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA5<- NormalizeData(IUA5, verbose = FALSE)
IUA5<- FindVariableFeatures(IUA5, selection.method = "vst", nfeatures = 2000)


#IUA6=======================================================================================
IUA6<- CreateSeuratObject(Read10X('E-20211115_IUA'),"IUA6")
IUA6$stim <- "IUA6"
IUA6[["percent.mt"]] <- PercentageFeatureSet(IUA6, pattern = "^mt-")
summary(IUA6$nCount_RNA)
summary(IUA6$nFeature_RNA)
VlnPlot(IUA6, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA6<- subset(IUA6, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA6, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA6<- NormalizeData(IUA6, verbose = FALSE)
IUA6<- FindVariableFeatures(IUA6, selection.method = "vst", nfeatures = 2000)


#IUA7=======================================================================================
IUA7<- CreateSeuratObject(Read10X('E-20211116'),"IUA7")
IUA7$stim <- "IUA7"
IUA7[["percent.mt"]] <- PercentageFeatureSet(IUA7 , pattern = "^mt-")
summary(IUA7$nCount_RNA)
summary(IUA7$nFeature_RNA)
VlnPlot(IUA7, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA7<- subset(IUA7, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA7, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA7<- NormalizeData(IUA7, verbose = FALSE)
IUA7<- FindVariableFeatures(IUA7, selection.method = "vst", nfeatures = 2000)


#IUA8=======================================================================================
IUA8<- CreateSeuratObject(Read10X('E-20211117_IUA'),"IUA8")
IUA8$stim <- "IUA8"
IUA8[["percent.mt"]] <- PercentageFeatureSet(IUA8, pattern = "^mt-")
summary(IUA8$nCount_RNA)
summary(IUA8$nFeature_RNA)
VlnPlot(IUA8, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA8<- subset(IUA8, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA8, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA8<- NormalizeData(IUA8, verbose = FALSE)
IUA8<- FindVariableFeatures(IUA8, selection.method = "vst", nfeatures = 2000)

#IUA9=======================================================================================
IUA9<- CreateSeuratObject(Read10X('E-20211201'),"IUA9")
IUA9$stim <- "IUA9"
IUA9[["percent.mt"]] <- PercentageFeatureSet(IUA9, pattern = "^mt-")
summary(IUA9$nCount_RNA)
summary(IUA9$nFeature_RNA)
VlnPlot(IUA9, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA9<- subset(IUA9, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(IUA9, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
IUA9<- NormalizeData(IUA9, verbose = FALSE)
IUA9<- FindVariableFeatures(IUA9, selection.method = "vst", nfeatures = 2000)

NE.anchors <- FindIntegrationAnchors(object.list = list(IUA1,IUA2,IUA3,IUA4,IUA5,IUA6,IUA7,IUA8,IUA9), dims = 1:20)
NE.combined <- IntegrateData(anchorset = NE.anchors, dims = 1:20)

###############################################################################################
DefaultAssay(NE.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
NE.combined <- ScaleData(NE.combined, verbose = FALSE)
NE.combined <- RunPCA(NE.combined, npcs = 30, verbose = FALSE)
ElbowPlot(NE.combined)
NE.combined$stim <- "IUA"

###############################################################################################

##load NP.Rata#############################################################################################

NP <- CreateSeuratObject(counts= NP.data, project = "NP", min.cells = 5)
NP$stim <- "NP"
NP[["percent.mt"]] <- PercentageFeatureSet(NP, pattern = "^MT-")
summary(NP$nCount_RNA)
summary(NP$nFeature_RNA)
VlnPlot(NP, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
NP <- subset(NP, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(NP, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
NP<- NormalizeData(NP, verbose = FALSE)
NP <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

##load NEN.combined#############################################################################################

#normal1=======================================================================================
mat <- Matrix::readMM("4861STDY7387181/matrix.mtx")
barcodes <- read.table("4861STDY7387181/4861STDY7387181_cells.tsv",header = T,sep = '\t',row.names = 1)
genes <- read.table("4861STDY7387181/genes.tsv", header= F, sep="\t",row.names = NULL)
colnames(mat) <- rownames(barcodes)
rownames(mat) <- genes$V1
Normal1 <- CreateSeuratObject(mat)
Normal1$stim <- "Normal1"
Normal1[["percent.mt"]] <- PercentageFeatureSet(Normal1, pattern = "^mt-")
summary(Normal1$nCount_RNA)
summary(Normal1$nFeature_RNA)
VlnPlot(Normal1, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal1<- subset(Normal1, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(Normal1, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal1<- NormalizeData(Normal1, verbose = FALSE)
Normal1<- FindVariableFeatures(Normal1, selection.method = "vst", nfeatures = 2000)

#normal2=======================================================================================
mat2 <- Matrix::readMM("4861STDY7387182/matrix.mtx")
barcodes <- read.table("4861STDY7387182/4861STDY7387182_cells.tsv",header = T,sep = '\t',row.names = 1)
genes <- read.table("4861STDY7387182/genes.tsv", header= F, sep="\t",row.names = NULL)
colnames(mat2) <- rownames(barcodes)
rownames(mat2) <- genes$V1
Normal2 <- CreateSeuratObject(mat2)
Normal2$stim <- "Normal2"
Normal2[["percent.mt"]] <- PercentageFeatureSet(Normal2, pattern = "^mt-")
summary(Normal2$nCount_RNA)
summary(Normal2$nFeature_RNA)
VlnPlot(Normal2, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal2<- subset(Normal2, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(Normal2, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal2<- NormalizeData(Normal2, verbose = FALSE)
Normal2<- FindVariableFeatures(Normal2, selection.method = "vst", nfeatures = 2000)

#normal9=======================================================================================
mat9 <- Matrix::readMM("MRCEndo8715415/matrix.mtx")
barcodes <- read.table("MRCEndo8715415/MRCEndo8715415_cells.tsv",header = T,sep = '\t',row.names = 1)
genes <- read.table("MRCEndo8715415/genes.tsv", header= F, sep="\t",row.names = NULL)
colnames(mat9) <- rownames(barcodes)
rownames(mat9) <- genes$V1
Normal9 <- CreateSeuratObject(mat9)
Normal9$stim <- "Normal9"
Normal9[["percent.mt"]] <- PercentageFeatureSet(Normal9, pattern = "^mt-")
summary(Normal9$nCount_RNA)
summary(Normal9$nFeature_RNA)
VlnPlot(Normal9, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal9<- subset(Normal9, subset = percent.mt < 50 & nCount_RNA < 60000 & nFeature_RNA <7500)
VlnPlot(Normal9, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal9<- NormalizeData(Normal9, verbose = FALSE)
Normal9<- FindVariableFeatures(Normal9, selection.method = "vst", nfeatures = 2000)

#normal10=======================================================================================
mat10 <- Matrix::readMM("MRCEndo8715416/matrix.mtx")
barcodes <- read.table("MRCEndo8715416/MRCEndo8715416_cells.tsv",header = T,sep = '\t',row.names = 1)
genes <- read.table("MRCEndo8715416/genes.tsv", header= F, sep="\t",row.names = NULL)
colnames(mat10) <- rownames(barcodes)
rownames(mat10) <- genes$V1
Normal10 <- CreateSeuratObject(mat10)
Normal10$stim <- "Normal10"
Normal10[["percent.mt"]] <- PercentageFeatureSet(Normal5, pattern = "^mt-")
summary(Normal10$nCount_RNA)
summary(Normal10$nFeature_RNA)
VlnPlot(Normal10, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal10<- subset(Normal10, subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(Normal10, features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
Normal10<- NormalizeData(Normal10, verbose = FALSE)
Normal10<- FindVariableFeatures(Normal10, selection.method = "vst", nfeatures = 2000)

#Integration=======================================================================================
NEN.anchors <- FindIntegrationAnchors(object.list = list(Normal1,Normal2,Normal9,Normal10), dims = 1:20)
save(NEN.anchors, file = "NEN.anchors.Rda")
NEN.combined <- IntegrateData(anchorset = NEN.anchors, dims = 1:20)
save(NEN.combined, file = "NEN_combined.Rda")
#=======================================================================================
#=======================================================================================

##load control-NM#############################################################################################
control <- CreateSeuratObject(counts= GSE111976_ct_endo_10x, project = "P1", min.cells = 5)
control $stim <- "control "
control [["percent.mt"]] <- PercentageFeatureSet(control , pattern = "^mt-")
summary(control $nCount_RNA)
summary(control $nFeature_RNA)
VlnPlot(control , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
control <- subset(control , subset = percent.mt < 10 & nCount_RNA < 40000 & nFeature_RNA <6000)
VlnPlot(control , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
control  <- NormalizeData(control , verbose = FALSE)
control  <- FindVariableFeatures(control , selection.method = "vst", nfeatures = 2000)
#=======================================================================================
#Integration=======================================================================================
IUA.anchors<- FindIntegrationAnchors(object.list = list(NE.combined,NP,control,NEN.combined), dims = 1:20)

All_IUA_proliferative <- IntegrateData(anchorset = IUA.anchors, dims = 1:20)

###############################################################################################
#scRNA-SEQ-annotation#########################################################################################

cellinfo<-read.table('All_IUA_proliferative_cell_annotation.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
All_IUA_proliferative<- AddMetaData(object = All_IUA_proliferative, metadata = cellinfo)


save(All_IUA_proliferative,file = "All_IUA_proliferative.Rdata")

ElbowPlot(All_IUA_proliferative)

All_IUA_proliferative <- RunUMAP(All_IUA_proliferative, reduction = "pca", dims = 1:8)
All_IUA_proliferative <- FindNeighbors(All_IUA_proliferative, reduction = "pca", dims = 1:8)
All_IUA_proliferative <- FindClusters(All_IUA_proliferative, resolution = 0.3)
# Visualization-UMAP
p1 <- DimPlot(All_IUA_proliferative, reduction = "umap", group.by = "stim" ,cols = c("#DC143C","#1E90FF"),pt.size = 0.1,raster =F)
#p2 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE, pt.size = 0.1,raster =F)
#p3 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE,group.by = "cell_type" ,pt.size = 0.1,raster =F)
p4 <- DimPlot(All_IUA_proliferative, reduction = "umap", group.by = "datasets" ,cols = c("#32CD32","#1E90FF","#DC143C"),pt.size = 0.1,raster =F)
p5 <- DimPlot(All_IUA_proliferative, reduction = "umap", group.by = "cell_types" ,cols = c("#20B2AA","#0000FF","#FFA500","#9370DB","#DC143C"),pt.size = 0.1,raster =F)
p6 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = F,group.by = "cell_cluster" ,cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                                                                                                        "#9370DB","#98FB98","#F08080","#1E90FF",
                                                                                                        "#7CFC00","#FFFF00","#808000","#FF00FF",
                                                                                                        "#FA8072","#7B68EE","#9400D3","#800080",
                                                                                                        "#A0522D","#D2B48C","#D2691E","#87CEEB",
                                                                                                        "#40E0D0","#5F9EA0","#FF1493","#0000CD",
                                                                                                        "#FFE4B5","#8A2BE2","#228B22"),pt.size = 0.1,raster =F)
#original_ID
p7 <- DimPlot(All_IUA_proliferative, reduction = "umap", group.by = "original_ID" ,cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                                                                                            "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                                                                                            "#808000","#FF00FF","#40E0D0","#7B68EE",
                                                                                            "#9400D3","#800080","#A0522D","#D2B48C"),pt.size = 0.1,raster =F)
plot_grid(p7)
plot_grid(p1)
#plot_grid(p2)
#plot_grid(p3)
plot_grid(p4)
plot_grid(p6)
plot_grid(p5)



markers.to.plot <- c("COL1A1","ECM1","ACTA2","PTPRC","CD68","CD3E","KRT8","KRT18","EPCAM","VWF","CLDN5")
DotPlot(All_IUA_proliferative, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_types5")+
  RotatedAxis()

###############################################
All_IUA_proliferative <- SetIdent(All_IUA_proliferative,cells = names(All_IUA_proliferative@meta.data$cell_types),value = All_IUA_proliferative@meta.data$cell_types)
STROMA_proliferative <- subset(All_IUA_proliferative,idents = c("STROMA"))
ENDOTHELIA_proliferative <- subset(All_IUA_proliferative,idents = c("ENDOTHELIA"))
EPITHELIA_proliferative <- subset(All_IUA_proliferative,idents = c("EPITHELIA"))
IMMUNE_proliferative <- subset(All_IUA_proliferative,idents = c("IMMUNE"))

#stroma==========================================================================
save(STROMA_proliferative,file = "STROMA_proliferative.Rdata")

STROMA_proliferative <- RunUMAP(STROMA_proliferative, reduction = "pca", dims = 1:11)
STROMA_proliferative <- FindNeighbors(STROMA_proliferative, reduction = "pca", dims = 1:11)
STROMA_proliferative <- FindClusters(STROMA_proliferative, resolution = 0.2)
# Visualization-UMAP

p1 <- DimPlot(STROMA_proliferative, reduction = "umap", group.by = "stim" ,cols = c("#DC143C","#1E90FF"),pt.size = 0.5)
p4 <- DimPlot(STROMA_proliferative, reduction = "umap", group.by = "original_ID" ,
              cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                       "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                       "#808000","#FF00FF","#40E0D0","#7B68EE",
                       "#9400D3","#800080","#A0522D","#D2B48C"),pt.size = 0.5)
p2 <- DimPlot(STROMA_proliferative, reduction = "umap", label = TRUE, pt.size = 0.5)
#p3 <- DimPlot(STROMA_proliferative, reduction = "umap", group.by = "cell_types" ,label = F,cols = c("#0000FF","#20B2AA","#FFA500","#DC143C",
#                                                                                             "#9370DB","#228B22"),pt.size = 0.5)

p5 <- DimPlot(STROMA_proliferative, reduction = "umap", group.by = "cell_cluster" ,label = F,cols = c("#0000FF","#20B2AA","#FFA500","#DC143C",
                                                                                             "#9370DB","#228B22"),pt.size = 0.5)


plot_grid(p1)
plot_grid(p4)
plot_grid(p2)
#plot_grid(p3)
plot_grid(p5)
###stroma2###############################################################
cellinfo3<-read.table('STROMA_proliferative-cell-cluster2.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
STROMA_proliferative<- AddMetaData(object = STROMA_proliferative, metadata = cellinfo3)


STROMA_proliferative <- RunUMAP(STROMA_proliferative, reduction = "pca", dims = 1:11)
STROMA_proliferative <- FindNeighbors(STROMA_proliferative, reduction = "pca", dims = 1:11)
STROMA_proliferative <- FindClusters(STROMA_proliferative, resolution = 0.2)
# Visualization-UMAP

p1 <- DimPlot(STROMA_proliferative, reduction = "umap", group.by = "stim" ,cols = c("#DC143C","#1E90FF"),pt.size = 0.5)
p4 <- DimPlot(STROMA_proliferative, reduction = "umap", group.by = "original_ID" ,
              cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                       "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                       "#808000","#FF00FF","#40E0D0","#7B68EE",
                       "#9400D3","#800080","#A0522D","#D2B48C"),pt.size = 0.5)
p2 <- DimPlot(STROMA_proliferative, reduction = "umap", label = TRUE, pt.size = 0.5)
p3 <- DimPlot(STROMA_proliferative, reduction = "umap", group.by = "cell_cluster" ,cols = c("#0000FF","#DC143C","#20B2AA","#FFA500",
                                                                                            "#9370DB","#228B22"),pt.size = 0.5)

plot_grid(p1)
plot_grid(p4)
plot_grid(p2)
plot_grid(p3)
###SCENIC########################################################################
library(Seurat) 
library(SCENIC)

Idents(object = STROMA_proliferative) <- 'cell_cluster'
# Downsample the number of cells per identity class
STROMA_proliferative1 <- subset(x = STROMA_proliferative, downsample = 100)

## Load data表达矩阵加上表型信息
exprMat  <-  as.matrix(STROMA_proliferative1@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  STROMA_proliferative1@meta.data
save(exprMat,cellInfo,file = "STROMA_proliferative1_scenic.Rdata")

### Initialize settings
org<-"hgnc"
dbDir<-"cisTarget_databases"

myDatasetTitle<-"SCENIC analysis on STROMA_proliferative"
data("defaultDbNames")
dbs<-defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

motifAnnotations_hgnc = motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
#创建一个对象：scenicOptions，后续数据会保存在当前目录的int文件夹
#输出结果会存储在output文件夹
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#将表达矩阵中未在cisTarget_databases收录的基因去除
library(RcisTarget)
dbFilePath<-getDatabases(scenicOptions)[[1]]
motifRanking<-importRankings(dbFilePath)
geneInDatabase<-colnames(getRanking(motifRanking))
genelefet_minCells<-rownames(exprMat)
length(genelefet_minCells)
genelefet_minCells_inDatabases<-genelefet_minCells[which(genelefet_minCells %in% geneInDatabase)]
length(genelefet_minCells_inDatabases)
genesKept<-genelefet_minCells_inDatabases
exprMat_filter<-exprMat[genesKept,]
dim(exprMat)
dim(exprMat_filter)

#计算相关性并对数据进行log处理
### Co-expression network
class(exprMat_filter)
runCorrelation(as.matrix(exprMat_filter), scenicOptions)
exprMat_filter <- log2(exprMat_filter+1) 
library(GENIE3)
runGenie3(as.matrix(exprMat_filter), scenicOptions)
save(exprMat_filter,scenicOptions,file="input_GENIE3_data.Rdata")


### Build the GRN
#load(stroma_scenic.Rdata)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
#closeAllConnections()
#scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget")) # Toy run settings
#library(doParallel)
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#结果可视化
library(AUCell)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)


#Regulators for known cell types or clusters

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cell_cluster3),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, col = col,name="Regulon activity")

###########################################################################
#stroma_proliferative-slingshot=============
library(RColorBrewer)
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(mclust)
library(tradeSeq)

library(slingshot)
library(Seurat)
library(devtools)
library(cowplot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tradeSeq)
library(RColorBrewer)
library(DelayedMatrixStats)
library(scales)
library(paletteer) 
library(viridis)
library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(reshape2)
#library(SeuratData)
library(patchwork)

MD_2_0<- as.SingleCellExperiment(STROMA_proliferative)
#save(MD_3, file = "MD_3_slingshot.Rda")
sce_slingshot1 <- slingshot(MD_2_0,      
                            reducedDim = 'UMAP',  
                            start.clus = 'PROLIFERATIVE_STROMA',
                            clusterLabels = STROMA_proliferative$cell_cluster,
                            end.clus =c('INFLAMMATORY_STROMA','MYOGENIC_STROMA'),
                            stretch = 2,
                            approx_points = 150)
SlingshotDataSet(sce_slingshot1) 


#定义颜色
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }}
#设置颜色并可视化
cell_colors <- cell_pal(sce_slingshot1$cell_cluster, brewer_pal("qual", "Set2"))
plot(reducedDims(sce_slingshot1)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.8)

#plot(reducedDims(sce_slingshot1)$UMAP, col = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB"), pch=16, asp = 1, cex = 0.8)

lines(SlingshotDataSet(sce_slingshot1),  lwd=2, col='black')

#计算celltype坐标位置，用于图中标记
celltype_label <- STROMA_proliferative@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(celltype = STROMA_proliferative@meta.data$cell_cluster) %>%
  group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP_1),
            UMAP2 = median(UMAP_2))

for (i in 1:8) {
  text(celltype_label$celltype[i], x=celltype_label$UMAP1[i]-1, y=celltype_label$UMAP2[i])
}

#也可以把slingshot结果放到单细胞seurat对象中-Metadata
pseudotime = slingPseudotime(sce_slingshot1)%>% as.data.frame() 
Lineages = colnames(pseudotime)
STROMA_proliferative <- AddMetaData(object = STROMA_proliferative,metadata = pseudotime)


#stroma_proliferative-slingshot-tradeseq=============

STROMA_proliferative_slingshot<- CreateSeuratObject(counts = STROMA_proliferative@assays$RNA@counts,meta.data=STROMA_proliferative@meta.data,project = "STROMA_proliferative_slingshot")
STROMA_proliferative_slingshot <- NormalizeData(STROMA_proliferative_slingshot, normalization.method = "LogNormalize", scale.factor = 10000)
STROMA_proliferative_slingshot <- FindVariableFeatures(STROMA_proliferative_slingshot, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(STROMA_proliferative_slingshot)
scale_gene<-VariableFeatures(STROMA_proliferative_slingshot)
#counts<-CE_MD_slingshot@assays$RNA@counts
#counts<counts[scale_gene,]
STROMA_proliferative_slingshot_HVG<-STROMA_proliferative_slingshot[scale_gene,]

#基因随轨迹表达变化
slingsce<-SlingshotDataSet(sce_slingshot1)
pseudotimeED <- as.matrix(unlist(slingPseudotime(slingsce, na = FALSE)))
cellWeightsED <-  as.matrix(unlist(slingCurveWeights(slingsce)))
#counts<-sce_slingshot1@assays@data@listData
#counts<- as.matrix(unlist(sce_slingshot1@assays@data@listData))
#counts<- as.matrix(unlist(MD_2@assays$RNA@counts))
counts<- as.matrix(unlist(STROMA_proliferative_slingshot_HVG@assays$RNA@counts))

#write.table(EPITHELIA_proliferative$cell_cluster,file="EPITHELIA_proliferative_cell_cluster.txt")

#write.table(pseudotime$Lineage1,file="Lineage1.txt")

#这里用HVG基因
system.time({ 
  sce_STROMA_proliferative_slinghot <- fitGAM(counts = counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T)
})


save(sce_STROMA_proliferative_slinghot, file = "sce_STROMA_proliferative_slinghot.Rda")


#探索基因表达与拟时序的相关性
assoRes <- associationTest(sce_STROMA_proliferative_slinghot)
head(assoRes)
write.table(assoRes,file="STROMA_proliferative_Pseudotime_association_GENE.txt")


#寻找与起止点相关性最高的基因
startRes <- startVsEndTest(sce_STROMA_proliferative_slinghot)
head(startRes)
write.table(startRes,file="STROMA_proliferative_Pseudotime_startRes_GENE.txt")

# 按相关性进行排序
oStart <- order(startRes$waldStat, decreasing = TRUE)

# 挑选相关性最强的基因，并可视化
sigGeneStart <- names(sce_STROMA_proliferative_slinghot)[oStart[1]]
plotSmoothers(sce_STROMA_proliferative_slinghot, counts, gene = "TPPP3")

#也可以用UMAP图展示
plotGeneCount(slingsce, counts, gene = "TPPP3")
plotGeneCount(slingsce, counts, gene = c("COL1A1","SFRP4","IGF1","ACTA2","NFKBIZ","MKI67","HOXA10","PDGFRB","ADIRF","CNN1","TAGLN"))
#=========================================================================

###EPITHELIA###########################################
save(EPITHELIA_proliferative, file = "EPITHELIA_proliferative.Rda")

cellinfo3<-read.table('EPITHELIA_proliferative-cluster.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
EPITHELIA_proliferative<- AddMetaData(object = EPITHELIA_proliferative, metadata = cellinfo3)

ElbowPlot(EPITHELIA_proliferative)

EPITHELIA_proliferative <- RunUMAP(EPITHELIA_proliferative, reduction = "pca", dims = 1:9)
EPITHELIA_proliferative <- FindNeighbors(EPITHELIA_proliferative, reduction = "pca", dims = 1:9)
EPITHELIA_proliferative <- FindClusters(EPITHELIA_proliferative, resolution = 0.2)
# Visualization-UMAP
p1 <- DimPlot(EPITHELIA_proliferative, reduction = "umap", group.by = "stim" ,cols = c("#DC143C","#1E90FF"),pt.size = 0.5)
p2 <- DimPlot(EPITHELIA_proliferative, reduction = "umap", group.by = "original_ID" ,cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                                                                                              "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                                                                                              "#808000","#FF00FF","#40E0D0","#7B68EE",
                                                                                              "#9400D3","#800080","#A0522D","#D2B48C"),pt.size = 0.5)

p3 <- DimPlot(EPITHELIA_proliferative, reduction = "umap", label = TRUE, pt.size = 0.5)
p4 <- DimPlot(EPITHELIA_proliferative, reduction = "umap", group.by = "cell_cluster" ,cols = c("#A0522D","#DC143C","#20B2AA","#FFA500",
                                                                                                "#9370DB","#9400D3","#0000FF"),label = T,pt.size = 0.5)

plot_grid(p1)
plot_grid(p2)
plot_grid(p3)

###SCENIC########################################################################
setwd("E:/IUA-revision/scenic_epi")
library(Seurat) 
library(SCENIC)

Idents(object = EPITHELIA_proliferative) <- 'cell_cluster'
# Downsample the number of cells per identity class
EPITHELIA_proliferative1 <- subset(x = EPITHELIA_proliferative, downsample = 100)

## Load data表达矩阵加上表型信息
exprMat  <-  as.matrix(EPITHELIA_proliferative1@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  EPITHELIA_proliferative1@meta.data
save(exprMat,cellInfo,file = "EPITHELIA_proliferative1_scenic.Rdata")

### Initialize settings
org<-"hgnc"
dbDir<-"cisTarget_databases"

myDatasetTitle<-"SCENIC analysis on EPITHELIA_proliferative"
data("defaultDbNames")
dbs<-defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

motifAnnotations_hgnc = motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
#创建一个对象：scenicOptions，后续数据会保存在当前目录的int文件夹
#输出结果会存储在output文件夹
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#将表达矩阵中未在cisTarget_databases收录的基因去除
library(RcisTarget)
dbFilePath<-getDatabases(scenicOptions)[[1]]
motifRanking<-importRankings(dbFilePath)
geneInDatabase<-colnames(getRanking(motifRanking))
genelefet_minCells<-rownames(exprMat)
length(genelefet_minCells)
genelefet_minCells_inDatabases<-genelefet_minCells[which(genelefet_minCells %in% geneInDatabase)]
length(genelefet_minCells_inDatabases)
genesKept<-genelefet_minCells_inDatabases
exprMat_filter<-exprMat[genesKept,]
dim(exprMat)
dim(exprMat_filter)

#计算相关性并对数据进行log处理
### Co-expression network
class(exprMat_filter)
runCorrelation(as.matrix(exprMat_filter), scenicOptions)
exprMat_filter <- log2(exprMat_filter+1) 
library(GENIE3)
runGenie3(as.matrix(exprMat_filter), scenicOptions)
save(exprMat_filter,scenicOptions,file="input_GENIE3_data.Rdata")


### Build the GRN
#load(stroma_scenic.Rdata)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
#closeAllConnections()
#scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget")) # Toy run settings
#library(doParallel)
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#结果可视化
library(AUCell)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)


#Regulators for known cell types or clusters

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cell_cluster3),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, col = col,name="Regulon activity")


###########################################################################
#EPITHELIA_proliferative-slingshot=============

library(RColorBrewer)
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(mclust)
library(tradeSeq)

library(slingshot)
library(Seurat)
library(devtools)
library(cowplot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tradeSeq)
library(RColorBrewer)
library(DelayedMatrixStats)
library(scales)
library(paletteer) 
library(viridis)
library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(reshape2)
#library(SeuratData)
library(patchwork)

MD_2_0<- as.SingleCellExperiment(EPITHELIA_proliferative)
#save(MD_3, file = "MD_3_slingshot.Rda")
sce_slingshot1 <- slingshot(MD_2_0,      #输入单细胞对象
                            reducedDim = 'UMAP',  #降维方式
                            start.clus = 'PROLIFERATIVE_EPI',       #轨迹起点,也可以不定义
                            clusterLabels = EPITHELIA_proliferative$cell_cluster,
                            end.clus =c('HEALING_EPI','SECRETORY_EPI','CILIATED_EPI','MYOGENIC_EPI'),
                            stretch = 4,
                            approx_points = 150)
SlingshotDataSet(sce_slingshot1) 


#定义颜色
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }}
#设置颜色并可视化
cell_colors <- cell_pal(sce_slingshot1$cell_cluster, brewer_pal("qual", "Set2"))
plot(reducedDims(sce_slingshot1)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.8)

#plot(reducedDims(sce_slingshot1)$UMAP, group.by = "sce_slingshot1$cell_cluster2" ,cols = c("#A0522D","#DC143C","#20B2AA","#FFA500","#9370DB","#9400D3","#0000FF"), pch=16, asp = 1, cex = 0.8)

lines(SlingshotDataSet(sce_slingshot1),  lwd=2, col='black')

#计算celltype坐标位置，用于图中标记
celltype_label <- EPITHELIA_proliferative@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(celltype = EPITHELIA_proliferative@meta.data$cell_cluster) %>%
  group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP_1),
            UMAP2 = median(UMAP_2))

for (i in 1:8) {
  text(celltype_label$celltype[i], x=celltype_label$UMAP1[i]-1, y=celltype_label$UMAP2[i])
}

#也可以把slingshot结果放到单细胞seurat对象中-Metadata
pseudotime = slingPseudotime(sce_slingshot1)%>% as.data.frame() 
Lineages = colnames(pseudotime)
EPITHELIA_proliferative <- AddMetaData(object = EPITHELIA_proliferative,metadata = pseudotime)


#EPITHELIA_proliferative-slingshot-tradeseq=============

EPITHELIA_proliferative_slingshot<- CreateSeuratObject(counts = EPITHELIA_proliferative@assays$RNA@counts,meta.data=EPITHELIA_proliferative@meta.data,project = "EPITHELIA_proliferative_slingshot")
EPITHELIA_proliferative_slingshot <- NormalizeData(EPITHELIA_proliferative_slingshot, normalization.method = "LogNormalize", scale.factor = 10000)
EPITHELIA_proliferative_slingshot <- FindVariableFeatures(EPITHELIA_proliferative_slingshot, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(EPITHELIA_proliferative_slingshot)
scale_gene<-VariableFeatures(EPITHELIA_proliferative_slingshot)
#counts<-CE_MD_slingshot@assays$RNA@counts
#counts<counts[scale_gene,]
EPITHELIA_proliferative_slingshot_HVG<-EPITHELIA_proliferative_slingshot[scale_gene,]

#基因随轨迹表达变化
slingsce<-SlingshotDataSet(sce_slingshot1)
pseudotimeED <- as.matrix(unlist(slingPseudotime(slingsce, na = FALSE)))
cellWeightsED <-  as.matrix(unlist(slingCurveWeights(slingsce)))
#counts<-sce_slingshot1@assays@data@listData
#counts<- as.matrix(unlist(sce_slingshot1@assays@data@listData))
#counts<- as.matrix(unlist(MD_2@assays$RNA@counts))
counts<- as.matrix(unlist(EPITHELIA_proliferative_slingshot_HVG@assays$RNA@counts))

#write.table(EPITHELIA_proliferative$cell_cluster,file="EPITHELIA_proliferative_cell_cluster.txt")

#write.table(pseudotime$Lineage1,file="Lineage1.txt")

#这里用HVG基因
system.time({ 
  sce_EPITHELIA_proliferative_slinghot <- fitGAM(counts = counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T)
})


save(sce_EPITHELIA_proliferative_slinghot, file = "sce_EPITHELIA_proliferative_slinghot202405.Rda")


#探索基因表达与拟时序的相关性
assoRes <- associationTest(sce_EPITHELIA_proliferative_slinghot)
head(assoRes)
write.table(assoRes,file="EPITHELIA_proliferative_Pseudotime_association_GENE202405.txt")


#寻找与起止点相关性最高的基因
startRes <- startVsEndTest(sce_EPITHELIA_proliferative_slinghot)
head(startRes)
write.table(startRes,file="EPITHELIA_proliferative_Pseudotime_startRes_GENE202405.txt")

# 按相关性进行排序
oStart <- order(startRes$waldStat, decreasing = TRUE)

# 挑选相关性最强的基因，并可视化
sigGeneStart <- names(sce_EPITHELIA_proliferative_slinghot)[oStart[1]]
plotSmoothers(sce_EPITHELIA_proliferative_slinghot, counts, gene = "TPPP3")


#也可以用UMAP图展示
plotGeneCount(slingsce, counts, gene = "TPPP3")
plotGeneCount(slingsce, counts, gene = c("UPK3B","MKI67","HMGB2","LHX2","COL1A1","COL3A1","COL1A2","ADAMTS9","RSPO3","HOXA10","BMP7","MMP7","NES"))
#=========================================================================
#https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
library("RColorBrewer")
library("pheatmap")

assogene<-read.table('EPITHELIA_proliferative_Pseudotime_association_GENE2.txt',header = T,row.names = 1)
topgenes <- rownames(assogene)[1:40]
pst.ord <- order(sce_EPITHELIA_proliferative_slinghot$crv, na.last = NA)
#heatdata <- CE_MD@assays$counts[topgenes, pst.ord]
heatdata<- as.matrix(unlist(EPITHELIA_proliferative@assays$RNA@counts))[topgenes, pst.ord]



heatmap(log1p(heatdata), Colv = NA,col = colorRampPalette(rev(brewer.pal(3, "RdGy")))(250))


heatclus <- unlist(CE_MD$cell_cluster)
heatmap(log1p(heatdata), Colv = NA,col = brewer.pal(3,"Reds"),
        ColSideColors = brewer.pal(3,"Reds")[heatclus])


mm=log1p(heatdata)
max_range = max(range(is.finite(mm)))
lim = c(-max_range, max_range)
library(pheatmap)
heatmap1 = pheatmap(mm, show_rownames= T, cluster_rows = TRUE,
                    cluster_cols = FALSE, show_colnames = FALSE, 
                    clustering_distance_rows = "euclidean",
                    clustering_method = "ward.D2",
                    treeheight_row = 10,
                    cutree_rows = 10, 
                    color = colorRampPalette(rev(brewer.pal(9, "RdGy")))(250),
                    breaks = seq(lim[1]/4, lim[2]/4, length.out = 251),
                    border_color = NA)

heatmap1 = pheatmap(mm, show_rownames= T, cluster_rows = TRUE,
                    cluster_cols = FALSE, show_colnames = FALSE, 
                    clustering_distance_rows = "correlation",
                    clustering_method = "manhattan",
                    treeheight_row = 20,
                    
                    color = colorRampPalette(rev(brewer.pal(9, "RdGy")))(250),
                    breaks = seq(lim[1]/4, lim[2]/4, length.out = 251),
                    border_color = NA)


heatmap1 = pheatmap(mm, show_rownames= T, cluster_rows = TRUE,
                    cluster_cols = FALSE, show_colnames = FALSE, 
                    
                    treeheight_row = 20,
                    cutree_rows = 10, 
                    color = colorRampPalette(rev(brewer.pal(9, "RdGy")))(250),
                    breaks = seq(lim[1]/4, lim[2]/4, length.out = 251),
                    border_color = NA)

#=========================================================================
##endothelia########################

save(ENDOTHELIA_proliferative, file = "ENDOTHELIA_proliferative.Rda")


cellinfo3<-read.table('ENDOTHELIA_proliferative-cluster.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
ENDOTHELIA_proliferative<- AddMetaData(object = ENDOTHELIA_proliferative, metadata = cellinfo3)

ElbowPlot(ENDOTHELIA_proliferative)
ENDOTHELIA_proliferative <- RunUMAP(ENDOTHELIA_proliferative, reduction = "pca", dims = 1:11)
ENDOTHELIA_proliferative <- FindNeighbors(ENDOTHELIA_proliferative, reduction = "pca", dims = 1:11)
ENDOTHELIA_proliferative <- FindClusters(ENDOTHELIA_proliferative, resolution = 0.2)
# Visualization-UMAP
p1 <- DimPlot(ENDOTHELIA_proliferative, reduction = "umap", group.by = "stim" ,cols = c("#DC143C","#1E90FF"),pt.size = 0.5)
p2 <- DimPlot(ENDOTHELIA_proliferative, reduction = "umap", label = T, pt.size = 1)
plot_grid(p1)
p3 <- DimPlot(ENDOTHELIA_proliferative, reduction = "umap", group.by = "original_ID" ,cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                                                                                               "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                                                                                               "#808000","#FF00FF","#40E0D0","#7B68EE",
                                                                                               "#9400D3","#800080","#A0522D","#D2B48C"),pt.size = 0.5)
plot_grid(p3)
p4 <- DimPlot(ENDOTHELIA_proliferative, reduction = "umap", group.by = "cell_cluster",cols = c("#A0522D","#DC143C","#20B2AA","#FFA500",
                                                                                                "#9370DB","#9400D3","#0000FF"),label = T,pt.size = 0.5)
plot_grid(p4)
plot_grid(p2)

###SCENIC########################################################################
setwd("E:/IUA-revision/scenic_endo")
library(Seurat) 
library(SCENIC)

Idents(object = ENDOTHELIA_proliferative) <- 'cell_cluster'
# Downsample the number of cells per identity class
ENDOTHELIA_proliferative1 <- subset(x = ENDOTHELIA_proliferative, downsample = 100)

## Load data表达矩阵加上表型信息
exprMat  <-  as.matrix(ENDOTHELIA_proliferative1@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  ENDOTHELIA_proliferative1@meta.data
save(exprMat,cellInfo,file = "ENDOTHELIA_proliferative1_scenic.Rdata")

### Initialize settings
org<-"hgnc"
dbDir<-"cisTarget_databases"

myDatasetTitle<-"SCENIC analysis on ENDOTHELIA_proliferative"
data("defaultDbNames")
dbs<-defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

motifAnnotations_hgnc = motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
#创建一个对象：scenicOptions，后续数据会保存在当前目录的int文件夹
#输出结果会存储在output文件夹
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#将表达矩阵中未在cisTarget_databases收录的基因去除
library(RcisTarget)
dbFilePath<-getDatabases(scenicOptions)[[1]]
motifRanking<-importRankings(dbFilePath)
geneInDatabase<-colnames(getRanking(motifRanking))
genelefet_minCells<-rownames(exprMat)
length(genelefet_minCells)
genelefet_minCells_inDatabases<-genelefet_minCells[which(genelefet_minCells %in% geneInDatabase)]
length(genelefet_minCells_inDatabases)
genesKept<-genelefet_minCells_inDatabases
exprMat_filter<-exprMat[genesKept,]
dim(exprMat)
dim(exprMat_filter)

#计算相关性并对数据进行log处理
### Co-expression network
class(exprMat_filter)
runCorrelation(as.matrix(exprMat_filter), scenicOptions)
exprMat_filter <- log2(exprMat_filter+1) 
library(GENIE3)
runGenie3(as.matrix(exprMat_filter), scenicOptions)
save(exprMat_filter,scenicOptions,file="input_GENIE3_data.Rdata")


### Build the GRN
#load(stroma_scenic.Rdata)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
#closeAllConnections()
#scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget")) # Toy run settings
#library(doParallel)
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#结果可视化
library(AUCell)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)


#Regulators for known cell types or clusters

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cell_cluster3),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, col = col,name="Regulon activity")


###########################################################################
#ENDOTHELIA_proliferative-slingshot=============

library(RColorBrewer)
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(mclust)
library(tradeSeq)

library(slingshot)
library(Seurat)
library(devtools)
library(cowplot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tradeSeq)
library(RColorBrewer)
library(DelayedMatrixStats)
library(scales)
library(paletteer) 
library(viridis)
library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(reshape2)
#library(SeuratData)
library(patchwork)

MD_2_0<- as.SingleCellExperiment(ENDOTHELIA_proliferative)
#save(MD_3, file = "MD_3_slingshot.Rda")
sce_slingshot1 <- slingshot(MD_2_0,      #输入单细胞对象
                            reducedDim = 'UMAP',  #降维方式
                            start.clus = 'PROLIFERATIVE_ENDO',       #轨迹起点,也可以不定义
                            clusterLabels = ENDOTHELIA_proliferative$cell_cluster,
                            end.clus =c('INFLAMMATORY_ENDO','JAG1_ENDO','MYOGENIC_ENDO'),
                            stretch = 3,
                            approx_points = 150)
SlingshotDataSet(sce_slingshot1) 


#定义颜色
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }}
#设置颜色并可视化
cell_colors <- cell_pal(sce_slingshot1$cell_cluster, brewer_pal("qual", "Set2"))
plot(reducedDims(sce_slingshot1)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.8)

#plot(reducedDims(sce_slingshot1)$UMAP, group.by = "sce_slingshot1$cell_cluster2" ,cols = c("#A0522D","#DC143C","#20B2AA","#FFA500","#9370DB","#9400D3","#0000FF"), pch=16, asp = 1, cex = 0.8)

lines(SlingshotDataSet(sce_slingshot1),  lwd=2, col='black')

#计算celltype坐标位置，用于图中标记
celltype_label <- ENDOTHELIA_proliferative@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(celltype = ENDOTHELIA_proliferative@meta.data$cell_cluster) %>%
  group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP_1),
            UMAP2 = median(UMAP_2))

for (i in 1:8) {
  text(celltype_label$celltype[i], x=celltype_label$UMAP1[i]-1, y=celltype_label$UMAP2[i])
}

#也可以把slingshot结果放到单细胞seurat对象中-Metadata
pseudotime = slingPseudotime(sce_slingshot1)%>% as.data.frame() 
Lineages = colnames(pseudotime)
ENDOTHELIA_proliferative <- AddMetaData(object = ENDOTHELIA_proliferative,metadata = pseudotime)


#ENDOTHELIA_proliferative-slingshot-tradeseq=============

ENDOTHELIA_proliferative_slingshot<- CreateSeuratObject(counts = ENDOTHELIA_proliferative@assays$RNA@counts,meta.data=ENDOTHELIA_proliferative@meta.data,project = "ENDOTHELIA_proliferative_slingshot")
ENDOTHELIA_proliferative_slingshot <- NormalizeData(ENDOTHELIA_proliferative_slingshot, normalization.method = "LogNormalize", scale.factor = 10000)
ENDOTHELIA_proliferative_slingshot <- FindVariableFeatures(ENDOTHELIA_proliferative_slingshot, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(ENDOTHELIA_proliferative_slingshot)
scale_gene<-VariableFeatures(ENDOTHELIA_proliferative_slingshot)
#counts<-CE_MD_slingshot@assays$RNA@counts
#counts<counts[scale_gene,]
ENDOTHELIA_proliferative_slingshot_HVG<-ENDOTHELIA_proliferative_slingshot[scale_gene,]

#基因随轨迹表达变化
slingsce<-SlingshotDataSet(sce_slingshot1)
pseudotimeED <- as.matrix(unlist(slingPseudotime(slingsce, na = FALSE)))
cellWeightsED <-  as.matrix(unlist(slingCurveWeights(slingsce)))
#counts<-sce_slingshot1@assays@data@listData
#counts<- as.matrix(unlist(sce_slingshot1@assays@data@listData))
#counts<- as.matrix(unlist(MD_2@assays$RNA@counts))
counts<- as.matrix(unlist(ENDOTHELIA_proliferative_slingshot_HVG@assays$RNA@counts))

#write.table(EPITHELIA_proliferative$cell_cluster,file="EPITHELIA_proliferative_cell_cluster.txt")

#write.table(pseudotime$Lineage1,file="Lineage1.txt")

#这里用HVG基因
system.time({ 
  sce_ENDOTHELIA_proliferative_slinghot <- fitGAM(counts = counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T)
})


save(sce_ENDOTHELIA_proliferative_slinghot, file = "sce_ENDOTHELIA_proliferative_slinghot202405.Rda")


#探索基因表达与拟时序的相关性
assoRes <- associationTest(sce_ENDOTHELIA_proliferative_slinghot)
head(assoRes)
write.table(assoRes,file="ENDOTHELIA_proliferative_Pseudotime_association_GENE202405.txt")


#寻找与起止点相关性最高的基因
startRes <- startVsEndTest(sce_ENDOTHELIA_proliferative_slinghot)
head(startRes)
write.table(startRes,file="ENDOTHELIA_proliferative_Pseudotime_startRes_GENE202405.txt")

# 按相关性进行排序
oStart <- order(startRes$waldStat, decreasing = TRUE)

# 挑选相关性最强的基因，并可视化
sigGeneStart <- names(sce_ENDOTHELIA_proliferative_slinghot)[oStart[1]]
plotSmoothers(sce_ENDOTHELIA_proliferative_slinghot, counts, gene = "TPPP3")


#也可以用UMAP图展示
plotGeneCount(slingsce, counts, gene = "TPPP3")
plotGeneCount(slingsce, counts, gene = c("UPK3B","MKI67","HMGB2","LHX2","COL1A1","COL3A1","COL1A2","ADAMTS9","RSPO3","HOXA10","BMP7","MMP7","NES"))
#=========================================================================
#immune=================================================================
save(IMMUNE.combined, file = "IMMUNE.combined.Rda")

#==================================================================
cellinfo3<-read.table('IMMUNE.combined_cell_cluster.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
IMMUNE.combined<- AddMetaData(object = IMMUNE.combined, metadata = cellinfo3)

ElbowPlot(IMMUNE.combined)

IMMUNE.combined <- RunUMAP(IMMUNE.combined, reduction = "pca", dims = 1:14)
IMMUNE.combined <- FindNeighbors(IMMUNE.combined, reduction = "pca", dims = 1:14)
IMMUNE.combined <- FindClusters(IMMUNE.combined, resolution = 0.5)
# Visualization-UMAP
p1 <- DimPlot(IMMUNE.combined, reduction = "umap", group.by = "stim" ,
              cols = c("#DC143C","#1E90FF"),
              pt.size = 0.5)

p4 <- DimPlot(IMMUNE.combined, reduction = "umap", group.by = "original_ID" ,
              cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                       "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                       "#808000","#FF00FF","#40E0D0","#7B68EE",
                       "#9400D3","#800080","#A0522D","#D2B48C"),pt.size = 0.5)

p6 <- DimPlot(IMMUNE.combined, reduction = "umap", group.by = "cell_cluster" , cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                                                                                         "#9370DB","#1E90FF","#FF00FF","#808000"),label = T,pt.size = 0.5)

plot_grid(p1)
plot_grid(p4)
plot_grid(p6)


##########################SCATAC-signac

library(Seurat) 
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)

#Pre-processing workflow
setwd("E:/IUA-ATAC/matrix")
counts_N <- Read10X('E458filtered_peak_bc_matrix',gene.column = 1)

metadata_N <- read.csv(
  file = "E_458_singlecell.csv",
  header = TRUE,
  row.names = 1)

chrom_assay_N <- CreateChromatinAssay(
  counts = counts_N,
  sep = c(":", "-"),
  fragments = 'E458fragments.tsv.gz',
  min.cells = 10,
  min.features = 200)

NORMAL_ATAC <- CreateSeuratObject(
  counts = chrom_assay_N,
  assay = "peaks",
  meta.data = metadata_N)

granges(NORMAL_ATAC)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(NORMAL_ATAC) <- annotations

# Computing QC Metrics-NORMAL_ATAC
# compute nucleosome signal score per cell
NORMAL_ATAC <- NucleosomeSignal(object = NORMAL_ATAC)
# compute TSS enrichment score per cell
NORMAL_ATAC <- TSSEnrichment(object = NORMAL_ATAC, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
NORMAL_ATAC$pct_reads_in_peaks <- NORMAL_ATAC$peak_region_fragments / NORMAL_ATAC$passed_filters * 100
NORMAL_ATAC$blacklist_ratio <- NORMAL_ATAC$blacklist_region_fragments / NORMAL_ATAC$peak_region_fragments


#DensityScatter(NORMAL_ATAC, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
NORMAL_ATAC$high.tss <- ifelse(NORMAL_ATAC$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(NORMAL_ATAC, group.by = 'high.tss') + NoLegend()
NORMAL_ATAC$nucleosome_group <- ifelse(NORMAL_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = NORMAL_ATAC, group.by = 'nucleosome_group')

VlnPlot(
  object = NORMAL_ATAC,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5)

NORMAL_ATAC <- subset(
  x = NORMAL_ATAC,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)

#Normalization and linear dimensional reduction
NORMAL_ATAC <- RunTFIDF(NORMAL_ATAC)
NORMAL_ATAC <- FindTopFeatures(NORMAL_ATAC, min.cutoff = 'q0')
NORMAL_ATAC <- RunSVD(NORMAL_ATAC)
saveRDS(NORMAL_ATAC,file="NORMAL_ATAC.rds")
#################################################################################

counts_IUA <- Read10X('E123679filtered_peak_bc_matrix',gene.column = 1)
metadata_IUA <- read.csv(
  file = "E_123679_singlecell.csv",
  header = TRUE,
  row.names = 1)
chrom_assay_IUA <- CreateChromatinAssay(
  counts = counts_IUA,
  sep = c(":", "-"),
  fragments = 'E123679fragments.tsv.gz',
  min.cells = 10,
  min.features = 200)

IUA_ATAC <- CreateSeuratObject(
  counts = chrom_assay_IUA,
  assay = "peaks",
  meta.data = metadata_IUA)

granges(IUA_ATAC)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(IUA_ATAC) <- annotations

# Computing QC Metrics-IUA_ATAC
# compute nucleosome signal score per cell
IUA_ATAC <- NucleosomeSignal(object = IUA_ATAC)
# compute TSS enrichment score per cell
IUA_ATAC <- TSSEnrichment(object = IUA_ATAC, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
IUA_ATAC$pct_reads_in_peaks <- IUA_ATAC$peak_region_fragments / IUA_ATAC$passed_filters * 100
IUA_ATAC$blacklist_ratio <- IUA_ATAC$blacklist_region_fragments / IUA_ATAC$peak_region_fragments

# DensityScatter(IUA_ATAC, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
IUA_ATAC$high.tss <- ifelse(IUA_ATAC$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(IUA_ATAC, group.by = 'high.tss') + NoLegend()
IUA_ATAC$nucleosome_group <- ifelse(IUA_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = IUA_ATAC, group.by = 'nucleosome_group')

VlnPlot(
  object = IUA_ATAC,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5)

IUA_ATAC <- subset(
  x = IUA_ATAC,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)

#Normalization and linear dimensional reduction
IUA_ATAC <- RunTFIDF(IUA_ATAC)
IUA_ATAC <- FindTopFeatures(IUA_ATAC, min.cutoff = 'q0')
IUA_ATAC <- RunSVD(IUA_ATAC)

saveRDS(IUA_ATAC,file="IUA_ATAC.rds")

#Non-linear dimension reduction and clustering
IUA_ATAC <- RunUMAP(object = IUA_ATAC, reduction = 'lsi', dims = 2:30)
IUA_ATAC <- FindNeighbors(object = IUA_ATAC, reduction = 'lsi', dims = 2:30)
IUA_ATAC <- FindClusters(object = IUA_ATAC, verbose = FALSE, algorithm = 3)
#DimPlot(object = IUA_ATAC, label = TRUE) + NoLegend()

#################################################################################
#################################################################################
gene.activities=GeneActivity(NORMAL_ATAC)
NORMAL_ATAC[["ACTIVITY"]]= CreateAssayObject(counts=gene.activities)
NORMAL_ATAC$dataset <- "Normal"

#NORMAL_ATAC <- NormalizeData(
#  object = NORMAL_ATAC,
#  assay = 'RNA',
#  normalization.method = 'LogNormalize',
#  scale.factor = median(NORMAL_ATAC$nCount_RNA))
#DefaultAssay(NORMAL_ATAC) <- 'RNA'

#FeaturePlot(
#  object = NORMAL_ATAC,
#  features = c('KRT8', 'SFRP4', 'COL1A1', 'IGF1', 'PTPRC', 'VWF'),
#  pt.size = 0.1,
#  max.cutoff = 'q95',
#  ncol = 3)

gene.activities2=GeneActivity(IUA_ATAC)
IUA_ATAC[["ACTIVITY"]]= CreateAssayObject(counts=gene.activities2)
IUA_ATAC$dataset <- "IUA"

IUA_ATAC <- NormalizeData(
  object = IUA_ATAC,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(IUA_ATAC$nCount_RNA))
# DefaultAssay(IUA_ATAC) <- 'RNA'

FeaturePlot(
  object = IUA_ATAC,
  features = c('KRT8', 'SFRP4', 'COL1A1', 'IGF1', 'PTPRC', 'VWF'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3)


# merge
atac.merged <- merge(NORMAL_ATAC, IUA_ATAC)

# process the combined dataset
# compute LSI
atac.merged <- FindTopFeatures(atac.merged, min.cutoff = 10)
atac.merged <- RunTFIDF(atac.merged)
atac.merged <- RunSVD(atac.merged)

atac.merged <- RunUMAP(atac.merged, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(atac.merged, group.by = "dataset")
p1 + ggtitle("ATAC_Merged")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(NORMAL_ATAC, IUA_ATAC),
  anchor.features = rownames(NORMAL_ATAC),
  dims = 2:30)

# integrate LSI embeddings
ATAC_integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = atac.merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30)

gene.activities2=GeneActivity(ATAC_integrated)
ATAC_integrated[["ACTIVITY"]]= CreateAssayObject(counts=gene.activities2)

# create a new UMAP using the integrated embeddings
ATAC_integrated <- RunUMAP(ATAC_integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(ATAC_integrated, group.by = "dataset")
p2 + ggtitle("ATAC_integrated")

#(p1 + ggtitle("ATAC_Merged")) | (p2 + ggtitle("ATAC_integrated"))

save(ATAC_integrated,file = "ATAC_integrated_IUA2.RData")

ElbowPlot(ATAC_integrated)

ATAC_integrated <- RunUMAP(object = ATAC_integrated, reduction = 'integrated_lsi', dims = 2:30)
ATAC_integrated <- FindNeighbors(object = ATAC_integrated, reduction = 'integrated_lsi', dims = 2:30)
ATAC_integrated <- FindClusters(object = ATAC_integrated, verbose = FALSE, algorithm = 3)
DimPlot(object = ATAC_integrated, label = TRUE) + NoLegend()

FeaturePlot(
  object = ATAC_integrated,
  features = c('KRT8','EPCAM', 'SFRP4', 'COL1A1', 'IGF1', 'VWF','ACTA2','CNN1','PTPRC','CD68','CD3E','NKG7'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3)

# compute LSI
ATAC_integrated <- FindTopFeatures(ATAC_integrated, min.cutoff = 10)
ATAC_integrated <- RunTFIDF(ATAC_integrated)
ATAC_integrated <- RunSVD(ATAC_integrated)

# create a new UMAP using the integrated embeddings
ATAC_integrated <- RunUMAP(ATAC_integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(ATAC_integrated, group.by = "dataset")
p2 + ggtitle("ATAC_integrated")


#Integrating with scRNA-seq data
#################################################################
transfer.anchors <- FindTransferAnchors(
  reference = All_IUA_proliferative,
  query = ATAC_integrated,
  features = VariableFeatures(object = All_IUA_proliferative), 
  reference.assay = "RNA", 
  query.assay = "ACTIVITY", 
  reduction = 'cca')

save(transfer.anchors,file = "transfer.anchors202406.RData")

#cell_types
#cell_cluster

# All_IUA_c.combined1102 (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  predicted.labels (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells 
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = All_IUA_proliferative$cell_cluster6,
  weight.reduction = ATAC_integrated[['integrated_lsi']],
  dims = 2:30)

save(predicted.labels,file = "predicted.labels_cell_types.Rdata")

save(predicted.labels,file = "predicted.labels_cell_cluster.Rdata")

Idents(object = All_IUA_proliferative) <- 'cell_types'

Idents(object = All_IUA_proliferative) <- 'cell_cluster'

ATAC_RNA_integrated <- AddMetaData(object = ATAC_integrated, metadata = predicted.labels)

save(ATAC_RNA_integrated,file = "ATAC_RNA_integrated_cell_types.Rdata")

save(ATAC_RNA_integrated,file = "ATAC_RNA_integrated_cell_cluster.Rdata")

plot1 <- DimPlot(
  object = All_IUA_proliferative,
  group.by = 'cell_types',
  label = TRUE,
  cols = c("#20B2AA","#0000FF","#FFA500","#9370DB","#DC143C"),
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = All_IUA_proliferative,
  group.by = 'cell_cluster',
  label = TRUE,
  cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
           "#9370DB","#98FB98","#F08080","#1E90FF",
           "#7CFC00","#FFFF00","#808000","#FF00FF",
           "#FA8072","#7B68EE","#9400D3","#800080",
           "#A0522D","#D2B48C","#D2691E","#87CEEB",
           "#40E0D0","#5F9EA0","#FF1493","#0000CD",
           "#FFE4B5","#8A2BE2","#228B22"),
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')


plot1 
plot2

plot1 <- DimPlot(
  object = ATAC_RNA_integrated,
  group.by = 'predicted.id',
  reduction = "umap",
  label = TRUE,
  raster =F,
  cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
           "#9370DB","#98FB98"          ,"#1E90FF",
           "#7CFC00","#FFFF00",          "#FF00FF",
           "#FA8072","#7B68EE",          "#800080",
           "#A0522D","#D2B48C","#D2691E","#87CEEB",
           "#40E0D0","#5F9EA0","#FF1493","#0000CD",
           "#008B8B","#FFE4B5","#8A2BE2","#228B22",
           "#E9967A"),pt.size = 1,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 

p4 <- DimPlot(All_IUA_proliferative, reduction = "umap", group.by = "dataset" ,cols = c("#20B2AA","#FFA500","#9370DB","#DC143C"),pt.size = 0.1,raster =F)
p6 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE,group.by = "cell_cluster" ,
              cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                       "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                       "#808000","#FF00FF","#FA8072","#7B68EE",
                       "#9400D3","#800080","#A0522D","#D2B48C",
                       "#D2691E","#87CEEB","#40E0D0","#5F9EA0",
                       "#FF1493","#0000CD","#008B8B","#FFE4B5",
                       "#8A2BE2","#228B22"),
              pt.size = 0.1,raster =F)
p7 <- DimPlot(All_IUA_proliferative, reduction = "umap", group.by = "original_ID" ,
              cols = c("#DC143C","#0000FF","#20B2AA","#FFA500",
                       "#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                       "#808000","#FF00FF","#FA8072","#7B68EE",
                       "#9400D3","#800080","#A0522D","#D2B48C"),
              pt.size = 0.1,raster =F)
plot_grid(p7)
plot_grid(p4)
plot_grid(p6)

plot1 <- DimPlot(
  object = All_IUA_proliferative,
  group.by = 'stim',
  label = TRUE,
  cols = c("#DC143C","#1E90FF"),
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = ATAC_RNA_integrated,
  group.by = 'dataset',
  label = TRUE,
  cols = c("#DC143C","#1E90FF"),
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 
plot2

save(ATAC_RNA_integrated,file = "ATAC_RNA_integrated2.RData")

#heatmap-proliferative
Idents(object = All_IUA_proliferative) <- 'cell_cluster'
All_IUA_proliferative.markers <- FindAllMarkers(All_IUA_proliferative, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
All_IUA_proliferative.markers %>%group_by(cluster) %>%top_n(n = 5, wt = avg_log2FC) -> top10

My_levels <- c('STROMA','EPITHELIA','ENDOTHELIA','IMMUNE','MUSCLE')
Idents(All_IUA_proliferative) <- factor(Idents(All_IUA_proliferative), levels= My_levels)

DoHeatmap(All_IUA_proliferative, features = top10$gene,label = F) + NoLegend()


DoHeatmap(All_IUA_proliferative,features = top10$gene,label = F,
          group.colors=c("#D2B48C","#008B8B","#8A2BE2","#FA8072",
                         "#0000CD","#7B68EE","#A0522D","#1E90FF",
                         "#FFFF00","#7CFC00","#D2691E","#20B2AA",
                         "#E9967A","#9370DB","#40E0D0","#FF1493",
                         "#FF00FF","#5F9EA0","#FFA500","#98FB98",
                         "#FFE4B5","#0000FF","#DC143C","#800080",
                         "#228B22","#87CEEB"))+scale_fill_gradientn (colors = c ("white","grey","firebrick3"))



markers.to.plot <- c("PTPRC","CD68","CD3E","COL1A1","ECM1","KRT8","KRT18","EPCAM","VWF","CLDN5","ACTA2")
DotPlot(All_IUA_proliferative, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_type2")+
  RotatedAxis()

markers.to.plot <- c("VWF","CLDN5","THY1","SLCO2A1","PTPRC","CD3E","COL1A1","ECM1","SFRP4","MKI67","KRT8","KRT18","EPCAM","GNLY","ACTA2","CD68","IL1B","JAG1","NFKBIA","CD1C","CLU","TPPP3","JCHAIN")
DotPlot(All_IUA_proliferative, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_cluster5")+
  RotatedAxis()

VlnPlot(All_IUA_proliferative, features = c("PTPRC","CD68","CD3E","COL1A1","ECM1","VIM","KRT8","KRT18","EPCAM","VWF","CLDN5","ACTA2"),pt.size = 0,ncol = 3)
FeaturePlot(All_IUA_proliferative, features = c("PTPRC","CD68","CD3E","COL1A1","ECM1","VIM","KRT8","KRT18","EPCAM","VWF","CLDN5","ACTA2"), min.cutoff = "q9")

#heatmap-ATAC-RNA-label
#https://blog.csdn.net/ZIGRA/article/details/132362313
Idents(object = ATAC_RNA_integrated) <- 'predicted.id'

ATAC_RNA_integrated.markers <- FindAllMarkers(ATAC_RNA_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ATAC_RNA_integrated.markers %>%group_by(cluster) %>%top_n(n = 5, wt = avg_log2FC) -> top10

DoHeatmap(ATAC_RNA_integrated, 
          features = top10$gene,slot = "data", label = F,
          group.colors =c("#FFA500","#DC143C","#0000FF","#9370DB","#20B2AA"))+
  scale_fill_gradientn (colors = c ("white","firebrick3"))


write.table(ATAC_RNA_integrated.markers,file="ATAC_RNA_integrated.markers_cell_types.txt")

write.table(ATAC_RNA_integrated.markers,file="ATAC_RNA_integrated.markers_cell_cluster.txt")

Idents(ATAC_RNA_integrated)="predicted.id"
EPITHELIA_ATAC<- subset(ATAC_RNA_integrated,idents = "EPITHELIA")
save(EPITHELIA_ATAC, file = "EPITHELIA_ATAC_202407.Rda")
STROMA_ATAC<- subset(ATAC_RNA_integrated,idents = "STROMA")
save(STROMA_ATAC, file = "STROMA_ATAC_202407.Rda")
ENDOTHELIA_ATAC<- subset(ATAC_RNA_integrated,idents = "ENDOTHELIA")
save(ENDOTHELIA_ATAC, file = "ENDOTHELIA_ATAC_202407.Rda")

#GeneActivity##################################################################
gene.activities <- GeneActivity(ATAC_RNA_integrated)
ATAC_RNA_integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
ATAC_RNA_integrated <- NormalizeData(
  object = ATAC_RNA_integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC_RNA_integrated$nCount_RNA))

DefaultAssay(ATAC_RNA_integrated) <- 'RNA'

FeaturePlot(
  object = ATAC_RNA_integrated,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3)

ATAC_RNA_integrated.markers <- FindAllMarkers(ATAC_RNA_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ATAC_RNA_integrated.markers %>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10

write.table(ATAC_RNA_integrated.markers,file="ATAC_RNA_integrated.ACTIVITY_markers_cell_type5.txt")

DoHeatmap(ATAC_RNA_integrated, 
          features = top10$gene,slot = "data", label = F,
          group.colors =c("#FFA500","#DC143C","#0000FF","#9370DB","#20B2AA"))+
  scale_fill_gradientn (colors = c ("grey","firebrick3"))


DoHeatmap(ATAC_RNA_integrated, 
          features = top10$gene,slot = "data", label = F, group.colors=c("#A0522D","#FF1493","#FFE4B5","#FF00FF",
                                                                         "#5F9EA0","#FA8072","#800080",         "#1E90FF",
                                                                         "#FFFF00","#7CFC00",         "#D2B48C",
                                                                         "#20B2AA"         ,"#228B22","#9370DB",
                                                                         "#87CEEB","#40E0D0","#FFA500","#98FB98",
                                                                         "#0000CD","#0000FF","#DC143C","#9400D3",
                                                                         "#D2691E"          ))+scale_fill_gradientn (colors = c ("grey","firebrick3"))

#Integrating with scRNA-seq data
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca')

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

# replace each label with its most likely prediction
for(i in levels(pbmc)) {
  cells_to_reid <- WhichCells(pbmc, idents = i)
  newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
  Idents(pbmc, cells = cells_to_reid) <- newid
}
markers.to.plot <- c("PTPRC","CD68","CD3E","COL1A1","ECM1","VIM","KRT8","KRT18","EPCAM","VWF","CLDN5","ACTA2")
DotPlot(All_IUA_proliferative, features = markers.to.plot, cols =  c("blue", "red"), dot.scale = 8, split.by = "stim", group.by ="cluster")+
  RotatedAxis()

markers.to.plot <- c("PTPRC","CD68","CD3E","COL1A1","ECM1","VIM","KRT8","KRT18","EPCAM","VWF","CLDN5","ACTA2")
DotPlot(ATAC_RNA_integrated, features = markers.to.plot, cols =  c("purple", "yellow"), dot.scale = 8, split.by = "stim", group.by ="cluster")+
  RotatedAxis()

CoveragePlot(
  object = ATAC_RNA_integrated,
  region = c("PTPRC","CD68","CD3E","COL1A1","ECM1","VIM","KRT8","KRT18","EPCAM","VWF","CLDN5","ACTA2"),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 3)


cols = c("#20B2AA","#0000FF","#FFA500","#9370DB","#DC143C")
p <- CoveragePlot(
  object = ATAC_RNA_integrated,group.by = 'predicted.id',
  region = c("PTPRC","CD68","CD3E"),
  extend.upstream = 1000,
  extend.downstream = 1000,ncol = 3)

p & scale_fill_manual(values =cols)
###############################################


#cellchat##########################################################################
All_IUA_proliferative <- SetIdent(All_IUA_proliferative,cells = names(All_IUA_proliferative@meta.data$cell_cluster5),value = All_IUA_proliferative@meta.data$cell_cluster5)
All_IUA_proliferative_cellchat <- subset(All_IUA_proliferative,idents = c("B_CELLS_1","B_CELLS_2","DC","NK_1","NK_2","M1_MACROPHAGE","M2_MACROPHAGE","PROLIFERATIVE_M","T_CELLS_1","T_CELLS_2","MYOGENIC_ENDO","MYOGENIC_EPI","MYOGENIC_STROMA"))

All_IUA_proliferative_cellchat <- SetIdent(All_IUA_proliferative,cells = names(All_IUA_proliferative@meta.data$stim),value = All_IUA_proliferative@meta.data$stim)
All_IUA_proliferative_cellchat_N <- subset(All_IUA_proliferative_cellchat,idents = c("Normal"))

All_IUA_proliferative_cellchat <- SetIdent(All_IUA_proliferative,cells = names(All_IUA_proliferative@meta.data$stim),value = All_IUA_proliferative@meta.data$stim)
All_IUA_proliferative_cellchat_IUA <- subset(All_IUA_proliferative_cellchat,idents = c("IUA"))

#normal-myogenic-immune-All_IUA_proliferative_cellchat_N================================================================
All_IUA_proliferative_cellchat_N <- RunUMAP(All_IUA_proliferative_cellchat_N, reduction = "pca", dims = 1:8)
All_IUA_proliferative_cellchat_N <- FindNeighbors(All_IUA_proliferative_cellchat_N, reduction = "pca", dims = 1:8)
All_IUA_proliferative_cellchat_N <- FindClusters(All_IUA_proliferative_cellchat_N, resolution = 0.3)
# Visualization-UMAP
p1 <- DimPlot(All_IUA_proliferative_cellchat_N, reduction = "umap", group.by = "stim" ,pt.size = 0.5,raster =F)
#p2 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE, pt.size = 0.1,raster =F)
#p3 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE,group.by = "cell_type" ,pt.size = 0.1,raster =F)
p4 <- DimPlot(All_IUA_proliferative_cellchat_N, reduction = "umap", group.by = "datasets" ,pt.size = 0.5,raster =F)
p5 <- DimPlot(All_IUA_proliferative_cellchat_N, reduction = "umap", group.by = "cell_types5" ,pt.size = 0.5,raster =F)
#p6 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE,group.by = "cell_cluster" ,pt.size = 0.1,raster =F)
#original_ID
p7 <- DimPlot(All_IUA_proliferative_cellchat_N, reduction = "umap", group.by = "original_ID" ,pt.size = 0.5,raster =F)
p8 <- DimPlot(All_IUA_proliferative_cellchat_N, reduction = "umap", label = F,group.by = "cell_cluster6" ,pt.size = 0.5,raster =F)

plot_grid(p7)
plot_grid(p1)
#plot_grid(p2)
#plot_grid(p3)
plot_grid(p4)
plot_grid(p5)
#plot_grid(p6)
plot_grid(p8)

library(CellChat)
#Extract the CellChat input files from a Seurat object

#建立cellchat对象
#https://www.jianshu.com/p/b3d26ac51c5a
#cellchat <- createCellChat(object = regeneration,group.by="cell_cluster")
cellchat_all_N <- createCellChat(All_IUA_proliferative_cellchat_N@assays$RNA@data,meta=All_IUA_proliferative_cellchat_N@meta.data,group.by="cell_cluster6")
groupSize<-as.numeric(table(cellchat_all_N@idents))
groupSize

#导入配体受体数据库
CellChatDB<-CellChatDB.human
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# set the used database in the object
cellchat_all_N@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_all_N <- subsetData(cellchat_all_N) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat_all_N <- identifyOverExpressedGenes(cellchat_all_N)
cellchat_all_N <- identifyOverExpressedInteractions(cellchat_all_N)

# Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
cellchat_all_N <- computeCommunProb(cellchat_all_N)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_all_N <- filterCommunication(cellchat_all_N, min.cells = 10)
df.net<-subsetCommunication(cellchat_all_N)
write.csv(df.net,"net_LR-All_IUA_proliferative_N3.csv")

#Infer the cell-cell communication at a signaling pathway level
cellchat_all_N <- computeCommunProbPathway(cellchat_all_N)
df.netp<-subsetCommunication(cellchat_all_N,slot.name = "netP")
write.csv(df.netp,"net_pathway_All_IUA_proliferative_N3.csv")

#Calculate the aggregated cell-cell communication network
#所有细胞总体互作数量和强度
cellchat_all_N <- aggregateNet(cellchat_all_N)
groupSize <- as.numeric(table(cellchat_all_N@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_all_N@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_all_N@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#每一种细胞的互作
mat <- cellchat_all_N@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}


pathways.show <- c("SPP1") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_all_N, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_all_N, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_all_N, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_all_N, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat_all_N, sources.use =1:13, targets.use = 6:8, signaling =pathways.show, remove.isolate = FALSE)

netVisual_bubble(cellchat_all_N,sources.use =c(3,10,11,16,17,20,26,27), targets.use =12:14,signaling = c("ANGPT","ANGPTL","LIFR","GDF","TRAIL","IL6","IL10","NT","ANNEXIN","BAFF","BAG","BMP","CALCR","CCL","CHEMERIN","COMPLEMENT","CSF","CXCL","EDN","EGF","FGF","GALECTIN",
                                                                                                         "GAS","GRN","IFN-II","IGF","IL16","IL2","LIGHT","LT","MIF","MK","ncWNT","OSM",
                                                                                                        "PARs","PDGF","PERIOSTIN","PROS","PTN","SEMA3","SPP1","TGFb","TNF","TWEAK", "VEGF","VISFATIN","WNT"), remove.isolate = FALSE)
save(cellchat_all_N, file = "cellchat_all_N3.Rda")

#==============================================================================================================================
#IUA-myogenic-immune-All_IUA_proliferative_cellchat_IUA================================================================
All_IUA_proliferative_cellchat_IUA <- RunUMAP(All_IUA_proliferative_cellchat_IUA, reduction = "pca", dims = 1:8)
All_IUA_proliferative_cellchat_IUA <- FindNeighbors(All_IUA_proliferative_cellchat_IUA, reduction = "pca", dims = 1:8)
All_IUA_proliferative_cellchat_IUA <- FindClusters(All_IUA_proliferative_cellchat_IUA, resolution = 0.3)
# Visualization-UMAP
p1 <- DimPlot(All_IUA_proliferative_cellchat_IUA, reduction = "umap", group.by = "stim" ,pt.size = 0.5,raster =F)
#p2 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE, pt.size = 0.1,raster =F)
#p3 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE,group.by = "cell_type" ,pt.size = 0.1,raster =F)
p4 <- DimPlot(All_IUA_proliferative_cellchat_IUA, reduction = "umap", group.by = "datasets" ,pt.size = 0.5,raster =F)
p5 <- DimPlot(All_IUA_proliferative_cellchat_IUA, reduction = "umap", group.by = "cell_types5" ,pt.size = 0.5,raster =F)
#p6 <- DimPlot(All_IUA_proliferative, reduction = "umap", label = TRUE,group.by = "cell_cluster" ,pt.size = 0.1,raster =F)
#original_ID
p7 <- DimPlot(All_IUA_proliferative_cellchat_IUA, reduction = "umap", group.by = "original_ID" ,pt.size = 0.5,raster =F)
p8 <- DimPlot(All_IUA_proliferative_cellchat_IUA, reduction = "umap", label = F,group.by = "cell_cluster6" ,pt.size = 0.5,raster =F)

plot_grid(p7)
plot_grid(p1)
#plot_grid(p2)
#plot_grid(p3)
plot_grid(p4)
plot_grid(p5)
#plot_grid(p6)
plot_grid(p8)

library(CellChat)
#Extract the CellChat input files from a Seurat object
#建立cellchat对象
cellchat_all_IUA <- createCellChat(All_IUA_proliferative_cellchat_IUA@assays$RNA@data,meta=All_IUA_proliferative_cellchat_IUA@meta.data,group.by="cell_cluster6")
groupSize<-as.numeric(table(cellchat_all_IUA@idents))
groupSize

#导入配体受体数据库
CellChatDB<-CellChatDB.human
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# set the used database in the object
cellchat_all_IUA@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_all_IUA <- subsetData(cellchat_all_IUA) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat_all_IUA <- identifyOverExpressedGenes(cellchat_all_IUA)
cellchat_all_IUA <- identifyOverExpressedInteractions(cellchat_all_IUA)

# Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
cellchat_all_IUA <- computeCommunProb(cellchat_all_IUA)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_all_IUA <- filterCommunication(cellchat_all_IUA, min.cells = 10)
df.net<-subsetCommunication(cellchat_all_IUA)
write.csv(df.net,"net_LR-All_IUA_proliferative_cellchat_IUA3.csv")

#Infer the cell-cell communication at a signaling pathway level
cellchat_all_IUA <- computeCommunProbPathway(cellchat_all_IUA)
df.netp<-subsetCommunication(cellchat_all_IUA,slot.name = "netP")
write.csv(df.netp,"net_pathway_All_IUA_proliferative_cellchat_IUA3.csv")

#Calculate the aggregated cell-cell communication network
#所有细胞总体互作数量和强度
cellchat_all_IUA <- aggregateNet(cellchat_all_IUA)
groupSize <- as.numeric(table(cellchat_all_IUA@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_all_IUA@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_all_IUA@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#每一种细胞的互作
mat <- cellchat_all_IUA@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}


pathways.show <- c("SPP1") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_all_IUA, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_all_IUA, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_all_IUA, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_all_IUA, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat_all_IUA,sources.use =c(1,2,5,12,13,18,19,22,28,29), targets.use =14:16,signaling = c("ANGPT","ANGPTL","LIFR","GDF","TRAIL","IL6","IL10","NT","ANNEXIN","BAFF","BAG","BMP","CALCR","CCL","CHEMERIN","COMPLEMENT","CSF","CXCL","EDN","EGF","FGF","GALECTIN",
                                                                                                               "GAS","GRN","IFN-II","IGF","IL16","IL2","LIGHT","LT","MIF","MK","ncWNT","OSM",
                                                                                                               "PARs","PDGF","PERIOSTIN","PROS","PTN","SEMA3","SPP1","TGFb","TNF","TWEAK", "VEGF","VISFATIN","WNT"), remove.isolate = FALSE)


netVisual_bubble(cellchat_all_IUA,sources.use =12, targets.use =c(14,15,16), remove.isolate = FALSE)
netVisual_bubble(cellchat_all_N,sources.use =12, targets.use =c(14,15,16), remove.isolate = FALSE)
c(1,2,5,12,13,18,19,22,28,29)
c(3,4,6:10,14:16,19:26,29)

library(NMF)
library(ggalluvial)
# show all the interactions sending from Inflam.FIB
selectK(cellchat_all_IUA, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat_all_IUA, pattern = "outgoing", k = nPatterns)

selectK(cellchat_all_IUA, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat_all_IUA, pattern = "incoming", k = nPatterns)

save(cellchat_all_IUA, file = "cellchat_all_IUA3.Rda")

#merge#############################################################
# Compute the network centrality scores
#https://cloud.tencent.com/developer/article/2146523
cellchat_all_N_2 <- netAnalysis_computeCentrality(cellchat_all_N, slot.name = "netP")
cellchat_all_IUA_2 <- netAnalysis_computeCentrality(cellchat_all_IUA, slot.name = "netP")

#https://cloud.tencent.com/developer/article/1935670
object.list <- list(NORMAL = cellchat_all_N_2, IUA = cellchat_all_IUA_2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#比较交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#比较不同细胞群之间的相互作用数量和强度
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))}

#不同细胞类型之间相互作用或交互强度的差异
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

#显示每个数据集中任意两个细胞类型之间的交互次数或交互强度。
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

#比较 2D 空间中的主要来源和目标
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#识别保守和环境特异的信号通路
#根据信号组的功能相似性识别信号组
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#基于结构相似性识别信号组
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2

#计算和可视化通路距离
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2

#比较每个信号通路的整体信息流
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 12, height = 13)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1],width = 12, height = 13)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 12, height = 13, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 12, height = 13, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 12, height = 13, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 12, height = 13, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#识别上调和下调的信号配体对
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IUA", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IUA", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

#使用层次结构图、圆图或和弦图可视比较细胞-细胞通信
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))}

pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))}
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

#比较不同数据集之间的信号基因表达分布
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NORMAL", "IUA")) # set factor level
plotGeneExpression(cellchat, signaling = c("SPP1","TGFb","GALECTIN"),split.by = "datasets", colors.ggplot = T)
#> The default behaviour of split.by has changed.
#> Separate violin plots are now plotted side-by-side.
#> To restore the old behaviour of a single split violin,



################################################################
################################################################
##rat-iua################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(reshape2)
setwd("E:/figureX_fibrosis/iua_rat")

#regeneration_mice=======================================================================================
sham <- CreateSeuratObject(Read10X('Control'),"control")
sham$stim <- "control"
sham[["percent.mt"]] <- PercentageFeatureSet(sham , pattern = "^mt-")
summary(sham $nCount_RNA)
summary(sham $nFeature_RNA)
VlnPlot(sham , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
sham  <- subset(sham , subset = nCount_RNA < 50000 & nFeature_RNA <5000 & percent.mt<15)
VlnPlot(sham , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
sham  <- NormalizeData(sham , verbose = FALSE)
sham  <- FindVariableFeatures(sham , selection.method = "vst", nfeatures = 2000)

injury <- CreateSeuratObject(Read10X('AS'),"IUA")
injury$stim <- "IUA"
injury[["percent.mt"]] <- PercentageFeatureSet(injury , pattern = "^mt-")
summary(injury $nCount_RNA)
summary(injury $nFeature_RNA)
VlnPlot(injury , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
injury  <- subset(injury , subset = nCount_RNA < 50000 & nFeature_RNA <5000 & percent.mt<15)
VlnPlot(injury , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
injury  <- NormalizeData(injury , verbose = FALSE)
injury  <- FindVariableFeatures(injury , selection.method = "vst", nfeatures = 2000)

regeneration.anchors <- FindIntegrationAnchors(object.list = list(sham, injury), dims = 1:20)
regeneration <- IntegrateData(anchorset = regeneration.anchors, dims = 1:20)


cellinfo<-read.table('all_cluster_markers2.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
regeneration<- AddMetaData(object = regeneration, metadata = cellinfo)


VlnPlot(regeneration , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)
regeneration  <- subset(regeneration , subset =nCount_RNA < 50000 & nFeature_RNA>650 & nFeature_RNA <7500 & percent.mt<5)
VlnPlot(regeneration , features = c("nCount_RNA","nFeature_RNA",  "percent.mt"), ncol = 3)

###############################################################################################
###############################################################################################

regeneration <- ScaleData(regeneration, verbose = FALSE)
regeneration <- RunPCA(regeneration, npcs = 30, verbose = FALSE)
ElbowPlot(regeneration)

regeneration <- RunUMAP(regeneration, reduction = "pca", dims = 1:13)
regeneration <- FindNeighbors(regeneration, reduction = "pca", dims = 1:13)
regeneration <- FindClusters(regeneration, resolution = 0.3)
# Visualization-UMAP
#p1 <- DimPlot(regeneration, reduction = "umap", group.by = "stim" ,cols = c("#DC143C","#98FB98","#1E90FF","#9400D3"),pt.size = 0.5)
p1 <- DimPlot(regeneration, reduction = "umap", group.by = "group" ,cols = c("#DC143C","#1E90FF"),pt.size = 0.5)
p2 <- DimPlot(regeneration, reduction = "umap", label = TRUE, pt.size = 0.5)
#p3 <- DimPlot(regeneration, reduction = "umap", label = TRUE, group.by = "cell_cluster" ,pt.size = 0.5)
p4 <- DimPlot(regeneration, reduction = "umap", group.by = "cell_type" ,cols = c("#20B2AA","#0000FF","#FFA500","#9370DB","#DC143C"),pt.size = 0.5)
#p5 <- DimPlot(regeneration, reduction = "umap", label = TRUE,group.by = "cell_cluster" ,pt.size = 0.5)
p6 <- DimPlot(regeneration, reduction = "umap", label = TRUE,group.by = "cell_cluster" ,cols =c("#DC143C","#0000FF","#20B2AA","#FFA500",
                                                                                                "#FF00FF","#98FB98","#9370DB","#F08080",
                                                                                                "#1E90FF","#7CFC00","#FFFF00","#808000",
                                                                                                "#9400D3","#800080","#A0522D",
                                                                                                "#D2691E","#87CEEB","#40E0D0",
                                                                                                "#FF1493","#008B8B","#FFE4B5",
                                                                                                "#8A2BE2","#228B22"),pt.size = 0.5)


#DimPlot(regeneration, reduction = "umap", group.by = "cell_type" ,split.by = "group",cols = c("#FFA500","#0000FF","#20B2AA","#DC143C","#9370DB"),pt.size = 0.5)

plot_grid(p1)
plot_grid(p2)
plot_grid(p6)
plot_grid(p4)
#plot_grid(p5)
#plot_grid(p6)
#heatmap
Idents(object = regeneration) <- 'cell_type'

regeneration.markers <- FindAllMarkers(regeneration, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
regeneration.markers %>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(regeneration, features = top10$gene) + NoLegend()

DoHeatmap(regeneration, 
          features = top10$gene,label = F,
          group.colors =c("#9370DB","#DC143C","#0000FF","#FFA500","#20B2AA"))+
  scale_fill_gradientn (colors = c ("white","grey","firebrick3"))



###############################################################################################
#epithelia##############################################################################################
epithelia <- ScaleData(epithelia, verbose = FALSE)
epithelia <- RunPCA(epithelia, npcs = 30, verbose = FALSE)

cellinfo3<-read.table('epithelia_cellcluster.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
epithelia<- AddMetaData(object = epithelia, metadata = cellinfo3)

ElbowPlot(epithelia)

epithelia <- RunUMAP(epithelia, reduction = "pca", dims = 1:20)
epithelia <- FindNeighbors(epithelia, reduction = "pca", dims = 1:20)
epithelia <- FindClusters(epithelia, resolution = 0.6)
# Visualization-UMAP
p1 <- DimPlot(epithelia, reduction = "umap", group.by = "group" ,cols = c("#DC143C","#1E90FF"),pt.size = 1)
#p2 <- DimPlot(epithelia, reduction = "umap", label = T, pt.size = 1)
p3 <- DimPlot(epithelia, reduction = "umap", group.by = "cell_cluster" ,cols = c("#DC143C","#20B2AA","#9370DB","#0000FF","#FFA500"
),pt.size = 1)

plot_grid(p1)
#plot_grid(p2)
plot_grid(p3)

DimPlot(epithelia, reduction = "umap", split.by = "group")

#heatmap
epithelia.markers <- FindAllMarkers(epithelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
epithelia.markers %>%group_by(cluster) %>%top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(epithelia, features = top10$gene,label = F) 

DoHeatmap(epithelia, features = top10$gene,label = T) 

Idents(object = epithelia) <- 'cell_cluster'

write.table(epithelia.markers,file="epithelia.markers2024.txt")

VlnPlot(epithelia, features = c("Acta2","Vim","Col1a1","Lrp2","Mmp7","Mki67","Krt8","Epcam","Sox9","Foxa2"),
        cols = c("#9370DB","#DC143C","#20B2AA","#0000FF","#FFA500"),pt.size = 0,ncol = 3)

FeaturePlot(epithelia, features = c("Krt8","Epcam","Mmp7","Sox9","Col1a1","Mki67","Vim","Acta2","Foxa2","Lrp2"), min.cutoff = "q9")

FeaturePlot(epithelia, features = c("Krt8","Wfdc2","Mmp7","Sox9","Foxa2","Col1a1","Spp1","Mki67","Vim","Csf1"), min.cutoff = "q9")

FeaturePlot(epithelia, features = c("Krt8","Wfdc2","Ncam1","Mmp7","Sox9","Foxa2","Mki67","Vim","Col1a1"), min.cutoff = "q9")
FeaturePlot(epithelia, features = c("Esr1","Pgr","Ncam1","Mmp7","Sox9","Foxa2","Mki67","Vim","Col1a1"), min.cutoff = "q9")
epithelia <- RenameIdents(epithelia, `0` = "0", `1` = "0", `2` = "2", 
                          `3` = "3", `4`= "4", `5` = "5")
DimPlot(epithelia, label = T,pt.size = 0.5)


write.table(epithelia$seurat_clusters,file="epithelia—cluster2.txt")
write.table(epithelia$stim,file="epithelia—stim.txt")

save(epithelia, file = "epithelia.Rda")
###############################################################################################
#STROMA##############################################################################################

cellinfo3<-read.table('stroma_CELLcluster_2024.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
stroma<- AddMetaData(object = stroma, metadata = cellinfo3)

ElbowPlot(stroma)

stroma <- RunUMAP(stroma, reduction = "pca", dims = 1:15)
stroma <- FindNeighbors(stroma, reduction = "pca", dims = 1:15)
stroma <- FindClusters(stroma, resolution = 0.2)
# Visualization-UMAP
p1 <- DimPlot(stroma, reduction = "umap", group.by = "group" ,cols = c("#DC143C","#1E90FF"),pt.size = 1)
#p2 <- DimPlot(stroma, reduction = "umap", label = T, pt.size = 1)
p3 <- DimPlot(stroma, reduction = "umap", group.by = "cell_cluster" ,cols = c("#0000FF","#20B2AA","#FFA500","#DC143C","#9370DB"),pt.size = 1)


plot_grid(p1)
#plot_grid(p2)
plot_grid(p3)

DimPlot(stroma, reduction = "umap", group.by = "cell_cluster" ,cols = c("#DC143C","#0000FF","#808000","#FFA500",
                                                                        "#9370DB","#20B2AA"),pt.size = 1)


DimPlot(stroma, reduction = "umap", split.by = "stim")

FeaturePlot(stroma, features = c("Col1a1","Esr1","Pgr","Sfrp4","Ngfr","Mki67","Igf1","Wnt5a","Hoxa10","Piezo2","Wnt4",
                                 "Mmp3","Cxcl2","Cxcl12","Cxcl14","Ccl4","Acta2"), min.cutoff = "q9")

VlnPlot(stroma, features = c("Col1a1","Vegfd","Esr1","Pgr","Sfrp4","Ngfr","Mki67","Igf1","Fst","Nog","Hoxa10","Piezo1","Piezo2","Wnt4","Wnt5a",
                             "Tnc","Saa3","Cxcl2","Il1b","Acta2"),pt.size = 0.1,ncol = 3)

VlnPlot(stroma, features = c("Col1a1","Sfrp4","Cd74","Mki67","Acta2","Wnt4"),pt.size = 0,cols = c("#DC143C","#0000FF","#FFA500","#20B2AA","#9370DB"),ncol = 3)


markers.to.plot <- c("Sfrp4","Igf1","Fst","Nog","Hoxa10","Piezo1","Piezo2","Col1a1","Mki67","Wnt4","Tnc","Saa3","Acta2","Cxcl2","Il1b","Esr1","Pgr","Ngfr","Wnt5a","Vegfd")
DotPlot(stroma, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_cluster")+
  RotatedAxis()

#TF
markers.to.plot <- c("Ets1","Pbx1","Sox4","Hoxa10","Cbx5","Ezh2","Gmnn","H2afx","Hells","Nucks1","Tead1","Top2a","Gata6","Klf6","Wt1","Atf3","Nlrp3","Spi1","Esr1","Pgr","Fosb","Hand2","Nr2f2","Nr4a1","Osr1","Irf8","NfiB")
DotPlot(stroma, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_cluster")+
  RotatedAxis()

#heatmap
Idents(object = stroma) <- 'cell_cluster'
stroma.markers <- FindAllMarkers(stroma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
stroma.markers %>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(stroma, features = top10$gene)

DoHeatmap(stroma, 
          features = top10$gene,label = F,
          group.colors =c("#DC143C","#0000FF","#808000","#FFA500","#9370DB","#20B2AA"))+
  scale_fill_gradientn (colors = c ("white","grey","firebrick3"))


write.table(stroma.markers,file="STROMA_.markers2024.txt")


write.table(stroma$seurat_clusters,file="stroma—cluster_2024.txt")
write.table(stroma$stim,file="stroma—stim.txt")

save(stroma, file = "stroma.Rda")


#endothelia##############################################################################################

cellinfo3<-read.table('endothelia_cellcluster_2024.txt',header = T,row.names = 1)
## Basic QC and selecting cells for further analysis
endothelia<- AddMetaData(object = endothelia, metadata = cellinfo3)

ElbowPlot(endothelia)

endothelia <- RunUMAP(endothelia, reduction = "pca", dims = 1:15)
endothelia <- FindNeighbors(endothelia, reduction = "pca", dims = 1:15)
endothelia <- FindClusters(endothelia, resolution = 0.2)
# Visualization-UMAP
p1 <- DimPlot(endothelia, reduction = "umap", group.by = "group" ,cols = c("#DC143C","#1E90FF"),pt.size = 1)
#p2 <- DimPlot(endothelia, reduction = "umap", label = T, pt.size = 1)
p3 <- DimPlot(endothelia, reduction = "umap", group.by = "cell_cluster" ,cols = c("#DC143C","#20B2AA","#FFA500",
                                                                                  "#9370DB"),pt.size = 1)

plot_grid(p1)
#plot_grid(p2)
plot_grid(p3)


Idents(object = endothelia) <- 'cell_cluster'

DimPlot(endothelia, reduction = "umap", group.by = "cell_cluster" ,split.by = "group",cols = c("#DC143C","#20B2AA","#FFA500",
                                                                                               "#9370DB"),pt.size = 1)


DimPlot(endothelia, reduction = "umap", split.by = "stim")

FeaturePlot(endothelia, features = c("Col1a1","Esr1","Pgr","Vwf","Mki67","Igf1","Cldn3","Cldn5",
                                     "Mmp3","Cxcl2","Cxcl12","Cxcl14","Ccl4","Acta2"), min.cutoff = "q9")

VlnPlot(endothelia, features = c("Vwf","Cldn5",
                                 "Cxcl12","Col1a1","Acta2","Mki67"),pt.size = 0,cols = c("#FFA500","#20B2AA","#DC143C",
                                                                                         "#9370DB"),ncol = 3)

markers.to.plot <- c("Sfrp4","Igf1","Fst","Nog","Hoxa10","Piezo1","Piezo2","Col1a1","Mki67","Wnt4","Tnc","Saa3","Acta2","Cxcl2","Il1b","Esr1","Pgr","Ngfr","Wnt5a","Vegfd")
DotPlot(endothelia, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_cluster")+
  RotatedAxis()

#TF
markers.to.plot <- c("Nr2f2","Hnf1b","Bptf","Ikzf1","Nfkb2","Ascl2",
                     "Pbx1","Dbp","Znf581","Runx1","Pou2f2","Pgr",
                     "Mitf","Fosl1","Maff","Hmga1","Tcf12","Klf9",
                     "Egr1","Tbx21","Mafb","Irf1","Fosb","Runx3",
                     "Fos","Tfe3","Foxo1","Smad4","Smad1","Pbx3",
                     "Spib","Stat1","Tcf7l2","Znf419")
DotPlot(endothelia, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_cluster")+
  RotatedAxis()



#heatmap
endothelia.markers <- FindAllMarkers(endothelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
endothelia.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <-endothelia.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(endothelia, features = top10$gene,label = F) 



Idents(object = endothelia) <- 'cell_cluster'
endothelia.markers <- FindAllMarkers(endothelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
endothelia.markers %>%group_by(cluster) %>%top_n(n = 20, wt = avg_log2FC) -> top10
DoHeatmap(endothelia, features = top10$gene,label = F) 


DoHeatmap(endothelia, 
          features = top10$gene,label = F,
          group.colors =c("#DC143C","#0000FF","#808000","#FFA500","#9370DB","#20B2AA"))+
  scale_fill_gradientn (colors = c ("white","grey","firebrick3"))


write.table(endothelia.markers,file="endothelia_.markers2024.txt")

write.table(endothelia$seurat_clusters,file="endothelia—cluster_2024.txt")
write.table(endothelia$stim,file="endothelia—stim.txt")

save(endothelia, file = "endothelia.Rda")

#IMMUNE##############################################################################################

regeneration <- SetIdent(regeneration,cells = names(regeneration@meta.data$cell_type),value = regeneration@meta.data$cell_type)
immune <- subset(regeneration,idents = c("IMMUNE_CELL"))
save(immune,file = "immune.Rdata")


ElbowPlot(immune)
#==================================================================
immune <- RunUMAP(immune, reduction = "pca", dims = 1:14)
immune <- FindNeighbors(immune, reduction = "pca", dims = 1:14)
immune <- FindClusters(immune, resolution = 0.3)
# Visualization-UMAP
p1 <- DimPlot(immune, reduction = "umap", group.by = "group" ,pt.size = 0.5)
p2 <- DimPlot(immune, reduction = "umap", label = T, pt.size = 0.5)
p3 <- DimPlot(immune, reduction = "umap", group.by = "cell_cluster", label = T, pt.size = 0.5)

plot_grid(p1)
#plot_grid(p2)
plot_grid(p3)

p1 <- DimPlot(immune, reduction = "umap", group.by = "group" ,cols = c("#DC143C","#1E90FF"),pt.size = 0.5)
#p2 <- DimPlot(immune, reduction = "umap", label = T, pt.size = 0.5)
p3 <- DimPlot(immune, reduction = "umap", group.by = "cell_cluster" ,cols = c("#FFA500","#0000FF","#DC143C",
                                                                              "#9370DB"),pt.size = 0.5)

plot_grid(p1)
#plot_grid(p2)
plot_grid(p3)

DimPlot(immune, reduction = "umap", group.by = "cell_cluster" ,split.by = "group",cols = c("#808000","#DC143C","#0000FF","#FFA500",
                                                                                           "#9370DB","#20B2AA"),pt.size = 0.5)





#heatmap
IMMUNE.markers <- FindAllMarkers(immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
IMMUNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <-IMMUNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(immune, features = top10$gene,label = F) 



Idents(object = immune) <- 'cell_cluster'
IMMUNE.markers <- FindAllMarkers(immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
IMMUNE.markers %>%group_by(cluster) %>%top_n(n = 20, wt = avg_log2FC) -> top10
DoHeatmap(immune, features = top10$gene) 

VlnPlot(immune, features = c("Hbb","Cd68","Csf1r","Cxcr2","Mrc1","Spp1","Cd3e","Nkg7","Cd14","Ccl3","Cxcl2"),
        cols = c("#9370DB","#0000FF","#DC143C","#FFA500"),pt.size = 0,ncol = 3)

VlnPlot(immune, features = c("Hbb","Cxcr2","Nkg7","Cd68","Mrc1","Spp1"),
        cols = c("#9370DB","#0000FF","#DC143C","#FFA500"),pt.size = 0,ncol = 3)

plots <- VlnPlot(immune, features = c("Spp1"), split.by = "group", group.by = "cell_cluster", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(immune, features = c("Cd3e","Il7r","Nkg7","Cd19","Cd14","S100a8","Cd68","Csf1r","Csf3r","Il1b","Mrc1","Spp1","Igkc","Mki67","Cst3"), min.cutoff = "q9")

DimPlot(immune, reduction = "umap", split.by = "stim")


write.table(immune$seurat_clusters,file="immune——cluster2024.txt")
write.table(immune$group,file="immune—group.txt")



FeaturePlot(immune, features = c("Cd3e","Il7r","Nkg7","Mki67","Cd14","S100a8","Csf3r","Il1b","Cd68","Csf1r","Mrc1","Spp1","Cst3","Igkc","Igha"), min.cutoff = "q9")

VlnPlot(immune, features = c("Cd3e","Il7r","Nkg7","Mki67","Cd14","S100a8","Csf3r","Il1b","Cd68","Csf1r","Mrc1","Spp1","Cst3","Igkc","Igha"),cols = c("#20B2AA","#DC143C","#0000FF","#FFA500","#9370DB","#808000"),pt.size = 0,ncol = 3)

markers.to.plot <- c("Cd3e","Il7r","Nkg7","Mki67","Cd14","S100a8","Csf3r","Il1b","Cd68","Csf1r","Mrc1","Spp1","Cst3","Igkc","Igha")
DotPlot(immune, features = markers.to.plot, cols =  c("grey", "firebrick3"), dot.scale = 8, group.by ="cell_cluster")+
  RotatedAxis()


#heatmap
Idents(object = immune) <- 'cell_cluster'
immune.markers <- FindAllMarkers(stroma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(immune, features = top10$gene) + NoLegend()

DoHeatmap(immune, 
          features = top10$gene,label = F,
          group.colors =c("#20B2AA","#DC143C","#0000FF","#FFA500","#9370DB","#808000"))+
  scale_fill_gradientn (colors = c ("white","grey","firebrick3"))


write.table(IMMUNE.markers,file="immune_.markers2024.txt")

save(immune, file = "immune_2024.Rda")

###############################################################################################

##SCNIC#############################################################################################
rm(list = ls()) 
regeneration <- SetIdent(regeneration,cells = names(regeneration@meta.data$cell_type),value = regeneration@meta.data$cell_type)
stroma <- subset(regeneration,idents = c("stroma"))
save(stroma,file = "stroma.Rdata")
rm(regeneration)
rm(stroma)


library(Seurat) 
library(SCENIC)
## Load data表达矩阵加上表型信息
exprMat  <-  as.matrix(endothelia@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  endothelia@meta.data
save(exprMat,cellInfo,file = "endothelia_scenic.Rdata")

### Initialize settings
org<-"mgi"
dbDir<-"cisTarget_databases"

myDatasetTitle<-"SCENIC analysis on regeneration stroma cells"
data("defaultDbNames")
dbs<-defaultDbNames[[org]]
#scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

motifAnnotations_mgi = motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
#创建一个对象：scenicOptions，后续数据会保存在当前目录的int文件夹
#输出结果会存储在output文件夹
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#将表达矩阵中未在cisTarget_databases收录的基因去除
library(RcisTarget)
dbFilePath<-getDatabases(scenicOptions)[[1]]
motifRanking<-importRankings(dbFilePath)
geneInDatabase<-colnames(getRanking(motifRanking))
genelefet_minCells<-rownames(exprMat)
length(genelefet_minCells)
genelefet_minCells_inDatabases<-genelefet_minCells[which(genelefet_minCells %in% geneInDatabase)]
length(genelefet_minCells_inDatabases)
genesKept<-genelefet_minCells_inDatabases
exprMat_filter<-exprMat[genesKept,]
dim(exprMat)
dim(exprMat_filter)

#计算相关性并对数据进行log处理
### Co-expression network
class(exprMat_filter)
runCorrelation(as.matrix(exprMat_filter), scenicOptions)
exprMat_filter <- log2(exprMat_filter+1) 
library(GENIE3)
runGenie3(as.matrix(exprMat_filter), scenicOptions)
save(exprMat_filter,scenicOptions,file="endothelia_input_GENIE3_data.Rdata")


### Build the GRN
#load(stroma_scenic.Rdata)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
closeAllConnections()
scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#结果可视化
library(AUCell)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

#Regulators for known cell types or clusters

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cell_cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, col = col,name="Regulon activity")

#############################################################################

library(Seurat) 
library(SCENIC)
## Load data表达矩阵加上表型信息
exprMat  <-  as.matrix(stroma@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  stroma@meta.data
save(exprMat,cellInfo,file = "stroma_scenic.Rdata")

### Initialize settings
org<-"mgi"
dbDir<-"cisTarget_databases"

myDatasetTitle<-"SCENIC analysis on regeneration stroma cells"
data("defaultDbNames")
dbs<-defaultDbNames[[org]]
#scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

motifAnnotations_mgi = motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
#创建一个对象：scenicOptions，后续数据会保存在当前目录的int文件夹
#输出结果会存储在output文件夹
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#将表达矩阵中未在cisTarget_databases收录的基因去除
library(RcisTarget)
dbFilePath<-getDatabases(scenicOptions)[[1]]
motifRanking<-importRankings(dbFilePath)
geneInDatabase<-colnames(getRanking(motifRanking))
genelefet_minCells<-rownames(exprMat)
length(genelefet_minCells)
genelefet_minCells_inDatabases<-genelefet_minCells[which(genelefet_minCells %in% geneInDatabase)]
length(genelefet_minCells_inDatabases)
genesKept<-genelefet_minCells_inDatabases
exprMat_filter<-exprMat[genesKept,]
dim(exprMat)
dim(exprMat_filter)

#计算相关性并对数据进行log处理
### Co-expression network
class(exprMat_filter)
runCorrelation(as.matrix(exprMat_filter), scenicOptions)
exprMat_filter <- log2(exprMat_filter+1) 
library(GENIE3)
runGenie3(as.matrix(exprMat_filter), scenicOptions)
save(exprMat_filter,scenicOptions,file="stroma_input_GENIE3_data.Rdata")


### Build the GRN
#load(stroma_scenic.Rdata)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
closeAllConnections()
scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#结果可视化
library(AUCell)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

#Regulators for known cell types or clusters

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cell_cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, col = col,name="Regulon activity")

####################################################################################################

library(Seurat) 
library(SCENIC)
## Load data表达矩阵加上表型信息
exprMat  <-  as.matrix(epithelia@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  epithelia@meta.data
save(exprMat,cellInfo,file = "epithelia_scenic.Rdata")

### Initialize settings
org<-"mgi"
dbDir<-"cisTarget_databases"

myDatasetTitle<-"SCENIC analysis on regeneration epithelia cells"
data("defaultDbNames")
dbs<-defaultDbNames[[org]]
#scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

motifAnnotations_mgi = motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,dbs=dbs['10kb'],datasetTitle =myDatasetTitle,  nCores=10)
#创建一个对象：scenicOptions，后续数据会保存在当前目录的int文件夹
#输出结果会存储在output文件夹
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#将表达矩阵中未在cisTarget_databases收录的基因去除
library(RcisTarget)
dbFilePath<-getDatabases(scenicOptions)[[1]]
motifRanking<-importRankings(dbFilePath)
geneInDatabase<-colnames(getRanking(motifRanking))
genelefet_minCells<-rownames(exprMat)
length(genelefet_minCells)
genelefet_minCells_inDatabases<-genelefet_minCells[which(genelefet_minCells %in% geneInDatabase)]
length(genelefet_minCells_inDatabases)
genesKept<-genelefet_minCells_inDatabases
exprMat_filter<-exprMat[genesKept,]
dim(exprMat)
dim(exprMat_filter)

#计算相关性并对数据进行log处理
### Co-expression network
class(exprMat_filter)
runCorrelation(as.matrix(exprMat_filter), scenicOptions)
exprMat_filter <- log2(exprMat_filter+1) 
library(GENIE3)
runGenie3(as.matrix(exprMat_filter), scenicOptions)
save(exprMat_filter,scenicOptions,file="epithelia_input_GENIE3_data.Rdata")


### Build the GRN
#load(epithelia_scenic.Rdata)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
closeAllConnections()
scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#结果可视化
library(AUCell)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

#Regulators for known cell types or clusters

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cell_cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")




