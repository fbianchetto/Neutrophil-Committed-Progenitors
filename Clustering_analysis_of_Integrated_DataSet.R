#####################################################################
### 2. Pipeline Clustering analysis of integrated dataset       #####
#####################################################################

#Load packages
source("/media/pg/Volume2_4T/RNA_SEQ_our/SingleCells/SingleCellsUtils.R")
options(future.globals.maxSize = 20000 * 1024^2)
library(Seurat)
library(ggplot2)
library(scater)
library(SingleR)
library(RColorBrewer)

MyUMAPPlot <- function(seurat,label=NULL,donors=NULL){
  library(ggthemes)
  umap = as.data.frame(Embeddings(object = seurat, reduction = "umap"))
  umap$Cells= seurat$Sample_Name
  umap$donor= paste("D",seurat$donor,sep="")
  data = subset(umap,umap$Cells==label)
  data = subset(data,data$donor%in%donors)
  p1 <- ggplot(umap, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(colour = "gray", size = 0.1) + geom_rangeframe() + theme_tufte() + geom_point(data=data,colour="darkblue", size = 0.1)+ ggtitle(label,donors)
  return(p1)

}

########################################
### Intergrated seurat object     ######
########################################

BMs.integrated = readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/precursors.integrated_BD_p1_p99_Mt_CC_genes.rds")
BMs.integrated <- RunPCA(BMs.integrated, verbose = FALSE,npcs = 50)
#### Use elbow plot to define number of dimensions
ElbowPlot(BMs.integrated,ndims = 50)
print(BMs.integrated[["pca"]], dims = 1:50, nfeatures = 20)
DimHeatmap(BMs.integrated, dims = 1:12, cells = 500, balanced = TRUE)
#### Run UMAP
BMs.integrated <- RunUMAP(BMs.integrated,dims = 1:50,reduction.name = "umap_learn50PC",reduction.key = "UMAPLearn50PC_",umap.method="umap-learn",metric = "correlation")

#### Plot UMAP
p1 <- DimPlot(object = BMs.integrated,
              reduction = 'umap_learn50PC',label = FALSE,
              group.by = "donor")+ scale_color_manual(values = color_clusters) + ggtitle("Group by donor,50PC");p1

#### Find neighbours and calculate clusters
ElbowPlot(BMs.integrated,ndims = 50)
PCA = 50
usepcs = seq(1,PCA,by =1)

#### Function to find clusters at several resolution
OptimalClusterRes <- function(seurat,res){
  for(i in res){
    seurat <- FindClusters(seurat, verbose = FALSE,resolution = i,algorithm = 1)
  }
  return(seurat)
}
res = seq(from=0,to=0.5,by=0.1)
#### Findneighbours
BMs.integrated <- FindNeighbors(BMs.integrated, dims = usepcs, verbose = FALSE, force.recalc = FALSE)
BMs.integrated = OptimalClusterRes(BMs.integrated,res)
library(clustree)
cl1 <- clustree(BMs.integrated, prefix = "integrated_snn_res.")
cl1

#### Visual inspection of clusters at several resolution
reduction = "umap_learn50PC"
R0.1 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                group.by="integrated_snn_res.0.1",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("integrated_snn_res.0.1") + theme(plot.title = element_text(size = 10, face = "bold"))

R0.2 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                group.by="integrated_snn_res.0.2",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("integrated_snn_res.0.2") + theme(plot.title = element_text(size = 10, face = "bold"))

R0.3 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                group.by="integrated_snn_res.0.3",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("integrated_snn_res.0.3") + theme(plot.title = element_text(size = 10, face = "bold"))

R0.4 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                group.by="integrated_snn_res.0.4",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("integrated_snn_res.0.4") + theme(plot.title = element_text(size = 10, face = "bold"))

R0.5 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                group.by="integrated_snn_res.0.5",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("integrated_snn_res.0.5") + theme(plot.title = element_text(size = 10, face = "bold"))

multiplot(R0.1,R0.2,R0.3,R0.4,R0.5,cols = 2)
#############################################################################################################
res = seq(from=0.3,to=0.5,by=0.1)
k = c(10,20,30)

### Function to find clusters at several resolution and k.param 
OptimalClusterRes <- function(seurat,res,k){
  for(i in k){
    for(j in res){
      seurat <- FindNeighbors(seurat, dims = usepcs, verbose = TRUE, force.recalc = FALSE,k.param=i)
      seurat <- FindClusters(seurat, verbose = TRUE,resolution = j,algorithm = 1)
      ClusterNames = paste("Clusters","k",i,"res",j,sep="_")
      print(ClusterNames)
      seurat <- StashIdent(object = seurat, save.name = ClusterNames)
    }
  }
  return(seurat)
}

DefaultAssay(BMs.integrated) <- "integrated"
BMs.integrated = OptimalClusterRes(BMs.integrated,res,k)
table(BMs.integrated$integrated_snn_res.0.3,BMs.integrated$Clusters_k_20_res_0.3)

control = data.frame(BMs.integrated$integrated_snn_res.0.3,BMs.integrated$Clusters_k_20_res_0.3)
colnames(control) <- gsub("BMs.integrated.","",colnames(control))
#####
reduction = "umap_learn50PC"
K_10_res0.3 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                        group.by="Clusters_k_10_res_0.3",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_10_res_0.3") + theme(plot.title = element_text(size = 10, face = "bold"))

K_10_res0.4 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                        group.by="Clusters_k_10_res_0.4",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_10_res_0.4") + theme(plot.title = element_text(size = 10, face = "bold"))

K_10_res0.5 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                        group.by="Clusters_k_10_res_0.5",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_10_res_0.5") + theme(plot.title = element_text(size = 10, face = "bold"))

#######
K_20_res0.3 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                       group.by="Clusters_k_20_res_0.3",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_20_res_0.3") + theme(plot.title = element_text(size = 10, face = "bold"))

K_20_res0.4 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                       group.by="Clusters_k_20_res_0.4",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_20_res_0.4") + theme(plot.title = element_text(size = 10, face = "bold"))

K_20_res0.5 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                       group.by="Clusters_k_20_res_0.5",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_20_res_0.5") + theme(plot.title = element_text(size = 10, face = "bold"))
#########
K_30_res0.3 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                       group.by="Clusters_k_30_res_0.3",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_30_res_0.3") + theme(plot.title = element_text(size = 10, face = "bold"))

K_30_res0.4 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                       group.by="Clusters_k_30_res_0.4",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_30_res_0.4") + theme(plot.title = element_text(size = 10, face = "bold"))

K_30_res0.5 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                       group.by="Clusters_k_30_res_0.5",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_30_res_0.5") + theme(plot.title = element_text(size = 10, face = "bold"))

p1 <- multiplot(K_10_res0.3,K_10_res0.4,K_10_res0.5,
                K_20_res0.3,K_20_res0.4,K_20_res0.5,
                K_30_res0.3,K_30_res0.4,K_30_res0.5,cols = 3)


##### Re- run OptimalClusterRes function
OptimalClusterRes <- function(seurat,res){
  for(i in res){
    seurat <- FindClusters(seurat, verbose = FALSE,resolution = i,algorithm = 1)
  }
  return(seurat)
}
res = seq(from=0,to=0.5,by=0.1)

# Findneighbours
BMs.integrated <- FindNeighbors(BMs.integrated, dims = usepcs, verbose = FALSE, force.recalc = FALSE)
BMs.integrated = OptimalClusterRes(BMs.integrated,res)

p1 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
              group.by="Clusters_k_20_res_0.3",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("Clusters_k_20_res_0.3") + theme(plot.title = element_text(size = 10, face = "bold"))

p2 <- DimPlot(object = BMs.integrated,reduction = reduction,pt.size = 0.01,
                      group.by="integrated_snn_res.0.3",label = TRUE) + scale_color_manual(values=color_clusters) + NoLegend() +
  ggtitle("integrated_snn_res.0.3") + theme(plot.title = element_text(size = 10, face = "bold"))

multiplot(p1,p2,cols = 2)

#### Set identity of cells with the best parameters for clustering 
Idents(BMs.integrated) <-"Clusters_k_20_res_0.3"
reduction = "umap_learn50PC"

MyUMAPPlot_cluster <- function(seurat,cluster=NULL,donors=NULL){
  library(ggthemes)
  umap = as.data.frame(Embeddings(object = seurat, reduction = reduction))
  colnames(umap) <- c("UMAP_1","UMAP_2")
  umap$Cells= seurat$Sample_Name
  umap$clusters= Idents(seurat)
  data = subset(umap,umap$clusters==cluster)
  data$density <- get_density(data$UMAP_1, data$UMAP_2, n = 100)
  p1 <- ggplot(umap, aes(x=UMAP_1, y=UMAP_2)) + geom_point(colour = "gray", size = 0.01) + geom_rangeframe() + theme_tufte() +
    geom_point(data=data,aes(x=UMAP_1, y=UMAP_2, color = density),size = 0.1) + scale_color_viridis() + ggtitle(cluster,donors)
    return(p1)
}

cl0 = MyUMAPPlot_cluster(BMs.integrated, cluster="0")
cl1 = MyUMAPPlot_cluster(BMs.integrated, cluster="1")
cl2 = MyUMAPPlot_cluster(BMs.integrated, cluster="2")
cl3 = MyUMAPPlot_cluster(BMs.integrated, cluster="3")
cl4 = MyUMAPPlot_cluster(BMs.integrated, cluster="4")
cl5 = MyUMAPPlot_cluster(BMs.integrated, cluster="5")
cl6 = MyUMAPPlot_cluster(BMs.integrated, cluster="6")
cl7 = MyUMAPPlot_cluster(BMs.integrated, cluster="7")
cl8 = MyUMAPPlot_cluster(BMs.integrated, cluster="8")

p1 <- multiplot(cl0,cl1,cl2,
                cl3,cl4,cl5,
                cl6,cl7,
                cl8,cols = 3)

########################################
Idents(BMs.integrated) <-"Clusters_k_20_res_0.3"
#saveRDS(BMs.integrated, file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.50PC_Mt_CC.rds")
BMs.integrated <- readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.50PC_Mt_CC.rds")

#### Figure 6A
reduction = "umap_learn50PC"
cluster_col <- c("0" = "#E41A1C",#
                 "1" = "#00008B",#
                 "2" = "#006400",#
                 "3" = "#FDB462",#
                 "4" = "#87CEEB",
                 "5" = "#7CFC00",#
                 "6" = "#000000",
                 "7" = "#F781BF",
                 "8" =  "#999999")

R0.3 <- DimPlot(object = BMs.integrated,reduction = reduction,label = TRUE) + scale_color_manual(values=cluster_col) +
  ggtitle("integrated_snn_res.0.3") + theme(plot.title = element_text(size = 10, face = "bold")) + xlim(-5,12) + ylim(-5,8)

cluster_tag <- c("cMOP" = "#006400",
                 "N0" = "#E41A1C",
                 "N1" = "#FDB462",
                 "N2" = "#87CEEB",
                 "N3" = "#00008B")

p6 <- DimPlot(object = BMs.integrated,
              reduction = reduction,label = FALSE,
              group.by = "Sample_Name") + ggtitle("Group by donor,50PC") + xlim(-5,12) + ylim(-5,8) + scale_color_manual(values = cluster_tag);p6

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig6_B_UMAP_mt_CC_without_contaminats_new_colors_all_clusters.eps", width = 20, height = 10)
scater::multiplot(p6,R0.3,cols = 2)
dev.off()

######### Figure S7A
MyUMAPPlot <- function(seurat,label=NULL,donors=NULL){
  library(ggthemes)
  umap = as.data.frame(Embeddings(object = seurat, reduction = "umap_learn50PC"))
  colnames(umap) <- c("UMAP_1","UMAP_2")
  umap$Cells= seurat$Sample_Name
  data = subset(umap,umap$Cells==label)
  data$density <- get_density(data$UMAP_1, data$UMAP_2, n = 100)
  p1 <- ggplot(umap, aes(x=UMAP_1, y=UMAP_2)) + geom_point(colour = "gray", size = 0.01) + geom_rangeframe() + theme_tufte() +
    geom_point(data=data,aes(x=UMAP_1, y=UMAP_2, color = density),size = 0.1) + scale_color_viridis() + ggtitle(label,donors) + xlim(-5,12) + ylim(-5,7.5)
  return(p1)
}

N0 = MyUMAPPlot(BMs.integrated, label="N0")
N1 = MyUMAPPlot(BMs.integrated, label="N1")
N2 = MyUMAPPlot(BMs.integrated, label="N2")
N3 = MyUMAPPlot(BMs.integrated, label="N3")
cMOPs = MyUMAPPlot(BMs.integrated, label="cMOP")

p1 <- multiplot(N0,N1,N2, N3,cMOPs,cols = 2)

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/FigS7B_A_UMAP_all_genes.eps", width = 8, height = 10)
multiplot(N0,N1,N2, N3,cMOPs,cols = 2)
dev.off()

#############################################################################
####### Analysis of cells contamination
#############################################################################
DefaultAssay(BMs.integrated) = "RNA" ####
BMs.integrated <- NormalizeData(BMs.integrated)

eosinphlis = c("IL5RA","CLC","EPX")
erythocytes = c("HBB")
platelet = c("PF4")
B = c("IGHG3")

VlnPlot(BMs.integrated,features = eosinphlis,pt.size = 0.1)
VlnPlot(BMs.integrated,features = erythocytes,pt.size = 0)
VlnPlot(BMs.integrated,features = platelet,pt.size = 0)
VlnPlot(BMs.integrated,features = B,pt.size = 0.1)

eosinphlis1 = WhichCells(object = BMs.integrated, cells=NULL ,expression = EPX > 0)
eosinphlis2 = WhichCells(object = BMs.integrated, cells=NULL ,expression = CLC > 0)

B =  WhichCells(object = BMs.integrated, cells=NULL ,expression = IGHG3 > 0)
platelet = WhichCells(object = BMs.integrated, cells=NULL ,expression = PF4 > 0)
erythocytes = WhichCells(object = BMs.integrated, cells=NULL ,expression = HBB > 0)

### Eosinophils
BMs.integrated = subset(BMs.integrated, subset = EPX > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = IL5RA > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = CLC > 0, invert = TRUE)
### Platelet
BMs.integrated = subset(BMs.integrated, subset = PF4 > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = ITGA2B > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = LTBP1 > 0, invert = TRUE)
### B
BMs.integrated = subset(BMs.integrated, subset = IGHG1 > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = IGHG3 > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = IGHG4 > 0, invert = TRUE)
### Eritrocytes
BMs.integrated = subset(BMs.integrated, subset = HBB > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = HBA1 > 0, invert = TRUE)
BMs.integrated = subset(BMs.integrated, subset = CA1 > 0, invert = TRUE)
### T cells
BMs.integrated = subset(BMs.integrated, subset = TRBC1 > 0, invert = TRUE)
##########################################################################
BMs.integrated = subset(BMs.integrated, subset = UMAPLearn50PC_2 > 4.6, invert = TRUE)
BM_filtered = subset(x = BM_filtered, subset = UMAP_1 < 11.5 & UMAP_2 < 1.3, invert = TRUE)

eosinphlis = c("IL5RA","CLC","EPX")
erythocytes = c("HBB","HBA1","HBD","CA1")
platelet = c("PF4","ITGA2B","LTBP1")
B = c("IGHG1","IGHG3","IGHG4")
VlnPlot(BMs.integrated,features = eosinphlis,pt.size = 0.1)
VlnPlot(BMs.integrated,features = platelet,pt.size = 1)
VlnPlot(BMs.integrated,features = B,pt.size = 0.1)
VlnPlot(BMs.integrated,features = erythocytes,pt.size = 0.1)

saveRDS(BMs.integrated, file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.clean_50PC_Mt_CC.rds")
##########################################################################
###### Cleaned data
##########################################################################
BMs.integrated <- readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.clean_50PC_Mt_CC.rds")
##########################################################################
DefaultAssay(BMs.integrated) = "RNA" ####"integrated"
BMs.integrated <- NormalizeData(BMs.integrated)
markers <- FindAllMarkers(object = BMs.integrated, only.pos = TRUE, min.pct = 0.25, min.diff.pct = -Inf, logfc.threshold = 0.25)
markers <- subset(markers,markers$p_val_adj < 0.01)
data.frame(markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC))
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

library(WriteXLS)
WriteXLS::WriteXLS(top10,ExcelFileName = "/home/pg/Documents/Shared-Win10/Dropbox/Precursors/DEGs_all_cluters_Neups_and_cMop.xls")

#######
BMs.integrated <- readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.clean_50PC_Mt_CC.rds")
DefaultAssay(BMs.integrated) <- "RNA"

BMs.integrated = subset(BMs.integrated,idents=c('6'), invert = TRUE)
BMs.integrated = NormalizeData(BMs.integrated)
DimPlot(BMs.integrated,reduction="umap_learn50PC")

###### Plots Figure 6C
genes <- c("ELANE","MPO","CSF3R","CSF1R","IRF8","CEBPE","CD34","GFI1","AZU1","PRTN3","SRGN","CTSG")

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/ViolinPlot_Neu_cMOP_1.eps", width = 8, height = 5)
VlnPlot(BMs.integrated,features = genes[1:6],pt.size = 0)
dev.off()

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/ViolinPlot_Neu_cMOP_2.eps", width = 8, height = 5)
VlnPlot(BMs.integrated,features = genes[7:12],pt.size = 0)
dev.off()

######

genes <- c("IGLL1","SMIM24","CPA3","SUCNR1","MS4A3")
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/ViolinPlot_Neu_cMOP.eps", width = 8, height = 5)
VlnPlot(BMs.integrated,features = genes,pt.size = 0)
dev.off()

genes = c("IRF8","SAMHD1","ANXA2")
# ##### function from https://github.com/satijalab/seurat/issues/2475
#
# draw median for a single plot
MedianVioSeu <- function() stat_summary(fun = median, geom='point', size = 1, colour = "black")

# for a list of ggplots
apply.MedianVioSeu.and.combine <- function(ls.ggplots=pl$A) {
  for (i in 1:length(ls.ggplots)) ls.ggplots[[i]] = ls.ggplots[[i]]+MedianVioSeu()
  return(CombinePlots(ls.ggplots, legend = 'none'))
}
#
C <- c("IRF8","ANXA2","CSF1R","LY86","HIST1H4C")
pl=NULL
pl$C <- VlnPlot(object = BMs.integrated, features = C, group.by =  NULL, pt.size = 0, assay = "RNA", log = F, combine = F,ncol=3)
pl$C <- apply.MedianVioSeu.and.combine(pl$C)
pl$C

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/ViolinPlots_contaminants.eps", width = 10, height = 10)
pl$C
dev.off()


####### Find gene expression markers for clusters 
markers <- FindAllMarkers(object = BMs.integrated, only.pos = TRUE, min.pct = 0.25, min.diff.pct = -Inf, logfc.threshold = 0.25)
markers <- subset(markers,markers$p_val_adj < 0.01)
data.frame(markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC))
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)


### Funtion to get the average gene expression of markers for each cluster

SignGenesTableCounts = function(seurat, markers){
  clustL = length(unique(Idents(seurat)))  +1
  data <- c()
  for(i in unique(markers$cluster)){
    df = cluster_AverageExpression(seurat,markers,clust=i)#[,1:clustL]
    data = rbind(data,df)
    data
  }
  data = data[!duplicated(data$gene),]
  return(data)
}

data = SignGenesTableCounts(BMs.integrated,markers=markers)

library(WriteXLS)
WriteXLS(data,ExcelFileName = "/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/DEGs_Neups_and_cMop.xls")

rownames(data) = data$gene; dim(data)

##### Perform hierarchical clustering dendrogram for pseudo-bulk and plot seriate_dendrogram
data = data[,2:7]
data[grepl("BEX1",rownames(data)),]
#
d <- dist(t(data),method = "euclidean")
hc <- hclust(d,method = "ward.D")
library(ggdendro)
#dend = ggdendrogram(hc, rotate = FALSE, size = 2);dend
dend <- as.dendrogram(hc)
plot(dend)

dend2 <- dendextend::seriate_dendrogram(dend, d,method = c("OLO"))
plot(dend2)

###### Figure S7B
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/HCl_euclidean_ward.D.eps", width = 8, height = 3)
plot(dend2)
dev.off()


# ####### Pseudobulk analysis
# #### https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
#
BMs.integrated <- readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.clean_50PC_Mt_CC.rds")

BMs.c0 = subset(BMs.integrated,idents=c('0'), invert = FALSE)
BMs.c1 = subset(BMs.integrated,idents=c('1'), invert = FALSE)
BMs.c3 = subset(BMs.integrated,idents=c('3'), invert = FALSE)
BMs.c4 = subset(BMs.integrated,idents=c('4'), invert = FALSE)

GetHighestExpressedGenes <- function(seurat){
  # Extract raw counts and metadata to create SingleCellExperiment object
  counts <- seurat@assays$RNA@counts
  counts <- counts[!grepl("MT-",rownames(counts)),] ### removing mitochondrial genes
  counts <- counts[!grepl("RPS|RPL",rownames(counts)),] ### removing ribosolmal genes
  metadata <- seurat@meta.data
  # Set up metadata as desired for aggregation and DE analysis
  metadata$cluster_id <- factor(seurat@active.ident)
  # Create single cell experiment object
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = metadata)
  sce <- addPerCellQC(sce,subsets=list(Mito=grep("MT-", rownames(sce))))
  p <- plotHighestExprs(sce, exprs_values = "counts")
  return(p)
}

c0 <- GetHighestExpressedGenes(BMs.c0)
c1 <- GetHighestExpressedGenes(BMs.c1)
c3 <- GetHighestExpressedGenes(BMs.c3)
c4 <- GetHighestExpressedGenes(BMs.c4)

########
GetAverageExpression <- function(seurat){
  aveExp = AverageExpression(seurat,features = NULL, assays=DefaultAssay(seurat),return.seurat = TRUE)
  counts <- GetAssayData(aveExp)
  counts <- counts[!grepl("MT-",rownames(counts)),] ### removing mitochondrial genes
  counts <- counts[!grepl("RPS|RPL",rownames(counts)),] ### removing ribosolmal genes
  require(dplyr)
  require(biomaRt)
  ensembl <- useMart("ensembl", host = "http://oct2018.archive.ensembl.org",dataset = "hsapiens_gene_ensembl")
  annotation = getBM(attributes=c('hgnc_symbol','gene_biotype','description'),mart = ensembl,filters = 'external_gene_name',values=rownames(counts))
  colnames(annotation) = c("SYMBOL","biotype","name")
  annotation$name = gsub('\\[.*','',annotation$name)
  annotation = annotation[!duplicated(annotation$SYMBOL),]
  aveExp_cluster = as.data.frame(subset(counts,rownames(counts)%in%annotation$SYMBOL))
  annotation = subset(annotation,annotation$SYMBOL%in%rownames(aveExp_cluster))
  aveExp_cluster$gene = rownames(aveExp_cluster)
  data = merge(aveExp_cluster,annotation,by.x="gene", by.y="SYMBOL",sort = FALSE,all=TRUE)
  return(data)
}

aveExp = GetAverageExpression(BMs.integrated)
WriteXLS::WriteXLS(aveExp,ExcelFileName = "/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/allgenes_Neups_and_cMop.xls")



