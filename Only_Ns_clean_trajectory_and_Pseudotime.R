################################################################
#### 3. Pipeline of analisis of NMs population             #####
###############################################################
### Set colors
cluster_col <- c("0" = "#E41A1C",#
                 "1" = "#00008B",#
                 "2" = "#006400",#
                 "3" = "#FDB462",#
                 "4" = "#87CEEB",
                 "5" = "#7CFC00",#
                 "6" = "#000000",
                 "7" = "#F781BF",
                 "8" =  "#999999")
####################################

source("/media/pg/Volume2_4T/RNA_SEQ_our/SingleCells/SingleCellsUtils.R")
library(Seurat)
library(destiny)
library(RColorBrewer)

colors = c("#737373","#969696","#BDBDBD","#D9D9D9","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D","#8C2D04") ### grey and reds scale

options(future.globals.maxSize = 20000 * 1024^2)
source("/media/pg/Volume2_4T/RNA_SEQ_our/SingleCells/SingleCellsUtils.R")
Ns = readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.Ns.clean_50PC_Mt_CC.rds")
DimPlot(Ns,label=T,reduction = "umap_learn50PC")
length(Ns@assays$integrated@var.features)

my_levels <- c("0","3","4","1");length(my_levels) ### 0.4
Ns@active.ident <- factor(x=Ns@active.ident,levels = my_levels)

Marco_barplot(Ns)

Ns_new = Ns
Ns_new$cluster = Idents(Ns_new)
Idents(Ns_new) <- "Sample_Name"
Ns_new = subset(Ns_new,idents=c('cMOP'), invert = TRUE)
Idents(Ns_new) <- "cluster"

#### Figure 7 A
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig6C_barplot_without_contaminants_only_Ns.eps", width = 8, height = 5)
Marco_barplot(Ns_new)
dev.off()

DefaultAssay(Ns) <- "RNA"
Ns <- NormalizeData(Ns,verbose = FALSE)
VlnPlot(Ns,features = c("TOP2A","TUBB"),pt.size = 0)
FeaturePlot(Ns,features = c("TOP2A","TUBB"),reduction = "umap_learn30PC")

Primary_granules = list(c("ACPP","AZU1","BPI","CD63","CECAM6","CTSG","DEFA3","DEFA4","ELANE","GLA","GLB1","GNS","HEXA","MAN2B1","MPO","PRTN3","VNN1"))
Secondary_granules = list(c("ITGAM","LTF","C3AR1","CECAM1","CLE5CA","GGH","OLFM4","PTX3","TSPAN14"))
Tertiary_granules = list(c("ADAM8","ARSB","ASAH1","ATP11A","CAMP","CD58","CECR1","CR1","CTSB","CTSD","CTSH","CTSS","FCER1G","FCN1",
                           "ITGAX","LAMP2","LILRB2","MGAM","MMP8","MMP9","ORM1","PLAUR","PRCP","SIGLEC5","SIRPB1","SLC11A1","TIMP2","TNFAIP6",
                           "TNFSF14","VNN1"))

Ns <- AddModuleScore(Ns,features = Primary_granules,name="Primary_granules",search = TRUE);summary(Ns@meta.data$Primary_granules1)
v1 = VlnPlot(Ns,features = "Primary_granules1",pt.size = 0) + ylim(-0.5,1.9) ;v1
fix.sc <- scale_color_gradientn( colours = colors,  limits = c(-0.5,1.9))
G_AG = FeaturePlot(Ns,features = "Primary_granules1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
pG_AG <- lapply(G_AG, function (x) x + fix.sc);pG_AG

Ns <- AddModuleScore(Ns,features = Secondary_granules,name="Secondary_granules",search = TRUE);summary(Ns@meta.data$Secondary_granules1)
v2 = VlnPlot(Ns,features = "Secondary_granules1",pt.size = 0) + ylim(-0.5,1.9);v2
fix.sc <- scale_color_gradientn( colours =  colors,  limits = c(-0.5,1.9))
G_SG = FeaturePlot(Ns,features = "Secondary_granules1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
pG_SG <- lapply(G_SG, function (x) x + fix.sc);pG_SG

Ns <- AddModuleScore(Ns,features = Tertiary_granules,name="Tertiary_granules",search = TRUE);summary(Ns@meta.data$Tertiary_granules1)
v3 = VlnPlot(Ns,features = "Tertiary_granules1",pt.size = 0) + ylim(-0.5,1.9);v3
fix.sc <- scale_color_gradientn( colours = colors,  limits = c(-0.5,1.9))
G_GG = FeaturePlot(Ns,features = "Tertiary_granules1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
pG_GG <- lapply(G_GG, function (x) x + fix.sc);pG_GG
# 
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/FigS7_Evrad_UMAP_violin_original UMAP.eps", width = 30, height = 20)
multiplot(pG_AG,v1,pG_SG,v2,pG_GG,v3,cols = 3)
dev.off()
# 
library(clipr)
library(biomaRt)
# ####### Table S6. /home/pg/Documents/Shared-Win10/PMN_precursors/Grassi/Granule Protein Assignment, Related to Figure 6. .xlsx
grassi_EnsemblID = read_clip_tbl(x=read_clip())
ensembl <- useMart("ensembl", host = "www.ensembl.org",dataset = "hsapiens_gene_ensembl") #### version 80
attributes = listAttributes(ensembl)###
Granule_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name','description'),mart = ensembl,filters = "ensembl_gene_id",values=grassi_EnsemblID$ENSEMBL_ID)
Granule_genes$description = gsub('\\[.*','',Granule_genes$description)
data_grassi = merge(x=Granule_genes,y=grassi_EnsemblID,by.x="ensembl_gene_id",by.y="ENSEMBL_ID",all=FALSE)

AG = list(subset(data_grassi,data_grassi$Peak_in_Fraction=="AG")$external_gene_name)
GG = list(subset(data_grassi,data_grassi$Peak_in_Fraction=="GG")$external_gene_name)
CM = list(subset(data_grassi,data_grassi$Peak_in_Fraction=="CM")$external_gene_name)
SG = list(subset(data_grassi,data_grassi$Peak_in_Fraction=="SG")$external_gene_name)
FG = list(subset(data_grassi,data_grassi$Peak_in_Fraction=="FG")$external_gene_name)
SV = list(subset(data_grassi,data_grassi$Peak_in_Fraction=="SV")$external_gene_name)

Ns <- AddModuleScore(Ns,features = AG,name="AG",search = TRUE);summary(Ns@meta.data$AG1)
v1 = VlnPlot(Ns,features = "AG1",pt.size = 0) + ylim(-0.5,1.5);v1
fix.sc <- scale_color_gradientn( colours = colors,  limits = c(-0.5,1.5))
G_AG = FeaturePlot(Ns,features = "AG1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
pG_AG <- lapply(G_AG, function (x) x + fix.sc);pG_AG

Ns <- AddModuleScore(Ns,features = GG,name="GG",search = TRUE);summary(Ns@meta.data$GG1)
v2 = VlnPlot(Ns,features = "GG1",pt.size = 0)+ ylim(-0.5,1.5);v2
fix.sc <- scale_color_gradientn( colours = colors,  limits = c(0,1.5))
G_GG = FeaturePlot(Ns,features = "GG1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
pG_GG <- lapply(G_GG, function (x) x + fix.sc);pG_GG

Ns <- AddModuleScore(Ns,features = CM,name="CM",search = TRUE);summary(Ns@meta.data$CM1)
v3 = VlnPlot(Ns,features = "CM1",pt.size = 0) + ylim(-0.5,1.5) ;v3
fix.sc <- scale_color_gradientn( colours = colors,  limits = c(-0.5,1.5))
G_CM = FeaturePlot(Ns,features = "CM1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
pG_CM <- lapply(G_CM, function (x) x + fix.sc);pG_CM

multiplot(pG_AG,v1,pG_GG,v2,pG_CM,v3,cols = 3)

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/FigS7_Grassi_UMAP_violin_original UMAP.eps", width = 30, height = 20)
multiplot(pG_AG,v1,pG_SG,v2,pG_GG,v3,cols = 3)
dev.off()
 
cc.genes <- readLines(con = "/media/pg/Volume2_4T/RNA_SEQ_our/SingleCells/10XGenomics/regev_lab_cell_cycle_genes.txt")
s.genes <- list(cc.genes[1:43])
g2m.genes <- list(cc.genes[44:97])

Ns <- AddModuleScore(Ns,features = s.genes,name="s.genes",search = TRUE);summary(Ns@meta.data$s.genes1)
v1 = VlnPlot(Ns,features = "s.genes1",pt.size = 0) + ylim(-0.4,0.6) ;v1
fix.sc <- scale_color_gradientn( colours = colors,  limits = c(-0.5,0.6))
p.s.genes = FeaturePlot(Ns,features = "s.genes1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
p.s.genes.p <- lapply(p.s.genes, function (x) x + fix.sc);p.s.genes.p

Ns <- AddModuleScore(Ns,features = g2m.genes,name="g2m.genes",search = TRUE);summary(Ns@meta.data$g2m.genes1)
v2 = VlnPlot(Ns,features = "g2m.genes1",pt.size = 0) + ylim(-0.3,1.1) ;v2
fix.sc <- scale_color_gradientn( colours = colors,  limits = c(-0.3,1.1))
p.g2m.genes = FeaturePlot(Ns,features = "g2m.genes1",pt.size=0.1,combine = F,reduction = "umap_learn50PC")
p.g2m.genes.p <- lapply(p.g2m.genes, function (x) x + fix.sc);p.g2m.genes.p
# 
multiplot(p.s.genes.p,v1,p.g2m.genes.p,v2,cols = 2)
# 
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/FigS7_CC_UMAP_violin_original UMAP.eps", width = 30, height = 20)
multiplot(p.s.genes.p,v1,p.g2m.genes.p,v2,cols = 2)
dev.off()

####### Find Markers
DefaultAssay(Ns) <- "RNA"
Ns<- NormalizeData(Ns)
markers <- FindAllMarkers(object = Ns, only.pos = TRUE, min.pct = 0.25, min.diff.pct = -Inf, logfc.threshold = 0.25)
markers <- subset(markers,markers$p_val_adj < 0.01)
dim(subset(markers,markers$cluster=="3"))
data.frame(markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC))
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

DotPlot(Ns,features = rev(unique(top10$gene)),dot.scale = 4)+
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + coord_flip() +
  theme(axis.text.y = element_text(color = "grey20", size = 8, angle = 0, hjust = 1, vjust = 0, face = "plain"))

#######################################################################################################################
library(clusterProfiler)
universe = rownames(Ns)

{
onto = "BP"
cl0_genes = subset(markers,cluster=="0")$gene
cl0_GO <- enrichGO(cl0_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
head(cl0_GO)
###
# cl0_GO = clusterProfiler::simplify(cl0_GO)
# dotplot(cl0_GO, showCategory=30) + ggtitle("up ")
# TopGO = data.frame(cl0_GO)
# tableGO = cbind(TopGO,cluster="Clust_0")

tmp <- rep(1, ncol(cl0_GO) )# Batman!
cl0_GO <- data.frame(cl0_GO)
df <- rbind(cl0_GO, tmp )
colnames(df) <- colnames(cl0_GO)
tableGO = cbind(df,cluster="Clust_0")

cl1_genes = subset(markers,cluster=="1")$gene
cl1_GO <- enrichGO(cl1_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
head(cl1_GO)
cl1_GO = clusterProfiler::simplify(cl1_GO)
dotplot(cl1_GO, showCategory=30) + ggtitle("Cluster 1")
dataGO = cbind(data.frame(cl1_GO),cluster="Clust_1")
tableGO = rbind(tableGO,dataGO);tableGO

cl3_genes = subset(markers,cluster=="3")$gene
cl3_GO <- enrichGO(cl3_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
head(cl3_GO)
cl3_GO = clusterProfiler::simplify(cl3_GO)
dotplot(cl3_GO, showCategory=30) + ggtitle("Cluster 3")
dataGO = cbind(data.frame(cl3_GO),cluster="Clust_3")
tableGO = rbind(tableGO,dataGO);tableGO

cl4_genes = subset(markers,cluster=="4")$gene
cl4_GO <- enrichGO(cl4_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
head(cl4_GO)
cl4_GO = clusterProfiler::simplify(cl4_GO)
dotplot(cl4_GO, showCategory=30) + ggtitle("Cluster 4")
dataGO = head(cbind(data.frame(cl4_GO),cluster="Clust_4"),10)
tableGO = rbind(tableGO,dataGO);tableGO

tableGO2 = tableGO[,c("ID","Description","GeneRatio","p.adjust","cluster","geneID")]
tableGO2$nGene = as.numeric(gsub('/.*','',tableGO2$GeneRatio))
tableGO2$TGene = as.numeric(gsub('.*/','',tableGO2$GeneRatio))
tableGO2$GeneRatio2 = round(tableGO2$nGene/tableGO2$TGene,2)
tableGO2$logPvalue = -log10(tableGO2$p.adjust)

tableGO2$cluster2 = as.numeric(gsub('Clust_','',tableGO2$cluster))+1
tableGO2$cluster_order <- tableGO2$cluster2
tableGO2$cluster_order <- gsub('1','1',tableGO2$cluster_order)
tableGO2$cluster_order <- gsub('2','3',tableGO2$cluster_order)
tableGO2$cluster_order <- gsub('4','2',tableGO2$cluster_order)
tableGO2$cluster_order <- gsub('5','4',tableGO2$cluster_order)
tableGO2$cluster_order <- as.numeric(tableGO2$cluster_order)
}

BP = tableGO2
BP$Process = c("BP")

{
  onto = "MF"
  cl0_genes = subset(markers,cluster=="0")$gene
  cl0_GO <- enrichGO(cl0_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl0_GO)
  
  cl0_GO = clusterProfiler::simplify(cl0_GO)
  dotplot(cl0_GO, showCategory=30) + ggtitle("up ")
  TopGO = data.frame(cl0_GO)
  tableGO = cbind(TopGO,cluster="Clust_0")
  
  cl1_genes = subset(markers,cluster=="1")$gene
  cl1_GO <- enrichGO(cl1_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl1_GO)
  cl1_GO = clusterProfiler::simplify(cl1_GO)
  dotplot(cl1_GO, showCategory=30) + ggtitle("Cluster 1")
  dataGO = cbind(data.frame(cl1_GO),cluster="Clust_1")
  tableGO = rbind(tableGO,dataGO);tableGO
  
  cl3_genes = subset(markers,cluster=="3")$gene
  cl3_GO <- enrichGO(cl3_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl3_GO)
  cl3_GO = clusterProfiler::simplify(cl3_GO)
  dotplot(cl3_GO, showCategory=30) + ggtitle("Cluster 3")
  dataGO = cbind(data.frame(cl3_GO),cluster="Clust_3")
  tableGO = rbind(tableGO,dataGO);tableGO
  
  cl4_genes = subset(markers,cluster=="4")$gene
  cl4_GO <- enrichGO(cl4_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl4_GO)
  cl4_GO = clusterProfiler::simplify(cl4_GO)
  dotplot(cl4_GO, showCategory=30) + ggtitle("Cluster 4")
  dataGO = head(cbind(data.frame(cl4_GO),cluster="Clust_4"),10)
  tableGO = rbind(tableGO,dataGO);tableGO
  
  tableGO2 = tableGO[,c("ID","Description","GeneRatio","p.adjust","cluster","geneID")]
  tableGO2$nGene = as.numeric(gsub('/.*','',tableGO2$GeneRatio))
  tableGO2$TGene = as.numeric(gsub('.*/','',tableGO2$GeneRatio))
  tableGO2$GeneRatio2 = round(tableGO2$nGene/tableGO2$TGene,2)
  tableGO2$logPvalue = -log10(tableGO2$p.adjust)
  
  tableGO2$cluster2 = as.numeric(gsub('Clust_','',tableGO2$cluster))+1
  tableGO2$cluster_order <- tableGO2$cluster2
  tableGO2$cluster_order <- gsub('1','1',tableGO2$cluster_order)
  tableGO2$cluster_order <- gsub('2','3',tableGO2$cluster_order)
  tableGO2$cluster_order <- gsub('4','2',tableGO2$cluster_order)
  tableGO2$cluster_order <- gsub('5','4',tableGO2$cluster_order)
  tableGO2$cluster_order <- as.numeric(tableGO2$cluster_order)
}

MF = tableGO2
MF$Process = c("MF")

{
  onto = "CC"
  cl0_genes = subset(markers,cluster=="0")$gene
  cl0_GO <- enrichGO(cl0_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl0_GO)
  
  cl0_GO = clusterProfiler::simplify(cl0_GO)
  dotplot(cl0_GO, showCategory=30) + ggtitle("up ")
  TopGO = data.frame(cl0_GO)
  tableGO = cbind(TopGO,cluster="Clust_0")
  
  cl1_genes = subset(markers,cluster=="1")$gene
  cl1_GO <- enrichGO(cl1_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl1_GO)
  cl1_GO = clusterProfiler::simplify(cl1_GO)
  dotplot(cl1_GO, showCategory=30) + ggtitle("Cluster 1")
  dataGO = cbind(data.frame(cl1_GO),cluster="Clust_1")
  tableGO = rbind(tableGO,dataGO);tableGO
  
  cl3_genes = subset(markers,cluster=="3")$gene
  cl3_GO <- enrichGO(cl3_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl3_GO)
  cl3_GO = clusterProfiler::simplify(cl3_GO)
  dotplot(cl3_GO, showCategory=30) + ggtitle("Cluster 3")
  dataGO = cbind(data.frame(cl3_GO),cluster="Clust_3")
  tableGO = rbind(tableGO,dataGO);tableGO
  
  cl4_genes = subset(markers,cluster=="4")$gene
  cl4_GO <- enrichGO(cl4_genes, 'org.Hs.eg.db', ont=onto, pvalueCutoff=0.01,keyType='SYMBOL',universe=universe)
  head(cl4_GO)
  cl4_GO = clusterProfiler::simplify(cl4_GO)
  dotplot(cl4_GO, showCategory=30) + ggtitle("Cluster 4")
  dataGO = head(cbind(data.frame(cl4_GO),cluster="Clust_4"),10)
  tableGO = rbind(tableGO,dataGO);tableGO
  
  tableGO2 = tableGO[,c("ID","Description","GeneRatio","p.adjust","cluster","geneID")]
  tableGO2$nGene = as.numeric(gsub('/.*','',tableGO2$GeneRatio))
  tableGO2$TGene = as.numeric(gsub('.*/','',tableGO2$GeneRatio))
  tableGO2$GeneRatio2 = round(tableGO2$nGene/tableGO2$TGene,2)
  tableGO2$logPvalue = -log10(tableGO2$p.adjust)
  
  tableGO2$cluster2 = as.numeric(gsub('Clust_','',tableGO2$cluster))+1
  tableGO2$cluster_order <- tableGO2$cluster2
  tableGO2$cluster_order <- gsub('1','1',tableGO2$cluster_order)
  tableGO2$cluster_order <- gsub('2','3',tableGO2$cluster_order)
  tableGO2$cluster_order <- gsub('4','2',tableGO2$cluster_order)
  tableGO2$cluster_order <- gsub('5','4',tableGO2$cluster_order)
  tableGO2$cluster_order <- as.numeric(tableGO2$cluster_order)
}

CC = tableGO2
CC$Process = c("CC")

data_GO2 = rbind(BP,MF,CC)
WriteXLS(data_GO2,ExcelFileName = "/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/GO_Ns.xls")

data_GO2 <-  xlsx::read.xlsx(file = "/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/GO_Ns.xls",1)
data_GO2 = data_GO2[order(data_GO2$cluster_order,decreasing = T),]

max(data_GO2$logPvalue)
top5GO = data_GO2 %>% group_by(cluster2) %>% top_n(n = 5, wt = logPvalue)
top10GO = data_GO2 %>% group_by(cluster2) %>% top_n(n = 10, wt = logPvalue)
table(data_GO2$p.adjust)

data <- as.data.frame(top10GO)

WriteXLS(top5GO,ExcelFileName = "/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/GO_Top5_Ns.xls")


library(plyr)
p <- ggplot(top10GO, # you can replace the numbers to the row number of pathway of your interest
            aes(x = reorder(cluster,cluster_order), y = reorder(Description,desc(cluster_order)))) + geom_point(aes(size = GeneRatio2, color = logPvalue)) +theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 43), low="blue",high = "red") + ylab(NULL) +
  ggtitle("Gene Ontology enrichment")
p

##### Figure S7 E
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/top10_all_GO_Ns.eps", width = 10, height = 20)
p
dev.off()

####################################################
######## Re-intetration of NMs population     ######
source("/media/pg/Volume2_4T/RNA_SEQ_our/SingleCells/SingleCellsUtils.R")
Ns = readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.Ns.clean_50PC_Mt_CC.rds")

DefaultAssay(Ns) <- "RNA"
Ns <- DietSeurat(Ns, assays = "RNA")
meta.data = Ns@meta.data[,1:17]
Ns@meta.data <- meta.data

options(future.globals.maxSize = 20000 * 1024^2)
variables = c("percent.mt","CC.Difference")
Ns.list <- SplitObject(Ns, split.by = "donor")
Ns.list <- Ns.list[c("1","2","3")]
for (i in 1:length(Ns.list)) {
  Ns.list[[i]] <- SCTransform(Ns.list[[i]], verbose = TRUE,vars.to.regress=variables)
}

Ns.features <- SelectIntegrationFeatures(object.list = Ns.list, nfeatures = 3000)
Ns.list <- PrepSCTIntegration(object.list = Ns.list, anchor.features = Ns.features, verbose = FALSE)

reference_dataset <- which(names(Ns.list) == "2")
Ns.anchors <- FindIntegrationAnchors(object.list = Ns.list, normalization.method = "SCT",
                                     anchor.features = Ns.features, reference = reference_dataset,dims = 1:30) ##
all_features = lapply(Ns.list,row.names) %>% Reduce(intersect, .)
Ns.integrated <- IntegrateData(anchorset = Ns.anchors, normalization.method = "SCT", verbose = FALSE,features.to.integrate=all_features, dims = 1:30)###
saveRDS(Ns.integrated, file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/precursors.BMs.integrated_contaminants_clean_Ns.var_Mt_CC_allGenes_50PC_Genes.rds")

###########################################################
##### Read Integrated DataSet
Ns = readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/precursors.BMs.integrated_contaminants_clean_Ns.var_Mt_CC_allGenes_50PC_Genes.rds")
head(Ns@meta.data)
Idents(Ns)

DefaultAssay(Ns) = "RNA" ####"integrated"
Idents(Ns)
Ns<- NormalizeData(Ns)
markers <- FindAllMarkers(object = Ns, only.pos = TRUE, min.pct = 0.25, min.diff.pct = -Inf, logfc.threshold = 0.25)
markers <- subset(markers,markers$p_val_adj < 0.01)
top10_Ns <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#### Scale variable genes for heatmap
Ns <- ScaleData(Ns,verbose = F,features=markers$gene,assay = "RNA",vars.to.regress = c("percent.mt","CC.Difference","donor"))

reorder_for_heatmap <- function(markers,old_cluster=old,new_cluster=new){
  markers = as.data.frame(markers)
  markers$cluster2 = markers$cluster
  markers$cluster2 <- qdap::mgsub(old_cluster,new_cluster,markers$cluster2,fixed = FALSE)
  markers = markers[order(markers$cluster2),]
  return(markers)
}

my_levels <- c("0","3","1","4");length(my_levels) ### 0.4
Ns@active.ident <- factor(x=Ns@active.ident,levels = my_levels)

old = as.character(c("0","1","3","4")) ## 3
new = as.character(c("a","c","b","d")) ## 3
top10_Ns = reorder_for_heatmap(top10_Ns)
markers = reorder_for_heatmap(markers)

reorder_for_heatmap <- function(markers,old_cluster=old,new_cluster=new){
  markers = as.data.frame(markers)
  markers$cluster2 = markers$cluster
  markers$cluster2 <- qdap::mgsub(old_cluster,new_cluster,markers$cluster2,fixed = FALSE)
  markers = markers[order(markers$cluster2),]
  return(markers)
}


My_seurat_heatmap <- function(seurat,sigMarkers,old_cluster=NULL,new_cluster,top_genes=NULL){
  sigMarkers <- sigMarkers[order(sigMarkers$gene, sigMarkers$p_val_adj, decreasing=TRUE),]
  sigMarkers = sigMarkers[!duplicated(sigMarkers$gene),]
  seurat.scale.data <- as.matrix(GetAssayData(object = seurat,slot="scale.data"))
  colnames(seurat.scale.data) <-Idents(seurat)
  seurat.scale.data = seurat.scale.data[,order(colnames(seurat.scale.data))]
  colnames(seurat.scale.data) <- qdap::mgsub(old_cluster,new_cluster,colnames(seurat.scale.data),fixed = FALSE)
  seurat.scale.data= seurat.scale.data[,order(colnames(seurat.scale.data))]
  seurat.scale.data.top = subset(seurat.scale.data,rownames(seurat.scale.data)%in%sigMarkers$gene)
  
  library(ComplexHeatmap)
  library(circlize)
  
  ##### Annotation
  sigMarkers$cluster = qdap::mgsub(old_cluster,new_cluster,sigMarkers$cluster,fixed = FALSE)
  sigMarkers = sigMarkers[order(sigMarkers$cluster),]
  # ##### Clusters
  #
  cluster = sigMarkers$cluster
  names(cluster) = unique(sigMarkers$gene)
  cluster <- cluster[!is.na(names(cluster))]
  seurat.scale.data.top = seurat.scale.data.top[order(match(rownames(seurat.scale.data.top), names(cluster))),]
  cluster_anno <-  colnames(seurat.scale.data)
  #names(split_row) <- NULL
  top <- sigMarkers %>% group_by(cluster) %>% top_n(top_genes, avg_logFC)
  at = which(rownames(seurat.scale.data.top)%in%top$gene)
  anno = rowAnnotation(link= anno_mark(at=at,labels = rownames(seurat.scale.data.top)[at],
                                       labels_gp = gpar(fontsize=8),padding = 1,
                                       link_width=unit(1,"cm"),link_height=unit(1,"cm")))
  cluster = subset(cluster,names(cluster)%in%rownames(seurat.scale.data.top))
  split_row = factor(cluster)
  col_fun = circlize::colorRamp2(c(-3, 0, 3), c("#FF00FF", "black", "#FFFF00"))
  ht1 = Heatmap(seurat.scale.data.top, name = "Expression",
                column_split = factor(cluster_anno),
                row_split = split_row,
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                cluster_column_slices = TRUE,
                cluster_row_slices = TRUE,
                column_title_gp = gpar(fontsize = 8),
                column_gap = unit(0.5, "mm"),
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                col = col_fun,
                row_names_gp = gpar(fontsize = 4),
                column_title_rot = 90,
                top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                show_column_names = FALSE,
                use_raster = TRUE,
                raster_quality = 4)
  
  anno2 = Heatmap(cluster, name = "cluster", show_row_names = FALSE, width = unit(2, "mm"),row_split=NULL,row_gap = unit(0.5, "mm"))
  p <- anno2 + ht1 #+ anno
  p
  
  #split_row
  
}

DefaultAssay(Ns) <- "integrated"
p = My_seurat_heatmap(seurat = Ns,sigMarkers = markers,old_cluster = old,new_cluster = new,top_genes = 10 )
p

###### Figure 7E
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Heatmap_markers_3000_genes.eps", width = 10, height = 5)
p
dev.off()

#######################################################################################################
#### Pseudotime and trajectory analysis
######### Destiny
DefaultAssay(Ns)  <- "integrated"
Ns$Cells = Idents(Ns)
#Ns = FindVariableFeatures(Ns)
length(VariableFeatures(Ns))
top10 <- head(VariableFeatures(Ns), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Ns)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

Log_Ns_SC = GetAssayData(object = Ns, slot = "data")
Log_Ns_SC = Log_Ns_SC[VariableFeatures(Ns)[1:2000],]

VarGenes <- VariableFeatures(Ns)
VarGenes[grepl("IRF8",VarGenes)]

sc_norm = t(as.matrix(Log_Ns_SC))
cellLabels <- Idents(Ns)
colnames(Log_Ns_SC) <- cellLabels
dm <- destiny::DiffusionMap(sc_norm)
meta.data = Ns@meta.data

saveRDS(dm, file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/precursors.dm.BMs.integrated_contaminants_clean_Ns.var_Mt_CC_allGenes_50PC_Genes.rds")
dm = readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/precursors.dm.BMs.integrated_contaminants_clean_Ns.var_Mt_CC_allGenes_50PC_Genes.rds")

idx = which(rownames(meta.data)=="342260_1")

head(Ns@meta.data)
DimPlot(Ns,reduction = "umap_learn30PC")
dpt <- destiny::DPT(dm,tips=idx)
saveRDS(dpt, file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/precursors.dpt.BMs.integrated_contaminants_clean_Ns.var_Mt_CC_allGenes_50PC_Genes.rds")
dpt = readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/precursors.dpt.BMs.integrated_contaminants_clean_Ns.var_Mt_CC_allGenes_50PC_Genes.rds")

# Plot DC1 vs DC2 and color the cells by their inferred diffusion pseudotime.
# We can accesss diffusion pseudotime via dpt$dpt.

df <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2],
                 dptval = dpt$dpt, cell_type2 = Idents(Ns))
p1 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = dptval))
p2 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = cell_type2))
p <- cowplot::plot_grid(p1, p2)
p

###### from https://github.com/satijalab/seurat/issues/1475
tmp <- data.matrix(data.frame(DC1 = eigenvectors(dm)[, 1],
                              DC2 = eigenvectors(dm)[, 2]))
tmp2 = tmp[order(tmp$DC2),]

### load the original Ns seurat object
Ns = readRDS(file = "/home/pg/RNA_SEQ_our/SingleCells/Rhapsody/precursors/IntegrationSamples/RDS_files/BMs.integrated.Ns.clean_50PC_Mt_CC.rds")


rownames(tmp) <- colnames(Ns)
Ns[["dm"]] <- CreateDimReducObject(embeddings = tmp, key="DC_", assay=DefaultAssay(Ns))
DimPlot(Ns, reduction="dm")

Ns$pseudotime_dpt <- rank(dpt$dpt)
v <- VlnPlot(Ns,features = "pseudotime_dpt",pt.size = 0) +
  stat_summary(fun.y = median, geom='point', size = 25, colour = "black", shape = 95)#  + coord_flip()
v

my_levels <- c("0","3","4","1");length(my_levels) ### 0.4
Ns@active.ident <- factor(x=Ns@active.ident,levels = my_levels)

median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out)
}

v <- VlnPlot(Ns,features = "pseudotime_dpt",pt.size = 0)  +
  stat_summary(fun.y = median.stat, geom='point', size = 2, colour = "black") + coord_flip()
v

######  Figura 7C
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig6_violin_destiny_reintegrations_Ns_UMAP_3000_genes_30PC.eps", width = 8, height = 4)
v
dev.off()

Ns <- ProjectDim(Ns, reduction = "dm",dims.print=1:2,assay="integrated")

DimHeatmap(Ns, reduction = "dm", dims = 1:2, cells = 500, projected = TRUE, balanced = TRUE,assays="integrated")

v1 = VlnPlot(Ns, features = "DC_1",pt.size = 0)
v2 = VlnPlot(Ns, features = "DC_2",pt.size = 0)
multiplot(v1,v2,cols = 1)

# load("destiny_outpu_local_sigma.RData")
p <- plot(dpt)
p0 <- plot(dpt, col_by='ELANE', pal=viridis::magma)
p2 <- plot(dpt, col_by='CPA3', pal=viridis::magma)
p3 <- plot(dpt, col_by='TUBB', pal=viridis::magma)
p4 <- plot(dpt, col_by='IFI6', pal=viridis::magma)
p5 <- plot(dpt, col_by='BEX1', pal=viridis::magma)
p6 <- plot(dpt, col_by='CENPF', pal=viridis::magma)

gridExtra::grid.arrange(p,p6,p3,p2,p4,p5,p0,ncol=4)

col <- rgl::dataset_get_feature(dataset(dpt), 'ELANE')
exp = data.frame(DC1 = eigenvectors(dm)[,1],
                 DC2 = eigenvectors(dm)[,2])

DefaultAssay(Ns) <- "RNA"
Ns <- NormalizeData(Ns)
Ns$pseudotime_dpt <- rank(dpt$dpt)
p1 <- FeaturePlot(Ns,features = "pseudotime_dpt",reduction = "umap_learn30PC")

colors = c("#D9D9D9","#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#99000D","#8C2D04") ### grey and reds scale
exp = data.frame(FetchData(object = Ns, vars = c("CD34","ELANE","CTSG","ident"),slot = "data"),
                 umap = Embeddings(object = Ns, reduction = "umap_learn50PC"))
colnames(exp) <- c("CD34","ELANE","CTSG","ident","UMAP1","UMAP2")

library(ggthemes)
brewer.pal(n = 8, name = "Blues")

#### Figure 7D
colors = c("#D9D9D9","#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#99000D","#8C2D04")##
data <- subset(exp,exp$CD34>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
data$CD34[data$CD34 >= 2 ] <- 2
e1 <- ggplot(data2, aes(x=UMAP1, y=UMAP2)) + geom_point(colour = "#D9D9D9", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=UMAP1, y=UMAP2, color = CD34),size = 0.1) + scale_colour_gradientn(colours=colors) +
  xlab("UMAP1") +
  ylab("UMAP2") + theme_classic();e1

#### Figure 7D
colors = c("#D9D9D9","#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#99000D","#8C2D04")##
data <- subset(exp,exp$CD34>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
data$CD34[data$CD34 >= 1.5 ] <- 1.5
e2 <- ggplot(data2, aes(x=UMAP1, y=UMAP2)) + geom_point(colour = "#D9D9D9", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=UMAP1, y=UMAP2, color = CD34),size = 0.1) + scale_colour_gradientn(colours=colors) +
  xlab("UMAP1") +
  ylab("UMAP2") + theme_classic();e2

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/UMAP_CD34_reintegrations_Ns_UMAP_3000_genes.eps", width = 20, height = 10)
multiplot(e1,e2,cols=2 )
dev.off()

data <- subset(exp,exp$ELANE>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
e1 <- ggplot(data2, aes(x=UMAP1, y=UMAP2)) + geom_point(colour = "#D9D9D9", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=UMAP1, y=UMAP2, color = ELANE),size = 0.1) + scale_colour_gradientn(colours=colors) +
  xlab("UMAP1") +
  ylab("UMAP2") + theme_classic();e1

data <- subset(exp,exp$CTSG>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
e2 <- ggplot(data2, aes(x=UMAP1, y=UMAP2)) + geom_point(colour = "#D9D9D9", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=UMAP1, y=UMAP2, color = CTSG),size = 0.1) + scale_colour_gradientn(colours=colors) +
  xlab("UMAP1") +
  ylab("UMAP2") + theme_classic();e2

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/UMAP_ELANE_CTSG_reintegrations_Ns_UMAP_3000_genes.eps", width = 20, height = 10)
multiplot(e1,e2,cols=2 )
dev.off()

Ns$pseudotime_dpt <- rank(dpt$dpt)
p2 <-FeaturePlot(Ns,features = "pseudotime_dpt",reduction = "umap_learn50PC",pt.size=0.1)

### Figura 7B
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/UMAP_pseudotime.eps", width = 8, height = 4)
multiplot(p2, cols=1)
dev.off()

##### function from https://github.com/satijalab/seurat/issues/2475
# draw median for a single plot
MedianVioSeu <- function() stat_summary(fun = median, geom='point', size = 1, colour = "black")

# for a list of ggplots
apply.MedianVioSeu.and.combine <- function(ls.ggplots=pl$A) {
  for (i in 1:length(ls.ggplots)) ls.ggplots[[i]] = ls.ggplots[[i]]+MedianVioSeu()
  return(CombinePlots(ls.ggplots, legend = 'none'))
}

clust_4 = subset(markers,markers$cluster=="0")
clust_4.ss = subset(clust_4,clust_4$gene%in%unlist(s.genes))
clust_4.gm = subset(clust_4,clust_4$gene%in%unlist(g2m.genes))
clust_4[order(clust_4$avg_logFC),]
data.frame(markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC))
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

C <- c("TUBB","CPA3","IGLL1","TYMS","HIST1H4C")
pl=NULL
pl$C <- VlnPlot(object = Ns, features = C, group.by =  NULL, pt.size = 0, assay = "RNA", log = F, combine = F,ncol=3)
pl$C <- apply.MedianVioSeu.and.combine(pl$C)
pl$C

D <- c("MKI67","TOP2A","UBE2C","CENPF")
pl=NULL
pl$D <- VlnPlot(object = Ns, features = D, group.by =  NULL, pt.size = 0, assay = "RNA", log = F, combine = F)
pl$D <- apply.MedianVioSeu.and.combine(pl$D)
pl$D

E <- c("ELANE","AZU1","PRTN3","CTSG")
pl$E <- VlnPlot(object = Ns, features = E, group.by =  NULL, pt.size = 0, assay = "RNA", log = F, combine = F)
pl$E <- apply.MedianVioSeu.and.combine(pl$E)
pl$E

clust_3 = subset(markers,markers$cluster=="3")

G <- c("BEX1","IFI6","ISG15","IFIT3")
pl$G <- VlnPlot(object = Ns, features = G, group.by =  NULL, pt.size = 0, assay = "RNA", log = F, combine = F)
pl$G <- apply.MedianVioSeu.and.combine(pl$G)
pl$G

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig7S_ViolinPlot_C_integrated_Ns.eps", width = 10, height = 10)
pl$C
dev.off()

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig7S_ViolinPlot_D_integrated_Ns.eps", width = 10, height = 10)
pl$D
dev.off()

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig7S_ViolinPlot_E_integrated_Ns.eps", width = 10, height = 10)
pl$E
dev.off()

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig7S_ViolinPlot_G_integrated_Ns.eps", width = 10, height = 10)
pl$G
dev.off()

library(RColorBrewer)
display.brewer.all()

cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Colors.eps", width = 20, height = 30)
display.brewer.all()
dev.off()

colors = c("#D9D9D9","#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#99000D","#8C2D04") ### grey and reds scale

DefaultAssay(Ns) = "RNA" ####"integrated"
Idents(Ns)
Ns <- NormalizeData(Ns)

exp = data.frame(FetchData(object = Ns, vars = c("TOP2A","IFI6","ELANE","BEX1","IFIT1",
                                                 "IFIT3","CTSG","MPO","CD34","ident"),slot = "data"),
                 DC1 = eigenvectors(dm)[,1],
                 DC2 = eigenvectors(dm)[,2])
library(ggthemes)

data <- subset(exp,exp$TOP2A>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
data$TOP2A[data$TOP2A >= 2 ] <- 2
e1 <- ggplot(data2, aes(x=DC1, y=DC2)) + geom_point(colour = "#D9D9D9", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=DC1, y=DC2, color = TOP2A),size = 0.1) + scale_colour_gradientn(colours=colors) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05);e1

data <- subset(exp,exp$IFI6 > 0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
data$IFI6[data$IFI6 >= 4 ] <- 4
e2 <- ggplot(data2 , aes(x=DC1, y=DC2)) + geom_point(colour = "#737373", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=DC1, y=DC2, color = IFI6),size = 0.1) + scale_colour_gradientn(colours=colors) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05);e2

data <- subset(exp,exp$ELANE >0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
e3 <- ggplot(data2, aes(x=DC1, y=DC2)) + geom_point(colour = "#737373", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=DC1, y=DC2, color = ELANE),size = 0.1) + scale_colour_gradientn(colours=colors)  +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05);e3

data <- subset(exp,exp$BEX1>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
max(data$BEX1)

e4 <- ggplot(data2, aes(x=DC1, y=DC2)) + geom_point(colour = "#D9D9D9", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=DC1, y=DC2, color = BEX1),size = 1) + scale_colour_gradientn(colours=colors)  +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05);e4

data <- subset(exp,exp$CTSG>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
max(data$CTSG)
#data$CTSG[data$CTSG>6.5] <- 6.5
e5 <- ggplot(data2, aes(x=DC1, y=DC2)) + geom_point(colour = "#737373", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=DC1, y=DC2, color = CTSG),size = 0.1) + scale_colour_gradientn(colours=colors)  +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05);e5

data <- subset(exp,exp$IFIT3>0)
data2 = subset(exp,!rownames(exp)%in%rownames(data))
e6 <- ggplot(data2, aes(x=DC1, y=DC2)) + geom_point(colour = "#737373", size = 0.1) + geom_rangeframe() + theme_tufte() +
  geom_point(data=data,aes(x=DC1, y=DC2, color = IFIT3),size = 0.1) + scale_colour_gradientn(colours=colors) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05);e6

#### Figura 7H
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig6_destiny_genes_reintegrations_Ns_UMAP_3000_genes_news.eps", width = 20, height = 30)
multiplot(e1,e2,e3,e4,e5,e6,cols=2 )
dev.off()

########
tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Clusters = Idents(Ns),
                  pseudotime_dpt= Ns$pseudotime_dpt)

library(ggplot2);library(ggthemes)
cluster_col <- c("0" = "#E41A1C",#
                 "1" = "#00008B",#
                 "2" = "#006400",#
                 "3" = "#FDB462",#
                 "4" = "#87CEEB",
                 "5" = "#7CFC00",#
                 "6" = "#000000",
                 "7" = "#F781BF",
                 "8" =  "#999999")

p1 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = Clusters)) +
  geom_point(size = 0.01) + scale_color_manual(values=cluster_col) + xlab("Diffusion component 1") +
  ylab("Diffusion component 2") + theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05)

p2 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = pseudotime_dpt)) +
  geom_point(size = 0.01) + scale_color_gradient(low = "lightgrey", high = "blue") +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic() + xlim(-0.05,0.05) + ylim(-0.05,0.05)

p2

#### Figura 7F e G
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig6_DimPlot_destiny_reintegrations_Ns_UMAP_3000_genes.eps", width = 20, height = 8)
multiplot(p1,p2, cols = 2)
dev.off()

my_destinyPlot <- function(df,label=NULL){
  library(ggthemes)
  library(ashr)
  data = subset(df,df$Clusters==label)
  data$density <- get_density(data$DC1, data$DC2, n = 100)
  p1 <- ggplot(df, aes(x=DC1, y=DC2)) + geom_point(colour = "gray", size = 0.1) + geom_rangeframe() + theme_tufte() +
    geom_point(data=data,aes(x=DC1, y=DC2, color = density),size = 0.1) + xlim(-0.05,0.05) + ylim(-0.05,0.05) + scale_color_viridis() + ggtitle(label)
  #return(data)
  return(p1)
}

my_destinyPlot(tmp,label = "0")
table(tmp$Clusters)

p1 <- my_destinyPlot(tmp,label = "0") # 6
p2 <- my_destinyPlot(tmp,label = "1") # 3
p3 <- my_destinyPlot(tmp,label = "3") # 4
p4 <- my_destinyPlot(tmp,label = "4") # 2

multiplot(p1,p2,p3,p4,cols = 2)

#### Figura S7
cairo_ps(file="/home/pg/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/Fig6_densityPlot_destiny_reintegrations_Ns_UMAP_3000_genes.eps", width = 10, height = 10)
multiplot(p1,p2,p3,p4,cols = 2)
dev.off()
