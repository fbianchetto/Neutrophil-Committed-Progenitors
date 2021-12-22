library(DESeq2)
library(pcaExplorer)
library(dplyr)
source("/home/simple/Documents/Shared-Win10/github_repository/precursors/Utilis.R")
ddsTxi <- readRDS(file = "/home/simple/Documents/Shared-Win10/github_repository/precursors_NCPs/ddsTxi_allSamples.rds")
SampleTable.csv = colData(ddsTxi)

### Annotation of ensembl ID

ddsTxi = Coding_Filtering(ddsTxi,version="http://apr2019.archive.ensembl.org") #### 27674  genes, 48 samples
names = as.vector(colData(ddsTxi)$Cells)

### Excluding low expressed genes
ddsTxi_expressed = ExpressedGenes_by_fpkm_dds(ddsTxi,name=names,threshold = 1);ddsTxi_expressed #### 12520 genes
ddsTxi_expressed = DESeq(ddsTxi_expressed,parallel=T,test="LRT",full= ~ batch + Cells,reduce= ~ batch)

### Getting differential expressed genes
res = results(ddsTxi_expressed,alpha = 0.01)
res$SYMBOL = mcols(ddsTxi_expressed)$SYMBOL
res$gene_biotype = mcols(ddsTxi_expressed)$gene_biotype
summary(res)

res_sig = subset(res,padj <=0.01) #### 9834 genes
res_sig =res_sig[order(res_sig$padj),]
select = which(rownames(ddsTxi_expressed)%in%rownames(res_sig));length(select) ####  10614
de <- rownames(res_sig)

######## PCA analysis
vsd = vst(ddsTxi, blind=T) ###### set blind TRUE
mat <- assay(vsd)
batch = vsd$batch
mod <- model.matrix(~ Cells, colData(ddsTxi))
mat <- limma::removeBatchEffect(mat, batch=batch,design=mod)
assay(vsd) <- mat

### Figure 5A
cairo_ps(file="/home/pg/Documents/Shared-Win10/PMN_precursors/figures/PCA_PMNs.eps", width = 10, height = 10)
pcaplot(vsd[de,], intgroup = c("Cells"),ellipse=F,pcX = 1, pcY = 2,ntop = nrow(vsd[de,]),text_labels = T)
dev.off()

colnames(mat) = colData(ddsTxi_expressed)$Cells
rv <- rowVars(mat)
ntop=nrow(mat)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
mat <- mat[select, ]
pca <- prcomp(t(mat),scale. =F)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
group <- factor(c(colnames(mat)))
names = factor(c(colnames(mat)))

d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,names = names)

pcX = c(1)
pcY = c(2)
d <- data.frame(PC1 = pca$x[, pcX], PC2 = pca$x[, pcY], group = group, names = names)
colnames(d)[1] <- paste0("PC", pcX)
colnames(d)[2] <- paste0("PC", pcY)

cells_colors = c("BC"="#228B22",
                 "CMP" = "#7F5AA2",
                 "PM" = "#E9967A",
                 "MY" = "#8A2BE2",
                 "HSC"="#EE83EE",
                 "MM"="#00008B",
                 "NCP1"="#FF8C00",
                 "NCP2"="#7FFF00",
                 "NCP3"="#FF1493",
                 "NCP4"="#00FFFF",
                 "PMN"="#8B0000",
                 "SN"="#6495ED")
summary(d$PC1)
summary(d$PC2)
g <- ggplot(data = d, aes_string(x = paste0("PC", pcX), y = paste0("PC", pcY), color = "group")) + geom_point(size = 3) + 
  xlab(paste0("PC", pcX, ": ", round(percentVar[pcX] * 100, digits = 2), "% variance")) + 
  ylab(paste0("PC", pcY, ": ", round(percentVar[pcY] * 100, digits = 2), "% variance")) + theme_bw() 

g <- g + geom_text(label=names,nudge_x = 0.25, nudge_y = 0.25, check_overlap = T)# + xlim(-200,130) + ylim(-95,85) + theme(legend.position = "none") + coord_fixed()
p2 <- g + geom_mark_ellipse()+ scale_color_manual(values = cells_colors,aesthetics = c("colour", "fill"));p2

#cairo_ps(file="/home/simple/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/PCA_PMNs_with_CMP.eps", width = 10, height = 10)
#p2
#dev.off()

################################################################################################àà
##### Holo dendrogram
###################################################
HCL_OLO = function(dds,genes,ntop,exclude=NULL){
  dds <- dds[, !(colnames(dds) %in% exclude)]
  dds <- dds[genes,]
  dds$Cells = droplevels(dds$Cells)
  dds$batch = droplevels(dds$batch)
  dds$donor = droplevels(dds$donor)
  rld  <- vst(dds, blind=T)
  mod <- model.matrix(~ Cells, colData(dds))
  mat <- assay(rld)
  batch = rld$batch
  mat <- limma::removeBatchEffect(mat, batch=batch,batch2 = NULL,design=mod)
  assay(rld) <- mat
  colnames(rld) = dds$SampleNames
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  mat <- assay(rld)[select, ]
  return(mat)
}

library(dendextend)
data1 = HCL_OLO(ddsTxi_expressed,rownames(res_sig),ntop=nrow(res_sig),exclude=NULL) ### Fig 
data1 = data1[,order(colnames(data1))]

d <- dist(t(data1),method = "euclidean")
hc <- hclust(d,method = "ward.D2")
dend <- as.dendrogram(hc)
plot(dend)
dend2 <- seriate_dendrogram(dend, d,method = c("OLO"))
#### Figure 5B
cairo_ps(file="/home/simple/Documents/Shared-Win10/PMN_precursors/figures/olo_dendrogram_Fig6C.eps",width = 14, height = 5)
plot(dend2,hang = -1)
dev.off()
#####
#### List from ~/Documents/Shared-Win10/PMN_precursors/Grassi/'Table S6. Granule Protein Assignment, Related to Figure 6.xlsx'
library(clipr)
grassi_EnsemblID = read_clip_tbl(x=read_clip())
############### Granule proteins
grassi_EnsemblID = grassi_EnsemblID[!duplicated(grassi_EnsemblID$ENSEMBL_ID),]
{ensembl <- useMart("ensembl", host = "http://may2015.archive.ensembl.org",dataset = "hsapiens_gene_ensembl") #### version 80
  attributes = listAttributes(ensembl)###
  Granule_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name','description'),mart = ensembl,filters = "ensembl_gene_id",values=grassi_EnsemblID$ENSEMBL_ID,useCache = FALSE)
  Granule_genes$description = gsub('\\[.*','',Granule_genes$description)
  grassi_EnsemblID = subset(grassi_EnsemblID,grassi_EnsemblID$ENSEMBL_ID%in%Granule_genes$ensembl_gene_id)
  data_grassi <- data.frame(grassi_EnsemblID[order(grassi_EnsemblID$ENSEMBL_ID),],Granule_genes[order(Granule_genes$ensembl_gene_id),])}

Median_values <- function(dds,groups){
  fpkm = DESeq2::fpkm(dds)
  rownames(fpkm) = mcols(dds)$SYMBOL
  data <- fpkm[order(rownames(fpkm)),,drop=TRUE] 
  data = data[!duplicated(rownames(data)),]
  data = as.data.frame(round(medianMatrix(data,groups),2))
  data = data[,c("HSC","CMP","NCP1","NCP2","NCP3","NCP4","PM","MY","MM","BC","SC","PMN")]
  return(data)
}

groups = colData(ddsTxi_expressed)$Cells
fpkm = Median_values(ddsTxi_expressed,groups=groups)

fpkm$gene = rownames(fpkm)
data <-  reshape2::melt(fpkm, id=c("gene"))
colnames(data) <- c("geneID","cells","expr")
table(data$cells)
data$log2Exp= log2(data$expr+1)
table(data$order)

GranulePlot <- function(granules,data){
  for(i in unique(granules$Peak_in_Fraction)){
    df <- subset(granules,granules$Peak_in_Fraction==i)
    print(dim(df))
    exp <- subset(data,data$geneID%in%df$external_gene_name)
    tm =theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p <- ggplot(exp, aes(x=cells, y=log2Exp)) +
      geom_boxplot(fill = "white", color = "black",notch = FALSE,outlier.shape=NA) +
      geom_smooth(method = "loess", se=TRUE, aes(group=1)) +
      theme_bw() + ggtitle(i) + ylim(0,15) + tm + theme(text = element_text(size=8)) +
      theme(axis.text.x = element_text(angle = 45))
    p
    svg(file=paste0("boxplot_with_HSC_", i, ".svg"), width = 2.5, height =2)
    print(p)
    dev.off()
  }
}

####  Figura 5D
setwd("/home/simple/Documents/Shared-Win10/PMN_precursors/Grassi")
GranulePlot(data_grassi,data)

##############################################################
#### Evrad GO analysis
#### list from /home/simple/Documents/Shared-Win10/PMN_precursors/Evrad
Evrad_GO = read_clip_tbl(x=read_clip())
geneSets = list("Phagocytosis"=toupper(Evrad_GO$Phagocytosis),
                "ROS biosynthetic process"=toupper(Evrad_GO$ROS.biosynthetic.process),
                "ROS biosynthetic process"=toupper(Evrad_GO$ROS.biosynthetic.process),
                "Chemotaxis"=toupper(Evrad_GO$Chemotaxis),
                "cell.cycle"=toupper(Evrad_GO$cell.cycle))
geneSets = lapply(geneSets, function(x) x[!is.na(x)])

GOPlot <- function(GO,data){
  for(i in 1:length(GO)){
    a = GO[i]
    b = unlist(a)
    names(b) <- NULL
    exp <- subset(data,data$geneID%in%b)
    tm =theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    p <- ggplot(exp, aes(x=cells, y=log2Exp)) +
      geom_boxplot(fill = "white", color = "black",notch = FALSE,outlier.shape=NA) +
      geom_smooth(method = "loess", se=TRUE, aes(group=1)) +
      theme_bw() + ggtitle(names(a)) + tm + theme(text = element_text(size=8)) +
      theme(axis.text.x = element_text(angle = 45))
    p
    svg(file=paste0("boxplot_with_HSC_", names(a), ".svg"), width = 2.5, height =2)
    print(p)
    dev.off()
  }
  p
}

#### Figura 5D e Figura S6 B,C,D
GOPlot(geneSets,data)
##################################################
normalization <-function(x){
  dimm=dim(x)
  for(i in 1:dimm[1]){
    x[i,]=(x[i,]-mean(x[i,]))/sd(x[i,]) ### zscore
  }
  return(x)
}
########## top 20% genes
PlotMediansHeatmap <- function(dds,genes,groups,res){
  fpkm = DESeq2::fpkm(dds)
  fpkm = as.data.frame(round(medianMatrix(fpkm,groups),2))
  fpkm = fpkm[,c("HSC","CMP","NCP1","NCP2","NCP3","NCP4","PM","MY","MM","BC","SC","PMN")]
  max <- apply(fpkm, 1, max)
  dataLog = as.matrix(log2(fpkm[, 1:12] + 1))
  #z <- t(scale(t(dataLog), scale=TRUE, center=TRUE)) #Z transformation
  z = as.matrix(normalization(dataLog)) #Z transformation
  df = data.frame(z,mean= rowMeans(z),sd = rowSds(z),max=max)
  df  = df [,c("HSC","CMP","NCP1","NCP2","NCP3","NCP4","PM","MY","MM","BC","SC","PMN","mean","sd","max")]
  data = data.frame(df,fdr=res$padj,SYMBOL=mcols(dds)$SYMBOL,gene_name=mcols(dds)$gene_name,gene_biotype=mcols(dds)$gene_biotype,fpkm)
  data = subset(data,rownames(data)%in%genes)
  data <- data[order(data$SYMBOL, abs(data$fdr) ), ] ### sort first
  data = data[!duplicated(data$SYMBOL),] ### Keep lowest fdr
  return(data)
}

GetTopVariableGenes = function(dds,ntop){
  rld  <- vst(dds, blind=F)
  colnames(rld) = colData(ddsTxi_expressed)$Cells
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  mat <- assay(rld)[select, ]
  return(mat)
}

normalization <-function(x){
  dimm=dim(x)
  for(i in 1:dimm[1]){
    x[i,]=(x[i,]-mean(x[i,]))/sd(x[i,]) ### zscore
  }
  return(x)
}

ProteinCoding = res_sig[res_sig$gene_biotype=="protein_coding",]
ProteinCoding <- ProteinCoding[order(ProteinCoding$SYMBOL, abs(ProteinCoding$padj) ), ] ### sort first
ProteinCoding = ProteinCoding[order(ProteinCoding$padj),] ### Keep lowest fdr
ProteinCoding = ProteinCoding[!duplicated(ProteinCoding$SYMBOL),] ### Keep lowest fdr
ddsTxi_expressed_sig = ddsTxi_expressed[rownames(ddsTxi_expressed)%in%rownames(ProteinCoding),] ###
selectedGenes = GetTopVariableGenes(ddsTxi_expressed_sig,ntop=1900) #### 2789.8 genes
groups = colData(ddsTxi_expressed)$Cells

Res_Sig = data.frame(res_sig,order=1:nrow(res_sig))
Res_Sig[match(geni_Interesanty,Res_Sig$SYMBOL),] %>% arrange(order)

data = PlotMediansHeatmap(ddsTxi_expressed,genes=rownames(selectedGenes),groups=groups,res=res)
data = data %>% drop_na()
data = droplevels(data)
data[match("CD34",data$SYMBOL),]
data[match("ELANE",data$SYMBOL),]

str(data)

x = data[,1:12]
rownames(x) = data$SYMBOL
# Determining the number od clusters by compute gap statistic
# see https://uc-r.github.io/kmeans_clustering#gap
# see https://stackoverflow.com/questions/36240695/connect-ggplot-boxplots-using-lines-and-multiple-factor
# see https://www.datanovia.com/en/blog/clustering-example-4-steps-you-should-know
#nboot >= 500 is recommended
nboot = 500
library(factoextra)
res.kmDEG <- eclust(x, FUNcluster = "kmeans", nstart = 25,k.max = 20,nboot=nboot,seed=123,k=10)
sort(table(res.kmDEG$cluster))
# res.km  = res.kmDEG
#Save an object to a file
#saveRDS(res.kmDEG, file = "/home/simple/Documents/Shared-Win10/PMN_precursors/R_Object/res.kmDEG_1900.rds")
# Restore the object
res.km = readRDS(file = "/home/simple/Documents/Shared-Win10/PMN_precursors/R_Object/res.kmDEG_1900.rds")
sort(table(res.km$cluster))
res.km$nbclust ### 15
res.km$clust_plot
res.km$size
fviz_silhouette(res.km)
# Silhouette width of observation
sil <- res.km$silinfo$widths[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
neg = sil[neg_sil_index, , drop = FALSE]

pheatmap::pheatmap(res.km$centers,cluster_cols = F)
d <- dist(res.km$centers,method = "euclidean")
hc <- hclust(d,method = "ward.D2")
dend <- as.dendrogram(hc)
plot(dend)
dend2 <-  dendextend::seriate_dendrogram(dend, d,method = c("OLO"))
plot(dend2)
split = factor(res.km$cluster, levels=c(10,4,2,7,1,6,9,8,3,5)) ## 10
table(res.km$cluster)
all(data$SYMBOL==rownames(res.km$data))
data2 = data.frame(cluster=res.km$cluster,genes=names(res.km$cluster),res.km$data)
WriteXLS::WriteXLS(data2,ExcelFileName = "/home/simple/Documents/Shared-Win10/PMN_precursors/TopDEG_1900genes_new.xls")

### /home/pg/Dropbox/Precursors/Analisi RNAseq/ClustersAnalysis
geni_Interesanty <- sort(c("MYC","CD34","SOX4","KIT","HLA-DRB1","HLA-DRB5","RPS3A","RPS6","RPL3","RPL7","HOXA9",
                           "ATP5F1A","IDH2","SSBP1","TOMM6","CEBPE","GFI1","FUT4","LYZ","LTF","LCN2","CEACAM8","MPO","AZU1",
                           "PRTN3","ELANE","CXCR1","CXCR2","FCGR3B","IFITM1","IRF1","IFIT2","ARG1","MMP9","MMP25","ICAM3","ANXA3",
                           "CD177","OASL","ISG15","RSAD2","IFIT1","IFIT3","CSF3R","CXCL8","MX1","MME","CEBPB","CYBA","CYBB",
                           "DEFA1","DEFA3"))
################################################################################################
library(ComplexHeatmap)
x = res.km$data
annotation = HeatmapAnnotation(df = data.frame(cells = colnames(x[,1:12])))
####### save clustering object 
km=10
mat = as.matrix(res.km$data)
mat[grepl("CD34",rownames(mat),)]
library(circlize)
f1 = colorRamp2(seq(-2, 2, length = 3), c("blue","white","red"))
colors = c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(6, "Set3"))
reorder.hmap <- Heatmap(mat, col= f1 ,name ="zscore",split=split,cluster_rows = T,cluster_columns = F, 
                        show_row_names = F,show_column_names = F,height = unit(180, "mm"),width = unit(50, "mm"))
cluster = Heatmap(split,col = structure(colors, names = as.character(1:km)), show_row_names = FALSE,cluster_rows = T,,show_column_names = F,
                  show_heatmap_legend = FALSE, name ="cluster", width = unit(2, "mm"),show_row_dend=F,split=split,height = unit(180,"mm"))
at = which(rownames(mat)%in%geni_Interesanty)
anno = rowAnnotation(link= anno_mark(at=at,labels = rownames(mat)[at],
                                     labels_gp = gpar(fontsize=6),padding = 1))

p <- cluster + reorder.hmap + anno
p

cairo_ps(file="/home/simple/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/heatmap_PMNs_with_CMP.eps", width = 5, height = 10)
p
dev.off()

library(clipr)
grassi_EnsemblID = read_clip_tbl(x=read_clip()) ##### /home/simple/Documents/Shared-Win10/PMN_precursors/Grassi

grassi_EnsemblID = grassi_EnsemblID[!duplicated(grassi_EnsemblID$ENSEMBL_ID),]
{ensembl <- useMart("ensembl", host = "http://may2015.archive.ensembl.org",dataset = "hsapiens_gene_ensembl") #### version 80
  attributes = listAttributes(ensembl)###
  Granule_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name','description'),mart = ensembl,filters = "ensembl_gene_id",values=grassi_EnsemblID$ENSEMBL_ID,useCache = FALSE)
  Granule_genes$description = gsub('\\[.*','',Granule_genes$description)
  grassi_EnsemblID = subset(grassi_EnsemblID,grassi_EnsemblID$ENSEMBL_ID%in%Granule_genes$ensembl_gene_id)
  data_grassi <- data.frame(grassi_EnsemblID[order(grassi_EnsemblID$ENSEMBL_ID),],Granule_genes[order(Granule_genes$ensembl_gene_id),])}

sort(data_grassi$external_gene_name)
at = which(rownames(mat)%in%data_grassi$external_gene_name)
anno = rowAnnotation(link= anno_mark(at=at,labels = rownames(mat)[at],
                                     labels_gp = gpar(fontsize=6),padding = 1))

p <- cluster + reorder.hmap + anno
p

cairo_ps(file="/home/simple/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/heatmap_PMNs_with_CMP_grassi.eps", width = 5, height = 10)
p
dev.off()

######################################################################################################
data = data2 
####### All ontology
cluterProfilerGO.symbol = function(dds,sigGenes){
  onts = c( "MF", "BP", "CC")
  tab = as.list(onts)
  names(tab) = onts
  for(i in 1:3){
    ego <- enrichGO(gene=sigGenes,'org.Hs.eg.db',universe= mcols(dds)$SYMBOL,minGSSize = 10,
                    keyType = 'SYMBOL',ont=onts[i],pAdjustMethod="BH",pvalueCutoff=0.05,
                    qvalueCutoff=0.2,readable=FALSE)
    tab[[i]] = ego
  }
  return(tab)
}


{
  cl1 = rownames(subset(data2,data2$cluster==10))
  cl1_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl1) #### expressed
  BP = cl1_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  tableGO = cbind(TopGO,cluster="c1")
  
  ########################################
  cl2 = rownames(subset(data2,data2$cluster==4))
  cl2_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl2) #### expressed
  BP = cl2_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c2")
  tableGO = rbind(tableGO,dataGO)
  
  ########################################
  cl3 = rownames(subset(data2,data2$cluster==2))
  cl3_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl3) #### expressed
  BP = cl3_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c3")
  tableGO = rbind(tableGO,dataGO)
  
  ########################################
  cl4 = rownames(subset(data2,data2$cluster==7))
  cl4_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl4) #### expressed
  BP = cl4_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c4")
  tableGO = rbind(tableGO,dataGO)
  
  ########################################
  cl5 = rownames(subset(data2,data2$cluster==1))
  cl5_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl5) #### expressed
  BP = cl5_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c5")
  tableGO = rbind(tableGO,dataGO)
  
  ###########################################################
  cl6 = rownames(subset(data2,data2$cluster==6))
  cl6_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl6) #### expressed
  BP = cl6_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c6")
  tableGO = rbind(tableGO,dataGO)
  
  ###########################################################
  cl7 = rownames(subset(data2,data2$cluster==9))
  cl7_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl7) #### expressed
  BP = cl7_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c7")
  tableGO = rbind(tableGO,dataGO)
  
  ###########################################################
  cl8 = rownames(subset(data2,data2$cluster==8))
  cl8_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl8) #### expressed
  BP = cl8_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c8")
  tableGO = rbind(tableGO,dataGO)
  
  ###########################################################
  cl9 = rownames(subset(data2,data2$cluster==3))
  cl9_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl9) #### expressed
  BP = cl9_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c9")
  tableGO = rbind(tableGO,dataGO)
  ###########################################################
  cl10 = rownames(subset(data2,data2$cluster==5))
  cl10_GO = cluterProfilerGO.symbol(dds = ddsTxi_expressed,sigGenes=cl10) #### expressed
  BP = cl10_GO[["BP"]]
  head(BP,10)
  BPs = clusterProfiler::simplify(BP,cutoff = 0.5)
  head(BPs,10)
  TopGO = data.frame(BPs)
  dataGO = cbind(TopGO,cluster="c10")
  tableGO = rbind(tableGO,dataGO)
}

table(tableGO$cluster)
tableGO2 = tableGO[,c("ID","Description","GeneRatio","pvalue","p.adjust","cluster","geneID")]
tableGO2$nGene = as.numeric(gsub('/.*','',tableGO2$GeneRatio))
tableGO2$TGene = as.numeric(gsub('.*/','',tableGO2$GeneRatio))
#tableGO2$GeneRatio2 = round(tableGO2$nGene/tableGO2$TGene,2)
tableGO2$GeneRatio2 = tableGO2$nGene/tableGO2$TGene
tableGO2$logPvalue = -log10(tableGO2$p.adjust)
tableGO2$cluster2 = as.numeric(gsub('c','',tableGO2$cluster))
tableGO2 = tableGO2[order(tableGO2$cluster2,decreasing = F),]

filter_GO_Table <- function(table,top=top){
  df <- NULL
  cluster = unique(table$cluster)
  for(i in cluster){
    data = subset(table,table$cluster==i)
    data = data[order(data$pvalue),]
    data = head(data,top)
    df <- rbind(df,data)
    
  }
  df
}

max(tableGO2$logPvalue)
top10GO <- filter_GO_Table(tableGO2,top=10)
table(top10GO$cluster)

p <- ggplot(top10GO, # you can replace the numbers to the row number of pathway of your interest
            aes(x = reorder(cluster,cluster2), y = reorder(Description,desc(cluster2)))) + geom_point(aes(size = GeneRatio2, color = logPvalue)) +theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 75), low="blue",high = "red") + ylab(NULL) +
  ggtitle("Gene Ontology enrichment") 

p 

top5GO <- filter_GO_Table(tableGO2,top=5)
table(top5GO$cluster)

top5GO$Description <- as.character(top5GO$Description)
#Then turn it back into a factor with the levels in the correct order
top5GO$Description  <- factor(top5GO$Description, levels=unique(top5GO$Description))

p2 <- ggplot(top5GO, # you can replace the numbers to the row number of pathway of your interest
             aes(x = reorder(cluster,cluster2), y = reorder(Description,desc(cluster2)))) + geom_point(aes(size = GeneRatio2, color = logPvalue)) +theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, max(top5GO$logPvalue)), low="blue",high = "red") + ylab(NULL) +
  ggtitle("Gene Ontology enrichment")
p2 + coord_fixed(ratio = 0.8)
max(top5GO$logPvalue)
cairo_ps(file="/home/simple/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/heatmap_PMNs_with_CMP.eps", width = 5, height = 10)

cairo_ps(filename="/home/simple/Dropbox/Precursors/CoralDrawFigures/NewCoreldraw/heatmap_PMNs_with_CMP_GO2.ps",width=10,height = 6.85*2,pointsize=10)
p2  + theme(axis.title = element_blank(),
            plot.title = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "none",
            axis.text.y = element_text(color = "black", size = 12 ,angle = 0, hjust = 1, vjust = 0, face = "plain"))
dev.off()

