library(DESeq2)
ddsTxi <- readRDS(file = "/home/pg/Documents/Shared-Win10/github_repository/precursors/ddsTxi_allSamples.rds")
### Annotation of ensembl ID

ddsTxi = Coding_Filtering(ddsTxi,version="http://apr2019.archive.ensembl.org") #### 27674  genes, 48 samples
names = as.vector(colData(ddsTxi)$Cell)

### Excluding low expressed genes
ddsTxi_expressed = ExpressedGenes_by_fpkm_dds(ddsTxi,name=names,threshold = 1);ddsTxi_expressed #### 12520 genes
ddsTxi_expressed = DESeq(ddsTxi_expressed,parallel=T,test="LRT",full= ~ batch + Cell,reduce= ~ batch)

##### GEO
fpkm = as.data.frame(fpkm(ddsTxi_expressed))
colnames(fpkm) <- colData(ddsTxi_expressed)$new.name
fpkm$genes <- mcols(ddsTxi_expressed)$SYMBOL
fpkm$gene_name <- mcols(ddsTxi_expressed)$gene_name
View(as.data.frame(colData(ddsTxi_expressed)))

WriteXLS::WriteXLS(fpkm,ExcelFileName="/home/pg/Documents/Shared-Win10/GEO_Data_Submission/Neutrophils_precursors/smartseq/FPKM.xls",row.names="TRUE")

### Getting differential expressed genes
res = results(ddsTxi_expressed)
res$SYMBOL = mcols(ddsTxi_expressed)$SYMBOL
summary(res)

res_sig = subset(res,padj <=0.01) #### 9834 genes
res_sig =res_sig[order(res_sig$padj),]
select = which(rownames(ddsTxi_expressed)%in%rownames(res_sig));length(select) ####  10614
de <- rownames(res_sig)

######## PCA analysis
vsd = vst(ddsTxi, blind=T) ###### set blind TRUE
mat <- assay(vsd)
batch = vsd$RUN
mod <- model.matrix(~ Cell, colData(ddsTxi))
mat <- limma::removeBatchEffect(mat, batch=batch,design=mod)
assay(vsd) <- mat

### Figure 5A
cairo_ps(file="/home/pg/Documents/Shared-Win10/PMN_precursors/figures/PCA_PMNs.eps", width = 10, height = 10)
pcaplot(vsd[de,], intgroup = c("Cell"),ellipse=F,pcX = 1, pcY = 2,ntop = nrow(vsd[de,]),text_labels = F)
dev.off()
################################################################################################àà
##### Holo dendrogram
###################################################
HCL_OLO = function(dds,genes,ntop,exclude=NULL){
  dds <- dds[, !(colnames(dds) %in% exclude)]
  dds <- dds[genes,]
  dds$Cell = droplevels(dds$Cell)
  dds$RUN = droplevels(dds$RUN)
  dds$Donor = droplevels(dds$Donor)
  rld  <- vst(dds, blind=T)
  mod <- model.matrix(~ Cell, colData(dds))
  mat <- assay(rld)
  batch = rld$RUN
  mat <- limma::removeBatchEffect(mat, batch=batch,batch2 = NULL,design=mod)
  assay(rld) <- mat
  colnames(rld) = paste(colData(rld)$Cell,colData(rld)$order,"#",colData(rld)$Donor,sep=" ")
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
cairo_ps(file="/home/pg/Documents/Shared-Win10/PMN_precursors/figures/olo_dendrogram_Fig6C.eps",width = 14, height = 5)
plot(dend2,hang = -1)
dev.off()
#####
#### List from /home/pg/Documents/Shared-Win10/PMN_precursors/Grassi 'Table S6. Granule Protein Assignment, Related to Figure 6. .xlsx'
library(clipr)
grassi_EnsemblID = read_clip_tbl(x=read_clip())
############### Granule proteins
grassi_EnsemblID = grassi_EnsemblID[!duplicated(grassi_EnsemblID$ENSEMBL_ID),]
{ensembl <- useMart("ensembl", host = "http://may2015.archive.ensembl.org",dataset = "hsapiens_gene_ensembl") #### version 80
  attributes = listAttributes(ensembl)###
  Granule_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name','description'),mart = ensembl,filters = "ensembl_gene_id",values=grassi_EnsemblID$ENSEMBL_ID)
  Granule_genes$description = gsub('\\[.*','',Granule_genes$description)
  grassi_EnsemblID = subset(grassi_EnsemblID,grassi_EnsemblID$ENSEMBL_ID%in%Granule_genes$ensembl_gene_id)
  data_grassi <- data.frame(grassi_EnsemblID[order(grassi_EnsemblID$ENSEMBL_ID),],Granule_genes[order(Granule_genes$ensembl_gene_id),])}

Median_values <- function(dds,groups){
  fpkm = DESeq2::fpkm(dds)
  rownames(fpkm) = mcols(ddsTxi_expressed)$SYMBOL
  data <- fpkm[order(rownames(fpkm)), ] ### sort first
  data = data[!duplicated(rownames(fpkm)),] ### Keep lowest fdr
  data = as.data.frame(round(medianMatrix(data,groups),2))
  data = data[,c("HSC","N0","N1","N2","N3","PM","M","MM","BC","SC","PMN")]
  return(data)
}

groups = colData(ddsTxi_expressed)$Cell
fpkm = Median_values(ddsTxi_expressed,groups=groups)
fpkm$gene = rownames(fpkm)
data <-  reshape2::melt(fpkm, id=c("gene"))
colnames(data) <- c("geneID","cells","expr")
table(data$cells)
data$order =as.numeric(data$cells)
data$order = as.numeric(gsub('HSC','1',data$order))
data$order = as.numeric(gsub('N0','2',data$order))
data$order = as.numeric(gsub('N1','3',data$order))
data$order = as.numeric(gsub('N2','4',data$order))
data$order = as.numeric(gsub('N3','5',data$order))
data$order = as.numeric(gsub('PMN','11',data$order))
data$order = as.numeric(gsub('MM','8',data$order))
data$order = as.numeric(gsub('PM','9',data$order))
data$order = as.numeric(gsub('M','7',data$order))
data$order = as.numeric(gsub('BC','9',data$order))
data$order = as.numeric(gsub('SC','10',data$order))
data$log2Exp= log2(data$expr+1)
table(data$order)

GranulePlot <- function(granules,data){
  for(i in unique(granules$Peak_in_Fraction)){
    df <- subset(granules,granules$Peak_in_Fraction==i)
    print(dim(df))
    exp <- subset(data,data$geneID%in%df$external_gene_name)
    tm =theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p <- ggplot(exp, aes(x=reorder(cells,order), y=log2Exp)) +
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
setwd("/home/pg/Documents/Shared-Win10/PMN_precursors/Grassi")
GranulePlot(data_grassi,data)

##############################################################
#### Evrad GO analysis
#### list from /home/pg/Documents/Shared-Win10/PMN_precursors/Evrad
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
    
    p <- ggplot(exp, aes(x=reorder(cells,order), y=log2Exp)) +
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
  fpkm = fpkm[,c("HSC","N0","N1","N2","N3","PM","M","MM","BC","SC","PMN")]
  max <- apply(fpkm, 1, max)
  dataLog = as.matrix(log2(fpkm[, 1:11] + 1))
  #z <- t(scale(t(dataLog), scale=TRUE, center=TRUE)) #Z transformation
  z = as.matrix(normalization(dataLog)) #Z transformation
  df = data.frame(z,mean= rowMeans(z),sd = rowSds(z),max=max)
  df  = df [,c("HSC","N0","N1","N2","N3","PM","M","MM","BC","SC","PMN","mean","sd","max")]
  data = data.frame(df,fdr=res$padj,SYMBOL=mcols(dds)$SYMBOL,gene_name=mcols(dds)$gene_name,gene_biotype=mcols(dds)$gene_biotype,fpkm)
  data = subset(data,rownames(data)%in%genes)
  data <- data[order(data$SYMBOL, abs(data$fdr) ), ] ### sort first
  data = data[!duplicated(data$SYMBOL),] ### Keep lowest fdr
  return(data)
}

GetTopVariableGenes = function(dds,ntop){
  rld  <- vst(dds, blind=F)
  colnames(rld) = paste(colData(ddsTxi_expressed)$Cell,colData(ddsTxi_expressed)$order,sep=" ")
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

ddsTxi_expressed_sig = ddsTxi_expressed[rownames(ddsTxi_expressed)%in%rownames(res_sig),]
selectedGenes = GetTopVariableGenes(ddsTxi_expressed_sig,ntop=2800) #### 2789.8 genes
groups = colData(ddsTxi_expressed)$Cell
# geni_Interesanty <- c("MYC","FLT3","CD34","SOX4","KIT","HLA-DRB1","HLA-DRB5","IRF8","RPS3A","RPS6","RPL3","RPL7",
#                       "ATP5F1A","IDH2","SSBP1","TOMM6","CEBPE","GFI1","FUT4","LYZ","LTF","LCN2","CEACAM8","MPO","AZU1",
#                       "PRTN3","ELANE","CXCR1","CXCR2","FCGR3B","IFITM1","IRF1","IFIT2","ARG1","MMP9","MMP25",
#                       "CD177","OASL","ISG15","RSAD2","IFIT1","IFIT3","CSF3R","CXCL8","MX1","MME","CEBPB","CYBA","CYBB",
#                       "HLA-DRA","HLA-DPB1","HOXA9","TOP2A", "MKI67","NFKBIA", "NFKB2","NFKBIZ", "CXCL8","IL1B", "LTB")

geni_Interesanty <- c("MYC","FLT3","CD34","SOX4","KIT","HLA-DRB1","HLA-DRB5","IRF8","RPS3A","RPS6","RPL3",
                      "TOMM6","CEBPE","GFI1","LYZ","LTF","LCN2","CEACAM8","MPO","AZU1",
                      "ELANE","CXCR1","CXCR2","FCGR3B","IFITM1","ARG1","MMP9","MMP25",
                      "CD177","ISG15","IFIT1","CSF3R","CXCL8","MX1","MME","CEBPB","CYBA","CYBB",
                      "HLA-DRA","HLA-DPB1","HOXA9","TOP2A", "MKI67","NFKBIA", "NFKB2","CXCL8","IL1B",
                      "MRPL1","MTCH2")
sort(geni_Interesanty)

Res_Sig = data.frame(res_sig,order=1:nrow(res_sig))
Res_Sig[match(geni_Interesanty,Res_Sig$SYMBOL),] %>% arrange(order)

data = PlotMediansHeatmap(ddsTxi_expressed,genes=rownames(selectedGenes),groups=groups,res=res)

data[match(geni_Interesanty,data$SYMBOL),]

x = data[,1:11]
rownames(x) = data$SYMBOL

nboot = 100
library(factoextra)
res.kmDEG <- eclust(x, FUNcluster = "kmeans", nstart = 25,k.max = 10,nboot=nboot,seed=123,k=10)
#Save an object to a file
#saveRDS(res.kmDEG, file = "/home/pg/Documents/Shared-Win10/PMN_precursors/R_Object/res.kmDEG.rds")
# Restore the object
res.km = readRDS(file = "/home/pg/Documents/Shared-Win10/PMN_precursors/R_Object/res.kmDEG.rds")
res.km= res.km
res.km$nbclust ### 10 clusters
res.km$clust_plot
res.km$size

split = factor(res.km$cluster, levels=c(7,2,3,9,1,4,8,6,10,5))
data2 = data.frame(data,cluster=res.km$cluster)
WriteXLS::WriteXLS(data2,ExcelFileName = "/home/pg/Documents/Shared-Win10/PMN_precursors/TopDEG_2800genes_new.xls")
### /home/pg/Dropbox/Precursors/Analisi RNAseq/ClustersAnalysis
################################################################################################
library(ComplexHeatmap)
x = res.km$data
annotation = HeatmapAnnotation(df = data.frame(cells = colnames(x[,1:11])))
####### save clustering object 
km=10
mat = res.km$data
library(circlize)
f1 = colorRamp2(seq(-2, 2, length = 3), c("blue","white","red"))
colors = c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(6, "Set3"))
reorder.hmap <- Heatmap(mat, col= f1 ,name ="zscore",split=split,cluster_rows = F,cluster_columns = F, 
                        show_row_names = F,show_column_names = F,height = unit(180, "mm"),width = unit(40, "mm"))
cluster = Heatmap(split,col = structure(colors, names = as.character(1:km)), show_row_names = FALSE,show_column_names = F,
                  show_heatmap_legend = FALSE, name ="cluster", width = unit(2, "mm"),show_row_dend=F,split=split,height = unit(180, "mm"))
at = which(rownames(mat)%in%geni_Interesanty)
anno = rowAnnotation(link= anno_mark(at=at,labels = rownames(mat)[at],
                                     labels_gp = gpar(fontsize=6),padding = 1))
p <- cluster + reorder.hmap + anno
p
#### Figura 5C
cairo_ps(file="/home/pg/Documents/Shared-Win10/Dropbox/Precursors/CoralDrawFigures/Heatmap_2800_DEG_genes_fig6D_new.eps", width = 5, height = 35)
p
dev.off()
######################################################################################################
data = data2 
####### All ontology
cluterProfilerGO <-function(DESeq2Res,sigGenes){
  onts = c( "MF", "BP", "CC")
  tab = as.list(onts)
  names(tab) = onts
  for(i in 1:3){
    ego <- enrichGO(gene=sigGenes,'org.Hs.eg.db',universe= rownames(DESeq2Res),minGSSize = 10,
                    keyType = 'ENSEMBL',ont=onts[i],pAdjustMethod="BH",pvalueCutoff=0.05,
                    qvalueCutoff=0.2,readable=TRUE)
    tab[[i]] = ego
  }
  return(tab)
}

rm(tableGO)

cl1 = subset(data,data$cluster==1)
cl1_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl1)) #### expressed

MFs = clusterProfiler::simplify(cl1_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl1_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl1_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
tableGO = cbind(TopGO,cluster="c5")

cl2 = subset(data,data$cluster==2)
cl2_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl2)) #### expressed

MFs = clusterProfiler::simplify(cl2_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl2_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl2_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c2")
tableGO = rbind(tableGO,dataGO)

cl3 = subset(data,data$cluster==3)
cl3_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl3)) #### expressed

MFs = clusterProfiler::simplify(cl3_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl3_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl3_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c3")
tableGO = rbind(tableGO,dataGO)

cl4 = subset(data,data$cluster==4)
cl4_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl4)) #### expressed

MFs = clusterProfiler::simplify(cl4_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl4_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl4_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c6")
tableGO = rbind(tableGO,dataGO)

cl5 = subset(data,data$cluster==5)
cl5_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl5)) #### expressed

MFs = clusterProfiler::simplify(cl5_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl5_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl5_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c10")
tableGO = rbind(tableGO,dataGO)

cl6 = subset(data,data$cluster==6)
cl6_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl6)) #### expressed

MFs = clusterProfiler::simplify(cl6_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl6_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl6_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c8")
tableGO = rbind(tableGO,dataGO)

cl7 = subset(data,data$cluster==7)
cl7_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl7))

MFs = clusterProfiler::simplify(cl7_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl7_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl7_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c1")
tableGO = rbind(tableGO,dataGO)

cl8 = subset(data,data$cluster==8)
cl8_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl8)) #### 

MFs = clusterProfiler::simplify(cl8_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl8_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl8_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c7")
tableGO = rbind(tableGO,dataGO)

cl9 = subset(data,data$cluster==9)
cl9_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl9)) ####

MFs = clusterProfiler::simplify(cl9_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl9_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl9_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c4")
tableGO = rbind(tableGO,dataGO)

cl10 = subset(data,data$cluster==10)
cl10_GO = cluterProfilerGO(DESeq2Res= ddsTxi,sigGenes=rownames(cl10)) #### expressed

MFs = clusterProfiler::simplify(cl10_GO[[1]],cutoff = 0.5) ### MF
Plot_GO_ClusterProfiler(head(MFs,10))
dotplot(MFs, showCategory=30) + ggtitle("cluters 1")
data_MFs = as.data.frame(MFs)
data_MFs$onto = "MF"

BPs = clusterProfiler::simplify(cl10_GO[[2]],cutoff = 0.5) ### BPs
Plot_GO_ClusterProfiler(head(BPs,10))
dotplot(BPs, showCategory=30) + ggtitle("cluters 1")
data_BPs = as.data.frame(BPs)
data_BPs$onto = "BF"

CCs = clusterProfiler::simplify(cl10_GO[[3]],cutoff = 0.5) ### CC
Plot_GO_ClusterProfiler(head(CCs,10))
dotplot(CCs, showCategory=30) + ggtitle("cluters 1")
data_CCs = as.data.frame(CCs)
data_CCs$onto = "CC"

TopGO = rbind(data_MFs,data_BPs,data_CCs)
dataGO = cbind(TopGO,cluster="c9")
tableGO = rbind(tableGO,dataGO)

table(tableGO$cluster)
tableGO2 = tableGO[,c("ID","Description","GeneRatio","pvalue","p.adjust","cluster","geneID","onto")]
tableGO2$nGene = as.numeric(gsub('/.*','',tableGO2$GeneRatio))
tableGO2$TGene = as.numeric(gsub('.*/','',tableGO2$GeneRatio))
tableGO2$GeneRatio2 = round(tableGO2$nGene/tableGO2$TGene,2)
tableGO2$logPvalue = -log10(tableGO2$p.adjust)
tableGO2$cluster2 = as.numeric(gsub('c','',tableGO2$cluster))
tableGO2 = tableGO2[order(tableGO2$cluster2,decreasing = F),]

WriteXLS::WriteXLS(tableGO2,ExcelFileName = "/home/pg/Documents/Shared-Win10/PMN_precursors/TopDEG_2800_GO2_all_annotate_allGO.xls")
tableGO2 = xlsx::read.xlsx(file = "/home/pg/Documents/Shared-Win10/PMN_precursors/TopDEG_2800_GO2_all_annotate_allGO.xls",1)

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

top5GO <- filter_GO_Table(tableGO2,top=5)
table(top5GO$cluster)
max(top5GO$logPvalue)

top5GO$Description <- as.character(top5GO$Description)
#Then turn it back into a factor with the levels in the correct order
top5GO$Description  <- factor(top5GO$Description, levels=unique(top5GO$Description))

p2 <- ggplot(top5GO, # you can replace the numbers to the row number of pathway of your interest
             aes(x = reorder(cluster,cluster2), y = reorder(Description,desc(cluster2)))) + geom_point(aes(size = GeneRatio2, color = logPvalue)) +theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 94), low="blue",high = "red") + ylab(NULL) +
  ggtitle("Gene Ontology enrichment") 
p2

cairo_ps(filename="/home/pg/Documents/Shared-Win10/PMN_precursors/figures/precursors_10_kmeans_heatmap_GO.ps",width=12,height = 10,pointsize=12)
p2
dev.off()
