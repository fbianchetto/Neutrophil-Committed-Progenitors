#############################################################################################################
##### 1 .Pipeline for reading SevenBriges Genomics output and integration of Single cells from 3 Donors #####
#############################################################################################################
#Load packages
source("/media/simple/Volume2_4T/RNA_SEQ_our/SingleCells/SingleCells_functions/SingleCellsUtils.R")
library(Seurat)
library(ggplot2)
library(scater)
library(SingleR)
library(RColorBrewer)

#### Creating a Seurat object from SevenBridges ouput  
CreateSeuratObjectRhapsody <- function(umi,sample=NULL,batch=NULL,donor=NULL,metadata=meta.data){
  seurat <- CreateSeuratObject(counts = t(umi), min.cells = 10, min.features = 200,project = sample,names.delim="-",meta.data=meta.data)
  seurat[["percent.ribosomal"]] <- PercentageFeatureSet(object = seurat, pattern = "^RPS|^RPL")
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat[["mitoRatio"]] <- seurat@meta.data$percent.mt / 100
  seurat[["log10GenesPerUMI"]] <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)
  seurat@meta.data$donor <- donor
  cc.genes <- readLines(con = "/media/simple/Volume2_4T/RNA_SEQ_our/SingleCells/10XGenomics/regev_lab_cell_cycle_genes.txt")
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]
  Idents(object = seurat) <- "Cells"
  seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
  seurat@meta.data$CC.Difference <- seurat@meta.data$S.Score - seurat@meta.data$G2M.Score
  Idents(object = seurat) <- seurat$Cells
  return(seurat)
}
##### donor1
setwd("/home/simple/Documents/Shared-Win10/GEO_Data_Submission/GSE_164687_new_submision/SingleCells")
donor1_NCP1  <- read.table("donor1_SampleTag03_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor1_NCP2  <- read.table("donor1_SampleTag04_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor1_NCP3  <- read.table("donor1_SampleTag05_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor1_NCP4  <- read.table("donor1_SampleTag06_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor1_cMOP <- read.table("donor1_SampleTag07_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")

donor1_NCP1$Cells = "NCP1"
donor1_NCP2$Cells = "NCP2"
donor1_NCP3$Cells = "NCP3"
donor1_NCP4$Cells = "NCP4"
donor1_cMOP$Cells = "cMOP"
donor1 <- rbind(donor1_NCP1,donor1_NCP2,donor1_NCP3,donor1_NCP4,donor1_cMOP)

meta.data = data.frame(barcorde=rownames(donor1),Cell_Index=donor1$Cell_Index,Cells = donor1$Cells)

rownames(meta.data) <- meta.data$Cell_Index
UMI = donor1[,3:ncol(donor1)-1]
rownames(UMI) <- donor1$Cell_Index
UMI <- as.sparse(UMI)

donor1_seurat <- CreateSeuratObjectRhapsody(umi=UMI,metadata=meta.data.donor1,donor=1,sample="26022020")
VlnPlot(donor1_seurat ,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size = 0)
saveRDS(donor1_seurat, file = "donor1_seurat.rds")

#### donor2
donor2_NCP1  <- read.table("donor2_SampleTag08_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor2_NCP2  <- read.table("donor2_SampleTag09_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor2_NCP3  <- read.table("donor2_SampleTag10_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor2_NCP4  <- read.table("donor2_SampleTag11_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor2_cMOP <- read.table("donor2_SampleTag12_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")

donor2_NCP1$Cells = "NCP1"
donor2_NCP2$Cells = "NCP2"
donor2_NCP3$Cells = "NCP3"
donor2_NCP4$Cells = "NCP4"
donor2_cMOP$Cells = "cMOP"
donor2 <- rbind(donor2_NCP1,donor2_NCP2,donor2_NCP3,donor2_NCP4,donor2_cMOP)

meta.data = data.frame(barcorde=rownames(donor2),Cell_Index=donor2$Cell_Index,Cells = donor2$Cells)
rownames(meta.data) <- meta.data$Cell_Index
UMI = donor2[,3:ncol(donor2)-1]
rownames(UMI) <- donor2$Cell_Index
UMI <- as.sparse(UMI)

donor2_seurat <- CreateSeuratObjectRhapsody(umi=UMI,metadata=meta.data,donor=2,sample="30032020")
saveRDS(donor2_seurat, file = "donor2_seurat.rds")

##### donor3
donor3_NCP1  <- read.table("donor3_SampleTag01_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor3_NCP2  <- read.table("donor3_SampleTag02_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor3_NCP3  <- read.table("donor3_SampleTag03_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor3_NCP4  <- read.table("donor3_SampleTag04_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")
donor3_cMOP <- read.table("donor3_SampleTag05_hs_RSEC_MolsPerCell_with_barcode.csv.gz", sep = ",", header = TRUE, row.names = 1,check.names=FALSE,comment.char = "#")

donor3_NCP1$Cells = "NCP1"
donor3_NCP2$Cells = "NCP2"
donor3_NCP3$Cells = "NCP3"
donor3_NCP4$Cells = "NCP4"
donor3_cMOP$Cells = "cMOP"
donor3 <- rbind(donor3_NCP1,donor3_NCP2,donor3_NCP3,donor3_NCP4,donor3_cMOP)

meta.data = data.frame(barcorde=rownames(donor3),Cell_Index=donor3$Cell_Index,Cells = donor3$Cells)
rownames(meta.data) <- meta.data$Cell_Index
UMI = donor3[,3:ncol(donor3)-1]
rownames(UMI) <- donor3$Cell_Index
UMI <- as.sparse(UMI)

donor3_seurat <- CreateSeuratObjectRhapsody(umi=UMI,metadata=meta.data,donor=3,sample="28072020")
saveRDS(donor3_seurat, file = "donor3_seurat.rds")

########################### Integration of data
#### function to calculate percentile
function_percentile <- function(seurat){
  x = seurat@meta.data$nFeature_RNA
  p <- quantile(x, c(.01, .99))
  print(p)
}
#########
BM.integrated = merge(x=donor1_seurat,y=c(donor2_seurat,donor3_seurat))
Idents(BM.integrated) <- "orig.ident"
VlnPlot(BM.integrated,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size = 0)

rm(donor1_seurat);rm(donor2_seurat);rm(donor3_seurat)

function_percentile(BM.integrated)  %>% kbl() %>%  kable_styling()
BM.integrated <- subset(BM.integrated, subset = nFeature_RNA > 805 & percent.mt < 25 & nFeature_RNA < 4762)
BM.integrated$donor = as.factor(BM.integrated$donor)
options(future.globals.maxSize = 20000 * 1024^2)

variables = c("percent.mt","CC.Difference")
BM.list <- SplitObject(BM.integrated, split.by = "donor")
BM.list <- BM.list[c("1","2","3")]
for (i in 1:length(BM.list)) {
  BM.list[[i]] <- SCTransform(BM.list[[i]], verbose = TRUE,vars.to.regress=variables)
}

BM.features <- SelectIntegrationFeatures(object.list = BM.list, nfeatures = 3000)
BM.list <- PrepSCTIntegration(object.list = BM.list, anchor.features = BM.features, verbose = FALSE)
reference_dataset <- which(names(BM.list) == "2")
BM.anchors <- FindIntegrationAnchors(object.list = BM.list, normalization.method = "SCT",
                                     anchor.features = BM.features, reference = reference_dataset)
BM.integrated <- IntegrateData(anchorset = BM.anchors, normalization.method = "SCT", verbose = FALSE)
saveRDS(BM.integrated, file = "precursors.integrated_BD_p1_p99_Mt_CC_genes.rds")

BM.integrated <- readRDS(file = "precursors.integrated_BD_p1_p99_Mt_CC_genes.rds")
