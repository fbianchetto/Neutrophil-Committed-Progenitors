fpkm_to_tpm = function(fpkm){
  exp(log(fpkm)-log(colSums(fpkm))+log(1e6))
}

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

G_sequenced <- function(genes,ensembl,DESeqOb) {
  data_genes <- c()
  for (i in 1:length(genes)) {
    tmp = subset(featureData,featureData$gene_biotype==genes[i])
    tmp = DESeqOb[rownames(DESeqOb)%in%tmp$rownames.ddsTxi.,]
    print(paste(nrow(tmp),genes[i],sep=" "))
    mapped = colSums(assay(tmp))
    data_genes = rbind(data_genes,mapped)
    
  }
  return(data_genes)
  rownames(data_genes) = genes
}

##### Library
library(gplots)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(plyr)
colors <- brewer.pal(12,"Set3")
#####

plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                 panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))

my_pcaplot <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE,title=NULL,
                        pcX = 1, pcY = 2,text_labels=TRUE,point_size=3,
                        ellipse=TRUE,ellipse.prob=0.95) # customized principal components
{
  
  
  colourCount =length(unique(colData(x)$group)) 
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(assay(x)[select,]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(x)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, pcX], PC2 = pca$x[, pcY], group = group,
                  intgroup.df, names = colData(x)$group)
  colnames(d)[1] <- paste0("PC",pcX)
  colnames(d)[2] <- paste0("PC",pcY)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  # clever way of positioning the labels - worked good, then no need with ggrepel
  d$hjust <- ifelse((sign(d[,paste0("PC",pcX)])==1),0.9,0.1)# (1 + varname.adjust * sign(PC1))/2)
  
  g <- ggplot(data = d, aes_string(x = paste0("PC",pcX), y = paste0("PC",pcY), color="group")) +
    geom_point(size = point_size) +
    xlab(paste0("PC",pcX,": ", round(percentVar[pcX] * 100,digits = 2), "% variance")) +
    ylab(paste0("PC",pcY,": ", round(percentVar[pcY] * 100,digits = 2), "% variance")) 
  ## plot confidence ellipse
  # credit to vince vu, author of ggbiplot
  if(ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(d, 'group', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x[[paste0("PC",pcX)]], x[[paste0("PC",pcY)]]))
      mu <- c(mean(x[[paste0("PC",pcX)]]), mean(x[[paste0("PC",pcY)]]))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$group[1])
    })
    # names(ell)[1:2] <- c('xvar', 'yvar')
    if(nrow(ell)>0) {
      g <-g + geom_path(data = ell, aes_string(x="X1",y="X2",color = "groups", group = "groups"))
    }
  }
  
  if(text_labels)
    g <- g + geom_label_repel(mapping = aes_string(label="names",fill="group"),
                              color="white", show.legend = TRUE) 
  if(!is.null(title)) g <- g + ggtitle(title)
  g <- g + theme_bw()
  g + plotTheme 
}


TF_function = function(data1,data2){
    TF = subset(data1,data1$ensembl_gene_id%in%data2$ensembl_gene_id)
    TF$is_TF = "SI"
    NO = subset(data1,!data1$ensembl_gene_id%in%data2$ensembl_gene_id)
    NO$is_TF = "NO"
    data_genes = rbind(TF,NO)
    data_genes = data_genes[order(data_genes$ensembl_gene_id),]
    return(data_genes)
}

GO_2_genes <- function(data) {
  GO = data$geneID
  names(GO) = rownames(data)
  tab = as.list(GO)
  for (i in 1:length(tab)) {
    genes = tab[[i]]
    tab[[i]] = unlist(strsplit(genes, "/"))
  }
  return(tab)
}

GenesFromGO <- function(data) {
  GO = data$geneID
  tab = as.list(GO)
  for (i in 1:length(tab)) {
    genes = tab[[i]]
    tab[[i]] = unlist(strsplit(genes, "/"))
    genes = unique(unlist(tab))
  }
  return(genes)
}


GO_heatmap = function(data){
  color = c("white","darkgreen")
  list= GO_2_genes(data)
  genes = GenesFromGO(data)
  d <- c()
  for (i in genes){
    for (j in 1:length(list)){
      UT = genes_GO[[j]][which(genes_GO[[j]]==i)]
      if(length(UT)>0){UT = c(1)} else {UT = c(0)}
      d = rbind(d,UT)
    }
  }
  d
  e = t(matrix(d, nrow=length(genes),length(list)))
  c = t(e)
  rownames(c) = genes
  colnames(c) = names(list)
  pheatmap(t(c),color = color,legend_breaks = c(0,1))
#  return(c)
}

dir_work = "/media/simple/Volume2_4T/RNA_SEQ_our/Smart-seq2/Patrizia"
setwd(dir_work)
file = file.path(paste(dir_work,"protein_class_Predicted.tsv",sep="/"))

Protein = read.delim(file=file, header=TRUE,sep="\t")
ProteinDataBase = Protein[,c(3,8,13)]

Protein_function = function(dataFrame,Protein){
  DataProtein = Protein[order(Protein$Ensembl),]
  SI = merge(dataFrame,DataProtein,by.x="ensembl_gene_id",by.y="Ensembl",all=FALSE)
  SI$isProteinMembrane="SI"
  NO = subset(dataFrame,!(dataFrame$ensembl_gene_id)%in%Protein$Ensembl)
  NO$Evidence = "NO"
  NO$Subcellular.location = "NO"
  NO$isProteinMembrane="NO"
  data_genes = rbind(SI,NO)
  data_genes = data_genes[order(data_genes$ensembl_gene_id),]
  return(data_genes)
}

####
GO_Profile_Analysis <- function(DESeq2Res,sigGenes) {
  fg = mapIds(org.Hs.eg.db,keys=sigGenes,keytype = "ENSEMBL",column = "SYMBOL")
  bg = mapIds(org.Hs.eg.db,keys=rownames(DESeq2Res),keytype = "ENSEMBL",column = "SYMBOL")
  GO =  gprofiler(query=fg, organism = "hsapiens",ordered_query=T,custom_bg= bg,hier_filtering="strong")
  GO = GO[order(GO$p.value),]
  return(GO)
}

library(topGO)

GoTopAnalysis <- function(DESeq2Res,sigGenes){
  overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
  
  sig_idx <- match(sigGenes, rownames(overallBaseMean))
  
  backG <- c()
  
  for(i in sig_idx){
    ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
    
  }
  
  backG <- unique(backG)
  backG <- rownames(overallBaseMean)[backG]
  
  backG <- setdiff(backG,  sigGenes)
  length(backG)
  
  onts = c( "MF", "BP", "CC" )
  
  geneIDs = rownames(overallBaseMean)
  inUniverse = geneIDs %in% c(sigGenes,  backG) 
  inSelection =  geneIDs %in%sigGenes 
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- geneIDs[inUniverse]
  
  tab = as.list(onts)
  names(tab) = onts
  for(i in 1:3){
    
    ## prepare data
    tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
                annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
    
    ## run tests
    resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
    
    showSigOfNodes(tgd,score(resultTopGO.classic),firstSigNodes = 10, useInfo = 'all')
    
    ## look at results
    tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                          Fisher.classic = resultTopGO.classic,
                          orderBy = "Fisher.classic" , topNodes = 10)
    
  }
  return(tab)
}

library(clusterProfiler)
#cluterProfilerGO <-function(DESeq2Res,sigGenes){
#  onts = c( "MF", "BP", "CC")
# tab = as.list(onts)
#  names(tab) = onts
#  for(i in 1:3){
#    ego <- enrichGO(gene=sigGenes,'org.Hs.eg.db',universe= rownames(DESeq2Res),minGSSize = 10,
#                    keyType = 'ENSEMBL',ont=onts[i],pAdjustMethod="BH",pvalueCutoff=0.1,
#                    qvalueCutoff=0.2,readable=TRUE)
#    tab[[i]] = head(data.frame(ego),10)
#    tab[[i]] = data.frame(ego)
#  }
#  return(tab)
#}


cluterProfilerGO_BP <-function(DESeq2Res,sigGenes){
	ego <- enrichGO(gene=sigGenes,'org.Hs.eg.db',universe= rownames(DESeq2Res),minGSSize = 10,
                    keyType = 'ENSEMBL',ont="BP",pAdjustMethod = "BH",pvalueCutoff=0.05,
                    qvalueCutoff=0.2,readable=TRUE,pool = FALSE)
	return(ego)
}

#Plot_GO_ClusterProfiler <- function(df){
#  df = df[order(df$p.adjust),]
#  log10(df$p.adjust)
  
#  plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                   panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))
  
#  p <- ggplot(data=df,aes(x=reorder(Description,-log10(pvalue)),y=Count)) +
#    coord_flip() + scale_fill_gradient(low="blue",high="red",guide="colourbar") +
#    geom_bar(stat="identity",aes(fill =-log10(p.adjust))) + plotTheme +
#    theme(axis.text = element_text(size=12,colour="black")) 
#  return(p)
#}


Plot_GO_ClusterProfiler <- function(df,labels=NULL,breaks=NULL){
  
  df = as.data.frame(df)
  df = subset(df,df$p.adjust<0.05)		
  df = df[order(df$p.adjust),]
  log10(df$p.adjust)
  
  plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))
  
  p <- ggplot(data=df,aes(x=reorder(Description,-log10(p.adjust)),y=Count)) +
    coord_flip() + scale_fill_gradient(low="blue",high="red",guide="colourbar",labels=labels,breaks=breaks) +
    geom_bar(stat="identity",aes(fill =-log10(p.adjust))) + plotTheme +
    theme(axis.text = element_text(size=12,colour="black")) 
  return(p)
}

Plot_GO_TopGO <- function(df){
  colnames(df) = c("GO.ID","Term","Annotated","counts","Expected","Rank.in.Fisher.classic","Fisher.elim","Fisher.classic")
  df = df[order(df$Rank.in.Fisher.classic),]
  Fisher.classic = as.numeric(df$Fisher.classic)
  df$pvalue = as.numeric(Fisher.classic)
  df = subset(df,df$pvalue < 0.05)	
  
  plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))
  
  p <- ggplot(data=df,aes(x=reorder(Term,-log10(pvalue)),y=counts)) + 
    coord_flip() + plotTheme + theme(axis.text = element_text(size=12,colour="black")) +
    geom_bar(stat="identity",aes(fill =-log10(pvalue))) + plotTheme + xlab("-log10 Pvalue") + ylab("number of genes")
  return(p)
}

Plot_GO_GoProfiler <- function(df){
  df = df[order(df$p.value),]
  
  plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))
  
  p <- ggplot(data=df,aes(x=reorder(term.name,-log10(p.value)),y=overlap.size)) +
    coord_flip() + plotTheme + theme(axis.text = element_text(size=12,colour="black")) +
    geom_bar(stat="identity",aes(fill =-log10(p.value))) + plotTheme + ylab("Counts")
  return(p)
}

########
my_pcaplot <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE,title=NULL,
                        pcX = 1, pcY = 2,text_labels=TRUE,point_size=3,
                        ellipse=TRUE,ellipse.prob=0.95) # customized principal components
{
  colourCount =length(unique(colData(x)$group)) 
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(assay(x)[select,]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(x)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, pcX], PC2 = pca$x[, pcY], group = group,
                  intgroup.df, names = colData(x)$ALT_name)
  colnames(d)[1] <- paste0("PC",pcX)
  colnames(d)[2] <- paste0("PC",pcY)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  # clever way of positioning the labels - worked good, then no need with ggrepel
  d$hjust <- ifelse((sign(d[,paste0("PC",pcX)])==1),0.9,0.1)# (1 + varname.adjust * sign(PC1))/2)
  
  g <- ggplot(data = d, aes_string(x = paste0("PC",pcX), y = paste0("PC",pcY), color="group")) +
    geom_point(size = point_size) +
    xlab(paste0("PC",pcX,": ", round(percentVar[pcX] * 100,digits = 2), "% variance")) +
    ylab(paste0("PC",pcY,": ", round(percentVar[pcY] * 100,digits = 2), "% variance")) 
  ## plot confidence ellipse
  # credit to vince vu, author of ggbiplot
  if(ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(d, 'group', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x[[paste0("PC",pcX)]], x[[paste0("PC",pcY)]]))
      mu <- c(mean(x[[paste0("PC",pcX)]]), mean(x[[paste0("PC",pcY)]]))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$group[1])
    })
    # names(ell)[1:2] <- c('xvar', 'yvar')
    if(nrow(ell)>0) {
      g <-g + geom_path(data = ell, aes_string(x="X1",y="X2",color = "groups", group = "groups"))
    }
  }
  
  if(text_labels)
    g <- g + geom_label_repel(mapping = aes_string(label="names",fill="group"),
                              color="white", show.legend = TRUE) 
  if(!is.null(title)) g <- g + ggtitle(title)
  g <- g + theme_bw()
  g + plotTheme 
}

######

GOfuncR_converterID <- function(genes,res){
  candi_gene_ids = mapIds(org.Hs.eg.db,keys=genes,keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
  candi_gene_ids <- candi_gene_ids[!is.na(candi_gene_ids)]
  background = mapIds(org.Hs.eg.db,keys=rownames(res),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
  bg_gene_ids <- background[!is.na(background)]
  bg_gene_ids = setdiff(bg_gene_ids,candi_gene_ids)
  is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
  genes = data.frame(gene_ids=c(candi_gene_ids, bg_gene_ids), is_candidate)
  return(genes)
}

#################################################################
#library(biomaRt)
# See available resources in Ensembl release 84
#listMarts(host='mar2016.archive.ensembl.org')
# Connect to the Ensembl Genes annotation release 84 for Bos taurus
#ensembl84 = useMart(host='http://mar2016.archive.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
## Download all the Ensembl gene annotations (no filtering)
#attributes = listAttributes(ensembl84 )
#allgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'),mart=ensembl84)
# Rename the gene identifier column to 'gene_id'
# This allows GOexpress to treat microarray and RNA-seq data identically
#colnames(allgenes.Ensembl)[1] = 'gene_id'
## Download all the gene ontology annotations (no filtering)
#allGO.Ensembl = getBM(attributes=c('go_id', 'name_1006', 'namespace_1003'),mart=ensembl84)
## Download all the mapping between gene and gene ontology identifiers
#GOgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'go_id'),mart=ensembl84)
# Rename the gene identifier column to 'gene_id'
#colnames(GOgenes.Ensembl)[1] = 'gene_id'
# Cleanup: remove some blank fields often found in both columns
#GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$go_id != '',]
#GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '',]

#save(GOgenes.Ensembl, file='/media/simple/Volume2_4T/RNA_SEQ_our/Smart-seq2/Giulia/R_scripts/GOgenes.Ensembl84.rda')
#save(allGO.Ensembl, file='/media/simple/Volume2_4T/RNA_SEQ_our/Smart-seq2/Giulia/R_scripts/allGO.Ensembl84.rda')
#save(allgenes.Ensembl, file='/media/simple/Volume2_4T/RNA_SEQ_our/Smart-seq2/Giulia/R_scripts/allgenes.Ensembl84.rda')

load(file='/media/simple/Volume2_4T/RNA_SEQ_our/Smart-seq2/Giulia/R_scripts/GOgenes.Ensembl84.rda')
load(file='/media/simple/Volume2_4T/RNA_SEQ_our/Smart-seq2/Giulia/R_scripts/allGO.Ensembl84.rda')
load(file='/media/simple/Volume2_4T/RNA_SEQ_our/Smart-seq2/Giulia/R_scripts/allgenes.Ensembl84.rda')

##############################
#' @title Chain topGO and REViGO analyses to produce treemaps
#'
#' @description This package takes a list of genes, a map file with the
#' correspondance between gene name and GO annotation and a prefix for the
#' output. It will then do a topGO analysis and send the results to the
#' REViGO (http://revigo.irb.hr/) website, to summarize the list of GO and
#' produce a treemap.
#'
#' @details The goal is from a list of genes and the corresponding GO map, to
#' be able to produce an enriched list of GO annotations and a treemap to
#' easily visualize. By default, the biological process is outputted, but
#' it can also be the Cellular Component or the Molecular Function. One can use
#' a installed db or a map file.
#' @param geneList The gene list must be a csv file without column name, each
#' line consisting of the gene name and a 1 or 0, separated by a ",". The 1 or
#' 0 corresponds to the fact that the gene is respectively selected or not.
#' What I mean by that is that the gene has previously been recognized by the
#' user as interesting, like belonging to a cluster, or selected by any other
#' way. It corresponds to the "Predefined list of interesting genes" from the
#' topGO vignette
#' (http://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf).
#' If you don't have a predefined list, please do the previous steps of the
#' topGO vignette (before 4.4).
#' @param prefix A prefix for the outputs.
#' @param mapFile A file containing the correspondance between the gene name and
#' a GO name, 1 per line, in the format :
#' GeneName<tabulation>GOName
#' It can also be a Db name, thus one needs to change the option mapOrDb.
#' @param ontology "BP", "CC" or "MF" for Biological Process, Cellular Component
#' or Molecular Function, this is the GO categories outputted by REViGO.
#' @param mapOrDb map if a map file is used, db if a database name is provided.
#' @param p the p-value for the weight Fisher test to accept the GO. The Ficher
#' test based on weights takes multiple testing into account directly without
#' further need to apply a fdr or bonferroni correction.
#' @return A csv file containing the enriched GO terms and a treemap pdf file
#' containing the image.
#' @import topGO
#' @examples
#' library(hgu133a.db)
#' selGenes <- sample(ls(hgu133aGO), 50)
#' allGenes <-  factor(as.integer(ls(hgu133aGO) %in% selGenes))
#' names(allGenes) <- ls(hgu133aGO)
#' topReviGO(allGenes, "toto", "hgu133a", mapOrDb = "db", p = 0.01)
#' @export
topReviGO <- function(geneList, prefix, mapFile, ontology = "BP",
                      mapOrDb = "map", p = 0.01){
  # Check that the geneList and prefix are provided
  if (missing(geneList) | missing(prefix) | missing(mapFile)){
    stop("geneList, prefix or mapFile is missing")
  }
  if (mapOrDb == "map"){
    # Loading the Potri map
    geneID2GO <- topGO::readMappings(file=mapFile)

    # Creation of the GOdata object
    GOdata <- methods::new("topGOdata", ontology = ontology, allGenes = geneList,
                  annot = topGO::annFUN.gene2GO, gene2GO = geneID2GO)
  } else if (mapOrDb == "db"){
    # Loading the db
    affyLib <- paste(mapFile, "db", sep = ".")
    library(package = affyLib, character.only = TRUE)

    # Creation of the GOdata object
    GOdata <- methods::new("topGOdata", ontology = ontology, allGenes = geneList,
                  annot = topGO::annFUN.db, affyLib = affyLib)

#    GOdata <- methods::new("topGOdata", ontology = ontology, allGenes = geneList,
 #                 annot = topGO::annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl")
 
  }  else {stop('mapFile option must be "map" or "db"')}
  # Calculation of the Fisher test weightCount
  test.stat <- methods::new("weightCount", testStatistic = topGO::GOFisherTest,
                   name = "Fisher test", sigRatio = "ratio")
  resultWeight <- topGO::getSigGroups(GOdata, test.stat)

  ## Showing test stats
  topGO::geneData(resultWeight)

  # Creation of the allRes object
  allRes <- topGO::GenTable(GOdata, weightFisher = resultWeight,
                     orderBy = "weightFisher", ranksOf = "weightFisher",
                     topNodes = length(topGO::score(resultWeight)))
  # Output the pvalues for the test
  utils::write.table(allRes, file=paste0(prefix, "_weightFisher.csv"),
                     quote=F, row.names=F, col.names=T)
  allResInf1 <- allRes[allRes$weightFisher < p,]
  if(nrow(allResInf1) == 0){
    stop("No gene has been found above the p value. Stopping now.
         Please try to increase the required p-value (option p).")
  }
  # Localization of the revigoDownload.py script
  revigoDownloadLocation <- paste(system.file(package="topReviGO"),
                                  "revigoDownload.py",
                                  sep="/")
  # separator <- if (.Platform$OS.type == "windows") ";" else ":"
  # paths <- unlist(strsplit(Sys.getenv("PATH"), separator))
  # whereRD <- sapply(paths,
  #                   function(x) file.exists(paste0(x, "/revigoDownload.py")))
  # revigoDownloadLocation <- paste0(paths[whereRD], "/revigoDownload.py")
  # Incorporation of the revigoDownload.py script
  aRevigorer = "aRevigorer.txt"
  utils::write.table(allResInf1[,c("GO.ID", "weightFisher")], file=aRevigorer,
              quote=F, row.names=F, col.names=F)
  system(command = paste0("python ", revigoDownloadLocation, " -tsap ",
                          prefix, " ", aRevigorer))
  file.remove(aRevigorer)
  revigo.data <- utils::read.csv(paste0(prefix, "_treemap.csv"), skip = 4)
  revigo.data$abslog10pvalue <- abs(as.numeric(as.character(
    revigo.data$log10pvalue)))
  revigo.data$freqInDbPercent <- as.numeric(gsub("%", "",
                                                 revigo.data$frequencyInDb))
  treemap::treemap(revigo.data, index = c("representative","description"),
          vSize = "abslog10pvalue", type = "categorical",
          vColor = "representative",
          title = paste0("REVIGO Gene Ontology treemap - ", prefix,
                         " - n=", sum(as.integer(geneList)-1)),
          inflate.labels = TRUE, lowerbound.cex.labels = 0,
          bg.labels = "#CCCCCCAA", position.legend = "none")
  grDevices::pdf(file=paste0(prefix, "_treemap.pdf"))
  treemap::treemap(revigo.data, index = c("representative","description"),
          vSize = "abslog10pvalue", type = "categorical",
          vColor = "representative",
          title = paste0("REVIGO Gene Ontology treemap - ", prefix,
                         " - n=", sum(as.integer(geneList)-1)),
          inflate.labels = TRUE, lowerbound.cex.labels = 0,
          bg.labels = "#CCCCCCAA", position.legend = "none")
  grDevices::dev.off()
}

#' csvFilePreparationForTopReviGo
#' @param csvFile The csv file containing the genes, that has to be prepared.
#' @return The genes list
#' @export
csvFilePreparationForTopReviGo <- function(csvFile){
  tmpList <- utils::read.csv(csvFile, header=F, row.names=1)
  geneList <- as.factor(tmpList$V2)
  names(geneList) <- rownames(tmpList)
  return(geneList)
}

#############
###Plot GO output in GOplot
### ClusterProfile


ClusterProfile_GO_ToGOplot = function(BP,MF,CC,ont){
## Biological process  
  BP$Category = rep("BP",nrow(BP))
  BP = BP[,c("Category","ID","Description","geneID","p.adjust")]
  colnames(BP) = c("Category","ID","Term","Genes","adj_pval")
## Molecular fuction
  MF$Category = rep("MF",nrow(MF))
  MF = MF[,c("Category","ID","Description","geneID","p.adjust")]
  colnames(MF) = c("Category","ID","Term","Genes","adj_pval")
## Cellular component
  CC$Category = rep("CC",nrow(CC))
  CC = CC[,c("Category","ID","Description","geneID","p.adjust")]
  colnames(CC) = c("Category","ID","Term","Genes","adj_pval")
  GO = rbind(BP,MF,CC)
  GO$Genes = as.factor(gsub("/",",",GO$Genes))
  GO = subset(GO,GO$Category%in%ont)
  return(GO)
}

ToGoPlot <- function(res,GO,genes,dds){
  onts = c("GO","genelist")
  tab = as.list(onts)
  names(tab) = onts
  data = data.frame(res,mcols(dds)$gene)
  data = subset(data,rownames(data)%in%genes)
  rownames(data) = data$gene
  colnames(data) = c("AveExpr","logFC","B","t","P.Value","adj.P.Val","ID")
  data = data[,(c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B"))]
  tab$GO = GO
  tab$genelist = data
  return(tab)
}

###
TopGOtoRegivo <- function(path,name,GO){
  GO <- GO[,c(1,8)]
  write.table(GO,file=paste(path,paste(name,"txt",sep="."),sep="/"),row.names = F,col.names = F,quote = F,sep="\t")
  return(GO)
}

GO_ProfiletoRegivo <- function(path,name,GO){
  GO <- GO[,c(9,3)]
  write.table(GO,file=paste(path,paste(name,"txt",sep="."),sep="/"),row.names = F,col.names = F,quote = F,sep="\t")
  return(GO)
}

##############################################################
go_enrichOutput <- function(go_enrich){
  onts = c("BP","MF","CC")
  tab = as.list(onts)
  names(tab) = onts
  data = go_enrich$results
  data= data[,c("ontology","node_id","node_name","raw_p_overrep","FWER_overrep")]
  data = subset(data,data$FWER_overrep<=0.05)
  BP = subset(data,data$ontology=="biological_process")
  tab$BP = BP[order(BP$FWER_overrep),]
  MF = subset(data,data$ontology=="molecular_function")
  tab$MF = MF[order(MF$FWER_overrep),]
  CC = subset(data,data$ontology=="cellular_component")
  tab$CC = CC[order(CC$FWER_overrep),]
  return(tab)
}

###################################
Plot_go_enrich <- function(df){
  colnames(df) = c("ontology","node_id","term.name","p.value","FWER_overrep")
  df = df[order(df$FWER_overrep),]
  plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))
  p <- ggplot(data=df,aes(x=reorder(term.name,-log10(p.value)),y=-log10(p.value))) +
    coord_flip() + plotTheme + theme(axis.text = element_text(size=12,colour="black")) +
    geom_bar(stat="identity",aes(fill =-log10(p.value))) + plotTheme
  return(p)
}

Plot_GOrilla <- function(df){
  df = df[,c(1:5,9)]
  colnames(df) = c("GO Term","Description","p.value","FDR q-value","Enrichment","Count")
  df = df[order(df$p.value),]
  
  plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))
  
  p <- ggplot(data=df,aes(x=reorder(Description,-log10(p.value)),y=Count)) +
    coord_flip() + scale_fill_gradient(low="blue",high="red",guide="colourbar") +
    geom_bar(stat="identity",aes(fill =-log10(p.value))) + plotTheme +
    theme(axis.text = element_text(size=12,colour="black")) 
  return(p)
}

################################################
cluterProfilerGOtoRevigo <- function(path,name,GO){
  GO <- GO[,c(1,5)]
  write.table(GO,file=paste(path,paste(name,"txt",sep="."),sep="/"),row.names = F,col.names = F,quote = F,sep="\t")
  return(GO)
}
#####################################################################################

GOplotdatGorilla <- function(GO_gorilla) {
  Category <- rep("BP",nrow(GO_gorilla))
  ID <- GO_gorilla$GO.Term
  # Delete everything after character for ID
  ID <- gsub("~.*", "", ID)
  Term <- GO_gorilla$Description
  # Delete everything before character for term
  Term <- gsub(".*~", "", Term)
  Count <- GO_gorilla$b
  Genes <- GO_gorilla$Genes
#  adj_pval <- GO_gorilla$FDR.q.value
  adj_pval <- GO_gorilla$P.value
  gopdat <- data.frame(Category, ID, Term,Count, Genes, adj_pval)
  colnames(gopdat) = c('Category','ID','Term','Count','Genes','adj_pval)')
#  gopdat <- subset(gopdat, adj_pval <= 0.05)
  return(gopdat)
}

GorillaGsub = function(GO){
  GO$Genes = gsub("[[:lower:]]","",GO$Genes)
  GO$Genes = gsub("[[:punct:]]","",GO$Genes)
  GO$Genes = gsub("\\b\\d+\\b","",GO$Genes)
  GO$Genes = gsub("( )\\1+", "\\1,",GO$Genes)
  GO$Genes = gsub(",$", "", GO$Genes)
  GO$Genes = gsub(" ,", ",",GO$Genes)
  GO$Genes = gsub(" ", "",GO$Genes)
  GO$Genes = as.factor(GO$Genes)
  #  rownames(GO) = GO$ID
  GO
}

ALIASToSYMBOL_Gorilla = function(GO,SigGenes) {
  ensembl = SigGenes
  SYMBOL <- mapIds(org.Hs.eg.db,keys=ensembl,keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
  SYMBOL <- SYMBOL[!is.na(SYMBOL)]
  names(SYMBOL) <- NULL
  GO = GO 
  data <- c()
  symbol <- c()
  for(i in 1:nrow(GO)) {
    genes =  GO$Genes
    tgenes <- as.character(unlist(strsplit(as.vector(genes[i]), ",")))
    S = mapIds(org.Hs.eg.db,keys=tgenes,keytype = "ALIAS",column = "SYMBOL",multiVals = "first")
    sig_genes = subset(SYMBOL,SYMBOL%in%S)
    d = cbind(GO[i,],SYMBOL=as.character(paste(sig_genes,collapse=",")))
    data = rbind(data,d)
    }
  colnames(data) = c("Category","ID","Term","Count","Genes","adj_pval","SYMBOL")
  data = data[,c("Category","ID","Term","Count","SYMBOL","adj_pval")]
  colnames(data) = c("Category","ID","Term","Count","Genes","adj_pval")
  return(data)
}

###### David
EnsemblToSYMBOL_david = function(GO,SigGenes) {
   GO = GO 
   data <- c()
   symbol <- c()
   for(i in 1:nrow(GO)) {
      genes =  GO$Genes
      tgenes <- as.character(unlist(strsplit(as.vector(genes[i]), ", ")))
      sig_genes = subset(SigGenes,SigGenes%in%tgenes)
      S = mapIds(org.Hs.eg.db,keys=sig_genes,keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
      d = cbind(GO[i,],as.character(paste(S,collapse=",")))
      data = rbind(data,d)
      }
   colnames(data) = c("Category","ID","Term","Count","Genes","adj_pval","SYMBOL")
   data = data[,c("Category","ID","Term","Count","SYMBOL","adj_pval")]
   colnames(data) = c("Category","ID","Term","Count","Genes","adj_pval")
   return(data)
}

GOplotdatDavid	<- function(david)	{
  Category	<- david$Category
  ID	<- david$Term
  #	Delete	everything after character for	ID
  ID	<- gsub("~.*","",ID)
  Term	<- david$Term
  #	Delete	everything before character for	term
  Term	<- gsub(".*~","",Term)
  Count <- david$Count
  Genes	<- david$Genes
  adj_pval <- david$Benjamini
  gopdat <- data.frame(Category,ID,Term,Count,Genes,adj_pval)
  return(gopdat)
}

GOplotdatGeneList <- function(res,genes,dds){
  dds = subset(dds,rownames(dds)%in%genes)	
  data = data.frame(res,mcols(dds)$gene)
  data = subset(data,rownames(data)%in%genes)
  rownames(data) =  NULL
  colnames(data) = c("AveExpr","logFC","B","t","P.Value","adj.P.Val","ID")
  data = data[,(c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B"))]
  data$AveExpr = log2(data$AveExpr)
  return(data)
}

Plot_David <- function(df){
  plotTheme =theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.background = element_blank(), panel.border = element_rect(colour="black",fill=NA,size=0.5))

  p <- ggplot(data=df,aes(x=reorder(Term,-log10(adj_pval)),y=Count)) +
    coord_flip() + scale_fill_gradient(low="blue",high="red",guide="colourbar") +
    geom_bar(stat="identity",aes(fill =-log10(adj_pval))) + plotTheme +
    theme(axis.text = element_text(size=12,colour="black")) 
  return(p)
}


########
mycircle_dat = function (terms, genes) 
{
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ", ")
  if (length(tgenes[[1]]) == 1) 
    tgenes <- strsplit(as.vector(terms$genes), ",")
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, 
                                                                genes$ID)])
  if (class(logFC) == "factor") {
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1
  zsc <- c()
  for (c in 1:length(count)) {
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    value[is.na(value)] <-1
    zsc <- c(zsc, sum(value)/sqrt(count[c]))
    s <- e + 1
  }
  if (is.null(terms$id)) {
    df <- data.frame(category = rep(as.character(terms$category), 
                                    count), term = rep(as.character(terms$term), count), 
                     count = rep(count, count), genes = as.character(unlist(tgenes)), 
                     logFC = logFC, adj_pval = rep(terms$adj_pval, count), 
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  else {
    df <- data.frame(category = rep(as.character(terms$category), 
                                    count), ID = rep(as.character(terms$id), count), 
                     term = rep(as.character(terms$term), count), count = rep(count, 
                                                                              count), genes = as.character(unlist(tgenes)), 
                     logFC = logFC, adj_pval = rep(terms$adj_pval, count), 
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}


Get_ontology_mart <- function(GO_term){ ####### https://support.bioconductor.org/p/90214/
  require(biomaRt)
  ensembl <- useMart("ensembl", host = "www.ensembl.org",dataset = "hsapiens_gene_ensembl") #### version 94
  attributes = listAttributes(ensembl)###
  GO_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','go_id','name_1006','namespace_1003'),mart = ensembl,filters = "go_parent_term",values=GO_term)
  print("Filtering genes with protein coding features")
  GO_genes = subset(GO_genes,GO_genes$gene_biotype=="protein_coding")
  print("Filtering GO term with biological process")
  GO_genes = subset(GO_genes,GO_genes$namespace_1003=="biological_process")
  return(GO_genes)
}

Get_TF <- function(GO_term){ ####### https://support.bioconductor.org/p/90214/
  require(biomaRt)
  ensembl <- useMart("ensembl", host = "www.ensembl.org",dataset = "hsapiens_gene_ensembl") #### version 94
  attributes = listAttributes(ensembl)###
  GO_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','go_id','name_1006','namespace_1003'),mart = ensembl,filters = "go",values=GO_term)
  print("Filtering genes with protein coding features")
  GO_genes = subset(GO_genes,GO_genes$gene_biotype=="protein_coding")
  print("Filtering GO term with biological process")
  GO_genes = subset(GO_genes,GO_genes$namespace_1003=="biological_process")
  return(GO_genes)
}

Coding_Filtering <- function(dds,version){
  dds = estimateSizeFactors(dds)
  v96 = "http://apr2019.archive.ensembl.org"
  v84 = "http://mar2016.archive.ensembl.org"
  v74 = "http://feb2014.archive.ensembl.org"
  us  = "uswest.ensembl.org"

  require(biomaRt)
  ensembl <- useMart("ensembl", host = version,dataset = "hsapiens_gene_ensembl") #### 
  attributes = listAttributes(ensembl)
  print("Getting annotations through biomart library...")
  if (version == v96) {
    print (paste("You are using ensembl version 96, URL",version,sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name'),mart = ensembl,useCache = FALSE)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  } else if (version == us){
    print (paste("You are using the emergency ensembl site, URL",version,sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name'),mart = ensembl,useCache = FALSE)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  } else if (version == v84){
    print (paste("You are using ensembl version 84, URL",version,sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name'),mart = ensembl,useCache = FALSE)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  } else {
    print (paste("You are using ensembl version 74, URL",version, sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_id','gene_biotype','description','chromosome_name'),mart = ensembl,useCache = FALSE)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  }
  annotation$name = gsub('\\[.*','',annotation$name)
  annotation = subset(annotation,annotation$ID%in%rownames(dds)) ### 57443
  annotation = annotation[order(annotation$ID),]
  dds = subset(dds,rownames(dds)%in%annotation$ID)
  featureData = data.frame(rownames(dds),annotation)
  mcols(dds) <- DataFrame(mcols(dds),SYMBOL=featureData$SYMBOL,gene_name=featureData$name,
                          chr=featureData$chr,gene_biotype=featureData$biotype)
  print("Filtering out mitochondrial and Y genes")
  autosome = subset(featureData,!featureData$chr%in%c("MT"))
  dds = dds[rownames(dds)%in%autosome$ID,]### exclude mithocondrial genes
  select = c("protein_coding","lincRNA")
  CodingGenes = subset(featureData,featureData$biotype%in%select); dim(CodingGenes) ### 27374 coding genes (lincRNA,protein_coding)
  print("keep only protein and lincRNA coding genes")
  dds = dds[rownames(dds)%in%CodingGenes$ID,]### keep only protein and lincRNA coding#
  return(dds)
}

############################################################################################
FilteringOutLowExpressedGenes = function(dds){
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
}
############################################################################################


Revigo = function(GO){
  GO_filtered = subset(GO,GO$qvalue <= 0.1)
  GO_filtered = data.frame(rownames(GO_filtered),GO_filtered$pvalue,row.names = NULL)
  clipr::write_clip(GO_filtered)
}

############################################################
ExpressedGenes_by_fpkm_dds <- function(dds,name,threshold=0.3){
  out = DESeq2::fpkm(dds)
  colnames(out) = name
  ###### from https://stackoverflow.com/questions/35205572/rowmeans-with-matching-column-names
  df = sapply(unique(colnames(out)), function(i) matrixStats::rowMedians(out[,colnames(out) == i]))
  rownames(df) = rownames(out)
  #### from https://stackoverflow.com/questions/38273643/r-keep-rows-with-at-least-one-column-greater-than-value
  df = df[apply(df>threshold,1,any),]
  keep = as.vector(rownames(df))
  dds = subset(dds,rownames(dds)%in%keep)
  return(dds)
}

######

ExpressedGenes_by_fpkm <- function(dds,name,threshold=0.3){
  out = DESeq2::fpkm(dds)
  colnames(out) = name
  #### from https://stackoverflow.com/questions/35205572/rowmeans-with-matching-column-names
  df = sapply(unique(colnames(out)), function(i) matrixStats::rowMedians(out[,colnames(out) == i]))
  rownames(df) = rownames(out)
  #### from https://stackoverflow.com/questions/38273643/r-keep-rows-with-at-least-one-column-greater-than-value
  df = df[apply(df>threshold,1,any),]
  keep = as.vector(rownames(df))
  dds = subset(dds,rownames(dds)%in%keep)
  out = subset(out,rownames(out)%in%keep)
  out = out[,order(colnames(out))]
  merged = data.frame(out,mcols(dds)[,1:4])
  fpkm = merged[order(rowSums(merged[,1:ncol(out)])),]
  return(fpkm)
}


mean_fpkm = function(dds,name,threshold=0.3){
  out = fpkm(dds)
  colnames(out) = name
  #### from https://stackoverflow.com/questions/35205572/rowmeans-with-matching-column-names
  df = sapply(unique(colnames(out)), function(i) matrixStats::rowMedians(out[,colnames(out) == i]))
  rownames(df) = rownames(out)
  #### from https://stackoverflow.com/questions/38273643/r-keep-rows-with-at-least-one-column-greater-than-value
  df = df[apply(df>threshold,1,any),]
  keep = as.vector(rownames(df))
  dds = subset(dds,rownames(dds)%in%keep)
  merged = data.frame(df,mcols(dds)[,1:4])
  fpkm = merged[order(rowSums(merged[,1:ncol(df)])),]
  return(fpkm)
}


mean_tpm = function(tpm,name,threshold=0.3,dds){
  out = tpm
  colnames(out) = name
  #### from https://stackoverflow.com/questions/35205572/rowmeans-with-matching-column-names
  df = sapply(unique(colnames(out)), function(i) matrixStats::rowMedians(out[,colnames(out) == i]))
  rownames(df) = rownames(out)
  #### from https://stackoverflow.com/questions/38273643/r-keep-rows-with-at-least-one-column-greater-than-value
  df = df[apply(df>threshold,1,any),]
  keep = as.vector(rownames(df))
  dds = subset(dds,rownames(dds)%in%keep)
  df = subset(df,rownames(df)%in%rownames(dds))
  merged = data.frame(df,mcols(dds)[,1:4])
  tpm = merged[order(rowSums(merged[,1:ncol(df)])),]
  return(tpm)
}

##############################################################################################################
Coding_Filtering_tpm <- function(tpm,version){
  v96 = "http://apr2019.archive.ensembl.org"
  v84 = "http://mar2016.archive.ensembl.org"
  v74 = "http://feb2014.archive.ensembl.org"
  us  = "uswest.ensembl.org"
  require(biomaRt)
  ensembl <- useMart("ensembl", host = version,dataset = "hsapiens_gene_ensembl") #### version 84
  attributes = listAttributes(ensembl)
  print("Getting annotations through biomart library...")
  if (version == v96) {
    print (paste("You are using ensembl version 96, URL",version,sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name'),mart = ensembl)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  } else if (version == us){
    print (paste("You are using the emergency ensembl site, URL",version,sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name'),mart = ensembl)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  } else if (version == v84){
    print (paste("You are using ensembl version 84, URL",version,sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name'),mart = ensembl)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  } else {
    print (paste("You are using ensembl version 74, URL",version, sep=" : "))
    annotation = getBM(attributes=c('ensembl_gene_id','external_gene_id','gene_biotype','description','chromosome_name'),mart = ensembl)
    colnames(annotation) = c("ID","SYMBOL","biotype","name","chr")
  }
  annotation$name = gsub('\\[.*','',annotation$name)
  annotation = subset(annotation,annotation$ID%in%rownames(tpm)) ### 57443
  annotation = annotation[order(annotation$ID),]
  tpm = subset(tpm,rownames(tpm)%in%annotation$ID)
  featureData = data.frame(rownames(tpm),annotation)
  print("Filtering out mitochondrial and Y genes")
  autosome = subset(featureData,!featureData$chr%in%c("MT","Y","X"))
  tpm = tpm[rownames(tpm)%in%autosome$ID,]### exclude mithocondrial genes
  select = c("protein_coding","lincRNA")
  CodingGenes = subset(featureData,featureData$biotype%in%select); dim(CodingGenes) ### 27374 coding genes (lincRNA,protein_coding)
  print("keep only protein and lincRNA coding genes")
  tpm = tpm[rownames(tpm)%in%CodingGenes$ID,]### keep only protein and lincRNA coding#
  CodingGenes = CodingGenes[CodingGenes$ID%in%rownames(tpm),]
  df = data.frame(tpm,CodingGenes)
  return(df)
}


medianMatrix = function(mat,groups) {
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if(sum(A)==1) {
      mat[,A]
    } else {
      rowMedians(mat[,A],na.rm=T)
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}


