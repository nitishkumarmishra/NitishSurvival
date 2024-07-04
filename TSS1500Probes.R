## This program developed by using getNearestTSS command from FDb.InfiniumMethylation.hg19
## Unlike Hui Shen program it's working for all gene, even if there is no CpG in genomica range of given gene
library('org.Hs.eg.db')
library(Homo.sapiens)
library(biomaRt)
library(FDb.InfiniumMethylation.hg19)
#######################################################
genes <- read.csv("TSS-HGNC-name.csv") ### TSS-HGNC.CSV is the list of TCGA genes in RNASeqV2 level-3 data
yy<-as.character(genes$hgnc, na.rm=FALSE)
keytypes(Homo.sapiens)
hm450 <- get450k()
TSS.nearest <- getNearestTSS(hm450)##  Make file for nearest TSS for each CpGs
TSS.nearest.uniq.gene <- unique(TSS.nearest$nearestGeneSymbol)
gene <- intersect(yy, TSS.nearest.uniq.gene)
#tmp <- TSS.nearest[grep("^A1CF$", TSS.nearest$nearestGeneSymbol),]
#rownames(tmp[tmp$distance <=1500,])

getProbes<-function(geneID){
  #tmp <- TSS.nearest[grep(geneID, TSS.nearest$nearestGeneSymbol, fixed = TRUE),]## This command have problem in NAT1, NAT10 and NAT14
  # In case of NAT1 it give list of all NAT1 (NAT1, NAT10, NAT14 etc)
  tmp <- TSS.nearest[which(TSS.nearest$nearestGeneSymbol==geneID),]
  if(is.null(tmp)){
    probes <- NULL}
    else{
  probes <- rownames(tmp[tmp$distance <=1500,])
    }
  return(probes)
}
#listprobes <- sapply(gene1, getProbes)
listprobes <- sapply(gene, getProbes)
# getExpr<-function(geneID){
#   tmp <- Expr[which(rownames(Expr)==geneID),]
#   return(tmp)
# }
# 
# getMeth<-function(probeID){
#   tmp <- Meth[which(rownames(Meth)==probeID),]
#   return(tmp)
# }
# 
# getMethVal<-function(geneID){
#   probes <- unlist(sapply(geneID, getProbes))
#   meth <- sapply(probes, getMeth)
#   return(meth)
# }
# 
# can <- cancor(t(matrix(unlist(tmp1), ncol = ncol(Expr))),t(matrix(unlist(tmp), ncol = ncol(Meth))))
# 
# getCanCorr <- function(geneID){
#   meth <- getMethVal(geneID)
#   meth <- t(matrix(unlist(meth), ncol = ncol(Meth)))
#   expr <- sapply(geneID, getExpr)
#   expr <- t(matrix(unlist(expr), ncol = ncol(Expr)))
#   can <- cancor(expr, meth)
#   return(can)
# }