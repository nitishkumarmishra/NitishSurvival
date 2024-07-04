library(ComplexHeatmap)
library(circlize)

SB <- read.csv("Spina Bifida Individual 206 targets Oct 11 2016.txt to Nitish.txt", header = TRUE, sep = "\t", row.names = 1)
rownames(SB) <- SB$TargetID
BetaSB.normal <- SB[, grep("Control.AVG_Beta", colnames(SB))]
BetaSB.affected <- SB[, grep("Affected.AVG_Beta", colnames(SB))]
Meth <- cbind(BetaSB.affected, BetaSB.normal)
Meth$TargetID <- rownames(Meth)

FC <- read.csv("SB_group analysis 239 Report Oct 5 2016 Final.txt", header = TRUE, sep = "\t")
rownames(FC) <- FC$TargetID
MergeData <- merge(Meth, FC, by="TargetID")
annotation <- MergeData$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
foldchange <- MergeData$Fold.chance
relation <- MergeData$RELATION_TO_UCSC_CPG_ISLAND
AUC <- MergeData$AUC
logP <- MergeData$LOG10p
foldchange <- MergeData$Fold.chance
direction <- ifelse(MergeData$Fold.chance >1,"Hyper","Hypo")
#delta <- abs(rowMeans(BetaSB.affected) - rowMeans(BetaSB.normal))
delta <- abs(MergeData$X..Methylation......Cases - MergeData$X..Methylation......Control)/100
Meth1 <- Meth
Meth1$TargetID <- NULL
colnames(Meth1) <- gsub(".AVG_Beta", "", colnames(Meth1))



set.seed(1110)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("SB", 24), rep("Control", 16))), col = list(type = c("SB" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","lightseagreen","aquamarine4", "violet","violetred1","violetred3","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(4, 2, 0), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta beta", col = colorRamp2(c(0, 0.1, 0.3), c("gray", "red", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "-log10 FDR", col = colorRamp2(c(1, 20, 40), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))
  #Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("SB_206.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()

save.image(file = "SB_New.RData")
