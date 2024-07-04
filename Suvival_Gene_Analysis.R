##  R code for the analysis of miR-22b, miR-551b targrt gene expression analysis
## Also hyper/hypomethylated CpGs in promoter region which are associated with survival of patients
setwd("F:/OneDrive - University of Nebraska Medical Center/Cholangiocarcinoma")
mir22 <- read.csv("miR-22 Target.txt", header = TRUE, sep = "\t")
diffExp <- read.csv("Supplem/Supplementary data/mRNA_edgeR_DEseq2.csv", header = TRUE, sep = ",", row.names = 1)
diffExp.mir22 <- diffExp[diffExp$symbol%in%mir22$Gene.Symbol,]
mir551b <- read.csv("miR-551b Target.txt", header = TRUE, sep = "\t")
grep("UPS6", diffExp$symbol)

hm450 <- read.table("F:/OneDrive - University of Nebraska Medical Center/Human450K-Extra/HuiShen/April 2018/hm450.hg38.manifest.gencode.v22.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names=5)
hm450 <- hm450[-c(1:5120),]
symbol <- strsplit(hm450$geneNames, split = ";")
type <- strsplit(hm450$transcriptTypes, split = ";")
transcript <- strsplit(hm450$transcriptIDs, split = ";")
positionTSS <- strsplit(hm450$distToTSS, split = ";")
feature <- strsplit(hm450$CGIposition, split = ";")
JHU.annotation <- data.frame(ProbeID = rep(rownames(hm450), sapply(positionTSS, length)), Chromosome=rep(hm450$CpG_chrm, sapply(positionTSS, length)),Strat=rep(hm450$CpG_beg, sapply(positionTSS, length)),End=rep(hm450$CpG_end, sapply(positionTSS, length)),Gene_Symbol = unlist(symbol),Transcript_ID = unlist(transcript), Gene_Type = unlist(type), Position_to_TSS=unlist(positionTSS), Feature_type = rep(hm450$CGIposition, sapply(positionTSS, length)))
JHU.annotation <- JHU.annotation[!is.na(JHU.annotation$Gene_Symbol),] ## Remove row which have <NA> in gene symbols
#tmp<- JHU.annotation[!is.na(JHU.annotation$Gene_Symbol),]
JHU.annotation$Position_to_TSS <-  varhandle::unfactor(JHU.annotation$Position_to_TSS)
JHU.annotation.TSS <-  JHU.annotation[((JHU.annotation$Position_to_TSS >= -1000) & (JHU.annotation$Position_to_TSS <= 500)),]
rm(symbol, type, transcript, positionTSS, feature, hm450)

Hyper <- list.files("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffMeth/DiffMeth_Hyper_ROC_Plot")
Hyper <- gsub("roc_", "", Hyper)
Hyper <- gsub(".png", "", Hyper)
Hyper.gene <- unique(JHU.annotation.TSS[JHU.annotation.TSS$ProbeID%in%Hyper,]$Gene_Symbol)
diffExp[diffExp$symbol%in%Hyper.gene,]

Hypo <- list.files("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffMeth/DiffMeth_Hypo_ROC_Plot")
Hypo <- gsub("roc_", "", Hypo)
Hypo <- gsub(".png", "", Hypo)
Hypo.gene <- unique(JHU.annotation.TSS[JHU.annotation.TSS$ProbeID%in%Hypo,]$Gene_Symbol)
diffExp[diffExp$symbol%in%Hypo.gene,]
save.image("Suvival_Gene_Analysis.RData")
