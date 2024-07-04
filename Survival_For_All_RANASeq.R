setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
load("PAAD.RData")
rm(list=setdiff(ls(), c("PAAD.Clinical", "PAAD.Meth", "PAAD.FPKM", "PAAD.FPKM.UQ", "PAAD.miRNA.RPM")))

####################### Expression data prepration ##########################
GDC.Harmonized <- read.csv("TCGA_GDC_Harmonided_Uniq_Status_GTF.txt", header = TRUE, sep = "\t", row.names = 1)
GDC.Harmonized.protein <- GDC.Harmonized[GDC.Harmonized$Gene_status=="protein_coding",]
GDC.Harmonized.lincRNA <- GDC.Harmonized[GDC.Harmonized$Gene_status=="lincRNA",]
#GDC.Harmonized.processed_pseudogene <- GDC.Harmonized[GDC.Harmonized$Gene_status=="processed_pseudogene",]

TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]

PAAD.FPKM.UQ.lncRNA <- PAAD.FPKM.UQ[rownames(GDC.Harmonized.lincRNA),]
PAAD.FPKM.UQ.protein <- PAAD.FPKM.UQ[rownames(GDC.Harmonized.protein),]

cancerID <- grep("01A", substr(colnames(PAAD.FPKM.UQ.lncRNA), 14, 16))
PAAD.FPKM.UQ.lncRNA <- PAAD.FPKM.UQ.lncRNA[,cancerID]

cancerID <- grep("01A", substr(colnames(PAAD.FPKM.UQ.protein), 14, 16))
PAAD.FPKM.UQ.protein <- PAAD.FPKM.UQ.protein[,cancerID]

cancerID <- grep("01A", substr(colnames(PAAD.miRNA.RPM), 14, 16))
PAAD.miRNA.RPM.cancer <- PAAD.miRNA.RPM[,cancerID]


colnames(PAAD.FPKM.UQ.lncRNA) <- substr(colnames(PAAD.FPKM.UQ.lncRNA), 1, 15)
colnames(PAAD.FPKM.UQ.protein) <- substr(colnames(PAAD.FPKM.UQ.protein), 1, 15)
colnames(PAAD.miRNA.RPM.cancer) <- substr(colnames(PAAD.miRNA.RPM.cancer), 1, 15)

PAAD.FPKM.UQ.lncRNA.ductal <-  PAAD.FPKM.UQ.lncRNA[,colnames(PAAD.FPKM.UQ.lncRNA)%in%TCGA_PAAD_select.ductal$Sample.ID]
PAAD.FPKM.UQ.protein.ductal <-  PAAD.FPKM.UQ.protein[,colnames(PAAD.FPKM.UQ.protein)%in%TCGA_PAAD_select.ductal$Sample.ID]
PAAD.miRNA.RPM.ductal <- PAAD.miRNA.RPM.cancer[,colnames(PAAD.miRNA.RPM.cancer)%in%TCGA_PAAD_select.ductal$Sample.ID]

## Remove genes which have more than 20% zeros
PAAD.FPKM.UQ.lncRNA.ductal <- PAAD.FPKM.UQ.lncRNA.ductal[apply(PAAD.FPKM.UQ.lncRNA.ductal == 0, 1, sum) <= ncol(PAAD.FPKM.UQ.lncRNA.ductal)*0.20, ]
PAAD.FPKM.UQ.protein.ductal <- PAAD.FPKM.UQ.protein.ductal[apply(PAAD.FPKM.UQ.protein.ductal == 0, 1, sum) <= ncol(PAAD.FPKM.UQ.protein.ductal)*0.20, ]
PAAD.miRNA.RPM.ductal <- PAAD.miRNA.RPM.ductal[apply(PAAD.miRNA.RPM.ductal ==0, 1, sum) <= ncol(PAAD.miRNA.RPM.ductal)*0.20,]

PAAD.FPKM.UQ.lncRNA.ductal <- log2(PAAD.FPKM.UQ.lncRNA.ductal + 1)
PAAD.FPKM.UQ.protein.ductal <- log2(PAAD.FPKM.UQ.protein.ductal + 1)
PAAD.miRNA.RPM.ductal <- log2(PAAD.miRNA.RPM.ductal + 1)

lncRNA.symbol <- subset(GDC.Harmonized.lincRNA, select = c(Gene_Symbol))
protein.symbol <- subset(GDC.Harmonized.protein, select = c(Gene_Symbol))

PAAD.FPKM.UQ.lncRNA.ductal <- merge(PAAD.FPKM.UQ.lncRNA.ductal, lncRNA.symbol, by="row.names")
rownames(PAAD.FPKM.UQ.lncRNA.ductal) <- PAAD.FPKM.UQ.lncRNA.ductal$Gene_Symbol
PAAD.FPKM.UQ.lncRNA.ductal$Row.names <- NULL

PAAD.FPKM.UQ.protein.ductal <- merge(PAAD.FPKM.UQ.protein.ductal, protein.symbol, by="row.names")
rownames(PAAD.FPKM.UQ.protein.ductal) <- PAAD.FPKM.UQ.protein.ductal$Gene_Symbol
PAAD.FPKM.UQ.protein.ductal$Row.names <- NULL

#############################################
############### For lncRNA ##################
source("plot_surv_RNAseq_function.R")

unlink("RNAseq_lncRNA_Surv/KM_Plot_*") 
results <- plot_surv(dir = "RNAseq_lncRNA_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.FPKM.UQ.lncRNA.ductal, Genelist = rownames(PAAD.FPKM.UQ.lncRNA.ductal), Survresult = FALSE, Median = TRUE, p.cut = 0.05)

results.DM.P.0.01 <- plot_surv(dir = "RNAseq_lncRNA_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.FPKM.UQ.lncRNA.ductal, Genelist = rownames(results), Survresult = TRUE, Median = TRUE, p.cut = 0.01)

# results.DM <- plot_surv(dir = "RNAseq_lncRNA_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.FPKM.UQ.lncRNA.ductal, Genelist = rownames(results), Survresult = FALSE, Median = TRUE, p.cut = 0.05)
# 
# results.DM.P.0.01 <- plot_surv(dir = "RNAseq_lncRNA_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.FPKM.UQ.lncRNA.ductal, Genelist = rownames(results.DM), Survresult = FALSE, Median = TRUE, p.cut = 0.01)
# 
# results.DM.P.0.01 <- plot_surv(dir = "RNAseq_lncRNA_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.FPKM.UQ.lncRNA.ductal, Genelist = rownames(results.DM.P.0.01), Survresult = TRUE, Median = TRUE, p.cut = 0.01)

#############################################
################ For RNAseq #################
unlink("RNAseq_Surv/KM_Plot_*") 
results.Gene <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.FPKM.UQ.protein.ductal, Genelist = rownames(PAAD.FPKM.UQ.protein.ductal), Survresult = FALSE, Median = TRUE, p.cut = 0.05)

results.Gene.DM.P.0.01 <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.FPKM.UQ.protein.ductal, Genelist = rownames(results.Gene), Survresult = TRUE, Median = TRUE, p.cut = 0.01)

#############################################
############# For miRNA #####################
unlink("RNAseq_miRNA_Surv/KM_Plot_*")
results.miRNA <- plot_surv(dir = "RNAseq_miRNA_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.miRNA.RPM.ductal, Genelist = rownames(PAAD.miRNA.RPM.ductal), Survresult = FALSE, Median = TRUE, p.cut = 0.05)

results.miRNA.DM.P.0.01 <- plot_surv(dir = "RNAseq_miRNA_Surv", clinical_patient = PAAD.Clinical, dataGE = PAAD.miRNA.RPM.ductal, Genelist = rownames(results.miRNA), Survresult = TRUE, Median = TRUE, p.cut = 0.05)

#############################################
save.image("Survival_For_All_RANASeq.RData")
