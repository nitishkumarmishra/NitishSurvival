setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
load("PAAD.RData")
rm(list=setdiff(ls(), c("PAAD.Clinical", "PAAD.Meth", "PAAD.FPKM", "PAAD.FPKM.UQ", "PAAD.miRNA.RPM")))
exp_meth_protein_coding <- read.csv("exp_meth_protein_coding.txt", header = TRUE, check.names = FALSE)
exp <- unique(exp_meth_protein_coding[, 5:129])
rownames(exp) <- exp$gene_id
exp$gene_id <- NULL
colnames(exp) <- paste0(colnames(exp), "-01A")

source("plot_surv_RNAseq_function.R")


unlink("RNAseq_Surv/KM_Plot_*") 
results <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(exp), Survresult = FALSE, Median = TRUE, p.cut = 0.05)


unlink("RNAseq_Surv/KM_Plot_*") 
results.DM <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(results), Survresult = FALSE, Median = TRUE, p.cut = 0.05)

unlink("RNAseq_Surv/KM_Plot_*") 
results.DM.P.0.01 <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(results.DM), Survresult = TRUE, Median = TRUE, p.cut = 0.01)


unlink("RNAseq_Surv/KM_Plot_*") 
results.DM.P.0.01 <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(results.DM.P.0.01), Survresult = TRUE, Median = TRUE, p.cut = 0.01)

#############################################
save.image("RNASeq_survival_analysis.RData")
