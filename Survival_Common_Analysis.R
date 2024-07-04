setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
diffMeth <- read.csv("DiffMeth_No_chrXYM_P.val_0.01.txt", header = TRUE, sep = ",", row.names = 1)
Survival.CpG_0.01 <- read.csv("Survival_CpGs_P_0.01.txt", header = TRUE, sep = ",", row.names = 1)
Spearman.Protein.Neg <- read.csv("spearman_Negative_Protein.txt", header = TRUE, sep = ",")
diffExp <- read.csv("diffExp_DESeq2.txt", header = TRUE, sep = ",", row.names = 1)
Survival.DEG_0.01 <- read.csv("Survival_Protein_Coding_P_0.01.txt", header = TRUE, sep = ",", row.names = 1)
## Below command will work multiple with rowname.
## But Spearman correlation matrix have no rownames, because CpG are repeated.
list_of_data = list(diffMeth, Survival.CpG_0.01)
common_names = Reduce(intersect, lapply(list_of_data, rownames))
Common_CpGs <- intersect(common_names, Spearman.Protein.Neg$Probe)
Common_CpGs.gene <- as.character(unique(Spearman.Protein.Neg[Spearman.Protein.Neg$Probe %in% Common_CpGs,]$gene))
intersect(Common_CpGs.gene, as.character(diffExp$Gene_Symbol)) ## BAIAP2  STX1A 
intersect(Common_CpGs.gene, rownames(Survival.DEG_0.01)) # IRX2  "ZFP28


#####################################################
Survival.CpG_0.05 <- read.csv("Survival_CpGs_P_0.05.txt", header = TRUE, sep = ",", row.names = 1)
Survival.DEG_0.05 <- read.csv("Survival_Protein_Coding_P_0.05.txt", header = TRUE, sep = ",", row.names = 1)

list_of_data.0.5 = list(diffMeth, Survival.CpG_0.05)
common_names.0.5 = Reduce(intersect, lapply(list_of_data.0.5, rownames))
Common_CpGs.0.5 <- intersect(common_names.0.5, Spearman.Protein.Neg$Probe)
Common_CpGs.gene.0.5 <- as.character(unique(Spearman.Protein.Neg[Spearman.Protein.Neg$Probe %in% Common_CpGs.0.5,]$gene))
intersect(Common_CpGs.gene.0.5, as.character(diffExp$Gene_Symbol)) #"BAIAP2" "MROH6"  "STX1A"
intersect(Common_CpGs.gene, rownames(Survival.DEG_0.05)) #"ASCL2"   "GNAL"    "IRX2"    "PCDHGB7" "ZFP28" 

###################################################
Survival.miRNA_0.01 <- read.csv("Survival_miRNA_P_0.01.txt", header = TRUE, sep = ",", row.names = 1)
diffExp.miRNA <- read.csv("diffExp_miRNA_DESeq2.txt", header = TRUE, sep = ",", row.names = 1)
intersect(rownames(diffExp.miRNA), rownames(Survival.miRNA_0.01)) #"hsa-mir-196b"

###################################################
save.image("Survival_Common_Analysis.RData")
