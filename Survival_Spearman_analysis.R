setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
survival.cpg.0.05 <- read.csv("Survival_CpGs_P_0.05.txt")
survival.cpg.0.01 <- read.csv("Survival_CpGs_P_0.01.txt")
survival.coding.0.05 <- read.csv("Survival_Protein_Coding_P_0.05.txt")
survival.coding.0.01 <- read.csv("Survival_Protein_Coding_P_0.01.txt")
survival.miRNA.0.01 <- read.csv("Survival_miRNA_P_0.01.txt")
survival.miRNA.0.05 <- read.csv("Survival_miRNA_P_0.05.txt")
survival.lncRNA.0.05 <- read.csv("Survival_lncRNA_P_0.05.txt")
survival.lncRNA.0.01 <- read.csv("Survival_lncRNA_P_0.01.txt")

DiffMeth <- read.csv("DiffMeth_No_chrXYM_P.val_0.01_UCSC.txt")
DiffMeth.FDR.0.005 <- DiffMeth[DiffMeth$adj.P.Val <= 0.005,]

Spearman.protein.pos <- read.csv("spearman_Postive_Protein_11_6_2018.txt")
Spearman.protein.neg <- read.csv("spearman_Negative_Protein_11_6_2018.txt")
spearman.protein.all <- rbind(Spearman.protein.neg, Spearman.protein.pos)

#spearman.protein.all <- rbind(Spearman.protein.neg, Spearman.protein.pos)
spearman.protein.all.BH.0.005 <- spearman.protein.all[spearman.protein.all$adjP < 0.005,]
spearman.protein.all.BH.0.005.pos <- spearman.protein.all.BH.0.005[spearman.protein.all.BH.0.005$cor.value > 0,]
spearman.protein.all.BH.0.005.neg <- spearman.protein.all.BH.0.005[spearman.protein.all.BH.0.005$cor.value < 0,]

spearman.protein.all.BH.0.005.rho0.5 <- spearman.protein.all.BH.0.005[abs(spearman.protein.all.BH.0.005$cor.value) >= 0.5,]
spearman.protein.all.BH.0.005.rho0.5.pos <- spearman.protein.all.BH.0.005.rho0.5[spearman.protein.all.BH.0.005.rho0.5$cor.value > 0,]
spearman.protein.all.BH.0.005.rho0.5.neg <- spearman.protein.all.BH.0.005.rho0.5[spearman.protein.all.BH.0.005.rho0.5$cor.value < 0,]


Spearman.lncRNA.pos <- read.csv("spearman_Postive_lincRNA_11_6_2018.txt")
Spearman.lncRNA.neg <- read.csv("spearman_Negative_lincRNA_11_6_2018.txt")
Spearman.lncRNA.all <- rbind(Spearman.lncRNA.neg, Spearman.lncRNA.pos)
Spearman.lncRNA.all.BH.0.005 <- Spearman.lncRNA.all[Spearman.lncRNA.all$adjP <= 0.005,]

Spearman.lncRNA.all.BH.0.005.neg <- Spearman.lncRNA.all.BH.0.005[Spearman.lncRNA.all.BH.0.005$cor.value < 0,]
Spearman.lncRNA.all.BH.0.005.pos <- Spearman.lncRNA.all.BH.0.005[Spearman.lncRNA.all.BH.0.005$cor.value > 0,]

Spearman.lncRNA.all.BH.0.005.rho.0.05 <- Spearman.lncRNA.all.BH.0.005[abs(Spearman.lncRNA.all.BH.0.005$cor.value)>= 0.5,]
Spearman.lncRNA.all.BH.0.005.rho.0.05.pos <- Spearman.lncRNA.all.BH.0.005[Spearman.lncRNA.all.BH.0.005$cor.value >= 0.5,]
Spearman.lncRNA.all.BH.0.005.rho.0.05.neg <- Spearman.lncRNA.all.BH.0.005[Spearman.lncRNA.all.BH.0.005$cor.value <= -0.5,]



diffExp.miRNA <- read.csv("diffExp_DESeq_miRNA_10_19_2018.txt")
diffExp.gene <- read.csv("Nitish_diffExp_DESeq_10_19_2018.csv")
#diffExp.miRNA <- read.csv("diffExp_miRNA_DESeq2.txt")


dim(spearman.protein.all.BH.0.005.neg[spearman.protein.all.BH.0.005.neg$Probe %in% intersect(spearman.protein.all.BH.0.005.neg$Probe, survival.cpg.0.05$X),])
dim(spearman.protein.all.BH.0.005.rho0.5.neg[spearman.protein.all.BH.0.005.rho0.5.neg$Probe %in% intersect(spearman.protein.all.BH.0.005.rho0.5.neg$Probe, survival.cpg.0.05$X),])
dim(spearman.protein.all.BH.0.005.rho0.5.neg[spearman.protein.all.BH.0.005.rho0.5.neg$Probe %in% intersect(spearman.protein.all.BH.0.005.rho0.5.neg$Probe, survival.cpg.0.01$X),])
intersect(spearman.protein.all.BH.0.005.rho0.5.neg$Probe, survival.cpg.0.05$X)
#spearman.protein.all.BH.0.005.rho0.5.neg[spearman.protein.all.BH.0.005.rho0.5.neg$Probe %in% c("cg26696870", "cg04956949", "cg05292954"),]
intersect(spearman.protein.all.BH.0.005.rho0.5.neg$Probe, survival.cpg.0.05$X)
intersect(spearman.protein.all.BH.0.005.rho0.5.neg$Probe, survival.cpg.0.01$X)


dim(Spearman.lncRNA.all.BH.0.005.neg[Spearman.lncRNA.all.BH.0.005.neg$Probe %in% intersect(Spearman.lncRNA.all.BH.0.005.neg$Probe, survival.cpg.0.05$X),])
dim(Spearman.lncRNA.all.BH.0.005.neg[Spearman.lncRNA.all.BH.0.005.neg$Probe %in% intersect(Spearman.lncRNA.all.BH.0.005.neg$Probe, survival.cpg.0.01$X),])
intersect(Spearman.lncRNA.all.BH.0.005.neg$Probe, survival.cpg.0.05$X)
intersect(Spearman.lncRNA.all.BH.0.005.neg$Probe, survival.cpg.0.01$X)
intersect(Spearman.lncRNA.all.BH.0.005.rho.0.05$Probe, survival.cpg.0.05$X)
intersect(Spearman.lncRNA.all.BH.0.005.rho.0.05$Probe, survival.cpg.0.01$X)
intersect(Spearman.lncRNA.all.BH.0.005.rho.0.05.neg$Probe, survival.cpg.0.01$X)
intersect(Spearman.lncRNA.all.BH.0.005.rho.0.05.pos, survival.cpg.0.01$X)


dim(survival.coding.0.05[survival.coding.0.05$X %in% diffExp.gene$Gene_Symbol,])
dim(survival.coding.0.01[survival.coding.0.01$X %in% diffExp.gene$Gene_Symbol,])
diffExp.miRNA[diffExp.miRNA$X %in% survival.miRNA.0.01$X,]
diffExp.miRNA[diffExp.miRNA$X %in% survival.miRNA.0.05$X,]

dim(survival.cpg.0.05[survival.cpg.0.05$X %in% intersect(DiffMeth.FDR.0.005$Row.names, spearman.protein.all.BH.0.005.neg$Probe),])
dim(survival.cpg.0.01[survival.cpg.0.01$X %in% intersect(DiffMeth.FDR.0.005$Row.names, spearman.protein.all.BH.0.005.rho0.5.neg$Probe),])

dim(survival.cpg.0.05[survival.cpg.0.05$X %in% intersect(DiffMeth.FDR.0.005$Row.names, Spearman.lncRNA.all.BH.0.005.neg$Probe),])
dim(survival.cpg.0.01[survival.cpg.0.01$X %in% intersect(DiffMeth.FDR.0.005$Row.names, spearman.protein.all.BH.0.005.rho0.5.neg$Probe),])

save.image("Spearman_Survival_Table.RData")
