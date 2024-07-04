#### R code for Kaplan-Meier survival curves for different patient groups
### We used median to make up and down genes
library(pacman)
p_load(limma, TCGAbiolinks, survival, survminer, edgeR)
#############################################################################################
clin <- GDCquery_clinic("TCGA-CHOL", type = "clinical", save.csv = FALSE)
cfu <- clin[clin[, "bcr_patient_barcode"] %in% substr(colnames(dataCancer), 1, 12), ]
### Next tow line is applicable if clin have colname "days_to_last_followup"
if ("days_to_last_followup" %in% colnames(cfu))
  colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <- "days_to_last_follow_up"
cfu <- as.data.frame(subset(cfu, select = c("bcr_patient_barcode",  "days_to_death", "days_to_last_follow_up", "vital_status")))
rownames(cfu) <- cfu$bcr_patient_barcode
notDead <- is.na(cfu$days_to_death)
if (any(notDead == TRUE)) {
  cfu[notDead, "days_to_death"] <- cfu[notDead, "days_to_last_follow_up"]
}
cfu$s <- grepl("dead|deceased", cfu$vital_status, ignore.case = TRUE)
cfu$vital_status <- ifelse(cfu$vital_status == "alive", 0, 1)

rna <- read.table('CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t', stringsAsFactors=FALSE, check.names = FALSE)
#rna <- rna[-1,] ## Use this line for normalized data and escape next line
rna <- rna[-c(1),which(rna[1,]=="raw_count")] ## Use this line for RSEM expected readcount data
rownames(rna) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(rna))
rna <- round(data.matrix(rna),0)
rna <- rna[-grep("\\?", rownames(rna)),] ## Remove genes which HGNC name is known i.e. ?

rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.25)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]

group1 <- TCGAquery_SampleTypes(colnames(rna), typesample = c("NT"))## Normal
group2 <- TCGAquery_SampleTypes(colnames(rna), typesample = c("TP"))## Tumor
rna <- cbind(rna[,group2], rna[,group1])

design <- model.matrix(~0 + factor(c(rep(2, length(group2)), rep(1, length(group1)))))
colnames(design) <- c("Normal", "Tumor")
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
cnts.dgelist <-  DGEList(rna, group=factor(c(rep("Tumor", length(group2)), rep("Normal", length(group1)))))
cnf <- calcNormFactors(cnts.dgelist, method = "TMM")
rna_vm <- voom(cnf, design, plot = FALSE)
rna_vm  <- rna_vm$E
rownames(rna_vm) <- sapply(strsplit(rownames(rna_vm),"\\|"),'[[',1)

Genelist <- rownames(rna_vm)
#Genelist <- intersect(rownames(rna_vm), Genelist)
dataCancer <- rna_vm[Genelist, group2]
dataNormal <- rna_vm[Genelist, group1]
colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)

tabSurv_Matrix <- matrix(0, nrow(rna_vm), 8)
colnames(tabSurv_Matrix) <- c("mRNA", "pvalue", "Cancer Deaths", 
                              "Cancer Deaths with Top", "Cancer Deaths with Down", 
                              "Mean Tumor Top", "Mean Tumor Down", "Mean Normal")
tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)
ngenes <- nrow(as.matrix(rownames(dataNormal)))
#for (i in 1:nrow(dataCancer)) {
for (i in 1:nrow(rna_vm)) { 
  cat(paste((ngenes - i), ".", sep = " "))
  mRNAselected <- as.matrix(rownames(dataCancer))[i]
  tabSurv_Matrix[i, "mRNA"] <- mRNAselected
  mRNAselected_values <- dataCancer[rownames(dataCancer) == mRNAselected, ]
  mRNAselected_values_normal <- dataNormal[rownames(dataNormal) == mRNAselected, ]
  
  mRNAselected_values_ordered <- sort(mRNAselected_values, decreasing = TRUE)
  numberOfSamples <- nrow(as.matrix(mRNAselected_values_ordered))  
  
  samples_top_mRNA_selected <- names(mRNAselected_values_ordered[mRNAselected_values_ordered > median(mRNAselected_values_ordered)])
  samples_down_mRNA_selected <- names(mRNAselected_values_ordered[mRNAselected_values_ordered < median(mRNAselected_values_ordered)])
  
  cfu_onlyTOP <- cfu[cfu[, "bcr_patient_barcode"] %in% samples_top_mRNA_selected, ]
  cfu_onlyDOWN <- cfu[cfu[, "bcr_patient_barcode"] %in% samples_down_mRNA_selected, ]
  
  cfu_ordered <- NULL
  cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
  #cfu <- cfu_ordered
  ttime <- as.numeric(cfu_ordered[, "days_to_death"])
  status <- cfu_ordered[, "vital_status"]
  #sum(status <- ttime > 0)
  deads_complete <- dim(cfu_ordered[cfu_ordered$vital_status > 0,])[1]
  deads_top <- dim(cfu_onlyTOP[cfu_onlyTOP$vital_status > 0,])[1]
  if (dim(cfu_onlyDOWN)[1] >= 1) {
    deads_down <- dim(cfu_onlyDOWN[cfu_onlyDOWN$vital_status > 0,])[1]
  }
  else {
    deads_down <- 0
  }
  
  tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
  tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
  tabSurv_Matrix[i, "Cancer Deaths with Down"] <- deads_down
  tabSurv_Matrix[i, "Mean Normal"] <- mean(mRNAselected_values_normal)
  dataCancer_onlyTop_sample <- dataCancer[, samples_top_mRNA_selected]
  dataCancer_onlyTop_sample_mRNASelected <- dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == mRNAselected, ]
  dataCancer_onlyDown_sample <- dataCancer[, samples_down_mRNA_selected]
  dataCancer_onlyDown_sample_mRNASelected <- dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == mRNAselected, ]
  tabSurv_Matrix[i, "Mean Tumor Top"] <- mean(dataCancer_onlyTop_sample_mRNASelected)
  tabSurv_Matrix[i, "Mean Tumor Down"] <- mean(dataCancer_onlyDown_sample_mRNASelected)
  ttime <- Surv(ttime, status)
  rownames(ttime) <- rownames(cfu_ordered)
  length(ttime)
  legendHigh <- paste(mRNAselected, "High")
  legendLow <- paste(mRNAselected, "Low")
  tabSurv <- survdiff(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN))))
  tabSurv_chis <- unlist(tabSurv)$chisq
  tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
  tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue
}

rownames(tabSurv_Matrix) <- tabSurv_Matrix$mRNA
Sig_KM_table <- tabSurv_Matrix[tabSurv_Matrix$pvalue < 0.05,]
################################################################################
save.image(file = "Sig_KM_table.RData")