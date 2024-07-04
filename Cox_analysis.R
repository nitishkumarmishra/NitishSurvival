library(survcomp)
####################################################
#### These lines are just example ##################
data(breastCancerData)
data(breastCancerData)
print(vdx7g)
print(head(phenoData(vdx7g)@data ))
dd <- phenoData(vdx7g)@data [, c("t.dmfs", "e.dmfs")]
colnames(dd) <- c("time", "event") ## Survival data
ge <- t(exprs(vdx7g)) ### expression data
dd <- cbind(dd, ge)
mm <- coxph(Surv(time, event) ~., data = data.frame(dd)) ##Multiple covariates analysis
summary(mm)
print(summary(mm))
tmp <- summary(mm)
tmp$coefficients

############################################
############  Real code from here ##########

clin <- GDCquery_clinic("TCGA-CHOL", type = "clinical", save.csv = FALSE)

rna <- read.table('CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t', stringsAsFactors=FALSE, check.names = FALSE)
rna <- rna[-1,]
#rna <- rna[-c(1),which(rna[1,]=="raw_count")]
rownames(rna) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(rna))
rna <- round(data.matrix(rna),0)
rna <- rna[-grep("\\?", rownames(rna)),] ## Remove genes which HGNC name is known i.e. ?
rownames(rna) <- sapply(strsplit(rownames(rna),"\\|"),'[[',1)
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
# rna_vm <- log2(rna+1)
Genelist <- rownames(rna_vm)
Genelist <- intersect(rownames(rna_vm), Genelist)
dataCancer <- rna_vm[Genelist, group2]
dataNormal <- rna_vm[Genelist, group1]
colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)

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

dd <- cfu[,c("days_to_death", "vital_status")]
colnames(dd) <- c("time", "event")
ge <- t(dataCancer)
dd <- cbind(dd, ge)
mm <- coxph(Surv(time, event) ~ ., data = data.frame(dd)) ##Multiple covariates analysis
summary(mm)
