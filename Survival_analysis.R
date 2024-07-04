library(limma)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(edgeR)
clin <- GDCquery_clinic("TCGA-CHOL", type = "clinical", save.csv = FALSE)

rna <- read.table('CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t', stringsAsFactors=FALSE, check.names = FALSE)
rna <- rna[-c(1),which(rna[1,]=="raw_count")]
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

#n_index <- which(substr(colnames(rna),14,14) == '1')
#t_index <- which(substr(colnames(rna),14,14) == '0')
#vm <- function(x){
#cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
#d <- model.matrix(~1+cond)
#x <- t(apply(x,1,as.numeric))
#ex <- voom(x,d,plot=F)
#return(ex$E)
#}
cancerID <- grep("01A", colnames(rna))
normalID <- grep("11A", colnames(rna))
rna <- cbind(rna[,cancerID], rna[,normalID])

design <- model.matrix(~0 + factor(c(rep(2, length(cancerID)), rep(1, length(normalID)))))
colnames(design) <- c("Normal", "Tumor")
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
factors <- factor(c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(rna, group=factors)
cnf <- calcNormFactors(cnts.dgelist, method = "TMM")
#cnf <- calcNormFactors(cnts.liver, method = "TMM") ## TMM normalization 
#v <- voom(cnts.liver, design, plot = TRUE, lib.size=colSums(cnts.liver) * cnf) ##Transform count data to log2-counts per million (logCPM)
v <- voom(cnf, design, plot = FALSE)

rna_vm  <- v$E
#colnames(rna) <- substr(colnames(rna), 1, 12)

scal <- function(x,y){
mean_n <- rowMeans(y)  # mean of normal
sd_n <- apply(y,1,sd)  # SD of normal
# z score as (value - mean normal)/SD normal
res <- matrix(nrow=nrow(x), ncol=ncol(x))
colnames(res) <- colnames(x)
rownames(res) <- rownames(x)
for(i in 1:dim(x)[1]){
for(j in 1:dim(x)[2]){
res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
}
}
return(res)
}
cancerID <- grep("01A", colnames(rna_vm))
normalID <- grep("11A", colnames(rna_vm))
z_rna <- scal(rna_vm[,cancerID],rna_vm[,normalID])
#rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,'\\|'))[[1]])
rm(rna_vm)


colnames(z_rna) <- substr(colnames(z_rna), 1, 12)
event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))
dim(event_rna)

data <- clin
rownames(data) <- data$submitter_id
notDead <- is.na(data$days_to_death)
if (any(notDead == TRUE)) {
  data[notDead, "days_to_death"] <- data[notDead, "days_to_last_follow_up"]
}
data$s <- grepl("dead|deceased", data$vital_status, ignore.case = TRUE)
data$type <- as.factor(data[, c("gender")])
data$vital_status <- ifelse(data$vital_status == "alive", 0, 1)
data <- data[, c("days_to_death", "s", "type", "vital_status")]

########################################
data1 <- data[rownames(data)%in%colnames(event_rna),]
t_event_rna <- t(event_rna)
t_event_rna <- t_event_rna[rownames(data1),]

ind_gene <- t_event_rna[,which(colnames(t_event_rna) == 'FOXC1')]


data1$type <- ifelse( ind_gene == "1", "up", "down")


f.m <- formula(Surv(as.numeric(data1$days_to_death), event = data1$vital_status) ~data1$type)
fit <- do.call(survfit, list(formula = f.m, data = data1))
label.add.n <- function(x) {
paste0(x, " (n = ", nrow(data1[data1[, "type"] == x, ]),
")")
}

#not necessary
if (is.null(labels)) {
d <- survminer::surv_summary(fit, data = data1)
order <- unname(sapply(levels(d$strata), function(x) unlist(strsplit(x,"="))[2]))
labels <- sapply(order, label.add.n)
}

#suppressWarnings({
#surv <- ggsurvplot(fit, risk.table = TRUE, pval = TRUE)})

suppressWarnings({
  surv <- ggsurvplot(fit, risk.table = TRUE, xlim = NULL, main="Kaplan-Meier Overall Survival Curves", ylab = "Probability of survival",
                     labels = NULL, legend = "top",legend.title="Expression", palette=c("blue4", "red4"), 
                     xlab = "Time since diagnosis (days)", pval = TRUE, pval.size = 5, pval.coord = c(1,0))})

ggsave(surv$plot, filename = "survival.pdf", width = 12, height = 8, dpi = 300)

##########################
save.image("Survival.RData")

### Top gene from TCGABiolinks
## In clin we replace "alive" < "Alive" and "dead" << "Dead"
Genelist <- rownames(rna_vm)
results <- TCGAanalyze_SurvivalKM(clin, dataGE = rna_vm, Genelist = Genelist, ThreshTop = 0.75, ThreshDown = 0.25, Survresult = TRUE, p.cut = 0.05, group1 = group1, group2 = group2)

Genelist1 <- head(rownames(results),15)
results1 <- plot_surv(clin, dataGE = rna, Genelist = Genelist1, ThreshTop = 0.6, ThreshDown = 0.4, Survresult = TRUE, p.cut = 0.01, group1 = group1, group2 = group2)

save.image("Survival_TCGAbiolinks.RData")
