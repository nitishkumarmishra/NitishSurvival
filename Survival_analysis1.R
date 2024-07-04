library(limma)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(edgeR)
library(matrixStats)
clin <- GDCquery_clinic("TCGA-CHOL", type = "clinical", save.csv = FALSE)

rna <- read.table('CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t', stringsAsFactors=FALSE, check.names = FALSE)
rna <- rna[-1,]
rownames(rna) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(rna))
#rna <- round(data.matrix(rna),0)
rna <- rna[-grep("\\?", rownames(rna)),] ## Remove genes which HGNC name is known i.e. ?
rownames(rna) <- sapply(strsplit(rownames(rna),"\\|"),'[[',1)
## Remove sample which have zero in more than 25% samples
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.25)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]
#################
n_index <- which(substr(colnames(rna),14,14) == '1')
t_index <- which(substr(colnames(rna),14,14) == '0')

### take log2
lg2 <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  res <- log2((x+1))
  return(res)
}
rna_lg <- lg2(rna)
colnames(rna_lg) <- substr(colnames(rna),1,12)
#n_index <- which(substr(colnames(rna),14,14) == '1')
#t_index <- which(substr(colnames(rna),14,14) == '0')
#vm <- function(x){
#cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
#d <- model.matrix(~1+cond)
#x <- t(apply(x,1,as.numeric))
#ex <- voom(x,d,plot=F)
#return(ex$E)
#}
#rna_vm  <- vm(rna)
#colnames(rna_vm) <- gsub('\\.','-',substr(colnames(rna),1,12))
#colnames(rna) <- substr(colnames(rna), 1, 12)

scal <- function(x,y){
  mean_n <- rowMeans2(y)  # mean of normal
  #sd_n <- apply(y,1,sd)  # SD of normal
  sd_n <- rowSds(y)
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
#cancerID <- grep("01A", colnames(rna_vm))
#normalID <- grep("11A", colnames(rna_vm))
z_rna <- scal(rna_lg[,t_index],rna_lg[,n_index])
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])
#rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,'\\|'))[[1]])
rm(rna_vm)

event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))
event_rna <- t(apply(z_rna, 1, function(x) ifelse(x > 1.5,2,ifelse(x < -1.5,1,0)))) 
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

ind_gene <- t_event_rna[,which(colnames(t_event_rna) == 'MUC13')]


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


