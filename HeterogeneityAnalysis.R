##########################################################

#### This script takes mQTL results and GWAS results from PGC and performs colocalisation analysis

#### Output is table with posterior probabilities for each probe within each region.


library(lme4)

### load genotype, methylation and covariates data for all four datasets
setwd("")

SNP_file_name = ""
expression_file_name = ""
meth.pfc<-read.table(expression_file_name, header = TRUE, row.names = 1)
meth.pfc<-as.matrix(meth.pfc)
geno.pfc<-read.table(SNP_file_name, header = TRUE, row.names = 1)
geno.pfc<-as.matrix(geno.pfc)
cov.pfc<-read.table("", header = TRUE, row.names = 1)
cov.pfc<-as.matrix(cov.pfc)

SNP_file_name = ""
expression_file_name = ""
meth.cer<-read.table(expression_file_name, header = TRUE, row.names = 1)
meth.cer<-as.matrix(meth.cer)
geno.cer<-read.table(SNP_file_name, header = TRUE, row.names = 1)
geno.cer<-as.matrix(geno.cer)
cov.cer<-read.table("", header = TRUE, row.names = 1)
cov.cer<-as.matrix(cov.cer)


SNP_file_name = ""
expression_file_name = ""
meth.str<-read.table(expression_file_name, header = TRUE, row.names = 1)
meth.str<-as.matrix(meth.str)
geno.str<-read.table(SNP_file_name, header = TRUE, row.names = 1)
geno.str<-as.matrix(geno.str)
cov.str<-read.table("", header = TRUE, row.names = 1)
cov.str<-as.matrix(cov.str)


geno.fetal<-read.table("", stringsAsFactors = FALSE)
geno.fetal<-as.matrix(geno.fetal)
meth.fetal<-read.table("", header = TRUE, row.names = 1)
meth.fetal<-as.matrix(meth.fetal)
cov.fetal<-read.table("", header = TRUE, row.names = 1)
cov.fetal<-as.matrix(cov.fetal)


### identify minor allele in fetal and adult datasets to check they are the same
adult.alleles<-unique(c(rownames(geno.pfc), rownames(geno.str), rownames(geno.cer)))
adult.alleles<-unlist(strsplit(unique(c(rownames(geno.pfc), rownames(geno.str), rownames(geno.cer))), "_"))[seq(from = 2, to = 2*length(adult.alleles), by = 2)]
names(adult.alleles)<-unlist(strsplit(unique(c(rownames(geno.pfc), rownames(geno.str), rownames(geno.cer))), "_"))[seq(from = 1, to = 2*length(adult.alleles), by = 2)]

fetal.alleles<-rownames(geno.fetal)
fetal.alleles<-unlist(strsplit(rownames(geno.fetal), "_"))[seq(from = 2, to = 2*length(fetal.alleles), by = 2)]
names(fetal.alleles)<-unlist(strsplit(rownames(geno.fetal), "_"))[seq(from = 1, to = 2*length(fetal.alleles), by = 2)]

rownames(geno.pfc)<-gsub("_.", "", rownames(geno.pfc))
rownames(geno.str)<-gsub("_.", "", rownames(geno.str))
rownames(geno.cer)<-gsub("_.", "", rownames(geno.cer))
rownames(geno.fetal)<-gsub("_.", "", rownames(geno.fetal))

## identify all mQTLs significant in fetal dataset but available to test in all
res<-read.csv("",  stringsAsFactors = FALSE)
res<-res[which(rowSums(cbind(is.na(res$BA9.P), is.na(res$STR.P), is.na(res$CER.P), is.na(res$Fetal.P))) == 0),] 
rownames(res)<-paste(res$SNP, res$Probe, sep = "|")

pairs<-cbind(res$SNP, res$Probe)
pairs<-unique(pairs)

### create variables for heterogeneity model
dev<-c(rep("Adult", ncol(meth.pfc)),rep("Adult", ncol(meth.str)),rep("Adult", ncol(meth.cer)), rep("Fetal", ncol(meth.fetal)))
br<-c(rep("PFC", ncol(meth.pfc)),rep("STR", ncol(meth.str)),rep("CER", ncol(meth.cer)), rep("Brain", ncol(meth.fetal)))
id<-c(colnames(cov.pfc), colnames(cov.str), colnames(cov.cer), colnames(cov.fetal))
cov<-cbind(cov.pfc[-3,], cov.str[-3,], cov.cer[-3,], cov.fetal)

### vectors to store output in 
hetP<-rep(NA, nrow(pairs))
hetP.nocov<-rep(NA, nrow(pairs))

for(i in 1:nrow(pairs)){
	meth<-c(meth.pfc[pairs[i,2],],meth.str[pairs[i,2],],meth.cer[pairs[i,2],], meth.fetal[pairs[i,2],])
	if(adult.alleles[pairs[i,1]] == fetal.alleles[pairs[i,1]]){
		geno<-c(geno.pfc[pairs[i,1],],geno.str[pairs[i,1],],geno.cer[pairs[i,1],], geno.fetal[pairs[i,1],])
	} else {
		geno<-c(geno.pfc[pairs[i,1],],geno.str[pairs[i,1],],geno.cer[pairs[i,1],], (2-geno.fetal[pairs[i,1],]))
	
	}
	nohet<-lmer(meth ~ geno + cov[1,] +cov[2,] + cov[3,] + cov[4,] + (1|id) + (1|br) + (1|dev), REML = FALSE)
	het<-lmer(meth ~ geno + cov[1,] +cov[2,] + cov[3,] + cov[4,] + (1|id) + (1|br) + (1|dev) + geno*dev, REML = FALSE)
	hetP[i]<-anova(het,nohet)[2,8]


	nohet<-lmer(meth ~ geno + (1|id)  + (1|br) + (1|dev), REML = FALSE)
	het<-lmer(meth ~ geno + (1|id) + (1|br) + (1|dev) + geno*dev, REML = FALSE)
	hetP.nocov[i]<-anova(het,nohet)[2,8]
	
}



hetP<-cbind(res[paste(pairs[,1], pairs[,2], sep = "|"),], hetP, hetP.nocov)
write.csv(hetP, "")
