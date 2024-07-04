##########################################################

#### This script takes mQTL results that have been clumped to establish quasi independent SNP set, calculates overlap with GWAS results at multiple thresholds, and performs simulation matched for allele frequency to test if overlap is more than expected by chance.
#### Simulations are run in parallel across multiple cores

#### Output is table with enrichment statistics

setwd("")
## Two output files from clumping as more stringent threshold used for MHC region.
qtl_results_1<-read.table("",header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
qtl_results_2<-read.table("",header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

qtl_results_1<-qtl_results_1[which(!is.na(qtl_results_1$CHR) == TRUE),]
qtl_results_2<-qtl_results_2[which(!is.na(qtl_results_2$CHR) == TRUE),]
qtl_snps<-unique(c(qtl_results_1$SNP, qtl_results_2$SNP)) ### independent SNP set

### load GWAS summary statistics
gwas_results<-read.table("", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
### load all SNPs passed QC
Map<-read.table("")
### idenitfy overlap
gwas_results<-gwas_results[match(Map[,1], gwas_results$QTL_SNP),]

### annotate with allele frequency
Freq<-read.table("", header = TRUE)
gwas_results<-cbind(gwas_results, Freq$MAF[match(gwas_results$QTL_SNP, Freq$SNP)])
colnames(gwas_results)[ncol(gwas_results)]<-"Data_Freq"

### filter to SNPs with GWAS results
gwas_results<-gwas_results[which(is.na(gwas_results$P) == FALSE),]

### pull out GWAS p value and SNP frequency for mQTL SNP set
qtl_snps<-gwas_results[match(intersect(qtl_snps, gwas_results$QTL_SNP), gwas_results$QTL_SNP),]

### match simulations on MAF:
### problems with exact pattern matching so added tolerance
freq_bins<-vector(length = length(seq(from = 0.02, to = 1, by = 0.02)))
thres<-seq(from = 0.02, to = 1, by = 0.02)
a<-NULL
for(i in 1:length(thres)){
	freq_bins[i]<-length(which(qtl_snps$Data_Freq <= thres[i]+ 1e-10 & qtl_snps$Data_Freq > (thres[i]-0.02)))
	a<-append(a, which(qtl_snps$Data_Freq <= thres[i] & qtl_snps$Data_Freq > (thres[i]-0.02)))

}

## create function to do one simulation
permutation<-function(pThres, thres, freq_bins, gwas_results){
tmp<-vector(length = length(pThres))
sim_snps<-c()
	for(i in 1:length(thres)){
		
		sim_snps<-append(sim_snps, sample(gwas_results$QTL_SNP[which(gwas_results$Data_Freq < thres[i]+ 1e-10 & gwas_results$Data_Freq >= (thres[i]-0.02))],freq_bins[i]))
	}
	gwas_sub<-gwas_results[match(sim_snps, gwas_results$QTL_SNP),]
	for(k in 1:length(pThres)){
		
		tmp[k]<-length(which(gwas_sub$P < pThres[k]))
	}
	return(tmp)
}


pThres<-c(5e-5, 5e-6,5e-7,5e-8)

### establish parallel environment
library(doParallel)
cl<-makeCluster(50)
registerDoParallel(cl)

nPerms<-100000000

perms.tmp<-foreach(i=1:nPerms, .combine=cbind) %dopar%{
		permutation(pThres, thres, freq_bins, gwas_results)
	}

true<-vector(length = 7)
emp<-vector(length = 7)
all<-vector(length = 7)
pdf("")
for(k in 1:length(pThres)){
	
	all[k]<-length(which(gwas_results$P < pThres[k]))
	true[k]<-length(which(qtl_snps$P < pThres[k]))
	emp[k]<-length(which(perms.tmp[k,] >= true[k]))/ncol(perms.tmp)
	hist(perms.tmp[k,], main = paste("GWAS P < ", pThres[k], sep = ""), xlim = c(min(c(perms.tmp[k,], true[k])), max(c(perms.tmp[k,], true[k]))))
	abline(v = true[k], col = "red")
}
qtl.rate<-true/nrow(qtl_snps)
all.rate<-all/nrow(gwas_results)
write.csv(rbind(pThres, true,qtl.rate, all.rate, emp), "")

dev.off()