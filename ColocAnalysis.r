##########################################################

#### This script takes mQTL results and GWAS results from PGC and performs colocalisation analysis

#### Output is table with posterior probabilities for each probe within each region.

library(coloc)
library(data.table)


removeAlleles<-function(text){
	tmp<-unlist(strsplit(unlist(strsplit(gsub("X", "", text), "\\.")), "_"))
	if(length(tmp) == 3){
		return(tmp)
	} else {
		return(c(tmp[1:2], paste(tmp[3:length(tmp)], collapse = "_")))
	}
}


### load GWAS results both SNP summary statistics and regions
gwas<-fread("daner_PGC_SCZ52_0513a.hq2")
regions<-read.table("scz2.anneal.108.txt", header = TRUE, stringsAsFactors = FALSE)

### load all DNAm sites that passed QC and probe annotation
setwd("")
meth<-read.table("", header = TRUE, row.names = 1)
load("")
probeAnnot<-probeAnnot[rownames(meth),]


### loop through chromosomes NOTE not all have a signifiacnt region
output<-NULL
for(chr in c(1:12, 14:20,22)){
	mQTL<-read.table(paste("", chr, "", sep = ""), stringsAsFactors = FALSE, header = TRUE) ## load cis mQTL results for each chromosome
	regions.sub<-regions[which(regions$hg19chrc == paste("chr", chr, sep = "")),] ## identify regions on this chr
	locs<-matrix(unlist(lapply(mQTL$SNP, removeAlleles)), ncol = 3, byrow = TRUE)
	mQTL<-cbind(mQTL, as.character(paste(locs[,1], locs[,2], sep = ":"))) ## match chromosome by chr:bp format
	
	freq<-read.table(paste("", chr, "", sep = ""), header = TRUE, stringsAsFactors = FALSE) ## load allele frequencies
	for(i in 1:nrow(regions.sub)){
		gwas.sub<-gwas[which(gwas$CHR == chr & gwas$BP <= regions.sub$anneal2[i]+500000 & gwas$BP >= regions.sub$anneal1[i]-500000),]
		dataset2<-list(beta= log(gwas.sub$OR), varbeta= (gwas.sub$SE^2), type = "cc", snp = paste(gwas.sub$CHR, gwas.sub$BP, sep = ":"), MAF = gwas.sub$FRQ_U_46839, N=82315) ### for each region SCZ GWAS results stay the same
		probes<-as.character(probeAnnot[which(probeAnnot$CHR == chr & probeAnnot$MAPINFO <= regions.sub$anneal2[i]+500000 & probeAnnot$MAPINFO >= regions.sub$anneal1[i]-500000),"NAME"])
			
		probes<-intersect(probes, rownames(meth))
		for(each in probes){
			mQTL.sub<-mQTL[which(mQTL$gene == each),]
			if(nrow(mQTL.sub) > 1){
				dataset1<-list(beta=mQTL.sub$beta,varbeta=(mQTL.sub$beta/mQTL.sub$t.stat)^2,type = "quant",snp = mQTL.sub[,7], MAF = freq$MAF[match(mQTL.sub[,7], freq$SNP)], N=166) ### mQTL results
				my.res<-coloc.abf(dataset1, dataset2)
				output<-rbind(output, c("region"=regions.sub$anneal.rank[i], "trait2"=each, my.res$summary))
			}
		}
	}
}

write.csv(output, "")

### recreate Fig 1 from http://hmg.oxfordjournals.org/content/early/2015/04/03/hmg.ddv077.full

d0<-density(as.numeric(as.character(output[,4])), from = 0, to = 1, bw = 0.01)
d1<-density(as.numeric(as.character(output[,5]))+as.numeric(as.character(output[,6])), from = 0, to = 1, bw = 0.01)
d2<-density(as.numeric(as.character(output[,7]))+as.numeric(as.character(output[,8])), from = 0, to = 1, bw = 0.01)

pdf("")
plot(d0)
lines(d1, col = "blue")
lines(d2, col = "red")
legend("topright", lty = 1, col = c("black", "blue", "red"), c("PP0", "PP1+PP2", "PP3+PP4"))
dev.off()