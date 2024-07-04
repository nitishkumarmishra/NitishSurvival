##########################################################

#### This script takes mQTL results and identify overlap of DNA methylation sites with regulatory regions from publically available data from ENCODE and epigenomics roadmap project.
#### These overlaps were then compared to randomly selected sets of DNA methylation sites

#### Output is table with enrichment statistics


load("/mnt/data1/reference_files/Roadmap/FetalBrain/450KProbesinFetalBrainRegulatoryRegions.rdata") ## annoMat matrix indicating overlap of DNA methylation sites in regions idenitified from chip-Ssq data from fetal brain in the ROadmap project
load("/mnt/data1/reference_files/SliekerAnnotation/Slieker_TFBindingSites-s17_AsListPerProbe.rdata")

checkTFAnnotation<-function(probes){
	return(probes.TF[match(probes, names(probes.TF))])

}

load("/mnt/data1/reference_files/SliekerAnnotation/Slieker_DHSSites-s18.rdata")
probes.DHS<-DHSsites$name

checkDHSAnnotation<-function(probes){
	DHS.indicator<-rep(0, length(probes))
	DHS.indicator[which(probes %in% probes.DHS)]<-1
	return(DHS.indicator)

}

load("/mnt/data1/Eilis/Projects/References/ENCODE/450KProbesinBrainDHSSites.rdata")

checkDHSCerebellum<-function(probes){
	DHS.indicator<-probe.cer.dhs[probes]
	return(DHS.indicator)

}

checkDHSFrontalCortex<-function(probes){
	DHS.indicator<-probe.front.dhs[probes]
	return(DHS.indicator)

}

checkDHSCerebrumFrontal<-function(probes){
	DHS.indicator<-probe.fc.dhs[probes]
	return(DHS.indicator)

}

### Illumina annotation for all DNAm sites
load("") ## probeAnnot

### take background list as all probes in fetal dataset
setwd("")
expression_file_name = "" ## input file to MatrixEQTL 
meth<-read.table(expression_file_name, header = TRUE)

probeAnnot<-probeAnnot[rownames(meth),]

probeAnnot<-cbind(probeAnnot, checkTFAnnotation(probeAnnot$TargetID), checkDHSAnnotation(probeAnnot$TargetID), checkDHSCerebellum(probeAnnot$TargetID), checkDHSCerebrumFrontal(probeAnnot$TargetID), checkDHSFrontalCortex(probeAnnot$TargetID))
colnames(probeAnnot)[(ncol(probeAnnot)-4):ncol(probeAnnot)]<-c("TranscriptionFactorBindingSiteAnnotation", "DHSAnnotation", "CerebellumDHSAnnotation","CerebrumFrontalDHSAnnotation","FrontalCortexDHSAnnotation")
probeAnnot$TranscriptionFactorBindingSiteAnnotation<-as.character(probeAnnot$TranscriptionFactorBindingSiteAnnotation)

### load mQTL results
test<-read.csv("", stringsAsFactors = FALSE)
test<-unique(test$NAME)
test.info<-probeAnnot[test,]
### randomly select nProbes
nProbes<-nrow(test.info)

### separate out transcription factors
allTF<-unique(unlist(strsplit(probeAnnot$TranscriptionFactorBindingSiteAnnotation, "\\|")))
allTF<-allTF[!is.na(allTF)]
toDo<-c(colnames(annoMat), "DHSAnnotation", "CerebellumDHSAnnotation","CerebrumFrontalDHSAnnotation","FrontalCortexDHSAnnotation","TranscriptionFactorBindingSiteAnnotation",  allTF)

### count true overlap for DNAm sites part of significant mQTLs
true<-c(colSums(annoMat[test,]), sum(test.info$DHSAnnotation), sum(test.info$CerebellumDHSAnnotation),sum(test.info$CerebrumFrontalDHSAnnotation), sum(test.info$FrontalCortexDHSAnnotation),length(which(test.info$TranscriptionFactorBindingSiteAnnotation != "")))
for(each in allTF){
	true<-append(true, length(which(unlist(strsplit(test.info$TranscriptionFactorBindingSiteAnnotation, "\\|")) == each)))
}
names(true)<-toDo

### count overlap with all background DNAm sites
total<-c(colSums(annoMat[rownames(probeAnnot),]), sum(probeAnnot$DHSAnnotation), sum(probeAnnot$CerebellumDHSAnnotation),sum(probeAnnot$CerebrumFrontalDHSAnnotation), sum(probeAnnot$FrontalCortexDHSAnnotation),length(which(probeAnnot$TranscriptionFactorBindingSiteAnnotation != "")))
for(each in allTF){
	total<-append(total, length(which(unlist(strsplit(probeAnnot$TranscriptionFactorBindingSiteAnnotation, "\\|")) == each)))
}
names(total)<-toDo


### for each set of regulatory regions, perform Fisher's test for enrichment
out.tmp<-matrix(data = NA, nrow = length(toDo), ncol = 6)
colnames(out.tmp)<-c("Total", "Overlap", "FC", "95%Lower", "95%Upper", "P")
rownames(out.tmp)<-toDo
out.tmp[,1]<-total
out.tmp[,2]<-true
for(i in 1:length(toDo)){
	a<-true[i]
	b<-nProbes-a
	d<-total[i]-a
	e<-nrow(probeAnnot)-d-b-a
	tmp<-fisher.test(matrix(c(a,b,d,e), nrow = 2))
	out.tmp[i,3]<-tmp$estimate
	out.tmp[i,4]<-tmp$conf.int[1]	
	out.tmp[i,5]<-tmp$conf.int[2]
	out.tmp[i,6]<-tmp$p.value		
}

write.csv(out.tmp, "")

