##########################################################

#### This script takes mQTL results and uses illumina annotation to generate a gene list for pathway analysis, it also creates a variable with the number of probes annotated to each each to control for potential bias. 
#### The pathway analysis is based on a logistic regression model across all genes, it then transform p values to one-tailed for over representation.
#### The initial pathway analysis is set up to be run in parallel
#### After performing the regression for each GO term, it takes the set of significant terms, and identifies were gene membership overlap with more significant terms is driving the association.

#### Output is two tables; 1 for each GO term with enrichment p value; 2 with results of the merging of overlapping pathways


### load background lists (all DNA methylation sites that passed QC)
load("") ## probeAnnot

setwd("")
expression_file_name = ""
meth<-read.table(expression_file_name, header = TRUE)
probeAnnot<-probeAnnot[rownames(meth),]

bg_probes<-rownames(probeAnnot)
bg_genes<-unique(unlist(strsplit(as.character(probeAnnot$UCSC_REFGENE_NAME), ";")))


##### load test lists
test<-read.csv("", stringsAsFactors = FALSE)
test_probes<-unique(test$NAME)

### load gene-pathway map 
gene_go<-read.csv("", stringsAsFactors = FALSE,)
# filter out genes not annotated to any pathway
gene_go<-gene_go[which(gene_go[,2] != ""),]

bg_gene_go<-gene_go[match(intersect(bg_genes, gene_go[,1]), gene_go[,1]),]
terms<-unique(unlist(strsplit(as.character(bg_gene_go[,2]), "\\|")))
names<-read.csv("", stringsAsFactors = FALSE) ### file matching GO ids to pathway names

### function to run over-representation analysis for each GO term.
apply_tests<-function(each){

		## logistic regression approach
		pathway<-vector(length = nrow(bg_gene_go))
		pathway[grep(each, bg_gene_go[,2])]<-1
		model<-glm(pathway ~ gene_test + gene_size)
		p1<-summary(model)$coefficients[c("gene_test"),c(4)]
		p2<-summary(model)$coefficients[c("gene_size"),c(4)]
		or<-exp(summary(model)$coefficients[c("gene_test"),c(1)])
		beta<-summary(model)$coefficients[c("gene_size"),c(1)]
		if(or > 1){
			p1<-p1/2
		} else {
			p1<-1-(p1/2)
		}
		return(unlist(c(each, names[match(each, names[,1]),2:3], length(which(pathway == 1)), length(which(pathway == 1 & gene_test)), p1,or,p2,beta, paste(bg_gene_go[which(pathway == 1 & gene_test == 1),1], collapse = "|"))))

		}

## create variable to control for nProbes annotated to each gene
probeAnnot$UCSC_REFGENE_NAME<-unlist(lapply(lapply(strsplit(as.character(probeAnnot$UCSC_REFGENE_NAME), ";"), unique), paste , collapse = ";"))
gene_size<-table(unlist(strsplit(as.character(probeAnnot$UCSC_REFGENE_NAME), ";")))[bg_gene_go[,1]]


## set up parallel environment
library(doParallel)
cl<-makeCluster(30)
registerDoParallel(cl)


test_genes<-probeAnnot[test_probes,]
test_genes<-unique(unlist(strsplit(as.character(test_genes$UCSC_REFGENE_NAME), ";")))

test_gene_go<-gene_go[match(intersect(test_genes, gene_go[,1]), gene_go[,1]),]

gene_test<-vector(length = nrow(bg_gene_go))
gene_test[match(test_genes, bg_gene_go[,1])]<-1

r<-foreach(i=1:length(terms), .combine=rbind, .export = c("gene_size")) %dopar%{
	apply_tests(terms[i])
	}
colnames(r)<-c("ID", "Name", "Type", "nGenesinPathway", "nTestListinPathway", "P:GenesinTestList", "OR", "P:GeneSize", "Beta:GeneSize", "GenesinTestListAndPathway")

### filter to terms with between 10 and 2000 genes
r<-r[which(as.numeric(r[,4]) < 2000 & as.numeric(r[,4]) > 9),]
	
write.csv(r, "")

### identify significant terms
r<-r[order(as.numeric(r[,6])),]
r<-r[which(as.numeric(r[,6]) < 0.05 & as.numeric(r[,7]) > 0),]
r.all<-r

### iterative procedure to merge terms where significance is due to overlapping gene membership

output<-c()
while(!is.null(r)){
	
	if( class(r) != "character"){
		
	### for all terms repeat analysis controlling for most significant terms
	best_term<-vector(length =  nrow(bg_gene_go))
	best_term[grep(r[1,1], bg_gene_go[,2])]<-1
	merge.id<-c()
	merge.name<-c()
	remove<-c()
		
	for(j in 2:nrow(r)){
		
		pathway<-vector(length = nrow(bg_gene_go))
		pathway[grep(r[j,1], bg_gene_go[,2])]<-1
		model<-glm(pathway ~ gene_test + gene_size + best_term)
		if(summary(model)$coefficients["gene_test", "Pr(>|t|)"] > 0.05){
			merge.id<-append(merge.id, r[j,1])
			merge.name<-append(merge.name, r[j,2])
			remove<-append(remove, j)
		}
		
	}
	merge.id<-paste(unlist(merge.id), collapse = "|")
	merge.name<-paste(unlist(merge.name), collapse = "|")
	
	output<-rbind(output, c(r[1,], merge.id, merge.name))
		
	r<-r[-c(1, remove),]
	} else {
		
		output<-rbind(output, c(r, "", ""))
		r<-NULL
	}
}


write.csv(output, "")

### create output table with p values for each merged term
r<-r.all

output_2<-NULL
for(i in 1:nrow(output)){
	if(output[i,11] != ""){
		tmp.1<-r[match(unlist(strsplit(output[i,11], "\\|")), r$ID),]
		tmp.2<-matrix(data = NA, nrow = (length(unlist(strsplit(output[i,11], "\\|")))-1), ncol = ncol(output))
		colnames(tmp.2)<-colnames(output)
		tmp.2<-rbind(output[i,], tmp.2)
		tmp.1<-tmp.1[order(tmp.1$P.GenesinTestList),]
		output_2<-rbind(output_2, cbind(tmp.1, tmp.2))
	} else {
		tmp.2<-output[i,]
		tmp.1<-rep(NA, ncol(r))
		names(tmp.1)<-colnames(r)
		output_2<-rbind(output_2, c(tmp.1, tmp.2))
		}

}

write.csv(output_2, "")
