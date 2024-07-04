##########################################################

#### This script takes genotypes, DNA methylation values and covariate data for all genome wide mQTL calculations.
#### Genotypes were filtered using PLINK for MAF, HWE, and missingness.
#### Files should be formatted as described in http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html

#### Output is text file of all mQTL above a certain threshold

setwd()
library(MatrixEQTL)
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR
SNP_file_name = ""
expression_file_name = ""
## if no covariates set to  character()
covariates_file_name = ""
output_file_name = ""

## threshold of results to save
pvOutputThreshold = 0.0001;
## set error covarriance. NOTE rarely used instead set to numeric()
errorCovariance = numeric();
	 
## load SNP data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

## load methylation data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name);

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name);

### run mqtls
me = Matrix_eQTL_engine(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	output_file_name = output_file_name,
	pvOutputThreshold = pvOutputThreshold,
	useModel = useModel, 
	errorCovariance = errorCovariance, 
	verbose = TRUE,
	pvalue.hist = "qqplot",
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = FALSE)	
