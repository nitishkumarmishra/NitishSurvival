library(TCGAbiolinks)
library(dplyr)
library(DT) # for datatable function
library(SummarizedExperiment) # for assay function

setwd("F:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/")
#The other fields (data.category, data.type, workflow.type, platform, file.type)
datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454",col_types = readr::cols()),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

dataClin <- GDCquery_clinic(project = "TARGET-OS","clinical") 
dataclin.select <- subset(dataClin,select=c(vital_status,  days_to_death, age_at_diagnosis, days_to_last_follow_up, primary_diagnosis, ethnicity, gender, site_of_resection_or_biopsy, tissue_or_organ_of_origin, submitter_id))
##Download HTSeq-Counts data
query.exp.hg38 <- GDCquery(project = "TARGET-OS", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
datatable(getResults(query.exp.hg38, cols = c("data_type","cases")),
          filter = 'top', options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), rownames = FALSE)
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "exp_HTseq_Counts.rda")
OS_HTSeq_Counts <- assay(expdat,"HTSeq - Counts")

##Download HTSeq-FPKM data
query.exp.hg38 <- GDCquery(project = "TARGET-OS", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "exp_HTseq_FPKM.rda")
OS_HTSeq_FPKM <- assay(expdat,"HTSeq - FPKM")

##Download HTSeq-FPKM-UQ data
query.exp.hg38 <- GDCquery(project = "TARGET-OS", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "exp_HTSeq_UQ.rda")
OS_HTSeq_FPKM_UQ <- assay(expdat,"HTSeq - FPKM-UQ")

rm(expdat, query.exp.hg38)


OS_Purity <- read.csv("TARGET_OS_ClinicalData_Discovery_and_Validation_Percent_Tumor_Supplement_20180624.txt", header = TRUE, sep = "\t")


####################################################
save.image("TARGET-OS-RNAseq.RData")
####################################################

load("F:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/OsteosarcomaData.RData")
## Common patients in DNA methylation and Gene Expression
intersect(substr(colnames(OS_HTSeq_Counts), 1, 19), substr(colnames(myNorm), 1, 19))


####################################################
save.image("TARGET-OS-RNAseq-Clin.RData")
####################################################