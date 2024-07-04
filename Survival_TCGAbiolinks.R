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
#############################################################################

plot_surv <- function (clinical_patient, dataGE, Genelist, Survresult, ThreshTop = 0.67, 
                       ThreshDown = 0.33, p.cut = 0.05, group1, group2) 
{
  Genelist <- intersect(rownames(dataGE), Genelist)
  dataCancer <- dataGE[Genelist, group2]
  dataNormal <- dataGE[Genelist, group1]
  colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)
  cfu <- clinical_patient[clinical_patient[, "bcr_patient_barcode"] %in% 
                            substr(colnames(dataCancer), 1, 12), ]
  if ("days_to_last_followup" %in% colnames(cfu)) 
    colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <- "days_to_last_follow_up"
  cfu <- as.data.frame(subset(cfu, select = c("bcr_patient_barcode", 
                                              "days_to_death", "days_to_last_follow_up", "vital_status")))
  cfu[which(cfu$vital_status == "Alive"), "days_to_death"] <- "-Inf"
  cfu[which(cfu$vital_status == "Dead"), "days_to_last_follow_up"] <- "-Inf"
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  followUpLevel <- FALSE
  #Survresult <- FALSE
  tabSurv_Matrix <- matrix(0, nrow(as.matrix(rownames(dataNormal))), 
                           8)
  colnames(tabSurv_Matrix) <- c("mRNA", "pvalue", "Cancer Deaths", 
                                "Cancer Deaths with Top", "Cancer Deaths with Down", 
                                "Mean Tumor Top", "Mean Tumor Down", "Mean Normal")
  tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)
  cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
  cfu$days_to_last_follow_up <- as.numeric(as.character(cfu$days_to_last_follow_up))
  rownames(cfu) <- cfu[, "bcr_patient_barcode"]
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  cfu_complete <- cfu
  ngenes <- nrow(as.matrix(rownames(dataNormal)))
  for (i in 1:nrow(as.matrix(rownames(dataNormal)))) {
    cat(paste((ngenes - i), ".", sep = ""))
    mRNAselected <- as.matrix(rownames(dataNormal))[i]
    tabSurv_Matrix[i, "mRNA"] <- mRNAselected
    mRNAselected_values <- dataCancer[rownames(dataCancer) == 
                                        mRNAselected, ]
    mRNAselected_values_normal <- dataNormal[rownames(dataNormal) == 
                                               mRNAselected, ]
    mRNAselected_values_ordered <- sort(mRNAselected_values, 
                                        decreasing = TRUE)
    mRNAselected_values_ordered_top <- as.numeric(quantile(mRNAselected_values_ordered, 
                                                           ThreshTop)[1])
    mRNAselected_values_ordered_down <- as.numeric(quantile(mRNAselected_values_ordered, 
                                                            ThreshDown)[1])
    mRNAselected_values_newvector <- mRNAselected_values
    if (is.na(mRNAselected_values_ordered_top) != 1) {
      numberOfSamples <- nrow(as.matrix(mRNAselected_values_ordered))
      lastelementTOP <- round(numberOfSamples/3)
      firstelementDOWN <- numberOfSamples - lastelementTOP
      samples_top_mRNA_selected <- rownames(as.matrix(mRNAselected_values_ordered[1:(lastelementTOP - 
                                                                                       1)]))
      samples_down_mRNA_selected <- rownames(as.matrix(mRNAselected_values_ordered[(firstelementDOWN + 
                                                                                      1):numberOfSamples]))
      samples_UNCHANGED_mRNA_selected <- rownames(as.matrix(which((mRNAselected_values_newvector) > 
                                                                    mRNAselected_values_ordered_down & mRNAselected_values_newvector < 
                                                                    mRNAselected_values_ordered_top)))
      cfu_onlyTOP <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% 
                                    samples_top_mRNA_selected, ]
      cfu_onlyDOWN <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% 
                                     samples_down_mRNA_selected, ]
      cfu_onlyUNCHANGED <- cfu_complete[cfu_complete[, 
                                                     "bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected, 
                                        ]
      cfu_ordered <- NULL
      cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
      cfu <- cfu_ordered
      ttime <- as.numeric(cfu[, "days_to_death"])
      sum(status <- ttime > 0)
      deads_complete <- sum(status <- ttime > 0)
      ttime_only_top <- cfu_onlyTOP[, "days_to_death"]
      deads_top <- sum(ttime_only_top > 0)
      if (dim(cfu_onlyDOWN)[1] >= 1) {
        ttime_only_down <- cfu_onlyDOWN[, "days_to_death"]
        deads_down <- sum(ttime_only_down > 0)
      }
      else {
        deads_down <- 0
      }
      tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
      tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
      tabSurv_Matrix[i, "Cancer Deaths with Down"] <- deads_down
      tabSurv_Matrix[i, "Mean Normal"] <- mean(mRNAselected_values_normal)
      dataCancer_onlyTop_sample <- dataCancer[, samples_top_mRNA_selected]
      dataCancer_onlyTop_sample_mRNASelected <- dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == 
                                                                            mRNAselected, ]
      dataCancer_onlyDown_sample <- dataCancer[, samples_down_mRNA_selected]
      dataCancer_onlyDown_sample_mRNASelected <- dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == 
                                                                              mRNAselected, ]
      tabSurv_Matrix[i, "Mean Tumor Top"] <- mean(dataCancer_onlyTop_sample_mRNASelected)
      tabSurv_Matrix[i, "Mean Tumor Down"] <- mean(dataCancer_onlyDown_sample_mRNASelected)
      ttime[!status] <- as.numeric(cfu[!status, "days_to_last_follow_up"])
      ttime[which(ttime == -Inf)] <- 0
      ttime <- Surv(ttime, status)
      rownames(ttime) <- rownames(cfu)
      length(ttime)
      legendHigh <- paste(mRNAselected, "High")
      legendLow <- paste(mRNAselected, "Low")
      tabSurv <- survdiff(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), 
                                    rep("down", nrow(cfu_onlyDOWN))))
      tabSurv_chis <- unlist(tabSurv)$chisq
      tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), 
                                              df = 1))
      tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue
      if (Survresult == TRUE) {
        titlePlot <- paste("Kaplan-Meier Survival analysis, pvalue=", 
                           tabSurv_pvalue)
        plot(survfit(ttime ~ c(rep("low", nrow(cfu_onlyTOP)), 
                               rep("high", nrow(cfu_onlyDOWN)))), col = c("green", 
                                                                          "red"), main = titlePlot, xlab = "Days", ylab = "Survival")
        legend(100, 1, legend = c(legendLow, legendHigh), 
               col = c("green", "red"), text.col = c("green", 
                                                     "red"), pch = 15)
        print(tabSurv)
      }
    }
  }
  tabSurv_Matrix[tabSurv_Matrix == "-Inf"] <- 0
  tabSurvKM <- tabSurv_Matrix
  tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < p.cut, ]
  tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA), ]
  rownames(tabSurvKM) <- tabSurvKM$mRNA
  tabSurvKM <- tabSurvKM[, -1]
  tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing = FALSE), 
                         ]
  colnames(tabSurvKM) <- gsub("Cancer", "Group2", colnames(tabSurvKM))
  colnames(tabSurvKM) <- gsub("Tumor", "Group2", colnames(tabSurvKM))
  colnames(tabSurvKM) <- gsub("Normal", "Group1", colnames(tabSurvKM))
  return(tabSurvKM)
}
############################################################################

Genelist <- rownames(rna_vm)
results <- TCGAanalyze_SurvivalKM(clin, dataGE = rna_vm, Genelist = Genelist, ThreshTop = 0.75, ThreshDown = 0.25, Survresult = TRUE, p.cut = 0.05, group1 = group1, group2 = group2)

results1 <- results[results$`Group2 Deaths with Top` >4 & results$`Group2 Deaths with Down` >4,]
Genelist1 <-rownames(results1)
results1 <- plot_surv(clin, dataGE = rna, Genelist = Genelist1, ThreshTop = 0.6, ThreshDown = 0.4, Survresult = TRUE, p.cut = 0.01, group1 = group1, group2 = group2)
save.image("Survival_TCGAbiolinks.RData")