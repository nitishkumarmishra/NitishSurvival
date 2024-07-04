library(TCGAbiolinks)
library(dplyr)
library(survival)
clinical_paad <- TCGAquery_clinic("paad","clinical_patient")
clinical_paad$new_death <- c()
for (i in 1:length(as.numeric(as.character(clinical_paad$days_to_death)))){
  clinical_paad$new_death[i] <- ifelse(is.na(as.numeric(as.character(clinical_paad$days_to_death))[i]),
                                       as.numeric(as.character(clinical_paad$days_to_last_followup))[i],as.numeric(as.character(clinical_paad$days_to_death))[i])
}
clinical_paad$death_event <- ifelse(clinical_paad$vital_status == "Alive", 0,1)
rownames(clinical_paad) <- clinical_paad$bcr_patient_barcode
## dplr command filter will select only those rows which have days_to_death >0
## I can do this without filter command also
#paad_survival <- filter(clinical_paad, days_to_death >0) ### Very good command, but at this point I am not using

## 1 Month = 30.4167 days
# paad_survival$bcr_patient_barcode # paad_survival$vital_status # paad_survival$histological_type
# paad_survival$gender # paad_survival$days_to_birth # paad_survival$prior_dx
# paad_survival$year_of_initial_pathologic_diagnosis # paad_survival$lymph_node_examined_count
# paad_survival$neoplasm_histologic_grade # paad_survival$maximum_tumor_dimension
# paad_survival$pathologic_stage # paad_survival$days_to_death # paad_survival$patient_death_reason
# paad_survival$age_at_initial_pathologic_diagnosis  # paad_survival$days_to_last_followup
# attach(paad_survival)
# coxph(formula = Surv(as.numeric(days_to_death)) ~ gender + ethnicity)
# # To get P-value for each variables
# tmp <- summary(coxph(formula = Surv(as.numeric(days_to_death)) ~ gender+pathologic_stage))
# p.value <- tmp$coefficients[,5]
## Use survival data in month by dividing 30.4167
coxph(formula = Surv(as.numeric(clinical_paad$new_death/30.4167), clinical_paad$death_event) ~ clinical_paad$gender)
plot(survfit(Surv(as.numeric(clinical_paad$new_death/30.4167), clinical_paad$death_event) ~ clinical_paad$neoplasm_histologic_grade))
#sample(clinical_paad$new_death) ## Randomly shuffle the new_death column

###### How to get P.val from survdiff log rank test ######
#https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
#http://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
#data.survdiff <- survdiff(Surv(time, status) ~ group)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
################################
tt <- survdiff(Surv(tmp,tmp1)~tmp3)
p.value <- 1 - pchisq(tt$chisq, length(tt$n) - 1)


########### This part of code is still not tested #########
## X ==Survival age, A== Censor, Y== MPS data for each gene, B=MPS-class level
myanalysis = function(X,Y){
  ntests = ncol(Y)
  rslts = as.data.frame(matrix(NA,nrow=ntests,ncol=2))
  names(rslts) = c("ID","pvalue")
  rslts[,"ID"] = 1:ntests
  for(i in 1:ntests){
    #fit = survdiff(Surv(as.numeric(clinical_paad$new_death/30.4167),clinical_paad$death_event)~tmp3)
    fit = survdiff(Surv(as.numeric(X/30.4167),censor)~ MPS+strata(MPS.class))
    p.value <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
    rslts[i,"pvalue"] = fit$p.value
  }
  return(rslts)
  return(rslts)
} # End myanalysis
## Generate observed results
obs = myanalysis(X,Y)
## Generate permuted results
perml = vector('list',nperm)
for(p_ in 1:nperm){
  X1 = X[order(runif(ncol_)),]
  ##I think in my case both "sample" and "order" command work better
  perml[[p_]] = myanalysis(X1,Y)
}
## FDR results table
fdrTbl(obs$pvalue,perml,"pvalue",ncol_,-2, -1.30103)
#log10(0.01)=-2,log10(0.05)=-1.30103-- We can change lower and upper cutoff
