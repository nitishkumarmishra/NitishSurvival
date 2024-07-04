## R code from Biostar
## https://www.biostars.org/p/169808/
## Read clinical data which can be downloaded from
## https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/skcm/bcr/biotab/clin/
clinical <- read.table("PATH/TO/nationwidechildrens.org_clinical_patient_skcm.txt", as.is = TRUE,
                       sep = "\t", header = TRUE)
clinical <- clinical[3:nrow(clinical), ]
rownames(clinical) <- clinical$bcr_patient_barcode
clinical <- clinical[, c("last_contact_days_to", "death_days_to", "ajcc_pathologic_tumor_stage")]

## Create Kaplan-Meier plot (survival plot)
library(survcomp)
set.seed(12345)
survival.time <- as.integer(ifelse(clinical$last_contact_days_to == "[Not Available]",
                                   clinical$death_days_to, clinical$last_contact_days_to))
censor <- ifelse(!is.na(as.integer(clinical$death_days_to)), 1, 0)

## Some processing to get the staging right
stage <- gsub("Stage |[A-D]|/.*", "", clinical$ajcc_pathologic_tumor_stage)
stage[stage == "I"] <- 1
stage[stage == "II"] <- 2
stage[stage == "III"] <- 3
stage[stage == "IV"] <- 4
stage[stage == "[Not vailable]"] <- NA
stage <- as.factor(stage)

df <- data.frame(survival.time, censor, stage)

## And here's the business guy
km.coxph.plot(formula.s = Surv(survival.time, censor) ~ stage, data.s = df, x.label = "Time (days)",
              y.label = "Probability of survival", main.title="",
              leg.text = c("0", "1", "2", "3", "4"), leg.pos = "topright", leg.inset = 0,
              .col = c("darkblue", "darkgreen", "darkred", "black", "orange"),
              .lty = c(1, 1, 1, 1, 1), show.n.risk = TRUE, n.risk.step = 1000, n.risk.cex = 0.85,
              verbose = FALSE)