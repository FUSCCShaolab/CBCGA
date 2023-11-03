rm(list=ls())

library(stringr)
library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(PMCMRplus)
library(pheatmap)
library(ggthemes)
library(export)
library(survcomp)
library(Hmisc)
library(survminer)

load("CBCGA_Cohort.Info_final.Rdata") # RFS time AND Status

Riskscore_train <- read.csv("RNA_Met_Path_IHC_Clin_Riskscore_Train.csv")
Riskscore_test <- read.csv("RNA_Met_Path_IHC_Clin_Riskscore_Test.csv")


Risk_train <- Riskscore_train[Riskscore_train$Algorithm=="avg",]
Risk_test  <- Riskscore_test[Riskscore_test$Algorithm=="avg",]

colnames(Risk_train)[5] <- "riskscore"
colnames(Risk_test) [5] <- "riskscore"

Risk_train$RFS_time <- CBCGA_Cohort.Info$`RFS time (month)`[match(Risk_train$PatientCode,CBCGA_Cohort.Info$PatientCode)]
Risk_train$RFS_status <- CBCGA_Cohort.Info$`RFS status`[match(Risk_train$PatientCode,CBCGA_Cohort.Info$PatientCode)]

Risk_test$RFS_time <- CBCGA_Cohort.Info$`RFS time (month)`[match(Risk_test$PatientCode,CBCGA_Cohort.Info$PatientCode)]
Risk_test$RFS_status <- CBCGA_Cohort.Info$`RFS status`[match(Risk_test$PatientCode,CBCGA_Cohort.Info$PatientCode)]

Risk_train$riskscore <- as.numeric(scale(Risk_train$riskscore))
Risk_test$riskscore <- as.numeric(scale(Risk_test$riskscore))

Cutoff = 0.5

Risk_train$risk <- ifelse(Risk_train$riskscore>=Cutoff,"high","low")
mypal <- c("#6281a5","#b33e38")


####### KM-plot in Multi Trainset

# Model Fit
surv_model_train <- survfit(Surv(RFS_time, RFS_status) ~ risk, data = Risk_train)
surv_pvalue(surv_model_train) # P-val = 1.704792e-23

# KMplot
survPlot_TMPIC_Train <- ggsurvplot(surv_model_train, 
                                   data = Risk_train,
                                   conf.int = FALSE,  
                                   palette = mypal, 
                                   pval = TRUE,
                                   pval.coord = c(0, 0.6),
                                   pval.size = 5,
                                   pval.method = TRUE,
                                   legend = c(0.8,0.2),
                                   legend.title = "Risk Score",
                                   linetype = "solid",
                                   size = 2,
                                   ylim = c(0,1),
                                   xlab = "Follow up (Months)",
                                   ylab = "RFS Probability",
                                   break.time.by = 12,
                                   risk.table = FALSE,
                                   censor = TRUE,
                                   censor.shape = 124,
                                   censor.size = 5
)

survPlot_TMPIC_Train


###### Risk in Multi Testset


Risk_test$risk <- ifelse(Risk_test$riskscore>=Cutoff,"high","low")


####### KM-plot in Multi testset

# Model Fit
surv_model_test <- survfit(Surv(RFS_time, RFS_status) ~ risk, data = Risk_test)
surv_pvalue(surv_model_test) # P-val = 1.189564e-05

# KMplot
survPlot_TMPIC_Test <- ggsurvplot(surv_model_test, 
                                  data = Risk_test,
                                  conf.int = FALSE,  
                                  palette = mypal, 
                                  pval = TRUE,
                                  pval.coord = c(0, 0.6),
                                  pval.size = 5,
                                  pval.method = TRUE,
                                  legend = c(0.8,0.2),
                                  legend.title = "Risk Score",
                                  linetype = "solid",
                                  size = 2,
                                  ylim = c(0,1),
                                  xlab = "Follow up (Months)",
                                  ylab = "RFS Probability",
                                  break.time.by = 12,
                                  risk.table = FALSE,
                                  censor = TRUE,
                                  censor.shape = 124,
                                  censor.size = 5
)

survPlot_TMPIC_Test

