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

load("CBCGA.Extended_MergedData_V2.6_220802.Rdata")
Riskscore_train <- read.csv("RNA_Met_Path_IHC_Clin_Riskscore_Train.csv")
Riskscore_test <- read.csv("RNA_Met_Path_IHC_Clin_Riskscore_Test.csv")

Cutoff = 0.5
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
Risk_train$risk <- ifelse(Risk_train$riskscore>=Cutoff,"high","low")
Risk_test$risk <- ifelse(Risk_test$riskscore>=Cutoff,"high","low")


########### Feature Heatmap ############

Risk_train$dataset <- "Train"
Risk_test$dataset <- "Test"

Risk_train <- Risk_train[,-c(2:4)]
Risk_test <- Risk_test[,-c(2:4)]

Risk_Table <- rbind(Risk_train,Risk_test)


Feat_file <- read.csv("Candidate_Feat.csv")


Feat_Metab <- Feat_file$RNA_Met_Path_IHC_Clin[1:5]
Feat_Path <- Feat_file$RNA_Met_Path_IHC_Clin[6:10]
Feat_RNA <- Feat_file$RNA_Met_Path_IHC_Clin[11:15]

CBCGA_metabolomics <- rbind(CBCGA.Extended_pol,CBCGA.Extended_lip)
CBCGA_metabolomics_T <- CBCGA_metabolomics[,str_detect(colnames(CBCGA_metabolomics),"_T")]
colnames(CBCGA_metabolomics_T) <- substring(colnames(CBCGA_metabolomics_T),1,4)
tmp <- as.data.frame(t(CBCGA_metabolomics_T))
Matrix_Metab <- tmp[Risk_Table$PatientCode,Feat_Metab]

CBCGA_RNA_T <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,str_detect(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),"_T")]
colnames(CBCGA_RNA_T) <- substring(colnames(CBCGA_RNA_T),1,4)
Matrix_RNA <- as.data.frame(t(CBCGA_RNA_T[Feat_RNA,Risk_Table$PatientCode]))

CBCGA_Path <- read.csv("CBCGA_Path_95_Feat.csv",row.names = 1)
Matrix_Path <- CBCGA_Path[Risk_Table$PatientCode,Feat_Path]

Risk_Multiomic <- cbind(Risk_Table,Matrix_Metab,Matrix_RNA,Matrix_Path)
table(is.na(Risk_Multiomic))
Risk_Multiomic_scale <- data.frame(Risk_Multiomic,row.names = 1)
Risk_Multiomic_scale <- Risk_Multiomic_scale[order(Risk_Multiomic_scale$riskscore,decreasing = F),]

for (i in c(6:20)) {
  Risk_Multiomic_scale[,i] <- as.numeric(scale(Risk_Multiomic_scale[,i]))
}


colnames(Risk_Multiomic_scale)

CBCGA_Stage_AJCC <- read.csv("Clin_Stage_AJCC_Clin_Subtype.csv")
load("SNF_C4.Rdata")
SNF_Cluster <- as.data.frame(SNF_Cluster)
SNF_Cluster$SNF_Cluster <- paste0("C",SNF_Cluster$SNF_Cluster)

anno_Clin <- Risk_Multiomic_scale[,c("riskscore","RFS_status","dataset", "risk")]
anno_Clin$IHC_subtype <- CBCGA_Cohort.Info$`Clinical subtype`[match(rownames(anno_Clin),CBCGA_Cohort.Info$PatientCode)]
anno_Clin$PAM50_subtype <- CBCGA_Cohort.Info$`Intrinsic subtype (PAM50)`[match(rownames(anno_Clin),CBCGA_Cohort.Info$PatientCode)]
anno_Clin$PAM50_subtype[is.na(anno_Clin$PAM50_subtype)] <- "Unknown"
anno_Clin$AJCC_Stage <- CBCGA_Stage_AJCC$AJCC_Stage[match(rownames(anno_Clin),CBCGA_Stage_AJCC$ID)]
anno_Clin$Refined_subtype <- SNF_Cluster$SNF_Cluster[match(rownames(anno_Clin),rownames(SNF_Cluster))]
anno_Clin <- anno_Clin[,c("IHC_subtype","AJCC_Stage","Refined_subtype","PAM50_subtype","riskscore","risk","RFS_status","dataset")]

str(anno_Clin)


cor_test <- c()

for (i in 6:20) {
  a <- cor.test(Risk_Multiomic_scale[,1],Risk_Multiomic_scale[,i], method = "spearman")
  cor_test <- rbind(cor_test,c(a$estimate,a$p.value))
}

cor_test <- as.data.frame(cor_test)
colnames(cor_test) <- c("Rho","p")
rownames(cor_test) <- colnames(Risk_Multiomic_scale)[c(6:20)]

Hm_feat_order <- rownames(cor_test)[order(cor_test$Rho,decreasing = T)]

tmp <- as.data.frame(t(Risk_Multiomic_scale[,c(6:20)]))
Hm_matrix <- tmp[Hm_feat_order,]


anno_Feat <- cor_test
anno_Feat$Feat_Cate <- c(rep("Metab",5),rep("RNA",5),rep("Path",5))
anno_Feat <- anno_Feat[,-2]


#### 5.3 Heatmap

range(Hm_matrix)

negcut <- 200
poscut <- 200

mybreaks1 <- seq(-2  , 0, by = 0.01)
mybreaks2 <- seq(0.01, 2, by = 0.01)


mybreaks <- c(mybreaks1, mybreaks2)
mycolors <- c(colorRampPalette(c("#83BCD9", "white"))(negcut), colorRampPalette(c("white", "#CB4A42"))(poscut))


anno_color <- list(
  
  AJCC_Stage = c("I"=brewer.pal(3,"Blues")[1], "II"=brewer.pal(3,"Blues")[2], "III"=brewer.pal(3,"Blues")[3]),
  IHC_subtype = c("HR+HER2-" ="#0085c4", "HR+HER2+" ="#7ab801", "HR-HER2+" ="#f2af01","TNBC"="#dc5035"),
  PAM50_subtype = c("LumA"="#0077c9", "LumB"="#74d2e8", "Her2"="#7552cd", "Basal"="#e4002c", "Normal"="#cecece","Unknown"="#f2f2f2"),
  riskscore = c("#edf8fb", "#2ca25f"),
  risk = c("low" = "#6281a5", "high" = "#b33e38"),
  dataset = c("Train"="#beced2", "Test"="#6281a5"),
  RFS_status = c("0" ="#cecece", "1" ="black"),
  Feat_Cate = c("Metab" = "#6281a5", "RNA" = "#b33e38", "Path" = "#ddb529"),
  Refined_subtype = c("C1" = "#4caf50", "C2" = "#c62828", "C3" = "#ffab00", "C4" ="#1565c0"),
  Rho = mycolors
  
)

## pheatmap

pheatmap(as.matrix(Hm_matrix), color = mycolors, breaks = mybreaks,
         annotation_col = anno_Clin, annotation_row = anno_Feat, annotation_color = anno_color, 
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE)
