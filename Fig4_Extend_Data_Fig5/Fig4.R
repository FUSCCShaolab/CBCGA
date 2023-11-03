#####################################################################################################################
## Code for Fig. 4
## Project: CBCGA; 
## Metabolic Figure

## Contents:
## 0. Preparation and Load data 
## 1. Data filtering and cleaning
## 2. Fig. 4b
## 3. Fig. 4c
## 4. Fig. 4e&f
## 4. Fig. 4g


# 0. Preparation and Load data  ---------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)
set.seed(123)

setwd("/Users/yinnyxu/Desktop/CBCGA_Meta/")


library(ggpubr)
library(ggplot2)
library(dplyr)

load("CBCGA.Extended_MergedData_V2.3_220714.Rdata")
load("2021AS_MPI.Rdata")
mypro <- intersect(unique(AS_gene$b),rownames(CBCGA.Extended_PRO_normalized))

Add_cohort.info <- read.csv("Supplementary Table S4a.csv")
colnames(Add_cohort.info) <- Add_cohort.info[1,]
Add_cohort.info <- Add_cohort.info[-1,]


# 1. Data filtering and cleaning  ---------------------------------------------------------

CBCGA_ID <- rownames(CBCGA.Extended_Cohort.Info[CBCGA.Extended_Cohort.Info$CBCGA_Cohort == "Yes",]) # 773

myclinicaldata <- Add_cohort.info[,c(1,5)]
colnames(myclinicaldata) <- c("PatientCode","PAM50_classifier")
myclinicaldata <- rbind(myclinicaldata,CBCGA.Extended_Cohort.Info[CBCGA_ID,c(1,8)])
rownames(myclinicaldata) <- myclinicaldata$PatientCode

Basal_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "Basal",1] # 173
Her2_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "Her2",1] # 181
LumA_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "LumA",1] # 244
LumB_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "LumB",1] # 242
Normal_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "Normal",1] # 75

Basal_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(Basal_ID,"_PRO_T",sep=""))]) # 59
Her2_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(Her2_ID,"_PRO_T",sep=""))]) # 59
LumA_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(LumA_ID,"_PRO_T",sep=""))]) # 56
LumB_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(LumB_ID,"_PRO_T",sep=""))]) # 77
Normal_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(Normal_ID,"_PRO_T",sep=""))]) # 20
PT_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(myclinicaldata$PatientCode,"_PRO_N",sep=""))]) # 86

Basal_PRO <- Basal_PRO[mypro,] # 59
colnames(Basal_PRO) <- substr(colnames(Basal_PRO),1,4)
Her2_PRO <- Her2_PRO[mypro,] # 59
colnames(Her2_PRO) <- substr(colnames(Her2_PRO),1,4)
LumA_PRO <- LumA_PRO[mypro,] # 56
colnames(LumA_PRO) <- substr(colnames(LumA_PRO),1,4)
LumB_PRO <- LumB_PRO[mypro,] # 77
colnames(LumB_PRO) <- substr(colnames(LumB_PRO),1,4)
Normal_PRO <- Normal_PRO[mypro,] # 20
colnames(Normal_PRO) <- substr(colnames(Normal_PRO),1,4)

TT_PRO <- cbind(Basal_PRO,Her2_PRO,LumA_PRO,LumB_PRO,Normal_PRO) # 271
PT_PRO <- PT_PRO[mypro,] # 86
colnames(PT_PRO) <- substr(colnames(PT_PRO),1,4)


Basal_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(Basal_ID,"_pol_T",sep=""))] # 92
colnames(Basal_pol) <- substr(colnames(Basal_pol),1,4)
Her2_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(Her2_ID,"_pol_T",sep=""))] # 110
colnames(Her2_pol) <- substr(colnames(Her2_pol),1,4)
LumA_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(LumA_ID,"_pol_T",sep=""))] # 120
colnames(LumA_pol) <- substr(colnames(LumA_pol),1,4)
LumB_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(LumB_ID,"_pol_T",sep=""))] # 144
colnames(LumB_pol) <- substr(colnames(LumB_pol),1,4)
Normal_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(Normal_ID,"_pol_T",sep=""))] # 35
colnames(Normal_pol) <- substr(colnames(Normal_pol),1,4)

TT_pol <- cbind(Basal_pol,Her2_pol,LumA_pol,LumB_pol,Normal_pol) # 501
PT_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste(myclinicaldata$PatientCode,"_pol_N",sep=""))] # 76
colnames(PT_pol) <- substr(colnames(PT_pol),1,4)


Basal_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(Basal_ID,"_pol_T",sep=""))] # 92
colnames(Basal_lip) <- substr(colnames(Basal_lip),1,4)
Her2_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(Her2_ID,"_pol_T",sep=""))] # 110
colnames(Her2_lip) <- substr(colnames(Her2_lip),1,4)
LumA_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(LumA_ID,"_pol_T",sep=""))] # 120
colnames(LumA_lip) <- substr(colnames(LumA_lip),1,4)
LumB_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(LumB_ID,"_pol_T",sep=""))] # 144
colnames(LumB_lip) <- substr(colnames(LumB_lip),1,4)
Normal_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(Normal_ID,"_pol_T",sep=""))] # 35
colnames(Normal_lip) <- substr(colnames(Normal_lip),1,4)

TT_lip <- cbind(Basal_lip,Her2_lip,LumA_lip,LumB_lip,Normal_lip) # 501
PT_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste(myclinicaldata$PatientCode,"_pol_N",sep=""))] # 76
colnames(PT_lip) <- substr(colnames(PT_lip),1,4)


# 2. Fig. 4b  ---------------------------------------------------------

exprSet <- CBCGA.Extended_pol[,c(paste(colnames(TT_pol),"_pol_T",sep = ""),# 501
                                 paste(colnames(PT_pol),"_pol_N",sep = ""))]# 76

grouplist <- as.data.frame(matrix(ncol = 2,nrow=577))
rownames(grouplist) <- rownames(t(exprSet))
colnames(grouplist) <- c("ID","Sample_Type")
grouplist$ID <- substr(rownames(grouplist),1,4)
grouplist$Sample_Type <- ifelse(substr(rownames(grouplist),10,10)=="T","TT","PT")

grouplist <- grouplist[rownames(t(exprSet)),]

set.seed(12345)
output <- TRUE
if (output){
  dat <- exprSet
  dat[1:4,1:4]
  dat <- t(dat)
  dat <- as.data.frame(dat)
  
  tsne1 <- Rtsne(dat, dimS=6, perplexity = 10, pca=T,theta=0.5, verbose=TRUE,
                 max_iter=1000,normalize =F,num_threads = 4)
  
  lis <- grouplist$Sample_Type
  b <- data.frame(component1=tsne1$Y[,1],component2=tsne1$Y[,2],PAM50=factor(lis))
  
  summary(b$component1)
  summary(b$component2)
  
  p <- ggplot(data=b,aes(x=component2,y=component1,fill=PAM50,cluster=PAM50))+
    geom_point(shape=19,size=1.8,aes(color=PAM50))+
    scale_color_manual(values = c("#80c687","#e7a839"))+
    theme_bw()+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(0.16, "cm"),
          text = element_text(size = 12,),
          axis.text = element_text(size = 12))+
    labs(x = "tsne2",y = "tsne1")
}
p
topptx(p, filename = "Fig4b_left.pptx", width = 4.5, height = 3)


exprSet <- CBCGA.Extended_lip[,c(paste(colnames(TT_lip),"_pol_T",sep = ""),# 501
                                 paste(colnames(PT_lip),"_pol_N",sep = ""))]# 76

grouplist <- as.data.frame(matrix(ncol = 3,nrow=577))
rownames(grouplist) <- rownames(t(exprSet))
colnames(grouplist) <- c("ID","Sample_Type","PAM50")
grouplist$ID <- substr(rownames(grouplist),1,4)
grouplist$PAM50 <- myclinicaldata[grouplist$ID,2]
grouplist$Sample_Type <- ifelse(substr(rownames(grouplist),10,10)=="T","TT","PT")

grouplist <- grouplist[rownames(t(exprSet)),]

set.seed(12345)
output <- TRUE
if (output){
  dat <- exprSet
  dat[1:4,1:4]
  dat <- t(dat)
  dat <- as.data.frame(dat)
  
  tsne1 <- Rtsne(dat, dimS=6, perplexity = 10, pca=T,theta=0.5, verbose=TRUE,
                 max_iter=1000,normalize =F,num_threads = 4)
  
  lis <- grouplist$Sample_Type
  b <- data.frame(component1=tsne1$Y[,1],component2=tsne1$Y[,2],PAM50=factor(lis))
  
  summary(b$component1)
  summary(b$component2)
  
  p <- ggplot(data=b,aes(x=component2,y=component1,fill=PAM50,cluster=PAM50))+
    geom_point(shape=19,size=1.8,aes(color=PAM50))+
    scale_color_manual(values = c("#80c687","#e7a839"))+
    theme_bw()+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(0.16, "cm"),
          text = element_text(size = 12,),
          axis.text = element_text(size = 12))+
    labs(x = "tsne2",y = "tsne1")
}
p
topptx(p, filename = "Fig4b_right.pptx", width = 4.5, height = 3)


# 3. Fig. 4c  ---------------------------------------------------------

## Metabolic protein correlation network construction
pro_Cor.Res <- matrix(nrow=nrow(TT_PRO),ncol=nrow(TT_PRO))
pro_P.val <- matrix(nrow=nrow(TT_PRO),ncol=nrow(TT_PRO))
colnames(pro_Cor.Res) <- rownames(TT_PRO)
rownames(pro_Cor.Res) <- rownames(TT_PRO)
colnames(pro_P.val) <- rownames(TT_PRO)
rownames(pro_P.val) <- rownames(TT_PRO)

TT_PRO <- t(TT_PRO)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(TT_PRO)){
  for (j in 1:ncol(TT_PRO)){
    TEMP_Res <- cor.test(TT_PRO[,i],TT_PRO[,j], method = "spearman")
    pro_Cor.Res[i,j] <- TEMP_Res$estimate
    pro_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(TT_PRO))
}

pro_FDR <- matrix(p.adjust(pro_P.val,method="fdr"),ncol=2252)
colnames(pro_FDR) <- colnames(pro_P.val)
rownames(pro_FDR) <- rownames(pro_P.val)

save(pro_Cor.Res,pro_P.val,pro_FDR,file="pro_Cor_spearman.Rdata")


## Polar metabolite correlation network construction
pol_Cor.Res <- matrix(nrow=nrow(TT_pol),ncol=nrow(TT_pol))
pol_P.val <- matrix(nrow=nrow(TT_pol),ncol=nrow(TT_pol))
colnames(pol_Cor.Res) <- rownames(TT_pol)
rownames(pol_Cor.Res) <- rownames(TT_pol)
colnames(pol_P.val) <- rownames(TT_pol)
rownames(pol_P.val) <- rownames(TT_pol)

TT_pol <- t(TT_pol)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(TT_pol)){
  for (j in 1:ncol(TT_pol)){
    TEMP_Res <- cor.test(TT_pol[,i],TT_pol[,j], method = "spearman")
    pol_Cor.Res[i,j] <- TEMP_Res$estimate
    pol_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(TT_pol))
}

pol_FDR <- matrix(p.adjust(pol_P.val,method="fdr"),ncol=669)
colnames(pol_FDR) <- colnames(pol_P.val)
rownames(pol_FDR) <- rownames(pol_P.val)

save(pol_Cor.Res,pol_P.val,pol_FDR,file="pol_Cor_spearman.Rdata")


## Metabolic protein and Polar metabolite correlation network construction
CBCGA_PRO_POL_ID <- intersect(rownames(TT_PRO),rownames(TT_pol)) # 236

pro_pol_matrix <- cbind(TT_PRO[CBCGA_PRO_POL_ID,],TT_pol[CBCGA_PRO_POL_ID,]) #236*2921
dim(pro_pol_matrix)

pro_pol_Cor.Res <- matrix(nrow=ncol(pro_pol_matrix),ncol=ncol(pro_pol_matrix))
pro_pol_P.val <- matrix(nrow=ncol(pro_pol_matrix),ncol=ncol(pro_pol_matrix))
colnames(pro_pol_Cor.Res) <- colnames(pro_pol_matrix)
rownames(pro_pol_Cor.Res) <- colnames(pro_pol_matrix)
colnames(pro_pol_P.val) <- colnames(pro_pol_matrix)
rownames(pro_pol_P.val) <- colnames(pro_pol_matrix)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(pro_pol_matrix)){
  for (j in 1:ncol(pro_pol_matrix)){
    TEMP_Res <- cor.test(pro_pol_matrix[,i],pro_pol_matrix[,j], method = "spearman")
    pro_pol_Cor.Res[i,j] <- TEMP_Res$estimate
    pro_pol_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(pro_pol_matrix))
}

pro_pol_FDR <- matrix(p.adjust(pro_pol_P.val,method="fdr"),ncol=2921)
colnames(pro_pol_FDR) <- colnames(pro_pol_P.val)
rownames(pro_pol_FDR) <- rownames(pro_pol_P.val)

save(pro_pol_Cor.Res,pro_pol_P.val,pro_pol_FDR,file="pro_pol_Cor_spearman.Rdata")


## Refer to Gephi


# 4. Fig. 4e&f  ---------------------------------------------------------

CBCGA_Cohort_META <- myclinicaldata
row.names(CBCGA_Cohort_META) <- paste0(myclinicaldata$PatientCode,"_pol_T","")
lipid_data_supplementary <- CBCGA.Extended_lip
TT_ID_supple <- colnames(lipid_data_supplementary[, !str_detect(colnames(lipid_data_supplementary), fixed("_N"))])
lipid_data_supplementary_TT <- lipid_data_supplementary[,TT_ID_supple]

Basal_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="Basal",])
LumA_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="LumA",])
LumB_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="LumB",])
Her2_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="Her2",])
Normal_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="Normal",])

Basal_TT_ID <- intersect(TT_ID_supple,Basal_TT_ID)
LumA_TT_ID <- intersect(TT_ID_supple,LumA_TT_ID)
LumB_TT_ID <- intersect(TT_ID_supple,LumB_TT_ID)
Her2_TT_ID  <- intersect(TT_ID_supple,Her2_TT_ID )
Normal_TT_ID  <- intersect(TT_ID_supple,Normal_TT_ID )


lipid_TT_data <- lipid_data_supplementary_TT
TT_ID_with_meta <- colnames(lipid_TT_data)
Basal_TT_ID_meta <- intersect(Basal_TT_ID,TT_ID_with_meta)#92
LumA_TT_ID_meta <- intersect(LumA_TT_ID,TT_ID_with_meta)#120
LumB_TT_ID_meta <- intersect(LumB_TT_ID,TT_ID_with_meta)#144
Her2_TT_ID_meta <- intersect(Her2_TT_ID,TT_ID_with_meta)#110
Normal_TT_ID_meta <- intersect(Normal_TT_ID,TT_ID_with_meta)#35
nonbasal_TT_ID_meta <- union(LumA_TT_ID_meta,LumB_TT_ID_meta)
nonbasal_TT_ID_meta <- union(nonbasal_TT_ID_meta,Her2_TT_ID_meta)
nonbasal_TT_ID_meta <- union(nonbasal_TT_ID_meta,Normal_TT_ID_meta)#409

####map/PAM50
#matrix basal vs non-basal
comparison_matrix5 <- matrix(ncol=5,nrow= nrow(lipid_TT_data))
colnames(comparison_matrix5) <- c("Basal_mean","LumA_mean","LumB_mean","Her2_mean","Normal_mean")
row.names(comparison_matrix5) <- row.names(lipid_TT_data)

for (i in rownames(comparison_matrix5)){
  comparison_matrix5[i,1] <- mean(as.numeric(lipid_TT_data [i,Basal_TT_ID_meta]))
  comparison_matrix5[i,2] <- mean(as.numeric(lipid_TT_data [i,LumA_TT_ID_meta]))
  comparison_matrix5[i,3] <- mean(as.numeric(lipid_TT_data [i,LumB_TT_ID_meta]))
  comparison_matrix5[i,4] <- mean(as.numeric(lipid_TT_data [i,Her2_TT_ID_meta]))
  comparison_matrix5[i,5] <- mean(as.numeric(lipid_TT_data [i,Normal_TT_ID_meta]))
}
comparison_matrix5 <- as.data.frame(comparison_matrix5)
CBCGA_lip_anno$Subclass
id_oxpi <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPI"),])
id_oxpc <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPC"),])
id_oxpe <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPE"),])
id_oxpg <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPG"),])
id_oxps <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPS"),])
id_ox_lip <- union(id_oxpi,id_oxpc)
id_ox_lip <- union(id_ox_lip,id_oxpe)
id_ox_lip <- union(id_ox_lip,id_oxpg)
id_ox_lip <- union(id_ox_lip,id_oxps)

lipid_mapping_new_order <- CBCGA_lip_anno[id_ox_lip,]
lipid_id_new_order <- row.names(lipid_mapping_new_order)
lipid_comparison_new_order <- comparison_matrix5[lipid_id_new_order,]

Plot_Data <- lipid_comparison_new_order
Plot_Data <- lipid_comparison_new_order[, c("LumA_mean","LumB_mean","Her2_mean","Basal_mean","Normal_mean")]
cluster_data <- Plot_Data
for (i in rownames(cluster_data)){
  cluster_data[i,] <- scale(as.numeric(cluster_data[i,]), center = T, scale = T)
}


library(pheatmap)
Plot_Data       <- cluster_data
Plot_Data       <- data.frame(Plot_Data,stringsAsFactors = F) 
datexpr2  <- as.data.frame(lapply(Plot_Data,as.numeric))
row.names(datexpr2) <- row.names(Plot_Data)
Plot_Data <- datexpr2

library(RColorBrewer)
pos      <- max(Plot_Data)+0.1 ; neg <- min(Plot_Data)-0.1
poscut   <- 100             ; negcut <- 100
mybreaks1<- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^ 2
mybreaks2<- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^ 2
mybreaks <- c(mybreaks1, mybreaks2)
mycolors <- c(colorRampPalette(c("#053061", "white"))(negcut), colorRampPalette(c("white", "#67001F"))(poscut))
#plot 

library(ggplot2)
library(pheatmap)
library(eoffice)
library(gplots)
mycolors <- c(colorRampPalette(c("#0077C9", "white"))(negcut), colorRampPalette(c("white", "#E4002C"))(poscut))

P_ferro_lip_V3 <- pheatmap(Plot_Data, color = mycolors, breaks = mybreaks ,
                           fontsize_row = 8,fontsize_col = 4, cluster_rows = F,cluster_cols=F,
                           clustering_distance_cols = "euclidean", clustering_method = "ward.D", treeheight_col = 50
                           #         show_colnames = T, show_rownames = T, labels_col <- anno$metabolite_mapping_name
)
P_ferro_lip_V3
ggsave(P_ferro_lip_V3,file="Fig4f_right_ox_lipids.pdf")
dev.off()


lipid_anno_new <- CBCGA_lip_anno
lipid_TT_test<- lipid_TT_data[row.names(CBCGA_lip_anno),]

lipid_TT_test <- cbind(lipid_TT_test,lipid_anno_new$Subclass)
colnames(lipid_TT_test)[604] <- 'Abbreviation'
test_value_1 <- tapply(lipid_TT_test[,1], lipid_TT_test$Abbreviation, mean)
test_value_2 <- tapply(lipid_TT_test[,2], lipid_TT_test$Abbreviation, mean)
#######matrix
inputdata_for_matrix <- lipid_TT_test
comparison_matrix <- matrix(ncol=604,nrow= 46)
colnames(comparison_matrix) <- colnames(inputdata_for_matrix)

test_matrix <- aggregate(lipid_TT_test[,c("KPIY_pol_T","KFBP_pol_T")],by=list(lipid_TT_test$Abbreviation),mean)

comparison_matrix[,1] <- test_value_1

test_matrix_2 <- as.data.frame(test_value_1)

for (i in 1:603){
  comparison_matrix[,i] <- tapply(lipid_TT_test[,i], lipid_TT_test$Abbreviation, mean)
}
comparison_matrix <- as.data.frame(comparison_matrix)
rownames(comparison_matrix) <- rownames(test_matrix_2)
write.csv(comparison_matrix,"comparison_matrix_46class_mean_per_sample_20220806.csv")


transfer_46class_mean_per_sample <- t(comparison_matrix)
transfer_46class_mean_per_sample <- transfer_46class_mean_per_sample[-604,]


CBCGA_Cohort_META_ID <- row.names(CBCGA_Cohort_META)
CBCGA_Cohort_META_TT_ID <- intersect(TT_ID_supple,CBCGA_Cohort_META_ID)
CBCGA_Cohort_META_TT <- CBCGA_Cohort_META[CBCGA_Cohort_META_TT_ID,]
CBCGA_Cohort_META_TT_PAM50_ID <- intersect(CBCGA_Cohort_META_TT_ID,row.names(transfer_46class_mean_per_sample))
transfer_46class_mean_per_sample_PAM50 <- transfer_46class_mean_per_sample[CBCGA_Cohort_META_TT_PAM50_ID,]
CBCGA_Cohort_META_TT <- CBCGA_Cohort_META_TT[row.names(transfer_46class_mean_per_sample_PAM50),]

transfer_46class_mean_per_sample <- cbind(transfer_46class_mean_per_sample_PAM50,CBCGA_Cohort_META_TT$PAM50_classifier)
transfer_46class_mean_per_sample <- as.data.frame(transfer_46class_mean_per_sample) 
colnames(transfer_46class_mean_per_sample)[47] <- 'PAM50'
testdata <- transfer_46class_mean_per_sample
library(ggplot2)
testdata <- na.omit(testdata)
#construct matrix3

inputdata_for_matrix <- t(testdata[,1:46])

comparison_matrix <- matrix(ncol=5,nrow=nrow(inputdata_for_matrix))
rownames(comparison_matrix) <- rownames(inputdata_for_matrix)
LumA_TT_ID_meta

for (i in rownames(comparison_matrix)){
  a <- as.numeric(inputdata_for_matrix[i,LumA_TT_ID_meta])
  b <- as.numeric(inputdata_for_matrix[i,LumB_TT_ID_meta])
  c <- as.numeric(inputdata_for_matrix[i,Her2_TT_ID_meta])
  d <- as.numeric(inputdata_for_matrix[i,Basal_TT_ID_meta ])
  e <- as.numeric(inputdata_for_matrix[i,Normal_TT_ID_meta])
  comparison_matrix[i,1] <- mean(a)
  comparison_matrix[i,2] <- mean(b)
  comparison_matrix[i,3] <- mean(c)
  comparison_matrix[i,4] <- mean(d)
  comparison_matrix[i,5] <- mean(e)
}
comparison_matrix <- as.data.frame(comparison_matrix)
colnames(comparison_matrix) <- c("LumA_TT","LumB_TT","HER2E_TT","Basal_TT","Normal_TT")
comparison_matrix9 <- comparison_matrix
write.csv(comparison_matrix,"comparison_matrix9_pam50_46class.csv")

Plot_Data <- comparison_matrix9
cluster_data <- Plot_Data


for (i in rownames(cluster_data)){
  cluster_data[i,] <- scale(as.numeric(cluster_data[i,]), center = T, scale = T)
}

########plot
library(pheatmap)
Plot_Data       <- cluster_data
Plot_Data       <- data.frame(Plot_Data,stringsAsFactors = F) 
datexpr2  <- as.data.frame(lapply(Plot_Data,as.numeric))
row.names(datexpr2) <- row.names(Plot_Data)
Plot_Data <- datexpr2
library(RColorBrewer)
pos      <- max(Plot_Data)+0.1 ; neg <- min(Plot_Data)-0.1
poscut   <- 100             ; negcut <- 100
mybreaks1<- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^ 2
mybreaks2<- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^ 2
mybreaks <- c(mybreaks1, mybreaks2)
library(ggplot2)
library(pheatmap)
mycolors <- c(colorRampPalette(c("#0077C9", "white"))(negcut), colorRampPalette(c("white", "#E4002C"))(poscut))
p3 <- pheatmap(Plot_Data, color = mycolors, breaks = mybreaks ,
               fontsize_row = 8,fontsize_col = 4, cluster_rows = T,cluster_cols=F,
               clustering_distance_cols = "euclidean", clustering_method = "ward.D", treeheight_col = 50
               #         show_colnames = T, show_rownames = T, labels_col <- anno$metabolite_mapping_name
)

p3
ggsave(p3,file="Fig4e_left.pdf")
dev.off()


# 5. Fig. 4g  ---------------------------------------------------------

####comparison_matrix10 volcano_46class
testdata <- transfer_46class_mean_per_sample
testdata <- na.omit(testdata)
inputdata_for_matrix <- t(testdata[,1:46])
comparison_matrix10 <- matrix(ncol=4,nrow=nrow(inputdata_for_matrix))
rownames(comparison_matrix10) <- rownames(inputdata_for_matrix)

for (i in rownames(comparison_matrix10)){
  a <- as.numeric(inputdata_for_matrix[i,Basal_TT_ID_meta])
  b <- as.numeric(inputdata_for_matrix[i,nonbasal_TT_ID_meta])
  comparison_matrix10[i,1] <- mean(a)
  comparison_matrix10[i,2] <- mean(b)
  c<-t.test(a,b,paired = FALSE)
  comparison_matrix10[i,3]<-c$p.value
}
comparison_matrix10[,4]<-p.adjust(c(comparison_matrix10[,3]),method = "fdr")
comparison_matrix10 <- as.data.frame(comparison_matrix10)
colnames(comparison_matrix10) <- c("basal_mean","nonbasal_mean","p","fdr")

comparison_matrix10$logFC <- comparison_matrix10[,1]-comparison_matrix10[,2]
## valcano plot

Expr             <- as.data.frame(comparison_matrix10)
Expr$label       <- rownames(Expr)
selected         <- matrix(nrow = 5,ncol = 1)
selected[,1] <- c("OxPC","OxPS","OxPI","OxPE","OxPG")
colnames(selected)[1] <- "gene"
row.names(selected) <- selected[,1]
Expr$gene        <- row.names(Expr)
selectgenes      <- merge(selected, Expr, by = "gene")
Expr$logFC <- Expr$basal_mean-Expr$nonbasal_mean
row.names(selectgenes) <- selectgenes$gene
selectgenes$logFC <- selectgenes$basal_mean-selectgenes$nonbasal_mean
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
logFCcut         <- 0.5
adjPCut          <- 0.05
xmin             <- (range(Expr$logFC)[1]- (range(Expr$logFC)[1]+ 1.5))
xmax             <- (range(Expr$logFC)[1]+ (1.5-range(Expr$logFC)[1]))
ymin             <- 0
ymax             <- -log10(Expr$p)[3] * 4
Expr$color       <- ifelse((Expr$p < adjPCut & Expr$logFC > logFCcut), "red", 
                           ifelse((Expr$p < adjPCut &Expr$logFC < -logFCcut), "blue","grey"))
size             <- ifelse((Expr$p < adjPCut & abs(Expr$logFC) > logFCcut), 2, 1)

Volcano          <- ggplot(data=Expr, aes( logFC, -log10(fdr), label = label)) +
  geom_point( size = 8, colour = Expr$color) +
  labs(x=bquote(~Log[2]~"(fold change)-basalvsnon"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(breaks = c(-1, -logFCcut, 0, logFCcut, 1),        
                     labels = c(-1, -logFCcut, 0, logFCcut, 1),
                     limits = c(-2, 2)) +
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", linetype="longdash", lwd = 0.6) +
  geom_hline(yintercept = -log10(adjPCut), color="grey40", linetype="longdash", lwd = 0.6) +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  geom_point(data = selectgenes, size= 8, shape = 1, color = "black") +
  guides(color=guide_legend(title = NULL))
Volcano 

ggsave(Volcano,file="Fig4g_Volcano.pdf")
dev.off()





