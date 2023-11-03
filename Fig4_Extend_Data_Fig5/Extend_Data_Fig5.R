#####################################################################################################################
## Code for Supplementary Fig. 5
## Project: CBCGA; 
## Metabolic Figure

## Contents:
## 0. Preparation and Load data 
## 1. Data filtering and cleaning
## 2. Suppl Fig. 5c 
## 3. Suppl Fig. 5d 
## 4. Suppl Fig. 5e
## 5. Suppl Fig. 5f 


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


# 2. Suppl Fig. 5c  ---------------------------------------------------------

numbers_polar <- c(101,117,240,51,55,53,26,26)
names_polar <- c("Amino acid (15.10%)","Carbohydrates (17.49%)",
                 "Lipid (35.87%)","Nucleotide (7.62%)",
                 "Other (8.22%)","Peptide (7.92%)",
                 "Vitamins and Cofactors (3.89%)","Xenobiotics (3.89%)")

numbers_lipid <- c(105,270,440,495,2)
names_lipid <- c("FA (8.00%)","SP (20.58%)","GP (33.54%)","GL (37.83%)","ST (0.15%)")


pdf(file = "Suppl Fig5c_left.pdf", height=5,width=8)
pie(numbers_polar,labels=names_polar,
    col =c("Amino acid"="#A6CEE3","Carbohydrates"="#B2DF8A",
           "Lipid"="#E31A1C","Nucleotide"="#6A3D9A",
           "Other"="#CAB2D6","Peptide"="#1F78B4",
           "Vitamins and Cofactors"="#FDBF6F","Xenobiotics"="#969696"),
    main = "Polar metabolites (n=669)",border = FALSE)
dev.off()


pdf(file = "Suppl Fig5c_right.pdf",height=5,width=8)
pie(numbers_lipid,labels=names_lipid,
    col = c("FA"="#E41A1C","SP"="#984EA3",
            "GP"="#4DAF4A","GL"="#377EB8",
            "ST"="#FF7F00"),
    main = "Lipids (n=1312)",border = FALSE)
dev.off()


# 3. Suppl Fig. 5d  ---------------------------------------------------------

group_A_matrix <- rbind(TT_pol)
group_B_matrix <- rbind(PT_pol)

#construct the comparison matrix
data_type <- "TT_vs_PT"
comparison_matrix <- matrix(ncol=6,nrow=nrow(group_A_matrix))
colnames(comparison_matrix) <- c("Mean_TT","Mean_PT","log2FC","P","FDR","name")
rownames(comparison_matrix) <- rownames(group_A_matrix)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,1] <- mean(as.numeric(group_A_matrix[i,]))
  comparison_matrix[i,2] <- mean(as.numeric(group_B_matrix[i,]))
  comparison_matrix[i,3] <- comparison_matrix[i,1]-comparison_matrix[i,2]
  a <- wilcox.test(as.numeric(group_A_matrix[i,]),as.numeric(group_B_matrix[i,]))
  comparison_matrix[i,4] <- a$p.value
}
comparison_matrix[,5] <- p.adjust(comparison_matrix[,4],method = "fdr")
comparison_matrix[,6] <- c(CBCGA_pol_anno$Putative_metabolite_name)

Peak_name <-rownames(comparison_matrix) 
comparison_matrix <- cbind(Peak_name,comparison_matrix)
#write.table(comparison_matrix,file=paste(data_type,"comparison.csv",sep="_"),sep=",",row.names=F,quote = F)
comparison_matrix<-as.data.frame(comparison_matrix,stringsAsFactors = F)
#########################       volcano plot    ########################################################
#volcano plot for polar metabolite
comparison_matrix_cat <- as.data.frame(comparison_matrix[1:nrow(CBCGA_pol_anno),1:7],stringsAsFactors = F)
comparison_matrix_cat$catagory <- as.character(CBCGA_pol_anno$Metabolite_class)
comparison_matrix_cat$catagory_sig <- comparison_matrix_cat$catagory
for (i in 2:6){
  comparison_matrix_cat[,i] <- as.numeric(comparison_matrix_cat[,i])
}
comparison_matrix_cat$log2FC <- as.numeric(comparison_matrix_cat$log2FC)

comparison_matrix_cat<-pol_comparison_matrix_cat
comparison_matrix_cat<-lip_comparison_matrix_cat


log2FC <- 5
FDR <- 0.05

for (i in rownames(comparison_matrix_cat)){
  if (abs(comparison_matrix_cat[i,"log2FC"])<log2FC | comparison_matrix_cat[i,"FDR"] >FDR) 
    comparison_matrix_cat[i,"catagory_sig"] <- "Z_insignificant"
}

table(comparison_matrix_cat$catagory_sig)

comparison_matrix_cat$label<-ifelse(comparison_matrix_cat$FDR < FDR & abs(comparison_matrix_cat$log2FC)>log2FC, comparison_matrix_cat$name,"")
comparison_matrix_cat[which(comparison_matrix_cat$label != ""),"name"]

labelmet<-c("L-Kynurenine","Valyl-Phenylalanine","L-Cystine","Leucyl-Methionine","PC(16:1(9Z)/14:0)",
            "PE(16:0/14:1(9Z))","PC(20:4(8Z,11Z,14Z,17Z)/14:0)","PC(14:1(9Z)/14:0)",
            "Aspartyl-Glutamate","Prolyl-Gamma-glutamate","Phenylalanyl-Methionine","(<U+00C2>±)-Tryptophan")
comparison_matrix_cat$label<-ifelse(comparison_matrix_cat$label %in% labelmet,comparison_matrix_cat$name,"")

library(ggplot2)

color_polar_metabolite =c("Amino acid"="#A6CEE3","Carbohydrates"="#B2DF8A",
                          "Lipid"="#E31A1C","Nucleotide"="#6A3D9A",
                          "Other"="#CAB2D6","Peptide"="#1F78B4",
                          "Vitamins and Cofactors"="#FDBF6F","Xenobiotics"="#969696",
                          "Insignificant"="#E5E5E5")

p<-ggplot(comparison_matrix_cat, aes(log2FC, -log10(FDR)))+
  geom_point(aes(col=comparison_matrix_cat$catagory_sig))+
  scale_color_manual(values=c("#A6CEE3","#B2DF8A","#E31A1C","#6A3D9A","#CAB2D6","#1F78B4","#FDBF6F","#969696","#E5E5E5"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-log2FC,log2FC), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(FDR),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(FDR)")+
  #theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  theme_all()+theme(legend.position = "none")

p<-ggplot(data = comparison_matrix_cat, aes(x=comparison_matrix_cat$log2FC,y=-log10(comparison_matrix_cat$FDR))) +
  geom_point(aes(colour=comparison_matrix_cat$catagory_sig),size=2)+
  scale_color_manual(values=c("#A6CEE3","#B2DF8A","#E31A1C","#6A3D9A","#CAB2D6","#1F78B4","#FDBF6F","#969696","#E5E5E5"))+
  theme_bw(base_size = 12, base_family = "Times") +
  xlim(-5,5) + ylim(0,10)+
  labs(x="log2FC",y="log10FDR")
p+geom_text_repel(comparison_matrix_cat, mapping=aes(log2FC,y=-log10(FDR),label = label),
                  size = 1, box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.3, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE, max.overlaps = 20)

#volcano plot for lipid (4-5 categories)
comparison_matrix_cat <- as.data.frame(comparison_matrix[(nrow(polar_metabolite_MS2_mapping)+1):nrow(comparison_matrix),1:6],stringsAsFactors = F)
comparison_matrix_cat$catagory <- lipid_MS2_mapping$Categories
comparison_matrix_cat$catagory_sig <- comparison_matrix_cat$catagory
for (i in 2:6){
  comparison_matrix_cat[,i] <- as.numeric(comparison_matrix_cat[,i])
}
comparison_matrix_cat$log2FC <- as.numeric(comparison_matrix_cat$log2FC)

#please set cutoff value for log2FC and P/FDR
log2FC <- 3
FDR <- 0.05

for (i in rownames(comparison_matrix_cat)){
  if (abs(comparison_matrix_cat[i,"log2FC"])<log2FC | comparison_matrix_cat[i,"FDR"] >FDR) 
    comparison_matrix_cat[i,"catagory_sig"] <- "Z_insignificant"
}

table(comparison_matrix_cat$catagory_sig)

library(ggplot2)

#please select optimal color and xlim, ylim

color_lipid_main <- c("FA"="#E41A1C",
                      "GL"="#377EB8",
                      "GP"="#4DAF4A",
                      "SP"="#984EA3",
                      "ST"="#FF7F00")

ggplot(comparison_matrix_cat, aes(log2FC, -log10(FDR)))+
  geom_point(aes(col=comparison_matrix_cat$catagory_sig))+
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#E5E5E5"))+
  labs(title = " ")+
  xlim(-11,11)+
  geom_vline(xintercept=c(-log2FC,log2FC), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(FDR),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(FDR)")+
  #theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  theme_all()+theme(legend.position = "none")

p+geom_text_repel(comparison_matrix_cat, mapping=aes(log2FC,y=-log10(FDR),label = label),
                  size = 1, box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.3, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE, max.overlaps = 20)


# 4. Suppl Fig. 5e  ---------------------------------------------------------

#boxplot of mean log2FC of lipid categories
boxplot(log2FC~catagory,data=comparison_matrix_cat,
        at=c(1:8),
        #xlim=c(0,6),ylim=c(-2,2),
        col=c("#A6CEE3","#B2DF8A","#E31A1C","#6A3D9A","#CAB2D6","#1F78B4","#FDBF6F","#969696"),
        #names=names(table(comparison_TT_PT_matrix_cat$catagory_whole)),
        outline=F,
        main="Fold change",
        xaxt="n",xlab=""
)
abline(h=0,lty=4,col="red",lwd=2)
axis(side=1, at=1:8, labels=FALSE)
text(names(table(comparison_matrix_cat$catagory)), x=0.5:7.5, y=-6, xpd=T, srt=30,cex=0.8)

names(table(comparison_matrix_cat$catagory))


# 5. Suppl Fig. 5f  ---------------------------------------------------------

#set test data
tested_TT_data <- rbind(TT_pol,TT_lip)
tested_PT_data <- rbind(PT_pol,PT_lip)


######################       annotate peak with KEGG and lipid abbreviation       ######################
peak_anno <- matrix(ncol=3,nrow=(nrow(CBCGA_pol_anno)+nrow(CBCGA_lip_anno)))
colnames(peak_anno) <- c("peak_name","KEGG_ID","Lipid_abbreviation")
rownames(peak_anno) <- c(rownames(CBCGA_pol_anno),rownames(CBCGA_lip_anno))
peak_anno[,1] <- rownames(peak_anno)
peak_anno[,2] <- c(CBCGA_pol_anno$KEGG,CBCGA_lip_anno$KEGG)
peak_anno[,3] <- c(rep(NA,nrow(CBCGA_pol_anno)),CBCGA_lip_anno$Subclass)
peak_anno <- as.data.frame(peak_anno,stringsAsFactors = F)

######################       construct differential abundance matrix       #############################
DA_score_matrix <- matrix(ncol=13,nrow=nrow(KEGG_annotation))
colnames(DA_score_matrix) <- c("Path_Name","Path_Class","Path_KEGG_Num","Anno_KEGG_Num","Anno_Lip_Num","Anno_Num",
                               "Up_KEGG_Num","Down_KEGG_Num","Up_Lip_Num","Down_Lip_Num","Up_Num","Down_Num",
                               "DA_Score")
rownames(DA_score_matrix) <- KEGG_annotation$Pathway
DA_score_matrix[,1:2] <- as.matrix(KEGG_annotation[,2:3])
DA_score_matrix <- as.data.frame(DA_score_matrix,stringsAsFactors=F)


########################      select differential metabolite     #######################################
## select the peak names of differential metabolites (please refer to previously estabilished tables)
#example:
comparison_matrix <- read.table("D:/研究方向/CBCGA计划/代谢组/分析结果/差异分析/subtype_specific/LumA_vs_LumB_comparison.csv",header=T,sep=",",stringsAsFactors=F)
rownames(comparison_matrix) <- comparison_matrix$Peak_name
name_1 <- rownames(comparison_matrix[which(comparison_matrix$log2FC >1 & comparison_matrix$FDR<0.01),])
name_2 <- rownames(comparison_matrix[which(comparison_matrix$log2FC < -1 & comparison_matrix$FDR<0.01),])

Up_peak_name <- name_1
Down_peak_name <- name_2

########################      calculate differential abundance score     #######################################
## if the annotated number of metabolites is less than 3, the DA score was set as NA
library(stringr)
for (i in 1:nrow(DA_score_matrix)){
  a <- strsplit(KEGG_annotation$Metabolites_1[i],";")[[1]]
  DA_score_matrix[i,"Path_KEGG_Num"] <- as.numeric(length(a))
  DA_score_matrix[i,"Anno_KEGG_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a]))
  b <- strsplit(KEGG_annotation$Metabolites_2[i],";")[[1]]
  DA_score_matrix[i,"Anno_Lip_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$Lipid_abbreviation%in%b]))
  DA_score_matrix[i,"Anno_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a | peak_anno$Lipid_abbreviation%in%b]))
  DA_score_matrix[i,"Up_KEGG_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a & rownames(peak_anno)%in%Up_peak_name]))
  DA_score_matrix[i,"Down_KEGG_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a & rownames(peak_anno)%in%Down_peak_name]))
  DA_score_matrix[i,"Up_Lip_Num"] <- as.numeric(length(rownames(peak_anno)[is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b & rownames(peak_anno)%in%Up_peak_name]))
  DA_score_matrix[i,"Down_Lip_Num"] <- as.numeric(length(rownames(peak_anno)[is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b & rownames(peak_anno)%in%Down_peak_name]))
  DA_score_matrix[i,"Up_Num"] <- as.numeric(length(rownames(peak_anno)[rownames(peak_anno)%in%Up_peak_name & (peak_anno$KEGG_ID%in%a | (is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b))]))
  DA_score_matrix[i,"Down_Num"] <- as.numeric(length(rownames(peak_anno)[rownames(peak_anno)%in%Down_peak_name & (peak_anno$KEGG_ID%in%a | (is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b))]))
  if (as.numeric(DA_score_matrix[i,"Anno_Num"])<=3) DA_score_matrix[i,"DA_Score"] <- NA 
  else DA_score_matrix[i,"DA_Score"] <- (as.numeric(DA_score_matrix[i,"Up_Num"])-as.numeric(DA_score_matrix[i,"Down_Num"]))/as.numeric(DA_score_matrix[i,"Anno_Num"])
}


for(i in 3:ncol(DA_score_matrix)){
  DA_score_matrix[,i] <- as.numeric(DA_score_matrix[,i])
}

#adjust the number and the order
DA_score_matrix <- DA_score_matrix[is.na(DA_score_matrix$DA_Score)==F,]
DA_score_matrix$Pathway_order <- c(1:nrow(DA_score_matrix))
DA_score_matrix$Path_Name_factor <- factor(DA_score_matrix$Path_Name,levels=DA_score_matrix$Path_Name)

#write.table(DA_score_matrix,file="DA_score_matrix.txt",sep="\t")


########################      DA plot    #######################################
library(ggplot2)
library(RColorBrewer)

#please select optimal color and xlim

ggplot(data=DA_score_matrix,aes(x=DA_Score,y=Path_Name_factor))+
  geom_point(aes(size=log2(DA_score_matrix$Anno_Num),color=DA_score_matrix$Path_Class)) +
  geom_segment(aes(color=DA_score_matrix$Path_Class),x=rep(0,nrow(DA_score_matrix)),xend=DA_score_matrix$DA_Score,y=DA_score_matrix$Pathway_order,yend=DA_score_matrix$Pathway_order)+
  theme_light(base_size = 12, base_family = "") +
  scale_color_manual(values=brewer.pal(11,"Paired"))+
  xlim(c(-0.5,1))






