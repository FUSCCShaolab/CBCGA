rm(list=ls())

load("CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("chromosomal_location.RData")
# load CNV-RNA/CNV-PRO correlation matrix of each PAM50 subtype

subtype <- "Normal"

FDR_mRNA_filter <- 0.05
FDR_Protein_filter <- 0.05

R_mRNA_filter <- 0.4
R_Protein_filter <- 0.4 

P_mRNA_filter <- 0.05
P_Protein_filter <- 0.05

AMP_filter_alldata <- 0.5 
AMP_filter_thre <- c(2,"2")
amp_percentage_filter <- 0.1

DEL_filter_alldata <- -0.5
DEL_filter_thre <- c(-2,"-2")
del_percentage_filter <- 0.05 

# sample selection
CBCGA_clinical    <-  CBCGA_Cohort.Info

patient <-  CBCGA_clinical[CBCGA_clinical$`Intrinsic subtype (PAM50)` == subtype,1]

CBCGA_cna <- CBCGA_GISTICgene.alldata[,which(colnames(CBCGA_GISTICgene.alldata)%in%patient)] 

RNA_matrix <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,which(substring(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),10,10) == "T")]
colnames(RNA_matrix) <- substring(colnames(RNA_matrix),1,4)
CBCGA_mRNA <- RNA_matrix[,which(colnames(RNA_matrix)%in%patient)] 

pro_matrix <- CBCGA.Extended_PRO_normalized[,which(substring(colnames(CBCGA.Extended_PRO_normalized),10,10) == "T")]
colnames(pro_matrix) <- substring(colnames(pro_matrix),1,4)
CBCGA_pro <- pro_matrix[,which(colnames(pro_matrix)%in%patient)] 

sample <- Reduce(intersect,list(colnames(CBCGA_cna),colnames(CBCGA_mRNA),colnames(CBCGA_pro)))
CNV <- CBCGA_GISTICgene.thre[colnames(cna_mRNA_P.Vals),sample]

RP_matrix_all <- matrix(nrow = ncol(cna_mRNA_P.Vals),ncol =13)
rownames(RP_matrix_all) <- colnames(cna_mRNA_P.Vals)
colnames(RP_matrix_all) <- c("Chro_pos",
                             "mRNA_R","mRNA_P","mRNA_FDR",
                             "Protein_R","Protein_P","Protein_FDR",
                             "AMP_percent","AMP_check","DEL_percent","DEL_check",
                             "AMP_peak","DEL_peak")

RP_matrix_all <- as.data.frame(RP_matrix_all)

RP_matrix_all$mRNA_R <- diag(cna_mRNA_Cor.Res)
RP_matrix_all$mRNA_P <- diag(cna_mRNA_P.Vals)
RP_matrix_all$mRNA_FDR <- p.adjust(RP_matrix_all$mRNA_P,method = "fdr")

RP_matrix_all$Protein_R <- diag(cna_pro_Cor.Res)
RP_matrix_all$Protein_P <- diag(cna_pro_P.Vals)
RP_matrix_all$Protein_FDR <- p.adjust(RP_matrix_all$Protein_P,method = "fdr")

amp_fun <-function(x){
  y <-length(x)-table(as.numeric(x) %in% c("2",2))[1]
  return(y)
}

RP_matrix_all$AMP_percent <- round(apply(CNV,1,amp_fun)/ncol(CNV),3)

del_fun<-function(x){
  y<-length(x)-table(as.numeric(x) %in% c("-2",-2))[1]
  return(y)
}

RP_matrix_all$DEL_percent <- round(apply(CNV,1,del_fun)/ncol(CNV),3)

RP_matrix_all$Chro_pos <- CBCGA_GISTICgene_Cytoband[rownames(RP_matrix_all)]

RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_amp_peak),"AMP_peak"]<-Genes_in_amp_peak[rownames(RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_amp_peak),]),2] 
RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_del_peak),"DEL_peak"]<-Genes_in_del_peak[rownames(RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_del_peak),]),2] 

# gene selection
RP_matrix_all[,"AMP_check"] <- RP_matrix_all$AMP_percent>amp_percentage_filter
RP_matrix_all[,"DEL_check"] <- RP_matrix_all$DEL_percent>del_percentage_filter

# Table S3
write.csv(RP_matrix_all, file = paste("TableS3_",subtype,".csv",sep = ""))

check1 <- RP_matrix_all$mRNA_FDR < FDR_mRNA_filter
check2 <- RP_matrix_all$Protein_FDR < FDR_Protein_filter
check3 <- RP_matrix_all$mRNA_R>0 
check4 <- RP_matrix_all$Protein_R>0

check_FDR <- check1&check2&check3&check4&(RP_matrix_all[,"AMP_check"] | RP_matrix_all[,"DEL_check"])

RP_matrix_FDR <- RP_matrix_all[check_FDR,]

write.csv(RP_matrix_FDR, file = paste("Cis_",subtype,".csv",sep = ""))

# genes in amp peaks

read_index <- c("LumA","LumB","Her2","Basal","Normal")

cis <- list()

for (i in read_index) {
  temp <- read.csv(paste("Cis_",i,".csv",sep=""),row.names=1)
  temp <- temp[!is.na(temp$AMP_peak),]
  cis[[i]] <- temp
}
common_peak <- unique(c(cis[[1]][,12],cis[[2]][,12],cis[[3]][,12],cis[[4]][,12],cis[[5]][,12]))

genes_in_common_peak <- as.data.frame(Genes_in_amp_peak[Genes_in_amp_peak[,2]%in%common_peak,])
genes_in_common_peak <- genes_in_common_peak[which(genes_in_common_peak$V2!="16p13.3"),1]

temp <- read.csv(paste("TableS3_LumA.csv",sep=""),row.names=1)
genes_in_common_peak <- rownames(temp)[rownames(temp)%in%genes_in_common_peak]

FDR_plot <- list()
FDR_save <- list()

for (i in read_index) {
  temp <- read.csv(paste("TableS3_",i,".csv",sep=""),row.names=1)
  
  FDR <- temp[genes_in_common_peak,c("mRNA_FDR","Protein_FDR")]
  FDR2 <- temp[genes_in_common_peak,c("mRNA_FDR","Protein_FDR","mRNA_R","Protein_R","AMP_percent")]
  colnames(FDR)<-paste(i,colnames(FDR),sep="_")
  colnames(FDR2)<-paste(i,colnames(FDR2),sep="_")
  
  FDR_plot[[i]] <- FDR
  FDR_save[[i]] <- FDR2
}

FDR_plot <-cbind(FDR_plot[[1]],FDR_plot[[2]],FDR_plot[[3]],FDR_plot[[4]],FDR_plot[[5]])
FDR_save <-cbind(FDR_save[[1]],FDR_save[[2]],FDR_save[[3]],FDR_save[[4]],FDR_save[[5]])

FDR_check <- function(x){
  y <- sum(x<0.05)
  y <- y>=1
  return(y)
}

FDR_plot_check <- apply(FDR_plot,1,FDR_check)

FDR_plot <- FDR_plot[FDR_plot_check,]
FDR_save <- FDR_save[FDR_plot_check,]

annotation_row <- as.data.frame(Genes_in_amp_peak[rownames(FDR_save),2])
colnames(annotation_row) <- "peaks"

cut_col <- c(rep("A",2),
             rep("B",2),
             rep("C",2))

exprSet<-FDR_plot

exprSet_loged<-exprSet

exprSet_loged[exprSet_loged == 0.000000e+00] <- 10^(-4)

exprSet_loged<- -log10(exprSet_loged)

library(pheatmap)
# heatmap-FDR

pdf('Figure3D_FDR.pdf',width = 3,height = 10)

choose_matrix = exprSet_loged

choose_matrix[choose_matrix > 4] = 4

annotation_row2= data.frame(Chro_peaks =annotation_row[,1])
rownames(annotation_row2)<-rownames(annotation_row)

annotation_row2[,1]<-paste('Chr',annotation_row2[,1],sep = "")

choose_matrix=as.matrix(choose_matrix)

pheatmap::pheatmap( fontsize = 6, choose_matrix, annotation_row = annotation_row2,
                    cutree_cols = 4,
                    show_rownames = T, show_colnames = T, breaks = NA,
                    annotation_legend = T, cluster_cols = F,cluster_rows = F,
                    color = colorRampPalette(c( "white", "firebrick3"))(50),border=FALSE,
                    display_numbers = matrix(ifelse(choose_matrix > -log10(0.05), "*", ""), nrow(choose_matrix)),number_color="Black",fontsize_number=8,
                    legend_breaks = c(0:4), legend_labels = c("0","1","2","3",">=4"))

dev.off()

# heatmap-amp rate

Amp <- FDR_save[,c(5,10,15,20,25)]

choose_matrix = Amp
choose_matrix[choose_matrix > 0.20] = 0.20
choose_matrix[choose_matrix < 0.10] = 0.10

annotation_row2= data.frame(Chro_peaks =annotation_row[,1])
rownames(annotation_row2)<-rownames(annotation_row)

annotation_row2[,1]<-paste('Chr',annotation_row2[,1],sep = "")

choose_matrix=as.matrix(choose_matrix)

pdf('Figure3D_amp.pdf',width = 3,height = 10)
pheatmap::pheatmap( fontsize = 6, choose_matrix, annotation_row = annotation_row2,
                    cutree_cols = 4,
                    show_rownames = T, show_colnames = T, breaks = NA,
                    annotation_legend = T, cluster_cols = F,cluster_rows = F,
                    color = colorRampPalette(c("#fff7ed","#fdaf61"))(200),border=FALSE,
                    legend_breaks = c(0.1:0.2), legend_labels = c("<=0.1"))

dev.off()
