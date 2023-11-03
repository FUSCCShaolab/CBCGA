rm(list = ls())

load("CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("chromosomal_location.RData")

# matrix prepare

patient <-  CBCGA_Cohort.Info$PatientCode

CBCGA_cna <- CBCGA_GISTICgene.alldata[,which(colnames(CBCGA_GISTICgene.alldata)%in%patient)] 

RNA_matrix <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,which(substring(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),10,10) == "T")]
colnames(RNA_matrix) <- substring(colnames(RNA_matrix),1,4)
CBCGA_mRNA <- RNA_matrix[,which(colnames(RNA_matrix)%in%patient)] 
CBCGA_mRNA <- log2(CBCGA_mRNA+1)

pro_matrix <- CBCGA.Extended_PRO_normalized[,which(substring(colnames(CBCGA.Extended_PRO_normalized),10,10) == "T")]
colnames(pro_matrix) <- substring(colnames(pro_matrix),1,4)
CBCGA_pro <- pro_matrix[,which(colnames(pro_matrix)%in%patient)] 

common_sample <- Reduce(intersect,list(colnames(CBCGA_cna),colnames(CBCGA_mRNA),colnames(CBCGA_pro)))
common_gene <- Reduce(intersect,list(rownames(CBCGA_cna),rownames(CBCGA_mRNA),rownames(CBCGA_pro)))

CBCGA_cna <- CBCGA_cna[common_gene,common_sample]
CBCGA_mRNA <- CBCGA_mRNA[common_gene,common_sample]
CBCGA_pro <- CBCGA_pro[common_gene,common_sample]

table(colnames(CBCGA_cna)==colnames(CBCGA_mRNA))
table(colnames(CBCGA_cna)==colnames(CBCGA_pro))
table(colnames(CBCGA_mRNA)==colnames(CBCGA_pro))

gene_order <- gene_order_duplicated$gene_id
gene_order <- gene_order[gene_order%in%common_gene]

CBCGA_cna <- as.matrix(CBCGA_cna[gene_order,])
CBCGA_mRNA <- as.matrix(CBCGA_mRNA[gene_order,])
CBCGA_pro <- as.matrix(CBCGA_pro[gene_order,])

table(rownames(CBCGA_cna)==rownames(CBCGA_mRNA))
table(rownames(CBCGA_cna)==rownames(CBCGA_pro))
table(rownames(CBCGA_mRNA)==rownames(CBCGA_pro))

# CNA-RNA correlation

cna_mRNA_Cor.Res <- matrix(ncol = length(gene_order),nrow = length(gene_order))
cna_mRNA_P.Vals <- matrix(ncol = length(gene_order),nrow = length(gene_order))

rownames(cna_mRNA_Cor.Res) <- gene_order
rownames(cna_mRNA_P.Vals) <- gene_order
colnames(cna_mRNA_Cor.Res) <- gene_order
colnames(cna_mRNA_P.Vals) <- gene_order

for(i in 1:length(gene_order)) {
  
  TEMP_CNA <- CBCGA_cna[i,]
  
  for(j in 1:length(gene_order)) {
    
    TEMP_Exp <- CBCGA_mRNA[j, ]
    
    cna_mRNA_Cor.Res[i,j] <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$estimate
    cna_mRNA_P.Vals[i,j]  <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$p.value
    
    if(j/500 == round(j/500)) {print(100*round(j/length(gene_order),2))}
  }
  print(paste("i = ", i, sep = ""))
}

cna_mRNA_Cor.Res  <-  t(cna_mRNA_Cor.Res)
cna_mRNA_P.Vals   <-  t(cna_mRNA_P.Vals)

p_matrix <- cna_mRNA_P.Vals
r_matirx <- cna_mRNA_Cor.Res

cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")

cna_mRNA_P.Vals_FDR <- p_matrix

for (i in 1:ncol(cna_mRNA_P.Vals_FDR)) {
  
  cna_mRNA_P.Vals_FDR[,i] <- p.adjust(cna_mRNA_P.Vals_FDR[,i],method = "fdr")
  
  print(i)
}
diag(cna_mRNA_P.Vals_FDR)  <- cis_fdr

# chromosomal location

anno_col <- data.frame(chro=as.factor(gene_order_duplicated[gene_order,"seqnames"]))
rownames(anno_col) <- gene_order

color_24=rep(rainbow(7)[-4],5)[1:24]
names(color_24)<-unique(anno_col$chro)
ann_colors = list(chro = color_24)

# heatmap
library(pheatmap)

data <- r_matirx

for (i in 1:ncol(data)){
  
  for (j in 1:nrow(data)){
    
    if (cna_mRNA_P.Vals_FDR[j,i]>=0.05 & !is.na(cna_mRNA_P.Vals_FDR[j,i])){ data[j,i]<-NA }
  }
  print(i)
}

data[data==0] <- NA
data[data>0] <- 1
data[data<0] <- -1

data <- data[rev(1:nrow(data)),]

png('Fig3A_left.png',width = 1500, height =  1500,res = 800)

pheatmap(data, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("#3653a5", "white", "#e21a21"))(15),na_col = "White",
         annotation_col = anno_col,annotation_colors = ann_colors,annotation_legend = F,
         show_rownames = F,show_colnames = F,legend = F)

dev.off()

# CNA-PRO

cna_pro_Cor.Res <- matrix(ncol = length(gene_order),nrow = length(gene_order))
cna_pro_P.Vals <- matrix(ncol = length(gene_order),nrow = length(gene_order))

rownames(cna_pro_Cor.Res) <- gene_order
rownames(cna_pro_P.Vals) <- gene_order
colnames(cna_pro_Cor.Res) <- gene_order
colnames(cna_pro_P.Vals) <- gene_order

for(i in 1:length(gene_order)) {
  
  TEMP_CNA <- CBCGA_cna[i,]
  
  for(j in 1:length(gene_order)) {
    
    TEMP_Exp <- CBCGA_pro[j,]
    
    cna_pro_Cor.Res[i,j] <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$estimate
    cna_pro_P.Vals[i,j]  <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$p.value
    
    if(j/500 == round(j/500)) {print(100*round(j/length(gene_order),2))}
  }
  print(paste("i = ", i, sep = ""))
}

cna_pro_Cor.Res  <-  t(cna_pro_Cor.Res)
cna_pro_P.Vals   <-  t(cna_pro_P.Vals)

p_matrix <- cna_pro_P.Vals
r_matirx <- cna_pro_Cor.Res

cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")

cna_pro_P.Vals_FDR <- p_matrix

for (i in 1:ncol(cna_pro_P.Vals_FDR)) {
  
  cna_pro_P.Vals_FDR[,i] <- p.adjust(cna_pro_P.Vals_FDR[,i],method = "fdr")
  
  print(i)
}
diag(cna_pro_P.Vals_FDR)  <- cis_fdr

# chromosomal location 

anno_col <- data.frame(chro=as.factor(gene_order_duplicated[gene_order,"seqnames"]))
rownames(anno_col) <- gene_order

color_24=rep(rainbow(7)[-4],5)[1:24]
names(color_24)<-unique(anno_col$chro)
ann_colors = list(chro = color_24)

# heatmap

data <- cna_pro_Cor.Res

for (i in 1:ncol(data)){
  
  for (j in 1:nrow(data)){
    
    if (cna_pro_P.Vals_FDR[j,i]>=0.05 & !is.na(cna_pro_P.Vals_FDR[j,i])){ data[j,i]<-NA }
  }
  print(i)
}

data[data==0]   <- NA
data[data>0]  <- 1
data[data<0]  <- -1

data <- data[rev(1:nrow(data)),]

png('Figure3A_right.png',width = 1500, height =  1500,res = 800)

pheatmap(data, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("#3653a5", "white", "#e21a21"))(15),na_col = "White",
         annotation_col = anno_col,annotation_colors = ann_colors,annotation_legend = F,
         show_rownames = F,show_colnames = F,legend = F)

dev.off()

# cis histogram
library(ggplot2)

FDR_filter_for_cis_p <- 0.05

# CNA-RNA/pro

# p_matrix <- cna_mRNA_P.Vals
# r_matirx <- cna_mRNA_Cor.Res

p_matrix <- cna_pro_P.Vals
r_matirx <- cna_pro_Cor.Res

cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")
cis_fdr <- -log10(cis_fdr)

cis_r <- diag(r_matirx)

gene_list <- colnames(p_matrix)
index <- 1:ncol(p_matrix)

chr_info <- gene_order_duplicated[gene_list,"seqnames"]

plot_dataframe<-data.frame(gene_id = gene_list, chr_c = chr_info, chr_n = chr_info, chr_info_num_total = chr_info,
                           log_cis_fdr = cis_fdr,cis_r = cis_r,gene_order = index)

plot_dataframe[which(plot_dataframe$cis_r<0),"log_cis_fdr"] <- NA

gene_label <- unique(plot_dataframe$chr_n)
gene_label <-sort(as.numeric(gene_label),decreasing = F)

gene_label_order <- rep(0,length(gene_label))

for (i in gene_label){
  gene_label_order[i] <- round(median(plot_dataframe[plot_dataframe$chr_info_num_total %in% i,"gene_order"]))
}

gene_label[gene_label %in% 23] <- "X"

chr_end_order <- plot_dataframe[!duplicated(plot_dataframe$chr_info_num_total),"gene_order"]
chr_end_order <- c(chr_end_order,max(index))

# histogram

# pdf(file = "Fig3B_upper_left.pdf",width = 9,height = 1.5)
pdf(file = "Fig3B_upper_right.pdf",width = 9,height = 1.5)

ggplot(plot_dataframe, aes(x = gene_order, y = log_cis_fdr),fill=gene_order) +
  geom_bar(stat = "identity",  size = 0.1, colour = "#F17162") +
  scale_fill_manual(values = c("#F17162","#F17162"))+
  scale_x_continuous(breaks = gene_label_order,labels = gene_label)+
  scale_y_continuous(limits = c(0,30),breaks = c(0,10,20,30),labels = c(0,10,20,30))+
  geom_vline(xintercept = chr_end_order,lty=1,lwd=0.6,alpha=1,color="#9FA0A0")+
  geom_hline(yintercept =c(-log10(FDR_filter_for_cis_p)),lty=1,lwd=0.6,alpha=1,color="Red")+
  labs(x="Chromosome",y="-log10(FDR)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()

# trans histogram

FDR_filter_for_trans_R_number <- 0.05

# CNA-RNA/pro

# p_matrix <- cna_mRNA_P.Vals
# r_matirx <- cna_mRNA_Cor.Res

p_matrix <- cna_pro_P.Vals
r_matirx <- cna_pro_Cor.Res

gene_list <- colnames(p_matrix)
index <- 1:ncol(p_matrix)

cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")

FDR <- p_matrix

for (i in 1:ncol(FDR)) {
  
  FDR[,i] <- p.adjust(FDR[,i],method = "fdr")
  
  print(i)
}
diag(FDR) <- cis_fdr

trans_r_number <- rep(NA,ncol(p_matrix))

for (i in index){  
  
  trans_r_number[i] <- table(FDR[,i] < FDR_filter_for_trans_R_number)[2]  
  trans_r_number[i] <- trans_r_number[i]-as.numeric(FDR[i,i] < FDR_filter_for_trans_R_number)
  print(i)
  
}

chr_info <- gene_order_duplicated[gene_list,"seqnames"]
plot_dataframe <- data.frame(gene_id = gene_list,chr_c = chr_info, chr_n = chr_info, chr_info_num_total = chr_info,
                             trans_r_number = trans_r_number,gene_order = index)

trans_r_number<-trans_r_number/max(index)
trans_r_number<-trans_r_number*100

plot_dataframe <- data.frame(gene_id = gene_list,chr_c = chr_info, chr_n = chr_info, chr_info_num_total = chr_info,
                             trans_r_number = trans_r_number,gene_order = index)

gene_label <- unique(plot_dataframe$chr_n)
gene_label <-sort(as.numeric(gene_label),decreasing = F)

gene_label_order <- rep(0,length(gene_label))

for (i in gene_label){
  gene_label_order[i] <- round(median(plot_dataframe[plot_dataframe$chr_info_num_total %in% i,"gene_order"]))
}

gene_label[gene_label %in% 23] <- "X"

chr_end_order <- plot_dataframe[!duplicated(plot_dataframe$chr_info_num_total),"gene_order"]
chr_end_order <- c(chr_end_order,max(index))

# histogram

# pdf(file = "Fig3B_lower_left.pdf",width = 9,height = 1.5)
pdf(file = "Fig3B_lower_right.pdf",width = 9,height = 1.5)

ggplot(plot_dataframe, aes(x = gene_order, y = trans_r_number),fill=gene_order) +
  geom_bar(stat = "identity",  size = 0.1, colour = "#317EB8") +
  scale_fill_manual(values = c("#317EB8","#317EB8"))+
  scale_x_continuous(breaks = gene_label_order,labels = gene_label)+
  scale_y_continuous(limits = c(0,40),breaks = c(0,10,20,30,40),labels = c(0,10,20,30,40))+   #Fig4B
  geom_vline(xintercept = chr_end_order,lty=1,lwd=0.6,alpha=1,color="#9FA0A0")+
  labs(x="Chromosome",y="Frequency(%)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()
