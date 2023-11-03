rm(list = ls())

load("CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("chromosomal_location.RData")

# matrix prepare

# patient <-  CBCGA_Cohort.Info[CBCGA_Cohort.Info$`Intrinsic subtype (PAM50)`=="LumA",1]
# patient <-  CBCGA_Cohort.Info[CBCGA_Cohort.Info$`Intrinsic subtype (PAM50)`=="LumB",1]
# patient <-  CBCGA_Cohort.Info[CBCGA_Cohort.Info$`Intrinsic subtype (PAM50)`=="Her2",1]
# patient <-  CBCGA_Cohort.Info[CBCGA_Cohort.Info$`Intrinsic subtype (PAM50)`=="Basal",1]
patient <-  CBCGA_Cohort.Info[CBCGA_Cohort.Info$`Intrinsic subtype (PAM50)`=="Normal",1]

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

# CNA-RNA

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

# venn diagram
library(grid)
library(futile.logger)
library(VennDiagram)

cis_p  <- diag(cna_mRNA_P.Vals)
cis_fdr  <- p.adjust(cis_p,method = "fdr")
cis_r <- diag(cna_mRNA_Cor.Res)
RNA_result <- cbind(cis_fdr,cis_r)

cis_sig_CNV_RNA <- rownames(RNA_result[RNA_result[,"cis_fdr"]<0.05&RNA_result[,"cis_r"]>0,])

cis_p  <- diag(cna_pro_P.Vals)
cis_fdr  <- p.adjust(cis_p,method = "fdr")
cis_r <- diag(cna_pro_Cor.Res)
pro_result <- cbind(cis_fdr,cis_r)

cis_sig_CNV_pro <- rownames(pro_result[pro_result[,"cis_fdr"]<0.05&pro_result[,"cis_r"]>0,])

venn.plot <- venn.diagram(
  x = list('cis_RNA' = unique(cis_sig_CNV_RNA[is.na(cis_sig_CNV_RNA)==F]),'cis_pro' = unique(cis_sig_CNV_pro[is.na(cis_sig_CNV_pro)==F])),filename = NULL
)

pdf(file="Figure3C_Normal.pdf")
grid.draw(venn.plot)
dev.off()
