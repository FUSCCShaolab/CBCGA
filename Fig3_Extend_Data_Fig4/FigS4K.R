rm(list = ls())

library(stringr)
library(tidyr)
library(GSVA)
library(ggplot2)
library(ggpubr)

load("CBCGA.Extended_MergedData_V2.5_220722.Rdata")

# load SNF cluster assignments in Fig3E

# immune related signature in ISPY2

exp <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode
exp <- exp[,substring(colnames(exp),10,10)=="T"]
colnames(exp) <- substring(colnames(exp),1,4)
exp <- exp[,rownames(pheno)] %>% as.matrix()

geneset <- read.csv(file = "ISPY2_immune_signature.csv")
geneset <- split(geneset$Symbol,geneset$signature)

ssGSEA <- gsva(expr = exp , gset.idx.list = geneset, min.sz=1, max.sz=500,method="ssgsea",kcdf='Gaussian')

# Logistic regression

exp <- t(ssGSEA) %>% scale()

pheno <- CBCGA_Cohort.Info[rownames(exp),c(10,44)]
colnames(pheno) <- c("PAM50","Clinical_subtype")

SNF <- as.data.frame(SNF_Cluster)
pheno <- cbind(SNF,pheno,exp)
pheno$SNF_Cluster <- paste("SNF",pheno$SNF_Cluster,sep = "")

# reg <- pheno[pheno$Clinical_subtype=="HR+HER2-",]

reg[reg$SNF_Cluster=="SNF1","SNF_Cluster"] <- 1
reg[reg$SNF_Cluster!=1,"SNF_Cluster"] <- 0
reg$SNF_Cluster <- as.factor(reg$SNF_Cluster)

# reg[reg$PAM50=="LumA","PAM50"] <- 1
# reg[reg$PAM50!="1","PAM50"] <- 0
# reg$PAM50 <- as.factor(reg$PAM50)

gene <- colnames(reg)[4:24]

result <- matrix(nrow = length(gene),ncol = 4)
rownames(result) <- gene
colnames(result) <- c("OR","CI_low","CI_high","p")

for (i in gene) {
  
  model <- glm(SNF_Cluster~reg[,i], data = reg, family = binomial()) 
  
  result[i,1:3] <- exp(cbind('OR' = coef(model), confint(model)))["reg[, i]",]
  result[i,4] <- (summary(model)$coefficients)["reg[, i]","Pr(>|z|)"]
  
  print(which(rownames(result)==i))
}

result <- as.data.frame(result)
result$fdr <- p.adjust(result$p,method = "fdr")

Logistic_ISPY2_signature_PAM50_LumA_in_All <- result

# visualization

# p value

data <- cbind(Logistic_ISPY2_signature_PAM50_LumA[,5],Logistic_ISPY2_signature_PAM50_LumB[,5],Logistic_ISPY2_signature_PAM50_Her2[,5],Logistic_ISPY2_signature_PAM50_Basal[,5],Logistic_ISPY2_signature_PAM50_Normal[,5],
              Logistic_ISPY2_signature_refined_c1[,5],Logistic_ISPY2_signature_refined_c2[,5],Logistic_ISPY2_signature_refined_c3[,5],Logistic_ISPY2_signature_refined_c4[,5])
rownames(data) <- rownames(Logistic_ISPY2_signature_PAM50_LumA)
colnames(data) <- c("LumA","LumB","Her2","Basal","Normal","refined_C1","refined_C2","refined_C3","refined_C4")

# data <- cbind(Logistic_ISPY2_signature_PAM50_LumA_in_HRpHER2n[,5],Logistic_ISPY2_signature_PAM50_LumB_in_HRpHER2n[,5],Logistic_ISPY2_signature_PAM50_Her2_in_HRpHER2n[,5],Logistic_ISPY2_signature_PAM50_Basal_in_HRpHER2n[,5],Logistic_ISPY2_signature_PAM50_Normal_in_HRpHER2n[,5],
#               Logistic_ISPY2_signature_refined_c1_in_HRpHER2n[,5],Logistic_ISPY2_signature_refined_c2_in_HRpHER2n[,5],Logistic_ISPY2_signature_refined_c3_in_HRpHER2n[,5],Logistic_ISPY2_signature_refined_c4_in_HRpHER2n[,5])
# rownames(data) <- rownames(Logistic_ISPY2_signature_PAM50_LumA_in_HRpHER2n)
# colnames(data) <- c("LumA","LumB","Her2","Basal","Normal","refined_C1","refined_C2","refined_C3","refined_C4")

data <- -log10(data)
data <- as.data.frame(as.table(data))
colnames(data) <- c("signature","subtype","p")

# OR
data1 <- cbind(Logistic_ISPY2_signature_PAM50_LumA[,1],Logistic_ISPY2_signature_PAM50_LumB[,1],Logistic_ISPY2_signature_PAM50_Her2[,1],Logistic_ISPY2_signature_PAM50_Basal[,1],Logistic_ISPY2_signature_PAM50_Normal[,1],
               Logistic_ISPY2_signature_refined_c1[,1],Logistic_ISPY2_signature_refined_c2[,1],Logistic_ISPY2_signature_refined_c3[,1],Logistic_ISPY2_signature_refined_c4[,1])
rownames(data1) <- rownames(Logistic_ISPY2_signature_PAM50_LumA)
colnames(data1) <- c("LumA","LumB","Her2","Basal","Normal","refined_C1","refined_C2","refined_C3","refined_C4")

# data1 <- cbind(Logistic_ISPY2_signature_PAM50_LumA_in_HRpHER2n[,1],Logistic_ISPY2_signature_PAM50_LumB_in_HRpHER2n[,1],Logistic_ISPY2_signature_PAM50_Her2_in_HRpHER2n[,1],Logistic_ISPY2_signature_PAM50_Basal_in_HRpHER2n[,1],Logistic_ISPY2_signature_PAM50_Normal_in_HRpHER2n[,1],
#               Logistic_ISPY2_signature_refined_c1_in_HRpHER2n[,1],Logistic_ISPY2_signature_refined_c2_in_HRpHER2n[,1],Logistic_ISPY2_signature_refined_c3_in_HRpHER2n[,1],Logistic_ISPY2_signature_refined_c4_in_HRpHER2n[,1])
# rownames(data1) <- rownames(Logistic_ISPY2_signature_PAM50_LumA_in_HRpHER2n)
# colnames(data1) <- c("LumA","LumB","Her2","Basal","Normal","refined_C1","refined_C2","refined_C3","refined_C4")

data1 <- as.data.frame(as.table(data1))
colnames(data1) <- c("signature","subtype","OR")

plot <- cbind(data,data1[,3])
colnames(plot)[4] <- "OR"

pdf(file = "FigS4K.pdf",width = 15,height = 8)
ggplot(plot,aes(x = subtype, y = signature, size = p ,fill = OR)) + 
  geom_point(data=subset(plot,plot$p < -log10(0.05)),aes(x = subtype, y = signature),size = 10, color = "grey",fill = "grey",shape =15)+
  geom_point(color="black",shape=21)+
  geom_point(data=subset(plot,plot$p > -log10(0.05)),aes(x = subtype, y = signature),size = 10, color = "black", shape =0)+
  scale_fill_gradient2(low = "#3b53a4", mid = "#f9fafa", high =  "#ee2125",midpoint = 1)+
  coord_fixed(ratio = 13/21)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()    
