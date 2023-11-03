rm(list = ls())

# calculated the CNA-RNA/CNA-PRO correlation as in Fig3A

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

# pdf(file="FigureS4C_CPTAC_2020.pdf")
pdf(file="FigureS4D_CPTAC_2020.pdf")
grid.draw(venn.plot)
dev.off()
