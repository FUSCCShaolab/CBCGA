rm(list = ls())
load("CBCGA.Extended_MergedData_V2.5_220722.Rdata")

# gene <- "WWP1"
gene <- "CCND1"

CBCGA_clinical    <-  CBCGA_Cohort.Info
patient <- CBCGA_clinical$PatientCode

CNV <- as.data.frame(t(CBCGA_GISTICgene.thre[gene,which(colnames(CBCGA_GISTICgene.thre)%in%patient)]))
CNV$PatientCode <- rownames(CNV)

RNA <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,which(substring(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),10,10) == "T")]
RNA <- log2(RNA+1)
colnames(RNA) <- substring(colnames(RNA),1,4)
RNA <- as.data.frame(t(RNA[gene,which(colnames(RNA)%in%patient)]))
RNA$PatientCode <- rownames(RNA)

pro <- CBCGA.Extended_PRO_normalized[,which(substring(colnames(CBCGA.Extended_PRO_normalized),10,10) == "T")]
colnames(pro) <- substring(colnames(pro),1,4)
pro <- as.data.frame(t(pro[gene,which(colnames(pro)%in%patient)])) 
pro$PatientCode <- rownames(pro)

data <- merge(CNV,RNA,all.x = T,all.y = T,by = "PatientCode")
data <- merge(data,pro,all.x = T,all.y = T,by = "PatientCode")
data <- merge(data,CBCGA_clinical[,c(1,10)],by = "PatientCode")

data <- data[(data[,2]==0|data[,2]==1|data[,2]==2)&is.na(data[,2])==F & is.na(data[,5])==F,]

data[,2] <- as.factor(data[,2])

a <- data[data$`Intrinsic subtype (PAM50)`=="Normal",]
table(is.na(a$CCND1.y))
table(is.na(a$CCND1))

library(ggplot2)

# pdf(file = "FigureS4E_upper.pdf",width = 6,height = 4)
pdf(file = "FigureS4F_upper.pdf",width = 6,height = 4)
ggplot(data, aes(x = data[,5],y = data[,3],color = data[,2]))+
  geom_boxplot(width=0.8,outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.3))+
  scale_color_manual(values = c("#34495E","#E67F22","#D60000"))+
  scale_x_discrete(limits= c("LumA","LumB","Her2","Basal","Normal"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  xlab(NULL)+ylab("mRNA level")
dev.off()

# pdf(file = "FigureS4E_lower.pdf",width = 6,height = 4)
pdf(file = "FigureS4F_lower.pdf",width = 6,height = 4)
ggplot(data, aes(x = data[,5],y = data[,4],color = data[,2]))+
  geom_boxplot(width=0.8,outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.3))+
  scale_color_manual(values = c("#34495E","#E67F22","#D60000"))+
  scale_x_discrete(limits= c("LumA","LumB","Her2","Basal","Normal"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  xlab(NULL)+ylab("Protein abundance")

dev.off()


