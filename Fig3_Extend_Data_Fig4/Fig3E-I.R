rm(list = ls())
load("CBCGA.Extended_MergedData_V2.5_220722.Rdata")

# SNF cluster

library(SNFtool)
library(tidyr)
library(CancerSubtypes)
library(paletteer) 
library(export)
library(ggpubr)
library(stringr)

# matrix preparation

EXP_Tumor <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,str_detect(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),fixed("_T"))]
colnames(EXP_Tumor) <- substr(colnames(EXP_Tumor),1,4)

PRO_Tumor <- CBCGA.Extended_PRO_normalized
PRO_Tumor <- PRO_Tumor[,str_detect(colnames(PRO_Tumor),fixed("_T"))]
colnames(PRO_Tumor) <- substr(colnames(PRO_Tumor),1,4)

Tes_Cases <- CBCGA_Cohort.Info$PatientCode
Tes_Cases <- Reduce(intersect,list(Tes_Cases,colnames(EXP_Tumor),colnames(PRO_Tumor)))

Tes_EXP <- EXP_Tumor[,Tes_Cases]
Tes_PRO <- PRO_Tumor[,Tes_Cases]


exp.fpkm.TT.log <- log2(Tes_EXP+1)
exp.fpkm.TT.log.SD2k <- exp.fpkm.TT.log[order(apply(exp.fpkm.TT.log,1,sd),decreasing = T),][1:2000,]

pro <- Tes_PRO[order(apply(Tes_PRO,1,sd),decreasing = T),][1:2000,]

exp.fpkm.TT.log.SD2k <- as.data.frame(t(exp.fpkm.TT.log.SD2k))
pro <- as.data.frame(t(pro))

param <- "logFPKM_pro"
K = 5
alpha = 0.8
T = 5  

exp.fpkm.m.n <- standardNormalization(exp.fpkm.TT.log.SD2k)
pro.m.n <- standardNormalization(pro)

dist_exp <- (SNFtool::dist2(as.matrix(exp.fpkm.m.n),as.matrix(exp.fpkm.m.n)))
dist_pro <- (SNFtool::dist2(as.matrix(pro.m.n),as.matrix(pro.m.n)))

W_exp <- affinityMatrix(dist_exp,K,alpha)
W_pro <- affinityMatrix(dist_pro,K,alpha)

W <- SNF(list(W_exp,W_pro),K,T)

# estimateNumberOfClustersGivenGraph

C <- 4
set.seed(1234)
group <- spectralClustering(W, C)

names(group) <- row.names(exp.fpkm.m.n)
SNF_Cluster  <- group

# heatmap
library(pheatmap)

cluster <- as.data.frame(SNF_Cluster)
cluster$PatientCode <- rownames(cluster)
cluster <- cluster[order(cluster$SNF_Cluster),]

# DEG and DEP
sample <- cluster$PatientCode
names(sample) <- cluster$SNF_Cluster

# data <- Tes_EXP
data <- Tes_PRO
data_A <- data[,sample[names(sample)==1]]
data_B <- data[,sample[names(sample)!=1]]

result <- matrix(ncol = 5,nrow = nrow(data))
rownames(result) <- rownames(data)
colnames(result) <- c("median_A","median_B","log2FC_A/B","p","fdr")

for (i in rownames(result)) {
  result[i,1] <- median(na.omit(as.numeric(data_A[i,])))
  result[i,2] <- median(na.omit(as.numeric(data_B[i,])))
  result[i,3] <- result[i,1]-result[i,2]
  result[i,4] <- wilcox.test(as.numeric(data_A[i,]),as.numeric(data_B[i,]))$p.value
  print(which(rownames(result)==i))
}
result <- as.data.frame(result)
result$fdr <- p.adjust(result$p,method = "fdr")

# DEG_4 <- rownames(result[result$`log2FC_A/B`>1 & result$fdr<0.05,])
DEP_1 <- rownames(result[result$`log2FC_A/B`>0.5 & result$fdr<0.05,])

data_RNA <- Tes_EXP[c(DEG_4,DEG_1,DEG_3,DEG_2),c(sample[names(sample)==4],sample[names(sample)==1],sample[names(sample)==3],sample[names(sample)==2])] %>% t() %>% scale() %>% t()
data_pro <- Tes_PRO[c(DEP_4,DEP_1,DEP_3,DEP_2),c(sample[names(sample)==4],sample[names(sample)==1],sample[names(sample)==3],sample[names(sample)==2])] %>% t() %>% scale() %>% t()

pos <- 2.5
neg <- -2.5
poscut <- 100
negcut <- 100
mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^ 2
mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^ 2
mybreaks <- c(mybreaks1, mybreaks2)
mycolors <- c(colorRampPalette(c("green","black"))(negcut), colorRampPalette(c("black","red"))(poscut))

anno <- CBCGA_Cohort.Info[sample,c(1,10,44)]
colnames(anno)[2:3] <- c("PAM50","Clinical_Subtype")
anno <- merge(anno,cluster,by = "PatientCode")
rownames(anno) <- anno$PatientCode
anno <- anno[,-1]

anno[anno=="NA"] <- NA

anno_color <-list(PAM50 = c("LumA"="#0077c9","LumB"="#74d2e8","Her2"="#7552cd","Basal"="#e4002c","Normal"="#cecece"),
                  Clinical_Subtype = c("HR+HER2-"="#0085c4","HR+HER2+"="#7ab801","HR-HER2+"="#f2af01","TNBC"="#dc5035"))

png(file = "Fig3E.png",width = 900,height = 600)
pheatmap(rbind(data_RNA,data_pro),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_col = anno,annotation_colors = anno_color,breaks = mybreaks,color = mycolors)
dev.off()

# K-M

library(survminer)
library(survival)

anno$PatientCode <- rownames(anno)
surv <- merge(anno,CBCGA_Cohort.Info[,c(1,50:55)],all.x = T,by = "PatientCode")
colnames(surv)[5:10] <- c("OS_status","OS_months","RFS_status","RFS_months","DMFS_status","DMFS_months")

# surv <- surv[surv$Clinical_Subtype=="HR+HER2-",]

fit_OS <- survfit(Surv(OS_months, OS_status) ~ SNF_Cluster ,surv)
fit_RFS <- survfit(Surv(RFS_months, RFS_status) ~ SNF_Cluster ,surv)
fit_DMFS <- survfit(Surv(DMFS_months, DMFS_status) ~ SNF_Cluster ,surv)

plot <-  ggsurvplot(fit_RFS,pval = TRUE, pval.coord = c(0,0.03), conf.int = F,
                    legend.title="expression", 
                    title=paste("RFS"),
                    ylab="Survival probability",xlab = " Time (months)", 
                    censor.shape = 124,censor.size = 1,
                    risk.table = TRUE, 
                    tables.height = 0.2,
                    tables.theme = theme_cleantable(),
                    risk.table.col = "strata", 
                    linetype = 1, 
                    surv.median.line = "hv", 
                    ggtheme = theme_bw())+guides(colour = guide_legend(nrow = 2))


pdf(file = "Fig3F.pdf",width = 4,height = 4)
# pdf(file = "Fig3G.pdf",width = 4,height = 4)
print(plot$plot)
dev.off()

# multivariate Cox regression

library(tableone) 
library(forestplot)
library(stringr)

data <- as.data.frame(SNF_Cluster)

# rename according to figure
data[data$SNF_Cluster==1,"Multi_omics_cluster"] <- "cluster2" 
data[data$SNF_Cluster==2,"Multi_omics_cluster"] <- "cluster4" 
data[data$SNF_Cluster==4,"Multi_omics_cluster"] <- "cluster1" 
data[data$SNF_Cluster==3,"Multi_omics_cluster"] <- "cluster3" 

data$PatientCode <- rownames(data)
data <- merge(data,CBCGA_Cohort.Info[,c(1,31,33,44,52,53)],all.x = T,by = "PatientCode")
colnames(data)[4:8] <- c("size","LN","Clinical_Subtype","RFS_status","RFS_months")

# data <- data[data$Clinical_Subtype=="HR+HER2-",]

data[which(data$size<=2),"size"] <- 0
data[which(data$size>2),"size"] <- 1
data$size <- as.factor(data$size)

data[which(data$LN==0),"LN"] <- 0
data[which(data$LN!=0),"LN"] <- 1
data$LN <- as.factor(data$LN)

data$Multi_omics_cluster <- as.factor(data$Multi_omics_cluster)
data$Multi_omics_cluster <- relevel(data$Multi_omics_cluster,ref = "cluster1")

for (i in 7:8) {
  data[,i] <- as.numeric(data[,i])
}

res.cox_RFS <- coxph(Surv(RFS_months, RFS_status) ~ Multi_omics_cluster + size + LN, data = data)

# forest plot
mul_cox1 <- summary(res.cox_RFS)
colnames(mul_cox1$conf.int)
multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
multi2<-ShowRegTable(res.cox_RFS,exp=TRUE,digits=2,pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)

result <-cbind(multi1,multi2)
result<-tibble::rownames_to_column(result, var = "Characteristics")

ins <- function(x) {c(x, rep(NA,ncol(result)-1))}

for(i in 5:6) {
  result[,i] <- as.character(result[,i])
}

result<-rbind(c("Characteristics", NA, NA, NA, "HR(95%CI)","p"),
              ins("Multi_omics_cluster"),
              ins("Cluster1"),  
              result[c(1:3),], 
              ins("Size"),
              ins("¡Ü2cm"),
              result[4,],
              ins("Lymph Node Status"),
              ins("Negative"),
              result[5,],
              c(NA, NA, NA, NA, NA,NA))

for(i in 2:4) {
  result[,i] <- as.numeric(result[,i])
}

myVars <- c("Multi_omics_cluster","Size","Lymph Node Status")
catVars <-  c("Multi_omics_cluster","Size","Lymph Node Status")

colnames(data)[3:5] <- c("Multi_omics_cluster","Size","Lymph Node Status")
table1<- print(CreateTableOne(vars = myVars,
                              data = data[,c(3:5)],
                              factorVars = catVars),
               showAllLevels=TRUE)

N<-rbind(c(NA,NA),table1[c(1:5),],c(NA,NA),table1[6:7,],c(NA,NA),table1[8:9,],c(NA,NA))       
N<-as.data.frame(N[,-1])

result1<-cbind(result,N)
result1<-result1[,c(1,7,2:6)]

for(i in 2:7) {
  result1[,i] <-  as.character(result1[,i])
}

result1<-rbind(c("Characteristics","Number(%)",NA,NA,NA,"HR (95%CI)","P.value"),result1[2:nrow(result1),])

for(i in 3:5) {
  result1[,i] <- as.numeric(result1[,i])
}

# pdf(file="Fig3H.pdf",width = 15,height = 6)   
pdf(file="Fig3I.pdf",width = 15,height = 6)  
forestplot(result1[,c(1,2,6,7)],mean=result1[,3],lower=result1[,4],upper=result1[,5],zero=1,boxsize=0.2,graph.pos= "right",
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),"2" = gpar(lty=2),"13"=  gpar(lwd=2,lty=1,columns=c(1:4))),
           graphwidth = unit(.25,"npc"), xlab="better survival  worse survival",xticks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13), 
           is.summary=c(T,T,F,F,F,F,F,T,F,F,T,F,F), 
           txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1), xlab=gpar(cex=1),title=gpar(cex=2)),lwd.zero=2,
           lwd.ci=2,lwd.xaxis=2,lty.ci=1,ci.vertices =T,ci.vertices.height=0.1,clip=c(0,13),lineheight=unit(16,'mm'),   
           line.margin=unit(16,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI",title=NULL,
           col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black"))
dev.off()

