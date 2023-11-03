rm(list = ls())

library(tidyr)
# load TME ssGSEA score in Fig5
# load SNF assignments in Fig3E

TME <- CBCGA_Microenvironment_ssGSEA[names(SNF_Cluster),] %>% scale() %>% t()
TME <- TME[,c(names(SNF_Cluster[SNF_Cluster==3]),names(SNF_Cluster[SNF_Cluster!=3]))]

pos <- 2.5
neg <- -2.5
poscut <- 100
negcut <- 100
mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^ 2
mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^ 2
mybreaks <- c(mybreaks1, mybreaks2)
mycolors <- c(colorRampPalette(c("green","black"))(negcut), colorRampPalette(c("black","red"))(poscut))

anno <- as.data.frame(SNF_Cluster)

png(file = "Fig3E.png",width = 900,height = 600)
pheatmap(TME,cluster_rows = T,cluster_cols = F,show_rownames = T,show_colnames = F,annotation_col = anno,breaks = mybreaks,color = mycolors)
dev.off()
