rm(list = ls())

# run GSEA software according to SNF cluster assignment in Fig3E

# visualization

library(ggplot2)

GSEA <- read.csv("GSEA_SNF_cluster.csv")
GSEA$Pathway <- paste(GSEA$Pathway,GSEA$annotation)

nameorder <- GSEA$Pathway[order(GSEA$annotation, GSEA$NES)] 
GSEA$Pathway <- factor(GSEA$Pathway, levels=nameorder)

pdf(file = "FigS4I.pdf",width = 8,height = 6)
ggplot(GSEA, aes(x=NES, y=Pathway))+ 
  geom_segment(aes(yend=Pathway), xend=0, colour="grey50")+
  geom_point(size=3, aes(colour=annotation))+
  scale_colour_brewer(palette="Set1", limits=c("cluster1","cluster2","cluster3","cluster4"),guide = F)+
  # scale_x_continuous(breaks = c(0,0.5,1,1,5,2,2,5))+
  theme_bw() + 
  theme(panel.grid.major.y = element_blank()) +
  facet_grid(annotation~., scales="free_y", space="free_y")
dev.off()
