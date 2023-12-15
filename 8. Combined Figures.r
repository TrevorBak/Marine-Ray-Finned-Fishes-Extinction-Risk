

library(patchwork)

rm(list = ls())


setwd("Your_WD")

#Using patchwork to combine together previously generated figures into multipanel figures

#Figure 6 Predictor variable and threat type

#N is 4915 for all 4 panels

load("Harvesting_Coef.RData")

load("Pollution_Coef.RData")

load("Climate_Change_Coef.RData")

load("Development_Coef.RData")





patch2<- (Harvesting_Coef  | Pollution_Coef ) / (Development_Coef | Climate_Change_Coef)

patch2.2<- patch2 + plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 20))



png("Best_Model_Threat_Types.png", width = 10, height = 8, units = "in", res = 600)
patch2.2 
dev.off()




#Figure 5 Threats by Threatened Non Threatened.

load("Threat_vs_non_threat.RData")

load("bar_total_threats.RData")

#Combined threats figure

patch3<-bar_total | Threat_Non_Threat




png("Threats.png", width = 12, height = 7, units = "in", res = 600)
patch3 +  plot_annotation(tag_levels = 'A')& theme(plot.tag = element_text(size = 24))
dev.off()





#Figure 2 Total Length
load("TL_Count.rdata")

load("TL_Percent.rdata")

patch4 <- TL_Count + TL_Percent + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
 


png("Total_Length.png", width = 12, height = 7, units = "in", res = 600)
patch4 + plot_annotation(tag_levels = 'A')& theme(plot.tag = element_text(size = 24))
dev.off()





#Figure 3 Critically endangered figures

load("Threat_Boxplot.rdata")

load("Eury_percent.rdata")


patch5 <- Threat_Boxplot | Eury_percent

patch5

#Save image
png("CR.png", width = 14, height = 7, units = "in", res = 600)
patch5 + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 28))
dev.off()



