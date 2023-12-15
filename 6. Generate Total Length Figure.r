

library(tidyverse)


rm(list = ls())

setwd("Your_WD")


Master_Dataset <- read.csv("Master_Dataset.csv", header = TRUE, sep = ",") 


Master_Dataset <- Master_Dataset %>%
  mutate(Risk = case_when(
    IUCN_Status %in% c("LC","NT") ~ "Non-Threatened", 
    IUCN_Status %in% c("VU", "EN", "CR", "EW","EX") ~ "Threatened", 
    #More conditions
  ))


Master_Dataset <- Master_Dataset %>% drop_na(Risk)  

Master_Dataset <- Master_Dataset %>% drop_na(tl_mm)  


Master_Dataset$tl_log=log(Master_Dataset$tl_mm, 10)


#Calculate outliers (using formula: outlier > Q3 + 1.5 x IQR or outlier < Q1 - 1.5 x IQR  )

#Calculate quartiles
quart <- quantile(Master_Dataset$tl_log)

#Calculate interquartile range
IQ<-IQR(Master_Dataset$tl_log)

L_Out<-quart[2] - 1.5 * IQ

U_Out<-quart[4] + 1.5 * IQ

#remove outliers
Master_Dataset<-subset(Master_Dataset,tl_log<U_Out)
Master_Dataset<-subset(Master_Dataset,tl_log>L_Out)

#Count histogram

TL_Count<-ggplot(Master_Dataset, aes(x=tl_log, fill=Risk)) + geom_histogram(alpha=0.5, binwidth = .075, position="identity")+
  scale_fill_manual(values=c("blue","red")) + theme_bw()+xlab("Total Length (log10)mm")+ylab("Number of Species")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_y_continuous(limits = c(0,400), expand = c(0, 0))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "bottom") +
  theme(legend.title=element_blank())+ 
  theme(legend.text=element_text(size=20))


TL_Count

save(TL_Count, file = "TL_Count.rdata")



#Percent histogram



TL_Percent<-ggplot(Master_Dataset,aes(tl_log,fill=Risk))+
  scale_fill_manual(values=c("blue","red"),name="Legend",labels=c("Non-Threatened","Threatened"))+
  geom_histogram(alpha=0.5,binwidth=0.075,position="fill")+
  xlab("Total Length (log10)mm")+ylab("Percent of Species") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_y_continuous(labels=scales::percent,limits = c(0,1), expand = c(0, 0))+
  scale_x_continuous(expand=expansion(mult=c(0,0)))+
  theme(text = element_text(size = 20)) +
  theme(legend.position = "bottom") +
  theme(legend.title=element_blank())+ 
  theme(legend.text=element_text(size=20))



TL_Percent

save(TL_Percent, file = "TL_Percent.rdata")








