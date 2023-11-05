# 3. - create threats figures, 


#Generate 5A figure, bar graph of threat numbers

library(tidyverse)

rm(list = ls())


setwd("Your_WD")


Master <- read.csv("Master_Dataset.csv")


#I want to sum across threat type columns to horizontally and vertically
#horizontal gets number of threats faced, vertical gets most common threat types

#Replacing NA's with 0's.  
Master <- Master %>% mutate_at(vars(Harvesting,
                                  Pollution, 
                                  Energy_Production_and_Mining,
                                  Natural_System_Modification,
                                  Invasive_Species_and_Disease,
                                  Climate_Change,
                                  Development,
                                  Human_Disturbance,
                                  Transportation,
                                  Aquaculture_and_Agriculture,
                                  Other,
                                  Geological_Events), ~replace_na(., 0))



Master$Threat_Row_Sums = rowSums(Master[,c("Harvesting",
                                              "Pollution", 
                                              "Energy_Production_and_Mining",
                                              "Natural_System_Modification",
                                              "Invasive_Species_and_Disease",
                                              "Climate_Change",
                                              "Development",
                                              "Human_Disturbance",
                                              "Transportation",
                                              "Aquaculture_and_Agriculture",
                                              "Other",
                                              "Geological_Events")])


Master <- Master %>%
  mutate(Risk = case_when(
    IUCN_Status %in% c("LC","NT") ~ "Not_Threatened", 
    IUCN_Status %in% c("VU", "EN", "CR") ~ "Threatened", 
    #More conditions
  ))

#Dropping out Data Deficient species.
Master <- Master %>% drop_na(Risk)


#Creating figure 5a.

Threat_Non_Threat <-ggplot(Master, aes(x=Risk, y=Threat_Row_Sums, fill=Risk)) +geom_violin(adjust=4,alpha=.5)+
  scale_y_continuous(breaks=seq(0,10,1))+
  labs(y="Number of Threats")+ scale_x_discrete(labels=c("Non-Threatened","Threatened"))+xlab("")+
  scale_fill_manual(values=c("Blue","Red"))+ theme_bw()+
  theme(legend.position = "none")+  
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 20)) +
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 3,
    # shape = 24,
    fill = "red")
  

Threat_Non_Threat

save(Threat_Non_Threat, file = "Threat_vs_non_threat.RData")


res <- wilcox.test(Threat_Row_Sums ~ Risk, paired = FALSE, exact = FALSE, conf.int = TRUE, data = Master)
res

aggregate(Master$Threat_Row_Sums~Master$Risk, FUN=median)



#create dataframe with count of columns and threat

Threat_Count <- data.frame(sum(Master$Harvesting), 
                           sum(Master$Pollution),
                           sum(Master$Energy_Production_and_Mining),
                           sum(Master$Natural_System_Modification),
                           sum(Master$Invasive_Species_and_Disease),
                           sum(Master$Climate_Change),
                           sum(Master$Development),
                           sum(Master$Human_Disturbance),
                           sum(Master$Transportation),
                           sum(Master$Aquaculture_and_Agriculture),
                           sum(Master$Other),
                           sum(Master$Geological_Events))
                           

names(Threat_Count)[1:12] <- c("Harvesting",
  "Pollution", 
  "Energy Production and Mining",
  "Natural System Modification",
  "Invasive Species and Disease",
  "Climate Change",
  "Development",
  "Human Disturbance",
  "Transportation",
  "Aquaculture and Agriculture",
  "Other",
  "Geological Events")


#Pivot longer to graph
Threat_Count <- Threat_Count %>% pivot_longer(1:12,
                                      names_to = "Threat", values_to = "Count")


#Create figure 5B.

#Plot figure
bar_total<-ggplot(data = Threat_Count, aes(x=reorder(Threat, -Count), y=Count)) + geom_col(fill="Red",alpha=.5)+
  labs(fill='Extinction Drivers') + 
  labs(y="Number of Species", x = "Threat Type")+ theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 20))



bar_total

save(bar_total, file = "bar_total_threats.RData")
  





