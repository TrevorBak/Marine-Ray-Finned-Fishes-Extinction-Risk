#CR subset

library(tidyverse)
library(rstatix)
library(FSA)


rm(list = ls())

setwd("Your_WD")

Master_Dataset <- read.csv("Master_Dataset.csv", header = TRUE, sep = ",") 

Master_Dataset$tl_log=log(Master_Dataset$tl_mm, 10)

#Create Euryhaline status
Master_Dataset <- Master_Dataset %>%
  mutate(Euryhaline_Status = case_when(
    is.na(brackish) & is.na(freshwater) ~ "Marine Only",
    !is.na(brackish) & !is.na(freshwater)
    ~ "Marine Brackish and Freshwater",
    !is.na(brackish) ~ "Marine and Brackish",
    !is.na(freshwater) & !is.na(marine)~ "Marine and Freshwater"))    

#Create dataset with only variables we are interested in for ease of use
Master_Dataset_Sub <- Master_Dataset %>% select(
  scientific_name, tl_mm, tl_log, Euryhaline_Status, IUCN_Status)


#Collapse LC and NT together into a non threatened category

Master_Dataset_Sub <- Master_Dataset_Sub %>%
  mutate(IUCN_Status2 = case_when(
    IUCN_Status %in% c("LC","NT") ~ "Non-Threatened", 
  )) %>% 
  mutate(IUCN_Status = coalesce(IUCN_Status2, IUCN_Status)) %>%
  #remove temporary column
  select(!IUCN_Status2)


#Reorder levels of factor to graph in order from least threatened to most threatened

Master_Dataset_Sub$IUCN_Status <- factor(Master_Dataset_Sub$IUCN_Status,     # Reorder factor levels
                         c("Non-Threatened", "VU", "EN", "CR"))



#Total Length----

#Exploratory analysis

CR <- filter(Master_Dataset_Sub, IUCN_Status == "CR")
EN <- filter(Master_Dataset_Sub, IUCN_Status == "EN")
VU <- filter(Master_Dataset_Sub, IUCN_Status == "VU")
Not_Threat <- filter(Master_Dataset_Sub, IUCN_Status == "Non-Threatened")

hist(CR$tl_log)
hist(EN$tl_log)
hist(VU$tl_log)
hist(Not_Threat$tl_log)


#Create boxplot figure


Master_Dataset_Sub2 <- Master_Dataset_Sub %>% drop_na(IUCN_Status)  

Master_Dataset_Sub2 <- Master_Dataset_Sub2 %>% drop_na(tl_log) 


Threat_Boxplot <- ggplot(Master_Dataset_Sub2, aes(x = IUCN_Status,y = tl_log)) + 
  geom_boxplot(fill = "red", alpha = .5) +
  theme_bw()+
  ylab("Total Length (log-10)mm" ) +
  xlab("IUCN Status") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20))


Threat_Boxplot


save(Threat_Boxplot, file = "Threat_Boxplot.rdata")



#Perform Kruskal Wallis test due to distributions being non-normal

kruskal.test(tl_log ~ IUCN_Status,  data = Master_Dataset_Sub)

#Perform Dunn test to see what pairs are significant

dunnTest(tl_log ~ IUCN_Status,  data = Master_Dataset_Sub, method = "bonferroni")






#Euryhaline Status----


#Drop NAs for graph
Master_Dataset_Sub3 <- Master_Dataset_Sub %>% drop_na(Euryhaline_Status)
Master_Dataset_Sub3 <- Master_Dataset_Sub3 %>% drop_na(IUCN_Status)  

#Eury_percent <- 
Eury_percent <- ggplot(Master_Dataset_Sub3, aes(IUCN_Status, fill = Euryhaline_Status))+
  scale_fill_manual(values=c("red","purple","blue"), name="Euryhaline Status") +
  geom_bar(position = "fill", alpha = .5) +
  ylab("Percent of Species") +
  xlab("") + 
  theme_bw() +
  scale_y_continuous(labels=scales::percent,limits = c(0,1), expand = c(0, 0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.position=c(0.5, -.25), legend.direction = "vertical") +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=15)) +
  theme(plot.margin=unit(c(1,1,1.15,1), 'cm'))

Eury_percent


save(Eury_percent, file = "Eury_percent.rdata")


#Chi Square Test
chi <- chisq.test(Master_Dataset_Sub$IUCN_Status, Master_Dataset_Sub$Euryhaline_Status)


chi$expected

#Standardized residuals for supporting info
chi$stdres


#Outputting as table to help visualize
count_eury <- data.frame(table(Master_Dataset_Sub2$Euryhaline_Status, Master_Dataset_Sub2$IUCN_Status))

count_eury <- count_eury %>% pivot_wider(names_from = Var2, values_from = Freq) %>%
  rename("Euryhaline_Status" = Var1)
 
write.csv(count_eury, "Count Euryhaline Status by IUCN Status.csv")
