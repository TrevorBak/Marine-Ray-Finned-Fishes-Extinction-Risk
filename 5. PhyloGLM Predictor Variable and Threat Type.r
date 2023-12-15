#PhyloGLM Best Model Association with Threat Types


#Packages----
library(ape) #for reading phylogenetic trees
library(caper) #for loops that match species in phylo tree to species in dataset
library(phylolm) #runs PhyloGLM
library(qpcR) # for akaike.weights function
library(rr2) # for R2 function
library(vioplot) #for violin plots
library(tidyverse)

rm(list = ls())

#Load in Data----
setwd("Your_WD")

#read in phylogenetic data and dataset
temp=list.files(pattern="*.tre$")

mytree=lapply(temp,read.tree)

Master_Dataset <- read.csv("Master_Dataset.csv", 
                           header = TRUE, sep = ",") 


Master_Dataset$tl_mm <- as.numeric(Master_Dataset$tl_mm)

#Modify Dataset----
#building binary risk categorization
Master_Dataset <- Master_Dataset %>%
  mutate(Risk = case_when(
    IUCN_Status %in% c("LC","NT") ~ "0", 
    IUCN_Status %in% c("VU", "EN", "CR") ~ "1", 
    #More conditions
  ))


#Log transforming length
Master_Dataset$tl_mmlog=log(Master_Dataset$tl_mm, 10)
#Create tiering category

#Create Euryhaline status
Master_Dataset <- Master_Dataset %>%
  mutate(Euryhaline_Status = case_when(
    is.na(brackish) & is.na(freshwater) ~ "Marine Only",
    !is.na(brackish) & !is.na(freshwater)
    ~ "Marine Brackish and Freshwater",
    !is.na(brackish) ~ "Marine and Brackish",
    !is.na(freshwater) & !is.na(marine)~ "Marine and Freshwater"))    


# substitutes the space between the genus and species with an underscore to match the species names in the phylogenetic tree
Master_Dataset$scientific_name<-sub(" ", "_", Master_Dataset$scientific_name) 

#I want to replace NA's with 0's so they are pulled into the phyloGLM.  0 just means
#threat type has not been assigned.

# select columns defined in vars(col1, col2, ...):
Master_Dataset <- Master_Dataset %>%
  mutate_at(vars(Harvesting, Pollution, Climate_Change, Development), ~replace_na(.,0))


#builds subset that will be used for analysis with only relevant variables and 
#ensures variables are correct format 
eco.species_wFactors<-cbind.data.frame(Master_Dataset[,c("scientific_name","IUCN_Status")], 
                                       TL_Log=as.numeric(Master_Dataset[,"tl_mmlog"]),
                                       Euryhaline_Status=as.factor(Master_Dataset[,"Euryhaline_Status"]),
                                       Harvesting=as.factor(Master_Dataset[,"Harvesting"]),
                                       Pollution=as.factor(Master_Dataset[,"Pollution"]),
                                       Climate_Change=as.factor(Master_Dataset[,"Climate_Change"]),
                                       Development=as.factor(Master_Dataset[,"Development"]),
                                       risk=as.factor(Master_Dataset[,"Risk"]))

#Match Phylo and Dataset Species----

# make sure the species names between the eco data and tree match 
# matches the names in the eco.species dataset to the names in the tree tip.label. If there is no match, an NA is returned
sp.matching<-vector(mode = "list", length (mytree))
for(j in 1:length(mytree)){
  sp.matching[[j]]<-match(eco.species_wFactors$scientific_name, mytree[[j]]$tip.label)
}

# keeps only the species names for which there was a match found. Others will be dropped because they had an NA
keep.names<-vector(mode = "list", length = 100)
for(j in 1:length(mytree)){
  keep.names[[j]]<-sp.matching[[j]][!is.na(sp.matching[[j]])]
}

# a vector of spp names we want to keep for the analyses
mam.meet<-vector(mode = "list", length = 100)
for(j in 1:length(mytree)){
  mam.meet[[j]]<-mytree[[j]]$tip.label[keep.names[[j]]]
}

# creates a list of spp we want to remove from the tree. It looks at our list of keepers (keep.names).
#If a species is not in our "keep" list it is now in a vector of names to drop
droppers<-vector(mode = "list", length = 100)
for(j in 1:length(mytree)){
  droppers[[j]]<-mytree[[j]]$tip.label[!mytree[[j]]$tip.label %in% mam.meet[[j]]]
}


# drops all of the spp from the "droppers" list (made above) from the tree 
Trimmedtrees<-vector(mode = "list", length = 100)
for(j in 1:length(mytree)){
  Trimmedtrees[[j]]<-drop.tip(mytree[[j]], droppers[[j]])
}

# matches the ecology data with the species from the tree to make a comparative data set
treedat<-vector(mode = "list", length = 100)
for(j in 1:length(mytree)){
  treedat[[j]]<-comparative.data(Trimmedtrees[[j]], eco.species_wFactors, "scientific_name", na.omit=FALSE)
}

phydat<-vector(mode = "list", length = 100)
for(j in 1:length(mytree)){
  #phydat[[j]]<-data.frame(treedat[[j]]$data)
  phydat[[j]]<- treedat[[j]]$data # don't want to make this a data.frame again -> will destroy the "factor" coding
}


for(j in 1:length(mytree)){
  phydat[[j]]$Euryhaline_Status<-relevel(phydat[[j]]$Euryhaline_Status, ref="Marine Only")}



#Start PhyloGLM Here----

model1 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model2 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model3 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model4 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results

##name the models to be averaged

Modnames <- c("Harvesting", "Pollution", "Climate_Change", "Development")
# 
# empty list to catch the 100 AIC results
aicResults <- vector(mode = "list", length = 100)
modelResults <- vector(mode = "list", length = 100)
bestModel_coef <- vector(mode = "list", length = 100)
aicWeights <- vector(mode = "list", length = 100)

for(j in 1: length(treedat)){
  one_tree <- treedat[[j]]$phy
  df <- phydat[[j]] 
  model1[[j]] <- phyloglm(Harvesting~TL_Log + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  model2[[j]] <- phyloglm(Pollution~TL_Log + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  model3[[j]] <- phyloglm(Climate_Change~TL_Log + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  model4[[j]] <- phyloglm(Development~TL_Log + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  
  allModels_j <- list(model1[[j]], model2[[j]], model3[[j]], model4[[j]])
  
  allSummaries_j <- lapply(allModels_j, summary)
  allAIC_j <- lapply(allModels_j, AIC)
  allR2_j <- as.matrix(unlist(c(NA, R2(model2[[j]],model1[[j]]), R2(model3[[j]],model1[[j]]), R2(model4[[j]],model1[[j]]))))
  
  bestModel_j <- which( allAIC_j==min(unlist(allAIC_j)) ) # get the best model of 4 (lowest AIC)
  deltaAIC_j <- unlist(allAIC_j)-allAIC_j[[bestModel_j]] # get the delta AIC values
  rank_j <- order(deltaAIC_j) # get the ranks of the 4 models
  
  # aic weights for all 8 models 
  aicWeights[[j]] <- akaike.weights( unlist(allAIC_j) )$weights
  
  ##model selection table based on AICc (now includes aicWeights too)
  aicResults[[j]] <- cbind.data.frame( aicWeight = round(aicWeights[[j]],6), deltaAIC = deltaAIC_j, rank =  rank_j, model = Modnames, tree = rep(j, length(Modnames)), R2 = allR2_j)
  modelResults[[j]] <- allModels_j # save the model results
  
  #summary results for the best model (model8[[j]])
  # getting coefficients for MODEL 5
  bestModel_coef[[j]] <- model1[[j]]$coefficients

  
}



#Saving each model's coefficients to it's own object for graphing


bestModel_coef1 <- vector(mode = "list", length = 100)
for(j in 1: length(treedat)){
  
  bestModel_coef1[[j]] <- model1[[j]]$coefficients
}

bestModel_coef2 <- vector(mode = "list", length = 100)
for(j in 1: length(treedat)){
  
bestModel_coef2[[j]] <- model2[[j]]$coefficients
}

bestModel_coef3 <- vector(mode = "list", length = 100)
for(j in 1: length(treedat)){
  
  bestModel_coef3[[j]] <- model3[[j]]$coefficients
}

bestModel_coef4 <- vector(mode = "list", length = 100)
for(j in 1: length(treedat)){
  
  bestModel_coef4[[j]] <- model4[[j]]$coefficients
}



#Writing to csv and bringing back in helps get the data in the form needed for graphing
write.csv(bestModel_coef1,file = "summaryResults.csv")

write.csv(bestModel_coef2,file = "summaryResults2.csv")

write.csv(bestModel_coef3,file = "summaryResults3.csv")

write.csv(bestModel_coef4,file = "summaryResults4.csv")



coef_df_Harvesting <- read.csv("summaryResults.csv")

coef_df_Harvesting2 <- coef_df_Harvesting %>%pivot_longer(!X, names_to = "col_head", values_to = "coef")

coef_df_Harvesting2<-subset(coef_df_Harvesting2, subset = !X %in% "(Intercept)")


Harvesting_Coef<-ggplot(coef_df_Harvesting2, aes(x=X, y=coef)) + geom_violin(fill="#31AE3B")+
  scale_y_continuous(breaks=seq(0,2.5,.5))+
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 3) +
  xlab("")+ylab("Extinction Risk (log-odds)")+
  scale_x_discrete(labels=c("","",""))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  theme_classic() + geom_hline(yintercept = 0,linetype="dotted")+
  theme(text = element_text(size = 14))+
  ggtitle("Harvesting")+ 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

Harvesting_Coef
#save as object

save(Harvesting_Coef, file = "Harvesting_Coef.RData")

#Calculating CI's as another way of getting at variance, do not report
#In the manuscript.

results_Harvesting <- group_by(coef_df_Harvesting2, X) %>% summarise(Coef = median(coef), Var = var(coef)) %>% 
  rename(Variable_Level=X)

results_Harvesting$SE <- sqrt(results_Harvesting$Var)

results_Harvesting$L95 <- results_Harvesting$Coef - (1.96 * results_Harvesting$SE)
results_Harvesting$U95 <- results_Harvesting$Coef + (1.96 * results_Harvesting$SE)

write.csv(results_Harvesting, "Harvesting_Model_Results.csv")





#Pollution
coef_df_Pollution <- read.csv("summaryResults2.csv")
#coef_df$trait <- factor(coef_df$X , levels=c("Fossorial","Cave-dwelling","Arboreal","Ground-dwelling","Pedal","Jumping","Brachiating"))


coef_df_Pollution2 <- coef_df_Pollution %>%pivot_longer(!X, names_to = "col_head", values_to = "coef")

coef_df_Pollution2<-subset(coef_df_Pollution2, subset = !X %in% "(Intercept)")


Pollution_Coef<-ggplot(coef_df_Pollution2, aes(x=X, y=coef)) + geom_violin(fill="#31AE3B")+
  scale_y_continuous(breaks=seq(-.5,2.5,.5),limits = c(-.5, 2.5))  +
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 3) +
  xlab("")+ylab("")+
  scale_x_discrete(labels=c("","",""))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  theme_classic() + geom_hline(yintercept = 0,linetype="dotted")+
  theme(text = element_text(size = 14))+
  ggtitle("Pollution")+ 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

Pollution_Coef
#save as object

save(Pollution_Coef, file = "Pollution_Coef.RData")


results_Pollution <- group_by(coef_df_Pollution2, X) %>% summarise(Coef = median(coef), Var = var(coef)) %>% 
  rename(Variable_Level=X)

results_Pollution$SE <- sqrt(results_Pollution$Var)

results_Pollution$L95 <- results_Pollution$Coef - (1.96 * results_Pollution$SE)
results_Pollution$U95 <- results_Pollution$Coef + (1.96 * results_Pollution$SE)

write.csv(results_Pollution, "Pollution_Model_Results.csv")

#Climate Change
coef_df_Climate_Change <- read.csv("summaryResults3.csv")
#coef_df$trait <- factor(coef_df$X , levels=c("Fossorial","Cave-dwelling","Arboreal","Ground-dwelling","Pedal","Jumping","Brachiating"))


coef_df_Climate_Change2 <- coef_df_Climate_Change %>%pivot_longer(!X, names_to = "col_head", values_to = "coef")

coef_df_Climate_Change2<-subset(coef_df_Climate_Change2, subset = !X %in% "(Intercept)")


Climate_Change_Coef<-ggplot(coef_df_Climate_Change2, aes(x=X, y=coef)) + geom_violin(fill="#31AE3B")+
  #scale_y_continuous(breaks=seq(-.5,2.5,.5),limits = c(-.5, 2.5))  +
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 3) +
  xlab("")+ylab("")+
  scale_x_discrete(labels=c("Marine and Brackish","Marine, Brackish, and Freshwater","Total Length (log10)mm"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  theme_classic() + geom_hline(yintercept = 0,linetype="dotted")+
  theme(text = element_text(size = 14))+
  ggtitle("Climate Change")+ 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

Climate_Change_Coef
#save as object

save(Climate_Change_Coef, file = "Climate_Change_Coef.RData")



results_Climate_change <- group_by(coef_df_Climate_Change2, X) %>% summarise(Coef = median(coef), Var = var(coef)) %>% 
  rename(Variable_Level=X)

results_Climate_change$SE <- sqrt(results_Climate_change$Var)

results_Climate_change$L95 <- results_Climate_change$Coef - (1.96 * results_Climate_change$SE)
results_Climate_change$U95 <- results_Climate_change$Coef + (1.96 * results_Climate_change$SE)

write.csv(results_Climate_change, "Climate_Change_Model_Results.csv")




#Development
coef_df_Development <- read.csv("summaryResults4.csv")
#coef_df$trait <- factor(coef_df$X , levels=c("Fossorial","Cave-dwelling","Arboreal","Ground-dwelling","Pedal","Jumping","Brachiating"))


coef_df_Development2 <- coef_df_Development %>%pivot_longer(!X, names_to = "col_head", values_to = "coef")

coef_df_Development2<-subset(coef_df_Development2, subset = !X %in% "(Intercept)")


Development_Coef<-ggplot(coef_df_Development2, aes(x=X, y=coef)) + geom_violin(fill="#31AE3B")+
  scale_y_continuous(breaks=seq(-1.5,2.5,.5),limits = c(-1.5, 2.5))  +
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 3) +
  xlab("")+ylab("Extinction Risk (log-odds)")+
  scale_x_discrete(labels=c("Marine and Brackish","Marine, Brackish, and Freshwater","Total Length (log10)mm"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  theme_classic() + geom_hline(yintercept = 0,linetype="dotted")+
  theme(text = element_text(size = 14))+
  ggtitle("Development")+ 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 


Development_Coef
#save as object

save(Development_Coef, file = "Development_Coef.RData")



results_Developent<- group_by(coef_df_Development2, X) %>% summarise(Coef = median(coef), Var = var(coef)) %>% 
  rename(Variable_Level=X)

results_Developent$SE <- sqrt(results_Developent$Var)

results_Developent$L95 <- results_Developent$Coef - (1.96 * results_Developent$SE)
results_Developent$U95 <- results_Developent$Coef + (1.96 * results_Developent$SE)

write.csv(results_Developent, "Development_Model_Results.csv")



