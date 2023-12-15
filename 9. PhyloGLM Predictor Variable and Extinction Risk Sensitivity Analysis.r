#PhyloGLM Model Selection and Results


#Packages----
library(ape) #for reading phylogenetic trees
library(caper) #for loops that match species in phylo tree to species in dataset
library(phylolm) #runs PhyloGLM
library(qpcR) # for akaike.weights function
library(rr2) # for R2 function
library(vioplot) #for violin plots
library(tidyverse)
library(moments)

rm(list = ls())

#Load in Data----

setwd("Your_WD")


#read in phylogenetic data and dataset
temp=list.files(pattern="*.tre$")

mytree=lapply(temp,read.tree)

Master_Dataset <- read.csv("Master_Dataset.csv", header = TRUE, sep = ",") 


Master_Dataset$tl_mm <- as.numeric(Master_Dataset$tl_mm)

Master_Dataset <- Master_Dataset %>% 
  mutate_all(~ifelse(. == "", NA, .))

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

Master_Dataset <- Master_Dataset %>%
  mutate(Tiering = case_when(
    #Create Benthic
    !is.na(benthic) | !is.na(demersal) ~ "Benthic",
    #create pelagic
    !is.na(bathypelagic) | !is.na(pelagic.neritic) | !is.na(pelagic) | !is.na(pelagic.oceanic)
    ~ "Pelagic",
  ))    

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


#builds subset that will be used for analysis with only relevant variables and 
#ensures variables are correct format 
eco.species_wFactors<-cbind.data.frame(Master_Dataset[,c("scientific_name","IUCN_Status")], 
                                       TL_Log=as.numeric(Master_Dataset[,"tl_mmlog"]),
                                       Trophic_Level=as.numeric(Master_Dataset[,"troph_level"]),
                                       Min_Doubling=as.factor(Master_Dataset[,"resil"]),
                                       Euryhaline_Status=as.factor(Master_Dataset[,"Euryhaline_Status"]),
                                       Tiering=as.factor(Master_Dataset[,"Tiering"]),
                                       risk=as.factor(Master_Dataset[,"Risk"]))

#eco.species_wFactors <- eco.species_wFactors %>% mutate_all(na_if,"")
#omit NA's to restrict analysis to complete cases only.
eco.species_wFactors <- na.omit(eco.species_wFactors)


# 
# levels(eco.species_wFactors$Min_Doubling)[levels(eco.species_wFactors$Min_Doubling) == ""] <- "Not specified"
# 
# eco.species_wFactors <- eco.species_wFactors %>% filter(!Min_Doubling == "Not specified")

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

class(phydat[[1]]$risk) 

#Double checking NAs not being pulled in for this variable
hold<-data.frame(table(eco.species_wFactors$Min_Doubling))

#Setting reference levels
for(j in 1:length(mytree)){
  phydat[[j]]$Tiering<-relevel(phydat[[j]]$Tiering, ref="Pelagic")}

for(j in 1:length(mytree)){
  phydat[[j]]$Euryhaline_Status<-relevel(phydat[[j]]$Euryhaline_Status, ref="Marine Only")}

for(j in 1:length(mytree)){
  phydat[[j]]$Min_Doubling<-relevel(phydat[[j]]$Min_Doubling, ref="less than 15 months")}



#Start PhyloGLM Here----

model1 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model2 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model3 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model4 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model5 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model6 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model7 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model8 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model9 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model10 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results

##name the models

Modnames <- c("Null", "TL_Log", "Trophic_Level","Min_Doubling","Tiering","Euryhaline_Status","All","Three_Var","Four_Var","Best_Guess")
# 
# empty list to catch the 100 AIC results
aicResults <- vector(mode = "list", length = 100)
modelResults <- vector(mode = "list", length = 100)
bestModel_coef <- vector(mode = "list", length = 100)
aicWeights <- vector(mode = "list", length = 100)

for(j in 1: length(treedat)){
  one_tree <- treedat[[j]]$phy
  df <- phydat[[j]] 
  model1[[j]] <- phyloglm(risk~1, phy=one_tree, data=df, method= "logistic_MPLE") 
  model2[[j]] <- phyloglm(risk~TL_Log, phy=one_tree, data=df, method= "logistic_MPLE")
  model3[[j]] <- phyloglm(risk~Trophic_Level, phy=one_tree, data=df, method= "logistic_MPLE") 
  model4[[j]] <- phyloglm(risk~Min_Doubling, phy=one_tree, data=df, method= "logistic_MPLE")
  model5[[j]] <- phyloglm(risk~Tiering, phy=one_tree, data=df, method= "logistic_MPLE") 
  model6[[j]] <- phyloglm(risk~Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  model7[[j]] <- phyloglm(risk~TL_Log + Trophic_Level + Min_Doubling + Tiering + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  model8[[j]] <- phyloglm(risk~TL_Log + Euryhaline_Status + Min_Doubling, phy=one_tree, data=df, method= "logistic_MPLE")
  model9[[j]] <- phyloglm(risk~TL_Log + Tiering + Euryhaline_Status + Min_Doubling, phy=one_tree, data=df, method= "logistic_MPLE")
  model10[[j]] <- phyloglm(risk~TL_Log + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  
  allModels_j <- list(model1[[j]], model2[[j]], model3[[j]], model4[[j]], 
                      model5[[j]], model6[[j]], model7[[j]], model8[[j]], model9[[j]], model10[[j]])
  
  allSummaries_j <- lapply(allModels_j, summary)
  allAIC_j <- lapply(allModels_j, AIC)
  allR2_j <- as.matrix(unlist(c(NA, R2(model2[[j]],model1[[j]]), R2(model3[[j]],model1[[j]]), R2(model4[[j]],model1[[j]]), 
                                R2(model5[[j]],model1[[j]]), R2(model6[[j]],model1[[j]]), R2(model7[[j]],model1[[j]]),
                                R2(model8[[j]],model1[[j]]),R2(model9[[j]],model1[[j]]),R2(model10[[j]],model1[[j]]))))
  
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
  bestModel_coef[[j]] <- model10[[j]]$coefficients

  
}

# unpack the list of AIC results:
ALL_aicResults <- do.call(rbind, aicResults)

# <<< Now you should have a dataframe of 4 COLUMNS * (4*100) ROWS, with the model deltaAIC, rank, modelName, and tree denoted
# Then you an SUMMARIZE further to get the 95% confidence interval for delta AIC per model
# something like:
delAIC95 <- list(); delAIC50 <- list(); numberBest <- c(); R2_95 <- list()
for( j in 1:length(Modnames)){
  modelRes <- ALL_aicResults[which(ALL_aicResults$model==Modnames[j]),] # gets all rows for a given model
  delAIC95[[j]] <- quantile(modelRes$deltaAIC, c(0.5, 0.025, 0.975)) # gets the median and 95% CI
  #delAIC50[[j]] <- quantile(modelRes$deltaAIC, c(0.5, 0.25, 0.75)) # gets the median and 50% CI (interquartile)
  numberBest[j] <- length( modelRes[which(modelRes$deltaAIC==0),1] )
  R2_95[[j]] <- quantile( na.omit(modelRes$R2), c(0.5, 0.025, 0.975)) # gets the median and 95% CI -- for R2 values
}

sum_R2_95 <- do.call(rbind, R2_95)
colnames(sum_R2_95) <- c("R2_median", "R2_low95", "R2_up95")
sum_delAIC95 <- do.call(rbind, delAIC95)
colnames(sum_delAIC95) <- c("delAIC_median", "delAIC_low95", "delAIC_up95")

ALL_delAIC95 <- cbind.data.frame( sum_delAIC95, numberBest = numberBest, sum_R2_95)
rownames(ALL_delAIC95) <- Modnames
write.csv(ALL_delAIC95, file=paste0("deltaAIC_100trees.csv"))

# save the 100 tree results of AIC weight, deltaAIC, modelRank, and R2
write.csv(ALL_aicResults, file="modelCompare_aicWeight_delAIC_R2_100trees.csv")

write.csv(bestModel_coef,file = "summaryResults.csv")



coef_df <- read.csv("summaryResults.csv")
#coef_df$trait <- factor(coef_df$X , levels=c("Fossorial","Cave-dwelling","Arboreal","Ground-dwelling","Pedal","Jumping","Brachiating"))


coef_df2 <- coef_df %>%pivot_longer(!X, names_to = "col_head", values_to = "coef")

coef_df2<-subset(coef_df2, subset = !X %in% "(Intercept)")


coef_df2 <- coef_df2 %>%
  mutate(color_group = case_when(
    X %in% c("TL_Log") ~ "One", 
    X %in% c("Euryhaline_StatusMarine and Brackish") ~ "Two", 
    X %in% c("Euryhaline_StatusMarine Brackish and Freshwater") ~ "Two", 
    #More conditions
  ))

# Create color pallete to match other figures:
cPalette <- c("#FF7F7F",  "#7F7FFF")

#Violin Plot----


#Generate Figure 1
Best_Model_Coef<-ggplot(coef_df2, aes(x=X, y=coef)) + 
  geom_violin(aes(fill=color_group), show.legend = FALSE) +
  scale_y_continuous(breaks=seq(-.5,2,.5),limits = c(-.75, 2))  +
  scale_fill_manual(values=cPalette) +
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 3) +
  xlab("")+ylab("Extinction Risk (log-odds)")+
  scale_x_discrete(labels=c("Marine and\nBrackish",
                            "Marine, Brackish, \nand Freshwater",
                            "Total Length\n(log10)mm"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  theme_classic() + geom_hline(yintercept = 0,linetype="dotted")+
  theme(axis.text.x = element_text(angle = 0)) +
  theme(axis.text.y = element_text(color="Black",size = 15)) +
  theme(axis.title.y = element_text(color="Black",size = 20)) +
  theme(axis.text.x = element_text(color="Black",size = 15)) +
  theme(axis.title.x = element_text(color="Black",size = 20))
  

Best_Model_Coef


save(Best_Model_Coef, file = "Predictor_by_Extinction_Risk1.RData")


#Saving figure
png("Best_Model_Coef.png", width = 10, height = 10, units = "in", res = 600)
Best_Model_Coef
dev.off()


skewness(model10[[j]]$coefficients)
#data not normally distributed, so used median instead of mean > ABS 1. 



#Median coefficient and 95% intervals 
results <- group_by(coef_df2, X) %>% summarize(
  Coef = median(coef), Var_Trees = var(coef), 
  Low95_Trees = quantile(coef, probs = .025), Up95_Trees = quantile(coef, probs = .975)) %>%
  rename(Variable_Level=X)


results$SE_Trees <- sqrt(results$Var_Trees)



#Average SE Across All Trees----
#Above calculates the variance and SE around the phylogenetic models 
#Below is to calculate the variance and SE WITHIN each phylogentic model

#Ie see this output
summary(model10[[1]])

#This is the results for the first phylogenetic tree and includes an estimate
#of std error within the model around the estimate.  This std error is 
#different from the std error calculated above (std error BETWEEN models)

#Below we will calculate the averaged SE from within each model.
#pulling out only SEs
SE <- lapply(model10, '[[',3)

#Turn to data frame, add in column identifying what each row is
SE <- data.frame(SE)
SE2 <- SE %>%
  add_column(Category = c("Intercept", 
                          "Total Length", "Marine and Brackish",
                          "Marine Brackish and Freshwater"),
             .before = "c.0.664099740252251..0.261262376840463..0.343301801383354..0.24940386927388")

#Pivot
SE2 <- SE2 %>% pivot_longer(!Category, names_to = "col_head", values_to = "SE")%>% select(Category, SE)

#remove intercept, not interested in
#SE2 <- filter(SE2, !Category == "Intercept")

#TL_SE <- filter(SE2, Category == "Total Length") %>% select(Category, SE)
#mmkay setup, now need to calculate pooled variance across dataframes
#could break into seperate dataframes if needed

#Get variance by squaring SE
SE3 <- group_by(SE2, Category)%>% summarize( SE = SE, Var = SE^2)

#Total Length Averaged SE
TL_SE <-  filter(SE3, Category == "Total Length") 

TL_Avg_Variance <- mean(TL_SE$Var)

TL_SE_Avg <- sqrt(TL_Avg_Variance)

#Marine and Brackish Averaged SE
MB_SE <-  filter(SE3, Category == "Marine and Brackish") 

MB_Avg_Variance <- mean(MB_SE$Var)

MB_SE_Avg <- sqrt(MB_Avg_Variance)

#Marine Brackish and Freshwater Averaged SE
MBF_SE <-  filter(SE3, Category == "Marine Brackish and Freshwater") 

MBF_Avg_Variance <- mean(MBF_SE$Var)

MBF_SE_Avg <- sqrt(MBF_Avg_Variance)

#Append to table

Hold <- data.frame(rbind(TL_SE_Avg, MB_SE_Avg, MBF_SE_Avg))

Hold <- Hold %>% rename(SE = rbind.TL_SE_Avg..MB_SE_Avg..MBF_SE_Avg.)


results2 <- cbind(results, Hold["SE"])


results2$L95 <- results2$Coef - (results2$SE * 1.96)
results2$U95 <- results2$Coef + (results2$SE * 1.96)


#Median AIC of best model----

AIC <- data.frame(lapply(model10, '[[',11))

#Adding in dummy column to make pivot longer work
AIC <- AIC %>% add_column(Model = "Best_Model")

AIC2 <- AIC %>% pivot_longer(!Model, names_to = "col_head", values_to = "AIC")

median(AIC2$AIC)

save.image("Sensitivity_AnalysisMar14.RData")

#Sensitivity Analysis----
#Vulnerable as Non-Threatened
#label everything as 2. Stick to best model
#Load in Data----

rm(list = ls())


setwd("Your_WD")


#read in phylogenetic data and dataset
temp=list.files(pattern="*.tre$")

mytree=lapply(temp,read.tree)


Master_Dataset_Sen <- read.csv("Master_Dataset.csv", header = TRUE, sep = ",") 


Master_Dataset_Sen$tl_mm <- as.numeric(Master_Dataset_Sen$tl_mm)

Master_Dataset_Sen <- Master_Dataset_Sen %>% mutate_all(na_if,"")

#Modify Dataset----
#building binary risk categorization
Master_Dataset_Sen <- Master_Dataset_Sen %>%
  mutate(Risk = case_when(
    IUCN_Status %in% c("LC","NT", "VU") ~ "0", 
    IUCN_Status %in% c("EN", "CR") ~ "1", 
    #More conditions
  ))


#Log transforming length
Master_Dataset_Sen$tl_mmlog=log(Master_Dataset_Sen$tl_mm, 10)
#Create tiering category

Master_Dataset_Sen <- Master_Dataset_Sen %>%
  mutate(Tiering = case_when(
    #Create Benthic
    !is.na(benthic) | !is.na(demersal) ~ "Benthic",
    #create pelagic
    !is.na(bathypelagic) | !is.na(pelagic.neritic) | !is.na(pelagic) | !is.na(pelagic.oceanic)
    ~ "Pelagic",
  ))    

#Create Euryhaline status
Master_Dataset_Sen <- Master_Dataset_Sen %>%
  mutate(Euryhaline_Status = case_when(
    is.na(brackish) & is.na(freshwater) ~ "Marine Only",
    !is.na(brackish) & !is.na(freshwater)
    ~ "Marine Brackish and Freshwater",
    !is.na(brackish) ~ "Marine and Brackish",
    !is.na(freshwater) & !is.na(marine)~ "Marine and Freshwater"))    


# substitutes the space between the genus and species with an underscore to match the species names in the phylogenetic tree
Master_Dataset_Sen$scientific_name<-sub(" ", "_", Master_Dataset_Sen$scientific_name) 


#builds subset that will be used for analysis with only relevant variables and 
#ensures variables are correct format 
eco.species_wFactors_Sen<-cbind.data.frame(Master_Dataset_Sen[,c("scientific_name","IUCN_Status")], 
                                       TL_Log=as.numeric(Master_Dataset_Sen[,"tl_mmlog"]),
                                       Trophic_Level=as.numeric(Master_Dataset_Sen[,"troph_level"]),
                                       Min_Doubling=as.factor(Master_Dataset_Sen[,"resil"]),
                                       Euryhaline_Status=as.factor(Master_Dataset_Sen[,"Euryhaline_Status"]),
                                       Tiering=as.factor(Master_Dataset_Sen[,"Tiering"]),
                                       risk=as.factor(Master_Dataset_Sen[,"Risk"]))

#eco.species_wFactors <- eco.species_wFactors %>% mutate_all(na_if,"")
#omit NA's to restrict analysis to complete cases only.
eco.species_wFactors_Sen <- na.omit(eco.species_wFactors_Sen)


# 
# levels(eco.species_wFactors$Min_Doubling)[levels(eco.species_wFactors$Min_Doubling) == ""] <- "Not specified"
# 
# eco.species_wFactors <- eco.species_wFactors %>% filter(!Min_Doubling == "Not specified")

hold<-data.frame(table(eco.species_wFactors_Sen$Min_Doubling))
#Match Phylo and Dataset Species----

# make sure the species names between the eco data and tree match 
# matches the names in the eco.species dataset to the names in the tree tip.label. If there is no match, an NA is returned
sp.matching<-vector(mode = "list", length (mytree))
for(j in 1:length(mytree)){
  sp.matching[[j]]<-match(eco.species_wFactors_Sen$scientific_name, mytree[[j]]$tip.label)
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
  treedat[[j]]<-comparative.data(Trimmedtrees[[j]], eco.species_wFactors_Sen, "scientific_name", na.omit=FALSE)
}

phydat<-vector(mode = "list", length = 100)
for(j in 1:length(mytree)){
  #phydat[[j]]<-data.frame(treedat[[j]]$data)
  phydat[[j]]<- treedat[[j]]$data # don't want to make this a data.frame again -> will destroy the "factor" coding
}

class(phydat[[1]]$risk) 

#Double checking NAs not being pulled in for this variable
hold<-data.frame(table(eco.species_wFactors_Sen$Min_Doubling))

#Setting reference levels
for(j in 1:length(mytree)){
  phydat[[j]]$Tiering<-relevel(phydat[[j]]$Tiering, ref="Pelagic")}

for(j in 1:length(mytree)){
  phydat[[j]]$Euryhaline_Status<-relevel(phydat[[j]]$Euryhaline_Status, ref="Marine Only")}

for(j in 1:length(mytree)){
  phydat[[j]]$Min_Doubling<-relevel(phydat[[j]]$Min_Doubling, ref="less than 15 months")}



#Start PhyloGLM Here----

model1 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model2 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model3 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model4 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model5 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model6 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model7 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model8 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model9 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results
model10 <- vector(mode = "list", length = 100) # make an empty list of length 100 to catch results

##name the models to be averaged

Modnames <- c("Null", "TL_Log", "Trophic_Level","Min_Doubling","Tiering","Euryhaline_Status","All","Three_Var","Four_Var","Best_Guess")
# 
# empty list to catch the 100 AIC results
aicResults <- vector(mode = "list", length = 100)
modelResults <- vector(mode = "list", length = 100)
bestModel_coef <- vector(mode = "list", length = 100)
aicWeights <- vector(mode = "list", length = 100)

for(j in 1: length(treedat)){
  one_tree <- treedat[[j]]$phy
  df <- phydat[[j]] 
  model1[[j]] <- phyloglm(risk~1, phy=one_tree, data=df, method= "logistic_MPLE") 
  model2[[j]] <- phyloglm(risk~TL_Log, phy=one_tree, data=df, method= "logistic_MPLE")
  model3[[j]] <- phyloglm(risk~Trophic_Level, phy=one_tree, data=df, method= "logistic_MPLE") 
  model4[[j]] <- phyloglm(risk~Min_Doubling, phy=one_tree, data=df, method= "logistic_MPLE")
  model5[[j]] <- phyloglm(risk~Tiering, phy=one_tree, data=df, method= "logistic_MPLE") 
  model6[[j]] <- phyloglm(risk~Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  model7[[j]] <- phyloglm(risk~TL_Log + Trophic_Level + Min_Doubling + Tiering + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  model8[[j]] <- phyloglm(risk~TL_Log + Euryhaline_Status + Min_Doubling, phy=one_tree, data=df, method= "logistic_MPLE")
  model9[[j]] <- phyloglm(risk~TL_Log + Tiering + Euryhaline_Status + Min_Doubling, phy=one_tree, data=df, method= "logistic_MPLE")
  model10[[j]] <- phyloglm(risk~TL_Log + Euryhaline_Status, phy=one_tree, data=df, method= "logistic_MPLE")
  
  allModels_j <- list(model1[[j]], model2[[j]], model3[[j]], model4[[j]], 
                      model5[[j]], model6[[j]], model7[[j]], model8[[j]], model9[[j]], model10[[j]])
  
  allSummaries_j <- lapply(allModels_j, summary)
  allAIC_j <- lapply(allModels_j, AIC)
  allR2_j <- as.matrix(unlist(c(NA, R2(model2[[j]],model1[[j]]), R2(model3[[j]],model1[[j]]), R2(model4[[j]],model1[[j]]), 
                                R2(model5[[j]],model1[[j]]), R2(model6[[j]],model1[[j]]), R2(model7[[j]],model1[[j]]),
                                R2(model8[[j]],model1[[j]]),R2(model9[[j]],model1[[j]]),R2(model10[[j]],model1[[j]]))))
  
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
  bestModel_coef[[j]] <- model10[[j]]$coefficients
  
  
}

# unpack the list of AIC results:
ALL_aicResults <- do.call(rbind, aicResults)

# <<< Now you should have a dataframe of 4 COLUMNS * (4*100) ROWS, with the model deltaAIC, rank, modelName, and tree denoted
# Then you an SUMMARIZE further to get the 95% confidence interval for delta AIC per model
# something like:
delAIC95 <- list(); delAIC50 <- list(); numberBest <- c(); R2_95 <- list()
for( j in 1:length(Modnames)){
  modelRes <- ALL_aicResults[which(ALL_aicResults$model==Modnames[j]),] # gets all rows for a given model
  delAIC95[[j]] <- quantile(modelRes$deltaAIC, c(0.5, 0.025, 0.975)) # gets the median and 95% CI
  #delAIC50[[j]] <- quantile(modelRes$deltaAIC, c(0.5, 0.25, 0.75)) # gets the median and 50% CI (interquartile)
  numberBest[j] <- length( modelRes[which(modelRes$deltaAIC==0),1] )
  R2_95[[j]] <- quantile( na.omit(modelRes$R2), c(0.5, 0.025, 0.975)) # gets the median and 95% CI -- for R2 values
}

sum_R2_95 <- do.call(rbind, R2_95)
colnames(sum_R2_95) <- c("R2_median", "R2_low95", "R2_up95")
sum_delAIC95 <- do.call(rbind, delAIC95)
colnames(sum_delAIC95) <- c("delAIC_median", "delAIC_low95", "delAIC_up95")

ALL_delAIC95 <- cbind.data.frame( sum_delAIC95, numberBest = numberBest, sum_R2_95)
rownames(ALL_delAIC95) <- Modnames
write.csv(ALL_delAIC95, file=paste0("deltaAIC_100trees_Sen.csv"))

# save the 100 tree results of AIC weight, deltaAIC, modelRank, and R2
write.csv(ALL_aicResults, file="modelCompare_aicWeight_delAIC_R2_100trees_Sen.csv")

write.csv(bestModel_coef,file = "summaryResults_Sen.csv")



coef_df <- read.csv("summaryResults_Sen.csv")
#coef_df$trait <- factor(coef_df$X , levels=c("Fossorial","Cave-dwelling","Arboreal","Ground-dwelling","Pedal","Jumping","Brachiating"))


coef_df2 <- coef_df %>%pivot_longer(!X, names_to = "col_head", values_to = "coef")

coef_df2<-subset(coef_df2, subset = !X %in% "(Intercept)")

#Violin Plot----

load("Sensitivity_AnalysisMar14.RData")
#Create column to group colorings by
coef_df2 <- coef_df2 %>%
  mutate(color_group = case_when(
    X %in% c("TL_Log") ~ "One", 
    X %in% c("Euryhaline_StatusMarine and Brackish") ~ "Two", 
    X %in% c("Euryhaline_StatusMarine Brackish and Freshwater") ~ "Two", 
    #More conditions
  ))

# Create color pallete to match other figures:
cPalette <- c("#FF7F7F",  "#7F7FFF")

#Generate Figure 1
Best_Model_Coef_Sensitivity<-ggplot(coef_df2, aes(x=X, y=coef)) + 
  geom_violin(aes(fill=color_group), show.legend = FALSE)+
  scale_y_continuous(breaks=seq(-.5,2,.5),limits = c(-.75, 2))  +
  scale_fill_manual(values=cPalette) +
  #scale_fill_brewer(palette="Set1") +
  #scale_fill_manual(values=c("#CC6666", "#CC6666", "#66CC99")) +
  stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 3) +
  xlab("")+
  ylab("")+
  scale_x_discrete(labels=c("Marine and\nBrackish",
                            "Marine, Brackish, \nand Freshwater",
                            "Total Length\n(log10)mm"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  theme_classic() + geom_hline(yintercept = 0,linetype="dotted")+
  theme(axis.text.x = element_text(angle = 0)) +
  theme(axis.text.y = element_text(color="Black",size = 15)) +
  theme(axis.title.y = element_text(color="Black",size = 20)) +
  theme(axis.text.x = element_text(color="Black",size = 15)) +
  theme(axis.title.x = element_text(color="Black",size = 20))


Best_Model_Coef_Sensitivity


save(Best_Model_Coef_Sensitivity, file = "Best_Model_Coef_Sensitivity1.rdata")



save.image("Sensitivity_AnalysisMar14.RData")


save(Best_Model_Coef_Sensitivity, file = "Best_Model_Coef_Sensitivity.rdata")
#Saving figure
png("Best_Model_Coef_Sen.png", width = 10, height = 10, units = "in", res = 600)
Best_Model_Coef
dev.off()


skewness(model10[[j]]$coefficients)
#data not normally distributed, so used median instead of mean > ABS 1. 



#Median coefficient and 95% intervals 
results <- group_by(coef_df2, X) %>% summarize(
  Coef = median(coef), Var_Trees = var(coef), 
  Low95_Trees = quantile(coef, probs = .025), Up95_Trees = quantile(coef, probs = .975)) %>%
  rename(Variable_Level=X)


results$SE_Trees <- sqrt(results$Var_Trees)



#Average SE Across All Trees----
#Above calculates the variance and SE around the phylogenetic models 
#Below is to calculate the variance and SE WITHIN each phylogentic model

#Ie see this output
summary(model10[[1]])

#This is the results for the first phylogenetic tree and includes an estimate
#of std error within the model around the estimate.  This std error is 
#different from the std error calculated above (std error BETWEEN models)

#Below we will calculate the averaged SE from within each model.
#pulling out only SEs
SE <- lapply(model10, '[[',3)

#Turn to data frame, add in column identifying what each row is
SE <- data.frame(SE)
SE2 <- SE %>%
  add_column(Category = c("Intercept", 
                          "Total Length", "Marine and Brackish",
                          "Marine Brackish and Freshwater"),
             .before = "c.1.03369162827528..0.377120662971914..0.470396766467823..0.384960420496807")

#Pivot
SE2 <- SE2 %>% pivot_longer(!Category, names_to = "col_head", values_to = "SE")%>% select(Category, SE)

#remove intercept, not interested in
#SE2 <- filter(SE2, !Category == "Intercept")

#TL_SE <- filter(SE2, Category == "Total Length") %>% select(Category, SE)
#mmkay setup, now need to calculate pooled variance across dataframes
#could break into seperate dataframes if needed

#Get variance by squaring SE
SE3 <- group_by(SE2, Category)%>% summarize( SE = SE, Var = SE^2)

#Total Length Averaged SE
TL_SE <-  filter(SE3, Category == "Total Length") 

TL_Avg_Variance <- mean(TL_SE$Var)

TL_SE_Avg <- sqrt(TL_Avg_Variance)

#Marine and Brackish Averaged SE
MB_SE <-  filter(SE3, Category == "Marine and Brackish") 

MB_Avg_Variance <- mean(MB_SE$Var)

MB_SE_Avg <- sqrt(MB_Avg_Variance)

#Marine Brackish and Freshwater Averaged SE
MBF_SE <-  filter(SE3, Category == "Marine Brackish and Freshwater") 

MBF_Avg_Variance <- mean(MBF_SE$Var)

MBF_SE_Avg <- sqrt(MBF_Avg_Variance)

#Append to table

Hold <- data.frame(rbind(TL_SE_Avg, MB_SE_Avg, MBF_SE_Avg))

Hold <- Hold %>% rename(SE = rbind.TL_SE_Avg..MB_SE_Avg..MBF_SE_Avg.)


results2 <- cbind(results, Hold["SE"])


results2$L95 <- results2$Coef - (results2$SE * 1.96)
results2$U95 <- results2$Coef + (results2$SE * 1.96)


write.csv(results2, "Results_Sen.csv")

#Median AIC of best model----

AIC <- data.frame(lapply(model10, '[[',11))

#Adding in dummy column to make pivot longer work
AIC <- AIC %>% add_column(Model = "Best_Model")

AIC2 <- AIC %>% pivot_longer(!Model, names_to = "col_head", values_to = "AIC")

median(AIC2$AIC)


