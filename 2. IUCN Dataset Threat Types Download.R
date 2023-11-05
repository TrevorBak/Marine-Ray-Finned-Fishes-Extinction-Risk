##Pull IUCN Threat Types Using API and Pivot to Wide Binary Format##



library(rredlist)
library(tidyverse)


rm(list = ls())


setwd("Your_WD")

#need to pull in a dataset that has all the species in it that you will want 
#threat type information for.  Code will query the IUCN database by each species.
#Species name should be in genus_SpecificEpithat (Homo_Sapien).

df<-read.csv("Master_Dataset.csv")


#creates a blank list the length of the number of species in the dataset.  In my
#dataset species column is titled scientific_name.

result <- vector('list', length(df$scientific_name))

#Loop that queries the dataset through the IUCN API, pulling threat types for each species
#in the dataset.  Key is a unique ID token assigned by the IUCN for a user that gives 
#them permission to access the dataset with the API.  Can be requested here:
#https://apiv3.iucnredlist.org/api/v3/token

for (i in df$scientific_name) {
  result[[i]] <- rl_threats(name=i, key = '1234', parse = TRUE)
}



#*Note that the code will need to run uninterrupted to successfully download
#threat types for all species.  I've had issues where, due to 
#drops in internet connection, the code will stop running and a partial list 
#will be generated.  Can mitigate this problem by running code in small chunks
#by reducing down dataset into say, quarters, and running each quarter at a time
#Can do with df<-df[1:2501,] which would cut a 10,000 species dataset
#into 2500.  2500 can be ran, then process repeated for each quarter.*


#A bit unwieldy but the 'result' list generated will be full of empty elements 
#for the first half and full of threat type information for the second half.
#So for say a dataset of 10,000 species the resulting list will look like:

#1-10,000 rows - empty
#10,001 - 20,000 - threat types

#So dataset needs to be cut to only contain the rows with threat type information
#This code manually does that with the 10,000 species example.  Adjust to size of your dataset.

result <- result[10018:19931]


#right now result is a nested list which is not very workable.  This converts
#result to a dataframe

df2 <- purrr::map_df(result, `[[`, 'result', .id = 'species')



##The following code is to to take the nested list output and convert it into 
#a binary set of variables indicating whether the species is impacted by each
#general threat (general threat type means the titles in bold here:
#https://www.iucnredlist.org/resources/threat-classification-scheme )

#Note the API code output also includes specific subthreats, severity,
#timing, etc. which could also be explored but is not explored in this code.




#We only want the general threat type number which is the first part of the 
#X.X.X sequence.  This line breaks apart the code by period.  We will then
#just use the first number.

df2<-within(df2, Gen_Code<-data.frame(do.call('rbind', strsplit(as.character(df2$code), '.', fixed=TRUE))))


#above code creates a dataframe within a dataframe, this code puts everything back into one dataframe
#and removes extra columns that don't contain the general code for clarity.
df2<-do.call(data.frame, df2)
df2 <- subset(df2, select = -c(Gen_Code.X2, Gen_Code.X3))
str(df2)

#reordering columns so general code is next to species name.  
df2 <- df2[, c(2,1,9,4,3,5,6,7,8)]

#issue - some species/general codes are duplicated.  This is because the API output
#can list sub-threat types along with more general, ie tiger 5.4 AND tiger 5.4.1 AND tiger 5.4.2
#When general code is extracted this leaves tiger 5 tiger 5 tiger 5.  Only want one
#instance of tiger 5.  This code does that and shrinks dataset to just species
#and general code.

df2 <- unique( df2[ , 2:3 ] )


#Start of conversion of this dataset to a binomial dataset with each threat type
#as a column. 

#Start by renaming general codes to their actual titles.
df2 <- df2%>%
  mutate(Threat_Type = case_when(
    Gen_Code.X1 %in% c("1") ~ "Development", 
    Gen_Code.X1 %in% c("2") ~ "Aquaculture_and_Agriculture",  
    Gen_Code.X1 %in% c("3") ~ "Energy_Production_and_Mining", 
    Gen_Code.X1 %in% c("4") ~ "Transportation", 
    Gen_Code.X1 %in% c("5") ~ "Harvesting", 
    Gen_Code.X1 %in% c("6") ~ "Human_Disturbance", 
    Gen_Code.X1 %in% c("7") ~ "Natural_System_Modification", 
    Gen_Code.X1 %in% c("8") ~ "Invasive_Species_and_Disease", 
    Gen_Code.X1 %in% c("9") ~ "Pollution", 
    Gen_Code.X1 %in% c("10") ~ "Geological_Events", 
    Gen_Code.X1 %in% c("11") ~ "Climate_Change", 
    Gen_Code.X1 %in% c("12") ~ "Other", 
  ))

#remove general code, will now only use titles

df2 <- subset(df2, select = -c(Gen_Code.X1))

#pivots to binary and wider format

df2<-df2 %>% mutate(val=1) %>% pivot_wider(names_from = Threat_Type,values_from = val) %>% 
  mutate(across(-species, ~replace_na(.x, 0))) %>%
  mutate(across(-species, ~ifelse(.x==1, 1,0)))

#issue - species listed on multiple rows, want 1 row per species.  This code collapses to 1 row

df2<-aggregate(.~species, df2, max)

#Now have a set of species with binary threat type in wide format.

write.csv(df2, "Threats.csv")

#Import master dataset to join together with threats
Master <- read.csv("Master_Dataset.csv")

df2 <- df2 %>% rename(scientific_name = species)

Master_Dataset <- left_join(Master, df2, by = "scientific_name")

write.csv(Master_Dataset, "Master_Dataset.csv")
