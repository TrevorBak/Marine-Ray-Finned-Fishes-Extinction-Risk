#Downloading list of species and IUCN extinction risk using IUCN API.

library(tidyverse)
library(jsonlite)
library(curl)


setwd("Your_WD")


#IUCN Download----
#Initialize data frame, redList will catch results from API request

redList <- data.frame()

#Token given out by IUCN as a password to access their API.  Request here: https://apiv3.iucnredlist.org/api/v3/token
token <- "1234"
spCount <- 100
i <- 0

while(spCount > 0) {
  species <- fromJSON(paste("http://apiv3.iucnredlist.org/api/v3/species/page/",i,"?token=",token, sep=""))	
  spCount <- species$count
  if(spCount > 0) {
    redList <- rbind(redList, species$result)
  }
  print(i)
  i <- i + 1
}

#Subset to animals then Actinopterygii
animalia <- subset(redList, kingdom_name == 'ANIMALIA')
Actinopterygii <- subset(animalia, class_name == 'ACTINOPTERYGII')


#remove excess columns not used in analyses for conciseness.

Actinopterygii <- Actinopterygii %>% select(scientific_name, category)

#Some duplicates of species that need to be removed


Actinopterygii <- Actinopterygii %>% filter(duplicated(scientific_name) == FALSE)



#Join with WoRMs----
#Join IUCN dataset with a dataset of taxonomy from the World Register of Marine Species (WoRMs)
#This will filter out any non-marine Actinopterygii species.

#CSV that I import was downloaded by requesting a copy of the database here: https://www.marinespecies.org/documents.php


WORMS <- read.csv("WoRMs.csv")

#Filter WORMS subset for species only in animal kingdom, documented at species level, and accepted taxonomy
#Also rename column label for species so that it matches column label for species in IUCN dataset 
WORMS <- WORMS %>% rename(scientific_name = scientificName) %>%
  filter(kingdom == "Animalia", taxonRank == "Species", taxonomicStatus == "accepted")


#remove excess columns not used in analyses for conciseness.

WORMS <- WORMS %>% select(kingdom, phylum, class, order, family, genus, scientific_name)

#Perform inner join to only preserve the species in both dataset (I.e. have IUCN status AND are marine)

IUCN_Marine <- inner_join(WORMS, Actinopterygii, by = "scientific_name")


#Combine with Fishbase.org data.  See PHP code for creating Fishbase data.  

#Join with Fishbase data----

FB <- read.csv("Fishbase.csv")

#Rename species column in Fishbase dataset to match species column in IUCN_Marine. Also remove
#taxon information as this is redundant with taxon information already in IUCN dataset.

FB <- FB %>% rename(scientific_name = taxon_name) %>%
  select(-c(phylum, class, order, family, genus, fb_taxon_id, taxon_uri))

#Perform inner join, want to keep only species in both dataset

Master_Dataset <- inner_join(IUCN_Marine, FB, by = "scientific_name")

#Check IUCN RedList assignments

table(Master_Dataset$category)

#LR/lc and LR/nt are outdated terminology for lower risk, least conern and lower risk, near threatend
#update to modern categorizations (for more see https://nc.iucnredlist.org/redlist/resources/files/1530881462-rl_criteria_1994_versus_2001.pdf)

Master_Dataset <- Master_Dataset %>% mutate(IUCN_Hold = case_when(
  category == "LR/lc" ~ "LC",
  category == "LR/nt" ~ "NT")) %>% 
  mutate(IUCN_Status = coalesce(IUCN_Hold, category)) %>%
  select(!c(category, IUCN_Hold))
  
#Filtering out Extinct and Extinct in the wild species

Master_Dataset <- Master_Dataset %>% filter(!c(IUCN_Status == "EW" | IUCN_Status == "EX"))

table(Master_Dataset$IUCN_Status)

#Write to csv for threat types script

write.csv(Master_Dataset, "Master_Dataset.csv")



