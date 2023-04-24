# Author: Mary Jewell
# Date: 3/6/2023
# Notes: 
# Last updated: 3/9/2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~

rm(list=ls())

## Libraries
library(dplyr)
library(janitor) #clean names
library(ggplot2) #data visualization
library(forcats) #sort barplot in ggplot2
library(ggpubr) #arrange multiple ggplot objects

setwd("C:/Users/marye/OneDrive/Desktop/Thesis/Data")
wd <- "C:/Users/marye/OneDrive/Desktop/Thesis/"

# Read data
# accessions <- read.csv("ecoli_wa_us_EnteroBase_strains.csv") #accession numbers: 1590
# accessions still needs to be moved over from the cluster
meta <- read.csv("enterobase_metadata.csv", na.strings = c("", "NA")) #metadata: 1590 isolates
mlst <- read.csv("mlst.csv", header = T) # MLST analysis results


##~~~~~~~~~~~~~~~~~~~~~~~
# Cleaning metadata ####
##~~~~~~~~~~~~~~~~~~~~~~~

meta <- clean_names(meta, "snake") #standardize names
names(meta)

### Create 4-category source variable
meta$source_cat[meta$source_niche == "Human"] <- "Human"
meta$source_cat[meta$source_niche == "Food"] <- "Food"
meta$source_cat[meta$source_niche == "Environment"] <- "Environment"
meta$source_cat[meta$source_niche %in% c("Animal Feed", "Companion Animal",
                                          "Livestock", "Poultry", "Wild animal", 
                                          "Wild Animal")] <- "Animal"
table(meta$source_cat, useNA = "always")


### Drop meaningless variables 
colMeans(is.na(meta))*100

meta$status <- NULL #useless assembly status info
meta$status_2 <- NULL #useless assembly status info
meta$collection_time <- NULL #completely missing
meta$serological_group <- NULL #completely missing
meta$species_1 <- NULL #redundant


### Clean location data
meta$district[meta$district == "King County"] <- "King"
table(meta$district)

meta$city[meta$city == "Seattle, WA"] <- "Seattle"
table(meta$city)

temp_set1 <- meta %>% filter(!is.na(latitude) & !is.na(longitude))
table(temp_set1$city, useNA = "always")
table(temp_set1$district, useNA = "always")
# there is one set of lat/long where city & district are NA


### Clean pathogen/disease data

meta <- transform(meta, 
                  disease_notes = ifelse(is.na(disease) | is.na(simple_patho), 
                                           coalesce(disease, simple_patho), 
                                           paste(disease, simple_patho)))
# very limited info on disease/pathogenicity

## COME BACK TO THIS ####

table(meta$disease_notes, useNA = "always")

temp_set1 <- meta %>% filter(!is.na(disease_notes)) #why do these accessions not start with SRR?
head(meta$accession)

temp_set2 <- filter(meta, !grepl("SRR", accession)) # 45 accessions don't start with SRR
temp_set3 <- filter(mlst, !grepl("SRR", SRR))

### Create shorter metadata set for analysis

meta_short <- meta %>% select(accession, sample_id, experiment_accession, species, 
                              region, district, city, collection_year, source_cat, 
                              serotype, path_nonpath, disease_notes)
meta_short_names <- c("SRR","SAMN", "SRX", "species", "state", "county", "city", "collection_year",
                      "source_cat", "serotype", "path", "disease_notes")
colnames(meta_short) <- meta_short_names

write.csv(meta_short, "meta_short.csv", row.names = F)


##~~~~~~~~~~~~~~~~~~~~~~~
# Cleaning MLST Data ####
##~~~~~~~~~~~~~~~~~~~~~~~

# Merge metadata and mlst data ####
full_df <- merge(meta_short, mlst, by = "SRR")

# Exclude all isolates marked as O157
full_df <- full_df[!grepl("O157", full_df$serotype),]
table(full_df$serotype, useNA = "always")

# Exclude all isolates with missing source_cat
full_df <- full_df[!is.na(full_df$source_cat),]
table(full_df$source_cat, useNA = "always")

write.csv(full_df, "meta_mlst.csv", row.names = F)


### Table of isolates by source_cat ####
table(full_df$source_cat)
prop.table(table(full_df$source_cat))

# Temporal diversity
table(full_df$collection_year, useNA = "always")
length(which(full_df$collection_year >= 2000)) # 934 isolates after 2000
length(which(full_df$collection_year < 2000)) # 97 isolates before 2000
length(which(is.na(full_df$collection_year))) # 418 isolates with no collection year

# Geographic location
table(full_df$county, useNA = "always")

# Serotype/pathogenicity
table(full_df$serotype, useNA = "always")
table(full_df$disease_notes, useNA = "always")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Housekeeping Gene Frequency ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

length(unique(full_df$ST))

### Number of alleles per housekeeping gene ####
length(unique(full_df$adk))
length(unique(full_df$fumC))
length(unique(full_df$gyrB))
length(unique(full_df$icd))
length(unique(full_df$mdh))
length(unique(full_df$purA))
length(unique(full_df$recA))

### Most frequent allele for each housekeeping gene ####
names(which.max(table(full_df$adk)))
length(which(full_df$adk == "adk(6)")) #343

names(which.max(table(full_df$fumC)))
length(which(full_df$fumC == "fumC(4)")) #182

names(which.max(table(full_df$gyrB)))
length(which(full_df$gyrB == "gyrB(19)")) #216

names(which.max(table(full_df$icd)))
length(which(full_df$icd == "icd(13)")) #238

names(which.max(table(full_df$mdh)))
length(which(full_df$mdh == "mdh(17)")) #262

names(which.max(table(full_df$purA)))
length(which(full_df$purA == "purA(8)")) #300

names(which.max(table(full_df$recA)))
length(which(full_df$recA == "recA(2)")) #198

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Most common ST types overall ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Overall
st_tab <- table(full_df$ST)
(st_tab_order <- as.data.frame(st_tab[order(st_tab, decreasing = T)])) 
# 66 are type "-" what does this mean?

table(st_tab_order$Freq)
sum(st_tab_order$Freq)

# Animal
st_tab_animal <- table(full_df$ST[full_df$source_cat == "Animal" &
                                  full_df$ST %in% c("131", "73","-","95", "12", "69")])
(st_tab_animal[order(st_tab_animal, decreasing = T)])

# Environmental
st_tab_env <- table(full_df$ST[full_df$source_cat == "Environmental" &
                               full_df$ST %in% c("131", "73","-","95", "12", "69")])
(st_tab_env[order(st_tab_env, decreasing = T)])

# Food
st_tab_food <- table(full_df$ST[full_df$source_cat == "Food" &
                                full_df$ST %in% c("131", "73","-","95", "12", "69")])
(st_tab_food[order(st_tab_food, decreasing = T)])

# Human
st_tab_human <- table(full_df$ST[full_df$source_cat == "Human" &
                                 full_df$ST %in% c("131", "73","-","95", "12", "69")])
(st_tab_human[order(st_tab_human, decreasing = T)])


## Most common ST types within each isolation source ####
# Animal
st_tab_animal <- table(full_df$ST[full_df$source_cat == "Animal"])
(st_tab_animal[order(st_tab_animal, decreasing = T)])

# Environmental
st_tab_environmental <- table(full_df$ST[full_df$source_cat == "Environment"])
(st_tab_environmental[order(st_tab_environmental, decreasing = T)])

# Food
st_tab_food <- table(full_df$ST[full_df$source_cat == "Food"])
(st_tab_food[order(st_tab_food, decreasing = T)])

# Human
st_tab_human <- table(full_df$ST[full_df$source_cat == "Human"])
(st_tab_human[order(st_tab_human, decreasing = T)])


##~~~~~~~~~~~~~~~~~~~~~~~
# Hypothesis testing ####
##~~~~~~~~~~~~~~~~~~~~~~~

# Select top six ST types?
top_sts <- full_df %>% filter(ST %in% c("131", "73","-","95", "12", "69"))

table(top_sts$ST, top_sts$source_cat)

# Null: there is no relationship between isolation source and ST type
# Chi-square is not appropriate: use Fisher's
chisq.test(top_sts$ST, top_sts$source_cat)$expected

fisher.test(top_sts$ST, top_sts$source_cat, simulate.p.value = TRUE)
# simulate.p.value = simulate p-value using Monte Carlo simulation in tables greater than 2x2
# p = 0.0004998

