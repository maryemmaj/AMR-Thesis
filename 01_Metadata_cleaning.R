# Author: Mary Jewell
# Date: 10/3/2022
# Notes: Cleaning NCBI metadata from txt form to produce useable metadata.
#       There is an easier way to get metadata from NCBI using rentrez 
#       (see the file biosample_extraction.R)
#       Rentrez produces an XML, which is hard to clean. This script cleans the text file.
# Last updated: 12/27/2022
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~

rm(list=ls())
library(stringr)
library(tidyverse)
library(naniar) #fixing NAs
library(zoo) #for na.locf function
library(janitor) #cleaning column names
library(lubridate) #cleaning date values
library(linelist) #please let this clean the dates
library(ggplot2) #data exploration/visualization

setwd("C:/Users/marye/OneDrive/Desktop/Thesis/Data")
wd <- "C:/Users/marye/OneDrive/Desktop/Thesis"

#Read data
biosample <- read.table(paste0(wd,"Data/biosample_search_result_full.txt"),
                        header = F, sep = "\t")

##~~~~~~~~~~~~~~~~~~~~~~~~
# Managing text file ####
##~~~~~~~~~~~~~~~~~~~~~~~~

#Split attribute rows
biosample[c('attribute', 'value')] <- str_split_fixed(biosample$V1, "=", 2)

#Split non-attribute rows
biosample[c('attribute2', 'value2')] <- str_split_fixed(biosample$attribute, ":", 2)

#Convert blank spaces to NA
biosample$value <- ifelse(biosample$value == "", NA, biosample$value)
biosample$value2 <- ifelse(biosample$value2 == "", NA, biosample$value2)

#Coalesce value and value2 into one column
biosample$value <- coalesce(biosample$value, biosample$value2)

biosample2 <- biosample %>% select(attribute2, value)

##~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to wide form ####
##~~~~~~~~~~~~~~~~~~~~~~~~~

#Generate unique ID column
biosample2$id <- as.numeric(biosample2$attribute2) #pull just the ID numbers
biosample2$id <- na.locf(biosample2$id) #fill in NAs with last non-NA value

#Clean attributes2 column
biosample3 <- subset(biosample2, attribute2 != "Attributes" & attribute2 != "Identifiers")
biosample3 <- subset(biosample3, nchar(attribute2) != 1 & nchar(attribute2) != 2 &
                       nchar(attribute2) != 3 & nchar(attribute2) != 4)

#Pivot-wider
biosample4 <- biosample3 %>% pivot_wider(
  id_cols = id,
  names_from = attribute2,
  values_from = value
)


#remove / from colnames
colnames(biosample4) <- gsub("/", "", colnames(biosample4))

names(biosample4)
biosample4 <- clean_names(biosample4, "snake") #clean & standardize names

bio5 <- biosample4 %>% select(1:27) #select the attributes that seem most useful/populated


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean covariate data types ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

names(bio5)

bio5 <- bio5 %>%
  mutate(across(everything(), as.character))

##~~~~~~~~~~~~~~~~~~~
# Geographic location
##~~~~~~~~~~~~~~~~~~~
table(bio5$geographic_location)
bio5$geographic_location <- "USA:WA" #standardize geo loc

##~~~~~~~~
# Organism
##~~~~~~~~
table(bio5$organism) #all e. coli, needs to be cleaned slightly
#trim leading/trailing whitespace
bio5$organism <- trimws(bio5$organism, which = c("both"), whitespace = "[ \t\r\n]") 
bio6 <- subset(bio5, bio5$organism == "Escherichia coli") #select non-O157 samples

##~~~~~~~~~~~~~~~
# Collection date
##~~~~~~~~~~~~~~~
table(bio6$collection_date) #some are from before Y2K...
#Three different date formats: year, ymd, ym:
#ymd
bio6$collection_year <- as.character(year(ymd(bio6$collection_date)))
table(bio6$collection_year)
#ym
bio6$year2 <- as.character(year(ym(bio6$collection_date)))
table(bio6$year2)
#year
bio6$year3 <- ifelse(nchar(bio6$collection_date) == 4, bio6$collection_date, NA)
table(bio6$year3)

#Create one year variable
bio6$collection_year <- coalesce(bio6$collection_year, bio6$year2, bio6$year3)
bio6$collection_year <- as.Date(bio6$collection_year) #error
table(bio6$collection_year, useNA = "always")

bio6$after2000 <- ifelse(bio6$collection_year < 2000, 0, 1)
cbind(head(bio6$after2000), head(bio6$collection_year)) #check that year matches
table(bio6$after2000)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolation source - basic cleaning
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(bio6$source_type)
table(bio6$isolation_source[is.null(bio6$source_type)]) 
#there are no iso source where source type is null

#trim leading/trailing whitespace
bio6$source_type <- trimws(bio6$source_type, which = c("both"), whitespace = "[ \t\r\n]")
#standardize text case
bio6$source_type <- tolower(bio6$source_type)
#standardize source categories
bio6$source_type[bio6$source_type == "animal feed"] <- "animal"
bio6$source_type[bio6$source_type == "Animal feed"] <- "animal"
bio6$source_type[bio6$source_type == 'c("animal", "animal")'] <- "animal"
bio6$source_type[bio6$source_type == 'c("food", "food")'] <- "food"

#Isolation source by year dichotomy
table(bio6$source_type, bio6$after2000)


write.csv(bio6, "cleaned-metadata.csv", row.names = F)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean isolation sources ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())

bio6 <- read.csv("cleaned-metadata.csv")
table(bio6$isolation_source[bio6$source_type == "animal"])
table(bio6$isolation_source[bio6$source_type == "environmental"])
table(bio6$isolation_source[bio6$source_type == "food"])
table(bio6$isolation_source[bio6$source_type == "human"])

table(bio6$source_type)


#If isolation source includes "animal" -> source type = "animal"
bio6$source_type[grepl("Animal", bio6$isolation_source)] <- "animal" #uppercase
bio6$source_type[grepl("animal", bio6$isolation_source)] <- "animal" #lowercase

#If isolation source includes specific animal, set source type "animal" 
bio6$source_type[grepl("Cow", bio6$isolation_source)] <- "animal" #cow
bio6$source_type[grepl("Canis", bio6$isolation_source)] <- "animal" #canine
bio6$source_type[grepl("Canine", bio6$isolation_source)] <- "animal" #canine

#wild animals?
bio6$source_type[grepl("ape", bio6$isolation_source)] <- "wild animal"
bio6$source_type[grepl("lion", bio6$isolation_source)] <- "wild animal"
bio6$source_type[grepl("gorilla", bio6$isolation_source)] <- "wild animal"
bio6$source_type[grepl("orangutan", bio6$isolation_source)] <- "wild animal"
bio6$source_type[grepl("cougar", bio6$isolation_source)] <- "wild animal"
bio6$source_type[grepl("marmoset", bio6$isolation_source)] <- "wild animal"
bio6$source_type[grepl("elephant", bio6$isolation_source)] <- "wild animal"
bio6$source_type[grepl("giraffe", bio6$isolation_source)] <- "wild animal"

#wtf is a cloacae?
miss_cloacae <- bio6 %>% filter(grepl("cloacae", bio6$isolation_source))
bio6$source_type[grepl("cloacae", bio6$isolation_source)] <- "wild animal" #duck cloacae


#Patient -> human
#filter rows with "patient"
miss_patient <- bio6 %>% filter(grepl("patient", bio6$isolation_source)) 
#all patients are Homo sapiens
table(miss_patient$host)
bio6$source_type[grepl("patient", bio6$isolation_source)] <- "human"


#Beef -> food
bio6$source_type[grepl("beef", bio6$isolation_source)] <- "food"
bio6$source_type[grepl("Beef", bio6$isolation_source)] <- "food"


#Misc environmental sources
bio6$source_type[grepl("water", bio6$isolation_source)] <- "environmental"
bio6$source_type[grepl("food sample", bio6$isolation_source)] <- "food"


#Categorize by host tags
bio6$source_type[grepl("Homo sapiens", bio6$host)] <- "human"

bio6$source_type[grepl("horse", bio6$host)] <- "animal"
bio6$source_type[grepl("cat", bio6$host)] <- "animal"
bio6$source_type[grepl("chicken", bio6$host)] <- "animal"
bio6$source_type[grepl("cattle", bio6$host)] <- "animal"
bio6$source_type[grepl("cow", bio6$host)] <- "animal"
bio6$source_type[grepl("dog", bio6$host)] <- "animal"
bio6$source_type[grepl("Bos taurus", bio6$host)] <- "animal"
bio6$source_type[grepl("pig", bio6$host)] <- "animal"
bio6$source_type[grepl("Gallus gallus domesticus", bio6$host)] <- "animal" #chicken

bio6$source_type[grepl("Duck", bio6$host)] <- "wild animal"
bio6$source_type[grepl("giraffe", bio6$host)] <- "wild animal"
bio6$source_type[grepl("marmoset", bio6$host)] <- "wild animal"
bio6$source_type[grepl("lion", bio6$host)] <- "wild animal"
bio6$source_type[grepl("leopard", bio6$host)] <- "wild animal"
bio6$source_type[grepl("orangutan", bio6$host)] <- "wild animal"
bio6$source_type[grepl("cougar", bio6$host)] <- "wild animal"
bio6$source_type[grepl("Neogale vison", bio6$host)] <- "wild animal" #mink


#End result: many isolates characterized
table(bio6$source_type)
table(bio6$isolation_source[bio6$source_type == "Unknown"])

miss <- bio6 %>% filter(source_type == "Unknown") 
# nothing I can see to do about these last 5


write.csv(bio6, "cleaned-metadata2.csv", row.names = F)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Trimming NA columns ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())

#Read in csv, converting na_strings to NA
bio6 <- read.csv("cleaned-metadata2.csv", na.strings = c("NULL", "null", "missing",
                                                         "Missing", "unknown", "Unknown",
                                                         "Not collected", "not applicable",
                                                         "Not applicable", "not collected"))

colSums(is.na(bio6))
table(bio6$source_type)

#Remove extraneous year columns
bio6$year2 <- NULL
bio6$year3 <- NULL

#There are 187 isolates with missing info for collection date
table(bio6$collection_year, useNA = "always")
table(bio6$collection_date)

#Remove columns with very high missingness, very little usefulness
bio6$package <- NULL
bio6$contact <- NULL
bio6$project_accession <- NULL
bio6$project_accession <- NULL
bio6$public_accession <- NULL
bio6$project_name <- NULL
bio6$ontological_term <- NULL
bio6$isolate <- NULL
bio6$latitude_and_longitude <- NULL

#Don't need genus or species - we have organism for every isolate
bio6$genus <- NULL
bio6$species <- NULL

###~~~~~~~~~~~
## REMOVE ALL O157 ISOLATES

table(bio6$serovar)

#Keep serovar values if not O157
bio7 <- subset(bio6, serovar != "O157:H7" | is.na(serovar))
#Check that it worked
table(bio7$serovar)
names(bio7)
table(bio7$source_type)

#Keep serotype values if not O157
bio7 <- subset(bio6, serotype != "O157:H7" | is.na(serotype))
#Check that it worked
table(bio7$serotype)
names(bio7)
table(bio7$source_type)


write.csv(bio7, "cleaned-metadata3.csv", row.names = F)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cleaning animal vs wild animal ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
bio8 <- read.csv("cleaned-metadata3.csv")

animal <- bio8 %>% filter(source_type == "animal" | source_type == "wild animal")

table(bio8$source_type, useNA = "always")

write.csv(animal, "animal_test2.csv", row.names = F)

# Convert some animal to wild animal
bio8$source_type[grepl("Odocoileus virginianus leucurus", bio8$host)] <- "wild animal" #deer
bio8$source_type[grepl("Accipitridae", bio8$host)] <- "wild animal" #hawk
bio8$source_type[grepl("mink", bio8$host)] <- "wild animal" #mink
bio8$source_type[grepl("duck", bio8$host)] <- "wild animal"
bio8$source_type[grepl("goose", bio8$host)] <- "wild animal"
bio8$source_type[grepl("Sus scrofa", bio8$isolation_source)] <- "wild animal" #wild boar


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Coalesce isolation source/host ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Adding common names to some host values
bio8$isolation_source[bio8$isolation_source == "Other - List in comments (Canis lupus)"] <-
  "Canis lupus"
bio8$host[bio8$host == "Neogale vison"] <- "Neogale vison (mink)"
bio8$host[bio8$host == "Myotis lucifugus"] <- "Myotis lucifugus (bat)"
bio8$host[bio8$host == "Odocoileus virginianus leucurus"] <- 
  "Odocoileus virginianus leucurus (deer)"
bio8$host[bio8$host == "Gallus gallus domesticus"] <- 
  "Gallus gallus domesticus (chicken)"


bio9 <- transform(bio8, 
                  isolation_notes = ifelse(is.na(isolation_source) | is.na(host), 
                                              coalesce(isolation_source, host), 
                                              paste(host, isolation_source)))

table(bio9$isolation_notes)

#Cleaning stray isolation_notes values
bio9$isolation_notes[bio9$isolation_notes == "Gorilla (gorilla) gorilla"] <- 
  "Gorilla"
bio9$isolation_notes[bio9$isolation_notes == "Macaca nigra (Celebes ape) celebese ape (Macaca nigra)" | 
                       bio9$isolation_notes == "Macaca nigra (Celebes ape) celebe ape (Macaca nigra)" |
                       bio9$isolation_notes == "Macaca nigra (Celebes ape) ape"] <- "Macaca nigra (ape)"

bio9$isolation_notes[bio9$isolation_notes == "mink mustelidae mink lung"] <- 
  "Mink lung"
bio9$isolation_notes[bio9$isolation_notes == "mink intestine mink mustelidae"] <- 
  "Mink intestine"

bio9$isolation_notes[bio9$isolation_notes == "Giraffa camelopardalis (giraffe) giraffa camelopardalis"] <- 
  "Giraffa camelopardalis (giraffe)"

write.csv(bio9, "cleaned-metadata4.csv", row.names = F)

##~~~~~~~~~~~~~~~~~~~~~~~
# Reordering columns ####
##~~~~~~~~~~~~~~~~~~~~~~~

names(bio9)

bio9$serotype2 <- coalesce(bio9$serotype, bio9$serovar)

#Select only the serotypes that are not O157
bio9 <- subset(bio9, serotype2 == "O103" | serotype2 == "O26" | is.na(serotype2))

table(bio9$serotype)
table(bio9$source_type)
names(bio9)

bio9$serotype <- NULL
bio9$serovar <- NULL


col_order <- c("id", "accession", "organism", "geographic_location", "source_type", #most important variables for ID
               "isolation_notes", "collection_year", "isolation_source", "host",
               "collected_by", "sequenced_by", "strain", "attribute_package",
               "isolate_name_alias", "ifsac_category", "serotype2", "host_disease",      
                "collection_date")

bio10 <- bio9[, col_order]
names(bio10)

table(bio10$serotype2, useNA = "always")

bio10$source_type[bio10$isolation_source == "cheddar cheese"] <- NA
bio10$isolation_notes[bio10$isolation_source == "cheddar cheese"] <- NA

write.csv(bio10, "cleaned-metadata5.csv", row.names = F)

table(bio10$source_type, useNA = "always")



##~~~~~~~~~~~~~~~~~~~~~~~
# Some exploration ####
##~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
bio10 <- read.csv("cleaned-metadata5.csv")

#isolates by year
plot1 <- ggplot(bio10, aes(x = collection_year)) + geom_histogram(bins = 30) +
  labs(title = "Number of Isolates Collected by Year", x = "Collection Year",
       y = "Number of Isolates") + theme(text = element_text(size = 20))
plot1

#isolates by year and source type
plot2 <- ggplot(bio10, 
                aes(x = collection_year, fill = source_type)) + 
  geom_histogram(bins = 30) + xlim(1983, 2022) +
  labs(title = "Number of Isolates Collected by Year", x = "Collection Year",
       y = "Number of Isolates") + theme(text = element_text(size = 20))
plot2

table(bio10$isolation_notes)
