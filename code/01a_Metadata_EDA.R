# Author: Mary Jewell
# Date: 12/27/2022
# Notes: Performing exploratory data analysis on cleaned metadata.
#         Uses cleaned-metadata5.csv
# Last updated: 1/27/2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~

rm(list=ls())
setwd("C:/Users/marye/OneDrive/Desktop/Thesis/Data")

library(tidyverse)
library(ggplot2)
library(Amelia) #describe missingness

##~~~~~~~~~~~~~~~~~~~~~~~
# Some exploration ####
##~~~~~~~~~~~~~~~~~~~~~~~

meta <- read.csv("C:/Users/marye/OneDrive/Desktop/Thesis/Data/cleaned-metadata-final.csv")
names(meta)
table(meta$source_type)

#isolates by year
plot1 <- ggplot(meta, aes(x = collection_year)) + geom_histogram(bins = 30) +
  labs(title = "Number of Isolates Collected by Year", x = "Collection Year",
       y = "Number of Isolates") + theme(text = element_text(size = 20))
plot1

#isolates by year and source type
plot2 <- ggplot(meta, 
                aes(x = collection_year, fill = source_type)) + 
  geom_histogram(bins = 30) + xlim(1983, 2022) +
  labs(title = "Number of Isolates Collected by Year", x = "Collection Year",
       y = "Number of Isolates", fill = "Isolation Source") + 
  theme(text = element_text(size = 20))
plot2

table(meta$isolation_notes)

#Describe and visualize missingness
missing_data <- meta %>% select("organism", "source_type", "isolation_notes", 
                                "collection_year","strain", "serotype2")
missmap(missing_data, main = "Missingness Map of Select Variables", legend = F)

sapply(meta, function(x) sum(is.na(x)))
