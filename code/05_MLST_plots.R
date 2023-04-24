# Author: Mary Jewell
# Date: 3/10/2023
# Notes: 
# Last updated: 3/10/2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~

rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("C:/Users/marye/OneDrive/Desktop/Thesis/Data")

full_df <- read.csv("meta_mlst.csv")

##~~~~~~~~~~~~~~~~
# MLST Plots ####
##~~~~~~~~~~~~~~~~

### Plot of top STs by isolation source
# Filter dataset for ST types with > 10 isolates
st_subset <- full_df %>% filter(ST %in% c("131","73","-","95","12","69","372","127",
                                          "10","21","11","227","58","48","101","17",
                                          "297","491","162","405","5337"))

# Create ancillary dataset of # of isolates per ST type
group_labels <- st_subset %>% group_by(ST) %>% summarise(n = n())

# Stacked barplot of top STs
st_plot2 <- ggplot(data = st_subset, aes(x = fct_infreq(ST), fill = source_cat)) + 
  geom_bar(stat = "count", position = "stack") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Number of Isolates of the Most Frequent ST Types") +
  geom_text(aes(ST, n, label = n, fill = NULL), data = group_labels, vjust = -0.5) + #labels above bars
  theme(axis.title = element_text(size = 20)) + #increase axis label size
  theme(plot.title = element_text(size = 25)) + #increase title size
  theme(legend.text = element_text(size = 15)) + #increase legend size
  theme(legend.title = element_text(size = 15))
st_plot2

length(which(full_df$ST == "-" & full_df$source_cat == "Environment"))


st_plot3 <- ggplot(data = st_subset, aes(x = fct_infreq(ST), fill = source_cat)) + 
  geom_bar(stat = "count", position = "stack") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Number of Isolates of the Most Frequent ST Types") +
  geom_text(stat = 'count', aes(label = after_stat(count)), #labels in each stacked bar
            position = position_stack(vjust = 0.5))
st_plot3


### Plot of top STs faceted by isolation source

## Set up source-specific dataframes
animal <- as.data.frame(st_tab_animal[order(st_tab_animal, decreasing = T)])
colnames(animal) <- c("ST", "Freq")

environment <- as.data.frame(st_tab_environmental[order(st_tab_environmental, decreasing = T)])
colnames(environment) <- c("ST", "Freq")

food <- as.data.frame(st_tab_food[order(st_tab_food, decreasing = T)])
colnames(food) <- c("ST", "Freq")

human <- as.data.frame(st_tab_human[order(st_tab_human, decreasing = T)])
colnames(human) <- c("ST", "Freq")



## Free axes:
# Animal
st_plot_animal <- ggplot(data = head(animal, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "salmon") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Animals")
st_plot_animal


# Environmental
st_plot_environment <- ggplot(data = head(environment, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "forestgreen") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Environmental Sources")
st_plot_environment

# Food
st_plot_food <- ggplot(data = head(food, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "dodgerblue") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Food")
st_plot_food

# Human
st_plot_human <- ggplot(data = head(human, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "darkorchid1") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Humans")
st_plot_human


facet_plot <- ggarrange(st_plot_animal, st_plot_environment, 
                        st_plot_food, st_plot_human, ncol = 2, nrow = 2)
facet_plot



## Scaled on the same axes:
# Animal scaled
st_plot_animal_scaled <- ggplot(data = head(animal, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "salmon") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Animals") + ylim(0, 150)
st_plot_animal_scaled


# Environmental scaled
st_plot_environment_scaled <- ggplot(data = head(environment, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "forestgreen") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Environmental Sources") + ylim(0, 150)
st_plot_environment_scaled

# Food scaled
st_plot_food_scaled <- ggplot(data = head(food, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "dodgerblue") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Food") + ylim(0, 150)
st_plot_food_scaled

# Human scaled
st_plot_human_scaled <- ggplot(data = head(human, 12), aes(x = ST, y = Freq)) + 
  geom_bar(stat = "identity", fill = "darkorchid1") + theme_minimal() +
  labs(x = "ST Type", y = "Number of Isolates", fill = "Isolation Source",
       title = "Isolates from Humans") + ylim(0, 150)
st_plot_human_scaled

# Arrange plots
facet_plot_scales <- ggarrange(st_plot_animal, st_plot_environment, 
                               st_plot_food, st_plot_human, ncol = 2, nrow = 2)
facet_plot_scales

