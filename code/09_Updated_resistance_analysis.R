# Author: Mary Jewell
# Date: 4/23/2023
# Notes: Resistance analysis, updated after Scott's feedback
# Last updated: 4/23/2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~~

rm(list = ls())
setwd("C:/Users/marye/OneDrive/Desktop/Thesis")

library(data.table)
library(dplyr)
library(ggplot2) #plotting
library(forcats) #reorder bars in ggplot
library(ggpubr) #ggarrange multiple plots

## Data ####
wide <- fread("Data/resfinder_wide.csv")
long <- read.csv("Data/resfinder_long_cleaned.csv")
full <- fread("Data/full_analytic_data.csv")

long <- long %>%
  filter(!is.na(class))

full <- full %>%
  filter(!is.na(source_cat))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Resistance exploratory analysis table ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Table 1: include numbers to show distribution of isolation year
# Create pre/post 2000 isolation year variable
full$year_cat <- ifelse(full$collection_year >= 2000, 1, 0)
table(full$year_cat, full$source_cat, useNA = "always")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Resistance category pie charts ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

table(full$source_cat, full$res_cat)


## Animal
animal <- full %>% filter(source_cat == "Animal")

animal_pie <- ggplot(data = animal, aes(x = source_cat, fill = as.factor(res_cat))) +
  geom_bar(position = "stack", width=1) + # create stacked bar plot
  coord_polar("y", start = 0) + # flip to pie chart with coord_polar
  labs(title = "Animal") + guides(fill="none") + #remove legend
  scale_fill_manual(values = c("#9db6ff", "#5c6ab4", "#1c266d" )) +
  theme_void() # remove all numbers & axes
animal_pie


## Environment
environment <- full %>% filter(source_cat == "Environment")

environment_pie <- ggplot(data = environment, aes(x = source_cat, fill = as.factor(res_cat))) +
  geom_bar(position = "stack", width=1) + # create stacked bar plot
  coord_polar("y", start = 0) + # flip to pie chart with coord_polar
  labs(title = "Environment") + guides(fill="none") + #remove legend
  scale_fill_manual(values = c("#9db6ff", "#5c6ab4", "#1c266d" )) +
  theme_void() # remove all numbers & axes
environment_pie


## Food
food <- full %>% filter(source_cat == "Food")

food_pie <- ggplot(data = food, aes(x = source_cat, fill = as.factor(res_cat))) +
  geom_bar(position = "stack", width=1) + # create stacked bar plot
  coord_polar("y", start = 0) + # flip to pie chart with coord_polar
  labs(title = "Food") + guides(fill="none") + #remove legend
  scale_fill_manual(values = c("#9db6ff", "#5c6ab4", "#1c266d" )) +
  theme_void() # remove all numbers & axes
food_pie


## Human
human <- full %>% filter(source_cat == "Human")

# Order factor levels to create ordered legens
human$res_cat <- factor(human$res_cat, ordered = T, levels = c(1:3),
                        labels = c("0-2 Antimicrobials",
                                   "3-9 Antimicrobials",
                                   "10+ Antimicrobials"))

human_pie <- ggplot(data = human, aes(x = source_cat, fill = as.factor(res_cat))) +
  geom_bar(position = "stack", width=1) + # create stacked bar plot
  coord_polar("y", start = 0) + # flip to pie chart with coord_polar
  labs(title = "Human", fill = "Predicted Resistance") + # Keep legend 
  scale_fill_manual(values = c("#9db6ff", "#5c6ab4", "#1c266d" )) +
  theme_void() # remove all numbers & axes
human_pie

# Arrange all four pie charts
ggarrange(animal_pie, environment_pie, food_pie, human_pie,
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom") #One legend for all plots


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bar plot alternative to pies ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create summary dataframe of proportions in each res_cat
bar_data <- full %>%
  group_by(source_cat) %>%
  summarize(
    prop_res_cat1 = mean(res_cat == 1) * 100,
    prop_res_cat2 = mean(res_cat == 2) * 100,
    prop_res_cat3 = mean(res_cat == 3) * 100
  )

# Melt to long form: row per source_cat and res_cat
bar_data <- melt(bar_data, id.vars = c("source_cat"), variable.name = "res_cat")

# Rename columns
names(bar_data) <- c("source_cat", "res_cat", "perc")

# Clean text from res_cat variable
bar_data$res_cat <- gsub("prop_res_cat", "", bar_data$res_cat)

# Order factor levels to create ordered legens
bar_data$res_cat <- factor(bar_data$res_cat, ordered = T, levels = c(3:1),
                           labels = c("10+ Antimicrobials",
                                      "3-9 Antimicrobials",
                                      "0-2 Antimicrobials"))

# Stacked bar plot of resistance categories in each source
ggplot(bar_data, aes(x = source_cat, y = perc, fill = res_cat)) + 
  geom_bar(stat = "identity") + theme_minimal() +
  labs(x = "Isolation Source", y = "Percentage of Isolates", 
       fill = "Predicted Resistance") +
  scale_fill_manual(values = c("#1c266d", "#5c6ab4", "#9db6ff"))



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Creating resistance profiles ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# amx/amp, cephtriaxone/cephtazidine, cipro, gentamycin, streptomycin, sxt, tetracycline

names(full)
alphnames <- sort(names(full))
alphnames

res_profile <- full %>% select(SRR, source_cat, collection_year, ST,
                               # These antibiotics will make up the profile:
                               amoxicillin, ampicillin, ceftazidime, 
                               ciprofloxacin, gentamicin, streptomycin, 
                               sulfamethoxazole, tetracycline,
                               # Also include other resistance data:
                               res_sum, res_binary, res_cat, res_any)

# Indicator 1 if isolate is resistant to either amoxicillin or ampicillin
res_profile$indicator1_amx_amp <- ifelse(res_profile$amoxicillin == "Resistant" |
                                        res_profile$ampicillin == "Resistant",
                                        1, 0)

# Indicator 2 if isolate is resistant to ceftazidime
res_profile$indicator2_caz <- ifelse(res_profile$ceftazidime == "Resistant", 1, 0)

# Indicator 3 if isolate is resistant to ciprofloxacin
res_profile$indicator3_cip <- ifelse(res_profile$ciprofloxacin == "Resistant", 1, 0)

# Indicator 4 if isolate is resistant to gentamicin
res_profile$indicator4_gene <- ifelse(res_profile$gentamicin == "Resistant", 1, 0)

# Indicator 5 if isolate is resistant to streptomycin
res_profile$indicator5_str <- ifelse(res_profile$streptomycin == "Resistant", 1, 0)

# Indicator 6 if isolate is resistant to sulfamethoxazole
res_profile$indicator6_sxt <- ifelse(res_profile$sulfamethoxazole == "Resistant", 1, 0)

# Indicator 7 if isolate is resistant to tetracycline
res_profile$indicator7_tet <- ifelse(res_profile$tetracycline == "Resistant", 1, 0)

# Create one 7-digit profile
res_profile$profile <- paste0(res_profile$indicator1_amx_amp,
                              res_profile$indicator2_caz,
                              res_profile$indicator3_cip,
                              res_profile$indicator4_gene,
                              res_profile$indicator5_str,
                              res_profile$indicator6_sxt,
                              res_profile$indicator7_tet)

length(unique(res_profile$profile)) # There are 68 unique profiles

# Table of profiles in decreasing order
profile_tab <- table(res_profile$profile)
profile_tab[order(profile_tab, decreasing = T)]


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualizing resistance profiles ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Create summary dataframe of proportions in each res_cat
summarized_profile_data <- res_profile %>%
  group_by(source_cat) %>%
  summarize( #come back to this
    prop_res_cat1 = mean(res_cat == 1) * 100,
    prop_res_cat2 = mean(res_cat == 2) * 100,
    prop_res_cat3 = mean(res_cat == 3) * 100
  )


##~~~~~~~~~
# Heatmap
##~~~~~~~~~

heatmap_data <- res_profile %>% select(SRR, indicator1_amx_amp:indicator7_tet)
names(heatmap_data) <- c("SRR", "Ampicillin/Amoxicillin",
                         "Ceftazidime", "Ciprofloxacin", "Gentamicin",
                         "Streptomycin", "Sulfamethoxazole", "Tetracycline")

# Melt the data frame to long format
heatmap_data <- melt(heatmap_data, id.vars = "SRR")


# Create a heatmap using ggplot2
ggplot(heatmap_data, aes(x = variable, y = SRR, fill = value)) +
  geom_tile() +   theme_bw() + coord_flip() + #flip axes to make readable
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(axis.text.x = element_blank(),  #remove x axis labels
        axis.ticks.x = element_blank(), #remove x axis ticks
        axis.title = element_blank(), #remove axis titles
        legend.position = "none") #remove legend


