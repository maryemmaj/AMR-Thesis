# Author: Mary Jewell
# Date: 3/15/2023
# Notes: Resistance analysis, including proportion tests and visualization of resistance data
# Last updated: 4/14/2023
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
library(tidyquant) #moving average in ggplot

## Data ####
wide <- fread("Data/resfinder_wide.csv")
long <- read.csv("Data/resfinder_long_cleaned.csv")
full <- fread("Data/full_analytic_data.csv")

long <- long %>%
  filter(!is.na(class))

full <- full %>%
  filter(!is.na(source_cat))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Resistance exploratory analysis tables ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Binary (none/any) resistance across isolation sources
table(full$source_cat, full$res_any)
prop.table(table(full$source_cat, full$res_any), margin = 1)

table(full$res_any)
prop.table(table(full$res_any))

# Binary (3+) resistance across isolation sources
(tab1 <- table(full$source_cat, full$res_binary))
prop.table(table(full$source_cat, full$res_binary), margin = 1)

table(full$res_binary)
prop.table(table(full$res_binary))

# Categorical (0-2, 3-9, 10+) resistance across isolation sources
(tab2 <- table(full$source_cat, full$res_cat))
prop.table(table(full$source_cat, full$res_cat), margin = 1)

table(full$res_cat)
prop.table(table(full$res_cat))


full$res_sum <- as.numeric(full$res_sum)
summary(full$res_sum[full$source_cat == "Animal"])
summary(full$res_sum[full$source_cat == "Environment"])
summary(full$res_sum[full$source_cat == "Food"])
summary(full$res_sum[full$source_cat == "Human"])


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hypothesis testing: Compare ARGs between isolation sources ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Proportion of 3+ resistance across categories
resist_vector <- tab1[,2] #select the total in each source cat with 3+ resistance
total_vector <- tab1[,1] + tab1[,2] #select the total in each source cat
prop.test(resist_vector, total_vector) # p-value < 2.2e-16

(cbind(resist_vector, total_vector))

# Proportion of 10+ resistance across categories
resist10_vector <- tab2[,3] #select the total in each source cat with 10+ resistance
prop.test(resist10_vector, total_vector) # p-value < 2.2e-16


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Hypothesis testing: Compare ARGs by isolation source within ST type ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Top 4 STs: ST131, ST73, Unidentified, ST95

# Difference in proportion of 3+ resistance across categories within each ST
## ST131
(tab3 <- table(full$source_cat[full$ST == 131], full$res_binary[full$ST == 131]))
fisher.test(tab3) # p = 1

## ST73
(tab4 <- table(full$source_cat[full$ST == 73], full$res_binary[full$ST == 73]))
fisher.test(tab4) # p = 0.03757

## Unidentified ST
(tab5 <- table(full$source_cat[full$ST == "-"], full$res_binary[full$ST == "-"]))
fisher.test(tab5) # p = 0.1995

## ST95
(tab6 <- table(full$source_cat[full$ST == 95], full$res_binary[full$ST == 95]))
fisher.test(tab6) # p = 0.4111


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots of resistance categories ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

binary_plot <- ggplot(data = full, aes(x = as.factor(res_binary), # bars are binary res categories
                                       fill = source_cat)) +  # colored by isolation source
  geom_bar(position = "stack") + theme_minimal() + 
  labs(x = "Predicted Resistance", y = "Number of Isolates", fill = "Isolation Source",
       title = "Distribution of Isolates Across Binary Resistance Categories") +
  scale_x_discrete(labels = c("< 3 Antimicrobials", "3+ Antimicrobials")) + #x-axis labels
  theme(axis.title = element_text(size = 15)) + #increase axis label size
  theme(plot.title = element_text(size = 20)) + #increase title size
  theme(axis.text = element_text(size = 15)) + #increase axis text size
  theme(legend.text = element_text(size = 15)) + #increase legend size
  theme(legend.title = element_text(size = 15))

binary_plot

cat_plot <- ggplot(data = full, aes(x = as.factor(res_cat), # bars are 3-category res categories
                                    fill = source_cat)) + # colored by isolation source
  geom_bar(position = "stack") + theme_minimal() + 
  labs(x = "Predicted Resistance", y = "Number of Isolates", fill = "Isolation Source",
       title = "Distribution of Isolates Across Three Resistance Categories") +
  scale_x_discrete(labels = c("< 3 Antimicrobials", "3-9 Antimicrobials", "10+ Antimicrobials")) +
  theme(axis.title = element_text(size = 15)) + #increase axis label size
  theme(plot.title = element_text(size = 20)) + #increase title size
  theme(axis.text = element_text(size = 15)) + #increase axis text size
  theme(legend.text = element_text(size = 15)) + #increase legend size
  theme(legend.title = element_text(size = 15))
cat_plot

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## By Number of antimicrobials ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Histogram of number of resistant isolates vs number of antimicrobials
res_histogram <- ggplot(data = full, aes(x = res_sum, # continuous x: number of resistance genes
                                         fill = source_cat)) + # colored by source
  geom_histogram(bins = 15) + theme_minimal() +
  labs(x = "Number of Antimicrobials", y = "Number of Isolates", fill = "Isolation Source") +
  theme(axis.title = element_text(size = 20)) + #increase axis label size
  theme(plot.title = element_text(size = 25)) + #increase title size
  theme(legend.text = element_text(size = 15)) + #increase legend size
  theme(legend.title = element_text(size = 15))
res_histogram


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## By Antimicrobial ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Number resistant per antimicrobial
antimicrobial_counts <- long %>%
  group_by(antimicrobial) %>% #get the number of isolates resistant to each antimicrobial
  summarize(num_resistant = sum(pheno_pred == "Resistant")) %>% 
  arrange(desc(num_resistant)) #sort descending by the number of isolates resistant

# Take the antimicrobials with > 100 isolates resistant
antimicrobial_highcounts <- antimicrobial_counts %>%
  filter(num_resistant > 100)

# add a variable of antibiotic abbreviations to use in plots
antimicrobial_highcounts$abv <- c("AMX", "AMP", "PIP", "TIC", "STR", "CEF", "SXT", 
                                  "CIP", "TMP", "DOX", "TET", "NAL", "SPT", "MIN", 
                                  "ERY", "AZM", "SPM", "TEL")

# Most frequent antimicrobials in general (those with > 100 isolates resistant)
bar1 <- ggplot(data = antimicrobial_highcounts, 
               aes(x = reorder(abv, -num_resistant), #order top antimicrobials in descending order
                   y = num_resistant)) + #plot counts of number of resistant isolates
  geom_bar(stat = "identity") + theme_minimal() +
  labs(x = "Antimicrobial", y = "Number of Resistant Isolates",
       title = "Antimicrobials with Greater than \n 100 Resistant Isolates") +
  theme(axis.title = element_text(size = 20)) + #increase axis label size
  theme(plot.title = element_text(size = 25)) + #increase title size
  theme(axis.text = element_text(size = 10)) #increase axis text size
bar1

# Get names of the top antimicrobials
top_antimicrobials <- head(antimicrobial_highcounts$antimicrobial, 20)

# Assign abbreviation for top antimicrobials
long$antimicrobial_abv[long$antimicrobial == "amoxicillin"] <- "AMX"
long$antimicrobial_abv[long$antimicrobial == "ampicillin"] <- "AMP"
long$antimicrobial_abv[long$antimicrobial == "piperacillin"] <- "PIP"
long$antimicrobial_abv[long$antimicrobial == "ticarcillin"] <- "TIC"
long$antimicrobial_abv[long$antimicrobial == "streptomycin"] <- "STR"
long$antimicrobial_abv[long$antimicrobial == "cephalothin"] <- "CEF"
long$antimicrobial_abv[long$antimicrobial == "sulfamethoxazole"] <- "SXT"
long$antimicrobial_abv[long$antimicrobial == "ciprofloxacin"] <- "CIP"
long$antimicrobial_abv[long$antimicrobial == "trimethoprim"] <- "TMP"
long$antimicrobial_abv[long$antimicrobial == "doxycycline"] <- "DOX"
long$antimicrobial_abv[long$antimicrobial == "tetracycline"] <- "TET"
long$antimicrobial_abv[long$antimicrobial == "nalidixic acid"] <- "NAL"
long$antimicrobial_abv[long$antimicrobial == "spectinomycin"] <- "SPT"
long$antimicrobial_abv[long$antimicrobial == "minocycline"] <- "MIN"
long$antimicrobial_abv[long$antimicrobial == "erythromycin"] <- "ERY"
long$antimicrobial_abv[long$antimicrobial == "azithromycin"] <- "AZM"
long$antimicrobial_abv[long$antimicrobial == "spiramycin"] <- "SPM"
long$antimicrobial_abv[long$antimicrobial == "telithromycin"] <- "TEL"

# Create ancillary dataset of # of isolates per antimicrobial for bar labels
group_labels <- long[long$pheno_pred == "Resistant" & long$antimicrobial %in% top_antimicrobials,] %>% 
  group_by(antimicrobial_abv) %>% summarise(n = n())

bar2 <- ggplot(data = long[long$pheno_pred == "Resistant" & #select isolates resistant to top antimicrobials
                           long$antimicrobial %in% top_antimicrobials,], 
               aes(x = fct_infreq(antimicrobial_abv), #sort in descending order by count of isolates
                   fill = source_cat)) + #color by isolation source
  geom_bar() + theme_minimal() +
  labs(x = "Antimicrobial", y = "Number of Isolates Predicted Resistant",
       title = "Most Common Antimicrobial Resistance Genes Detected",
       fill = "Isolation Source") +
  geom_text(aes(antimicrobial_abv, n, label = n, fill = NULL), 
            data = group_labels, vjust = -0.5) + #labels from ancillary counts dataset above bars
  theme(axis.title = element_text(size = 15)) + #increase axis label size
  theme(plot.title = element_text(size = 15)) + #increase title size
  theme(axis.text = element_text(size = 10)) #increase axis text size
bar2


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## By Antimicrobial Class ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create ancillary dataset of # of isolates per antimicrobial for bar labels
group_labels <- long[long$pheno_pred == "Resistant",] %>% 
  group_by(class) %>% summarise(n = n())

# Bar plot of resistance distribution across antimicrobial class
bar3 <- ggplot(data = long[long$pheno_pred == "Resistant",], #select resistant isolates
               aes(x = fct_infreq(class), #sort in descending order by count of isolates
                   fill = source_cat)) + #color by isolation source
  geom_bar() + theme_minimal() +
  labs(x = "Antimicrobial Class", y = "Number of Isolates Predicted Resistant",
       title = "Most Common Classes of Antimicrobial Resistance Gene Detected",
       fill = "Isolation Source") +
  geom_text(aes(class, n, label = n, fill = NULL), 
            data = group_labels, vjust = -0.5) + #labels above bars from ancillary counts dataset
  theme(axis.text.x = element_text(angle = 30, hjust=1)) + #adjust x-axis to be readable
  theme(axis.title = element_text(size = 15)) + #increase axis label size
  theme(plot.title = element_text(size = 15)) + #increase title size
  theme(axis.text = element_text(size = 10)) #increase axis text size
bar3


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Trends of resistance over time ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Based on feedback from Peter:

### Creating datasets by year ####
# Time series by isolation source: proportion resistant in each isolation source in each year
time_series <- full %>%
  group_by(collection_year, source_cat) %>% #by isolation source and year
  summarize(prop_resistant = sum(res_any == 1) / n(), #proportion of resistant isolates
            num_resistant = n()) #number of resistant isolates


# Time series overall: proportion resistant overall in each year
time_series_overall <- full %>%
  group_by(collection_year) %>% #by isolation source and year
  summarize(prop_resistant = sum(res_any == 1) / n(), #proportion of resistant isolates
            num_resistant = n()) #number of resistant isolates


### Plotting proportion resistant ####
# proportion of resistance over time by isolation source
timeplot_source_prop <- ggplot(data = time_series, 
                   aes(x = collection_year, y = prop_resistant, 
                       group = source_cat, color = source_cat)) + 
  geom_line(size = 1) + theme_minimal() +
  labs(x = "Collection Year", y = "Proportion of Resistant Isolates", color = "Isolation Source",
       title = "Proportion of Isolates With Any Antimicrobial Resistance Over Time")
timeplot_source_prop


timeplot_overall_prop <- ggplot(data = time_series_overall, 
                   aes(x = collection_year, y = prop_resistant)) + 
  geom_line(size = 1) + theme_minimal() +
  geom_ma(ma_fun = SMA, n = 5, color = "red", size = 2) + #5-year moving average
  labs(x = "Collection Year", y = "Proportion of Isolates Resistant", 
       title = "Proportion of Isolates With Any Antimicrobial Resistance Over Time",
       subtitle = "Plotted With Five Year Moving Average")
timeplot_overall_prop


### Plotting number resistant ####
# Overall
timeplot_overall_num <- ggplot(data = time_series_overall, 
                           aes(x = collection_year, y = num_resistant)) + 
  geom_line(size = 1) + theme_minimal() + xlim(1985, 2022) + #xlim to cut off sparse data before 1980
  geom_ma(ma_fun = SMA, n = 5, color = "red") + #five-year moving average
  labs(x = "Collection Year", y = "Number of Isolates Resistant", 
       title = "Number of Isolates With Any Antimicrobial Resistance Over Time",
       subtitle = "Plotted With Five Year Moving Average")
timeplot_overall_num

# By isolation source
timeplot_source_num <- ggplot(data = time_series, 
                   aes(x = collection_year, y = num_resistant, 
                       group = source_cat, color = source_cat)) + 
  geom_line(size = 1) + theme_minimal() + 
  xlim(1985, 2022) + #xlim to cut off sparse data before 1980
  ylim(0, 250) + #pick a more appropriate ylim for disaggregated data
  labs(x = "Collection Year", y = "Number of Resistant Isolates", color = "Isolation Source",
       title = "Number of Isolates With Any Antimicrobial Resistance Over Time")
timeplot_source_num


