# Author: Mary Jewell
# Date: 11/21/2022
# Notes: Power calculations for multiple proportions test and chi-square tests
# Last updated: 12/20/2022
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~

rm(list=ls())

library(pwr)
library(EnvStats)
setwd("C:/Users/marye/OneDrive/Desktop/Thesis")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prop test power calculations ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#I will use a proportion test to evaluate the differences in the proportion of samples 
# carrying ARGs from each isolation source. 

#this will be the number of cases (isolates carrying ARGs) from each source
# for now, I'm just making up a number...
case_vector <- c(40, 50, 25, 10) 
#this will be the total number of isolates from each source
total_vector <- c(167, 264, 125, 100)

#prop.test will compare the proportions across each source (df = 3)
prop.test(case_vector, total_vector)
# in this case, with made-up data, the p-value = 0.046 with chi-square value = 7.99 and df = 3


## Calculating power ####
# for prop tests:
## h = 0.2 is a small effect size, h = 0.5 is medium, and h = 0.8 is a large effect

# with power set at 0.80 & small effect size, n = 392
pwr.2p.test(h = 0.2, power = 0.80, sig.level = 0.05,
            alternative = "two.sided")

#with power set at 0.80 & medium effect size, n = 63
pwr.2p.test(h = 0.5, power = 0.80, sig.level = 0.05,
            alternative = "two.sided")

#with power set at 0.80 & large effect size, n = 25
pwr.2p.test(h = 0.8, power = 0.80, sig.level = 0.05,
            alternative = "two.sided")


## Trying again: calculate Cohen's h with fixed power:

#Known sample sizes
Animal <- 345
Environmental <- 263
Food <- 128
Human <- 500



#Minimum detectable differences: animal vs environmental
pwr.2p2n.test(h = NULL, 
              n1 = Animal, 
              n2 = Environmental, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: animal vs food
pwr.2p2n.test(h = NULL, 
              n1 = Animal, 
              n2 = Food, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: animal vs human
pwr.2p2n.test(h = NULL, 
              n1 = Animal, 
              n2 = Human, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: environmental vs food
pwr.2p2n.test(h = NULL, 
              n1 = Environmental, 
              n2 = Food, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: environmental vs human
pwr.2p2n.test(h = NULL, 
              n1 = Environmental, 
              n2 = Human, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: food vs human
pwr.2p2n.test(h = NULL, 
              n1 = Food, 
              n2 = Human, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ST-specific analysis power calculations ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Known sample sizes
Animal <- 345*0.20
Environmental <- 263*0.20
Food <- 128*0.20
Human <- 500*0.20



#Minimum detectable differences: animal vs environmental
pwr.2p2n.test(h = NULL, 
              n1 = Animal, 
              n2 = Environmental, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: animal vs food
pwr.2p2n.test(h = NULL, 
              n1 = Animal, 
              n2 = Food, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: animal vs human
pwr.2p2n.test(h = NULL, 
              n1 = Animal, 
              n2 = Human, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: environmental vs food
pwr.2p2n.test(h = NULL, 
              n1 = Environmental, 
              n2 = Food, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: environmental vs human
pwr.2p2n.test(h = NULL, 
              n1 = Environmental, 
              n2 = Human, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))

#Minimum detectable differences: food vs human
pwr.2p2n.test(h = NULL, 
              n1 = Food, 
              n2 = Human, 
              sig.level = 0.05, 
              power = 0.80,
              alternative = c("two.sided"))




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chi-square power calculations ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# for chi-square tests:
## effect size is small if w = 0.10, medium if w = 0.30, and large if w = 0.50

# with power set at 0.80 & small effect size, n = 1090
pwr.chisq.test(w = 0.10, df = 3, sig.level = 0.05, power = 0.80)

# with power set at 0.80 & medium effect size, n = 121
pwr.chisq.test(w = 0.30, df = 3, sig.level = 0.05, power = 0.80)

# with power set at 0.80 & large effect size, n = 44
pwr.chisq.test(w = 0.50, df = 3, sig.level = 0.05, power = 0.80)

#chi-square with set sample size
pwr.chisq.test(N = Animal+Environmental+Food+Human, df = 3, sig.level = 0.05, power = 0.80)
