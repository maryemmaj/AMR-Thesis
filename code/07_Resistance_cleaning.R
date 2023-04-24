# Author: Mary Jewell
# Date: 3/14/2023
# Notes: Clean ResFinder results and convert to wide form for analysis
# Last updated: 3/14/2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~~

rm(list = ls())
setwd("C:/Users/marye/OneDrive/Desktop/Thesis")

library(tidyr) #for reshaping to wide form
library(dplyr) #data manipulation: rename
library(ggplot2) #plotting

resfinder_long <- read.csv("Data/resfinder_long.csv")
resfinder_long$X <- NULL

meta_mlst <- read.csv("Data/meta_mlst.csv")
meta <- read.csv("Data/meta_short.csv")


##~~~~~~~~~~~~~~~~~~~~~~~~~~
# Long form data ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~

resfinder_long <- resfinder_long %>% rename("pheno_pred" = "WGS_pred_pheno")
names(resfinder_long)

# There are several values that belong in "gene" but are in antimicrobial or class instead
table(resfinder_long$antimicrobial)
misplaced_genes <- c("aac(3)-IIa (aac(3)-IIa_L22613), aac(3)-IIa (aac(3)-IIa_CP023555)",
                     "blaCMY-2 (blaCMY-2_X91840), blaTEM-1B (blaTEM-1B_AY458016)",
                     "blaTEM-1B (blaTEM-1B_AY458016), blaCTX-M-15 (blaCTX-M-15_AY044436)",
                     "blaTEM-1D (blaTEM-1D_AF188200), blaNDM-5 (blaNDM-5_JN104597)",
                     "dfrA14 (dfrA14_DQ388123), dfrA12 (dfrA12_AM040708), dfrA14 (dfrA14_AF393510)",
                     "dfrA14 (dfrA14_DQ388123), dfrA14 (dfrA14_AF393510)",
                     "dfrA14 (dfrA14_DQ388123), dfrA14 (dfrA14_KF921535)",
                     "erm(B) (erm(B)_JN899585), mph(A) (mph(A)_D16251)",
                     "erm(B) (erm(B)_M19270), mph(A) (mph(A)_D16251)",
                     "parC;;1;;CP084529.1_108_t_AA",
                     "qnrB10 (qnrB10_DQ631414), qnrB19 (qnrB19_EU432277), qnrB82 (qnrB82_KX372672)",
                     "qnrB5 (qnrB5_DQ303919), qnrB10 (qnrB10_HM439644), qnrB19 (qnrB19_EU432277), qnrB36 (qnrB36_JN173058), qnrB67 (qnrB67_KC580656), qnrB70 (qnrB70_KC580659), qnrB81 (qnrB81_KX372671), qnrB82 (qnrB82_KX372672)",
                     "sul1 (sul1_EU780013), sul2 (sul2_AY034138)",
                     "sul1 (sul1_U12338), sul2 (sul2_AY034138)",
                     "tet(A) (tet(A)_AF534183), tet(B) (tet(B)_AF326777)",
                     "tet(A) (tet(A)_AF534183), tet(D) (tet(D)_AF467077)",
                     "tet(B) (tet(B)_AP000342), tet(B) (tet(B)_AF326777)")

temp1 <- resfinder_long %>% 
  filter(antimicrobial %in% misplaced_genes)

temp2 <- resfinder_long %>%
  filter(class == "parC;;1;;CP084529.1")

# Append misplaced gene value from antimicrobial to gene
resfinder_long$gene <- ifelse(
  resfinder_long$antimicrobial %in% misplaced_genes, # if gene value is in antimicrobial column
  resfinder_long$antimicrobial, # assign value to gene
  resfinder_long$gene) # else keep existing gene value

# Replace the gene value with NA in antimicrobial column
resfinder_long$antimicrobial <- ifelse(
  resfinder_long$antimicrobial %in% misplaced_genes, # if gene is in antimicrobial
  NA, # replace with NA
  resfinder_long$antimicrobial) # else keep existing antimicrobial value


## Do the same with the misplaced genes in class column:

# Append misplaced gene value from class to gene
resfinder_long$gene <- ifelse(
  resfinder_long$class == "parC;;1;;CP084529.1", # if gene value is in class column
  paste(resfinder_long$gene, resfinder_long$class), # assign to gene column (paste existing value)
  resfinder_long$gene) # else keep class value

resfinder_long$class <- ifelse(
  resfinder_long$class == "parC;;1;;CP084529.1", # if gene value is in class column
  NA, # replace with NA
  resfinder_long$class) # else keep existing class value



# Check results:
temp1 <- resfinder_long %>% 
  filter(antimicrobial %in% misplaced_genes)

temp2 <- resfinder_long %>%
  filter(class == "parC;;1;;CP084529.1")

table(resfinder_long$antimicrobial)
table(resfinder_long$class)

## Merge metadata with long data ####

long_merged <- merge(resfinder_long, meta, by = "SRR", all = F)
names(long_merged)

write.csv(long_merged, "Data/resfinder_long_cleaned.csv", row.names = F)


##~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wide form data ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~

resfinder_wide <- resfinder_long %>% 
  pivot_wider(id_cols = SRR, names_from = antimicrobial, values_from = pheno_pred)

names(resfinder_wide)
resfinder_wide$`NA` <- NULL #removing rows where antimicrobial was == NA (??)


resfinder_wide$res_sum <- rowSums(resfinder_wide == "Resistant")
summary(resfinder_wide$res_sum)

# Binary resistance: less than 3 or 3+
resfinder_wide$res_binary <- ifelse(resfinder_wide$res_sum >= 3, 1, 0)
table(resfinder_wide$res_binary)
prop.table(table(resfinder_wide$res_binary))

# Categorical resistance: 0-2, 3-9, or 10+
resfinder_wide$res_cat[resfinder_wide$res_sum < 3] <- 1 #susceptible
resfinder_wide$res_cat[resfinder_wide$res_sum >= 3 &
                       resfinder_wide$res_sum < 10] <- 2 #resistant: 3-9 antimicrobials
resfinder_wide$res_cat[resfinder_wide$res_sum >= 10] <- 3 #highly resistant: 10+ antimicrobials
table(resfinder_wide$res_cat)

# Any resistance
resfinder_wide$res_any <- ifelse(resfinder_wide$res_sum >= 1, 1, 0)
table(resfinder_wide$res_any)
prop.table(table(resfinder_wide$res_any))



# Convert everything to character in order to write csv
resfinder_wide <- resfinder_wide %>%
  mutate(across(everything(), as.character))

write.csv(resfinder_wide, "Data/resfinder_wide.csv", row.names = F)



## Merge meta_mlst and resfinder data ####
meta_res <- merge(meta_mlst, resfinder_wide, by = "SRR", all.y = T)
names(meta_res)

write.csv(meta_res, "Data/full_analytic_data.csv", row.names = F)



