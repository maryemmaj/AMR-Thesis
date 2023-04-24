# Author: Mary Jewell
# Date: 3/14/2023
# Notes: Compile ResFinder results into a dataframe
#       Needs to be run on the cluster where the original ResFinder results live.
# Last updated: 3/14/2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~
# Setup ####
##~~~~~~~~~~~~

rm(list = ls())
setwd("/home/mj49/projects/Ecoli_AMR/data/fastq/resfinder_results")
getwd()

library(tidyr) #for reshaping to wide form
library(dplyr)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading SRR files in a loop ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

total_files <- list.files(pattern = "^SRR") # 1485 files

# create an empty data frame to hold the combined data
combined_df <- data.frame()

# loop over all folders starting with "SRR"
for (dir in list.files(pattern = "^SRR")) {
  
  # get the file path for the current file
  file_path <- file.path(dir, "pheno_table.txt")
  
  # get the SRR number from the directory name
  srr <- gsub("^SRR([0-9]+).*", "\\1", dir)
  
  # read the file into a data frame
  tryCatch({
    # read each individual file
    df <- read.table(file_path, header = F, comment.char = "#", 
                     skip = 5, sep = "\t", na.strings = c("", " "), fill = T,
                     col.names = c("antimicrobial", "class", "WGS_pred_pheno", "match", "gene"))
    # add SRR number
    df$SRR <- paste0("SRR", srr)
    
    # add the data frame to the combined data frame
    combined_df <- rbind(combined_df, df)
    
  }, error = function(e){
    # if there's an error, print a message and move on to the next file
    message(paste("Error reading file", file_path))
  })
  
}



# print the first few rows of the combined data frame
head(combined_df)

length(unique(combined_df$SRR)) # read in 1473 SRR files

# Missing files:
# Error reading file SRR11068088/pheno_table.txt
# Error reading file SRR1198928/pheno_table.txt
# Error reading file SRR1314215/pheno_table.txt
# Error reading file SRR1314232/pheno_table.txt
# Error reading file SRR1314233/pheno_table.txt
# Error reading file SRR1314358/pheno_table.txt
# Error reading file SRR1314560/pheno_table.txt
# Error reading file SRR1314598/pheno_table.txt
# Error reading file SRR13182993/pheno_table.txt
# Error reading file SRR1656097/pheno_table.txt
# Error reading file SRR1795857/pheno_table.txt
# Error reading file SRR2889879/pheno_table.txt
# Warning message:
#   In file(file, "rt") :
#   cannot open file 'SRR2889879/pheno_table.txt': No such file or directory

# Checking data quality
sum(is.na(combined_df$antimicrobial))
sum(is.na(combined_df$class))
sum(is.na(combined_df$WGS_pred_pheno))
sum(is.na(combined_df$match))
sum(is.na(combined_df$SRR))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fixing the skipped file errors ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compile a list of missed files:
files <- c("SRR11068088",
           "SRR1198928",
           "SRR1314215",
           "SRR1314232",
           "SRR1314233",
           "SRR1314358",
           "SRR1314560",
           "SRR1314598",
           "SRR13182993",
           "SRR1656097",
           "SRR1795857",
           "SRR2889879")


# create an empty data frame to hold the combined data
second_pass_df <- data.frame()

# loop over all folders starting with "SRR"
for (dir in files) {
  
  # get the file path for the current file
  file_path <- file.path(dir, "pheno_table.txt")
  
  # get the SRR number from the directory name
  srr <- gsub("^SRR([0-9]+).*", "\\1", dir)
  
  # read the file into a data frame
  tryCatch({
    # read each individual file
    df <- read.table(file_path, header = F, comment.char = "#", 
                     skip = 5, sep = "\t", na.strings = c("", " "), fill = T)
    
    # assign column names, accommodating stray column
    colnames(df) <- c("antimicrobial", "class", "WGS_pred_pheno", "match", "temp", "gene")
    
    # paste relevant information into gene column
    df$gene <- ifelse(!is.na(df$temp), paste(df$temp, df$gene), NA)
    
    # remove temp column
    df$temp <- NULL
    
    # add SRR number
    srr <- gsub("^SRR([0-9]+).*", "\\1", file_path)
    df$SRR <- paste0("SRR", srr)
    
    # add the data frame to the combined data frame
    second_pass_df <- rbind(second_pass_df, df)
  }, error = function(e){
    # if there's an error, print a message and move on to the next file
    message(paste("Error reading file", file_path))
  })
  
}



# Checking data quality
sum(is.na(second_pass_df$antimicrobial))
sum(is.na(second_pass_df$class))
sum(is.na(second_pass_df$WGS_pred_pheno))
sum(is.na(second_pass_df$match))
sum(is.na(second_pass_df$SRR))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Append complete dataframes & write long form ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

resfinder_long <- rbind(combined_df, second_pass_df)
length(unique(resfinder_long$SRR)) # 1484 complete files

write.csv(resfinder_long, "/home/mj49/projects/Ecoli_AMR/data/cleaned_data/resfinder_long.csv")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create wide form & write ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

resfinder_long <- rename(resfinder_long, "pheno_pred" = "WGS_pred_pheno")

resfinder_wide <- resfinder_long %>% 
  pivot_wider(id_cols = SRR, names_from = antimicrobial, values_from = pheno_pred)
names(resfinder_wide)
