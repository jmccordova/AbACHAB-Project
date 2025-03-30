# Part 1: Packages
# Bioconductor Installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = package_loc)
library(BiocManager, lib.loc = package_loc)

BiocManager::install(
  c(
    "withr", "ggplot2", "tidyverse", 
    
    "Amelia", "psych", "dyplr", "mice", "missForest", "mvnmle", "naniar", "misty", "Hmisc",
    
    "lmtest",
    
    "nlme"
  ), 
  force = TRUE, 
  dependencies = TRUE, 
  lib = package_loc
)

library(utf8, lib.loc = package_loc); library(tibble, lib.loc = package_loc); library(readxl, lib.loc = package_loc)
library(Amelia, lib.loc = package_loc); library(psych, lib.loc = package_loc); library(backports, lib.loc = package_loc); library(mice, lib.loc = package_loc); library(missForest, lib.loc = package_loc); library(withr, lib.loc = package_loc); library(ggplot2, lib.loc = package_loc); library(dplyr, lib.loc = package_loc); library(cowplot, lib.loc = package_loc); library(mvnmle, lib.loc = package_loc); library(naniar, lib.loc = package_loc); library(misty, lib.loc = package_loc); library(Hmisc, lib.loc = package_loc)
library(lmtest, lib.loc = package_loc); 
library(nlme, lib.loc = package_loc); 

# Part 2: Reading values
data <- read_excel(paste(datadir, "abachab final results.xlsx", sep = "/"), sheet = "Consolidated")

# Part 3: Separating dataset into their own bins
antibodies <- data[, c("1st Extraction (Day 0)", "2nd Extraction (Day 14)", "3rd Extraction (Day 30)", "4th Extraction (Day 90)", "5th Extraction (Day 180)")]
colnames(antibodies) <- c("Day 0", "Day 14", "Day 30", "Day 90", "Day 180")
antibodies <- mutate_all(antibodies, as.numeric)
perc_inhibition_orig_variant <- data[, c("1st Extraction (Day 0)2", "2nd Extraction (Day 14)3", "3rd Extraction (Day 30)4", "4th Extraction (Day 90)5", "5th Extraction (Day 180)2")]
colnames(perc_inhibition_orig_variant) <- c("Day 0", "Day 14", "Day 30", "Day 90", "Day 180")
perc_inhibition_orig_variant <- mutate_all(perc_inhibition_orig_variant, as.numeric)
perc_inhibition_delta_variant <- data[, c("1st Extraction DELTA", "2nd Extraction DELTA", "3rd Extraction DELTA", "4th Extraction DELTA", "5th Extraction DELTA")]
colnames(perc_inhibition_delta_variant) <- c("Day 0", "Day 14", "Day 30", "Day 90", "Day 180")
perc_inhibition_delta_variant <- mutate_all(perc_inhibition_delta_variant, as.numeric)

# Part 4: Describe datasets
summary(antibodies)
describe(antibodies)
summary(perc_inhibition_orig_variant)
describe(perc_inhibition_orig_variant)
summary(perc_inhibition_delta_variant)
describe(perc_inhibition_delta_variant)

# Part 5: Visualization
# Visualize outliers
draw_boxplot <- function(df0, df1, df2) {
  par(mfrow = c(3, 1))  # 3 rows, 1 column
  boxplot(df0)
  boxplot(df1)
  boxplot(df2)
  # Reset plotting area
  par(mfrow = c(1, 1))
}

# Visualize the distribution of each column
draw_histograms <- function(df, title, label) {
  par(mfrow = c(2, 3))  # 2 rows, 3 columns
  # Create histograms
  if (!grepl("Delta", title)) {
    hist(df$`Day 0`, main = paste(title,"(Day 0)"), xlab = label, col = "red") 
  }
  hist(df$`Day 14`, main = paste(title,"(Day 14)"), xlab = label, col = "orange")
  hist(df$`Day 30`, main = paste(title,"(Day 30)"), xlab = label, col = "yellow")
  hist(df$`Day 90`, main = paste(title,"(Day 90)"), xlab = label, col = "green")
  hist(df$`Day 180`, main = paste(title,"(Day 180)"), xlab = label, col = "blue")
  # Reset plotting area
  par(mfrow = c(1, 1))
}

#IgG Antibodies
draw_histograms(antibodies, "IgG antibodies", "IgG")
# %Inhibition original variant
draw_histograms(perc_inhibition_orig_variant, "%Inhibition Original Variant", "% inhibition")
# %Inhibition delta variant
draw_histograms(perc_inhibition_delta_variant, "%Inhibition Delta Variant", "% inhibition")

# Part 6: Get vaccine values
data.vaccine <- read_excel(paste(datadir, "AbACHAB MASTERLIST.xlsx", sep = "/"), sheet = "Complete (Cleaned)")
