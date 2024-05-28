# WinterEcology I experiment
# Script author: Emil A.S. Andersen
#
#
# Load Libraries
library(tidyverse) # For the loading data in tibble and organizing it
library(readxl) # For reading excel files with the "read_excel" command
library(writexl) # For writing excel files
#
#
# Soil extraction (SE and SEF) mg nitrogen per g DW
# [Only SE right now]
SE_N_DW <- read.csv2("Raw_data/Soil_extraction_mg_N_g_DW.csv",header=TRUE)
SE_N_DW <- SE_N_DW %>%
  mutate(across(mg_N_pr_g_DW, as.numeric))# %>%
#  mutate(across(MP, as.character))
#
# Order Month-yr in order as given in file
SE_N_DW$Month.yr <- factor(SE_N_DW$Month.yr, levels = unique(SE_N_DW$Month.yr))
#
# Plot SE based on Month-yr and site
box_mgN_gDW <- ggplot(SE_N_DW, aes(Month.yr, mg_N_pr_g_DW, colour = Site))
box_mgN_gDW + geom_boxplot()
#
