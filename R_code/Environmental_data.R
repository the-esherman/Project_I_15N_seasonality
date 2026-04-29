# WinterEcology I experiment
# Script author: Emil A.S. Andersen
#
# Load weather data
# From own data
#
# To get months in the right format (i.e. not in whatever local the computer has, e.g. Swedish)
Sys.setlocale("LC_ALL", 'en_GB.UTF-8')
#
library(tidyverse)
library(readr)
library(lubridate)
library(readxl)
library(viridis)
library(gridExtra)
library(cowplot)
#
#
#------- ### Load the data ### -------
#------- # Air temperature logger # -------|
Tair_cols <- c("date", "date", "text","date","text","text","text","date","text","text","text","date","text","text","text","text", "text")
Abisko_Tair <- read_xlsx("raw_data/Abisko_Tair.xlsx", skip = 5, col_names = TRUE, col_types = Tair_cols)
Vassijaure_Tair <- read_xlsx("raw_data/Vassijaure_Tair.xlsx", skip = 5, col_names = TRUE, col_types = Tair_cols)
#
# Soil temperature and moisture data and PAR sensor
EM50_cols <- c("date", "date", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "date", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "date", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "date", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "date", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
Abisko_EM50 <- read_xlsx("raw_data/allEMdata Abisko.xlsx", col_names = TRUE, col_types = EM50_cols)
Vassijaure_EM50 <- read_xlsx("raw_data/allEMdata Vassijaure.xlsx", col_names = TRUE, col_types = EM50_cols)
#
#
# SMHI weather station data
SMHI_Abisko_Tair <- read_csv2("raw_data/smhi_Temp_Abisko_new.csv", skip = 10, col_select = c(1:4))
SMHI_Katterjakk_Tair <- read_csv2("raw_data/smhi_Kattejakk_Temp.csv", skip = 9, col_select = c(1:4))
#
# Define the winter period as snow covered period
# First and last day of measuring snow on entire plots
# Abisko: 20191021 - 20200506; Vassijaure: 20191028-20200601
# Safe to assume Vassijaure also started earlier than when measured, thus assumed same day.
winterP <- data.frame(wstart = c(ymd(20191021), ymd(20191021)), wend = c(ymd(20200506), ymd(20200601)))
winterP_date <- data.frame(wstart = c(as.Date("2019-11-10"),as.Date("2019-11-12")), wend = c(as.Date("2020-05-06"),as.Date("2020-06-01")))
#
# Snow depth
snowData <- read_xlsx("raw_data/Snow_depth.xlsx", col_names = TRUE)
#
# Core data on soils and days of labelling and harvest as well as soil moisture (FW and DW) and snowdepth
coreData <- read_csv("clean_data/Core_data.csv", col_names = TRUE)
#
#
#
#=======  ###   Functions    ### =======
# From http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     max  = max    (xx[[col]], na.rm=na.rm),
                     min  = min    (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult 
  return(datac)
}
#
#
#
#------- ### Clean data ### -------
#------- # Air Temperature # -------
# Select specific rows
Abisko_Tair <- Abisko_Tair %>%
  dplyr::select(1:4, "Tair_C1", 8, "Tair_31", 12, "Tair_39_2", "Tair_A", "Tair_A2") %>%
  rename("Time_A39_1" = "...2",
         "Time_C1" = "...4",
         "Time_31" = "...8",
         "Time_39_2" = "...12") %>%
  mutate(across(1:10, as.character)) %>%
  mutate(across(where(is.character), ~na_if(.,"NaN"))) %>%
  mutate(across(c("Tair_A39_1", "Tair_C1", "Tair_31", "Tair_39_2", "Tair_A", "Tair_A2"), as.numeric)) %>%
  mutate(across(c(Date_A, Time_A39_1, Time_C1, Time_31, Time_39_2), ymd_hms))
#
Abisko_Tair_long <- Abisko_Tair %>%
  dplyr::select(1,3,5,7,9,11) %>%
  pivot_longer(cols = 2:6, names_to = "Sensor", values_to = "Air_temp")
#
Vassijaure_Tair <- Vassijaure_Tair %>%
  dplyr::select(1:4, "Tair_C7", 8, "Tair_38", 12, "Tair_C1", "Tair_V", "Tair_V2") %>%
  rename("Time_B3" = "...2",
         "Time_C7" = "...4",
         "Time_38" = "...8",
         "Time_C1" = "...12") %>%
  mutate(across(1:10, as.character)) %>%
  mutate(across(where(is.character), ~na_if(.,"NaN"))) %>%
  mutate(across(c("Tair_B3", "Tair_C7", "Tair_38", "Tair_C1", "Tair_V", "Tair_V2"), as.numeric)) %>%
  mutate(across(c(Date_V, Time_B3, Time_C7, Time_38, Time_C1), ymd_hms))
#
Vassijaure_Tair_long <- Vassijaure_Tair %>%
  dplyr::select(1,3,5,7,9,11) %>%
  pivot_longer(cols = 2:6, names_to = "Sensor", values_to = "Air_temp")
#
# Mean diel temperature
Abisko_avgTair <- Abisko_Tair %>%
  group_by(date(Date_A)) %>%
  summarise(Abisko = mean(Tair_A2, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Date_A)")
#
Vassijaure_avgTair <- Vassijaure_Tair %>%
  group_by(date(Date_V)) %>%
  summarise(Vassijaure = mean(Tair_V2, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Date_V)")
#
avgTair_wide <- left_join(Abisko_avgTair, Vassijaure_avgTair) %>%
  rename("Abisko_Tair" = Abisko,
         "Vassijaure_Tair" = Vassijaure)
# Combine Abisko and Vassijaure and pivot longer
avgTair_long <- left_join(Abisko_avgTair, Vassijaure_avgTair) %>%
  pivot_longer(2:3, names_to = "Site", values_to = "dielT_air")
#
#
#
#------- # EM50 loggers # -------
# Remove last lines without a date and extra header part. Set NaN to NA
Abisko_EM50 <- Abisko_EM50 %>%
  slice(5:n()) %>%
  dplyr::select(-c(NA...9, NA...18, NA...27, NA...36, NA...45)) %>%
  dplyr::filter(!(is.na(Date))) %>%
  mutate(across(where(is.character), ~na_if(.,"NaN")),
         across(where(is.numeric), ~na_if(.,NaN))) %>%
  mutate(across(c(Date, A1_Date, A2_Date, A3_Date, A4_Date, A5_Date), ymd_hms))
#
Vassijaure_EM50 <- Vassijaure_EM50 %>%
  slice(5:(n()-24)) %>% # No data for the last 24 logs
  dplyr::select(-c(NA...18, NA...27, NA...36, NA...37, NA...43, NA...44)) %>%
  filter(!(is.na(V_Date))) %>% # if the date is not there, we cannot use the data
  mutate(across(where(is.character), ~na_if(.,"NaN")),
         across(where(is.numeric), ~na_if(.,NaN))) %>%
  mutate(across(c(V_Date, V1_Date, V2_Date, V3_Date, V4_Date, V5_Date), ymd_hms))
#
#
#
#------- # Soil Temperature # -------
# Mean Soil temperature
Abisko_avgTsoil <- Abisko_EM50 %>%
  select(Date, A1N_Tsoil, A2N_Tsoil, A3N_Tsoil, A4N_Tsoil, A5N_Tsoil) %>%
  # Remove soil temperatures below -15 C
  mutate(A1N_Tsoil = replace(A1N_Tsoil, A1N_Tsoil < -15, NA), # Has very low values in July 2019 and during winter
         A2N_Tsoil = replace(A2N_Tsoil, A2N_Tsoil < -15, NA),
         A3N_Tsoil = replace(A3N_Tsoil, A3N_Tsoil < -15, NA), # Has a single measure very low measure in March 2020
         A4N_Tsoil = replace(A4N_Tsoil, A4N_Tsoil < -15, NA),
         A5N_Tsoil = replace(A5N_Tsoil, A5N_Tsoil < -15, NA)) %>%
  mutate(A1N_Tsoil = if_else(month(Date) == 7 & A1N_Tsoil < 5, NA, A1N_Tsoil),
         A2N_Tsoil = if_else(month(Date) == 7 & A2N_Tsoil < 0, NA, A2N_Tsoil),
         A3N_Tsoil = if_else(month(Date) == 7 & A3N_Tsoil < 0, NA, A3N_Tsoil),
         A4N_Tsoil = if_else(month(Date) == 7 & A4N_Tsoil < 0, NA, A4N_Tsoil),
         A5N_Tsoil = if_else(month(Date) == 7 & A5N_Tsoil < 0, NA, A5N_Tsoil)) %>%
  group_by(date(Date)) %>%
  summarise(A1N_Tsoil = mean(A1N_Tsoil, na.rm = TRUE),
            A2N_Tsoil = mean(A2N_Tsoil, na.rm = TRUE),
            A3N_Tsoil = mean(A3N_Tsoil, na.rm = TRUE),
            A4N_Tsoil = mean(A4N_Tsoil, na.rm = TRUE),
            A5N_Tsoil = mean(A5N_Tsoil, na.rm = TRUE),
            .groups = "keep") %>%
  rename("Date" = "date(Date)") %>%
  mutate(Abisko = mean(c(A1N_Tsoil, A2N_Tsoil, A3N_Tsoil, A4N_Tsoil, A5N_Tsoil), na.rm = TRUE))
#
Vassijaure_avgTsoil <- Vassijaure_EM50 %>%
  select(V_Date, V1N_Tsoil, V2N_Tsoil, V3N_Tsoil, V4N_Tsoil, V5N_Tsoil) %>%
  # Remove soil temperatures below -15 C
  mutate(V1N_Tsoil = replace(V1N_Tsoil, V1N_Tsoil < -15, NA), # Has very low values in July 2019 and during winter
         V2N_Tsoil = replace(V2N_Tsoil, V2N_Tsoil < -15, NA),
         V3N_Tsoil = replace(V3N_Tsoil, V3N_Tsoil < -15, NA), # Has a single measure very low measure in March 2020
         V4N_Tsoil = replace(V4N_Tsoil, V4N_Tsoil < -15, NA),
         V5N_Tsoil = replace(V5N_Tsoil, V5N_Tsoil < -15, NA)) %>%
  group_by(date(V_Date)) %>%
  summarise(V1N_Tsoil = mean(V1N_Tsoil, na.rm = TRUE),
            V2N_Tsoil = mean(V2N_Tsoil, na.rm = TRUE),
            V3N_Tsoil = mean(V3N_Tsoil, na.rm = TRUE),
            V4N_Tsoil = mean(V4N_Tsoil, na.rm = TRUE),
            V5N_Tsoil = mean(V5N_Tsoil, na.rm = TRUE),
            .groups = "keep") %>%
  rename("Date" = "date(V_Date)") %>%
  mutate(Vassijaure = mean(c(V1N_Tsoil, V2N_Tsoil, V3N_Tsoil, V4N_Tsoil, V5N_Tsoil), na.rm = TRUE))
#
# Average across plots
avgTsoil_wide <- left_join(Abisko_avgTsoil, Vassijaure_avgTsoil) %>%
  select(c(Date, Abisko, Vassijaure)) %>%
  rename("Abisko_Tsoil" = Abisko,
         "Vassijaure_Tsoil" = Vassijaure)
# Combine Abisko and Vassijaure and pivot longer
avgTsoil_long <- left_join(Abisko_avgTsoil, Vassijaure_avgTsoil) %>%
  dplyr::select(c(Date, Abisko, Vassijaure)) %>%
  pivot_longer(2:3, names_to = "Site", values_to = "dielT_soil")
#
#
# Combine all temperatures
avgT_wide <- left_join(avgTair_wide, avgTsoil_wide) %>% ungroup()
avgT_long <- left_join(avgTair_long, avgTsoil_long) %>% ungroup()
#
#
#
#------- # Soil Moisture # -------
#
#
# Mean Soil moisture
Abisko_avgVWC <- Abisko_EM50 %>%
  select(Date, A1N_VWC, A2N_VWC, A3N_VWC, A4N_VWC, A5N_VWC) %>%
  group_by(date(Date)) %>%
  summarise(A1N_VWC = mean(A1N_VWC, na.rm = TRUE),
            A2N_VWC = mean(A2N_VWC, na.rm = TRUE),
            A3N_VWC = mean(A3N_VWC, na.rm = TRUE),
            A4N_VWC = mean(A4N_VWC, na.rm = TRUE),
            A5N_VWC = mean(A5N_VWC, na.rm = TRUE),
            .groups = "keep") %>%
  rename("Date" = "date(Date)") %>%
  mutate(Abisko = mean(c(A1N_VWC, A2N_VWC, A3N_VWC, A4N_VWC, A5N_VWC), na.rm = TRUE)) %>%
  ungroup()
#
Vassijaure_avgVWC <- Vassijaure_EM50 %>%
  select(V_Date, V1N_VWC, V2N_VWC, V3N_VWC, V4N_VWC, V5N_VWC) %>%
  group_by(date(V_Date)) %>%
  summarise(V1N_VWC = mean(V1N_VWC, na.rm = TRUE),
            V2N_VWC = mean(V2N_VWC, na.rm = TRUE),
            V3N_VWC = mean(V3N_VWC, na.rm = TRUE),
            V4N_VWC = mean(V4N_VWC, na.rm = TRUE),
            V5N_VWC = mean(V5N_VWC, na.rm = TRUE),
            .groups = "keep") %>%
  rename("Date" = "date(V_Date)") %>%
  mutate(Vassijaure = mean(c(V1N_VWC, V2N_VWC, V3N_VWC, V4N_VWC, V5N_VWC), na.rm = TRUE)) %>% # Average across plots
  ungroup()
#
# Format wide: Abisko and Vassijaure separate columns
avgVWC_wide <- left_join(Abisko_avgVWC, Vassijaure_avgVWC, by = join_by(Date)) %>%
  select(c(Date, Abisko, Vassijaure)) %>%
  rename("Abisko_VWC" = Abisko,
         "Vassijaure_VWC" = Vassijaure) %>%
  mutate(Abisko_VWC = Abisko_VWC*100,
         Vassijaure_VWC = Vassijaure_VWC*100) # % by Volume
#
# Format long: Abisko and Vassijaure as Site factor
avgVWC_long <- left_join(Abisko_avgVWC, Vassijaure_avgVWC, by = join_by(Date)) %>%
  dplyr::select(c(Date, Abisko, Vassijaure)) %>%
  pivot_longer(2:3, names_to = "Site", values_to = "dielVWC_soil") %>%
  mutate(dielVWC_soil = dielVWC_soil*100) # % by Volume
#
# Remove dates during the snow-covered period
# For long dataset
avgVWC_long2 <- avgVWC_long %>%
  mutate(dielVWC_soil = case_when(Site == "Abisko" & ymd(Date) >= winterP$wstart[1] & ymd(Date) <= winterP$wend[1] ~ NA,
                                  Site == "Vassijaure" & ymd(Date) >= winterP$wstart[2] & ymd(Date) <= winterP$wend[2] ~NA,
                                  TRUE ~ dielVWC_soil))
#
# For wide dataset
avgVWC_wide2 <- avgVWC_wide %>%
  mutate(Abisko_VWC = case_when(ymd(Date) >= winterP$wstart[1] & ymd(Date) <= winterP$wend[1] ~ NA,
                                TRUE ~ Abisko_VWC),
         Vassijaure_VWC = case_when(ymd(Date) >= winterP$wstart[2] & ymd(Date) <= winterP$wend[2] ~ NA,
                                    TRUE ~ Vassijaure_VWC))
#
#
# Gravimetric water content (GWC) and estimated bulk density (BD)
# Taken from coreData
# Convert gravimetric to volumetric by estimating bulk density
coreGWC <- coreData %>%
  mutate(Soil_GWC = (SM_FW_g - SM_DW_g)/SM_DW_g*100,
         Soil_GWC_FW = (SM_FW_g - SM_DW_g)/SM_FW_g*100,
         Soil_RF_core_g_FW = Soil_RF_FW_g*Soil_mass_g,) %>%
  mutate(Soil_BD_est = Soil_RF_core_g_FW*(DW_FW_frac/Soil_vol_cm3),
         Soil_BD_est1 = case_when(Site == "Abisko" ~ 0.110854207987199, # Bulk density based on a few samples in the autumn 2019
                                  Site == "Vassijaure" ~ 0.101975726059452,
                                  TRUE ~ 0)) %>%
  mutate(Soil_VWC_est = Soil_GWC*Soil_BD_est,
         Soil_VWC_est1 = Soil_GWC*Soil_BD_est1) %>%
  select(1:3, Day_of_harvest, Soil_GWC, Soil_GWC_FW, Soil_VWC_est, Soil_VWC_est1) %>%
  rename("Date" = Day_of_harvest)
#
# Summarise with function
# Estimate of BD from soil mass and water content
coreGWC_1.1 <- summarySE(coreGWC, measurevar = "Soil_VWC_est", groupvars = c("Site", "MP", "Date")) %>%
  rename("N_VWC" = N,
         "sd_VWC" = sd,
         "max_VWC" = max,
         "min_VWC" = min,
         "se_VWC" = se,
         "ci_VWC" = ci)
# Estimate of BD from a few samples autumn 2019
coreGWC_1.2 <- summarySE(coreGWC, measurevar = "Soil_VWC_est1", groupvars = c("Site", "MP", "Date")) %>%
  rename("N1" = N,
         "sd1" = sd,
         "max1" = max,
         "min1" = min,
         "se1" = se,
         "ci1" = ci)
#
coreGWC_1.3 <- summarySE(coreGWC, measurevar = "Soil_GWC_FW", groupvars = c("Site", "MP", "Date")) %>%
  rename("N_GWC" = N,
         "sd_GWC" = sd,
         "max_GWC" = max,
         "min_GWC" = min,
         "se_GWC" = se,
         "ci_GWC" = ci)
#
# Combine both VWC estimates and GWC
coreGWC_1 <- left_join(coreGWC_1.1, coreGWC_1.2, by = join_by(Site, MP, Date)) %>%
  left_join(coreGWC_1.3, by = join_by(Site, MP, Date)) %>%
  select(Site, Date, 5:10, 12:17, 19:24)
#
#
# Combine VWC(sensor) and VWC(GWC)
# Wide data
avgVWC_wide3 <- avgVWC_wide2 %>%
  left_join(coreGWC_1, by = join_by(Date))
#
# Long data
avgVWC_long3 <- avgVWC_long2 %>%
  left_join(coreGWC_1, by = join_by(Date, Site))
#
#
#
#------- # PAR # -------
# Mean Soil PAR
# Abisko
Abisko_avgPAR <- Abisko_EM50 %>%
  select(Date, A1_PAR, A2_PAR, A3_PAR, A4_PAR, A5_PAR) %>%
  # Remove PAR above 3000 µmol pr m2 pr s1
  mutate(A1_PAR = replace(A1_PAR, A1_PAR > 3000, NA),
         A2_PAR = replace(A2_PAR, A2_PAR > 3000, NA),
         A3_PAR = replace(A3_PAR, A3_PAR > 3000, NA),
         A4_PAR = replace(A4_PAR, A4_PAR > 3000, NA),
         A5_PAR = replace(A5_PAR, A5_PAR > 3000, NA)) %>%
  group_by(date(Date)) %>%
  summarise(A1_PAR = mean(A1_PAR, na.rm = TRUE),
            A2_PAR = mean(A2_PAR, na.rm = TRUE),
            A3_PAR = mean(A3_PAR, na.rm = TRUE),
            A4_PAR = mean(A4_PAR, na.rm = TRUE),
            A5_PAR = mean(A5_PAR, na.rm = TRUE),
            .groups = "keep") %>%
  rename("Date" = "date(Date)") %>%
  mutate(Abisko = mean(c(A1_PAR, A2_PAR, A3_PAR, A4_PAR, A5_PAR), na.rm = TRUE)) %>%
  ungroup()
#
# Vassijaure
Vassijaure_avgPAR <- Vassijaure_EM50 %>%
  select(V_Date, V1_PAR, V2_PAR, V3_PAR, V4_PAR, V5_PAR) %>%
  group_by(date(V_Date)) %>%
  summarise(V1_PAR = mean(V1_PAR, na.rm = TRUE),
            V2_PAR = mean(V2_PAR, na.rm = TRUE),
            V3_PAR = mean(V3_PAR, na.rm = TRUE),
            V4_PAR = mean(V4_PAR, na.rm = TRUE),
            V5_PAR = mean(V5_PAR, na.rm = TRUE),
            .groups = "keep") %>%
  rename("Date" = "date(V_Date)") %>%
  mutate(Vassijaure = mean(c(V1_PAR, V2_PAR, V3_PAR, V4_PAR, V5_PAR), na.rm = TRUE)) %>%
  ungroup()
#
# Average across plots
avgPAR_wide <- left_join(Abisko_avgPAR, Vassijaure_avgPAR, by = join_by(Date)) %>%
  select(c(Date, Abisko, Vassijaure)) %>%
  rename("Abisko_PAR" = Abisko,
         "Vassijaure_PAR" = Vassijaure)
# Combine Abisko and Vassijaure and pivot longer
avgPAR_long <- left_join(Abisko_avgPAR, Vassijaure_avgPAR, by = join_by(Date)) %>%
  dplyr::select(c(Date, Abisko, Vassijaure)) %>%
  pivot_longer(2:3, names_to = "Site", values_to = "diel_PAR")
#
#
#
#------- # SMHI weather data # -------
# Abisko
SMHI_Abisko_Tair <- SMHI_Abisko_Tair %>%
  dplyr::select(1:3) %>%
  mutate(across(Datum, ymd), # If also for time: across("Tid (UTC)", hms))
         across(Lufttemperatur, as.numeric)) %>% 
  dplyr::filter(year(Datum) == 2019 | year(Datum) == 2020) %>%
  dplyr::filter(!(year(Datum) == 2020 & month(Datum) > 10),
         !(year(Datum) == 2019 & month(Datum) < 6))
#
SMHI_A_avgTair <- SMHI_Abisko_Tair %>%
  group_by(date(Datum)) %>%
  summarise(Abisko_Tair_SMHI = mean(Lufttemperatur, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Datum)") %>%
  ungroup()
#
# Katterjåkk (closest station to Vassijaure)
SMHI_Katterjakk_Tair <- SMHI_Katterjakk_Tair %>%
  dplyr::select(1:3) %>%
  mutate(across(Datum, ymd), # If also for time: across("Tid (UTC)", hms))
         across(Lufttemperatur, as.numeric)) %>% 
  dplyr::filter(year(Datum) == 2019 | year(Datum) == 2020) %>%
  dplyr::filter(!(year(Datum) == 2020 & month(Datum) > 10),
         !(year(Datum) == 2019 & month(Datum) < 6))
#
SMHI_K_avgTair <- SMHI_Katterjakk_Tair %>%
  group_by(date(Datum)) %>%
  summarise(Katterjakk_Tair_SMHI = mean(Lufttemperatur, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Datum)") %>%
  ungroup()
#
# Combine

# Combine all
avgT_wide2 <- left_join(SMHI_A_avgTair, avgT_wide) %>%
  left_join(SMHI_K_avgTair) %>%
  relocate(Abisko_Tair_SMHI, .before = Katterjakk_Tair_SMHI)
#
# Pivot long
avgT_long2 <- avgT_wide2 %>%
  pivot_longer(2:7, names_to = "Site", values_to = "dielT")
#
#
#
#------- ### Plot ### -------
#
# Snow depth on plot
snowData_2 <- summarySE(snowData, measurevar = "Snow_depth_cm", groupvars = c("Site", "MP", "Date"))
snowData_2 <- snowData_2 %>%
  mutate(Date = ymd(Date))
#
#
# Extract date for sample. Here using Day of harvest
DayOf <- coreData %>%
  select(Site, Round, Day_of_harvest) %>%
  distinct(Day_of_harvest, .keep_all = TRUE) %>%
  mutate(across(Day_of_harvest, ~ as.character(.x)))
DayOfL <- coreData %>%
  select(Site, Round, Day_of_label) %>%
  distinct(Day_of_label, .keep_all = TRUE) %>%
  mutate(across(Day_of_label, ~ as.character(.x)))
DayOfLH <- DayOf %>%
  left_join(DayOfL, by = join_by(Site, Round)) %>%
  filter(Site == "Abisko")
#
#
# <><><><><> MAIN PLOT - FIG 1 <><><><><>
#
#
# Start and end date shown
plotDayStart <- as.Date("2019-07-12")
plotDayEnd <- as.Date("2020-08-26")
#
#
# Air temperatures - all
airT_plot <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) + # Katterjåkk
  #geom_line(aes(x = Date, y = Vassijaure_Tair, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) + 
  geom_line(aes(x = Date, y = Abisko_Tair_SMHI, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  #geom_line(aes(x = Date, y = Abisko_Tair, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(plotDayStart,plotDayEnd)) +
  labs(x = NULL, y = "Air (°C)") + # x = "Time of year",  , title = "Air temperature" 
  #guides(lty = guide_legend(title = "Mean diel temperature"))+ #lty = guide_legend(title = "Mean diel temperature")) +
  theme_bw(base_size = 25) +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.text.y = element_text(size = 15))
#
airT_legend <- get_plot_component(airT_plot, "guide-box", return_all = TRUE)#[[4]]  # 1 is right, 2 is left, 3 is bottom, 4 is top
airT_legend
airT_plot.2 <- airT_plot + theme_bw(base_size = 25) + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 18))#, axis.title.y = element_text(size = 25)) 
#airT_plot <- airT_plot + guides(lty = NULL)
#
# Soil temperatures - all
soilT_plot <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  #annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -3, ymax = 0, fill = "white", alpha = 0.9) +
  #geom_hline(yintercept = -3, color = "#D55E00") +
  geom_line(aes(x = Date, y = Vassijaure_Tsoil, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) +
  geom_line(aes(x = Date, y = Abisko_Tsoil, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15))+#, minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%d-%b") +
  coord_cartesian(xlim = c(plotDayStart,plotDayEnd), ylim = c(-5,15)) +
  labs(x = NULL, y = "Soil (°C)") + # x = "Time of year", , title = "Soil temperature"
  guides(lty = "none") + # guide_legend(title = "Soil temperature")
  theme_bw(base_size = 25) +
  theme(legend.position = "top")#, axis.title.y = element_text(size = 25))
#
soilT_legend <- get_x_axis(soilT_plot) + theme_bw(base_size = 17)
soilT_plot.2 <- soilT_plot + theme_bw(base_size = 25) + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 18)) 
# #soilT_plot <- soilT_plot + guides(lty = NULL)
#
#
snowDepth_plot <- snowData_2 %>%
  ggplot(aes(x = Date, y = Snow_depth_cm, ymin = min, ymax = max, fill = Site, linetype = Site)) +
  geom_ribbon(alpha = 0.5) +
  geom_line(linewidth = 1) +
  scale_fill_grey() +
  #scale_fill_viridis_d() +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%b-%d") +
  coord_cartesian(xlim = c(plotDayStart,plotDayEnd)) +
  labs(x = "Time of year", y = "Snow cover (cm)") + #, title = "Snow cover measured over the entire plot (around all 15 patches)") +
  guides(fill = guide_legend(title = "Snow"), linetype = "none") +
  theme_bw(base_size = 25) +
  theme(legend.position = "bottom")
#
snowData_legend <- get_plot_component(snowDepth_plot, "guide-box", return_all = TRUE)#[[3]]  # 1 is right, 2 is left, 3 is bottom, 4 is top
snowDepth_plot.2 <- snowDepth_plot + theme_bw(base_size = 25) + theme(legend.position = "none")#, axis.title.y = element_text(size = 25)) 
#
# Plot
#grid.arrange(airT_plot, soilT_plot, snowDepth_plot)
#grid.arrange(airT_plot, soilT_plot, snowDepth_plot, top = grid::textGrob('Mean diel temperature', gp=grid::gpar(fontsize=20)))
#
# Fig. 1: Air temperature, Soil temperature, and Snow depth
grid.arrange(airT_legend, airT_plot.2, soilT_plot, snowDepth_plot.2, snowData_legend, ncol = 1, widths = c(2.7), heights = c(0.5, 3, 3, 3, 0.5))#, top = grid::textGrob('Mean diel temperature', gp=grid::gpar(fontsize=20)))
#
#
# Snow and moisture in one graph
snowDepth_plot_2 <- avgVWC_wide3 %>%
  # mutate(Site = case_when(Site == "Abisko" ~ "A",
  #                         Site == "Vassijaure" ~ "V",
  #                         TRUE ~ Site)) %>%
  rename("SiteGWC" = Site) %>%
  left_join(snowData_2, by = join_by(Date)) %>%
  ggplot() +
  #
  # GWC
  # geom_point(aes(x = Date, y = Soil_GWC_FW, shape = SiteGWC), na.rm = TRUE) +
  # geom_errorbar(aes(x = Date, y = Soil_GWC_FW, ymin=Soil_GWC_FW-se_GWC, ymax=Soil_GWC_FW+se_GWC), position=position_dodge(.9)) +
  # VWC converted from GWC
  geom_point(aes(x = Date, y = Soil_VWC_est, shape = SiteGWC), na.rm = TRUE, size = 2) +
  geom_errorbar(aes(x = Date, y = Soil_VWC_est, ymin=Soil_VWC_est-se_VWC, ymax=Soil_VWC_est+se_VWC), position=position_dodge(.9)) +
  #
  # Snow
  geom_line(aes(x = Date, y = Snow_depth_cm, linetype = Site), na.rm = TRUE, linewidth = 1) +
  geom_ribbon(aes(x = Date, y = Snow_depth_cm, ymin = min, ymax = max, fill = Site, linetype = Site), alpha = 0.5) +
  scale_fill_grey(na.translate = F) +
  #
  # 0-line
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  #
  # VWC from sensors
  geom_line(aes(x = Date, y = Abisko_VWC, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  geom_line(aes(x = Date, y = Vassijaure_VWC, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) +
  #
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Soil Moisture (% vol)"))+#, minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%d-%b") +
  #scale_shape(na.translate = F) +
  scale_shape_manual(values = c(19, 1), na.translate = F) +
  coord_cartesian(xlim = c(plotDayStart,plotDayEnd)) +
  labs(x = "Time of year", y = "Snow cover (cm)") +
  guides(fill = guide_legend(title = "", order = 1), shape = guide_legend(title = "Soil GWC", order = 2), linetype = "none") +
  theme_bw(base_size = 25) +
  theme(legend.position = "bottom")
#
snowData_legend_2 <- get_plot_component(snowDepth_plot_2, "guide-box", return_all = TRUE)#[[3]]  # 1 is right, 2 is left, 3 is bottom, 4 is top
snowData_yaxis_2 <- get_y_axis(snowDepth_plot_2)
snowDepth_plot_2.2 <- snowDepth_plot_2 + theme_bw(base_size = 25) + theme(legend.position = "none")#, axis.title.y = element_blank())
#
#
# Measuring periods
MP_plot <- avgT_wide2 %>%
  ggplot() +
  #annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  #annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[1]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[1]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[4]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[4]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[5]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[5]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[6]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[6]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[7]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[7]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[8]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[8]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[9]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[9]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[10]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[10]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[11]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[11]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[12]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[12]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[13]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[13]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[14]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[14]), yend = 10), linewidth = 5, fill = "#000000") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[15]), y = 10, xend = as.Date(DayOfLH$Day_of_harvest[15]), yend = 10), linewidth = 5, fill = "#000000") +
  scale_y_continuous(breaks = c(0, 10), minor_breaks = c(0, 10)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%d-%b") +
  coord_cartesian(xlim = c(plotDayStart,plotDayEnd)) +
  labs(x = element_blank(), y = "MP") +
  theme_bw(base_size = 25) +
  theme(axis.title.x = element_blank())
#
#
# Layout of final graph:
hlay <- rbind(c(1,1),
              c(2,NA),
              c(3,NA),
              c(4,NA),
              c(5,5),
              c(6,6))
#
# Final graph:
# Fig. 2: Air temperature, Soil temperature, and Snow depth with soil moisture as VWC and converted GWC
env_plot <- grid.arrange(airT_legend, airT_plot.2, soilT_plot.2, MP_plot, snowDepth_plot_2.2, snowData_legend_2, widths = c(2.8,0.15), heights = c(0.7, 4, 4, 1.3, 4.7, 0.7), layout_matrix=hlay)

ggsave("Fig_2_EnvironmentalData.4.png", plot = env_plot, path = "images", width = 50, height = 27, units = "cm", dpi = 1200, bg = "white")

#env_plot.flux <- plot_grid(airT_plot, soilT_plot.2, soilM_plot.2, PAR_plot.2, Flux_period_plot, align = "v", ncol = 1, rel_heights = c(3,3,2.5,3.5,1.5))
#plot_grid(env_plot.flux, soilT_legend, ncol = 1, rel_heights = c(9, 1))
#
#
# Save important data from environmental:
# write_tsv(avgT_wide2, "export/Temperature_Air_Soil.tsv")
# write_tsv(snowData, "export/Snow.tsv")
# write_tsv(snowData_2, "export/Snow_avg.tsv")
# write_tsv(avgVWC_wide3, "export/VWC_avg.tsv")
#
# For quick access to make the graphs above:
# avgT_wide2 <- read_tsv("export/Temperature_Air_Soil.tsv")
# snowData <- read_tsv("export/Snow.tsv")
# snowData_2 <- read_tsv("export/Snow_avg.tsv")
# avgVWC_wide3 <- read_tsv("export/VWC_avg.tsv")
#
#
# <><><><><> END FIG 1 <><><><><>
#
#
#
#
#------- The End -------