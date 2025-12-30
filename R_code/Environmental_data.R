# WinterEcology I experiment
# Script author: Emil A.S. Andersen
#
# Load weather data
# From own data
#
# To get months in the right format (i.e. not in whatever local the computer has, e.g. Swedish)
Sys.setlocale("LC_ALL", 'en_GB.UTF-8')
#
#library(plyr)
library(tidyverse)
library(readr)
library(lubridate)
library(readxl)
library(viridis)
#library(dygraphs)
#library(xts)
#library(plotly)
#library(hrbrthemes)
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
SMHI_Abisko_Tair <- read_csv2("raw_data/smhi_Temp_Abisko_new.csv", skip = 10)
SMHI_Katterjakk_Tair <- read_csv2("raw_data/smhi_Kattejakk_Temp.csv", skip = 9)
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
# How the temperature should be calculated:
test <- Abisko_Tair %>%
  rowwise() %>%
  mutate(Tair_A3 = mean(c(Tair_A39_1, Tair_C1, Tair_31, Tair_39_2), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(Tair_A3))
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
avgT_wide <- left_join(avgTair_wide, avgTsoil_wide)
avgT_long <- left_join(avgTair_long, avgTsoil_long)
#
# Add confidence interval
# 95% CI
test1 <- test %>%
  dplyr::select("Date_A", "Tair_A3") %>%
  mutate(Date_A = date(Date_A))

test2 <- summarySE(test1, measurevar = "Tair_A3", groupvars = c("Date_A"))

test2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Abisko snow
  geom_line(aes(x = Date_A, y = Tair_A3, lty = "Abisko"), na.rm = TRUE) + 
  geom_ribbon(aes(x = Date_A, y = Tair_A3, ymin = Tair_A3-ci, ymax = Tair_A3+ci), alpha = 0.5) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time", y = "Air temperature (°C)") + # x = "Time of year",  , title = "Air temperature" 
  guides(lty = guide_legend(title = "Mean diel temperature")) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")


# test3 <- avgT_wide2 %>%
#   dplyr::filter(Date >= winterP_date$wstart[1]) %>%
#   dplyr::filter(Date <= winterP_date$wend[2])

# avgT_wide2 %>%
#   ggplot(aes(x = Abisko_Tair, y = Abisko_Tsoil)) +
#   geom_point()
# test3 %>%
#   ggplot(aes(x = Vassijaure_Tair, y = Vassijaure_Tsoil)) +
#   geom_point()

#
#
#
#------- # Soil Moisture # -------
#
# Check all sensors
# Abisko
Abisko_EM50 %>%
  select(Date, A1N_VWC, A2N_VWC, A3N_VWC, A4N_VWC, A5N_VWC) %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "#999999") +
  geom_line(aes(x = date(Date), y = A1N_VWC, lty = "A1"), na.rm = TRUE) +
  geom_line(aes(x = date(Date), y = A2N_VWC, lty = "A2"), na.rm = TRUE) +
  geom_line(aes(x = date(Date), y = A3N_VWC, lty = "A3"), na.rm = TRUE) +
  geom_line(aes(x = date(Date), y = A4N_VWC, lty = "A4"), na.rm = TRUE) +
  geom_line(aes(x = date(Date), y = A5N_VWC, lty = "A5"), na.rm = TRUE) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time of year", y = expression("Soil VWC ("*m^3*" "*m^-3*")"), title = "Abisko soil VWC") + 
  guides(lty = "none") +
  theme_bw(base_size = 17) +
  theme(legend.position = "top")
#
# Vassijaure
Vassijaure_EM50 %>%
  select(V_Date, V1N_VWC, V2N_VWC, V3N_VWC, V4N_VWC, V5N_VWC) %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "#999999") +
  geom_line(aes(x = date(V_Date), y = V1N_VWC, lty = "V1"), na.rm = TRUE) +
  geom_line(aes(x = date(V_Date), y = V2N_VWC, lty = "V2"), na.rm = TRUE) +
  geom_line(aes(x = date(V_Date), y = V3N_VWC, lty = "V3"), na.rm = TRUE) +
  geom_line(aes(x = date(V_Date), y = V4N_VWC, lty = "V4"), na.rm = TRUE) +
  geom_line(aes(x = date(V_Date), y = V5N_VWC, lty = "V5"), na.rm = TRUE) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time of year", y = expression("Soil VWC ("*m^3*" "*m^-3*")"), title = "Vassijaure soil VWC") + 
  guides(lty = "none") +
  theme_bw(base_size = 17) +
  theme(legend.position = "top")
#
#
# Mean Soil temperature
Abisko_avgVWC <- Abisko_EM50 %>%
  select(Date, A1N_VWC, A2N_VWC, A3N_VWC, A4N_VWC, A5N_VWC) %>%
  # Remove soil temperatures below -15 C
  # mutate(A1N_VWC = replace(A1N_VWC, A1N_VWC < -15, NA), # Has very low values in July 2019 and during winter
  #        A2N_VWC = replace(A1N_VWC, A2N_VWC < -15, NA),
  #        A3N_VWC = replace(A1N_VWC, A3N_VWC < -15, NA), # Has a single measure very low measure in March 2020
  #        A4N_VWC = replace(A1N_VWC, A4N_VWC < -15, NA),
  #        A5N_VWC = replace(A1N_VWC, A5N_VWC < -15, NA)) %>%
  # mutate(A1N_VWC = if_else(month(Date) == 7 & A1N_VWC < 5, NA, A1N_VWC),
  #        A2N_VWC = if_else(month(Date) == 7 & A2N_VWC < 0, NA, A1N_VWC),
  #        A3N_VWC = if_else(month(Date) == 7 & A3N_VWC < 0, NA, A1N_VWC),
  #        A4N_VWC = if_else(month(Date) == 7 & A4N_VWC < 0, NA, A1N_VWC),
  #        A5N_VWC = if_else(month(Date) == 7 & A5N_VWC < 0, NA, A1N_VWC)) %>%
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
  # Remove soil temperatures below -15 C
  # mutate(V1N_VWC = replace(V1N_VWC, V1N_VWC < -15, NA), # Has very low values in July 2019 and during winter
  #        V2N_VWC = replace(V2N_VWC, V2N_VWC < -15, NA),
  #        V3N_VWC = replace(V3N_VWC, V3N_VWC < -15, NA), # Has a single measure very low measure in March 2020
  #        V4N_VWC = replace(V4N_VWC, V4N_VWC < -15, NA),
  #        V5N_VWC = replace(V5N_VWC, V5N_VWC < -15, NA)) %>%
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
# Format wide: Abisko and Vassijaureseparate columns
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







# First attempt at graphing

avgVWC_wide2 %>% 
  ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999") +
  geom_line(aes(x = Date, y = Vassijaure_VWC, lty = "Vassijaure"), na.rm = TRUE) +
  geom_line(aes(x = Date, y = Abisko_VWC, lty = "Abisko"), na.rm = TRUE) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15), sec.axis = sec_axis(~.*1, name = "Soil Moisture (% by volume)"))+#, minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = NULL, y = "Snow cover (cm)") +
  guides(lty = "none") +
  theme_bw(base_size = 17) +
  theme(legend.position = "top")



# Second attempt. Incl snow depth
# 
# avgVWC_wide2 %>% 
#   left_join(snowData_2, by = join_by(Date)) %>%
#   ggplot() +
#   geom_line(aes(x = Date, y = Snow_depth_cm, linetype = Site)) +
#   geom_ribbon(aes(x = Date, y = Snow_depth_cm, ymin = min, ymax = max, fill = Site, linetype = Site), alpha = 0.5) +
#   #geom_point(aes(x = Date, y = Snow_depth_cm, shape = Site)) +
#   #geom_errorbar(aes(x = Date, y = Snow_depth_cm, ymin=min, ymax=max), position=position_dodge(.9)) +
#   scale_fill_grey() +
#   geom_hline(yintercept = 0, color = "#999999") +
#   geom_line(aes(x = Date, y = Vassijaure_VWC, lty = "Vassijaure"), na.rm = TRUE) +
#   geom_line(aes(x = Date, y = Abisko_VWC, lty = "Abisko"), na.rm = TRUE) +
#   scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Soil Moisture (% by volume)"))+#, minor_breaks = c(-15, -5, 5, 15)) +
#   scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
#   coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
#   labs(x = "Time of year", y = "Snow cover (cm)") +
#   guides(fill = guide_legend(title = "Snow"), linetype = "none") +
#   theme_bw(base_size = 25) +
#   theme(legend.position = "bottom")








# avgVWC_wide3 %>% 
#   left_join(snowData_2, by = join_by(Date, Site)) %>%
#   ggplot() +
#   geom_line(aes(x = Date, y = Snow_depth_cm, linetype = Site)) +
#   geom_ribbon(aes(x = Date, y = Snow_depth_cm, ymin = min, ymax = max, fill = Site, linetype = Site), alpha = 0.5) +
#   #geom_point(aes(x = Date, y = Snow_depth_cm, shape = Site)) +
#   #geom_errorbar(aes(x = Date, y = Snow_depth_cm, ymin=min, ymax=max), position=position_dodge(.9)) +
#   scale_fill_grey() +
#   geom_hline(yintercept = 0, color = "#999999") +
#   geom_line(aes(x = Date, y = Vassijaure_VWC, lty = "Vassijaure"), na.rm = TRUE) +
#   geom_line(aes(x = Date, y = Abisko_VWC, lty = "Abisko"), na.rm = TRUE) +
#   scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Soil Moisture (% by volume)"))+#, minor_breaks = c(-15, -5, 5, 15)) +
#   scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
#   coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
#   labs(x = "Time of year", y = "Snow cover (cm)") +
#   guides(fill = guide_legend(title = "Snow"), linetype = "none") +
#   theme_bw(base_size = 25) +
#   theme(legend.position = "bottom")




# x <- left_join(coreGWC_1, snowData_2, by = join_by(Date, Site))




# avgVWC_wide3 %>%
#   mutate(Site = case_when(Site == "Abisko" ~ "A",
#                           Site == "Vassijaure" ~ "V",
#                           TRUE ~ Site)) %>%
#   rename("SiteGWC" = Site) %>%
#   left_join(snowData_2, by = join_by(Date)) %>%
#   ggplot() +
#   geom_point(aes(x = Date, y = Soil_VWC_est, shape = SiteGWC), na.rm = TRUE) +
#   geom_errorbar(aes(x = Date, y = Soil_VWC_est, ymin=Soil_VWC_est-se_VWC, ymax=Soil_VWC_est+se_VWC), position=position_dodge(.9)) +
#   geom_line(aes(x = Date, y = Snow_depth_cm, linetype = Site), na.rm = TRUE) +
#   geom_ribbon(aes(x = Date, y = Snow_depth_cm, ymin = min, ymax = max, fill = Site, linetype = Site), alpha = 0.5) +
#   #geom_point(aes(x = Date, y = Snow_depth_cm, shape = Site)) +
#   
#   scale_fill_grey() +
#   geom_hline(yintercept = 0, color = "#999999") +
#   #geom_line(aes(x = Date, y = dielVWC_soil, linetype = Site), na.rm = TRUE) +
#   geom_line(aes(x = Date, y = Abisko_VWC, lty = "Abisko"), na.rm = TRUE) +
#   geom_line(aes(x = Date, y = Vassijaure_VWC, lty = "Vassijaure"), na.rm = TRUE) +
#   #geom_point(aes(x = Date, y = Abisko_GWC, shape = "Abisko"), na.rm = TRUE) +
#   #geom_point(aes(x = Date, y = Vassijaure_GWC, shape = "Vassijaure"), na.rm = TRUE) +
#   scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Soil Moisture (% by volume)"))+#, minor_breaks = c(-15, -5, 5, 15)) +
#   scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
#   coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
#   labs(x = "Time of year", y = "Snow cover (cm)") +
#   guides(fill = guide_legend(title = "Snow"), shape = guide_legend(title = "Soil GWC"), linetype = "none") +
#   theme_bw(base_size = 25) +
#   theme(legend.position = "bottom")








# Potential graph update:



# avgVWC_wide3 %>%
#   mutate(Site = case_when(Site == "Abisko" ~ "A",
#                           Site == "Vassijaure" ~ "V",
#                           TRUE ~ Site)) %>%
#   rename("SiteGWC" = Site) %>%
#   left_join(snowData_2, by = join_by(Date)) %>%
#   ggplot() +
#   #
#   # GWC
#   # geom_point(aes(x = Date, y = Soil_GWC_FW, shape = SiteGWC), na.rm = TRUE) +
#   # geom_errorbar(aes(x = Date, y = Soil_GWC_FW, ymin=Soil_GWC_FW-se_GWC, ymax=Soil_GWC_FW+se_GWC), position=position_dodge(.9)) +
#   # VWC converted from GWC
#   geom_point(aes(x = Date, y = Soil_VWC_est, shape = SiteGWC), na.rm = TRUE) +
#   geom_errorbar(aes(x = Date, y = Soil_VWC_est, ymin=Soil_VWC_est-se_VWC, ymax=Soil_VWC_est+se_VWC), position=position_dodge(.9)) +
#   #
#   # Snow
#   geom_line(aes(x = Date, y = Snow_depth_cm, linetype = Site), na.rm = TRUE) +
#   geom_ribbon(aes(x = Date, y = Snow_depth_cm, ymin = min, ymax = max, fill = Site, linetype = Site), alpha = 0.5) +
#   scale_fill_grey(na.translate = F) +
#   #
#   # 0-line
#   geom_hline(yintercept = 0, color = "#999999") +
#   #
#   # VWC from sensors
#   geom_line(aes(x = Date, y = Abisko_VWC, lty = "Abisko"), na.rm = TRUE) +
#   geom_line(aes(x = Date, y = Vassijaure_VWC, lty = "Vassijaure"), na.rm = TRUE) +
#   scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Soil Moisture (% by volume)"))+#, minor_breaks = c(-15, -5, 5, 15)) +
#   scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
#   scale_shape(na.translate = F) +
#   coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
#   labs(x = "Time of year", y = "Snow cover (cm)") +
#   guides(fill = guide_legend(title = "Snow"), shape = guide_legend(title = "Soil GWC"), linetype = "none") +
#   theme_bw(base_size = 25) +
#   theme(legend.position = "bottom")













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
  mutate(Abisko = mean(c(A1_PAR, A2_PAR, A3_PAR, A4_PAR, A5_PAR), na.rm = TRUE))
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
  mutate(Vassijaure = mean(c(V1_PAR, V2_PAR, V3_PAR, V4_PAR, V5_PAR), na.rm = TRUE))

  
  # mutate(A1_PAR = if_else(month(Date) == 7 & A1_PAR < 5, NA, A1_PAR),
  #        A2_PAR = if_else(month(Date) == 7 & A2_PAR < 0, NA, A1_PAR),
  #        A3_PAR = if_else(month(Date) == 7 & A3_PAR < 0, NA, A1_PAR),
  #        A4_PAR = if_else(month(Date) == 7 & A4_PAR < 0, NA, A1_PAR),
  #        A5_PAR = if_else(month(Date) == 7 & A5_PAR < 0, NA, A1_PAR)) #%>%
  
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
         !(year(Datum) == 2019 & month(Datum) < 7))
#
SMHI_A_avgTair <- SMHI_Abisko_Tair %>%
  group_by(date(Datum)) %>%
  summarise(Abisko_Tair_SMHI = mean(Lufttemperatur, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Datum)")
#
# Katterjåkk (closest station to Vassijaure)
SMHI_Katterjakk_Tair <- SMHI_Katterjakk_Tair %>%
  dplyr::select(1:3) %>%
  mutate(across(Datum, ymd), # If also for time: across("Tid (UTC)", hms))
         across(Lufttemperatur, as.numeric)) %>% 
  dplyr::filter(year(Datum) == 2019 | year(Datum) == 2020) %>%
  dplyr::filter(!(year(Datum) == 2020 & month(Datum) > 10),
         !(year(Datum) == 2019 & month(Datum) < 7))
#
SMHI_K_avgTair <- SMHI_Katterjakk_Tair %>%
  group_by(date(Datum)) %>%
  summarise(Katterjakk_Tair_SMHI = mean(Lufttemperatur, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Datum)")
#
# Combine

# Combine all
avgT_wide2 <- left_join(avgT_wide, SMHI_A_avgTair) %>%
  left_join(SMHI_K_avgTair)
#
# Pivot long
avgT_long2 <- avgT_wide2 %>%
  pivot_longer(2:7, names_to = "Site", values_to = "dielT")
#
#
#
#------- ### Plot ### -------

# Plot it
avgTair_long %>%
  ggplot(aes(x = Date, y = dielT_air, colour = Site)) +
  geom_point() +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  xlab("Time") + ylab("Temperature C")
#
avgTsoil_long %>%
  ggplot(aes(x = Date, y = dielT_soil, colour = Site)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  xlab("Time") + ylab("Temperature C")
#

# For getting the labels, use aes of e.g. linetype/shape, color or fill
#
# ggplot() + geom_rect(winterP, aes(x=, ymin=-Inf, ymax=Inf), stat="identity", alpha = 0.5, fill = 'grey', inherit.aes = FALSE)
# Vassijaure air and soil
avgT_wide2 %>% ggplot() +
  geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Katterjåkk air temperature")) + 
  geom_line(aes(x = Date, y = Vassijaure_Tair, lty = "Vassijaure air temperature"), na.rm = TRUE) + 
  geom_point(aes(x = Date, y = Vassijaure_Tsoil, fill = "Vassijaure soil temperature"), shape = 6, na.rm = TRUE) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time of year", y = "Mean diel temperature °C", title = "Air and soil temperature") +
  guides(fill = guide_legend(title = "Soil temperature"), lty = guide_legend(title = "Air temperature")) +
  theme_bw(base_size = 15)
  #theme(axis.text.x=element_text(angle=10, hjust=0.5))
#
# Abisko air and soil
avgT_wide2 %>% ggplot() +
  geom_line(aes(x = Date, y = Abisko_Tair_SMHI, lty = "Abisko air temperature SMHI")) +
  geom_line(aes(x = Date, y = Abisko_Tair, lty = "Abisko air temperature")) +
  geom_point(aes(x = Date, y = Abisko_Tsoil, fill = "Abisko soil temperature"), shape = 6, na.rm = TRUE) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time of year", y = "Air temperature °C", title = "Air and soil temperature") +
  guides(fill = guide_legend(title = "Soil temperature"), lty = guide_legend(title = "Air temperature")) +
  theme_bw(base_size = 15)
#
#
# <><><><><> MAIN PLOT - FIG 1 <><><><><>
#
#
# Snow depth on plot
snowData_2 <- summarySE(snowData, measurevar = "Snow_depth_cm", groupvars = c("Site", "MP", "Date"))
snowData_2 <- snowData_2 %>%
  mutate(Date = ymd(Date))
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
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-08-26"))) +
  labs(x = NULL, y = "Air (°C)") + # x = "Time of year",  , title = "Air temperature" 
  guides(lty = guide_legend(title = "Mean diel temperature"))+ #lty = guide_legend(title = "Mean diel temperature")) +
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
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-08-26")), ylim = c(-5,15)) +
  labs(x = NULL, y = "Soil (°C)") + # x = "Time of year", , title = "Soil temperature"
  guides(lty = "none") + # guide_legend(title = "Soil temperature")
  theme_bw(base_size = 25) +
  theme(legend.position = "top")#, axis.title.y = element_text(size = 25))
#
# soilT_legend <- get_x_axis(soilT_plot) + theme_bw(base_size = 17)
# soilT_plot.2 <- soilT_plot + theme_bw(base_size = 20) + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 15)) 
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
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-08-26"))) +
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
  mutate(Site = case_when(Site == "Abisko" ~ "A",
                          Site == "Vassijaure" ~ "V",
                          TRUE ~ Site)) %>%
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
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-08-26"))) +
  labs(x = "Time of year", y = "Snow cover (cm)") +
  guides(fill = guide_legend(title = "Snow", order = 1), shape = guide_legend(title = "Soil GWC", order = 2), linetype = "none") +
  theme_bw(base_size = 25) +
  theme(legend.position = "bottom")
#
snowData_legend_2 <- get_plot_component(snowDepth_plot_2, "guide-box", return_all = TRUE)#[[3]]  # 1 is right, 2 is left, 3 is bottom, 4 is top
snowData_yaxis_2 <- get_y_axis(snowDepth_plot_2)
snowDepth_plot_2.2 <- snowDepth_plot_2 + theme_bw(base_size = 25) + theme(legend.position = "none")#, axis.title.y = element_blank()) 
#
#
# Layout of final graph:
hlay <- rbind(c(1,1),
              c(2,NA),
              c(3,NA),
              c(4,4),
              c(5,5))
#
# Final graph:
# Fig. 1: Air temperature, Soil temperature, and Snow depth with soil moisture as VWC and converted GWC
grid.arrange(airT_legend, airT_plot.2, soilT_plot, snowDepth_plot_2.2, snowData_legend_2, widths = c(2.8,0.15), heights = c(0.5, 3, 3, 3.5, 0.5), layout_matrix=hlay)
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



soilT_plot.MP <- soilT_plot +
  # From label to harvest
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[1]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[1]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[2]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[2]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[3]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[3]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[4]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[4]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[5]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[5]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[6]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[6]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[7]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[7]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[8]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[8]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[9]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[9]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[10]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[10]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[11]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[11]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[12]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[12]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[13]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[13]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[14]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[14]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[15]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[15]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3)

airT_plot.2.MP <- airT_plot.2 +
  # From label to harvest
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[1]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[1]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[2]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[2]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[3]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[3]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[4]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[4]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[5]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[5]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[6]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[6]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[7]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[7]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[8]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[8]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[9]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[9]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[10]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[10]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[11]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[11]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[12]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[12]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[13]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[13]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[14]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[14]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[15]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[15]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3)

snowDepth_plot_2.2.MP <- snowDepth_plot_2.2 +
  # From label to harvest
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[1]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[1]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[2]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[2]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[3]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[3]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[4]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[4]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[5]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[5]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[6]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[6]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[7]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[7]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[8]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[8]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[9]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[9]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[10]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[10]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[11]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[11]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[12]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[12]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[13]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[13]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[14]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[14]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3) +
  annotate("rect", xmin = as.Date(DayOfLH$Day_of_label[15]), ymin = -Inf, xmax = as.Date(DayOfLH$Day_of_harvest[15]), ymax = Inf, linewidth = 0.9, fill = "#D55E00", alpha = 0.3)


# Final graph 2, with dates of measurements:
# Fig. 1: Air temperature, Soil temperature, and Snow depth with soil moisture as VWC and converted GWC
grid.arrange(airT_legend, airT_plot.2.MP, soilT_plot.MP, snowDepth_plot_2.2, snowData_legend_2, widths = c(2.8,0.15), heights = c(0.5, 3, 3, 3.5, 0.5), layout_matrix=hlay)




# Air temperatures - all
airT_plot2 <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999") +
  #geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Katterjakk air temperature")) + 
  geom_line(aes(x = Date, y = Vassijaure_Tair, lty = "Vassijaure"), na.rm = TRUE) + 
  #geom_line(aes(x = Date, y = Abisko_Tair_SMHI, lty = "Abisko air temperature SMHI")) +
  geom_line(aes(x = Date, y = Abisko_Tair, lty = "Abisko"), na.rm = TRUE) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = NULL, y = "Air temperature (°C)") + # x = "Time of year",  , title = "Air temperature" 
  guides(lty = guide_legend(title = "Mean diel temperature")) +
  theme_classic(base_size = 15) +
  theme(legend.position = "top", axis.text.x = element_blank())
#
# Soil temperatures - all
soilT_plot2 <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999") +
  
  # Arrows from label to harvest
  #geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  #geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  
  # Label arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[1]), y = -3, xend = as.Date(DayOfLH$Day_of_label[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = -3, xend = as.Date(DayOfLH$Day_of_label[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = -3, xend = as.Date(DayOfLH$Day_of_label[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[4]), y = -3, xend = as.Date(DayOfLH$Day_of_label[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[5]), y = -3, xend = as.Date(DayOfLH$Day_of_label[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[6]), y = 3, xend = as.Date(DayOfLH$Day_of_label[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[7]), y = 3, xend = as.Date(DayOfLH$Day_of_label[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[8]), y = 3, xend = as.Date(DayOfLH$Day_of_label[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[9]), y = 3, xend = as.Date(DayOfLH$Day_of_label[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[10]), y = 3, xend = as.Date(DayOfLH$Day_of_label[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[11]), y = 3, xend = as.Date(DayOfLH$Day_of_label[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[12]), y = 3, xend = as.Date(DayOfLH$Day_of_label[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[13]), y = -3, xend = as.Date(DayOfLH$Day_of_label[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[14]), y = -3, xend = as.Date(DayOfLH$Day_of_label[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[15]), y = -3, xend = as.Date(DayOfLH$Day_of_label[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  
  # Harvest arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[1]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[2]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[3]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[4]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[5]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[6]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[7]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[8]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[9]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[10]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[11]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[12]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[13]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[14]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[15]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  
  
  geom_line(aes(x = Date, y = Vassijaure_Tsoil, lty = "Vassijaure"), na.rm = TRUE) +
  geom_line(aes(x = Date, y = Abisko_Tsoil, lty = "Abisko"), na.rm = TRUE) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15))+#, minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16")), ylim = c(-5,15)) +
  labs(x = NULL, y = "Soil temperature (°C)") + # x = "Time of year", , title = "Soil temperature"
  guides(lty = "none") + # guide_legend(title = "Soil temperature")
  theme_classic(base_size = 15) +
  theme(legend.position = "top")
#
# Snow depth on plot
snowDepth_plot2 <- snowData_2 %>%
  ggplot(aes(x = Date, y = Snow_depth_cm, ymin = Snow_depth_cm-ci, ymax = Snow_depth_cm+ci, fill = Site)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_grey() +
  #scale_fill_viridis_d() +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time of year", y = "Snow cover (cm)") + #, title = "Snow cover measured over the entire plot (around all 15 patches)") +
  guides(fill = guide_legend(title = "Snow")) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom")
#
# Plot
grid.arrange(airT_plot2, soilT_plot2, snowDepth_plot2)






# Air temperatures - all
airT_plot4 <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  #geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Katterjakk air temperature")) + 
  geom_line(aes(x = Date, y = Vassijaure_Tair, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) + 
  #geom_line(aes(x = Date, y = Abisko_Tair_SMHI, lty = "Abisko air temperature SMHI")) +
  geom_line(aes(x = Date, y = Abisko_Tair, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-07-12"),as.Date("2020-08-26"))) +
  labs(x = NULL, y = "Air (°C)") + # x = "Time of year",  , title = "Air temperature" 
  guides(lty = guide_legend(title = "Mean diel temperature"))+ #lty = guide_legend(title = "Mean diel temperature")) +
  theme_bw(base_size = 25) +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.text.y = element_text(size = 15))
#
airT_legend4 <- get_legend(airT_plot4)
airT_plot4.2 <- airT_plot4 + theme_bw(base_size = 25) + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 18))#, axis.title.y = element_text(size = 25)) 
#airT_plot <- airT_plot + guides(lty = NULL)




soilT_plot4 <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  #annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -3, ymax = 0, fill = "white", alpha = 0.9) +
  #geom_hline(yintercept = -3, color = "#D55E00") +
  
  # Arrows from label to harvest
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[1]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[1]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[4]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[4]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[5]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[5]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[6]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[6]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[7]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[7]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[8]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[8]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[9]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[9]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[10]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[10]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[11]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[11]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[12]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[12]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[13]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[13]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[14]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[14]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[15]), y = -5, xend = as.Date(DayOfLH$Day_of_harvest[15]), yend = -5), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  
  # Label arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[1]), y = -2, xend = as.Date(DayOfLH$Day_of_label[1]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = -2, xend = as.Date(DayOfLH$Day_of_label[2]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = -2, xend = as.Date(DayOfLH$Day_of_label[3]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[4]), y = -2, xend = as.Date(DayOfLH$Day_of_label[4]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[5]), y = -2, xend = as.Date(DayOfLH$Day_of_label[5]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[6]), y = -2, xend = as.Date(DayOfLH$Day_of_label[6]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[7]), y = -2, xend = as.Date(DayOfLH$Day_of_label[7]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[8]), y = -2, xend = as.Date(DayOfLH$Day_of_label[8]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[9]), y = -2, xend = as.Date(DayOfLH$Day_of_label[9]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[10]), y = -2, xend = as.Date(DayOfLH$Day_of_label[10]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[11]), y = -2, xend = as.Date(DayOfLH$Day_of_label[11]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[12]), y = -2, xend = as.Date(DayOfLH$Day_of_label[12]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[13]), y = -2, xend = as.Date(DayOfLH$Day_of_label[13]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[14]), y = -2, xend = as.Date(DayOfLH$Day_of_label[14]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[15]), y = -2, xend = as.Date(DayOfLH$Day_of_label[15]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  
  # Harvest arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[1]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[1]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[2]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[3]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[4]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[4]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[5]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[5]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[6]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[6]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[7]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[7]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[8]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[8]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[9]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[9]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[10]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[10]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[11]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[11]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[12]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[12]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[13]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[13]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[14]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[14]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[15]), y = -2, xend = as.Date(DayOfLH$Day_of_harvest[15]), yend = -5), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  
  geom_line(aes(x = Date, y = Vassijaure_Tsoil, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) +
  geom_line(aes(x = Date, y = Abisko_Tsoil, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15))+#, minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%d-%b") +
  coord_cartesian(xlim = c(as.Date("2019-07-12"),as.Date("2020-08-26")), ylim = c(-7,15)) +
  labs(x = NULL, y = "Soil (°C)") + # x = "Time of year", , title = "Soil temperature"
  guides(lty = "none") + # guide_legend(title = "Soil temperature")
  theme_bw(base_size = 25) +
  theme(legend.position = "top")#, axis.title.y = element_text(size = 25))

#
grid.arrange(airT_legend, soilT_plot4, widths = c(2.8), heights = c(0.5, 3))


grid.arrange(airT_legend4, airT_plot4.2, soilT_plot4, widths = c(2.8), heights = c(0.5, 3, 3))








#
library(car)
library(nlme)
#
snowData_stat <- snowData %>%
  add_row(Date = 20200603, MP = 13, Site = "Abisko", Plot = 1, Snow_depth_cm = 0) %>% # Add Abisko for MP13 where it is snow-free
  add_row(Date = 20200603, MP = 13, Site = "Abisko", Plot = 2, Snow_depth_cm = 0) %>%
  add_row(Date = 20200603, MP = 13, Site = "Abisko", Plot = 3, Snow_depth_cm = 0) %>%
  add_row(Date = 20200603, MP = 13, Site = "Abisko", Plot = 4, Snow_depth_cm = 0) %>%
  add_row(Date = 20200603, MP = 13, Site = "Abisko", Plot = 5, Snow_depth_cm = 0) %>%
  add_row(Date = 20200624, MP = 13, Site = "Abisko", Plot = 1, Snow_depth_cm = 0) %>%
  add_row(Date = 20200624, MP = 13, Site = "Abisko", Plot = 2, Snow_depth_cm = 0) %>%
  add_row(Date = 20200624, MP = 13, Site = "Abisko", Plot = 3, Snow_depth_cm = 0) %>%
  add_row(Date = 20200624, MP = 13, Site = "Abisko", Plot = 4, Snow_depth_cm = 0) %>%
  add_row(Date = 20200624, MP = 13, Site = "Abisko", Plot = 5, Snow_depth_cm = 0) %>%
  mutate(MP = case_when(Date == 20200504 | Date == 20200506 ~ 12.1,
                        Date == 20200525 | Date == 20200527 ~ 12.2,
                        Date == 20200601 | Date == 20200603 ~ 13.1,
                        Date == 20200622 | Date == 20200624 ~ 13.2,
                        TRUE ~ MP)) %>%
  mutate(Round = case_when(MP == 5 ~ "05_Nov_2019",
                           MP == 6 ~ "06_Dec_2019",
                           MP == 7 ~ "07_Jan_2020",
                           MP == 8 ~ "08_Feb_2020",
                           MP == 9 ~ "09_Mar_2020",
                           MP == 10 ~ "10_Apr_2020",
                           MP == 11 ~ "11_Apr_2020",
                           MP == 12.1 ~ "12_Maj_2020_injection",
                           MP == 12.2 ~ "12_Maj_2020_harvest",
                           MP == 13.1 ~ "13_Jun_2020_injection",
                           MP == 13.2 ~ "13_Jun_2020_harvest")) #%>%
  filter(MP != 13.1 & MP != 13.2) # Remove round 13 as Abisko i snow-free by then
lmeSnow<-lme(asin(sqrt(Snow_depth_cm/100)) ~ Site*Round,
             random = ~1|Plot/Site,
             data = snowData_stat, na.action = na.exclude , method = "REML")
#
#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lmeSnow), resid(lmeSnow), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lmeSnow), main = "Normally distributed?")                 
qqline(resid(lmeSnow), main = "Homogeneity of Variances?", col = 2) #OK
plot(lmeSnow)
par(mfrow = c(1,1))
#
#model output
Anova(lmeSnow, type=2)



# Everything - chaos
avgT_wide2 %>% ggplot() +
  geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Katterjakk air temperature")) + 
  geom_line(aes(x = Date, y = Vassijaure_Tair, lty = "Vassijaure air temperature")) + 
  geom_point(aes(x = Date, y = Vassijaure_Tsoil, fill = "Vassijaure soil temperature"), shape = 6) +
  geom_line(aes(x = Date, y = Abisko_Tair_SMHI, lty = "Abisko air temperature SMHI")) +
  geom_line(aes(x = Date, y = Abisko_Tair, lty = "Abisko air temperature")) +
  geom_point(aes(x = Date, y = Abisko_Tsoil, fill = "Abisko soil temperature"), shape = 4) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time of year", y = "Mean diel temperature °C", title = "Air and soil temperature") +
  guides(fill = guide_legend(title = "Soil temperature"), lty = guide_legend(title = "Air temperature")) +
  theme_bw(base_size = 15)


avgT_wide2 %>%
  ggplot() +
  geom_line(aes(x = Date, y = Abisko_Tair_SMHI), color = "#E69F00", linetype = "solid", size=1, alpha = 0.4) + #                                                                                                         # Orange
  geom_point(aes(x = Date, y = Abisko_Tair_SMHI), shape = 6) +
  geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Katterjakk"), color = "#009E73", linetype = "solid", size=1, alpha = 0.4) + #                                                                                  # Bluish green
  geom_point(aes(x = Date, y = Katterjakk_Tair_SMHI), shape = 4) +
  geom_line(aes(x = Date, y = Abisko_Tsoil), color = "#D55E00", linetype = "dashed", size=1) + #      # Vermilion
  geom_point(aes(x = Date, y = Abisko_Tsoil), color = "#D55E00", shape = 6, size=1) + #              # Vermilion
  geom_line(aes(x = Date, y = Vassijaure_Tsoil), color = "#0072B2", linetype = "dotdash", size=1) + ## Blue
  geom_point(aes(x = Date, y = Vassijaure_Tsoil), color = "#0072B2", shape = 4, size=1) + # Blue
  #scale_color_viridis(discrete = TRUE, option = "E") +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = "Time of year", y = "Mean diel temperature °C", title = "Air and soil temperature") +
  #xlab("Time") + ylab("Temperature C") +
  theme_ipsum() +
  theme_bw()
  
  legend("bottomleft", 
         legend = c("Group 1", "Group 2"), 
         col = c(rgb(0.2,0.4,0.1,0.7), 
                 rgb(0.8,0.4,0.1,0.7)), 
         pch = c(17,19), 
         bty = "n", 
         pt.cex = 2, 
         cex = 1.2, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.1, 0.1))
#
avgT_long2 %>%
  dplyr::filter(Site == "Abisko_Tair" | Site == "Abisko_Tsoil") %>%
  ggplot() +
  geom_line(aes(x = Date, y = dielT, color = Site))


avgT_long %>%
  #pivot_longer(3:4, names_to = "Temp", values_to = "dielT") %>%
  #group_by(Site) %>%
  ggplot() +
  geom_line(aes(x = Date, y = dielT_air, colour = Site),linetype = "solid", size = 1) +
  geom_line(aes(x = Date, y = dielT_soil, colour = Site), linetype = "dotdash", size = 1) +
  #scale_color_viridis(discrete = TRUE, option = "E") +
  xlab("Time") + ylab("Temperature C") +
  theme_bw()
#facet_wrap( ~ Site)
#
# ### PAR ###
avgPAR_long %>%
  ggplot(aes(x = Date, y = diel_PAR, colour = Site)) +
  geom_point() +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  xlab("Time") + ylab("Temperature C")
#
#
#
#------- # Outliers # -------
#
hist(avgT_wide2$Abisko_Tsoil, main = "Histogram")
hist(avgT_wide2$Vassijaure_Tsoil, main = "Histogram")
# Cleveland plot
dotchart(avgT_wide2$Abisko_Tsoil, 
         main="Cleveland plot", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12),
         gpch = 12, gcolor = 1)
dotchart(avgT_wide2$Vassijaure_Tsoil, 
         main="Cleveland plot", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12),
         gpch = 12, gcolor = 1)
#
# Abisko soil temperature sensors
hist(Abisko_EM50$A1N_Tsoil, main = "Histogram - A1N_Tsoil")
dotchart(Abisko_EM50$A1N_Tsoil,
         main="Cleveland plot - A1N_Tsoil", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12),
         gpch = 12, gcolor = 1)
hist(Abisko_EM50$A2N_Tsoil, main = "Histogram - A2N_Tsoil")
hist(Abisko_EM50$A3N_Tsoil, main = "Histogram - A3N_Tsoil")
hist(Abisko_EM50$A4N_Tsoil, main = "Histogram - A4N_Tsoil")
hist(Abisko_EM50$A5N_Tsoil, main = "Histogram - A5N_Tsoil")
#
Abisko_EM50 %>% ggplot(aes(x = Date, y = A1N_Tsoil)) + geom_point()
Abisko_EM50 %>% ggplot(aes(x = Date, y = A2N_Tsoil)) + geom_point()
Abisko_EM50 %>% ggplot(aes(x = Date, y = A3N_Tsoil)) + geom_point()
Abisko_EM50 %>% ggplot(aes(x = Date, y = A4N_Tsoil)) + geom_point()
Abisko_EM50 %>% ggplot(aes(x = Date, y = A5N_Tsoil)) + geom_point()
Abisko_EM50 %>% ggplot() + 
  geom_point(aes(x = date(Date), y = A1N_Tsoil, shape = "A1N")) +
  geom_point(aes(x = date(Date), y = A2N_Tsoil, shape = "A2N")) +
  geom_point(aes(x = date(Date), y = A3N_Tsoil, shape = "A3N")) +
  geom_point(aes(x = date(Date), y = A4N_Tsoil, shape = "A4N")) +
  geom_point(aes(x = date(Date), y = A5N_Tsoil, shape = "A5N")) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  #coord_cartesian(xlim = c(ymd("2019-07-01"),ymd("2020-09-28"))) +
  labs(x = "Time of year", y = "Measured temperature °C", title = "Soil temperature - Abisko outliers") +
  theme_bw()
#
#
# Vassijaure soil temperature sensors
hist(Vassijaure_EM50$V1N_Tsoil, main = "Histogram - V1N_Tsoil")
dotchart(Vassijaure_EM50$V1N_Tsoil,
         main="Cleveland plot - V1N_Tsoil", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12),
         gpch = 12, gcolor = 1)
hist(Vassijaure_EM50$V2N_Tsoil, main = "Histogram - V2N_Tsoil")
hist(Vassijaure_EM50$V3N_Tsoil, main = "Histogram - V3N_Tsoil")
hist(Vassijaure_EM50$V4N_Tsoil, main = "Histogram - V4N_Tsoil")
hist(Vassijaure_EM50$V5N_Tsoil, main = "Histogram - V5N_Tsoil")
#
Vassijaure_EM50 %>% ggplot(aes(x = V_Date, y = V1N_Tsoil)) + geom_point()
Vassijaure_EM50 %>% ggplot(aes(x = V_Date, y = V2N_Tsoil)) + geom_point()
Vassijaure_EM50 %>% ggplot(aes(x = V_Date, y = V3N_Tsoil)) + geom_point()
Vassijaure_EM50 %>% ggplot(aes(x = V_Date, y = V4N_Tsoil)) + geom_point()
Vassijaure_EM50 %>% ggplot(aes(x = V_Date, y = V5N_Tsoil)) + geom_point()
Vassijaure_EM50 %>% ggplot() + 
  geom_point(aes(x = date(V_Date), y = V1N_Tsoil, shape = "V1N")) +
  geom_point(aes(x = date(V_Date), y = V2N_Tsoil, shape = "V2N")) +
  geom_point(aes(x = date(V_Date), y = V3N_Tsoil, shape = "V3N")) +
  geom_point(aes(x = date(V_Date), y = V4N_Tsoil, shape = "V4N")) +
  geom_point(aes(x = date(V_Date), y = V5N_Tsoil, shape = "V5N")) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  #coord_cartesian(xlim = c(ymd("2019-07-01"),ymd("2020-09-28"))) +
  labs(x = "Time of year", y = "Measured temperature °C", title = "Soil temperature - Vassijaure outliers") +
  theme_bw()
#
#
# Air temperature
# Abisko
hist(Abisko_Tair$Tair_A2, main = "Histogram - Tair Abisko")
dotchart(Abisko_Tair$Tair_A2,
         main="Cleveland plot - Tair Abisko", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12),
         gpch = 12, gcolor = 1)
Abisko_Tair %>% ggplot(aes(x = date(Date_A), y = Tair_A2)) + geom_point() +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  labs(x = "Time of year", y = "Measured temperature °C", title = "Air temperature - Abisko") +
  theme_bw()
#
# Vassijaure
Vassijaure_Tair %>% ggplot(aes(x = date(Date_V), y = Tair_V2)) + geom_point() +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  labs(x = "Time of year", y = "Measured temperature °C", title = "Air temperature - Vassijaure") +
  theme_bw()

#
# ### PAR ###
Abisko_EM50 %>% ggplot() + 
  geom_point(aes(x = date(Date), y = A1_PAR, shape = "A1N")) +
  geom_point(aes(x = date(Date), y = A2_PAR, shape = "A2N")) +
  geom_point(aes(x = date(Date), y = A3_PAR, shape = "A3N")) +
  geom_point(aes(x = date(Date), y = A4_PAR, shape = "A4N")) +
  geom_point(aes(x = date(Date), y = A5_PAR, shape = "A5N")) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  #coord_cartesian(xlim = c(ymd("2019-07-01"),ymd("2020-09-28"))) +
  labs(x = "Time of year", y = expression("Measured PAR µmol "*m^-2*s^-1), title = "PAR - Abisko outliers") +
  theme_bw()
#
Vassijaure_EM50 %>% ggplot() + 
  geom_point(aes(x = date(V_Date), y = V1_PAR, shape = "V1N")) +
  geom_point(aes(x = date(V_Date), y = V2_PAR, shape = "V2N")) +
  geom_point(aes(x = date(V_Date), y = V3_PAR, shape = "V3N")) +
  geom_point(aes(x = date(V_Date), y = V4_PAR, shape = "V4N")) +
  geom_point(aes(x = date(V_Date), y = V5_PAR, shape = "V5N")) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  #coord_cartesian(xlim = c(ymd("2019-07-01"),ymd("2020-09-28"))) +
  labs(x = "Time of year", y = expression("Measured PAR µmol "*m^-2*s^-1), title = "PAR - Vassijaure outliers") +
  theme_bw()



#------- # Leftovers # -------

data <- data.frame(
  time=seq(from=Sys.Date()-40, to=Sys.Date(), by=1 ), 
  value1=runif(41), 
  value2=runif(41)+0.7
)
# Then you can create the xts format:
don=xts( x=data[,-1], order.by=data$time)

# Chart
p <- dygraph(don)
p

# Interesting, but does not let me zoom down to negative values
don=xts( x=avgT_long[,-1], order.by=avgT_long$Date)
p <- dygraph(don)
p

# Not really working as of now
#
# Abisko_Tair <- Abisko_Tair %>%
#   mutate(avgTair_A = rowMeans(do.call(rbind, list(Abisko_Tair$Tair_A39_1, Abisko_Tair$Tair_C1, Abisko_Tair$Tair_31, Abisko_Tair$Tair_39_2)), na.rm = TRUE))
# mean(c(Abisko_Tair$Tair_A39_1, Abisko_Tair$Tair_C1, Abisko_Tair$Tair_31, Abisko_Tair$Tair_39_2), na.rm = TRUE)
#

avgT_wide2 %>% 
  ggplot(aes(x = Date, y = Vassijaure_Tair)) + geom_smooth()
#
#
# Layout matrix
library(grid)
library(lattice)
gs <- lapply(1:3, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
lay <- rbind(c(2,2,1),
             c(3,3,1))
grid.arrange(grobs = gs, layout_matrix = lay)
#
#

# Soil temperatures - all
soilT_plot2 <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  #annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -3, ymax = 0, fill = "white", alpha = 0.9) +
  #geom_hline(yintercept = -3, color = "#D55E00") +
  
  # Label arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[1]), y = -3, xend = as.Date(DayOfLH$Day_of_label[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = -3, xend = as.Date(DayOfLH$Day_of_label[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = -3, xend = as.Date(DayOfLH$Day_of_label[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[4]), y = -3, xend = as.Date(DayOfLH$Day_of_label[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[5]), y = -3, xend = as.Date(DayOfLH$Day_of_label[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[6]), y = 3, xend = as.Date(DayOfLH$Day_of_label[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[7]), y = 3, xend = as.Date(DayOfLH$Day_of_label[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[8]), y = 3, xend = as.Date(DayOfLH$Day_of_label[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[9]), y = 3, xend = as.Date(DayOfLH$Day_of_label[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[10]), y = 3, xend = as.Date(DayOfLH$Day_of_label[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[11]), y = 3, xend = as.Date(DayOfLH$Day_of_label[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[12]), y = 3, xend = as.Date(DayOfLH$Day_of_label[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[13]), y = -3, xend = as.Date(DayOfLH$Day_of_label[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[14]), y = -3, xend = as.Date(DayOfLH$Day_of_label[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[15]), y = -3, xend = as.Date(DayOfLH$Day_of_label[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  
  # Harvest arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[1]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[2]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[3]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[4]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[5]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[6]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[7]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[8]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[9]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[10]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[11]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[12]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[13]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[14]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[15]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  
  geom_line(aes(x = Date, y = Vassijaure_Tsoil, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) +
  geom_line(aes(x = Date, y = Abisko_Tsoil, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15))+#, minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%d-%b") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-08-26")), ylim = c(-5,15)) +
  labs(x = NULL, y = "Soil (°C)") + # x = "Time of year", , title = "Soil temperature"
  guides(lty = "none") + # guide_legend(title = "Soil temperature")
  theme_bw(base_size = 25) +
  theme(legend.position = "top")#, axis.title.y = element_text(size = 25))


#
# Air temperatures - all
airT_plot3 <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  #geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Katterjakk air temperature")) + 
  geom_line(aes(x = Date, y = Vassijaure_Tair, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) + 
  #geom_line(aes(x = Date, y = Abisko_Tair_SMHI, lty = "Abisko air temperature SMHI")) +
  geom_line(aes(x = Date, y = Abisko_Tair, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  
  # Label arrows
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[1]), y = -5, xend = as.Date(DayOfL$Day_of_label[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[2]), y = -5, xend = as.Date(DayOfL$Day_of_label[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[3]), y = -5, xend = as.Date(DayOfL$Day_of_label[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[4]), y = -5, xend = as.Date(DayOfL$Day_of_label[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[5]), y = -5, xend = as.Date(DayOfL$Day_of_label[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[6]), y = -5, xend = as.Date(DayOfL$Day_of_label[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[7]), y = -5, xend = as.Date(DayOfL$Day_of_label[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[8]), y = -5, xend = as.Date(DayOfL$Day_of_label[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[9]), y = -5, xend = as.Date(DayOfL$Day_of_label[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[10]), y = -5, xend = as.Date(DayOfL$Day_of_label[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[11]), y = -5, xend = as.Date(DayOfL$Day_of_label[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[12]), y = -5, xend = as.Date(DayOfL$Day_of_label[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[13]), y = -5, xend = as.Date(DayOfL$Day_of_label[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[14]), y = -5, xend = as.Date(DayOfL$Day_of_label[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[15]), y = -5, xend = as.Date(DayOfL$Day_of_label[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[16]), y = -5, xend = as.Date(DayOfL$Day_of_label[16]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[17]), y = -5, xend = as.Date(DayOfL$Day_of_label[17]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[18]), y = -5, xend = as.Date(DayOfL$Day_of_label[18]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[19]), y = -5, xend = as.Date(DayOfL$Day_of_label[19]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[20]), y = -5, xend = as.Date(DayOfL$Day_of_label[20]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[21]), y = -5, xend = as.Date(DayOfL$Day_of_label[21]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[22]), y = -5, xend = as.Date(DayOfL$Day_of_label[22]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[23]), y = -5, xend = as.Date(DayOfL$Day_of_label[23]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[24]), y = -5, xend = as.Date(DayOfL$Day_of_label[24]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[25]), y = -5, xend = as.Date(DayOfL$Day_of_label[25]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[26]), y = -5, xend = as.Date(DayOfL$Day_of_label[26]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[27]), y = -5, xend = as.Date(DayOfL$Day_of_label[27]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[28]), y = -5, xend = as.Date(DayOfL$Day_of_label[28]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[29]), y = -5, xend = as.Date(DayOfL$Day_of_label[29]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfL$Day_of_label[30]), y = -5, xend = as.Date(DayOfL$Day_of_label[30]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
  
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[1]), y = -5, xend = as.Date(DayOf$Day_of_harvest[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[2]), y = -5, xend = as.Date(DayOf$Day_of_harvest[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[3]), y = -5, xend = as.Date(DayOf$Day_of_harvest[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[4]), y = -5, xend = as.Date(DayOf$Day_of_harvest[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[5]), y = -5, xend = as.Date(DayOf$Day_of_harvest[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[6]), y = -5, xend = as.Date(DayOf$Day_of_harvest[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[7]), y = -5, xend = as.Date(DayOf$Day_of_harvest[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[8]), y = -5, xend = as.Date(DayOf$Day_of_harvest[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[9]), y = -5, xend = as.Date(DayOf$Day_of_harvest[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[10]), y = -5, xend = as.Date(DayOf$Day_of_harvest[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[11]), y = -5, xend = as.Date(DayOf$Day_of_harvest[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[12]), y = -5, xend = as.Date(DayOf$Day_of_harvest[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[13]), y = -5, xend = as.Date(DayOf$Day_of_harvest[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[14]), y = -5, xend = as.Date(DayOf$Day_of_harvest[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[15]), y = -5, xend = as.Date(DayOf$Day_of_harvest[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[16]), y = -5, xend = as.Date(DayOf$Day_of_harvest[16]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[17]), y = -5, xend = as.Date(DayOf$Day_of_harvest[17]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[18]), y = -5, xend = as.Date(DayOf$Day_of_harvest[18]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[19]), y = -5, xend = as.Date(DayOf$Day_of_harvest[19]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[20]), y = -5, xend = as.Date(DayOf$Day_of_harvest[20]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[21]), y = -5, xend = as.Date(DayOf$Day_of_harvest[21]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[22]), y = -5, xend = as.Date(DayOf$Day_of_harvest[22]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[23]), y = -5, xend = as.Date(DayOf$Day_of_harvest[23]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[24]), y = -5, xend = as.Date(DayOf$Day_of_harvest[24]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[25]), y = -5, xend = as.Date(DayOf$Day_of_harvest[25]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[26]), y = -5, xend = as.Date(DayOf$Day_of_harvest[26]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[27]), y = -5, xend = as.Date(DayOf$Day_of_harvest[27]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[28]), y = -5, xend = as.Date(DayOf$Day_of_harvest[28]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[29]), y = -5, xend = as.Date(DayOf$Day_of_harvest[29]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[30]), y = -5, xend = as.Date(DayOf$Day_of_harvest[30]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
  
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%d-%b") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-08-26"))) +
  labs(x = NULL, y = "Air (°C)") + # x = "Time of year",  , title = "Air temperature" 
  guides(lty = guide_legend(title = "Mean diel temperature"))+ #lty = guide_legend(title = "Mean diel temperature")) +
  theme_bw(base_size = 25) +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.text.y = element_text(size = 15))
#
airT_legend3 <- get_legend(airT_plot3)
airT_plot3.2 <- airT_plot3 + theme_bw(base_size = 25) + theme(legend.position = "none", axis.text.y = element_text(size = 18))#, axis.title.y = element_text(size = 25)) 
#airT_plot <- airT_plot + guides(lty = NULL)


#
# Final graph with arrows:
# Fig. 1 with arrows: Air temperature, Soil temperature, and Snow depth with soil moisture as VWC and converted GWC
grid.arrange(airT_legend, airT_plot.2, soilT_plot2, widths = c(2.8), heights = c(0.5, 3, 3))
#
# Just Soil Temperature
grid.arrange(airT_legend, soilT_plot2, widths = c(2.8), heights = c(0.5, 3))
# just Air Temperature
grid.arrange(airT_legend3, airT_plot3.2, widths = c(2.8), heights = c(0.5, 3))










soilT_plot3 <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + # Abisko snow
  geom_hline(yintercept = 0, color = "#999999", linewidth = 1) +
  #annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -3, ymax = 0, fill = "white", alpha = 0.9) +
  #geom_hline(yintercept = -3, color = "#D55E00") +
  
  # Arrows from label to harvest
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[1]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[1]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[4]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[4]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[5]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[5]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[6]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[6]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[7]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[7]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[8]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[8]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[9]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[9]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[10]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[10]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[11]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[11]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[12]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[12]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[13]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[13]), yend = 0), linewidth = 0.9, curvature = 0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[14]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[14]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  geom_curve(aes(x = as.Date(DayOfLH$Day_of_label[15]), y = 0, xend = as.Date(DayOfLH$Day_of_harvest[15]), yend = 0), linewidth = 0.9, curvature = -0.5, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "#CC79A7") +
  
  # Label arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[1]), y = -3, xend = as.Date(DayOfLH$Day_of_label[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[2]), y = -3, xend = as.Date(DayOfLH$Day_of_label[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[3]), y = -3, xend = as.Date(DayOfLH$Day_of_label[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[4]), y = -3, xend = as.Date(DayOfLH$Day_of_label[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[5]), y = -3, xend = as.Date(DayOfLH$Day_of_label[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[6]), y = 3, xend = as.Date(DayOfLH$Day_of_label[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[7]), y = 3, xend = as.Date(DayOfLH$Day_of_label[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[8]), y = 3, xend = as.Date(DayOfLH$Day_of_label[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[9]), y = 3, xend = as.Date(DayOfLH$Day_of_label[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[10]), y = 3, xend = as.Date(DayOfLH$Day_of_label[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[11]), y = 3, xend = as.Date(DayOfLH$Day_of_label[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[12]), y = 3, xend = as.Date(DayOfLH$Day_of_label[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[13]), y = -3, xend = as.Date(DayOfLH$Day_of_label[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[14]), y = -3, xend = as.Date(DayOfLH$Day_of_label[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_label[15]), y = -3, xend = as.Date(DayOfLH$Day_of_label[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#56B4E9") +
  
  # Harvest arrows
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[1]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[2]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[3]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[4]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[5]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[6]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[7]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[8]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[9]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[10]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[11]), y = 3, xend = as.Date(DayOfLH$Day_of_harvest[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[12]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[13]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[14]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  geom_segment(aes(x = as.Date(DayOfLH$Day_of_harvest[15]), y = -3, xend = as.Date(DayOfLH$Day_of_harvest[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.3, "cm")), color = "#D55E00") +
  
  geom_line(aes(x = Date, y = Vassijaure_Tsoil, lty = "Vassijaure"), na.rm = TRUE, linewidth = 1) +
  geom_line(aes(x = Date, y = Abisko_Tsoil, lty = "Abisko"), na.rm = TRUE, linewidth = 1) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15))+#, minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day", date_labels = "%d-%b") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-08-26")), ylim = c(-5,15)) +
  labs(x = NULL, y = "Soil (°C)") + # x = "Time of year", , title = "Soil temperature"
  guides(lty = "none") + # guide_legend(title = "Soil temperature")
  theme_bw(base_size = 25) +
  theme(legend.position = "top")#, axis.title.y = element_text(size = 25))

#
grid.arrange(airT_legend, soilT_plot3, widths = c(2.8), heights = c(0.5, 3))




# # Many arrows. Each unique for each date
# geom_segment(aes(x = as.Date(DayOfL$Day_of_label[1]), y = -5, xend = as.Date(DayOfL$Day_of_label[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[2]), y = -5, xend = as.Date(DayOfL$Day_of_label[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[3]), y = -5, xend = as.Date(DayOfL$Day_of_label[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[4]), y = -5, xend = as.Date(DayOfL$Day_of_label[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[5]), y = -5, xend = as.Date(DayOfL$Day_of_label[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[6]), y = -5, xend = as.Date(DayOfL$Day_of_label[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[7]), y = -5, xend = as.Date(DayOfL$Day_of_label[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[8]), y = -5, xend = as.Date(DayOfL$Day_of_label[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[9]), y = -5, xend = as.Date(DayOfL$Day_of_label[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[10]), y = -5, xend = as.Date(DayOfL$Day_of_label[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[11]), y = -5, xend = as.Date(DayOfL$Day_of_label[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[12]), y = -5, xend = as.Date(DayOfL$Day_of_label[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[13]), y = -5, xend = as.Date(DayOfL$Day_of_label[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[14]), y = -5, xend = as.Date(DayOfL$Day_of_label[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[15]), y = -5, xend = as.Date(DayOfL$Day_of_label[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[16]), y = -5, xend = as.Date(DayOfL$Day_of_label[16]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[17]), y = -5, xend = as.Date(DayOfL$Day_of_label[17]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[18]), y = -5, xend = as.Date(DayOfL$Day_of_label[18]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[19]), y = -5, xend = as.Date(DayOfL$Day_of_label[19]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[20]), y = -5, xend = as.Date(DayOfL$Day_of_label[20]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[21]), y = -5, xend = as.Date(DayOfL$Day_of_label[21]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[22]), y = -5, xend = as.Date(DayOfL$Day_of_label[22]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[23]), y = -5, xend = as.Date(DayOfL$Day_of_label[23]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[24]), y = -5, xend = as.Date(DayOfL$Day_of_label[24]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[25]), y = -5, xend = as.Date(DayOfL$Day_of_label[25]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[26]), y = -5, xend = as.Date(DayOfL$Day_of_label[26]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[27]), y = -5, xend = as.Date(DayOfL$Day_of_label[27]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[28]), y = -5, xend = as.Date(DayOfL$Day_of_label[28]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[29]), y = -5, xend = as.Date(DayOfL$Day_of_label[29]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   geom_segment(aes(x = as.Date(DayOfL$Day_of_label[30]), y = -5, xend = as.Date(DayOfL$Day_of_label[30]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#56B4E9") +
#   
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[1]), y = -5, xend = as.Date(DayOf$Day_of_harvest[1]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[2]), y = -5, xend = as.Date(DayOf$Day_of_harvest[2]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[3]), y = -5, xend = as.Date(DayOf$Day_of_harvest[3]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[4]), y = -5, xend = as.Date(DayOf$Day_of_harvest[4]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[5]), y = -5, xend = as.Date(DayOf$Day_of_harvest[5]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[6]), y = -5, xend = as.Date(DayOf$Day_of_harvest[6]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[7]), y = -5, xend = as.Date(DayOf$Day_of_harvest[7]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[8]), y = -5, xend = as.Date(DayOf$Day_of_harvest[8]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[9]), y = -5, xend = as.Date(DayOf$Day_of_harvest[9]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[10]), y = -5, xend = as.Date(DayOf$Day_of_harvest[10]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[11]), y = -5, xend = as.Date(DayOf$Day_of_harvest[11]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[12]), y = -5, xend = as.Date(DayOf$Day_of_harvest[12]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[13]), y = -5, xend = as.Date(DayOf$Day_of_harvest[13]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[14]), y = -5, xend = as.Date(DayOf$Day_of_harvest[14]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[15]), y = -5, xend = as.Date(DayOf$Day_of_harvest[15]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[16]), y = -5, xend = as.Date(DayOf$Day_of_harvest[16]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[17]), y = -5, xend = as.Date(DayOf$Day_of_harvest[17]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[18]), y = -5, xend = as.Date(DayOf$Day_of_harvest[18]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[19]), y = -5, xend = as.Date(DayOf$Day_of_harvest[19]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[20]), y = -5, xend = as.Date(DayOf$Day_of_harvest[20]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[21]), y = -5, xend = as.Date(DayOf$Day_of_harvest[21]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[22]), y = -5, xend = as.Date(DayOf$Day_of_harvest[22]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[23]), y = -5, xend = as.Date(DayOf$Day_of_harvest[23]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[24]), y = -5, xend = as.Date(DayOf$Day_of_harvest[24]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[25]), y = -5, xend = as.Date(DayOf$Day_of_harvest[25]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[26]), y = -5, xend = as.Date(DayOf$Day_of_harvest[26]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[27]), y = -5, xend = as.Date(DayOf$Day_of_harvest[27]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[28]), y = -5, xend = as.Date(DayOf$Day_of_harvest[28]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[29]), y = -5, xend = as.Date(DayOf$Day_of_harvest[29]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +
#   geom_segment(aes(x = as.Date(DayOf$Day_of_harvest[30]), y = -5, xend = as.Date(DayOf$Day_of_harvest[30]), yend = 0), linewidth = 2, arrow = arrow(length = unit(0.5, "cm")), color = "#D55E00") +



#
#
#------- The End -------