# Load weather data
# From own data
#
library(tidyverse)
library(lubridate)
library(readxl)
library(viridis)
#
#
# Load the data
# Air temperature logger
Abisko_Tair <- read_xlsx("raw_data/Abisko_Tair.xlsx", skip = 5, col_names = TRUE, col_types = c("date", "date", "text","date","text","text","text","date","text","text","text","date","text","text","text","text", "text"))
Vassijaure_Tair <- read_xlsx("raw_data/Vassijaure_Tair.xlsx", skip = 5, col_names = TRUE, col_types = c("date", "date", "text","date","text","text","text","date","text","text","text","date","text","text","text","text", "text"))
#
# Soil temperature and moisture data and PAR sensor
Abisko_EM50 <- read_xlsx("raw_data/allEMdata Abisko.xlsx", col_names = TRUE)
Vassijaure_EM50 <- read_xlsx("raw_data/allEMdata Vassijaure.xlsx", col_names = TRUE)
#
#
#
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
# Mean diel temperature
Abisko_avgTair <- Abisko_Tair %>%
  group_by(date(Date_A)) %>%
  summarise(Abisko = mean(Tair_A2, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Date_A)")
#
Vassijuare_avgTair <- Vassijaure_Tair %>%
  group_by(date(Date_V)) %>%
  summarise(Vassijaure = mean(Tair_V2, na.rm = TRUE), .groups = "keep") %>%
  rename("Date" = "date(Date_V)")
#
avgTair <- left_join(Abisko_avgTair, Vassijuare_avgTair) %>%
  pivot_longer(2:3, names_to = "Site", values_to = "dielT")
# Plot it
avgTair %>%
  ggplot(aes(x = Date, y = dielT, colour = Site)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  xlab("Time") + ylab("Temperature C")





# Not really working as of now
#
# Abisko_Tair <- Abisko_Tair %>%
#   mutate(avgTair_A = rowMeans(do.call(rbind, list(Abisko_Tair$Tair_A39_1, Abisko_Tair$Tair_C1, Abisko_Tair$Tair_31, Abisko_Tair$Tair_39_2)), na.rm = TRUE))
# mean(c(Abisko_Tair$Tair_A39_1, Abisko_Tair$Tair_C1, Abisko_Tair$Tair_31, Abisko_Tair$Tair_39_2), na.rm = TRUE)
