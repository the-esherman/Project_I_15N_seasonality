# Load weather data
# From ANS station. Find annual temperature for the period 2010-2019
#
library(tidyverse)
library(lubridate)
#
# Data from SMHI
T_ANS <- read.csv2("raw_data/smhi_Abisko_Temp_old.csv", skip = 9) # Temperature
Precip_ANS <- read.csv2("raw_data/smhi_Abisko_precip.csv", skip = 9) # Precipitation
Snow_ANS <- read.csv2("raw_data/smhi_Abisko_snow.csv", skip = 9) # Snow depth
#
T_ANS <- T_ANS %>%
  select(1:3) %>%
  rename("Tid_UTC" = "Tid..UTC.") %>%
  mutate(across(c("Lufttemperatur"), as.numeric)) %>%
  mutate(across("Datum", ymd)) %>%
  mutate(across("Tid_UTC", hms))
#
T_ANS %>%
  filter(across(year(Datum) == 2019))
  group_by(across)
