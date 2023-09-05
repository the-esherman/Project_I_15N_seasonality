# Load weather data
# From own data
#
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
         A2N_Tsoil = replace(A1N_Tsoil, A2N_Tsoil < -15, NA),
         A3N_Tsoil = replace(A1N_Tsoil, A3N_Tsoil < -15, NA), # Has a single measure very low measure in March 2020
         A4N_Tsoil = replace(A1N_Tsoil, A4N_Tsoil < -15, NA),
         A5N_Tsoil = replace(A1N_Tsoil, A5N_Tsoil < -15, NA)) %>%
  mutate(A1N_Tsoil = if_else(month(Date) == 7 & A1N_Tsoil < 5, NA, A1N_Tsoil),
         A2N_Tsoil = if_else(month(Date) == 7 & A2N_Tsoil < 0, NA, A1N_Tsoil),
         A3N_Tsoil = if_else(month(Date) == 7 & A3N_Tsoil < 0, NA, A1N_Tsoil),
         A4N_Tsoil = if_else(month(Date) == 7 & A4N_Tsoil < 0, NA, A1N_Tsoil),
         A5N_Tsoil = if_else(month(Date) == 7 & A5N_Tsoil < 0, NA, A1N_Tsoil)) %>%
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
#
#
#------- # Soil Moisture # -------

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

  
  mutate(A1_PAR = if_else(month(Date) == 7 & A1_PAR < 5, NA, A1_PAR),
         A2_PAR = if_else(month(Date) == 7 & A2_PAR < 0, NA, A1_PAR),
         A3_PAR = if_else(month(Date) == 7 & A3_PAR < 0, NA, A1_PAR),
         A4_PAR = if_else(month(Date) == 7 & A4_PAR < 0, NA, A1_PAR),
         A5_PAR = if_else(month(Date) == 7 & A5_PAR < 0, NA, A1_PAR)) #%>%
  
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
ggplot() + geom_rect(winterP, aes(x=, ymin=-Inf, ymax=Inf), stat="identity", alpha = 0.5, fill = 'grey', inherit.aes = FALSE)
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

# Air temperatures - all
airT_plot <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Abisko snow
  #geom_line(aes(x = Date, y = Katterjakk_Tair_SMHI, lty = "Katterjakk air temperature")) + 
  geom_line(aes(x = Date, y = Vassijaure_Tair, lty = "Vassijaure"), na.rm = TRUE) + 
  #geom_line(aes(x = Date, y = Abisko_Tair_SMHI, lty = "Abisko air temperature SMHI")) +
  geom_line(aes(x = Date, y = Abisko_Tair, lty = "Abisko"), na.rm = TRUE) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16"))) +
  labs(x = NULL, y = "Air temperature (°C)") + # x = "Time of year",  , title = "Air temperature" 
  guides(lty = guide_legend(title = "Mean diel temperature")) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top", axis.text.x = element_blank())
#
# Soil temperatures - all
soilT_plot <- avgT_wide2 %>% ggplot() +
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Vassijaure snow
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) + # Abisko snow
  geom_line(aes(x = Date, y = Vassijaure_Tsoil, lty = "Vassijaure"), na.rm = TRUE) +
  geom_line(aes(x = Date, y = Abisko_Tsoil, lty = "Abisko"), na.rm = TRUE) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20), minor_breaks = c(-15, -5, 5, 15)) +
  scale_x_date(date_breaks = "30 day", date_minor_breaks = "5 day") +
  coord_cartesian(xlim = c(as.Date("2019-08-06"),as.Date("2020-09-16")), ylim = c(-10,20)) +
  labs(x = "Time of year", y = "Soil temperature (°C)") + # , title = "Soil temperature"
  guides(lty = "none") + # guide_legend(title = "Soil temperature")
  theme_bw(base_size = 15) +
  theme(legend.position = "top")
#
grid.arrange(airT_plot, soilT_plot, ncol = 1)


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
#
#------- The End -------