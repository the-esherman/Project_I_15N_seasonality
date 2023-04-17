# Load weather data
# From own data
#
library(tidyverse)
library(lubridate)
library(readxl)
library(viridis)
library(dygraphs)
library(xts)
library(plotly)
library(hrbrthemes)
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

#
#
#------- ### Clean data ### -------
#------- # Air Temperature # -------|
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
#------- # EM50 loggers # -------|
# Remove last lines without a date and extra header part. Set NaN to NA
Abisko_EM50 <- Abisko_EM50 %>%
  slice(5:n()) %>%
  select(-c(NA...9, NA...18, NA...27, NA...36, NA...45)) %>%
  filter(!(is.na(Date))) %>%
  mutate(across(where(is.character), ~na_if(.,"NaN")),
         across(where(is.numeric), ~na_if(.,NaN))) %>%
  mutate(across(c(Date, A1_Date, A2_Date, A3_Date, A4_Date, A5_Date), ymd_hms))
#
Vassijaure_EM50 <- Vassijaure_EM50 %>%
slice(5:(n()-24)) %>% # No data for the last 24 logs
  select(-c(NA...18, NA...27, NA...36, NA...37, NA...43, NA...44)) %>%
  filter(!(is.na(V_Date))) %>% # if the date is not there, we cannot use the data
  mutate(across(where(is.character), ~na_if(.,"NaN")),
         across(where(is.numeric), ~na_if(.,NaN))) %>%
  mutate(across(c(V_Date, V1_Date, V2_Date, V3_Date, V4_Date, V5_Date), ymd_hms))
#
# Mean Soil temperature
Abisko_avgTsoil <- Abisko_EM50 %>%
  select(Date, A1N_Tsoil, A2N_Tsoil, A3N_Tsoil, A4N_Tsoil, A5N_Tsoil) %>%
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
  select(c(Date, Abisko, Vassijaure)) %>%
  pivot_longer(2:3, names_to = "Site", values_to = "dielT_soil")
#
#
# Combine all temperatures
avgT_wide <- left_join(avgTair_wide, avgTsoil_wide)
avgT_long <- left_join(avgTair_long, avgTsoil_long)
#
#------- # Soil Temperature # -------|


#
#
#
#------- # Soil Moisture # -------|

#
#
#
#------- # PAR # -------|

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

avgT_wide %>%
  ggplot() +
  geom_line(aes(x = Date, y = Abisko_Tair), color = "#E69F00", linetype = "solid", size=1, alpha = 0.4) + # Orange
  geom_point(aes(x = Date, y = Abisko_Tair), shape = 6) +
  geom_line(aes(x = Date, y = Vassijaure_Tair), linetype = "solid", size=1, alpha = 0.4) + # Bluish green, color = "#009E73"
  geom_point(aes(x = Date, y = Vassijaure_Tair), shape = 4) +
  geom_line(aes(x = Date, y = Abisko_Tsoil), linetype = "dashed", size=1) + # Vermilion, color = "#D55E00"
  geom_line(aes(x = Date, y = Vassijaure_Tsoil), linetype = "dotdash", size=1) + # Blue, color = "#0072B2"
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
         inset = c(0.1, 0.1)) +
  

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
#
#
#------- The End -------