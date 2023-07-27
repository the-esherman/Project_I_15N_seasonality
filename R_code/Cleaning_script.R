# Cleaning script
#
# - 15N vegetation and root data
# By Emil A.S. Andersen
# 
#------- ### Libraries ### -------
library(ggpubr)
library(tidyverse)
library(readxl)
library(gridExtra)
#
#
#
#------- ### Load data ### -------
#
DataName <- "raw_data/15N vegetation and roots v0.38.xlsx"
#
# Core data; Snow, soil mass, DW/FW of soil
core15N <- read_xlsx(DataName, sheet = "Core", skip = 1, col_names = TRUE)
#
# Biomass, d15N, atom% 15N, and recovery ("R_") (discard later)
vegroot15N <- read_xlsx(DataName, sheet = "15N", skip = 1, col_names = TRUE)
# Natural abundance
vegrootsNatAbu <- read_xlsx(DataName, sheet = "NatAbu", col_names = TRUE)
#
# Microbial biomass, d15N, atom% 15N, and recovery ("R_") (discard later)
#Mic15N <- read_xlsx(DataName, sheet = "MBN", skip = 1, col_names = TRUE)
Mic15N <- read_csv("clean_data/Mic_15N_data.csv", skip = 1, col_names = TRUE)
#
#vegroot15Nlong <- read_xlsx(DataName, sheet = "Long", col_names = TRUE)
#
IRMS <- read_xlsx("raw_data/IRMS_data v0.9.xlsx", col_names = TRUE)
#
Nconc <- read_xlsx(DataName, sheet = "Nconc", skip = 1, col_names = TRUE)
#
# To format inorganic N
#inorgN <- read_xlsx("raw_data/Inorganic N v1.xlsx", sheet = "InorgN", col_names = TRUE, col_types = c("text", "text", "text","text", "numeric", "numeric"))
#soilExtr <- read_xlsx("raw_data/Inorganic N v1.xlsx", sheet = "Soil_extr", col_names = TRUE, col_types = c("text", "text", "text","text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
#
# Table data for N concentrations and snow
#table_dat_N_atom <- read_xlsx(DataName, sheet = "Table", skip = 1, col_names = TRUE)
#
#
# Inorganic N
inorgN <- read_xlsx("raw_data/Inorganic N v1.9.xlsx", sheet = "Soil_extr", col_names = TRUE, col_types = c("text", "text", "text","text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
Blanks <- read_xlsx("raw_data/Inorganic N v1.9.xlsx", sheet = "Blanks", col_names = TRUE, col_types = c("text", "text", "text","text", "numeric", "numeric"))
#
# Reference atmospheric Nitrogen
# Either 0.003676 or 1/272 (more decimals)
Nair_Rst = 1/272
#
# Added 15N; mg 15N pr patch
N_add <- 1.084
#
#
# Plant 15N data ----
# Transform to long format - 15N enriched data and natural abundance
#
# Biomass
vegroot15N_bioLong <- vegroot15N %>%
  dplyr::select(Site, Plot, MP, Round, 23:37) %>%
  dplyr::rename("ES_S" = "ESS",
                "DS_S" = "DSS",
                "G_S" = "GS",
                "O_S" = "OS",
                "U_S" = "US",
                "ES_CR" = "ESCR",
                "DS_CR" = "DSCR",
                "G_CR" = "GCR",
                "ES_FR" = "ESFR",
                "DS_FR" = "DSFR",
                "G_FR" = "GFR",
                "O_FR" = "OFR",
                "Root_CR" = "CR",
                "Root_FR" = "FR",
                "RootG_FR" = "RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "Biomass") %>%
  mutate(Biomass = na_if(Biomass, 0))
#
# d15N values
vegroot15N_dLong <- vegroot15N %>%
  dplyr::select(Site, Plot, MP, Round, 41:55) %>%
  dplyr::rename("ES_S" = "d15N_ESS",
                "DS_S" = "d15N_DSS",
                "G_S" = "d15N_GS",
                "O_S" = "d15N_OS",
                "U_S" = "d15N_US",
                "ES_CR" = "d15N_ESCR",
                "DS_CR" = "d15N_DSCR",
                "G_CR" = "d15N_GCR",
                "ES_FR" = "d15N_ESFR",
                "DS_FR" = "d15N_DSFR",
                "G_FR" = "d15N_GFR",
                "O_FR" = "d15N_OFR",
                "Root_CR" = "d15N_CR",
                "Root_FR" = "d15N_FR",
                "RootG_FR" = "d15N_RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "d15N")
#
# N concentration
vegroot15N_NLong <- vegroot15N %>%
  dplyr::select(Site, Plot, MP, Round, 86:100) %>%
  dplyr::rename("ES_S" = "Nconc_ESS",
                "DS_S" = "Nconc_DSS",
                "G_S" = "Nconc_GS",
                "O_S" = "Nconc_OS",
                "U_S" = "Nconc_US",
                "ES_CR" = "Nconc_ESCR",
                "DS_CR" = "Nconc_DSCR",
                "G_CR" = "Nconc_GCR",
                "ES_FR" = "Nconc_ESFR",
                "DS_FR" = "Nconc_DSFR",
                "G_FR" = "Nconc_GFR",
                "O_FR" = "Nconc_OFR",
                "Root_CR" = "Nconc_CR",
                "Root_FR" = "Nconc_FR",
                "RootG_FR" = "Nconc_RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "Nconc")
#
# atom% 15N
vegroot15N_atomLong <- vegroot15N %>%
  dplyr::select(Site, Plot, MP, Round, 71:85) %>%
  dplyr::rename("ES_S" = "atom_ESS",
                "DS_S" = "atom_DSS",
                "G_S" = "atom_GS",
                "O_S" = "atom_OS",
                "U_S" = "atom_US",
                "ES_CR" = "atom_ESCR",
                "DS_CR" = "atom_DSCR",
                "G_CR" = "atom_GCR",
                "ES_FR" = "atom_ESFR",
                "DS_FR" = "atom_DSFR",
                "G_FR" = "atom_GFR",
                "O_FR" = "atom_OFR",
                "Root_CR" = "atom_CR",
                "Root_FR" = "atom_FR",
                "RootG_FR" = "atom_RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "atom_pc")
#
# Recovery (Replace later)
vegroot15N_RLong <- vegroot15N %>%
  dplyr::select(Site, Plot, MP, Round, 56:70) %>%
  dplyr::rename("ES_S" = "R_ESS",
                "DS_S" = "R_DSS",
                "G_S" = "R_GS",
                "O_S" = "R_OS",
                "U_S" = "R_US",
                "ES_CR" = "R_ESCR",
                "DS_CR" = "R_DSCR",
                "G_CR" = "R_GCR",
                "ES_FR" = "R_ESFR",
                "DS_FR" = "R_DSFR",
                "G_FR" = "R_GFR",
                "O_FR" = "R_OFR",
                "Root_CR" = "R_CR",
                "Root_FR" = "R_FR",
                "RootG_FR" = "R_RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "Recovery")
#
#
# Combine to one file
vegroot15N_long <- vegroot15N_bioLong %>%
  left_join(vegroot15N_dLong, by = join_by(Site, Plot, MP, Round, Type)) %>%
  left_join(vegroot15N_atomLong, by = join_by(Site, Plot, MP, Round, Type)) %>%
  left_join(vegroot15N_NLong, by = join_by(Site, Plot, MP, Round, Type))
vegroot15N_long <- vegroot15N_long %>%
  add_column(Organ = str_split_fixed(vegroot15N_long$Type,"\\w+_",n=2)[,2], .after = "Type") %>%
  add_column(Species = str_split_fixed(vegroot15N_long$Type,"_\\w+",n=2)[,1], .after = "Type")
#
#
#
# Natural abundance atom% 15N
vegrootsNatAbu_atomLong <- vegrootsNatAbu %>%
  dplyr::select(Site, Plot, MP, Round, 20:34) %>%
  dplyr::rename("ES_S" = "atom_NatAb_ESS",
                "DS_S" = "atom_NatAb_DSS",
                "G_S" = "atom_NatAb_GS",
                "O_S" = "atom_NatAb_OS",
                "U_S" = "atom_NatAb_US",
                "ES_CR" = "atom_NatAb_ESCR",
                "DS_CR" = "atom_NatAb_DSCR",
                "G_CR" = "atom_NatAb_GCR",
                "ES_FR" = "atom_NatAb_ESFR",
                "DS_FR" = "atom_NatAb_DSFR",
                "G_FR" = "atom_NatAb_GFR",
                "O_FR" = "atom_NatAb_OFR",
                "Root_CR" = "atom_NatAb_CR",
                "Root_FR" = "atom_NatAb_FR",
                "RootG_FR" = "atom_NatAb_RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "atom_pc_NatAb")
#
# Natural abundance d15N
vegrootsNatAbu_dLong <- vegrootsNatAbu %>%
  dplyr::select(Site, Plot, MP, Round, 5:19) %>%
  dplyr::rename("ES_S" = "d15N_NatAb_ESS",
                "DS_S" = "d15N_NatAb_DSS",
                "G_S" = "d15N_NatAb_GS",
                "O_S" = "d15N_NatAb_OS",
                "U_S" = "d15N_NatAb_US",
                "ES_CR" = "d15N_NatAb_ESCR",
                "DS_CR" = "d15N_NatAb_DSCR",
                "G_CR" = "d15N_NatAb_GCR",
                "ES_FR" = "d15N_NatAb_ESFR",
                "DS_FR" = "d15N_NatAb_DSFR",
                "G_FR" = "d15N_NatAb_GFR",
                "O_FR" = "d15N_NatAb_OFR",
                "Root_CR" = "d15N_NatAb_CR",
                "Root_FR" = "d15N_NatAb_FR",
                "RootG_FR" = "d15N_NatAb_RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "d15N_NatAb")
#
# Natural abundance N concentration
vegrootsNatAbu_NLong <- vegrootsNatAbu %>%
  dplyr::select(Site, Plot, MP, Round, 35:49) %>%
  dplyr::rename("ES_S" = "Nconc_NatAb_ESS",
                "DS_S" = "Nconc_NatAb_DSS",
                "G_S" = "Nconc_NatAb_GS",
                "O_S" = "Nconc_NatAb_OS",
                "U_S" = "Nconc_NatAb_US",
                "ES_CR" = "Nconc_NatAb_ESCR",
                "DS_CR" = "Nconc_NatAb_DSCR",
                "G_CR" = "Nconc_NatAb_GCR",
                "ES_FR" = "Nconc_NatAb_ESFR",
                "DS_FR" = "Nconc_NatAb_DSFR",
                "G_FR" = "Nconc_NatAb_GFR",
                "O_FR" = "Nconc_NatAb_OFR",
                "Root_CR" = "Nconc_NatAb_CR",
                "Root_FR" = "Nconc_NatAb_FR",
                "RootG_FR" = "Nconc_NatAb_RG"
  ) %>%
  pivot_longer(cols = 5:19, names_to = "Type", values_to = "Nconc_NatAb")
#
# Combine natural abundance
vegrootsNatAbu_long <- vegrootsNatAbu_dLong %>%
  left_join(vegrootsNatAbu_atomLong, by = join_by(Site, Plot, MP, Round, Type)) %>%
  left_join(vegrootsNatAbu_NLong, by = join_by(Site, Plot, MP, Round, Type))
vegrootsNatAbu_long <- vegrootsNatAbu_long %>%
  add_column(Organ = str_split_fixed(vegrootsNatAbu_long$Type,"\\w+_",n=2)[,2], .after = "Type") %>%
  add_column(Species = str_split_fixed(vegrootsNatAbu_long$Type,"_\\w+",n=2)[,1], .after = "Type")
#
#
#
# Combine both 15N enriched and natural abundance data
vegroot_long <- left_join(vegroot15N_long, vegrootsNatAbu_long, by = join_by(Site, Plot, MP, Round, Type, Species, Organ))
#
#
#
# Calculate recovery and add as column the other recovered R check that they match
vegroot_long <- vegroot_long %>%
  mutate(Recovery = ((atom_pc - atom_pc_NatAb)/100 * Nconc/100 * Biomass)/(N_add/1000) * 100) %>%
  mutate(Recovery = if_else(Recovery < 0, 0, Recovery)) %>%
  left_join(vegroot15N_RLong, by = join_by(Site, Plot, MP, Round, Type))
ggplot(vegroot_long, aes(Recovery.x, Recovery.y)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(label.y = 10, aes(label = after_stat(eq.label))) + stat_regline_equation(label.y = 9, aes(label = after_stat(rr.label))) # Perfect fit
#
vegroot_long <- vegroot_long %>%
  select(-(Recovery.y)) %>%
  rename(Recovery = Recovery.x)
#
#
#
# Save 15N enriched data and Natural abundance to a data file
write_csv(vegroot_long, "clean_data/Plant_15N_data.csv", na = "NA")
#
#
#
# 
# Microbial 15N data ----
# Transform to long format - 15N enriched data and natural abundance
#
# d15N values
Mic15N_d15N <- Mic15N %>%
  dplyr::select(Site, Plot, MP, Round, d15N_SE, d15N_SEF) %>%
  dplyr::rename("SE" = d15N_SE,
                "SEF" = d15N_SEF) %>%
  pivot_longer(cols = 5:6, names_to = "SE_SEF", values_to = "d15N_soil")
#
# atom% 15N
Mic15N_atom <- Mic15N %>%
  dplyr::select(Site, Plot, MP, Round, atom_pc_SE, atom_pc_SEF) %>%
  dplyr::rename("SE" = atom_pc_SE,
                "SEF" = atom_pc_SEF) %>%
  pivot_longer(cols = 5:6, names_to = "SE_SEF", values_to = "atom_pc_soil")
#
# N concentration
Mic15N_Nconc <- Mic15N %>%
  dplyr::select(Site, Plot, MP, Round, Nconc_SE, Nconc_SEF) %>%
  dplyr::rename("SE" = Nconc_SE,
                "SEF" = Nconc_SEF) %>%
  pivot_longer(cols = 5:6, names_to = "SE_SEF", values_to = "Nconc_soil")
#
# Also for natural abundance
# d15N values
Mic15N_d15N_NatAb <- Mic15N %>%
  dplyr::select(Site, Plot, MP, Round, d15N_SE_NatAb, d15N_SEF_NatAb) %>%
  dplyr::rename("SE" = d15N_SE_NatAb,
                "SEF" = d15N_SEF_NatAb) %>%
  pivot_longer(cols = 5:6, names_to = "SE_SEF", values_to = "d15N_soil_NatAb")
#
# atom% 15N
Mic15N_atom_NatAb <- Mic15N %>%
  dplyr::select(Site, Plot, MP, Round, atom_pc_SE_NatAb, atom_pc_SEF_NatAb) %>%
  dplyr::rename("SE" = atom_pc_SE_NatAb,
                "SEF" = atom_pc_SEF_NatAb) %>%
  pivot_longer(cols = 5:6, names_to = "SE_SEF", values_to = "atom_pc_soil_NatAb")
#
# N concentration
Mic15N_Nconc_NatAb <- Mic15N %>%
  dplyr::select(Site, Plot, MP, Round, Nconc_SE_NatAb, Nconc_SEF_NatAb) %>%
  dplyr::rename("SE" = Nconc_SE_NatAb,
                "SEF" = Nconc_SEF_NatAb) %>%
  pivot_longer(cols = 5:6, names_to = "SE_SEF", values_to = "Nconc_soil_NatAb")
#
# Combine
Mic15N_long <- Mic15N_d15N %>%
  left_join(Mic15N_atom, by = join_by(Site, Plot, MP, Round, SE_SEF)) %>%
  left_join(Mic15N_Nconc, by = join_by(Site, Plot, MP, Round, SE_SEF))
Mic15N_long <- Mic15N_long %>%
  left_join(Mic15N_d15N_NatAb, by = join_by(Site, Plot, MP, Round, SE_SEF)) %>%
  left_join(Mic15N_atom_NatAb, by = join_by(Site, Plot, MP, Round, SE_SEF)) %>%
  left_join(Mic15N_Nconc_NatAb, by = join_by(Site, Plot, MP, Round, SE_SEF))
#
# Save
write_csv(Mic15N_long, "clean_data/Soil_15N.csv", na = "NA")
#
#
#
# Inorganic N ----
#
# Blanks filtered from all data
Blanks <- Blanks %>%
  filter(Site == "Blank" | Site == "Water")
#
# Calculating averages of blanks
Blanks_avg <- Blanks %>%
  group_by(MP, SE_SEF) %>%
  summarise(across(c(NO3_µg_L, NH4_µg_L), ~ mean(.x, na.rm = TRUE))) %>%
  rename("Blank_NO3_1" = NO3_µg_L,
         "Blank_NH4_1" = NH4_µg_L) %>%
  ungroup()
#
# Join blanks to inorganic N values
inorgN_1 <- left_join(inorgN, Blanks_avg, by = join_by(MP, SE_SEF))
#
# Replace the missing SEF blanks with the average from SE
inorgN_1 <- inorgN_1 %>%
  mutate(Blank_NO3_1 = case_when(is.na(Blank_NO3_1) & SE_SEF == "SEF" & MP == 1 ~ inorgN_1$Blank_NO3_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 1 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NO3_1) & SE_SEF == "SEF" & MP == 2 ~ inorgN_1$Blank_NO3_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 2 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NO3_1) & SE_SEF == "SEF" & MP == 3 ~ inorgN_1$Blank_NO3_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 3 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NO3_1) & SE_SEF == "SEF" & MP == 4 ~ inorgN_1$Blank_NO3_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 4 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NO3_1) & SE_SEF == "SEF" & MP == 5 ~ inorgN_1$Blank_NO3_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 5 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 TRUE ~ Blank_NO3_1),
         Blank_NH4_1 = case_when(is.na(Blank_NH4_1) & SE_SEF == "SEF" & MP == 1 ~ inorgN_1$Blank_NH4_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 1 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NH4_1) & SE_SEF == "SEF" & MP == 2 ~ inorgN_1$Blank_NH4_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 2 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NH4_1) & SE_SEF == "SEF" & MP == 3 ~ inorgN_1$Blank_NH4_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 3 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NH4_1) & SE_SEF == "SEF" & MP == 4 ~ inorgN_1$Blank_NH4_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 4 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 is.na(Blank_NH4_1) & SE_SEF == "SEF" & MP == 5 ~ inorgN_1$Blank_NH4_1[which(inorgN_1$Site == "Abisko" & 
                                                                                                               inorgN_1$MP == 5 & 
                                                                                                               inorgN_1$Plot == 1 & 
                                                                                                               inorgN_1$SE_SEF == "SE")],
                                 TRUE ~ Blank_NH4_1))
#
# Correct values for blanks and calculate concentration per g DW
inorgN_2 <- inorgN_1 %>%
  mutate(NO3_µg_L_corr = if_else(NO3_µg_L - Blank_NO3_1 <= 0, NA, NO3_µg_L - Blank_NO3_1, missing = NO3_µg_L),
         NH4_µg_L_corr = if_else(NH4_µg_L - Blank_NH4_1 <= 0, NA, NH4_µg_L - Blank_NH4_1, missing = NH4_µg_L),
         Soil_extr_g_DW = DW_soil_m_g/FW_soil_m_g * Soil_extr_g_FW) %>%
  mutate(SW_mL = Soil_extr_g_FW - Soil_extr_g_DW) %>%
  mutate(NO3_µg_DW = (NO3_µg_L_corr/1000 * (SW_mL + 40)) / Soil_extr_g_DW,
         NH4_µg_DW = (NH4_µg_L_corr/1000 * (SW_mL + 40)) / Soil_extr_g_DW) %>% # 40mL MilliQ water added to the extraction
  mutate(NO3_µg_DW = replace_na(NO3_µg_DW, 0),
         NH4_µg_DW = replace_na(NH4_µg_DW, 0))
#
inorgN_3 <- inorgN_1 %>%
  mutate(NO3_µg_L_corr = if_else(NO3_µg_L <= 0, NA, NO3_µg_L, missing = NO3_µg_L),
         NH4_µg_L_corr = if_else(NH4_µg_L <= 0, NA, NH4_µg_L, missing = NH4_µg_L),
         Soil_extr_g_DW = DW_soil_m_g/FW_soil_m_g * Soil_extr_g_FW) %>%
  mutate(SW_mL = Soil_extr_g_FW - Soil_extr_g_DW) %>%
  mutate(NO3_µg_DW = (NO3_µg_L_corr/1000 * (SW_mL + 40)) / Soil_extr_g_DW,
         NH4_µg_DW = (NH4_µg_L_corr/1000 * (SW_mL + 40)) / Soil_extr_g_DW) %>% # 40mL MilliQ water added to the extraction
  mutate(NO3_µg_DW = replace_na(NO3_µg_DW, 0),
         NH4_µg_DW = replace_na(NH4_µg_DW, 0))


#
# Save inorganic concentrations
write_csv(inorgN_2, "clean_data/Soil_inorganic_N.csv", na = "NA")

#
#
#
# Plotting Blanks and results to find outliers or odd values
#
#
# Plot blanks as Cleveland plot
# Transform numbered months to another format
Month_yr <- tribble(~MP, ~Round,
                    01,	"01_July_19",
                    02,	"02_Aug_19",
                    03,	"03_Sep_19",
                    04,	"04_Oct_19",
                    05,	"05_Nov_19",
                    06,	"06_Dec_19",
                    07,	"07_Jan_20",
                    08,	"08_Feb_20",
                    09,	"09_Mar_20",
                    10,	"10_Apr_20",
                    11,	"11_Apr_20",
                    12,	"12_May_20",
                    13,	"13_Jun_20",
                    14,	"14_Jul_20",
                    15,	"15_Aug_20"
)
# Make MP character
Month_yr <- Month_yr %>%
  mutate(across(MP, as.character))

# combine to have blanks by Round
Blanks_plot <- Blanks %>%
  left_join(Month_yr, by = join_by(MP)) %>%
  relocate(Round, .after = MP) %>%
  mutate(Round = if_else(is.na(Round), MP, Round))

# Plot
bc1 <- Blanks_plot %>%
  #filter(SE_SEF == "SE") %>%
  ggplot() +
  geom_point(aes(y = Round,
                 x = NH4_µg_L)) +
  labs(x = "NH4 µg pr L in blank samples",
       y = "Measuring period (MP)")  #+ facet_wrap(vars(SE_SEF), scales = "free")

bc2 <- Blanks_plot %>%
  filter(SE_SEF == "SE") %>%
  ggplot() +
  geom_point(aes(y = Round,
                 x = NO3_µg_L)) +
  labs(x = "NO3 µg pr L in blank samples",
       y = "Measuring period (MP)")  #+ facet_wrap(vars(SE_SEF), scales = "free")


# Cleveland plot of uncorrected inorganic values

# combine to have blanks by Round
inorgN_plot <- inorgN_1 %>%
  left_join(Month_yr, by = join_by(MP)) %>%
  relocate(Round, .after = MP) %>%
  mutate(Round = if_else(is.na(Round), MP, Round))
inorgN_plot2 <- inorgN_2 %>%
  left_join(Month_yr, by = join_by(MP)) %>%
  relocate(Round, .after = MP) %>%
  mutate(Round = if_else(is.na(Round), MP, Round))


# Plot
nc1 <- inorgN_plot %>%
  filter(SE_SEF == "SE") %>%
  ggplot() +
  geom_point(aes(y = Round,
                 x = NH4_µg_L)) +
  labs(x = "NH4 µg pr L",
       y = "Measuring period (MP)") #+ facet_wrap(vars(SE_SEF), scales = "free")

nc2 <- inorgN_plot %>%
  filter(SE_SEF == "SE") %>%
  ggplot() +
  geom_point(aes(y = Round,
                 x = NO3_µg_L)) +
  labs(x = "NO3 µg pr L",
       y = "Measuring period (MP)") #+ facet_wrap(vars(SE_SEF), scales = "free")

# Combine
grid.arrange(bc1, bc2, nc1, nc2, ncol = 2)


bc1
nc1
inorgN_plot %>%
  filter(SE_SEF == "SE") %>%
  ggplot() +
  geom_point(aes(y = Round,
                 x = NH4_µg_L_corr)) +
  labs(x = "NH4 µg pr L corrected",
       y = "Measuring period (MP)") #+ facet_wrap(vars(SE_SEF), scales = "free")
inorgN_plot %>%
  filter(SE_SEF == "SE") %>%
  ggplot() +
  geom_point(aes(y = Round,
                 x = NH4_µg_DW)) +
  labs(x = "NH4 µg pr g DW",
       y = "Measuring period (MP)") #+ facet_wrap(vars(SE_SEF), scales = "free")
inorgN_plot2 %>%
  filter(SE_SEF == "SE") %>%
  ggplot() +
  geom_point(aes(y = Round,
                 x = NH4_µg_DW)) +
  labs(x = "NH4 µg pr g DW",
       y = "Measuring period (MP)") #+ facet_wrap(vars(SE_SEF), scales = "free")


dotchart(Blanks$NH4_µg_L,
         main="Cleveland plot - soil subsample", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = Blanks$MP, 
         groups = Blanks$Plot,
         gpch = 12, gcolor = 1)



#Plot Blanks and uncorrected inorganic values. Box-plots
b1 <- Blanks %>%
  ggplot(aes(MP, NO3_µg_L)) + geom_boxplot() + facet_wrap(vars(SE_SEF))
b2 <- Blanks %>%
  ggplot(aes(MP, NH4_µg_L)) + geom_boxplot() + facet_wrap(vars(SE_SEF))
n1 <- inorgN %>%
  ggplot(aes(MP, NO3_µg_L)) + geom_boxplot() + facet_wrap(vars(Site, SE_SEF), scales = "free")
n2 <- inorgN %>%
  ggplot(aes(MP, NH4_µg_L)) + geom_boxplot() + facet_wrap(vars(Site, SE_SEF), scales = "free")
grid.arrange(b1, b2, n1, n2, ncol = 2)


inorgN_1 %>%
  ggplot(aes(MP, NO3_µg_DW)) + geom_boxplot() + facet_wrap(vars(Site, SE_SEF), scales = "free") + labs(x = "Measuring period (MP)", y = expression(NO[3]*" "*µg^-1*" DW"), title = expression(NO[3]*" "*µg^-1*" DW with correction"))




soilExtr_2 <- soilExtr %>%
  mutate(DW_FW_frac = DW_soil_m_g/FW_soil_m_g) %>%
  mutate(NO3_µg_L)

soilExtr %>% 
  ggplot(aes(MP, NO3_µg_L)) + geom_boxplot() + facet_wrap(vars(Site,SE_SEF))
soilExtr %>% 
  ggplot(aes(MP, NH4_µg_L)) + geom_boxplot() + facet_wrap(vars(Site,SE_SEF))
soilExtr %>% 
  ggplot(aes(MP, NO3_µg_L_corr)) + geom_boxplot() + facet_wrap(vars(Site,SE_SEF))
soilExtr %>% 
  ggplot(aes(MP, NH4_µg_L_corr)) + geom_boxplot() + facet_wrap(vars(Site,SE_SEF))







#
# Combine data ----


# Combine all types of recovery into one:
vegroot15N_RLong_one <- vegroot15N_RLong %>%
  group_by(across(c("Site", "Plot", "Round"))) %>%
  summarise(TotalRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
#
# Partition only in organs
vegroot15N_RLong_Organ <- vegroot15N_RLong %>%
  group_by(across(c("Site","Plot", "Round", "Part"))) %>%
  summarise(OrganRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  dplyr::rename("Organ" = "Part") %>%
  ungroup()
vegroot15N_RLong_Organ_original <- vegroot15N_RLong_Organ
#


# Split to have species and organ
vegroot15N_bioLong <- vegroot15N_bioLong %>%
  

#
vegroot15N_d15N_Nconc_Long0 <- vegroot15N_d15N_Nconc_Long0 %>%
  add_column(Part = str_split_fixed(vegroot15N_d15N_Nconc_Long0$Type,"\\w+_",n=2)[,2]) %>%
  add_column(Species = str_split_fixed(vegroot15N_d15N_Nconc_Long0$Type,"_\\w+",n=2)[,1])
#
# Actual needed
vegroot15N_NLong1 <- vegroot15N_NLong %>%
  add_column(Part = str_split_fixed(vegroot15N_NLong$Type,"\\w+_",n=2)[,2]) %>%
  add_column(Species = str_split_fixed(vegroot15N_NLong$Type,"_\\w+",n=2)[,1])

# Get functional group and part
vegroot15N_RLong <- vegroot15N_RLong %>%
  add_column(Part = str_split_fixed(vegroot15N_RLong$Type,"\\w+_",n=2)[,2]) %>%
  add_column(Species = str_split_fixed(vegroot15N_RLong$Type,"_\\w+",n=2)[,1])


#
Rec15N <- Mic15N %>%
  left_join(vegroot15N_RLong_one) %>%
  mutate(sysRec = replace_na(R_SE, 0) + replace_na(R_MBN, 0) + replace_na(TotalRecovery, 0))
#