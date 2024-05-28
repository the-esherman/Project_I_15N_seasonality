# Cleaning script - Project I WinterEcology
# Version 2.0
#
# using all base information
#
# 15N vegetation, root, and soil data
# By Emil A.S. Andersen
# 
#------- ### Libraries ### -------
library(tidyverse)
library(readxl)
#
#
#
#------- ### Load data ### -------
#
DataName <- "raw_data/base data v1.3_R.xlsx"
#
# Core data:
# Harvest information, snow measurements, core diameter, soil depth and mass (and sub-mass), label values, soil moisture, sub-sample root mass
Core <- read_xlsx(DataName, sheet = "Core", skip = 1, col_names = TRUE, na = "NA")
#
# Vegetation biomass (shoots and roots)
vegrootMass <- read_xlsx(DataName, sheet = "Vegetation_data", skip = 1, col_names = TRUE, na="NA")
#
# Vegetation 15N data: delta15N, [N], Atom%
vegroot15N <- read_xlsx(DataName, sheet = "15N_data", skip = 1, col_names = TRUE, na = "NA")
#
# Soil extraction mass
extrMass <- read_xlsx(DataName, sheet = "Extracts", skip = 1, col_names = TRUE, na = "NA")
#
# Soil extraction 15N data: delta15N, [N], Atom%
extr15N <- read_xlsx(DataName, sheet = "Extr_15N", skip = 1, col_names = TRUE, na = "NA")
extr15N_last <- read_xlsx(DataName, sheet = "Extr_15N_last", skip = 1, col_names = TRUE, na = "")
#
# Soil extraction inorganic N concentrations and for blanks
extrInorgN <- read_xlsx(DataName, sheet = "Extr_inorgN_raw", skip = 1, col_names = TRUE, na = "NA")
extrInorgN_blanks <- read_xlsx(DataName, sheet = "Extr_inorgN_blanks", skip = 1, col_names = TRUE, na = "NA")
#
# Soil extractions TDN and blanks
extrTDN <- read_xlsx(DataName, sheet = "Extr_TDN", skip = 1, col_names = TRUE, na = "NA")
extrTDN_blanks <- read_xlsx(DataName, sheet = "Extr_TDN_blanks", skip = 1, col_names = TRUE, na = "NA")
#
# 
# Control values
# Vegetation biomass
vegrootMass_control <- read_xlsx(DataName, sheet = "Vegetation_control", skip = 1, col_names = TRUE, na = "NA")
#
# Vegetation 15N data: delta15N, [N], Atom%
vegroot15N_control <- read_xlsx(DataName, sheet = "15N_control", skip = 1, col_names = TRUE, na = "NA")
#
# Soil extraction 15N data: delta15N, [N], Atom%
extr15N_control <- read_xlsx(DataName, sheet = "Extr_15N_control", skip = 1, col_names = TRUE, na = "NA")
#
# Soil extraction inorganic N concentrations and for blanks
extrInorgN_control <- read_xlsx(DataName, sheet = "Extr_inorgN_control", skip = 1, col_names = TRUE, na = "NA")
#
#
#
#------- ### Cleaning  ### -------
# Label atomic data ----
# Added 15N; mg 15N pr patch
N_add <- 1.084
Label_atom_pc <- 0.987 # 98.7% double labelled 15N-NH4NO3
Atom_mass_14N_NH4NO3 <- 2*14.007+4*1.008+3*15.999 # The atomic mass of 14N NH4NO3
Atom_mass_15N_NH4NO3 <- 2*15+4*1.008+3*15.999 # The atomic mass of double 15N NH4NO3

#(2*15*Label_atom_pc)/(2*14.007*(1-Label_atom_pc)+2*15*Label_atom_pc+4*1.008+3*15.999)

Label_15N_frac <- (2*15*Label_atom_pc)/(Atom_mass_14N_NH4NO3*(1-Label_atom_pc)+Atom_mass_15N_NH4NO3*Label_atom_pc) # The atom mass of the 15N to the total label NH4NO3
Label_N_frac <- (2*15*Label_atom_pc+2*14*(1-Label_atom_pc))/(Atom_mass_14N_NH4NO3*(1-Label_atom_pc)+Atom_mass_15N_NH4NO3*Label_atom_pc)
# In the numerator: 2 15N per molecule, but only 98.7%
# In the denominator: The molecule's average mass, as 14N NH4NO3 is 1.3% and 15N NH4NO3 is 98.7%

#
#
#
# Core calculations ----
Core2 <- Core %>%
  mutate(Soil_vol_cm3 = (Soil_diameter_cm/2)^2*pi*Soil_depth_cm, # Soil core vol cm3
         Soil_area_cm2 = (Soil_diameter_cm/2)^2*pi, # Patch area cm2
         Soil_ratio = Soil_sub_mass_g/Soil_mass_g, # Subsample to full core ratio
         DW_FW_frac = SM_DW_g/SM_FW_g, # The soil Dry weight to Fresh weight ratio
         Injection_N_mg_pr_patch = Label_conc_mg/(Label_conc_mL/1000)*(19/1000)*Label_N_frac, # Actual added amount of N, assuming 19mL vol of solution per plot, 1mL per hole
         Injection_15N_mg_pr_patch = Label_conc_mg/(Label_conc_mL/1000)*(19/1000)*Label_15N_frac, # Actual added amount of labelled 15N, assuming 19mL vol of solution per plot, 1mL per hole
         Soil_RF_FW_g = Soil_subRF_mass_g/Soil_sub_mass_g) %>% # Mass of root free soil in FW
  mutate(Soil_RF_DW_g = Soil_RF_FW_g*Soil_mass_g*DW_FW_frac) # Mass of root free soil in DW

#
Core3 <- Core2 %>%
  filter(MP != "EX") %>%
  mutate(across(c(Day_of_label, Day_of_harvest), ymd)) %>%
  mutate(DaysLH = Day_of_harvest - Day_of_label) %>%
  mutate(across(c(MP, DaysLH), as.numeric))
# Calculate days between harvests
Core3_A <- Core3 %>%
  select(Site, Plot, MP) %>%
  filter(Site == "Abisko")
Core3_V <- Core3 %>%
  select(Site, Plot, MP) %>%
  filter(Site == "Vassijaure")
Core3_A["DaysHH"] <- NA
Core3_V["DaysHH"] <- NA
# Double check I can count
ct<- vector("double")
# Loop over each to calculate rate from one round to the next for the same plot and then around 7 days if when labeling happens
for (i in 2:15) {
  for (j in 1:5) {
    ct[5*(i-1)+j] <- 5*(i-1)+j
    Core3_A$DaysHH[5*(i-1)+j] <- (Core3$Day_of_harvest[which(Core3$Site == "Abisko" & 
                                                            Core3$Plot == j & 
                                                            Core3$MP == (i))] - 
                                    Core3$Day_of_harvest[which(Core3$Site == "Abisko" & 
                                                              Core3$Plot == j & 
                                                              Core3$MP == (i-1))])
    Core3_V$DaysHH[5*(i-1)+j] <- (Core3$Day_of_harvest[which(Core3$Site == "Vassijaure" & 
                                                            Core3$Plot == j & 
                                                            Core3$MP == i)] - 
                                    Core3$Day_of_harvest[which(Core3$Site == "Vassijaure" & 
                                                              Core3$Plot == j & 
                                                              Core3$MP == (i-1))])
  }
}

Core3_AV <- Core3_A %>%
  bind_rows(Core3_V)

Core4 <- Core3 %>%
  left_join(Core3_AV, by = join_by(Site, Plot, MP)) %>%
  relocate(c(DaysLH, DaysHH), .after = Day_of_harvest) %>%
  mutate(DaysHL = DaysHH - DaysLH)

# Save as csv
write_csv(Core4, "clean_data/Core_data.csv")

#
#
#
# Vegetation ----

# Pivot longer for bulk root samples
vegrootBulk <- Core2 %>%
  select(Site, Plot, MP, Soil_mass_g, Soil_sub_mass_g, Bulk_FR, Bulk_CR, BulkG_FR, Soil_ratio) %>%
  mutate(Bulk_FR = Bulk_FR/Soil_ratio,
         Bulk_CR = Bulk_CR/Soil_ratio,
         BulkG_FR = BulkG_FR/Soil_ratio) %>% # Scale subsample to full core
  select(Site, Plot, MP, Bulk_FR, Bulk_CR, BulkG_FR) %>%
  pivot_longer(cols = 4:6, names_to = "Organ", values_to = "Biomass_DW_g")

# Pivot longer for 15N label added
vegrootInject <- Core2 %>%
  select(Site, Plot, MP, Injection_15N_mg_pr_patch) %>%
  filter(MP != "EX") %>%
  mutate(across(MP, as.numeric))

# Combine bulk roots with rest of root mass
vegroot <- vegrootBulk %>%
  filter(MP != "EX") %>%
  mutate(across(MP, as.numeric)) %>%
  bind_rows(vegrootMass)

# Combine 15N data with biomass
vegroot <- vegroot %>%
  left_join(vegroot15N, by = join_by(Site, Plot, MP, Organ)) %>%
  left_join(vegrootInject, by = join_by(Site, Plot, MP))

# Combine labelled with natural abundance
vegroot15N_control_sum <- vegroot15N_control %>%
  group_by(Site, Organ) %>%
  summarise(NatAb_d15N = mean(d15N, na.rm = TRUE), 
            NatAb_atom = mean(Atom_pc, na.rm = TRUE), # Average with atom% instead of d15N ?
            .groups = "keep") %>%
  mutate(NatAb_atom_pc = 100/(1+272/(1+NatAb_d15N/1000)))

x <- vegroot15N_control_sum %>%
  mutate(diff = NatAb_atom - NatAb_atom_pc) # it makes no noticeable difference if atom% is calculated from an average d15N or from an average atom%

vegroot2 <- vegroot %>%
  left_join(vegroot15N_control_sum, by = join_by(Site, Organ)) %>%
  select(1:9, NatAb_atom_pc) %>%
  rename(Type = Organ)
vegroot2 <- vegroot2 %>%
  add_column(Organ = str_split_fixed(vegroot2$Type,"\\w+_",n=2)[,2], .after = "Type") %>%
  add_column(Species = str_split_fixed(vegroot2$Type,"_\\w+",n=2)[,1], .after = "Type")
vegroot2 <- vegroot2 %>%
  mutate(Recovery = ((Atom_pc - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/(N_add/1000) * 100) %>%
  mutate(Recovery2 = ((Atom_pc - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/(Injection_15N_mg_pr_patch/1000) * 100) %>%
  mutate(Recovery = if_else(Recovery < 0, 0, Recovery),
         Recovery2 = if_else(Recovery2 < 0, 0, Recovery2))
vegroot_output <- vegroot2 %>%
  select(1:6, Biomass_DW_g, Nconc_pc, d15N, Atom_pc, Injection_15N_mg_pr_patch, NatAb_atom_pc, Recovery2) %>%
  rename("Recovery" = Recovery2)

# write a csv with the data
write_csv(vegroot_output, "clean_data/Plant_15N_data_v2.csv", na = "NA")

#
#
#
# Soil N - TDN and inorganic N - ----
# » TDN « ----
#
# From the core dataset, extract Dry weight/fresh weight fraction
Core_extr <- Core2 %>%
  select(Site, Plot, MP, DW_FW_frac)

# Add DW/FW fraction to extraction mass (g FW)
extrMass_DW <- extrMass %>%
  left_join(Core_extr, by = join_by(Site, Plot, MP)) %>%
  mutate(Soil_DW_g = Soil_FW_g*DW_FW_frac,
         SW_mL = Soil_FW_g - Soil_DW_g)

# Remove control round
extrMass2 <- extrMass_DW %>%
  filter(MP != "EX") %>%
  mutate(across(MP, as.numeric))

# Add mass to extraction and calculate the [N] (TDN from IRMS)
extr15N_Nconc <- extr15N %>%
  select(1:10) %>%
  bind_rows(extr15N_last) %>% # Append the previously missing samples to the end
  left_join(extrMass2, by = join_by(Site, Plot, MP, Extr_type)) %>%
  mutate(Nconc_microg2 = ((Nconc_sample_microg/FreezeD_mL)*(MilliQ_mL+SW_mL))/Soil_DW_g)

# TDN (from quAAtro) values for MP 9, 13, and 15 Nconc (same as [N] or Nconc, but measured slightly differently and thus in a different sheet and unit!!!)
# Blanks
extrTDN_blanks_avg <- extrTDN_blanks %>%
  group_by(MP, Extr_type) %>%
  summarise(Blank_TDN = mean(TDN_sample_mg_pr_L, na.rm = TRUE), .groups = "keep") %>%
  ungroup()

# Add mass and blanks and calculate TDN as µg pr g DW
extrTDN_1 <- extrTDN %>%
  mutate(Plot = if_else(Site == "Abisko" & MP == "13" & Extr_type == "SEF" & TDN_sample_mg_pr_L >= 20, 1, Plot)) %>% #  # # # OBS!!!
  left_join(extrMass2, by = join_by(Site, Plot, MP, Extr_type)) %>%
  left_join(extrTDN_blanks_avg, by = join_by(MP, Extr_type)) %>%
  mutate(TDN_sample_mg_pr_L_corr = if_else(TDN_sample_mg_pr_L - Blank_TDN <= 0, NA, TDN_sample_mg_pr_L - Blank_TDN, missing = TDN_sample_mg_pr_L)) %>%
  mutate(TDN_mg_pr_gDW = (TDN_sample_mg_pr_L/1000 * (SW_mL + 40)) / Soil_DW_g) %>% # 40mL MilliQ water added to the extraction
  mutate(TDN_mg_pr_gDW = replace_na(TDN_mg_pr_gDW, 0)) %>%
  mutate(TDN_microg_pr_gDW = TDN_mg_pr_gDW*1000)

# Select only parts that match [N] (TDN from IRMS) and rename to fit
extrTDN_2 <- extrTDN_1 %>%
  select(1:4, Soil_FW_g, DW_FW_frac, Soil_DW_g, SW_mL, TDN_microg_pr_gDW) %>%
  rename(Nconc_microg3 = TDN_microg_pr_gDW)

# Combine all TDN [N]
extr15N_Nconc_1 <- extr15N_Nconc %>%
  left_join(extrTDN_2, by = join_by(Site, Plot, MP, Extr_type, Soil_FW_g, DW_FW_frac, Soil_DW_g, SW_mL)) %>%
  mutate(Nconc_microg2 = if_else(is.na(Nconc_microg), Nconc_microg3, Nconc_microg2)) %>% # We don't have exact values for the missing or last samples that were analysed with quAAtro instead of IRMS. So results are replaced there with TDN values
  select(1:15)

# Controls of TDN
# Average the control atom% by site and type of extraction
extr15N_control_avg <- extr15N_control %>%
  group_by(across(c("Site", "Extr_type"))) %>%
  summarise(Atom_pc_NatAb = mean(Atom_pc, na.rm = TRUE), .groups = "keep") %>%
  ungroup()

# Add control atom% to rest
extr15N_TDN <- extr15N_Nconc_1 %>%
  left_join(extr15N_control_avg, by = join_by(Site, Extr_type))

# Narrow down to what is essential
extr15N_TDN_essential <- extr15N_TDN %>%
  select(1:4, d15N, Atom_pc, Nconc_microg2, Atom_pc_NatAb) %>%
  relocate(Nconc_microg2, .before = d15N) %>%
  rename("Nconc_microg_pr_gDW" = Nconc_microg2)

# Save as csv
#write_csv(extr15N_output, "clean_data/Soil_15N_v2.csv", na = "NA")

#
#
# » Inorganic N « ----
#
# Inorganic N values - 
# raw data cleaned by changing certain wrong or missing names and removing SEF_old
extrInorgN_clean_0 <- extrInorgN %>%
  mutate(MP = case_when(Site == "Water" ~ "1",
                        MP == "EX (w)" ~ "EX",
                        TRUE ~ MP),
         Extr_type = case_when(Site == "Blank" & MP == "9" & is.na(Extr_type) ~ "SE",
                               Site == "Blank" & MP == "13" & is.na(Extr_type) ~ "SEF",
                               Site == "Blank" & MP == "15" & is.na(Extr_type) ~ "SE",
                               Site == "Water" ~ "SE",
                               TRUE ~ Extr_type)) %>%
  mutate(across(NO3_sample_microg_pr_L, ~ na_if(.x, "n.d."))) %>%
  mutate(across(NO3_sample_microg_pr_L, ~ as.numeric(.x))) %>%
  mutate(Site = if_else(Site == "Water", "Blank", Site)) %>%
  filter(Extr_type != "SEF_old")

# Check duplicates
extrInorgN_clean_1 <- extrInorgN_clean_0 %>%
  unite("ID", 1:4, remove = FALSE)
duplicated_IDS_inorgN <- extrInorgN_clean_1[duplicated(extrInorgN_clean_1$ID),]  
duplicate_subset_inorgN <- extrInorgN_clean_1[extrInorgN_clean_1$ID %in% duplicated_IDS_inorgN$ID, ]

test <- extrInorgN_clean_0 %>%
  filter(Site != "Blank")
test %>%
  ggplot(aes(MP)) + geom_bar()

# V_1_4_SE and V_4_6_SEF have true duplicates. An average value should be taken as their value (as values are very similar this should not change much)
# A_2_13_SEF is a duplicate where one is A_1_13_SEF
# The remaining duplicates are blanks. But with very different values since they were done one "Abisko" or "Vassijaure" samples.

# Average V_1_4_SE and V_4_6_SEF duplicates. Change one A_2_13_SEF to A_1_13_SEF
# Blanks have been removed
extrInorgN_clean_2 <- extrInorgN_clean_0 %>%
  filter(Site != "Blank") %>%
  mutate(Plot = if_else(Site == "Abisko" & MP == "13" & Extr_type == "SEF" & NH4_sample_microg_pr_L >= 8000, 1, Plot)) %>% #  # # # OBS!!!
  group_by(Site, Plot, MP, Extr_type) %>%
  mutate(NO3_sample_microg_pr_L = mean(NO3_sample_microg_pr_L),
         NH4_sample_microg_pr_L = mean(NH4_sample_microg_pr_L)) %>%
  slice(1) %>%
  ungroup()

dotchart(extrInorgN_clean_2$NH4_sample_microg_pr_L, 
         main="Cleveland plot - NH4", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_clean_2$MP, 
         groups = extrInorgN_clean_2$Plot,
         gpch = 12, gcolor = 1)
dotchart(extrInorgN_clean_2$NO3_sample_microg_pr_L, 
         main="Cleveland plot - NO3", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_clean_2$MP, 
         groups = extrInorgN_clean_2$Plot,
         gpch = 12, gcolor = 1)

# A few outliers of very high values, but uncorrected for blanks

# Blanks
extrInorgN_blanks <- extrInorgN_clean_0 %>%
  filter(Site == "Blank")

# Calculating averages of blanks
extrInorgN_blanks_avg <- extrInorgN_blanks %>%
  group_by(MP, Extr_type) %>%
  summarise(across(c(NO3_sample_microg_pr_L, NH4_sample_microg_pr_L), ~ mean(.x, na.rm = TRUE)), .groups = "keep") %>%
  rename("Blank_NO3" = NO3_sample_microg_pr_L,
         "Blank_NH4" = NH4_sample_microg_pr_L) %>%
  ungroup()

# Split into "labelled" blanks and "control" blanks
extrInorgN_blanks_avg_label <- extrInorgN_blanks_avg %>%
  filter(MP != "EX") %>%
  mutate(across(MP, as.numeric))
extrInorgN_blanks_avg_control <- extrInorgN_blanks_avg %>%
  filter(MP == "EX")



extrInorgN_blanks.1 <- extrInorgN_blanks %>%
  filter(Extr_type == "SE")
extrInorgN_blanks.2 <- extrInorgN_blanks %>%
  filter(Extr_type == "SEF")

dotchart(extrInorgN_blanks.1$NH4_sample_microg_pr_L, 
         main="Cleveland plot - NH4 SE blanks", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_blanks.1$MP, 
         groups = extrInorgN_blanks.1$Plot,
         gpch = 12, gcolor = 1)
dotchart(extrInorgN_blanks.2$NH4_sample_microg_pr_L, 
         main="Cleveland plot - NH4 SEF blanks", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_blanks.2$MP, 
         groups = extrInorgN_blanks.2$Plot,
         gpch = 12, gcolor = 1)


extrInorgN_clean_2.1 <- extrInorgN_clean_2 %>%
  filter(Extr_type == "SE")
extrInorgN_clean_2.2 <- extrInorgN_clean_2 %>%
  filter(Extr_type == "SEF")

dotchart(extrInorgN_clean_2.1$NH4_sample_microg_pr_L, 
         main="Cleveland plot - NH4 SE", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_clean_2.1$MP, 
         groups = extrInorgN_clean_2.1$Plot,
         gpch = 12, gcolor = 1)
dotchart(extrInorgN_clean_2.2$NH4_sample_microg_pr_L, 
         main="Cleveland plot - NH4 SE", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_clean_2.2$MP, 
         groups = extrInorgN_clean_2.2$Plot,
         gpch = 12, gcolor = 1)




# 
extrInorgN_clean <- extrInorgN_clean_2 %>%
  filter(MP != "EX") %>%
  mutate(across(MP, as.numeric))

# Join blanks to inorganic N values
extrInorgN_1 <- extrInorgN_clean %>%
  left_join(extrMass2, by = join_by(Site, Plot, MP, Extr_type)) %>%
  left_join(extrInorgN_blanks_avg_label, by = join_by(MP, Extr_type))

# Replace the missing SEF blanks with the average from SE. This is for the first few MPs when only one set of blanks were done
extrInorgN_1 <- extrInorgN_1 %>%
  mutate(Blank_NO3_1 = case_when(is.na(Blank_NO3) & Extr_type == "SEF" & MP == 1 ~ extrInorgN_1$Blank_NO3[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 1 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NO3) & Extr_type == "SEF" & MP == 2 ~ extrInorgN_1$Blank_NO3[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 2 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NO3) & Extr_type == "SEF" & MP == 3 ~ extrInorgN_1$Blank_NO3[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 3 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NO3) & Extr_type == "SEF" & MP == 4 ~ extrInorgN_1$Blank_NO3[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 4 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NO3) & Extr_type == "SEF" & MP == 5 ~ extrInorgN_1$Blank_NO3[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 5 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 TRUE ~ Blank_NO3),
         Blank_NH4_1 = case_when(is.na(Blank_NH4) & Extr_type == "SEF" & MP == 1 ~ extrInorgN_1$Blank_NH4[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 1 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NH4) & Extr_type == "SEF" & MP == 2 ~ extrInorgN_1$Blank_NH4[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 2 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NH4) & Extr_type == "SEF" & MP == 3 ~ extrInorgN_1$Blank_NH4[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 3 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NH4) & Extr_type == "SEF" & MP == 4 ~ extrInorgN_1$Blank_NH4[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 4 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 is.na(Blank_NH4) & Extr_type == "SEF" & MP == 5 ~ extrInorgN_1$Blank_NH4[which(extrInorgN_1$Site == "Abisko" & 
                                                                                                                      extrInorgN_1$MP == 5 & 
                                                                                                                      extrInorgN_1$Plot == 1 & 
                                                                                                                      extrInorgN_1$Extr_type == "SE")],
                                 TRUE ~ Blank_NH4))

# Correct values based on blanks, replacing negative values with NA. Further calculating the concentration to µg pr g DW
extrInorgN_1 <- extrInorgN_1 %>%
  mutate(NO3_sample_microg_pr_L_corr = if_else(NO3_sample_microg_pr_L - Blank_NO3 <= 0, NA, NO3_sample_microg_pr_L - Blank_NO3, missing = NO3_sample_microg_pr_L),
         NH4_sample_microg_pr_L_corr = if_else(NH4_sample_microg_pr_L - Blank_NH4 <= 0, NA, NH4_sample_microg_pr_L - Blank_NH4, missing = NH4_sample_microg_pr_L)) %>%
  mutate(NO3_microg_pr_gDW = (NO3_sample_microg_pr_L_corr/1000 * (SW_mL + 40)) / Soil_DW_g,
         NH4_microg_pr_gDW = (NH4_sample_microg_pr_L_corr/1000 * (SW_mL + 40)) / Soil_DW_g) %>% # 40mL MilliQ water added to the extraction
  mutate(NO3_microg_pr_gDW = replace_na(NO3_microg_pr_gDW, 0),
         NH4_microg_pr_gDW = replace_na(NH4_microg_pr_gDW, 0))
extrInorgN_essential <- extrInorgN_1 %>%
  select(1:4, NO3_microg_pr_gDW, NH4_microg_pr_gDW)

# Outliers
dotchart(extrInorgN_essential$NH4_microg_pr_gDW, 
         main="Cleveland plot - NH4", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_essential$MP, 
         groups = extrInorgN_essential$Plot,
         gpch = 12, gcolor = 1)
dotchart(extrInorgN_essential$NO3_microg_pr_gDW, 
         main="Cleveland plot - NO3", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = extrInorgN_essential$MP, 
         groups = extrInorgN_essential$Plot,
         gpch = 12, gcolor = 1)

# Outliers have changed for NH4, where only 1 value is a bit further out. The NO3 outlier is still present


# Save inorganic N to csv
#write_csv(extrInorgN_output, "clean_data/Soil_inorganic_N.csv", na = "NA")

#
#
# » TDN + inorganic N « ----
#
# Combine extract data in one file (inorganic [N] and TDN 15N data)
extr15N_all <- extr15N_TDN_essential %>%
  left_join(extrInorgN_essential, by = join_by(Site, Plot, MP, Extr_type))

#
#
# » Microbial mass and saving « ---- 
#
# Select the microbial mass, or root-free soil core dry weight
Mic_mass <- Core2 %>%
  select(Site, Plot, MP, Soil_RF_DW_g) %>%
  filter(MP != "EX") %>%
  mutate(across(MP, as.numeric))

# Combine with all else
extr15N_all <- extr15N_all %>%
  left_join(Mic_mass, by = join_by(Site, Plot, MP))

# Save to csv
write_csv(extr15N_all, "clean_data/Soil_N.csv", na = "NA")

#
#
#
#
#
# Graphs for checking calculations against old ----
# Checking if the old recovery matches the new calculations

# library(ggpubr)

firstVegroot <- read_csv("clean_data/Plant_15N_data.csv", col_names = TRUE)

vegroot3 <- vegroot2 %>%
  mutate(Species = case_when(Species == "Bulk" ~ "Root",
                             Species == "BulkG" ~"RootG",
                             TRUE ~ Species)) %>%
  left_join(firstVegroot, by = join_by(Site, Plot, MP, Species, Organ))

vegroot3 %>%
  ggplot(aes(Recovery.x, Recovery.y)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(label.y = 10, aes(label = after_stat(eq.label))) + stat_regline_equation(label.y = 9, aes(label = after_stat(rr.label)))
vegroot3 %>%
  ggplot(aes(Biomass_DW_g, Biomass)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(label.y = 10, aes(label = after_stat(eq.label))) + stat_regline_equation(label.y = 9, aes(label = after_stat(rr.label)))
# Some weird outliers because of mismatch in biomass. The other data set has outdated values somehow.


# Changed the value of mg 15N pr patch
vegroot3 %>%
  ggplot(aes(Recovery.x, Recovery2)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(label.y = 10, aes(label = after_stat(eq.label))) + stat_regline_equation(label.y = 9, aes(label = after_stat(rr.label)))
vegroot3 %>%
  ggplot(aes(Recovery.y, Recovery2)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(label.y = 10, aes(label = after_stat(eq.label))) + stat_regline_equation(label.y = 9, aes(label = after_stat(rr.label)))

# Show only where the difference between the two Recovery values is not 0
test <- vegroot3 %>%
  mutate(Diff = Recovery.x - Recovery.y) %>%
  filter(Diff != 0)


#
#
#

# Checking the [N] as calculated in excel and in R
# Essentially no difference, and a regression line close to y = x
extr15N_Nconc %>%
  ggplot(aes(Nconc_microg, Nconc_microg2)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(label.y = 750, aes(label = after_stat(eq.label))) + stat_regline_equation(label.y = 700, aes(label = after_stat(rr.label)))

# Show only where the difference between the two concentrations is not 0
Test <- extr15N_Nconc %>%
  mutate(Diff = Nconc_microg - Nconc_microg2) %>%
  filter(Diff != 0)


#
#
#
#
#
# Control 15N graphs ----

vegroot15N %>%
  ggplot(aes(Atom_pc, color = Organ)) + geom_histogram() + facet_wrap(vars(Site), scales = "free")
vegroot15N %>%
  ggplot(aes(d15N, color = Site)) + geom_histogram(binwidth = 1) + facet_wrap(vars(Organ), scales = "free") + coord_cartesian(xlim = c(-10,5))



vegroot15N_control %>%
  ggplot(aes(Atom_pc, color = Organ)) + geom_histogram() + facet_wrap(vars(Site), scales = "free")

vegroot15N_control %>%
  ggplot(aes(d15N, color = Site)) + geom_histogram(binwidth = 1) + facet_wrap(vars(Organ), scales = "free")



# Visualizing the difference in natural abundance occurrence in organs over time and space
# C vs EX (September vs March)
x <- vegroot15N_control %>%
  select(1:4, d15N) %>%
  mutate(MP = if_else(MP == "?", "C", MP)) %>%
  distinct(Site, Plot, MP, Organ, .keep_all = TRUE) %>%
  pivot_wider(names_from = MP, values_from = d15N) 
x %>%
  ggplot(aes(EX, C, color = Organ)) + geom_point() + facet_wrap(vars(Site), scales = "free")
x %>%
  ggplot(aes(EX, C, color = Site)) + geom_point() + facet_wrap(vars(Organ), scales = "free")
# Abisko vs Vassijaure
y <- vegroot15N_control %>%
  select(1:4, d15N) %>%
  mutate(MP = if_else(MP == "?", "C", MP)) %>%
  distinct(Site, Plot, MP, Organ, .keep_all = TRUE) %>%
  pivot_wider(names_from = Site, values_from = d15N) 
y %>%
  ggplot(aes(Abisko, Vassijaure, color = Organ)) + geom_point() + facet_wrap(vars(MP), scales = "free")
y %>%
  ggplot(aes(Abisko, Vassijaure, color = MP)) + geom_point() + facet_wrap(vars(Organ), scales = "free")
#
# Control times vs each other
vegroot15N_control %>%
  # Remove two duplicates from EX. Very similar values
  filter(Site != "Vassijaure" | Plot != 1 | MP != "EX" | Organ != "DS_S" | Nconc_pc <= 1) %>% # V_1_DS_S
  filter(Site != "Abisko" | Plot != 4 | MP != "EX" | Organ != "ES_CR" | d15N <= -7.9) %>% # A_4_ES_CR
  select(1:4, d15N) %>%
  filter(MP != "?") %>%
  complete(Site, Plot, MP, Organ) %>%
  pivot_wider(names_from = MP, values_from = d15N) %>%
  ggplot(aes(x = C, y = EX, color = Site)) + 
  geom_point() +
  #scale_colour_viridis_d() + # option = "H"
  labs(title = expression("Control "*delta^15*"N")) + 
  theme_bw(base_size = 20)
#
# Control sites vs each other
vegroot15N_control %>%
  # Remove two duplicates from EX. Very similar values
  filter(Site != "Vassijaure" | Plot != 1 | MP != "EX" | Organ != "DS_S" | Nconc_pc <= 1) %>% # V_1_DS_S
  filter(Site != "Abisko" | Plot != 4 | MP != "EX" | Organ != "ES_CR" | d15N <= -7.9) %>% # A_4_ES_CR
  select(1:4, d15N) %>%
  filter(MP != "?") %>%
  complete(Site, Plot, MP, Organ) %>%
  pivot_wider(names_from = Site, values_from = d15N) %>%
  # Plot
  ggplot(aes(x = Abisko, y = Vassijaure, color = MP)) + 
  geom_point() +
  scale_colour_viridis_d() + # option = "H"
  labs(title = expression("Control "*delta^15*"N")) + 
  theme_bw(base_size = 20)
