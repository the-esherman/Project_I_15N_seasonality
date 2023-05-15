# Cleaning script
#
# - 15N vegetation and root data
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
DataName <- "raw_data/15N vegetation and roots v0.35.xlsx"
#
# Biomass, d15N, atom% 15N, and recovery ("R_") (discard later)
vegroot15N <- read_xlsx(DataName, sheet = "15N", skip = 1, col_names = TRUE)
# Natural abundance
vegrootsNatAbu <- read_xlsx(DataName, sheet = "NatAbu", col_names = TRUE)
#
# Microbial biomass, d15N, atom% 15N, and recovery ("R_") (discard later)
Mic15N <- read_xlsx(DataName, sheet = "MBN", skip = 1, col_names = TRUE)
#
vegroot15Nlong <- read_xlsx(DataName, sheet = "Long", col_names = TRUE)
#
IRMS <- read_xlsx("raw_data/IRMS_data v0.9.xlsx", col_names = TRUE)
#
Nconc <- read_xlsx(DataName, sheet = "Nconc", skip = 1, col_names = TRUE)
inorgN <- read_xlsx(DataName, sheet = "InorgN", skip = 1, col_names = TRUE)
#
# Table data for N concentrations and snow
table_dat_N_atom <- read_xlsx(DataName, sheet = "Table", skip = 1, col_names = TRUE)
#
#
# Reference atmospheric Nitrogen
# Either 0.003676 or 1/272 (more decimals)
Nair_Rst = 1/272
#
# Added 15N; mg 15N pr patch
N_add <- 1.084
#
#
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
#
#
#
# Transform to long format - 15N enriched data 
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
  left_join(vegroot15N_RLong, by = join_by(Site, Plot, MP, Round, Type))
ggplot(vegroot_long, aes(Recovery.x, Recovery.y)) + geom_point()
#
#
#
# Save 15N enriched data and Natural abundance to a data file
write_csv(vegroot_long, "clean_data/Plant_15N_data.csv", na = "NA")






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