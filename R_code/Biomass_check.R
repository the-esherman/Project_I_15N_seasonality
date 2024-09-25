# WinterEcology I experiment
# Script author: Emil A.S. Andersen
#
# Biomass
# 
#=======  ###   Libraries    ### =======
library(plotly)
library(tidyverse)
#
#
#
#=======  ###   Load data    ### =======
#
# Core data on soils and days of labelling and harvest as well as snowdepth
coreData <- read_csv("clean_data/Core_data.csv", col_names = TRUE)
#
# Biomass, N concentration, d15N, atom% 15N for enriched and natural abundance samples, and Recovery
vegroot15N <- read_csv("clean_data/Plant_15N_data_v2.csv", col_names = TRUE)
#
#
#
#=======  ###   Biomass      ### =======
#
# The full number of unique possible rows would be:
# 15 measuring periods × 2 sites (locations) × 5 replicates × 15 categories of plant species and organs
15*2*5*15
# = 2250
# Dataset is complete
#
# A bunch of density plots
vegroot15N %>%
  mutate(across(MP, ~as.character(.x))) %>%
  filter(Site == "Vassijaure") %>%
  ggplot(aes(x = Biomass_DW_g, group = Organ, fill = Organ)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~MP)
#
# Very high values
High <- vegroot15N %>%
  filter(Biomass_DW_g > 10)
#
# Missing values, but has isotopic measures
Empty <- vegroot15N %>%
  filter(is.na(Biomass_DW_g) & !is.na(Nconc_pc))

Empty2 <- vegroot15N %>%
  filter(is.na(Biomass_DW_g) & !is.na(Atom_pc))
#
# Biomass but no isotopic value
Missing <- vegroot15N %>%
  filter(is.na(Nconc_pc) & !is.na(Biomass_DW_g)) %>%
  filter(Biomass_DW_g >= 0.05)
#
#
#
#=======  ###  { The End }   ### =======