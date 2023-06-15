# Testing of various parameters
# By Emil A.S. Andersen
# 
#=======  ###   Libraries    ### =======
library(ggpubr)
library(tidyverse)
library(readxl)
library(gridExtra)
#
#=======  ###   Load data    ### =======
#
DataName <- "raw_data/15N vegetation and roots v0.37_outlierTest.xlsx"
#
# Core data; Snow, soil mass, DW/FW of soil
core15N <- read_xlsx(DataName, sheet = "Core", skip = 1, col_names = TRUE)
#
# Vegetation data: Biomass, d15N, atom% 15N, and N concentration for enriched and natural abundance samples
vegroot15N <- read_csv("clean_data/Plant_15N_data.csv", col_names = TRUE)
#
#=======  ###    Plotting    ### =======
#
core15N %>%
  ggplot(aes(MP, Soil_mass_g)) + geom_point() + facet_wrap(vars(Site, Plot), scales = "free")
#
dotchart(core15N$Soil_mass_g, 
         main="Cleveland plot - soil mass", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = core15N$Round, 
         groups = core15N$Plot,
         gpch = 12, gcolor = 1)
#
dotchart(core15N$Soil_sub_mass_g, 
         main="Cleveland plot - soil subsample", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = core15N$Round, 
         groups = core15N$Plot,
         gpch = 12, gcolor = 1)
#
dotchart(core15N$Soil_subRF_mass_g, 
         main="Cleveland plot - soil subsample", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = core15N$Round, 
         groups = core15N$Plot,
         gpch = 12, gcolor = 1)
#
core15N <- core15N %>%
  mutate(Soil_sortMass = Soil_sub_mass_g - Soil_subRF_mass_g)
#
dotchart(core15N$Soil_sortMass, 
         main="Cleveland plot - soil subsample", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = core15N$Round, 
         groups = core15N$Plot,
         gpch = 12, gcolor = 1)
#
# Vegetation
#
vegroot15N <- vegroot15N %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round"), as.factor))
# Biomass of shoots
vegroot15N %>%
  filter(Organ == "S") %>%
  ggdotchart(x = "Round", y = "Biomass",
             group = "Round", color = "Species",
             rotate = TRUE,
             sorting = "descending",
             ggtheme = theme_bw()) + 
  facet_wrap(vars(Site)) + 
  labs(title = expression("Shoots"))
#
# Fine roots from bulk soil
vegroot15N %>%
  filter(Organ == "FR" & Species == "Root") %>%
  ggdotchart(x = "Round", y = "Biomass",
             group = "Round", color = "Plot",
             rotate = TRUE,
             sorting = "descending",
             ggtheme = theme_bw()) + 
  facet_wrap(vars(Site)) + 
  labs(title = expression("Bulk fine roots"))
#
# Coarse roots from bulk soil
vegroot15N %>%
  filter(Organ == "CR" & Species == "Root") %>%
  ggdotchart(x = "Round", y = "Biomass",
             group = "Round", color = "Plot",
             rotate = TRUE,
             sorting = "descending",
             ggtheme = theme_bw()) + 
  facet_wrap(vars(Site)) + 
  labs(title = expression("Bulk coarse roots"))
#
dotchart(vegroot15N$Biomass, 
         main="Cleveland plot - plant biomass", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = vegroot15N$Round, 
         groups = vegroot15N$Plot,
         gpch = 12, gcolor = 1)
#
# Recovery
vegroot15N %>%
  #filter(Species != "Root" & Species != "RootG") %>%
  filter(Recovery <= 0) %>%
  ggdotchart(x = "Round", y = "Recovery",
             group = "Round", color = "Species",
             rotate = TRUE,
             sorting = "descending",
             ggtheme = theme_bw()) + 
  facet_wrap(vars(Site)) + 
  labs(title = expression("Recovery"))

#
# [N]
vegroot15N %>%
  ggdotchart(x = "Round", y = "Nconc",
             group = "Round", color = "Species",
             rotate = TRUE,
             sorting = "descending",
             ggtheme = theme_bw()) + 
  facet_wrap(vars(Site)) + 
  labs(title = expression("N concentration"))

#
# APE (Atom% excess)
vegroot15N %>%
  mutate(APE = atom_pc - atom_pc_NatAb) %>%
  filter(APE <= 0) %>%
  ggdotchart(x = "Round", y = "APE",
             group = "Round", color = "Species",
             rotate = TRUE,
             sorting = "descending",
             ggtheme = theme_bw()) + 
  facet_wrap(vars(Site)) + 
  labs(title = expression("Atom% excess"))


#=======  ###   Rick roll    ### =======
# The original
https://www.youtube.com/watch?v=dQw4w9WgXcQ
#
#=======  ###  { The End }   ### =======