# Vegetation cover
#
library(tidyverse)
library(readxl)
#
#
VegInvName <- "raw_data/Vegetation inventarization_15N_updated.xlsx"
A_1 <- read_excel(VegInvName, sheet = "A_1", col_names = TRUE)
A_2 <- read_excel(VegInvName, sheet = "A_2", col_names = TRUE)
A_3 <- read_excel(VegInvName, sheet = "A_3", col_names = TRUE)
A_4 <- read_excel(VegInvName, sheet = "A_4", col_names = TRUE)
A_5 <- read_excel(VegInvName, sheet = "A_5", col_names = TRUE)
V_1 <- read_excel(VegInvName, sheet = "V_1", col_names = TRUE)
V_2 <- read_excel(VegInvName, sheet = "V_2", col_names = TRUE)
V_3 <- read_excel(VegInvName, sheet = "V_3", col_names = TRUE)
V_4 <- read_excel(VegInvName, sheet = "V_4", col_names = TRUE)
V_5 <- read_excel(VegInvName, sheet = "V_5", col_names = TRUE)
#
# Filter classes and total coverage (averaged in excel), and rename unnamed columns
A_1.x <- A_1 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("A_1_1_pres" = ...3,
         "A_1_2_pres" = ...5,
         "A_1_3_pres" = ...7,
         "A_1_4_pres" = ...9,
         "A_1_5_pres" = ...11,
         "A_1_6_pres" = ...13,
         "A_1_7_pres" = ...15,
         "A_1_8_pres" = ...17,
         "A_1_9_pres" = ...19,
         "A_1_10_pres" = ...21,
         "A_1_11_pres" = ...23,
         "A_1_12_pres" = ...25,
         "A_1_13_pres" = ...27,
         "A_1_14_pres" = ...29,
         "A_1_15_pres" = ...31)
A_2.x <- A_2 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("A_2_1_pres" = ...3,
         "A_2_2_pres" = ...5,
         "A_2_3_pres" = ...7,
         "A_2_4_pres" = ...9,
         "A_2_5_pres" = ...11,
         "A_2_6_pres" = ...13,
         "A_2_7_pres" = ...15,
         "A_2_8_pres" = ...17,
         "A_2_9_pres" = ...19,
         "A_2_10_pres" = ...21,
         "A_2_11_pres" = ...23,
         "A_2_12_pres" = ...25,
         "A_2_13_pres" = ...27,
         "A_2_14_pres" = ...29,
         "A_2_15_pres" = ...31)
A_3.x <- A_3 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("A_3_1_pres" = ...3,
         "A_3_2_pres" = ...5,
         "A_3_3_pres" = ...7,
         "A_3_4_pres" = ...9,
         "A_3_5_pres" = ...11,
         "A_3_6_pres" = ...13,
         "A_3_7_pres" = ...15,
         "A_3_8_pres" = ...17,
         "A_3_9_pres" = ...19,
         "A_3_10_pres" = ...21,
         "A_3_11_pres" = ...23,
         "A_3_12_pres" = ...25,
         "A_3_13_pres" = ...27,
         "A_3_14_pres" = ...29,
         "A_3_15_pres" = ...31)
A_4.x <- A_4 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("A_4_1_pres" = ...3,
         "A_4_2_pres" = ...5,
         "A_4_3_pres" = ...7,
         "A_4_4_pres" = ...9,
         "A_4_5_pres" = ...11,
         "A_4_6_pres" = ...13,
         "A_4_7_pres" = ...15,
         "A_4_8_pres" = ...17,
         "A_4_9_pres" = ...19,
         "A_4_10_pres" = ...21,
         "A_4_11_pres" = ...23,
         "A_4_12_pres" = ...25,
         "A_4_13_pres" = ...27,
         "A_4_14_pres" = ...29,
         "A_4_15_pres" = ...31)
A_5.x <- A_5 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("A_5_1_pres" = ...3,
         "A_5_2_pres" = ...5,
         "A_5_3_pres" = ...7,
         "A_5_4_pres" = ...9,
         "A_5_5_pres" = ...11,
         "A_5_6_pres" = ...13,
         "A_5_7_pres" = ...15,
         "A_5_8_pres" = ...17,
         "A_5_9_pres" = ...19,
         "A_5_10_pres" = ...21,
         "A_5_11_pres" = ...23,
         "A_5_12_pres" = ...25,
         "A_5_13_pres" = ...27,
         "A_5_14_pres" = ...29,
         "A_5_15_pres" = ...31)
#
# Merge all Abisko
Abisko <- A_1.x %>%
  left_join(A_2.x, by = join_by(Species)) %>%
  left_join(A_3.x, by = join_by(Species)) %>%
  left_join(A_4.x, by = join_by(Species)) %>%
  left_join(A_5.x, by = join_by(Species))
#
# Filter classes and total coverage (averaged in excel), and rename unnamed columns
V_1.x <- V_1 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("V_1_1_pres" = ...3,
         "V_1_2_pres" = ...5,
         "V_1_3_pres" = ...7,
         "V_1_4_pres" = ...9,
         "V_1_5_pres" = ...11,
         "V_1_6_pres" = ...13,
         "V_1_7_pres" = ...15,
         "V_1_8_pres" = ...17,
         "V_1_9_pres" = ...19,
         "V_1_10_pres" = ...21,
         "V_1_11_pres" = ...23,
         "V_1_12_pres" = ...25,
         "V_1_13_pres" = ...27,
         "V_1_14_pres" = ...29,
         "V_1_15_pres" = ...31)
V_2.x <- V_2 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("V_2_1_pres" = ...3,
         "V_2_2_pres" = ...5,
         "V_2_3_pres" = ...7,
         "V_2_4_pres" = ...9,
         "V_2_5_pres" = ...11,
         "V_2_6_pres" = ...13,
         "V_2_7_pres" = ...15,
         "V_2_8_pres" = ...17,
         "V_2_9_pres" = ...19,
         "V_2_10_pres" = ...21,
         "V_2_11_pres" = ...23,
         "V_2_12_pres" = ...25,
         "V_2_13_pres" = ...27,
         "V_2_14_pres" = ...29,
         "V_2_15_pres" = ...31)
V_3.x <- V_3 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("V_3_1_pres" = ...3,
         "V_3_2_pres" = ...5,
         "V_3_3_pres" = ...7,
         "V_3_4_pres" = ...9,
         "V_3_5_pres" = ...11,
         "V_3_6_pres" = ...13,
         "V_3_7_pres" = ...15,
         "V_3_8_pres" = ...17,
         "V_3_9_pres" = ...19,
         "V_3_10_pres" = ...21,
         "V_3_11_pres" = ...23,
         "V_3_12_pres" = ...25,
         "V_3_13_pres" = ...27,
         "V_3_14_pres" = ...29,
         "V_3_15_pres" = ...31)
V_4.x <- V_4 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("V_4_1_pres" = ...3,
         "V_4_2_pres" = ...5,
         "V_4_3_pres" = ...7,
         "V_4_4_pres" = ...9,
         "V_4_5_pres" = ...11,
         "V_4_6_pres" = ...13,
         "V_4_7_pres" = ...15,
         "V_4_8_pres" = ...17,
         "V_4_9_pres" = ...19,
         "V_4_10_pres" = ...21,
         "V_4_11_pres" = ...23,
         "V_4_12_pres" = ...25,
         "V_4_13_pres" = ...27,
         "V_4_14_pres" = ...29,
         "V_4_15_pres" = ...31)
V_5.x <- V_5 %>%
  filter(Species != "Large classes" & Species != "Evergreen shrubs" & Species != "Decideous shrubs" & Species != "Graminoids" & Species != "Forbs" & Species != "Ferns and Equistums") %>%
  select(!Total) %>%
  rename("V_5_1_pres" = ...3,
         "V_5_2_pres" = ...5,
         "V_5_3_pres" = ...7,
         "V_5_4_pres" = ...9,
         "V_5_5_pres" = ...11,
         "V_5_6_pres" = ...13,
         "V_5_7_pres" = ...15,
         "V_5_8_pres" = ...17,
         "V_5_9_pres" = ...19,
         "V_5_10_pres" = ...21,
         "V_5_11_pres" = ...23,
         "V_5_12_pres" = ...25,
         "V_5_13_pres" = ...27,
         "V_5_14_pres" = ...29,
         "V_5_15_pres" = ...31)
#
# Merge all Vassijaure
Vassijaure <- V_1.x %>%
  left_join(V_2.x, by = join_by(Species)) %>%
  left_join(V_3.x, by = join_by(Species)) %>%
  left_join(V_4.x, by = join_by(Species)) %>%
  left_join(V_5.x, by = join_by(Species))
#




#
Abisko_1 <- Abisko %>%
  select(where(~ is.na(.x)))
