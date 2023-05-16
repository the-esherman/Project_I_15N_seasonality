# 15N vegetation and root data
# By Emil A.S. Andersen
# 
#------- ### Libraries ### -------
library(tidyverse)
library(readxl)
library(gridExtra)
library(viridis)
#library(ggpubr)
#library(rstatix)
library(car)
library(pastecs)
library(nlme)
#library(multcomp)
#library(WRS)
#library(ez)
#
#
#
#------- ### Load data ### -------
#
DataName <- "raw_data/15N vegetation and roots v0.35.xlsx"
#
# Biomass, d15N, atom% 15N, and N concentration for enriched and natural abundance samples
vegroot15N <- read_csv("clean_data/Plant_15N_data.csv", col_names = TRUE)
#
# Microbial biomass, d15N, atom% 15N, and recovery ("R_")
Mic15N <- read_csv("clean_data/Mic_15N_data.csv", skip = 1, col_names = TRUE)
#
#
# Define the winter period as snow covered period
winterP <- data.frame(wstart = c(05, 12), wend = c(12, 13))
winterP2 <- data.frame(wstart = c("05_Nov_19", "05_Nov_19"), wend = c("12_May_20", "13_Jun_20"))
#
# List of Measuring periods as they should appear in graphs
measuringPeriod <- c("July-19",	"Aug-19",	"Sep-19",	"Oct-19",	"Nov-19",	"Dec-19",	"Jan-20",	"Feb-20",	"Mar-20",	"Apr-20",	"Apr-20",	"May-20",	"Jun-20",	"Jul-20",	"Aug-20")
#
# Reference atmospheric Nitrogen
# Either 0.003676 or 1/272 (more decimals)
Nair_Rst = 1/272
#
# Added 15N; mg 15N pr patch
N_add <- 1.084
#
# Extraction correction factor
K_EN = 0.4
# See https://climexhandbook.w.uib.no/2019/11/06/soil-microbial-biomass-c-n-and-p/ and UCPH bio lab protocol (where K_EN = 0.4)
#
#
#
#------- ### Main data ### -------
#
# Calculate recovery
#
# Transform numbered months to another format
Month_yr <- tribble(~MP, ~Round,
  01,	"01_Jul_19",
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
# Calculate recovery for plant partition
vegroot15N <- vegroot15N %>%
  mutate(Recovery = ((atom_pc - atom_pc_NatAb)/100 * Nconc/100 * Biomass)/(N_add/1000) * 100)
#
# Vegetation recovery for total core
vegroot15N_total_Plant <- vegroot15N %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
#
#
# Calculate recovery for microbial partition
Mic15N <- Mic15N %>%
  mutate(R_TDN = ((atom_pc_SE - atom_pc_SE_NatAb)/100 * Nconc_SE*10^(-6) * Mic_mass)/(N_add/1000)* 100) %>%
  mutate(R_MBN = (((atom_pc_SEF/100 * Nconc_SEF*10^(-6) - atom_pc_SE/100 * Nconc_SE*10^(-6)) - (atom_pc_SEF_NatAb/100 * Nconc_SEF_NatAb*10^(-6) - atom_pc_SE_NatAb/100 * Nconc_SE_NatAb*10^(-6)))/K_EN * Mic_mass)/(N_add/1000) * 100)
#
# Calculate recovery as a proportion of total recovered in each plot
Rec15N <- vegroot15N_total_Plant %>%
  left_join(Mic15N, by = join_by(Site, Plot, MP, Round)) %>%
  select(Site, Plot, MP, Round, PlantRecovery, R_TDN, R_MBN) %>%
  rowwise() %>%
  mutate(sysRec = sum(PlantRecovery, R_TDN, R_MBN, na.rm = TRUE)) %>%
  mutate(PlantR_frac = PlantRecovery/sysRec*100,
         R_TDN_frac = R_TDN/sysRec*100,
         R_MBN_frac = R_MBN/sysRec*100)
#
#
#
#-------  ### Statistics ### -------
#-------   ##     Q1     ## -------
#
# Model
# Response variable: whole plant recovery of 15N (% of combined functional groups and plant organs)
# Factors: Time, Site, Time*Site
# Most important factor: Time
#
vegroot15N_total_Plant <- vegroot15N_total_Plant %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round"), as.factor))
#
# Contrasts - whole plant recovery
# Month            (J, A, S, O, N, D, J, F, M, A, A, M, J, J, A) # Two times April
SummervsWinter <- c(1, 1, 0, 0,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1)
SpringvsAutumn <- c(0, 0, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1, 0, 0, 0)
SnowvsNot      <- c(-8,-8,-8,-8,7, 7, 7, 7, 7, 7, 7, 7,-8,-8,-8) # Snow from November to May, but June in Vassijaure!
JulvsJan       <- c(1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 1)
OctvsApr       <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0)
Summervs2      <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1) # Summer '19 vs summer '20
SpringChA      <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-2, 0, 0, 0) # Spring change in Abisko: April vs May
SpringChV      <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-2, 0, 0) # Spring change in Vassijaure: April/May vs June
AutumnCh       <- c(0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
WinterCh       <- c(0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0)
cont11         <- c(0, 0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
cont12         <- c(0, 0, 0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0)
cont13         <- c(0, 0, 0, 0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0, 0)
cont14         <- c(1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # Just to get everything to balance out
#
AvsV<-c(1,-1)
contrasts(vegroot15N_total_Plant$Site)<-AvsV
contrasts(vegroot15N_total_Plant$Round)<-cbind(SummervsWinter,SpringvsAutumn,SnowvsNot, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh, cont11, cont12, cont13, cont14)
#contrasts(vegroot15N_total_Plant$Round)<-contr.helmert # Contrasts that compare each new round with the previous ones.
#
# Check if contrasts work, by using a two-way ANOVA
PlantModel_alias <- aov(PlantRecovery ~ Round*Site, data = vegroot15N_total_Plant)
Anova(PlantModel_alias, type ="III")
# alias checks dependencies
alias(PlantModel_alias)
#
# transform data
vegroot15N_total_Plant <- vegroot15N_total_Plant %>%
  mutate(sqrtPlantRecovery = sqrt(PlantRecovery)) %>%
  mutate(invPlantRecovery = 1/PlantRecovery) %>%
  mutate(logPlantRecovery = log(PlantRecovery+1)) %>% # Works the best. Values are small, so even if percent act like log dist.
  mutate(arcPlantRecovery = asin(sqrt(PlantRecovery/100))) # Look into this for general percentages!!

#model:
lme1<-lme(logPlantRecovery ~ Round*Site,
          random = ~1|Plot/Site,
          data = vegroot15N_total_Plant, na.action = na.exclude, method = "REML")

#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1), resid(lme1), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1), main = "Normally distributed?")                 
qqline(resid(lme1), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1)
par(mfrow = c(1,1))

#model output
Anova(lme1, type=2)
summary(lme1)
# Significant for Round
#
#
#
#-------  ### Statistics ### -------
#-------   ##     Q1a    ##  -------
#
# Model
# Response variable: plant recovery of 15N (% of combined functional groups particitioned into organs)
# Factors: Time, Site, Organ, Time*Organ, Time*Site, Site*Organ, Time*Site*Organ
# Most important factors: Organ, Time*Organ, Time*Site*Organ
#
vegroot15N_Organ <- vegroot15N %>%
  select(1:4,6,7,15) %>%
  group_by(across(c("Site","Plot", "MP", "Round", "Organ"))) %>%
  summarise(OrganRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  ungroup() %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "Organ"), as.factor))
#
# Contrasts - plant organs
# For contrasts of Round see Q1
SvsR<-c(-1, -1, 2) # Shoots vs roots
CRvsFR<-c(1,-1,0) # Coarse roots vs fine roots
AvsV<-c(1,-1) # Abisko vs Vassijaure. Maybe not much point in specifying this, but not sure if I dare remove it
contrasts(vegroot15N_Organ$Site)<-AvsV
contrasts(vegroot15N_Organ$Round)<-cbind(SummervsWinter,SpringvsAutumn,SnowvsNot, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh, cont11, cont12, cont13, cont14) # Contrasts defined in Q1
#contrasts(vegroot15N_Organ$Round)<-contr.helmert # Contrasts that compare each new round with the previous ones.
contrasts(vegroot15N_Organ$Organ)<-cbind(SvsR,CRvsFR)
#
# Check if contrasts work, by using a two-way ANOVA
OrganModel_alias <- aov(OrganRecovery ~ Round*Site, data = vegroot15N_Organ)
Anova(OrganModel_alias, type ="III")
# alias checks dependencies
alias(OrganModel_alias)
#
#
#
# transform data
vegroot15N_Organ <- vegroot15N_Organ %>%
  #dplyr::filter(OrganRecovery > 0) %>% 
  mutate(sqrtOrganRecovery = sqrt(OrganRecovery+1)) %>%
  mutate(invOrganRecovery = 1/OrganRecovery+1) %>%
  mutate(logOrganRecovery = log(OrganRecovery+1)) %>%
  mutate(expOrganRecovery = exp(OrganRecovery+1)) %>%
  mutate(cubeOrganRecovery = (OrganRecovery^2)^(1/9)) %>%
  mutate(sqOrganRecovery = OrganRecovery^2) %>%
  mutate(arcOrganRecovery = asin(sqrt((OrganRecovery+1)/100)))
# BoxCox transformation?
#
#
# model:
lme1a<-lme(cubeOrganRecovery ~ Round*Site*Organ,
          random = ~1|Plot/Site,
          data = vegroot15N_Organ, na.action = na.exclude , method = "REML")

#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1a), resid(lme1a), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1a), main = "Normally distributed?")                 
qqline(resid(lme1a), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1a)
par(mfrow = c(1,1))

#model output
Anova(lme1a, type=3)
summary(lme1a)
# Highly significant for Round, Organ Round*Organ and significant for three-way interaction
#
#
#
#-------  ### Statistics ### -------
#-------   ##     Q2     ##  -------
#
# Model
# Response variable: Microbial recovery of 15N (%), MBN
# Factors: Time, Site, Time*Site
# Most important factor: Time
#
# Load data from excel instead of calculated combined
Mic15N_R <- Mic15N %>%
  select("Site", "Plot", "MP", "Round", "R_MBN")
Mic15N_R <- Mic15N_R %>%
  filter(!(is.na(R_MBN))) %>% # remove empty rows, MP9, 13 and 15 as they are missing from either SE or SEF
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round"), as.factor))
#
# Contrasts - MBN recovery
# Month                (J, A, S, O, N, D, J, F, M, A, A, M, J, J, A) # Two times April
# Month                (J, A, S, O, N, D, J, F,  , A, A, M,  , J,  ) # Missing March, June and August-20
# Month                (J, A, S, O, N, D, J, F, A, A, M, J)
SummervsWinter_Mic <- c(1, 1, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1)
SpringvsAutumn_Mic <- c(0, 0, 1, 1, 1, 0, 0, 0,-1,-1,-1, 0)
SnowvsNot_Mic      <- c(7, 7, 7, 7,-5,-5,-5,-5,-5,-5,-5, 7)
Cont4_Mic          <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1) # July vs July
Cont5_Mic          <- c(0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0) # 2 times april
Cont6_Mic          <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-2, 0) # Spring change
Cont7_Mic          <- c(0, 0, 0, 0, 2,-1,-1, 0, 0, 0, 0, 0) # November is highly different from the rest of winter?
Cont8_Mic          <- c(0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0, 0) # Or from autumn?
Cont9_Mic          <- c(1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # Just getting the last to work 
Cont10_Mic         <- c(1, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Cont11_Mic         <- c(1, 1, 1, 1, 1,-5, 0, 0, 0, 0, 0, 0)
AvsV<-c(1,-1)
contrasts(Mic15N_R$Site)<-AvsV
contrasts(Mic15N_R$Round)<-cbind(SummervsWinter_Mic,SpringvsAutumn_Mic,SnowvsNot_Mic, Cont4_Mic, Cont5_Mic, Cont6_Mic, Cont7_Mic, Cont8_Mic, Cont9_Mic, Cont10_Mic, Cont11_Mic) # Contrasts that compare each new round with the previous ones.
#contrasts(Mic15N_R$Round)<-contr.helmert
#
# Alternative: Two-way ANOVA
# have time as a factor in a two-way ANOVA, combined with Site. As each sampling is destructive, the samples are technically independent of each other, although it does not account for the block design
MicModel3 <- aov(R_MBN ~ Round*Site, data = Mic15N_R)
Anova(MicModel3, type ="III")
#
#
# Transform data
Mic15N_R <- Mic15N_R %>%
  mutate(sqrtR_MBN = sqrt(R_MBN+10)) %>%
  mutate(logR_MBN = log(R_MBN+10)) %>%
  mutate(cubeR_MBN = (R_MBN+10)^(1/3)) %>%
  mutate(arcR_MBN = asin(sqrt(((R_MBN+10)/max(R_MBN))/10))) # Standardised as proportion of max microbial recovery
#
#
#model:
lme2<-lme(arcR_MBN ~ Site*Round,
          random = ~1|Plot/Site,
          data = Mic15N_R, na.action = na.exclude , method = "REML")
#
#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme2), resid(lme2), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme2), main = "Normally distributed?")                 
qqline(resid(lme2), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme2)
par(mfrow = c(1,1))

#model output
Anova(lme2, type=2)
summary(lme2)
# Highly significant for Round
#
#
#
#
#
#------- ####   Plotting   ####-------
#-------   ## Plot Biomass ## -------
vegroot15N_bioLong_one <- vegroot15N_bioLong %>%
  group_by(across(c("Site", "Plot", "Round"))) %>%
  summarise(TotalBiomass = sum(Biomass, na.rm = TRUE), .groups = "keep") %>%
  ungroup()

# Plant biomass total
vegroot15N_bioLong_one %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgBiomass = mean(TotalBiomass, na.rm = TRUE), se = sd(TotalBiomass)/sqrt(length(TotalBiomass)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgBiomass, ymin=avgBiomass-se, ymax=avgBiomass+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgBiomass),color = "black") +
  #ylim(0,20) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("Biomass (g pr sample)"), title = expression("Plant biomass")) + #, title = "15N Biomass in plants") + #guides(x = guide_axis(n.dodge = 2)) + 
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Plant biomass by organ
vegroot15N_bioLong %>%
  group_by(across(c("Site", "Plot", "Round", "Part"))) %>%
  summarise(TotalBiomass = sum(Biomass, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  summarise(avgBiomass = mean(TotalBiomass, na.rm = TRUE), se = sd(TotalBiomass)/sqrt(length(TotalBiomass)), .groups = "keep") %>%
  mutate(avgBiomass = if_else(Part == "S", avgBiomass, -avgBiomass),
         se = if_else(Part == "S", se, -se),
         avgR_SE = if_else(Part == "CR", avgBiomass, 0)) %>%
  group_by(across(c("Site", "Round"))) %>%
  mutate(avgR_SE = if_else(Part == "FR", cumsum(avgR_SE)+avgBiomass, avgBiomass)) %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  #
  # Plot 
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgBiomass, fill = factor(Part, levels=c("S","FR","CR"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(-18,3)) +
  scale_fill_viridis_d() +
  #scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Biomass") +
  geom_errorbar(aes(x = Round, y = avgBiomass, ymin=avgR_SE, ymax=avgR_SE+se), position=position_dodge(.9)) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(-18, -15, -12, -9, -6, -3, 0, 3), labels = abs) +
  #scale_fill_discrete(labels = c("Shoots", "Fine Roots", "Course roots")) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("Biomass (g pr sample)"), title = expression("Plant biomass")) + #guides(x = guide_axis(n.dodge = 2)) + 
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
#
#-------   ##    Nconc     ## -------
#
Nconc_one <- Nconc %>%
  dplyr::select(Site, Plot, MP, Round, 6:8) %>%
  dplyr::rename("Plant_S" = "Plant_S_N",
                "Plant_CR" = "Plant_CR_N",
                "Plant_FR" = "Plant_FR_N"
  ) %>%
  pivot_longer(cols = 5:7, names_to = "Type", values_to = "Nconc")
Nconc_one <- Nconc_one %>%
  add_column(Part = str_split_fixed(Nconc_one$Type,"\\w+_",n=2)[,2]) %>%
  add_column(Species = str_split_fixed(Nconc_one$Type,"_\\w+",n=2)[,1])
#
#
# Total N
Nconc %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgNconc = mean(Plant_total_N, na.rm = TRUE), se = sd(Plant_total_N)/sqrt(length(Plant_total_N)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgNconc, ymin=avgNconc-se, ymax=avgNconc+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgNconc),color = "black") +
  #ylim(0,20) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("[N] % (g "*g^-1*" DW)"), title = expression("Plant N concentration")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Plant N by organ
Nconc_one %>%
  group_by(across(c("Site", "Plot", "Round", "Part"))) %>%
  summarise(TotalNconc = sum(Nconc, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  summarise(avgNconc = mean(TotalNconc, na.rm = TRUE), se = sd(TotalNconc)/sqrt(length(TotalNconc)), .groups = "keep") %>%
  mutate(avgNconc = if_else(Part == "S", avgNconc, -avgNconc),
         se = if_else(Part == "S", se, -se),
         avgR_SE = if_else(Part == "CR", avgNconc, 0)) %>%
  group_by(across(c("Site", "Round"))) %>%
  mutate(avgR_SE = if_else(Part == "FR", cumsum(avgR_SE)+avgNconc, avgNconc)) %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  #
  # Plot 
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgNconc, fill = factor(Part, levels=c("S","FR","CR"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(-1.5,1.5)) +
  scale_fill_viridis_d() +
  #scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Nconc") +
  geom_errorbar(aes(x = Round, y = avgNconc, ymin=avgR_SE, ymax=avgR_SE+se), position=position_dodge(.9)) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), labels = abs) +
  #scale_fill_discrete(labels = c("Shoots", "Fine Roots", "Course roots")) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("[N] % (g "*g^-1*" DW)"), title = expression("Plant N concentration")) +
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# For each species
# (Change the filter of species and the plot title)
vegroot15N_NLong1 %>%
  dplyr::filter(Species == "ES") %>%
  group_by(across(c("Site", "Plot", "Round", "Part"))) %>%
  summarise(TotalNconc = sum(Nconc, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  summarise(avgNconc = mean(TotalNconc, na.rm = TRUE), se = sd(TotalNconc)/sqrt(length(TotalNconc)), .groups = "keep") %>%
  mutate(avgNconc = if_else(Part == "S", avgNconc, -avgNconc),
         se = if_else(Part == "S", se, -se),
         avgR_SE = if_else(Part == "CR", avgNconc, 0)) %>%
  group_by(across(c("Site", "Round"))) %>%
  mutate(avgR_SE = if_else(Part == "FR", cumsum(avgR_SE)+avgNconc, avgNconc)) %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  #
  # Plot 
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgNconc, fill = factor(Part, levels=c("S","FR","CR"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(-2,2)) +
  scale_fill_viridis_d() +
  geom_errorbar(aes(x = Round, y = avgNconc, ymin=avgR_SE, ymax=avgR_SE+se), position=position_dodge(.9)) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), labels = abs) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("[N] % (g "*g^-1*" DW)"), 
       title = expression("Evergreen N concentration")) + # <--------- Title
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Inorganic N content
# NO3 in SE
inorgN %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgN = mean(NO3_SE, na.rm = TRUE), se = sd(NO3_SE)/sqrt(length(NO3_SE)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgN, ymin=avgN-se, ymax=avgN+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgN),color = "black") +
  coord_cartesian(ylim = c(0,250)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("["*NO[3]*"] (µg "*L^-1*")"), 
       title = expression(NO[3]*" concentration in extract, non-fumigated")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
# NO3 in SEF
inorgN %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgN = mean(NO3_SEF, na.rm = TRUE), se = sd(NO3_SEF)/sqrt(length(NO3_SEF)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgN, ymin=avgN-se, ymax=avgN+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgN),color = "black") +
  coord_cartesian(ylim = c(0,250)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("["*NO[3]*"] (µg "*L^-1*")"), 
       title = expression(NO[3]*" concentration in extract, fumigated")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
# NH4 in SE
inorgN %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgN = mean(NH4_SE, na.rm = TRUE), se = sd(NH4_SE)/sqrt(length(NH4_SE)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgN, ymin=avgN-se, ymax=avgN+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgN),color = "black") +
  coord_cartesian(ylim = c(0,1000)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("["*NH[4]*"] (µg "*L^-1*")"), 
       title = expression(NH[4]*" concentration in extract, non-fumigated")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
# NH4 in SEF
inorgN %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgN = mean(NH4_SEF, na.rm = TRUE), se = sd(NH4_SEF)/sqrt(length(NH4_SEF)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgN, ymin=avgN-se, ymax=avgN+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgN),color = "black") +
  coord_cartesian(ylim = c(0,6000)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("["*NH[4]*"] (µg "*L^-1*")"), 
       title = expression(NH[4]*" concentration in extract, fumigated")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
#-------   ##     d15N     ##-------
#

#-------   ##   Recovery   ## -------
#
# Sum recovery and calculate average
# This means combining
# Shoots: ES, DS, G, O, U
# CR: ES, DS, G, bulk
# FR: ES, DS, G, O, bulk
vegroot15N %>%
  group_by(across(c("Site", "Plot", "Round", "Organ"))) %>%
  summarise(TotalRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  summarise(avgRecovery = mean(TotalRecovery, na.rm = TRUE), se = sd(TotalRecovery)/sqrt(length(TotalRecovery)), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Organ == "S", avgRecovery, -avgRecovery),
         se = if_else(Organ == "S", se, -se),
         avgR_SE = if_else(Organ == "CR", avgRecovery, 0)) %>%
  group_by(across(c("Site", "Round"))) %>%
  mutate(avgR_SE = if_else(Organ == "FR", cumsum(avgR_SE)+avgRecovery, avgRecovery)) %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  #
  # Plot 
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery, fill = factor(Organ, levels=c("S","FR","CR"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(-15,3)) +
  scale_fill_viridis_d() +
#  scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgR_SE, ymax=avgR_SE+se), position=position_dodge(.9)) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(-15, -12, -9, -6, -3, 0, 3), labels = abs) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) +
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Each species or part separated. For a quick overview of where patterns might come from
vegroot15N %>%
  group_by(across(c("Site", "Plot", "Round", "Organ", "Species"))) %>%
  summarise(TotalRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Organ", "Species"))) %>%
  summarise(avgRecovery = mean(TotalRecovery, na.rm = TRUE), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Organ == "S", avgRecovery, -avgRecovery)) %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  #
  # Plot 
  ggplot() +
  #geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery, fill = factor(Organ, levels=c("S","FR","CR"))), position = "stack") +
  coord_cartesian(ylim = c(-8,3)) +
  scale_fill_viridis_d() +
  #scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") +
  #geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgR_SE-se, ymax=avgR_SE+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site + Species, ncol = 5) +
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) +
  guides(fill = guide_legend(title = "Plant organ")) + 
  theme_light(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1)) 
#
#
#
# Abisko and Vassijaure plant recovery faceted
vegroot15N_total_Plant %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(PlantRecovery, na.rm = TRUE), se = sd(PlantRecovery)/sqrt(length(PlantRecovery)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgRecovery),color = "black") +
  ylim(0,20) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) + #, title = "15N recovery in plants") + #guides(x = guide_axis(n.dodge = 2)) + 
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Proportional to total recovery
# Plant recover fraction
Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(PlantR_frac, na.rm = TRUE), se = sd(PlantR_frac)/sqrt(length(PlantR_frac)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery, ymax=avgRecovery+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgRecovery)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period", y = expression("% of total recovered "*{}^15*"N"), title = expression({}^15*"N recovery in plants, proportional to total recovery")) +
  theme_classic(base_size = 20)  +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
# Mic - MBN fraction
Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(R_MBN_frac, na.rm = TRUE), se = sd(R_MBN_frac)/sqrt(length(R_MBN_frac)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery)) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period", y = expression("% of total recovered "*{}^15*"N"), title = expression({}^15*"N recovery in plants, proportional to total recovery")) +
  theme_classic(base_size = 20)  +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1)) 
#
#
# Microbial and soil part
# MBN
Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgR = mean(R_MBN, na.rm = TRUE), 
            se = sd(R_MBN)/sqrt(length(R_MBN)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgR, ymin=avgR-se, ymax=avgR+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgR),color = "black") +
  coord_cartesian(ylim = c(0,120)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("% of added "*{}^15*"N"), 
       title = expression("Microbial "*{}^15*"N tracer recovery")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
# Soil
Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgR = mean(R_TDN, na.rm = TRUE), 
            se = sd(R_TDN)/sqrt(length(R_TDN)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgR, ymin=avgR-se, ymax=avgR+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgR),color = "black") +
  coord_cartesian(ylim = c(0,0.8)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("% of added "*{}^15*"N"), 
       title = expression("Soil "*{}^15*"N tracer recovery")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# System recovery
Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(sysRec, na.rm = TRUE), se = sd(sysRec)/sqrt(length(sysRec)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery)) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("% of added "*{}^15*"N"), 
       title = expression("Total "*{}^15*"N tracer recovery")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#


margin <- qt()


#
#
#
#------- ### Checking data ### -------
#
library("outliers") 
library("VIM") # For matrixplot
#
# Outliers
#
hist(log(Rec15N$TotalRecovery+1), main = "Histogram")
# Cleveland plot
dotchart(log(Rec15N$TotalRecovery+1), 
         main="Cleveland plot", xlab = "Observed values", 
         pch = 19, color = hcl.colors(12), 
         labels = Rec15N$Round, 
         groups = Rec15N$Plot,
         gdata = tapply(Rec15N$TotalRecovery, Rec15N$Plot, mean),
         gpch = 12, gcolor = 1)
#
# Missing values
aggr(vegroot15N[56:70])
matrixplot(vegroot15N[56:70])
aggr(Rec15N)
matrixplot(Rec15N)
#
# Outliers
# Grubb's test for single outliers
grubbs.test(log(Rec15N$TotalRecovery+1))
grubbs.test(log(Rec15N$TotalRecovery+1), opposite = TRUE)
# Data log transformed to meet normality
# Both highest and lowest values are outliers
#
# Rosner's test for multiple (k) outliers (generalized ESD many-outliers)
EnvStats::rosnerTest(log(Rec15N$TotalRecovery+1), k = 2)$all.stats
# No outliers!
#
# Duplicate IRMS numbers
IRMS_dupli <- IRMS %>%
  group_by(ID) %>%
  filter(n()>1)

#
#
#------- The End -------