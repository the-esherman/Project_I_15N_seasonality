# 15N vegetation and root data
# By Emil A.S. Andersen
# 
#------- ### Libraries ### -------
library(tidyverse)
library(readxl)
#library(gridExtra)
#library(viridis)
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
winterP2 <- data.frame(wstart = c("05_Nov-19", "05_Nov-19"), wend = c("12_May-20", "13_Jun-20"))
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
  group_by(across(c("Site", "Plot", "Round"))) %>%
  summarise(PlantRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
#
#
# Calculate recovery for microbial partition
Mic15N <- Mic15N %>%
  mutate(R_TDN = ((atom_pc_SE - atom_pc_SE_NatAb)/100 * Nconc_SE*10^(-6) * Mic_mass)/(N_add/1000)* 100) %>%
  mutate(R_MBN = (((atom_pc_SEF/100 * Nconc_SEF*10^(-6) - atom_pc_SE/100 * Nconc_SE*10^(-6)) - (atom_pc_SEF_NatAb/100 * Nconc_SEF_NatAb*10^(-6) - atom_pc_SE_NatAb/100 * Nconc_SE_NatAb*10^(-6)))/K_EN * Mic_mass)/(N_add/1000) * 100)
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
  left_join(Month_yr, by = join_by("Round")) %>%
  relocate("MP", .after = "Plot") %>%
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
vegroot15N_RLong_Organ <- vegroot15N_RLong_Organ_original %>%
  mutate(across(c("Plot"), as.character))%>%
  mutate(across(c("Site", "Round", "Organ"), as.factor))
#
# Contrasts - plant organs
SvsR<-c(-1, -1, 2) # Shoots vs roots
CRvsFR<-c(1,-1,0) # Coarse roots vs fine roots
AvsV<-c(1,-1) # Abisko vs Vassijaure. Maybe not much point in specifying this, but not sure if I dare remove it
contrasts(vegroot15N_RLong_Organ$Site)<-AvsV
contrasts(vegroot15N_RLong_Organ$Round)<-cbind(SummervsWinter,SpringvsAutumn,SnowvsNot, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh, cont11, cont12, cont13, cont14) # Contrasts defined in Q1
#contrasts(vegroot15N_RLong_Organ$Round)<-contr.helmert # Contrasts that compare each new round with the previous ones.
contrasts(vegroot15N_RLong_Organ$Organ)<-cbind(SvsR,CRvsFR)
#
# Check if contrasts work, by using a two-way ANOVA
OrganModel_alias <- aov(OrganRecovery ~ Round*Site, data = vegroot15N_RLong_Organ)
Anova(OrganModel_alias, type ="II")
# alias checks dependencies
alias(OrganModel_alias)
#
#
#
# Still not working. See comment above
organModel <-ezANOVA(data = vegroot15N_RLong_Organ, dv = .(OrganRecovery), wid = .(Plot), within = .(Organ, Round, Site), type = 3, detailed = TRUE)
#
# Test Homogeneity of variance
leveneTest(vegroot15N_RLong_Organ$OrganRecovery, vegroot15N_RLong_Organ$Organ, center = median) # Organ
leveneTest(vegroot15N_RLong_Organ$OrganRecovery, vegroot15N_RLong_Organ$Round, center = median) # Round
leveneTest(vegroot15N_RLong_Organ$OrganRecovery, vegroot15N_RLong_Organ$Site, center = median) # Site
leveneTest(vegroot15N_RLong_Organ$OrganRecovery, interaction(vegroot15N_RLong_Organ$Round, vegroot15N_RLong_Organ$Site), center = median) # interaction Round*Site
# Very much not equal variance, except for site
#
# transform data
vegroot15N_RLong_Organ <- vegroot15N_RLong_Organ %>%
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
# Check distribution
vegroot15N_RLong_Organ %>% 
  ggplot(aes(OrganRecovery), color = Site) + geom_histogram()
vegroot15N_RLong_Organ %>% 
  ggplot(aes(sqOrganRecovery), color = Site) + geom_histogram()
vegroot15N_RLong_Organ %>% 
  ggplot(aes(cubeOrganRecovery), color = Site) + geom_histogram()
hist(vegroot15N_RLong_Organ$cubeOrganRecovery, main = "Historgram of plant recovery")
hist(vegroot15N_RLong_Organ$arcOrganRecovery, main = "Historgram of plant recovery")
#
# Shapiro-Wilk Test
shapiro.test(vegroot15N_RLong_Organ$cubeOrganRecovery)
#
# Check qq-plot of transformations
# log transformation best
qqnorm(vegroot15N_RLong_Organ$OrganRecovery, main = "Normal Q-Q Plot")
qqline(vegroot15N_RLong_Organ$OrganRecovery)
qqnorm(vegroot15N_RLong_Organ$logOrganRecovery, main = "Log Q-Q Plot")
qqline(vegroot15N_RLong_Organ$logOrganRecovery)
qqnorm(vegroot15N_RLong_Organ$invOrganRecovery, main = "Inverse Q-Q Plot")
qqline(vegroot15N_RLong_Organ$invOrganRecovery)
qqnorm(vegroot15N_RLong_Organ$sqrtOrganRecovery, main = "sqrt Q-Q Plot")
qqline(vegroot15N_RLong_Organ$sqrtOrganRecovery)
qqnorm(vegroot15N_RLong_Organ$expOrganRecovery, main = "exponential Q-Q Plot")
qqline(vegroot15N_RLong_Organ$expOrganRecovery)
qqnorm(vegroot15N_RLong_Organ$cubeOrganRecovery, main = "Cube root Q-Q Plot")
qqline(vegroot15N_RLong_Organ$cubeOrganRecovery)
qqnorm(vegroot15N_RLong_Organ$arcOrganRecovery, main = "Angular transformation (arcsine of sqrt) Q-Q Plot")
qqline(vegroot15N_RLong_Organ$arcOrganRecovery)
#
# Check homogeneity of variance
# log is homogeneous
leveneTest(vegroot15N_RLong_Organ$cubeOrganRecovery, vegroot15N_RLong_Organ$Organ, center = median) # Organ
leveneTest(vegroot15N_RLong_Organ$cubeOrganRecovery, vegroot15N_RLong_Organ$Round, center = median) # Round
leveneTest(vegroot15N_RLong_Organ$cubeOrganRecovery, vegroot15N_RLong_Organ$Site, center = median) # Site
leveneTest(vegroot15N_RLong_Organ$cubeOrganRecovery, interaction(vegroot15N_RLong_Organ$Round, vegroot15N_RLong_Organ$Site), center = median) 
leveneTest(vegroot15N_RLong_Organ$cubeOrganRecovery, interaction(vegroot15N_RLong_Organ$Organ, vegroot15N_RLong_Organ$Site), center = median)
leveneTest(vegroot15N_RLong_Organ$cubeOrganRecovery, interaction(vegroot15N_RLong_Organ$Round, vegroot15N_RLong_Organ$Organ), center = median)
# cube root transformation not equal variance for Round or interaction
#
# Cube root is approx. normal-distributed, but variance not equal
#
# Two-way ANOVA
PlantModel3 <- aov(logOrganRecovery ~ Round*Site*Organ, data = vegroot15N_RLong_Organ)
Anova(PlantModel3, type ="III")

# Multilevel linear model approach
baseline_organ <- lme(logOrganRecovery ~ 1, random = ~1|Plot/Round/Site, data = vegroot15N_RLong_Organ, method = "ML")
organOrganModel <- update(baseline_organ, .~. + Organ)
organRoundModel <- update(organOrganModel, .~. + Round)
organSiteModel <- update(organRoundModel, .~. + Site)

organIntactModel <- update(organSiteModel, .~. + Round:Site)


#
anova(baseline_organ, organRoundModel, organSiteModel, organIntactModel)
anova(baseline_organ, organSiteModel, organRoundModel, organIntactModel)
#
summary(organIntactModel)

vegroot15N_RLong_Organ %>%
  ggplot(aes(x = reorder(Round, sort(as.numeric(Round))), OrganRecovery, colour = Site)) + stat_summary(fun.y = mean, geom = "point") + stat_summary(fun.y = mean, geom = "line", aes(group= Site)) + stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + labs(x = "Round", y = "Mean Recovery", colour = "Site")

#
#
#summary(aov(logPlantRecovery ~ Round*Site + Error(Plot/Round/Site), data = vegroot15N_total_Plant))
#
# Posthoc test
PostOrgan <- glht(organIntactModel, linfct = mcp(Round = "Tukey"))
summary(PostOrgan)
confint(PostOrgan)
#
OrganIntactModel_post <- lme(logOrganRecovery ~ 1 + Site_Round, random = ~1|Plot/Round/Site, data = vegroot15N_RLong_Organ, method = "ML")
PostOrgan2 <- glht(OrganIntactModel_post, linfct = mcp(Site_Round = "Tukey"))
summary(PostOrgan2)
confint(PostOrgan2)
#
#
#
# From Signe
#model:
lme1a<-lme(cubeOrganRecovery ~ Round*Site*Organ,
          random = ~1|Site/Plot,
          data = vegroot15N_RLong_Organ, na.action = na.exclude , method = "REML")

#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1a), resid(lme1a), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1a), main = "Normally distributed?")                 
qqline(resid(lme1a), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1a)
par(mfrow = c(1,1))

#model output
Anova(lme1a, type=2)
summary(lme1a)
# Highly significant for Round, Organm Round*Organ and significant for three-way interaction
#
#
#
# From: https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
# Check data
# Summary statistics
vegroot15N_RLong_Organ %>%
  group_by(Site, Organ, Round) %>%
  get_summary_stats(cubeOrganRecovery, type = "mean_sd")
#
#
OrganBoxP <- ggboxplot(vegroot15N_RLong_Organ, x = "Round", y = "OrganRecovery", color = "Site", palette = "jco", facet.by = "Organ", short.panel.labs = FALSE) + guides(x = guide_axis(n.dodge = 2))
OrganBoxP
#
# Identify outliers
vegroot15N_RLong_Organ %>%
  group_by(Site, Round, Organ) %>%
  identify_outliers(cubeOrganRecovery)
# Several extreme outliers
#
# Normality
print(
  vegroot15N_RLong_Organ %>%
    group_by(Site, Round, Organ) %>%
    shapiro_test(cubeOrganRecovery), n = 90)
ggqqplot(vegroot15N_RLong_Organ, "cubeOrganRecovery", ggtheme = theme_bw()) + facet_grid(Site + Round ~ Organ, labeller = "label_both")
# Not normally distributed at each Round-Site: at least 3 significantly different, log-transformed: 2, arcsin: 1 (always Vassi MP15)
# 
#
# ANOVA
organ_aov <- anova_test(data = vegroot15N_RLong_Organ, dv = "cubeOrganRecovery", wid = "Plot", within = c("Site", "Round", "Organ"))
organ_aov
get_anova_table(organ_aov)
#
# Post-hoc tests
# Two-way ANOVA at each diet level
organTwoWay <- vegroot15N_RLong_Organ %>%
  group_by(Organ) %>%
  anova_test(dv = cubeOrganRecovery, wid = Plot, within = c(Site, Round)) %>%
  adjust_pvalue(method = "bonferroni")
organTwoWay
# Extract anova table
get_anova_table(organTwoWay)
# One-way ANOVA
organOneWay <- vegroot15N_RLong_Organ %>%
  group_by(Organ, Site) %>%
  anova_test(dv = cubeOrganRecovery, wid = Plot, within = Round) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
organOneWay
#
# Pair-wise comparison: paired t-test # Unbalanced setup as missing some organs at some time-points
organPWC <- vegroot15N_RLong_Organ %>%
  group_by(Organ, Round) %>%
  pairwise_t_test(cubeOrganRecovery ~ Site, 
    paired = TRUE, p.adjust.method = "bonferroni"
  )
organPWC
#
# Visualization: box plots with p-values
organPWC_plot <- organPWC %>% add_xy_position(x = "Round", y.trans = function(x){(x)^5})
OrganBoxP + 
  stat_pvalue_manual(organPWC_plot, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(organ_aov, detailed = TRUE),
    caption = get_pwc_label(organPWC_plot)
  )


#
# Different plot # not working :(
#                                                     #1                                    #2                                   #3
OrganBoxP2 <- ggboxplot(vegroot15N_RLong_Organ, x = "Organ", y = "OrganRecovery", color = "Round", palette = "jco", facet.by = "Site", short.panel.labs = FALSE) + guides(x = guide_axis(n.dodge = 2))
OrganBoxP2
# Pair-wise
organPWC <- vegroot15N_RLong_Organ %>%
  #         #3     #1
  group_by(Site, Organ) %>% #          #2
  pairwise_t_test(cubeOrganRecovery ~ Round, 
                  paired = TRUE, p.adjust.method = "bonferroni"
  )
# Visualization: box plots with p-values
organPWC_plot2 <- organPWC %>% add_xy_position(x = "exercises", y.trans = function(x){(x)^5})
OrganBoxP2 + 
  stat_pvalue_manual(organPWC_plot2, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(organ_aov, detailed = TRUE),
    caption = get_pwc_label(organPWC_plot2)
  )
#


#
#
#
#
#
#
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
Mic15N_RLong <- Mic15N %>%
  dplyr::select("Site", "Plot", "MP", "Round", "R_MBN")
Mic15N_RLong <- Mic15N_RLong %>%
  dplyr::filter(!(is.na(R_MBN))) %>% # remove empty rows, MP9, 13 and 15 as they are missing from either SE or SEF
  mutate(across(c("Plot", "Round"), as.character))%>%
  mutate(across(c("Site", "Round"), as.factor))
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
contrasts(Mic15N_RLong$Site)<-AvsV
contrasts(Mic15N_RLong$Round)<-cbind(SummervsWinter_Mic,SpringvsAutumn_Mic,SnowvsNot_Mic, Cont4_Mic, Cont5_Mic, Cont6_Mic, Cont7_Mic, Cont8_Mic, Cont9_Mic, Cont10_Mic, Cont11_Mic) # Contrasts that compare each new round with the previous ones.
#contrasts(Mic15N_RLong$Round)<-contr.helmert
#
# Check for missing values in data set. See esp. Plot vs Round vs Site
ezDesign(Mic15N_RLong, x = Round, y = Plot)
ezDesign(Mic15N_RLong, x = Site, y = Plot)
ezDesign(Mic15N_RLong, x = Round, y = Site)
ezDesign(Mic15N_RLong, x = Plot, y = R_MBN)
ezDesign(Mic15N_RLong, x = Site, y = R_MBN)
#
# Still not working. See comment above
Mic_Model <-ezANOVA(data = Mic15N_RLong, dv = .(sqrtR_MBN), wid = .(Plot), within = .(Round, Site), type = 3, detailed = TRUE)
#
# Alternative: Two-way ANOVA
# have time as a factor in a two-way ANOVA, combined with Site. As each sampling is destructive, the samples are technically independent of each other, although it does not account for the block design
MicModel3 <- aov(R_MBN ~ Round*Site, data = Mic15N_RLong)
Anova(MicModel3, type ="III")

qqplot(PlantModel3, main = "Normal Q-Q Plot")

# Test Homogeneity of variance
leveneTest(Mic15N_RLong$R_MBN, Mic15N_RLong$Round, center = median)
leveneTest(Mic15N_RLong$R_MBN, Mic15N_RLong$Site, center = median)
leveneTest(Mic15N_RLong$R_MBN, interaction(Mic15N_RLong$Round, Mic15N_RLong$Site), center = median)
# No transformation necessary
#
# Check distribution
Mic15N_RLong %>% 
  ggplot(aes(R_MBN), color = Site) + geom_histogram()
qqnorm(Mic15N_RLong$R_MBN, main = "Normal Q-Q Plot")
qqline(Mic15N_RLong$R_MBN)
shapiro.test(Mic15N_RLong$R_MBN)
# Not normally distributed, but close
#
# transform data
Mic15N_RLong <- Mic15N_RLong %>%
  mutate(sqrtR_MBN = sqrt(R_MBN+10)) %>%
  mutate(logR_MBN = log(R_MBN+10)) %>%
  mutate(cubeR_MBN = (R_MBN+10)^(1/3)) %>%
  mutate(arcR_MBN = asin(sqrt(((R_MBN+10)/max(R_MBN))/100)))
  
#
# Check distribution of transformed data
Mic15N_RLong %>% 
  ggplot(aes(sqrtR_MBN), color = Site) + geom_histogram()
hist(Mic15N_RLong$sqrtR_MBN)
Mic15N_RLong %>% 
  ggplot(aes(logR_MBN), color = Site) + geom_histogram()
hist(Mic15N_RLong$cubeR_MBN)
hist(Mic15N_RLong$arcR_MBN)
#
# Check qq-plot of transformations
# sqrt transformation best
qqnorm(Mic15N_RLong$logR_MBN, main = "Log Q-Q Plot")
qqline(Mic15N_RLong$logR_MBN)
qqnorm(Mic15N_RLong$sqrtR_MBN, main = "sqrt Q-Q Plot")
qqline(Mic15N_RLong$sqrtR_MBN)
qqnorm(Mic15N_RLong$cubeR_MBN, main = "cube Q-Q Plot")
qqline(Mic15N_RLong$cubeR_MBN)
qqnorm(Mic15N_RLong$arcR_MBN, main = "angular transformation Q-Q Plot")
qqline(Mic15N_RLong$arcR_MBN)
#
# Shapiro-Wilk test for normal distribution
shapiro.test(Mic15N_RLong$sqrtR_MBN)
# sqrt transformation is normally distributed
#
# Check homogeneity of variance
leveneTest(Mic15N_RLong$sqrtR_MBN, Mic15N_RLong$Round, center = median)
leveneTest(Mic15N_RLong$sqrtR_MBN, Mic15N_RLong$Site, center = median)
leveneTest(Mic15N_RLong$sqrtR_MBN, interaction(Mic15N_RLong$Round, Mic15N_RLong$Site), center = median)
# sqrt is homogeneous
#
# sqrt is normal-distributed, and variance is equal
#
# Two-way ANOVA
MicModel3 <- aov(sqrtR_MBN ~ Round*Site, data = Mic15N_RLong)
Anova(MicModel3, type ="II")

alias(MicModel3)
# contrasts have no conflicts

# Multilevel linear model approach
baseline_Mic <- lme(sqrtR_MBN ~ 1, random = ~1|Plot/Round/Site, data = Mic15N_RLong, method = "ML")
micRoundModel <- update(baseline_Mic, .~. + Round)
micSiteModel <- update(micRoundModel, .~. + Site)
micIntactModel <- update(micSiteModel, .~. + Round:Site)
#
anova(baseline_Mic, micRoundModel, micSiteModel, micIntactModel)
anova(baseline_Mic, micSiteModel, micRoundModel, micIntactModel)
#

summary(micIntactModel)

Mic15N_RLong %>%
  ggplot(aes(x = reorder(Round, sort(as.numeric(Round))), R_MBN, colour = Site)) + stat_summary(fun.y = mean, geom = "point") + stat_summary(fun.y = mean, geom = "line", aes(group= Site)) + stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + labs(x = "Round", y = "Mean Recovery", colour = "Site")

#
#
#summary(aov(logR_MBN ~ Round*Site + Error(Plot/Round/Site), data = Mic15N_RLong))
#
# Posthoc test
PostPlant <- glht(plantIntactModel, linfct = mcp(Round = "Tukey"))
summary(PostPlant)
confint(PostPlant)
#
plantIntactModel_post <- lme(logR_MBN ~ 1 + Site_Round, random = ~1|Plot/Round/Site, data = Mic15N_RLong, method = "ML")
PostPlant2 <- glht(plantIntactModel_post, linfct = mcp(Site_Round = "Tukey"))
summary(PostPlant2)
confint(PostPlant2)
#

by(Mic15N_RLong$logR_MBN, list(Mic15N_RLong$Round, Mic15N_RLong$Site), stat.desc)
#
#
#
# From Signe
#model:
lme2<-lme(sqrtR_MBN ~ Site*Round,
          random = ~1|Plot/Site,
          data = Mic15N_RLong, na.action = na.exclude , method = "REML")

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
# To get a idea of nested vs crossed design:
# https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified
#
#
# From: https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
# Check data
# Summary statistics
Mic15N_RLong %>%
  group_by(Round, Site) %>%
  get_summary_stats(R_MBN, type = "mean_sd")
#
#
micBoxP <- ggboxplot(Mic15N_RLong, x = "Round", y = "R_MBN", color = "Site", palette = "jco")
micBoxP
micBoxP_sqrt <- ggboxplot(Mic15N_RLong, x = "Round", y = "sqrtR_MBN", color = "Site", palette = "jco")
micBoxP_sqrt
#
# Identify outliers
Mic15N_RLong %>%
  group_by(Round, Site) %>%
  identify_outliers(sqrtR_MBN)
# Several extreme outliers
#
# Normality
print(
  Mic15N_RLong %>%
    group_by(Round, Site) %>%
    shapiro_test(sqrtR_MBN), n = 30)
ggqqplot(Mic15N_RLong, "sqrtR_MBN", ggtheme = theme_bw()) + facet_grid(Round ~ Site, labeller = "label_both")
# Not normally distributed at each Round-Site
# 
#
# ANOVA
mic_aov <- anova_test(data = Mic15N_RLong, dv = "sqrtR_MBN", wid = "Plot", within = c("Site", "Round"))
mic_aov
get_anova_table(mic_aov)
#
# Post-hoc tests
# One-way ANOVA
micOneWay <- Mic15N_RLong %>%
  group_by(Round) %>%
  anova_test(dv = sqrtR_MBN, wid = Plot, within = Site) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
micOneWay
# Pair-wise comparison: paired t-test
micPWC <- Mic15N_RLong %>%
  group_by(Round) %>%
  pairwise_t_test(
    sqrtR_MBN ~ Site, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
micPWC
#
# Similarly for Round
micOneWay2 <- Mic15N_RLong %>%
  group_by(Site) %>%
  anova_test(dv = sqrtR_MBN, wid = Plot, within = Round) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
micOneWay2
# Pair-wise comparison: paired t-test
micPWC2 <- Mic15N_RLong %>%
  group_by(Site) %>%
  pairwise_t_test(
    sqrtR_MBN ~ Round, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
micPWC2
#
# When no interaction effect
# comparisons for Site variable
Mic15N_RLong %>%
  pairwise_t_test(
    sqrtR_MBN ~ Site, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )
# comparisons for Round variable
print(
  Mic15N_RLong %>%
    pairwise_t_test(
      sqrtR_MBN ~ Round, paired = TRUE, 
      p.adjust.method = "bonferroni"
    ), n = 105)
# Visualization: box plots with p-values
micPWC_plot <- micPWC %>% add_xy_position(x = "Round", y.trans = function(x){x^2-10})
micBoxP + 
  stat_pvalue_manual(micPWC_plot, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(mic_aov, detailed = TRUE),
    caption = get_pwc_label(micPWC_plot)
  )

#
#
#
#
# # Trying a few different things
# For this method check link:
# https://stats.stackexchange.com/questions/58435/repeated-measures-error-in-r-ezanova-using-more-levels-than-subjects-balanced-d
set.seed(123)  ## make reproducible
N  <- 5 #18       ## number of subjects
P  <- 2 #3        ## number of conditions
Q  <- 12 #29       ## number of sites
voltage <- matrix(round(rnorm(N*P*Q), 2), nrow=N)   ## (N x (PxQ))-matrix with voltages
fit  <- lm(voltage ~ 1)   ## between-subjects design (here: no between factors)
inDf <- expand.grid(channel=gl(P, 1), electrode=gl(Q, 1))  ## within design
library(car)              ## for Anova()
AnRes <- Anova(fit, idata=inDf, idesign=~channel*electrode)
summary(AnRes, multivariate=FALSE, univariate=TRUE)
# With sphericity test
anova(fit, M=~channel, X=~1, idata=inDf, test="Spherical")
anova(fit, M=~channel + electrode, X=~channel, idata=inDf, test="Spherical")
anova(fit, M=~channel + electrode + channel:electrode, X=~channel + electrode, idata=inDf, test="Spherical")
#
#
#
Mic15N_RLong_wideData <- Mic15N_RLong %>%
  dplyr::select("Site", "Plot", "Round", "R_MBN") %>%
  mutate(Site_Round = str_c(Site, Round, sep = "_")) %>%
  dplyr::select(!c("Site", "Round")) %>%
  pivot_wider(names_from = "Site_Round", values_from = "R_MBN")
Mic15N_RLong_wideData1 <- Mic15N_RLong_wideData %>%
  dplyr::select(!c("Plot"))
Mic15N_RLong_wideData2 <- Mic15N_RLong %>%
  dplyr::select("Site", "Round") %>%
  mutate(across(c("Site", "Round"), as.factor))
#
fitMic <- lm(Mic15N_RLong_wideData1 ~ 1)

#


#
# Trying the way described by Peter Dalgaard here:
# https://cran.r-project.org/doc/Rnews/Rnews_2007-2.pdf
MicModel4_base <- lm(R_MBN ~ Round, data = Mic15N_RLong)
mauchly.test(MicModel4_base, X=~1)
# But not working
#
reacttime <- matrix(c(420, 420, 480, 480, 600, 780, 420, 480, 480, 360, 480, 600, 480, 480, 540, 660, 780, 780, 420, 540, 540, 480, 780, 900, 540, 660, 540, 480, 660, 720, 360, 420, 360, 360, 480, 540, 480, 480, 600, 540, 720, 840, 480, 600, 660, 540, 720, 900, 540, 600, 540, 480, 720, 780, 480, 420, 540, 540, 660, 780), ncol = 6, byrow = TRUE, dimnames=list(subj=1:10, cond=c("deg0NA", "deg4NA", "deg8NA", "deg0NP", "deg4NP", "deg8NP")))
mlmfit <- lm(reacttime~1)
mauchly.test(mlmfit, X=~1)
#
#
#
#
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
# Plot average recovery of each component per site over each measuring period
# Evergreen shrubs
vegroot15N_RLong %>%
  filter(Species == "ES") %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  summarise(avgRecovery = mean(Recovery, na.rm = TRUE), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Part == "S", avgRecovery, -avgRecovery)) %>%
  ggplot(aes(Round, avgRecovery, fill = factor(Part, levels=c("S","FR","CR")))) + geom_col(position = "stack") + scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") + facet_wrap( ~ Site) + labs(x = "Measuring period", y = "% of added N", title = "15N recovery for Evergreen shrubs") + guides(x = guide_axis(n.dodge = 2))
#
# Deciduous shrubs
vegroot15N_RLong %>%
  filter(Species == "DS") %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  summarise(avgRecovery = mean(Recovery, na.rm = TRUE), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Part == "S", avgRecovery, -avgRecovery)) %>%
  ggplot(aes(Round, avgRecovery, fill = factor(Part, levels=c("S","FR","CR")))) + geom_col(position = "stack") + scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") + facet_wrap( ~ Site) + labs(x = "Measuring period", y = "% of added N", title = "15N recovery for Deciduous shrubs") + guides(x = guide_axis(n.dodge = 2))
#
# Graminoids
vegroot15N_RLong %>%
  filter(Species == "G") %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  summarise(avgRecovery = mean(Recovery, na.rm = TRUE), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Part == "S", avgRecovery, -avgRecovery)) %>%
  ggplot(aes(Round, avgRecovery, fill = factor(Part, levels=c("S","FR","CR")))) + geom_col(position = "stack") + scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") + facet_wrap( ~ Site) + labs(x = "Measuring period", y = "% of added N", title = "15N recovery for Graminoids") + guides(x = guide_axis(n.dodge = 2))
#
# Sum recovery and calculate average
# This means combining
# Shoots: ES, DS, G, O, U
# CR: ES, DS, G, bulk
# FR: ES, DS, G, O, bulk
# R_G is "G" and removed from here # not anymore, changed naming. It is now part of FR
vegroot15N_RLong %>%
#  mutate(Part = if_else(Part == "R", "FR", Part)) %>%
  group_by(across(c("Site", "Plot", "Round", "Part"))) %>%
  summarise(TotalRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
 # filter(Species != "RootG") %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  summarise(avgRecovery = mean(TotalRecovery, na.rm = TRUE), se = sd(TotalRecovery)/sqrt(length(TotalRecovery)), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Part == "S", avgRecovery, -avgRecovery),
         se = if_else(Part == "S", se, -se),
         avgR_SE = if_else(Part == "CR", avgRecovery, 0)) %>%
  group_by(across(c("Site", "Round"))) %>%
  mutate(avgR_SE = if_else(Part == "FR", cumsum(avgR_SE)+avgRecovery, avgRecovery)) %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  #
  # Plot 
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery, fill = factor(Part, levels=c("S","FR","CR"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(-15,3)) +
  scale_fill_viridis_d() +
  #scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgR_SE, ymax=avgR_SE+se), position=position_dodge(.9)) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(-15, -12, -9, -6, -3, 0, 3), labels = abs) +
  #scale_fill_discrete(labels = c("Shoots", "Fine Roots", "Course roots")) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) + #guides(x = guide_axis(n.dodge = 2)) + 
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))#,#legend.position=c(1,1),
        #legend.justification=c(1, 1))#, 
        #legend.key.width=unit(1, "lines"), 
        #legend.key.height=unit(1, "lines"), 
        #plot.margin = unit(c(5, 1, 0.5, 0.5), "lines")) 
#+ annotate("rect", xmin = winterP2$wstart, xmax = winterP2$wend, ymin = Inf, ymax = -Inf, fill = "grey", alpha = 0.5)
#
#
#
vegroot15N_RLong %>%
  group_by(across(c("Site", "Plot", "Round", "Part", "Species"))) %>%
  summarise(TotalRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Part", "Species"))) %>%
  summarise(avgRecovery = mean(TotalRecovery, na.rm = TRUE), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Part == "S", avgRecovery, -avgRecovery)) %>%
  group_by(across(c("Site", "Round", "Part"))) %>%
  #
  # Plot 
  ggplot() +
  #geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery, fill = factor(Part, levels=c("S","FR","CR"))), position = "stack") +
  coord_cartesian(ylim = c(-8,3)) +
  scale_fill_viridis_d() +
  #scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") +
  #geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgR_SE-se, ymax=avgR_SE+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site + Species, ncol = 7) +
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery in plants") + guides(x = guide_axis(n.dodge = 3)) + 
  theme_light() 



#
#
# Abisko Vassijaure plant recovery facet
vegroot15N_RLong_one %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(TotalRecovery, na.rm = TRUE), se = sd(TotalRecovery)/sqrt(length(TotalRecovery)), .groups = "keep") %>%
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
Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean((TotalRecovery/sysRec*100), na.rm = TRUE), se = sd((TotalRecovery/sysRec*100))/sqrt(length((TotalRecovery/sysRec*100))), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery, ymax=avgRecovery+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgRecovery)) +
  #ylim(0,30) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period", y = expression("% of added "*{}^15*"N"), title = expression({}^15*"N recovery in plants, proportional to total recovery")) + #guides(x = guide_axis(n.dodge = 2)) + 
  theme_classic(base_size = 20)  +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
# Mic - MBN
Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean((R_MBN/sysRec*100), na.rm = TRUE), se = sd((R_MBN/sysRec*100))/sqrt(length((R_MBN/sysRec*100))), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery)) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site, ncol = 2) + 
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery in MBN, proportional to total recovery") + guides(x = guide_axis(n.dodge = 2)) + 
  theme_light() 


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
  summarise(avgR = mean(R_SE, na.rm = TRUE), 
            se = sd(R_SE)/sqrt(length(R_SE)), .groups = "keep") %>%
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
# Soil 15N concentration
Rec15N %>%
  dplyr::rename("atomP_SE" = "atom%_SE") %>%
  group_by(across(c("Site", "Plot", "Round"))) %>%
  mutate(SE_15N = (Nconc_SE*10^-6)*(atomP_SE/100)*10^6) %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgR = mean(SE_15N, na.rm = TRUE), 
            se = sd(SE_15N)/sqrt(length(SE_15N)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = avgR, ymin=avgR-se, ymax=avgR+se), position=position_dodge(.9)) +
  geom_col(aes(Round, avgR),color = "black") +
  coord_cartesian(ylim = c(0,0.5)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", 
       y = expression("["*{}^15*"N] (µg "*{}^15*"N "*g^-1*" DW)"), 
       title = expression({}^15*"N concentration in soil extract")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))








Recov_FR <- ggplot(vegroot15N, aes(Round, R_FR, colour = Site))
Recov_FR + geom_boxplot()

Recov_DSS <- ggplot(vegroot15N, aes(Round, R_DSS, colour = Site))
Recov_DSS + geom_boxplot()
Recov_DSCR <- ggplot(vegroot15N, aes(Round, R_DSCR, colour = Site))
Recov_DSCR + geom_boxplot()
Recov_DSFR <- ggplot(vegroot15N, aes(Round, R_DSFR, colour = Site))
Recov_DSFR + geom_boxplot()


Recov_DS <- vegroot15Nlong %>%
  filter(Species == "DS") %>%
  mutate(Recovery = if_else(Part == "S", Recovery, -Recovery)) %>%
  ggplot(aes(MP, Recovery, colour = Site))
Recov_DS + geom_col(position = "stack")




# Microbial recovery
# SE - TDN
Mic15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(R_SE, na.rm = TRUE), se = sd(R_SE)/sqrt(length(R_SE)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery)) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site, ncol = 2) + 
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery in TDN") + guides(x = guide_axis(n.dodge = 2)) + 
  theme_light()
#
# SEF
Mic15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(R_SEF, na.rm = TRUE), se = sd(R_SEF)/sqrt(length(R_SEF)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery)) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site, ncol = 2) + 
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery SEF") + guides(x = guide_axis(n.dodge = 2)) + 
  theme_light() 
#
# Mic - MBN
Mic15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(R_MBN, na.rm = TRUE), se = sd(R_MBN)/sqrt(length(R_MBN)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery)) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site, ncol = 2) + 
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery in MBN") + guides(x = guide_axis(n.dodge = 2)) + 
  theme_light() 


Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(sysRec, na.rm = TRUE), se = sd(sysRec)/sqrt(length(sysRec)), .groups = "keep") %>%
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery)) +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site, ncol = 2) + 
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery total") + guides(x = guide_axis(n.dodge = 2)) + 
  theme_light() 
#
#
# Total recovery
Total_plot <- Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(sysRec, na.rm = TRUE), se = sd(sysRec)/sqrt(length(sysRec)), .groups = "keep") %>%
  ggplot(aes(y = avgRecovery, x = Round, fill = Site)) +
 # geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), stat="identity", alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_bar(position = position_dodge(), stat = "identity", alpha = 0) + # invisible bars to annotate on top of
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) +
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_bar(position = position_dodge(), stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  scale_fill_viridis_d(option= "plasma") +
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery total") + guides(x = guide_axis(n.dodge = 2)) + 
  theme_light()
#
# Microbial N recovery
MBN_plot <- Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(R_MBN, na.rm = TRUE), se = sd(R_MBN)/sqrt(length(R_MBN)), .groups = "keep") %>%
  ggplot(aes(y = avgRecovery, x = Round, fill = Site)) +
  # geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), stat="identity", alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_bar(position = position_dodge(), stat = "identity", alpha = 0) + # invisible bars to annotate on top of
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) +
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_bar(position = position_dodge(), stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  scale_fill_viridis_d(option= "plasma") +
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery microbial N")+ guides(x = guide_axis(n.dodge = 2))+ 
  theme_light()
#
# TDN
TDN_plot <- Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(R_SE, na.rm = TRUE), se = sd(R_SE)/sqrt(length(R_SE)), .groups = "keep") %>%
  ggplot(aes(y = avgRecovery, x = Round, fill = Site)) +
  # geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), stat="identity", alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_bar(position = position_dodge(), stat = "identity", alpha = 0) + # invisible bars to annotate on top of
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) +
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_bar(position = position_dodge(), stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  scale_fill_viridis_d(option= "plasma") +
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery TDN") + guides(x = guide_axis(n.dodge = 2))+ 
  theme_light()
#
# Plant
plant_plot <- Rec15N %>%
  group_by(across(c("Site", "Round"))) %>%
  summarise(avgRecovery = mean(TotalRecovery, na.rm = TRUE), se = sd(TotalRecovery)/sqrt(length(TotalRecovery)), .groups = "keep") %>%
  ggplot(aes(y = avgRecovery, x = Round, fill = Site)) +
  # geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), stat="identity", alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_bar(position = position_dodge(), stat = "identity", alpha = 0) + # invisible bars to annotate on top of
  annotate("rect", xmin = winterP$wstart[2], xmax = winterP$wend[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) +
  annotate("rect", xmin = winterP$wstart[1], xmax = winterP$wend[1], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_bar(position = position_dodge(), stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin=avgRecovery-se, ymax=avgRecovery+se), position=position_dodge(.9)) +
  scale_fill_viridis_d(option= "plasma") +
  labs(x = "Measuring period", y = "% of added N", title = "15N recovery plants") + guides(x = guide_axis(n.dodge = 2)) + 
  theme_light()
#
# combine figures
grid.arrange(Total_plot, MBN_plot, TDN_plot,plant_plot, ncol = 2)
#
#
#
#
# Combine all types of recovery into one:
vegroot15N_d15N_Nconc_Long1 <- vegroot15N_d15N_Nconc_Long0 %>%
  group_by(across(c("Site", "Plot", "Round"))) %>%
  summarise(Nconc = sum(Nconc, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegroot15N_d15N_Nconc_Long2 <- vegroot15N_d15N_Nconc_Long0 %>%
  group_by(across(c("Site", "Plot", "Round"))) %>%
  summarise(d15N = sum(d15N, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegroot15N_d15N_Nconc_Long <- vegroot15N_d15N_Nconc_Long1 %>% left_join(vegroot15N_d15N_Nconc_Long2)
#
vegroot15N_d15N_Nconc_Long0 %>%
  ggplot(aes(x = Nconc, y = d15N, color = Site)) +
  geom_point()
#
vegroot15N_d15N_Nconc_Long3 <- vegroot15N_d15N_Nconc_Long0 %>%
  group_by(across(c("Site","Plot", "Round", "Part"))) %>%
  summarise(Nconc = sum(Nconc, na.rm = TRUE), .groups = "keep") %>%
  dplyr::rename("Organ" = "Part") %>%
  ungroup()
vegroot15N_d15N_Nconc_Long4 <- vegroot15N_d15N_Nconc_Long0 %>%
  group_by(across(c("Site","Plot", "Round", "Part"))) %>%
  summarise(d15N = sum(d15N, na.rm = TRUE), .groups = "keep") %>% # You CANNOT sum delta15N values!!
  dplyr::rename("Organ" = "Part") %>%
  ungroup()
vegroot15N_d15N_Nconc_Long <- vegroot15N_d15N_Nconc_Long3 %>% left_join(vegroot15N_d15N_Nconc_Long4)
#
vegroot15N_d15N_Nconc_Long %>%
  ggplot(aes(x = log((Nconc)^2+1), y = log((d15N)^2+1), color = Organ)) +
  geom_point()
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