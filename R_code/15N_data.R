# 15N vegetation and root data
# By Emil A.S. Andersen
# 
#=======  ###   Libraries    ### =======
library(plyr)
library(tidyverse)
library(car)
library(nlme)
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
# Long formatted soil extraction data (SE/SEF): [N] (TDN), d15N, atom% 15N for enriched and natural abundance samples, inorganic [N] (NH4 and NO3), and root free soil mass
soil15N <- read_csv("clean_data/Soil_N.csv", col_names = TRUE)
#
#
# Define the winter period as snow covered period
winterP <- data.frame(wstart = c(05, 12), wend = c(12, 13))
winterP2 <- data.frame(wstart = c("05_Nov_19", "05_Nov_19"), wend = c("12_May_20", "13_Jun_20"))
winterP_date <- data.frame(wstart = c(as.Date("2019-11-10"),as.Date("2019-11-12")), wend = c(as.Date("2020-05-27"),as.Date("2020-06-22")))
#
# List of Measuring periods as they should appear in graphs
measuringPeriod <- c("July-19",	"Aug-19",	"Sep-19",	"Oct-19",	"Nov-19",	"Dec-19",	"Jan-20",	"Feb-20",	"Mar-20",	"Apr-20",	"Apr-20",	"May-20",	"Jun-20",	"Jul-20",	"Aug-20")
#
# Days between periods. Will need adjusting as not always 21 days from labeling to harvest and at some point Abisko and Vassijaure shifted which was done first = a few days difference!
# Days between labeling and harvest
dayLH <- 21
#
# Days between harvests
dayHH <- 28
dayHL <- dayHH - dayLH
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
#=======  ###   Functions    ### =======
# From http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}
#
# Plot recovery proportional to total recovery. Split in Abisko and Vassijaure
# dataF: dataframe
# plotvar: variable to plot on y-axis. Needs to be in format dataF$plotvar
# titleExp: the title of the plot
plot_prop_Recovery <- function(dataF=NULL, plotvar, titleExp){
  dataF %>%
    ggplot() + 
    # Plot the snow covered period. This is fixed for all samples
    geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
    # Plot errorbars as 95% confidence interval. The 
    geom_errorbar(aes(x = Round, y = plotvar, ymin=plotvar, ymax=plotvar+ci), position=position_dodge(.9)) +
    # Plot the columns
    geom_col(aes(Round, plotvar), color = "black") +
    # Limit the graph to the range 0-100
    #coord_cartesian(ylim=c(0,100)) +
    # Change the x-axis labels
    scale_x_discrete(labels = measuringPeriod) +
    # Split into Abisko and Vassijaure. Each with their own y-axis
    facet_wrap( ~ Site, ncol = 2, scales = "free") + 
    # Labeling of axis and title
    labs(x = "Measuring period (MP)", y = expression("% of total recovered "*{}^15*"N"), title = titleExp) + 
    # General size of text and lines
    theme_classic(base_size = 20) + 
    # Angle the x-axis and space the two graphs
    theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
}
#
#
#
#=======  ###   Main data    ### =======
#
# Add Round to vegetation data as well as other relevant information
vegroot15N <- vegroot15N %>%
  left_join(coreData, by = join_by(Site, Plot, MP)) %>%
  select(1:13, Round) %>%
  relocate(Round, .after = MP)


# Vegetation recovery for total core
vegroot15N_total_Plant <- vegroot15N %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
#
#
# Calculate recovery for microbial partition
# Pivot wide to get SE and SEF separated
Mic15N <- soil15N %>%
  pivot_wider(names_from = "Extr_type", values_from = c(Nconc_microg_pr_gDW, d15N, Atom_pc, NO3_microg_pr_gDW, NH4_microg_pr_gDW, Atom_pc_NatAb))

# Add Round to microbial data as well as other relevant information
Mic15N <- Mic15N %>%
  left_join(coreData, by = join_by(Site, Plot, MP, Soil_RF_DW_g)) %>%
  select(1:16, Round) %>%
  relocate(Round, .after = MP)

Mic15N <- Mic15N %>%
  mutate(R_TDN = ((Atom_pc_SE - Atom_pc_NatAb_SE)/100 * Nconc_microg_pr_gDW_SE *10^(-6) * Soil_RF_DW_g)/(N_add/1000)* 100) %>%
  mutate(R_TDN = if_else(R_TDN < 0, 0, R_TDN)) %>%
  mutate(atom_pc_MBN = ((Atom_pc_SEF/100 * Nconc_microg_pr_gDW_SEF - Atom_pc_SE/100 * Nconc_microg_pr_gDW_SE) - 
                          (Atom_pc_NatAb_SEF/100 * Nconc_microg_pr_gDW_SEF - Atom_pc_NatAb_SE/100 * Nconc_microg_pr_gDW_SE))/
           (Nconc_microg_pr_gDW_SEF - Nconc_microg_pr_gDW_SE)*100) %>%
  mutate(R_MBN = (((Atom_pc_SEF/100 * Nconc_microg_pr_gDW_SEF*10^(-6) - Atom_pc_SE/100 * Nconc_microg_pr_gDW_SE*10^(-6)) - 
                     (Atom_pc_NatAb_SEF/100 * Nconc_microg_pr_gDW_SEF*10^(-6) - Atom_pc_NatAb_SE/100 * Nconc_microg_pr_gDW_SE*10^(-6)))/
                    K_EN * Soil_RF_DW_g)/(N_add/1000) * 100) %>%
  mutate(R_MBN = if_else(R_MBN < 0, 0, R_MBN))
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
#=======  ### Mineralization ### =======
#
# To calculate mineralization and immobilization
# See Wessel & Tietema 1992, Soil Biol. Biochem. Vol. 24, No. 10, pp. 931-942
# As well as see Bengtson et al. 2005, OIKOS 111: 81-90
# For production
# p = (ln((f_t-k)/(f_0-k)) / ln(W_t/W_0)) * (W_0-W_t)/t
# For consumption
# c = (1 + ln((f_t-k)/(f_0-k)) / ln(W_t/W_0)) * (W_0-W_t)/t
#
# Where
# f is atom%
# k is atom% for natural pool/background
# W is N concentration, mg N pr g
# t is time
#
# except in cases where production and consumption is equal, then
# p = c = - W/t * ln((f_t-k)/(f_0-k))
#
# For mineralization and immobilization the atom% of the inorganic N is needed.
# Info needed:
# Calculate at every harvest point the inorganic atom% (or d15N):
# [N]_inorganic = [NH4] + [NO3]                                     inorgN$NH4_µg_DW & inorgN$NO3_µg_DW
# [N]_total (total of extract, found along Mic15N data)             soil15N$Nconc
# [N]_org = [N]_total - [N]_inorganic (organic N concentration)
# atom%_organic (Organic 15N ratio, equal to background)            soil15N$atom_pc_NatAb
# atom%_inorganic = ?
# atom%_total (the total extracted N's 15N ratio)                   soil15N$atom_pc
#
# Thus
# atom%_tot = (atom%_org * [N]_org + atom%_inorg * [N]_inorg) / [N]_tot
# =>
# atom%_tot * [N]_tot = atom%_org * [N]_org + atom%_inorg * [N]_inorg
# =>
# atom%_tot * [N]_tot - atom%_org * [N]_org = atom%_inorg * [N]_inorg
# =>
# atom%_inorg = (atom%_tot * [N]_tot - atom%_org * [N]_org) / [N]_inorg
#
# New data set with needed information

soil15N_2 <- soil15N %>%
  left_join(coreData, by = join_by(Site, Plot, MP)) %>%
  select(1:11, DaysLH, DaysHH, Round, Injection_mg_pr_patch) %>%
  mutate(DaysHL = DaysHH - DaysLH) %>%
  relocate(c(Round, DaysLH, DaysHH, DaysHL), .after = MP) %>%
  rename(Soil_RF_DW_g = Soil_RF_DW_g.x)

# soil15N <- soil15N %>%
#   mutate(across(c(Plot, MP), as.character))
#inorgN <- inorgN %>%
#  mutate(across(c(Plot, MP), as.character))
#
# mineral <- inorgN %>%
#   select(1:4, NH4_µg_DW, NO3_µg_DW) %>%
#   left_join(soil15N, by = join_by(Site, Plot, MP, SE_SEF)) %>%
#   relocate(Round, .after = MP)
#
# Calculate
mineral <- soil15N_2 %>%
  filter(Extr_type == "SE") %>%
  filter(MP != "EX (w)") %>%
  mutate(Nconc_inorg = NH4_microg_pr_gDW + NO3_microg_pr_gDW) %>% # Inorganic N concentration is equal to the sum of NH4 and NO3: [N]_in = [NH4] + [NO3]
  mutate(Nconc_soil = if_else(Nconc_microg_pr_gDW - Nconc_inorg < 0, Nconc_inorg, Nconc_microg_pr_gDW)) %>%
  mutate(Nconc_org = Nconc_soil - Nconc_inorg) %>% # Organic N concentration is the difference between TDN and inorganic N: [N]_org = [TDN] - [N]_in
  mutate(atom_pc_in_harvest_high = if_else(Nconc_inorg == 0, NA, (Atom_pc * Nconc_soil - Atom_pc_NatAb * Nconc_org)/Nconc_inorg)) %>% # Assuming organic atom% does not change considerably from natural abundance (!!!)
  mutate(atom_pc_in_harvest_low = if_else(Nconc_inorg == 0, NA, Atom_pc_NatAb),
         atom_pc_in_harvest_low2 = if_else(Nconc_inorg == 0, NA, (Atom_pc * Nconc_soil - Atom_pc * Nconc_org)/Nconc_inorg)) %>% # Since atom%_organic < atom%_total*[N]_total/[N]_organic when [N]_org < [N]_tot
  mutate(atom_pc_in_extreme = if_else(Nconc_inorg == 0, NA, (Atom_pc * Nconc_soil - (Atom_pc+0.1) * Nconc_org)/Nconc_inorg))
#
# arrange values in order of MP so that the following loop works correctly
mineral <- mineral %>%
  arrange(Plot) %>%
  arrange(Site) %>%
  arrange(MP)
#
# Extrapolate linearly the following for labeling point (1/4 from one harvest to the next):
# Make subset of data for Abisko and Vassijaure
mineral_A <- mineral %>%
  select(Site, Plot, MP, Round, Extr_type) %>%
  filter(Site == "Abisko")
mineral_A["Nconc_in0.A"] <- NA
mineral_V <- mineral %>%
  select(Site, Plot, MP, Round, Extr_type) %>%
  filter(Site == "Vassijaure")
mineral_V["Nconc_in0.V"] <- NA
# Double check I can count
mt<- vector("double")
# Loop over each to calculate rate from one round to the next for the same plot and then around 7 days if when labeling happens
for (i in 2:15) {
  for (j in 1:5) {
    mt[5*(i-1)+j] <- 5*(i-1)+j
    mineral_A$Nconc_in0.A[5*(i-1)+j] <- (mineral$Nconc_inorg[which(mineral$Site == "Abisko" & 
                                                                     mineral$Plot == j & 
                                                                     mineral$MP == (i))] - 
                                          mineral$Nconc_inorg[which(mineral$Site == "Abisko" & 
                                                                     mineral$Plot == j & 
                                                                     mineral$MP == (i-1))])/mineral$DaysHH[which(mineral$Site == "Abisko" & 
                                                                                                                   mineral$Plot == j &
                                                                                                                   mineral$MP == (i))] * mineral$DaysHL[which(mineral$Site == "Abisko" &
                                                                                                                                                                  mineral$Plot == j &
                                                                                                                                                                  mineral$MP == (i))] +
                                          mineral$Nconc_inorg[which(mineral$Site == "Abisko" & 
                                                                     mineral$Plot == j & 
                                                                     mineral$MP == (i-1))]
    mineral_V$Nconc_in0.V[5*(i-1)+j] <- (mineral$Nconc_inorg[which(mineral$Site == "Vassijaure" & 
                                                                     mineral$Plot == j & 
                                                                     mineral$MP == i)] - 
                                          mineral$Nconc_inorg[which(mineral$Site == "Abisko" & 
                                                                     mineral$Plot == j & 
                                                                     mineral$MP == (i-1))])/mineral$DaysHH[which(mineral$Site == "Abisko" & 
                                                                                                                   mineral$Plot == j &
                                                                                                                   mineral$MP == (i))] * mineral$DaysHL[which(mineral$Site == "Abisko" &
                                                                                                                                                                  mineral$Plot == j &
                                                                                                                                                                  mineral$MP == (i))] + 
                                          mineral$Nconc_inorg[which(mineral$Site == "Vassijaure" & 
                                                                     mineral$Plot == j & 
                                                                     mineral$MP == (i-1))]
  }
}


(mineral$Nconc_inorg[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == 6)] - mineral$Nconc_inorg[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == 5)])/mineral$DaysHH[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == (6-1))] * mineral$DaysHL[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == (6-1))] + mineral$Nconc_inorg[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == (6-1))]

(mineral$Nconc_inorg[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == 2)] - mineral$Nconc_inorg[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == 5)])/dayHH * dayHL + mineral$Nconc_inorg[which(mineral$Site == "Abisko" & mineral$Plot == 5 & mineral$MP == (2-1))]



test <- mineral %>%
  left_join(mineral_A, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  left_join(mineral_V, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  mutate(Nconc_in0 = if_else(!is.na(Nconc_in0.A), Nconc_in0.A, Nconc_in0.V)) %>%
  select(!c(Nconc_in0.A, Nconc_in0.V))

# Combine Abisko and Vassijaure and join the values
mineral_test <- mineral %>%
  left_join(mineral_A, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  left_join(mineral_V, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  mutate(Nconc_in0 = if_else(!is.na(Nconc_in0.A), Nconc_in0.A, Nconc_in0.V)) %>%
  select(!c(Nconc_in0.A, Nconc_in0.V))
#
# Add core mass in DW
mic_mass <- Mic15N %>%
  select(Site, Plot, MP, Soil_RF_DW_g)# %>%
  mutate(across(c(Plot, MP), as.character))
#
mineral_test <- mineral_test %>%
  left_join(mic_mass, by = join_by(Site, Plot, MP, Soil_RF_DW_g)) %>%
  relocate(Soil_RF_DW_g, .after = Extr_type)
#
# calculate how much 15N was added and add it to the inorganic pool at labeling
mineral_test <- mineral_test %>%
  mutate(inj_15N = (N_add*1000)/Soil_RF_DW_g) %>% # µg N added pr g DW
  mutate(Nconc_in0 = if_else(Nconc_in0 <= 0, 0, Nconc_in0)) %>%
  mutate(Nconc_in0_l = Nconc_in0 + inj_15N) %>%
#  mutate(Nconc_in0_l_2 = (Nconc_in0*Mic_mass + N_add*1000)/Mic_mass) %>% # Same as above but made explicit
  mutate(atom_pc_in0_l = (98.7*inj_15N + Atom_pc_NatAb*Nconc_in0)/Nconc_in0_l) #  NB!!! 98.7% atom% for label
#
# Calculate isotope ratio at injection and harvest. Then calculate average ratio over the 3 weeks and convert to average AP
mineral_test <- mineral_test %>%
  mutate(R_inj = if_else(is.na(atom_pc_in0_l), NA, atom_pc_in0_l/100 / (1 - atom_pc_in0_l/100)),
         R_harv_high = if_else(is.na(atom_pc_in_harvest_high), NA, atom_pc_in_harvest_high/100 / (1 - atom_pc_in_harvest_high/100)),
         R_harv_low = if_else(is.na(atom_pc_in_harvest_low), NA, atom_pc_in_harvest_low/100 / (1 - atom_pc_in_harvest_low/100))) %>%
  # If the inorganic N pool is depleted at harvest, it is assumed that the ratio stays constant
  mutate(R_avg_high = if_else(is.na(R_inj), NA, if_else(is.na(R_harv_high), R_inj, (R_inj - R_harv_high)/2)),
         R_avg_low = if_else(is.na(R_inj), NA, if_else(is.na(R_harv_low), R_inj, (R_inj - R_harv_low)/2))) %>%
  # Calculate atom% from the ratio
  mutate(AP_avg_high = if_else(is.na(R_avg_high), NA, R_avg_high/(1 + R_avg_high)*100),
         AP_avg_low = if_else(is.na(R_avg_low), NA, R_avg_low/(1 + R_avg_low)*100))





# Now the mineralization can be calculated
mineral_test <- mineral_test %>%
  mutate(gross_p = (log((atom_pc_in_harvest_low2 - Atom_pc_NatAb # ft - k
                         ) / (atom_pc_in0_l - Atom_pc_NatAb) # f0 - k
                        ) / log(Nconc_inorg / Nconc_in0_l) # log(Wt/W0)
                    ) * ((Nconc_in0_l - Nconc_inorg # W0 - Wt
                          ) / dayLH), 
         # gross production
         gross_c = (1 + log((atom_pc_in_harvest_low2 - Atom_pc_NatAb # ft - k
                           ) / (atom_pc_in0_l - Atom_pc_NatAb) # f0 - k
                          ) / log(Nconc_inorg / Nconc_in0_l) # log(Wt/W0)
                    ) * ((Nconc_in0_l-Nconc_inorg # W0 - Wt
                        ) / dayLH)
         ) # gross consumption
#
mineral_test %>%
  ggplot(aes(Round, gross_p)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
mineral_test %>%
  ggplot(aes(Round, gross_c)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")


mineral_test %>%
  ggplot(aes(Round, Nconc_inorg)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

# How does the 15N/14N ratio change
mineral_test %>%
  ggplot(aes(Round, R_avg_low)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
mineral_test %>%
  ggplot(aes(Round, R_avg_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")



# Testing on plant uptake
vegroot15N_format <- vegroot15N %>%
  #mutate(across(c("Plot", "MP"), as.character)) %>%
  rename("Atom_pc_plant" = Atom_pc)
vegMineral <- vegroot15N_format %>%
  left_join(mineral_test, by = join_by(Site, Plot, MP, Round)) %>%
  dplyr::select(1:15, 34, 39:42)

vegMineral <- vegMineral %>%
  mutate(ext14N_high = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc/100 * Biomass)/R_avg_high,
         ext14N_low = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc/100 * Biomass)/R_avg_low) %>%
  mutate(NtotRec_high = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc/100 * Biomass)/(AP_avg_high/100),
         NtotRec_low = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc/100 * Biomass)/(AP_avg_low)/100) %>%
  mutate(frac15Nstart = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc/100 * Biomass)/Nconc_in0_l,
         Ntot_stand_high = NtotRec_high/Nconc_in0_l,
         Ntot_stand_low = NtotRec_low/Nconc_in0_l,
         Recov_stand = Recovery/Nconc_in0_l)


vegMineral_total_0 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_1 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantExt14N_high = sum(ext14N_high, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_2 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantExt14N_low = sum(ext14N_low, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_3 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantN_high = sum(NtotRec_high, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_4 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantN_low = sum(NtotRec_low, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_5 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(Plantstd_high = sum(Ntot_stand_high, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_6 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(Plantstd_low = sum(Ntot_stand_low, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_7 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(Plantstd_15N = sum(frac15Nstart, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_8 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantRecovstd = sum(Recov_stand, na.rm = TRUE), .groups = "keep") %>%
  ungroup()

vegMineral_total <- vegMineral_total_0 %>%
  left_join(vegMineral_total_1, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_2, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_3, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_4, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_5, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_6, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_7, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_8, by = join_by(Site, Plot, MP, Round))


vegMineral_total %>%
  ggplot(aes(Round, PlantRecovery)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

vegMineral_total %>%
  ggplot(aes(Round, PlantExt14N_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(Round, PlantExt14N_low)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

vegMineral_total %>%
  ggplot(aes(Round, PlantN_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(Round, PlantN_low)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

vegMineral_total %>%
  ggplot(aes(Round, Plantstd_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(Round, Plantstd_low)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(Round, Plantstd_15N)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

vegMineral_total %>%
  ggplot(aes(Plantstd_low, Plantstd_high)) + geom_smooth() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(PlantRecovstd, Plantstd_low)) + geom_smooth() + facet_wrap(vars(Site), scales = "free")

vegMineral_total %>%
  ggplot(aes(PlantN_low, PlantN_high)) + geom_smooth() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(PlantN_low, PlantN_high)) + geom_point() + facet_wrap(vars(Site), scales = "free")


vegMineral_total %>%
  dplyr::select(1:4, 8,9) %>%
  pivot_longer(cols = 5:6, names_to = "Estimate", values_to = "PlantN_uptake") %>%
  mutate(across("MP", as.numeric)) %>%
  ggplot(aes(Round, PlantN_uptake, color = Site)) + geom_boxplot() + facet_wrap(vars(Estimate), scales = "free")

vegMineral_total %>%
  dplyr::select(1:4, 10,11) %>%
  pivot_longer(cols = 5:6, names_to = "Estimate", values_to = "PlantN_uptake_std") %>%
  mutate(across("MP", as.numeric)) %>%
  ggplot(aes(MP, PlantN_uptake_std, color = Site)) + geom_point() + facet_wrap(vars(Estimate), scales = "free")
vegMineral_total %>%
  dplyr::select(1:4, 10,11) %>%
  pivot_longer(cols = 5:6, names_to = "Estimate", values_to = "PlantN_uptake_std") %>%
  mutate(across("MP", as.numeric)) %>%
  ggplot(aes(Round, PlantN_uptake_std, color = Site)) + geom_boxplot() + facet_wrap(vars(Estimate), scales = "free")
vegMineral_total %>%
  dplyr::select(1:4, 10,11) %>%
  pivot_longer(cols = 5:6, names_to = "Estimate", values_to = "PlantN_uptake_std") %>%
  mutate(across("MP", as.numeric)) %>%
  ggplot(aes(Round, PlantN_uptake_std, color = Estimate)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")


vegMineral_total %>%
  ggplot(aes(PlantRecovery, PlantRecovstd)) + geom_smooth() + facet_wrap(vars(Site), scales = "free")


vegMineral_total %>%
  ggplot(aes(Round, PlantRecovstd)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")




vegroot15N_total_Plant %>%
  ggplot(aes(Round, PlantRecovery)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

#
#
#=======  ###   Statistics   ### =======
#-------   ##       Q1       ## -------
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
# Check contrasts are orthogonal
crossprod(cbind(SummervsWinter,SpringvsAutumn,SnowvsNot))#, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh, cont11, cont12, cont13, cont14))


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
#
#model:
lme1<-lme(arcPlantRecovery ~ Round*Site,
          random = ~1|Plot/Site,
          data = vegroot15N_total_Plant, na.action = na.exclude, method = "REML")
#
#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1), resid(lme1), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1), main = "Normally distributed?")                 
qqline(resid(lme1), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1)
par(mfrow = c(1,1))
#
#model output
Anova(lme1, type=2)
summary(lme1)
# Significant for Round
#
#
#
#=======  ###   Statistics   ### =======
#-------   ##       Q1a      ##  -------
#
# Model
# Response variable: plant recovery of 15N (% of combined functional groups partitioned into organs)
# Factors: Time, Site, Organ, Time*Organ, Time*Site, Site*Organ, Time*Site*Organ
# Most important factors: Organ, Time*Organ, Time*Site*Organ
#
vegroot15N_Organ <- vegroot15N %>%
  select(1:4,Species,Organ,Recovery) %>%
  group_by(across(c("Site","Plot", "MP", "Round", "Organ"))) %>%
  summarise(OrganRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  ungroup() %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "Organ"), as.factor)) %>%
  left_join(vegroot15N_total_Plant, by = join_by(Site, Plot, MP, Round)) %>%
  select(1:7) %>%
  mutate(OrganRecovery_ofTotal = OrganRecovery, 
         OrganRecovery = OrganRecovery/PlantRecovery*100) %>%
  select(1:5,OrganRecovery)
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
  # select(1:5,7) %>%
  # rename(OrganRecovery = OrganRecovery_frac) %>%
  #dplyr::filter(OrganRecovery > 0) %>% 
  mutate(sqrtOrganRecovery = sqrt(OrganRecovery)) %>%
  mutate(invOrganRecovery = 1/OrganRecovery) %>%
  mutate(logOrganRecovery = log(OrganRecovery+1)) %>%
  mutate(expOrganRecovery = exp(OrganRecovery+1)) %>%
  mutate(rootOrganRecovery = (OrganRecovery^2)^(1/9)) %>%
  mutate(sqOrganRecovery = OrganRecovery^2) %>%
  mutate(arcOrganRecovery = asin(sqrt((OrganRecovery)/100))) %>%
  mutate(logsqrtOrganRecovery = sqrt(log(OrganRecovery+1))) %>%
  mutate(logarcOrganRecovery = log(asin(sqrt((OrganRecovery)/100))+1)) %>%
  mutate(funOrganRecovery = (OrganRecovery/100)^(1/9))
# BoxCox transformation?
#
# model:
lme1a<-lme(arcOrganRecovery ~ Round*Site*Organ,
          random = ~1|Plot/Site,
          data = vegroot15N_Organ, na.action = na.exclude , method = "REML")
#
#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1a), resid(lme1a), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1a), main = "Normally distributed?")                 
qqline(resid(lme1a), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1a)
par(mfrow = c(1,1))
#
#model output
Anova(lme1a, type=2)
summary(lme1a)
# For organ recovery as fraction of whole plant recovery:
# Highly significant for organ and all interactions except Round:Site
#
#
#
#=======  ###   Statistics   ### =======
#-------   ##       Q2       ##  -------
#
# Model
# Response variable: Microbial recovery of 15N (%), MBN
# Factors: Time, Site, Time*Site
# Most important factor: Time
#
# Load data from excel instead of calculated combined
Mic15N_R <- Rec15N %>%
  select("Site", "Plot", "MP", "Round", "R_MBN")
Mic15N_R <- Mic15N_R %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round"), as.factor))
#
# Contrasts - MBN recovery
# Month           (J, A, S, O, N, D, J, F, M, A, A, M, J, J, A) # Two times April
# MP              (1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15)
Cont4_Mic     <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0) # July vs July
Cont5_Mic     <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0) # 2 times april
Cont6_Mic     <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-2, 0, 0, 0) # Spring change
Cont7_Mic     <- c(0, 0, 0, 0, 2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0) # November is highly different from the rest of winter?
Cont8_Mic     <- c(0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # Or from autumn?
Cont9_Mic     <- c(1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # Just getting the last to work 
Cont10_Mic    <- c(1, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Cont11_Mic    <- c(1, 1, 1, 1, 1,-5, 0, 0, 0, 0, 0, 0, 0, 0, 0)
AvsV<-c(1,-1)
contrasts(Mic15N_R$Site)<-AvsV
# Same contrasts as for plant recovery:
# SummervsWinter, SpringvsAutumn, SnowvsNot
contrasts(Mic15N_R$Round)<-cbind(SummervsWinter,SpringvsAutumn,SnowvsNot, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh, cont11, cont12, cont13, cont14)
#contrasts(Mic15N_R$Round)<-cbind(SummervsWinter, SpringvsAutumn, SnowvsNot)#, Cont4_Mic, Cont5_Mic, Cont6_Mic, Cont7_Mic, Cont8_Mic, Cont9_Mic, Cont10_Mic, Cont11_Mic) # Contrasts that compare each new round with the previous ones.
#contrasts(Mic15N_R$Round)<-contr.helmert
#
# Check if contrasts work
# Check contrasts are orthogonal
crossprod(cbind(SummervsWinter,SpringvsAutumn,SnowvsNot, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh))#, cont11, cont12, cont13, cont14))
#
# Alternative: Two-way ANOVA
# have time as a factor in a two-way ANOVA, combined with Site. As each sampling is destructive, the samples are technically independent of each other, although it does not account for the block design
MicModel3 <- aov(R_MBN ~ Round*Site, data = Mic15N_R)
Anova(MicModel3, type ="III")
#
#
# Transform data
Mic15N_R <- Mic15N_R %>%
  mutate(sqrtR_MBN = sqrt(R_MBN)) %>%
  mutate(logR_MBN = log(R_MBN)) %>%
  mutate(cubeR_MBN = (R_MBN)^(1/3)) %>%
  mutate(arcR_MBN = asin(sqrt(R_MBN/max(Mic15N_R$R_MBN))))
#
#
#model:
lme2<-lme(R_MBN ~ Site*Round,
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
#
#model output
Anova(lme2, type=2)
summary(lme2)
# Highly significant for Round
#
#
#
#
#
#=======  ###    Plotting    ### =======
#-------   ## Plant Biomass  ## -------
#
vegroot15N_bio <- vegroot15N %>%
  group_by(across(c("Site", "Plot", "Round"))) %>%
  summarise(TotalBiomass = sum(Biomass_DW_g, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
#
# Plant biomass total +/- SE
vegroot15N_bio %>%
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
# Plant biomass by organ +/- SE
vegroot15N %>%
  group_by(across(c("Site", "Plot", "Round", "Organ"))) %>%
  summarise(TotalBiomass = sum(Biomass_DW_g, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  summarise(avgBiomass = mean(TotalBiomass, na.rm = TRUE), se = sd(TotalBiomass)/sqrt(length(TotalBiomass)), .groups = "keep") %>%
  mutate(avgBiomass = if_else(Organ == "S", avgBiomass, -avgBiomass),
         se = if_else(Organ == "S", se, -se),
         avgR_SE = if_else(Organ == "CR", avgBiomass, 0)) %>%
  group_by(across(c("Site", "Round"))) %>%
  mutate(avgR_SE = if_else(Organ == "FR", cumsum(avgR_SE)+avgBiomass, avgBiomass)) %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  #
  # Plot 
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgBiomass, fill = factor(Organ, levels=c("S","FR","CR"))), position = "stack", color = "black") +
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
#-------   ##    Recovery    ## -------
#
# For species recovery, it is also possible to 
vegroot15N_Species <- vegroot15N %>%
  select(1:4,Species,Organ,Recovery) %>%
  #group_by(across(c("Site","Plot", "MP", "Round", "Organ"))) %>%
  #summarise(OrganRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  #ungroup() %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "Organ"), as.factor)) %>%
  left_join(vegroot15N_total_Plant, by = join_by(Site, Plot, MP, Round)) %>%
  select(1:8) %>%
  mutate(SpOrganRecovery = Recovery/PlantRecovery*100)
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
# Same as above, but Recovery proportional to total plant recovery
# organRecovery / PlantRecovery
vegroot15N_Organ %>%
  group_by(across(c("Site", "Plot", "Round", "Organ"))) %>%
  summarise(TotalRecovery = sum(OrganRecovery, na.rm = TRUE), .groups = "keep") %>%
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
  coord_cartesian(ylim = c(-110,75)) +
  scale_fill_viridis_d() +
  #  scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") +
  geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgR_SE, ymax=avgR_SE+se), position=position_dodge(.9)) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75), labels = abs) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of total plant recovered "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) +
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Each species or part separated. For a quick overview of where patterns might come from
vegroot15N %>%
  #group_by(across(c("Site", "Plot", "Round", "Organ", "Species"))) %>%
  #summarise(TotalRecovery = sum(Recovery, na.rm = TRUE), .groups = "keep") %>%
  group_by(across(c("Site", "Round", "Organ", "Species"))) %>%
  summarise(avgRecovery = mean(Recovery, na.rm = TRUE), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Organ == "S", avgRecovery, -avgRecovery)) %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  #
  # Plot 
  ggplot() +
  #geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery, fill = factor(Organ, levels=c("S","FR","CR"))), position = "stack") +
  coord_cartesian(ylim = c(-8,3)) +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = abs) +
  #scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") +
  #geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgR_SE-se, ymax=avgR_SE+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site + Species, ncol = 5) +
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery per species per organ")) +
  guides(fill = guide_legend(title = "Plant organ")) + 
  theme_light(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1)) 
#
# Species, excluding the belowground Bulk roots. Proportional to total recovered in plants
vegroot15N_Species %>%
  group_by(across(c("Site", "Round", "Organ", "Species"))) %>%
  summarise(avgRecovery = mean(SpOrganRecovery, na.rm = TRUE), .groups = "keep") %>%
  mutate(avgRecovery = if_else(Organ == "S", avgRecovery, -avgRecovery)) %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  filter(Species != "Bulk" & Species != "BulkG") %>%
  #
  # Plot 
  ggplot() +
  #geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, avgRecovery, fill = factor(Organ, levels=c("S","FR","CR"))), position = "stack") +
  coord_cartesian(ylim = c(-20,40)) +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = abs) +
  #scale_fill_manual(values = c("darkgreen", "navy", "brown"), name = "Recovery") +
  #geom_errorbar(aes(x = Round, y = avgRecovery, ymin=avgR_SE-se, ymax=avgR_SE+se), position=position_dodge(.9)) +
  facet_wrap( ~ Site + Species, ncol = 5) +
  labs(x = "Measuring period (MP)", y = expression("% of total plant recovered "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery per species per organ")) +
  guides(fill = guide_legend(title = "Plant organ")) + 
  theme_light(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1)) 
#
#
#
# Abisko and Vassijaure plant and microbial recovery faceted
# Calculate means and 95% CI
vegroot15N_total_Plant_sum <- summarySE(vegroot15N_total_Plant, measurevar="PlantRecovery", groupvars=c("Site", "Round"))
Mic15N_sum <- summarySE(Mic15N, measurevar="R_MBN", groupvars=c("Site", "Round"), na.rm=TRUE)
TDN15N_sum <- summarySE(Mic15N, measurevar="R_TDN", groupvars=c("Site", "Round"), na.rm=TRUE)
#
# Same calculations as with the summarySE function, but less flexible and would need to write code each place
# vegroot15N_total_Plant %>%
#   group_by(across(c("Site", "Round"))) %>%
#   summarise(avgRecovery = mean(PlantRecovery, na.rm = TRUE), 
#                    sd = sd(PlantRecovery),
#                    se = sd(PlantRecovery)/sqrt(length(PlantRecovery)), 
#                    N = length(PlantRecovery), # 5 replicates
#                    ci = qt(0.95/2 + .5, length(PlantRecovery)-1) * (sd(PlantRecovery)/sqrt(length(PlantRecovery))), # t-distribution of 95% (97.5% as 2.5% each end) for N-1 times standard error
#                    .groups = "keep")
#
# Plant total recovery +/- 95% CI
vegroot15N_total_Plant_sum %>%  
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = PlantRecovery, ymin=PlantRecovery-ci, ymax=PlantRecovery+ci), position=position_dodge(.9)) +
  #geom_point(aes(Round, PlantRecovery)) +
  geom_col(aes(Round, PlantRecovery),color = "black") +
  coord_cartesian(ylim=c(0,30)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) + 
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
vegroot15N_total_Plant %>%
  ggplot() +
  geom_boxplot(aes(Round, PlantRecovery))
#
# Microbial total recovery +/- 95% CI
Mic15N_sum %>%  
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = R_MBN, ymin=R_MBN, ymax=R_MBN+ci), position=position_dodge(.9)) +
  geom_col(aes(Round, R_MBN),color = "black") +
  #coord_cartesian(ylim=c(0,30)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Microbial "*{}^15*"N tracer recovery")) + 
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
# TDN total recovery +/- 95% CI 
TDN15N_sum %>%  
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = R_TDN, ymin=R_TDN, ymax=R_TDN+ci), position=position_dodge(.9)) +
  geom_col(aes(Round, R_TDN),color = "black") +
  #coord_cartesian(ylim=c(0,30)) +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("TDN "*{}^15*"N tracer recovery")) + 
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Proportional to total recovery
# Calculate means and 95% CI
Rec15N_Plant_sum <- summarySE(Rec15N, measurevar = "PlantR_frac", groupvars = c("Site", "Round"))
Rec15N_TDN_sum <- summarySE(Rec15N, measurevar = "R_TDN_frac", groupvars = c("Site", "Round"))
Rec15N_MBN_sum <- summarySE(Rec15N, measurevar = "R_MBN_frac", groupvars = c("Site", "Round"))
# Plot values using function
plot_prop_Recovery(Rec15N_Plant_sum, plotvar=Rec15N_Plant_sum$PlantR_frac, titleExp = expression("Plant "*{}^15*"N tracer recovery"))
plot_prop_Recovery(Rec15N_TDN_sum, plotvar=Rec15N_TDN_sum$R_TDN_frac, titleExp = expression("TDN "*{}^15*"N tracer recovery"))
plot_prop_Recovery(Rec15N_MBN_sum, plotvar=Rec15N_MBN_sum$R_MBN_frac, titleExp = expression("Microbial "*{}^15*"N tracer recovery"))
#
#
vegroot15N_Organ_sum <- summarySE(vegroot15N_Organ, measurevar = "OrganRecovery", groupvars = c("Site", "Round", "Organ"))
#
vegroot15N_Organ_sum %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  mutate(OrganRecovery = if_else(Organ == "S", OrganRecovery, -OrganRecovery),
         ci = if_else(Organ == "S", ci, -ci),
         avgR_CI = if_else(Organ == "CR", OrganRecovery, 0)) %>%
  group_by(across(c("Site", "Round"))) %>%
  mutate(avgR_CI = if_else(Organ == "FR", cumsum(avgR_CI)+OrganRecovery, OrganRecovery)) %>%
  group_by(across(c("Site", "Round", "Organ"))) %>%
  # Plot 
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, OrganRecovery, fill = factor(Organ, levels=c("S","FR","CR"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(-125,75)) +
  scale_fill_viridis_d() +
  geom_errorbar(aes(x = Round, y = OrganRecovery, ymin=avgR_CI, ymax=avgR_CI+ci), position=position_dodge(.9)) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(-125, -100, -75, -50, -25, 0, 25, 50, 75), labels = abs) +
  #scale_fill_discrete(labels = c("Shoots", "Fine Roots", "Course roots")) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of total plant recovered "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery per organ")) + #guides(x = guide_axis(n.dodge = 2)) + 
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
vegroot15N_Organ_sum %>%
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, OrganRecovery, fill = factor(Organ, levels=c("S","FR","CR"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(0,100)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), labels = abs) +
  #scale_fill_discrete(labels = c("Shoots", "Fine Roots", "Course roots")) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of total plant recovered "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery per organ")) + #guides(x = guide_axis(n.dodge = 2)) + 
  guides(fill = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
Rec15N_sum <- Rec15N_Plant_sum %>%
  left_join(Rec15N_MBN_sum, by = join_by(Site, Round)) %>%
  left_join(Rec15N_TDN_sum, by = join_by(Site, Round)) %>%
  select(Site, Round, PlantR_frac, ci.x, R_MBN_frac, ci.y, R_TDN_frac, ci) %>%
  rename("ci_plant" = ci.x,
         "ci_MBN" = ci.y,
         "ci_TDN" = ci)
Rec15N_sum <- Rec15N_sum %>%
  select(Site, Round, PlantR_frac, R_MBN_frac, R_TDN_frac) %>%
  pivot_longer(cols = 3:5, names_to = "Type", values_to = "Recovery")
#
Rec15N_sum %>%
  ggplot() +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_col(aes(Round, Recovery, fill = factor(Type, levels=c("PlantR_frac","R_MBN_frac","R_TDN_frac"))), position = "stack", color = "black") +
  coord_cartesian(ylim = c(0,100)) +
  scale_fill_viridis_d(labels = c("Plant", "Microbial", "TDN")) +
  scale_x_discrete(labels = measuringPeriod) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), labels = abs) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of total system recovered "*{}^15*"N"), title = expression("Plant, microbial, and TDN "*{}^15*"N tracer recovery per part of the system")) + #guides(x = guide_axis(n.dodge = 2)) + 
  guides(fill = guide_legend(title = "System part")) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"),axis.text.x=element_text(angle=60, hjust=1))
#
#
# Attempting to add the 95% CI as smooth line
test <- Rec15N %>%
  left_join(coreData, by = join_by(Site, Plot, MP, Round)) %>%
  select(1:11, Day_of_harvest) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x)))

# Extract date for sample. Here using Day of harvest
DayOf <- coreData %>%
  select(Site, Round, Day_of_harvest) %>%
  distinct(Day_of_harvest, .keep_all = TRUE)
#
# Plant recovery fraction
Rec15N_Plant_sum2 <- Rec15N_Plant_sum %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, PlantR_frac, ci, Day_of_harvest) %>%
  add_column(Type = "Plant_frac") %>%
  rename(Recov_frac = "PlantR_frac")
#
# Microbial recovery fraction
Rec15N_MBN_sum2 <- Rec15N_MBN_sum %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, R_MBN_frac, ci, Day_of_harvest) %>%
  add_column(Type = "MBN_frac") %>%
  rename(Recov_frac = "R_MBN_frac")
#
# TDN recovery fraction
Rec15N_TDN_sum2 <- Rec15N_TDN_sum %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, R_TDN_frac, ci, Day_of_harvest) %>%
  add_column(Type = "TDN_frac") %>%
  rename(Recov_frac = "R_TDN_frac")
#
# Combine recovery types to one file
Rec15N_sum2 <- Rec15N_Plant_sum2 %>%
  bind_rows(Rec15N_MBN_sum2) %>%
  bind_rows(Rec15N_TDN_sum2)
#
# Graph recovery by type and add 95% CI
Rec15N_sum2 %>%
  ggplot(aes(x = Day_of_harvest, y = Recov_frac, ymin = Recov_frac-ci, ymax = Recov_frac+ci, fill=Type, linetype=Type)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("Microbial", "Plant", "TDN")) +
  scale_linetype(labels = c("Microbial", "Plant", "TDN")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))+
  labs(x = "Time of harvest", y = expression("% of total system recovered "*{}^15*"N"), title = expression("Plant, microbial, and TDN "*{}^15*"N tracer recovery per part of the system"))+# to wrap the title properly around use atop() ))) +
  facet_wrap( ~ Site, ncol = 2)+#, scales = "free") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))

Rec15N_sum2 %>%
  ggplot(aes(x = Day_of_harvest, y = Recov_frac, ymin = Recov_frac-ci, ymax = Recov_frac+ci, fill=Type, linetype=Site)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("Microbial", "Plant", "TDN")) +
  scale_linetype(labels = c("Abisko", "Vassijaure")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  labs(x = "Time of harvest", y = expression("% of total system recovered "*{}^15*"N"), title = expression("Plant, microbial, and TDN "*{}^15*"N tracer recovery per part of the system")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))
#
# Recovery with all plots still available
Rec15N_2 <- Rec15N %>%
  pivot_longer(cols = c(PlantR_frac, R_MBN_frac, R_TDN_frac), names_to = "Type", values_to = "Recovery") %>%
  select(1:4, Type, Recovery) %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x)))
#
# Plot with geom_smooth
Rec15N_2 %>%
  ggplot(aes(x = Day_of_harvest, y = Recovery, color = Type)) +
  geom_point() +
  geom_smooth(se = TRUE) +
  facet_wrap( ~ Site)



#
#
#
#=======  ###  { The End }   ### =======