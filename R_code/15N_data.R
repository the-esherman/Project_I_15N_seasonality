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
winterP2 <- data.frame(wstart = c("05_Nov_19", "05_Nov_19"), wend = c("11_Apr_20", "12_May_20"))
winterP_date <- data.frame(wstart = c(as.Date("2019-11-10"),as.Date("2019-11-12")), wend = c(as.Date("2020-05-06"),as.Date("2020-06-01")))
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
Label_atom_pc <- 0.987
#
# Extraction correction factor
K_EN = 0.4
# See https://climexhandbook.w.uib.no/2019/11/06/soil-microbial-biomass-c-n-and-p/ and UCPH bio lab protocol (where K_EN = 0.4)
#
# Extract date for sample. Here using Day of harvest
DayOf <- coreData %>%
  select(Site, Round, Day_of_harvest) %>%
  distinct(Day_of_harvest, .keep_all = TRUE)
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
# Add Round to vegetation data
vegroot15N <- vegroot15N %>%
  left_join(coreData, by = join_by(Site, Plot, MP)) %>%
  select(1:13, Round) %>%
  relocate(Round, .after = MP)
#
#
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
#
# Add Round to microbial data
Mic15N <- Mic15N %>%
  left_join(coreData, by = join_by(Site, Plot, MP, Soil_RF_DW_g)) %>%
  select(1:16, Round, Injection_15N_mg_pr_patch) %>%
  relocate(Round, .after = MP)
#
# Calculate recovery in TDN (SE) and microbial parts
# For microbial recovery a few steps are needed:
# Calculate a new Atom% for the microbial part based on fumigated (SEF) minus non-fumigated (SE) extractions:
#    ((AP_SEF*[N]_SEF - AP*[N]_SE) - (AP_NatAb_SEF*[N]_SEF - AP_NatAb*[N]_SE)) / (([N]_SEF - [N]_SE)*100)
# Then calculate recovery
#    
Mic15N <- Mic15N %>%
  # Recovery of TDN
  mutate(R_TDN = ((Atom_pc_SE - Atom_pc_NatAb_SE)/100 * Nconc_microg_pr_gDW_SE *10^(-6) * Soil_RF_DW_g)/(N_add/1000)* 100) %>%
  mutate(R_TDN = if_else(R_TDN < 0, 0, R_TDN)) %>% # Recovery cannot be negative
  # Recovery of MBN
  # AP for MBN
  mutate(atom_pc_MBN = ((Atom_pc_SEF/100 * Nconc_microg_pr_gDW_SEF - Atom_pc_SE/100 * Nconc_microg_pr_gDW_SE) - 
                          (Atom_pc_NatAb_SEF/100 * Nconc_microg_pr_gDW_SEF - Atom_pc_NatAb_SE/100 * Nconc_microg_pr_gDW_SE))/
           (Nconc_microg_pr_gDW_SEF - Nconc_microg_pr_gDW_SE)*100) %>%
  # Recovery
  mutate(R_MBN = (((Atom_pc_SEF/100 * Nconc_microg_pr_gDW_SEF*10^(-6) - Atom_pc_SE/100 * Nconc_microg_pr_gDW_SE*10^(-6)) - 
                     (Atom_pc_NatAb_SEF/100 * Nconc_microg_pr_gDW_SEF*10^(-6) - Atom_pc_NatAb_SE/100 * Nconc_microg_pr_gDW_SE*10^(-6)))/
                    K_EN * Soil_RF_DW_g)/(Injection_15N_mg_pr_patch/1000) * 100) %>%
  mutate(R_MBN = if_else(R_MBN < 0, 0, R_MBN)) # Recovery cannot be negative
  # # Method 2:
  # # {15N_lab} & {15N_NatAb}
  # mutate( # {15N_lab}
  #         MBN_15N_lab = (Atom_pc_SEF/100 * Nconc_microg_pr_gDW_SEF*10^(-6) - Atom_pc_SE/100 * Nconc_microg_pr_gDW_SE*10^(-6))*Soil_RF_DW_g/K_EN,
  #        # {15N_NatAb}
  #        MBN_15N_NatAb = (Atom_pc_NatAb_SEF/100 * Nconc_microg_pr_gDW_SEF*10^(-6) - Atom_pc_NatAb_SE/100 * Nconc_microg_pr_gDW_SE*10^(-6))*Soil_RF_DW_g/K_EN) %>%
  # # Recovery:
  # # ({15N_lab} - {15N_NatAb})/15N_inj * 100 
  # mutate(R_MBN2 = ((MBN_15N_lab - MBN_15N_NatAb)*100)/(Injection_15N_mg_pr_patch/1000))
  # Gives exactly the same results
#
# Calculate recovery as a proportion of total recovered in each plot
Rec15N <- vegroot15N_total_Plant %>%
  left_join(Mic15N, by = join_by(Site, Plot, MP, Round)) %>%
  select(Site, Plot, MP, Round, PlantRecovery, R_TDN, R_MBN) %>%
  rowwise() %>%
  mutate(sysRec = sum(PlantRecovery, R_TDN, R_MBN, na.rm = TRUE)) %>%
  mutate(PlantR_frac = PlantRecovery/sysRec*100,
         R_TDN_frac = R_TDN/sysRec*100,
         R_MBN_frac = R_MBN/sysRec*100) %>%
  ungroup()
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
  left_join(coreData, by = join_by(Site, Plot, MP, Soil_RF_DW_g)) %>%
  select(1:11, DaysLH, DaysHH, Round, Injection_15N_mg_pr_patch, Injection_N_mg_pr_patch) %>%
  mutate(DaysHL = DaysHH - DaysLH) %>%
  relocate(c(Round, DaysLH, DaysHH, DaysHL), .after = MP)


test2 <- coreData %>%
  mutate(Soil_RF_FW_g_all = Soil_RF_FW_g*Soil_mass_g)
#
test3 <- summarySE(test2, measurevar = "Soil_RF_FW_g_all", groupvars = c("Site", "MP"))

#
# Calculate
mineral <- soil15N_2 %>%
  # Need only non-fumigated extraction (SE) - this is the soil N content
  filter(Extr_type == "SE") %>%
  # Remove control samples
  filter(MP != "EX (w)") %>%
  # Inorganic N concentration is equal to the sum of NH4 and NO3: [N]_in = [NH4] + [NO3] µg pr g DW
  mutate(Nconc_inorg = NH4_microg_pr_gDW + NO3_microg_pr_gDW) %>%
  #The total N concentration in soil cannot be less than the total inorganic concentration!
  # Three cases where this has to be corrected: A_4_4, A_3_4, V_5_8
  mutate(Nconc_soil = if_else(Nconc_microg_pr_gDW - Nconc_inorg < 0, Nconc_inorg, Nconc_microg_pr_gDW)) %>%
  # Organic N concentration is the difference between TDN and inorganic N: [N]_org = [TDN] - [N]_in µg pr g DW
  mutate(Nconc_org = Nconc_soil - Nconc_inorg) %>%
  # Two estimates of AP:
  # High: Assuming organic atom% does not change considerably from natural abundance (!!!)
  mutate(atom_pc_in_harvest_high = if_else(Nconc_inorg == 0, NA, (Atom_pc * Nconc_soil - Atom_pc_NatAb * Nconc_org)/Nconc_inorg),
         atom_pc_in_harvest_high2 = if_else(Nconc_inorg == 0, NA, Label_atom_pc*100)) %>% 
  # Low: Assuming inorganic atom% does not change considerably from natural abundance (!!!)
  # Or assume [N]_org < [N]_total, then AP_org < AP_total*[N]_total/[N]_org
  mutate(atom_pc_in_harvest_low = if_else(Nconc_inorg == 0, NA, Atom_pc_NatAb),
         atom_pc_in_harvest_low2 = if_else(Nconc_inorg == 0, NA, (Atom_pc * Nconc_soil - Atom_pc * Nconc_org)/Nconc_inorg)) %>% # Since atom%_organic < atom%_total*[N]_total/[N]_organic when [N]_org < [N]_tot
  # Third alternative:
  # How far can the atom% be taken? Can the organic N pool have a atom% higher than TDN?
  mutate(atom_pc_in_extreme = if_else(Nconc_inorg == 0, NA, (Atom_pc * Nconc_soil - (Atom_pc+0.000003) * Nconc_org)/Nconc_inorg)) %>%
  mutate(atom_pc_in_harvest_high = if_else(atom_pc_in_harvest_high > Label_atom_pc*100, Label_atom_pc*100, atom_pc_in_harvest_high),
         atom_pc_in_harvest_low2)
# test <- mineral %>%
#   mutate(atom_pc_in_extreme = if_else(Nconc_inorg == 0, NA, (Atom_pc * Nconc_soil - (Atom_pc+0.000003) * Nconc_org)/Nconc_inorg)) %>%
#   mutate(deltaAP_extreme = atom_pc_in_extreme - Atom_pc_NatAb) %>%
#   mutate(AP_plus = if_else(Nconc_org == 0, NA, (atom_pc_in_harvest_low*Nconc_inorg + Atom_pc * Nconc_soil)/Nconc_org - Atom_pc)) %>%
#   mutate(deltaNconc = Nconc_soil - Nconc_inorg)
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
mt<- vector("double") # the count should start 5* NA and then 6-75
MP_plot<- matrix(NA,nrow = 15*5,ncol = 2) # Making a matrix with MP and plot at the correct spot
#
# Loop over each to calculate rate from one round to the next for the same plot and then around 7 days if when labeling happens
for (i in 2:15) { # 15MP excluding first round
  for (j in 1:5) { # 5 replicates
    # Test if looping count is correct: Skip first 5 (as NA), continue counting from 6 to 75 (15 MP * 5 replicates)
    mt[5*(i-1)+j] <- 5*(i-1)+j
    MP_plot[5*(i-1)+j,1] <- i # MP
    MP_plot[5*(i-1)+j,2] <- j # Plot
    # Abisko initial inorganic N pool
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
    # Vassijaure initial inorganic N pool
    mineral_V$Nconc_in0.V[5*(i-1)+j] <- (mineral$Nconc_inorg[which(mineral$Site == "Vassijaure" & 
                                                                   mineral$Plot == j & 
                                                                   mineral$MP == i)] - 
                                         mineral$Nconc_inorg[which(mineral$Site == "Vassijaure" & 
                                                                   mineral$Plot == j & 
                                                                   mineral$MP == (i-1))])/mineral$DaysHH[which(mineral$Site == "Vassijaure" & 
                                                                                                               mineral$Plot == j &
                                                                                                               mineral$MP == (i))] * mineral$DaysHL[which(mineral$Site == "Vassijaure" &
                                                                                                                                                                  mineral$Plot == j &
                                                                                                                                                                  mineral$MP == (i))] + 
                                         mineral$Nconc_inorg[which(mineral$Site == "Vassijaure" & 
                                                                   mineral$Plot == j & 
                                                                   mineral$MP == (i-1))]
  }
}

#
#
#




test <- mineral %>%
  left_join(mineral_A, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  left_join(mineral_V, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  mutate(Nconc_in0 = if_else(!is.na(Nconc_in0.A), Nconc_in0.A, Nconc_in0.V)) %>%
  select(!c(Nconc_in0.A, Nconc_in0.V))

# Combine Abisko and Vassijaure and join the values
mineral_combined <- mineral %>%
  left_join(mineral_A, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  left_join(mineral_V, by = join_by(Site, Plot, MP, Round, Extr_type)) %>%
  # Values alternate between Abisko and Vassijaure and should be combined into one value
  mutate(Nconc_in0 = if_else(!is.na(Nconc_in0.A), Nconc_in0.A, Nconc_in0.V)) %>% 
  select(!c(Nconc_in0.A, Nconc_in0.V)) %>%
  relocate(Soil_RF_DW_g, .after = Extr_type)
#
# calculate how much 15N was added and add it to the inorganic pool at labeling
mineral_combined <- mineral_combined %>%
  # N label (98.7% 15N) injected per dry weight soil
  mutate(inj_15N = (Injection_N_mg_pr_patch*1000)/Soil_RF_DW_g) %>% # µg N added pr g DW
  # Inorganic N concentration at injection point should be set at min 0
  mutate(Nconc_in0 = if_else(Nconc_in0 <= 0, 0, Nconc_in0)) %>%
  # Total inorganic N pool is injected N and inorganic N, µg pr g DW
  mutate(Nconc_in0_l = Nconc_in0 + inj_15N) %>%
  # Atom% of the inorganic N pool at injection. Corresponds to the 15N added (in µg pr g DW) plus natural abundance 15N in pre-labelled pool divided by new pool size
  mutate(atom_pc_in0_l = ((Injection_15N_mg_pr_patch*1000)/Soil_RF_DW_g + (Atom_pc_NatAb/100)*Nconc_in0)/Nconc_in0_l*100) #  NB!!! 98.7% atom% for label
#
# How much of a fertilizer was the injection? 
N_fertilizer <- mineral_combined %>%
  mutate(deltaInj_pc = inj_15N/Nconc_in0_l*100) %>%
  select(1:4,inj_15N, Nconc_in0, deltaInj_pc)
N_fertilizer_sum <- summarySE(N_fertilizer, measurevar = "deltaInj_pc", groupvars = c("Site", "MP"))
N_fertilizer_sum %>%
  ggplot(aes(x = MP, y = deltaInj_pc, color = Site)) + 
  geom_line() + 
  scale_x_continuous(breaks = 1:15) +
  labs(x = "Time of harvest", y = expression("% = "*frac("[N]"*{}[label],"[N]"*{}[inorg] + "[N]"*{}[label])*symbol("\264")*100*" (µg N pr g DW)"), title = "Label [N] as % of total inorganic [N] at time of injection") +
  theme_classic(base_size = 20)
#
#
# Calculate isotope ratio (isoR) at injection and harvest. Then calculate average ratio over the 3 weeks and convert to average AP
mineral_combined <- mineral_combined %>%
  mutate(isoR_inj = if_else(is.na(atom_pc_in0_l), NA, atom_pc_in0_l/100 / (1 - atom_pc_in0_l/100)),
         isoR_harv_high = if_else(is.na(atom_pc_in_harvest_high), NA, atom_pc_in_harvest_high/100 / (1 - atom_pc_in_harvest_high/100)),
         isoR_harv_low = if_else(is.na(atom_pc_in_harvest_low), NA, atom_pc_in_harvest_low/100 / (1 - atom_pc_in_harvest_low/100)),
         isoR_harv_low2 = if_else(is.na(atom_pc_in_harvest_low2), NA, atom_pc_in_harvest_low2/100 / (1 - atom_pc_in_harvest_low2/100))) %>%
  # If the inorganic N pool is depleted at harvest, it is assumed that the ratio stays constant
  mutate(isoR_avg_high = if_else(is.na(isoR_inj), NA, if_else(is.na(isoR_harv_high), isoR_inj, (isoR_inj + isoR_harv_high)/2)),
         isoR_avg_low = if_else(is.na(isoR_inj), NA, if_else(is.na(isoR_harv_low), isoR_inj, (isoR_inj + isoR_harv_low)/2)),
         isoR_avg_low2 = if_else(is.na(isoR_inj), NA, if_else(is.na(isoR_harv_low2), isoR_inj, (isoR_inj + isoR_harv_low2)/2))) %>%
  # Calculate atom% from the ratio
  mutate(AP_avg_high = if_else(is.na(isoR_avg_high), NA, isoR_avg_high/(1 + isoR_avg_high)*100),
         AP_avg_low = if_else(is.na(isoR_avg_low), NA, isoR_avg_low/(1 + isoR_avg_low)*100),
         AP_avg_low2 = if_else(is.na(isoR_avg_low2), NA, isoR_avg_low2/(1 + isoR_avg_low2)*100)) %>%
  mutate(isoF_avg_high = if_else(is.na(atom_pc_in0_l), NA, if_else(is.na(atom_pc_in_harvest_high), atom_pc_in0_l/100, (atom_pc_in0_l/100 + atom_pc_in_harvest_high/100)/2)),
         isoF_avg_low = if_else(is.na(atom_pc_in0_l), NA, if_else(is.na(atom_pc_in_harvest_low), atom_pc_in0_l/100, (atom_pc_in0_l/100 + atom_pc_in_harvest_low/100)/2)),
         isoF_avg_low2 = if_else(is.na(atom_pc_in0_l), NA, if_else(is.na(atom_pc_in_harvest_low2), atom_pc_in0_l/100, (atom_pc_in0_l/100 + atom_pc_in_harvest_low2/100)/2))) %>%
  mutate(isoR_avg_high_v2 = isoF_avg_high/(1-isoF_avg_high),
         isoR_avg_low_v2 = isoF_avg_low/(1-isoF_avg_low))

test <- mineral_combined %>%
  mutate(delta_inj = Injection_N_mg_pr_patch - Injection_15N_mg_pr_patch) %>%
  relocate(delta_inj, .before = inj_15N) %>%
  mutate(delta_Nconc = Nconc_org/Nconc_soil) %>%
  relocate(delta_Nconc, .after = Nconc_org) %>%
  filter(delta_Nconc != 1)
  # mutate(delta_isoR = isoR_avg_high - isoR_avg_low) %>%
  # select(1:4, atom_pc_in_harvest_high, atom_pc_in_harvest_low, isoR_avg_high, isoR_avg_low, delta_isoR) %>%
  # filter(delta_isoR < 0)

test4 <- summarySE(test, measurevar = "delta_Nconc", groupvars = c("Site", "Round"))
ggplot(test4, aes(x = Round, y = delta_Nconc, fill = Site)) + geom_col()

# How much of the total N pool is in the organic fraction
ggplot(test, aes(x = Round, y = delta_Nconc, fill = Site)) + 
  geom_boxplot() +
  labs(x = "Time of harvest", y = expression("Fraction "*frac("[N]"*{}[organic],"[N]"*{}[total])), title = "Organic [N] as fraction of total inorganic [N] at time of harvest") +
  theme_classic(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))
#
#







# Is the 15R always lower in harvest than at injection?
# That is, is there more 15N per 14N at harvest than at injection?
# This would imply that there is fractionation, or a bias against using the 15N over using the 14N. Which is possible and likely, but does not take into account dilution from mineralization.
mineral_combined %>%
  select(1:3, isoR_inj, isoR_harv_high, isoR_harv_low, isoR_harv_low2) %>%
  filter(!is.na(isoR_inj)) %>%
  mutate(isoR_harv_high = if_else(is.na(isoR_harv_high), isoR_inj, isoR_harv_high),
         isoR_harv_low = if_else(is.na(isoR_harv_low), isoR_inj, isoR_harv_low),
         isoR_harv_low2 = if_else(is.na(isoR_harv_low2), isoR_inj, isoR_harv_low2)) %>%
  #group_by(across(c("Site", "MP"))) %>%
  summarise(across(c(isoR_inj, isoR_harv_high, isoR_harv_low, isoR_harv_low2), ~mean(.)), .by = c("Site", "MP")) %>%
  pivot_longer(cols = 3:6, names_to = "Time", values_to = "isoR") %>%
  mutate(Time = case_when(Time == "isoR_inj" ~ 1,
                          Time == "isoR_harv_high" ~ 3,
                          Time == "isoR_harv_low" ~ 2,
                          Time == "isoR_harv_low2" ~ 2.5)) %>%
  #filter(Time != 3) %>%
  ggplot(aes(x = Time, y = isoR, color = Site)) + geom_point() + facet_wrap(~ MP)
  


mineral_combined %>%
  ggplot(aes(x = MP, y = Nconc_inorg, color = Site)) + geom_point()

mineral_combined %>%
  filter(!is.na(isoR_harv_high),
         !is.na(isoR_inj)) %>%
  ggplot(aes(x = isoR_avg_high, y = isoR_avg_high_v2, color = Site)) + geom_point()

mineral_combined %>%
  filter(!is.na(isoR_harv_high),
         !is.na(isoR_inj)) %>%
  ggplot(aes(x = isoR_avg_high/(1+isoR_avg_high), y = isoF_avg_high, color = Site)) + geom_point()



mineral_combined %>%
  ggplot(aes(x = Round, y = atom_pc_in_harvest_high, fill = Site)) + geom_boxplot() + coord_cartesian(ylim = c(0,1))

mineral_combined %>%
  ggplot(aes(x = Round, y = atom_pc_in_harvest_low2, fill = Site)) + geom_boxplot() + coord_cartesian(ylim = c(0,1))


mineral_combined %>%
  mutate(deltaAP_high_low2 = atom_pc_in_harvest_high - atom_pc_in_harvest_low2) %>%
  ggplot(aes(x = Round, y = deltaAP_high_low2, fill = Site)) + geom_boxplot() + coord_cartesian(ylim = c(0,1))


test1 <- summarySE(mineral_combined, measurevar = "atom_pc_in_harvest_high", groupvars = c("Site", "Round"), na.rm = TRUE)
test2 <- summarySE(mineral_combined, measurevar = "atom_pc_in_harvest_low2", groupvars = c("Site", "Round"), na.rm = TRUE)
test_delta <- mineral_combined %>%
  mutate(deltaAP = atom_pc_in_harvest_low2 / atom_pc_in_harvest_high* 100)
test_delta <- summarySE(test_delta, measurevar = "deltaAP", groupvars = c("Site", "Round"), na.rm = TRUE)

test1 <- test1 %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, Day_of_harvest, atom_pc_in_harvest_high, ci) %>%
  rename("AP_high" = atom_pc_in_harvest_high,
         "ci_high" = ci)
test2 <- test2 %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, Day_of_harvest, atom_pc_in_harvest_low2, ci) %>%
  rename("AP_low" = atom_pc_in_harvest_low2,
         "ci_low" = ci)
test_delta <- test_delta %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, Day_of_harvest, deltaAP, ci) %>%
  rename("AP_delta" = deltaAP,
         "ci_delta" = ci)


test3 <- test1 %>%
  left_join(test2, by = join_by(Site, Round, Day_of_harvest)) %>%
  left_join(test_delta, by = join_by(Site, Round, Day_of_harvest))

test3_1 <- test3 %>%
  rename("High" = AP_high,
         "Low" = AP_low) %>%
  select(Site, Round, Day_of_harvest, High, Low) %>%
  pivot_longer(4:5, names_to = "Estimate", values_to = "AP")
test3_2 <- test3 %>%
  rename("High" = ci_high,
         "Low" = ci_low) %>%
  select(Site, Round, Day_of_harvest, High, Low) %>%
  pivot_longer(4:5, names_to = "Estimate", values_to = "ci")
test4 <- test3_1 %>%
  left_join(test3_2, by = join_by(Site, Round, Day_of_harvest, Estimate))

test3 %>%
  ggplot(aes(x = Day_of_harvest, y = AP_delta, ymin = AP_delta-ci_delta, ymax = AP_delta+ci_delta, fill = Site)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("Abisko", "Vassijaure")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  coord_cartesian(ylim = c(0,100)) +
  labs(x = "Measuring period (MP)", y = expression("Atom% "*{}^15*"N"), title = expression("Soil "*{}^15*"N atom% high estimate")) +
  theme_light(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))

test4 %>%
  ggplot(aes(x = Day_of_harvest, y = AP, ymin = AP-ci, ymax = AP+ci, fill = Estimate)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.4, fill = 'grey', inherit.aes = FALSE) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("High", "Low")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  coord_cartesian(ylim = c(0,100)) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") +
  labs(x = "Measuring period (MP)", y = expression("Atom% "*{}^15*"N"), title = expression("Soil "*{}^15*"N atom% high estimate")) +
  theme_light(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))





#  
#
# Now the mineralization can be calculated
mineral_combined <- mineral_combined %>%
  # Low estimate
  mutate(gross_p_low2 = (log((atom_pc_in_harvest_low2 - Atom_pc_NatAb # ft - k
                         ) / (atom_pc_in0_l - Atom_pc_NatAb) # f0 - k
                        ) / log(Nconc_inorg / Nconc_in0_l) # log(Wt/W0)
                    ) * ((Nconc_in0_l - Nconc_inorg # W0 - Wt
                          ) / dayLH), 
         # gross production
         gross_c_low2 = (1 + log((atom_pc_in_harvest_low2 - Atom_pc_NatAb # ft - k
                           ) / (atom_pc_in0_l - Atom_pc_NatAb) # f0 - k
                          ) / log(Nconc_inorg / Nconc_in0_l) # log(Wt/W0)
                    ) * ((Nconc_in0_l-Nconc_inorg # W0 - Wt
                        ) / dayLH)
         ) %>% # gross consumption
  # High estimate
  mutate(gross_p_high = (log((atom_pc_in_harvest_high - Atom_pc_NatAb # ft - k
                             ) / (atom_pc_in0_l - Atom_pc_NatAb) # f0 - k
                            ) / log(Nconc_inorg / Nconc_in0_l) # log(Wt/W0)
                        ) * ((Nconc_in0_l - Nconc_inorg # W0 - Wt
                              ) / dayLH), 
         # gross production
         gross_c_high = (1 + log((atom_pc_in_harvest_high - Atom_pc_NatAb # ft - k
                          ) / (atom_pc_in0_l - Atom_pc_NatAb) # f0 - k
                         ) / log(Nconc_inorg / Nconc_in0_l) # log(Wt/W0)
                 ) * ((Nconc_in0_l-Nconc_inorg # W0 - Wt
                       ) / dayLH)
         ) # gross consumption
#



mineral_combined %>%
  ggplot(aes(Round, gross_p_low2)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

mineral_combined %>%
  ggplot(aes(Round, gross_p_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

mineral_combined %>%
  ggplot(aes(Round, gross_c_low2)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

mineral_combined %>%
  ggplot(aes(Round, gross_c_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

mineral_combined %>%
  ggplot(aes(Round, Nconc_inorg)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

# How does the 15N/14N ratio change
mineral_combined %>%
  ggplot(aes(Round, isoR_avg_low)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
mineral_combined %>%
  ggplot(aes(Round, isoR_avg_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")



#
# Testing on plant uptake ----
vegroot15N_format <- vegroot15N %>%
  #mutate(across(c("Plot", "MP"), as.character)) %>%
  rename("Atom_pc_plant" = Atom_pc,
         "d15N_plant" = d15N,
         "R_plant" = Recovery)
vegMineral <- vegroot15N_format %>%
  left_join(mineral_combined, by = join_by(Site, Plot, MP, Round)) %>%
  select(1:7, Soil_RF_DW_g, Biomass_DW_g, Nconc_pc, d15N_plant, Atom_pc_plant, NatAb_atom_pc, R_plant, 20:50)

vegMineral <- vegMineral %>%
  mutate(ext14N_high = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/isoR_avg_high,
         ext14N_low = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/isoR_avg_low,
         ext14N_low2 = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/isoR_avg_low2) %>%
  mutate(NtotRec_high = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/(AP_avg_high/100),
         NtotRec_low = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/(AP_avg_low/100),
         NtotRec_low2 = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/(AP_avg_low2/100)) %>%
  mutate(frac15Nstart = ((Atom_pc_plant - NatAb_atom_pc)/100 * Nconc_pc/100 * Biomass_DW_g)/Nconc_in0_l,
         Ntot_stand_high = NtotRec_high/Nconc_in0_l,
         Ntot_stand_low = NtotRec_low/Nconc_in0_l,
         Ntot_stand_low2 = NtotRec_low2/Nconc_in0_l,
         Recov_stand = R_plant/Nconc_in0_l) %>%
  mutate(ext14N_high_pr15N = ext14N_high/(Injection_15N_mg_pr_patch/1000),
         ext14N_low_pr15N = ext14N_low/(Injection_15N_mg_pr_patch/1000),
         ext14N_low2_pr15N = ext14N_low2/(Injection_15N_mg_pr_patch/1000),
         NtotRec_high_pr15N = NtotRec_high/(Injection_15N_mg_pr_patch/1000),
         NtotRec_low_pr15N = NtotRec_low/(Injection_15N_mg_pr_patch/1000),
         NtotRec_low2_pr15N = NtotRec_low2/(Injection_15N_mg_pr_patch/1000))



vegMineral %>%
  ggplot(aes(x = ext14N_high, y = ext14N_low, color = Site)) + geom_point()
vegMineral %>%
  ggplot(aes(x = ext14N_high_pr15N, y = ext14N_low_pr15N, color = Site)) + geom_point()
vegMineral %>%
  ggplot(aes(x = ext14N_low2_pr15N, y = ext14N_low_pr15N, color = Site)) + geom_point()
vegMineral %>%
  ggplot(aes(x = NtotRec_high_pr15N, y = NtotRec_low_pr15N, color = Site)) + geom_point()
vegMineral %>%
  ggplot(aes(x = NtotRec_high, y = NtotRec_low, color = Site)) + geom_point()



vegMineral_total_0 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantRecovery = sum(R_plant, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_1 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantExt14N_high = sum(ext14N_high_pr15N, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_2 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantExt14N_low = sum(ext14N_low_pr15N, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_2_5 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantExt14N_low2 = sum(ext14N_low2_pr15N, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_3 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantN_high = sum(NtotRec_high_pr15N, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_4 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantN_low = sum(NtotRec_low_pr15N, na.rm = TRUE), .groups = "keep") %>%
  ungroup()
vegMineral_total_4_5 <- vegMineral %>%
  group_by(across(c("Site", "Plot", "MP", "Round"))) %>%
  summarise(PlantN_low2 = sum(NtotRec_low2_pr15N, na.rm = TRUE), .groups = "keep") %>%
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
  left_join(vegMineral_total_2_5, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_3, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_4, by = join_by(Site, Plot, MP, Round)) %>%
  left_join(vegMineral_total_4_5, by = join_by(Site, Plot, MP, Round)) %>%
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
  ggplot(aes(Round, PlantExt14N_low2)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

vegMineral_total %>%
  ggplot(aes(Round, PlantN_high)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(Round, PlantN_low)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")
vegMineral_total %>%
  ggplot(aes(Round, PlantN_low2)) + geom_boxplot() + facet_wrap(vars(Site), scales = "free")

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
#-------   ##   Contrasts    ## -------
#
# Abisko vs Vassijaure
AvsV<-c(1,-1)
#
#
# Measuring period - Snow vs Snow-free and in between
#
# Abisko
# Month             ( J, A, S, O, N, D, J, F, M, A, A, M, J, J, A) # Two times April
# MP                ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15)
Abi_SnowvsNot   <- c(-8,-8,-8,-8, 7, 7, 7, 7, 7, 7, 7, 7,-8,-8,-8) # Snow found on plots from November to May
Abi_SnowCvsW    <- c( 0, 0, 0, 0, 3, 3, 3, 3, 3,-5,-5,-5, 0, 0, 0) # Cold vs warm part of the snow covered period
Abi_freeWvsC    <- c( 2, 2,-5,-5, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2) # Warm vs cold part of the snow-free period
# The rest of the contrasts: Necessary for balanced test, but not interesting
Abi_Cont4       <- c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 1, 1)
Abi_Cont5       <- c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1) # Summer 2019 vs 2020
Abi_Cont6       <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-2, 0, 0, 0) # April(x2) vs May
Abi_Cont7       <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0)
Abi_Cont8       <- c( 0, 0, 0, 0, 3, 3,-2,-2,-2, 0, 0, 0, 0, 0, 0) # Nov-Dec vs Jan-Mar
Abi_Cont9       <- c( 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Abi_Cont10      <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1)
Abi_Cont11      <- c( 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Abi_Cont12      <- c( 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Abi_Cont13      <- c( 0, 0, 0, 0, 0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0)
Abi_Cont14      <- c( 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0)
#
Contr_Abisko_MP <- cbind(Abi_SnowvsNot, Abi_SnowCvsW, Abi_freeWvsC, Abi_Cont4, Abi_Cont5, Abi_Cont6, Abi_Cont7, Abi_Cont8, Abi_Cont9, Abi_Cont10, Abi_Cont11, Abi_Cont12, Abi_Cont13, Abi_Cont14)
#
# Check contrasts are orthogonal
crossprod(Contr_Abisko_MP)
#
# Vassijaure
# Month             ( J, A, S, O, N, D, J, F, M, A, A, M, J, J, A) # Two times April
# MP                ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15)
Vas_SnowvsNot   <- c(-7,-7,-7,-7, 8, 8, 8, 8, 8, 8, 8, 8, 8,-7,-7) # Snow found on plots from November to June
Vas_SnowCvsW    <- c( 0, 0, 0, 0, 4, 4, 4, 4, 4,-5,-5,-5,-5, 0, 0) # Cold vs warm part of the snow covered period
Vas_freeWvsC    <- c( 2, 2,-4,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2) # Warm vs cold part of the snow-free period
# The rest of the contrasts: Necessary for balanced test, but not interesting
Vas_Cont4       <- c( 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Vas_Cont5       <- c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1) # Summer 2019 vs 2020
Vas_Cont6       <- c( 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Vas_Cont7       <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1)
Vas_Cont8       <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0) # April(x2) vs May & June
Vas_Cont9       <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0)
Vas_Cont10      <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0)
Vas_Cont11      <- c( 0, 0, 0, 0, 3, 3,-2,-2,-2, 0, 0, 0, 0, 0, 0) # Nov-Dec vs Jan-Mar
Vas_Cont12      <- c( 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Vas_Cont13      <- c( 0, 0, 0, 0, 0, 0, 1, 1,-2, 0, 0, 0, 0, 0, 0)
Vas_Cont14      <- c( 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0)
#
Contr_Vassijaure_MP <- cbind(Vas_SnowvsNot, Vas_SnowCvsW, Vas_freeWvsC, Vas_Cont4, Vas_Cont5, Vas_Cont6, Vas_Cont7, Vas_Cont8, Vas_Cont9, Vas_Cont10, Vas_Cont11, Vas_Cont12, Vas_Cont13, Vas_Cont14)
#
# Check contrasts are orthogonal
crossprod(Contr_Vassijaure_MP)
#
#
# Plant organs
SvsR<-c(-1, -1, 2) # Shoots vs roots
CRvsFR<-c(1,-1,0) # Coarse roots vs fine roots
#
# Check contrasts are orthogonal
crossprod(cbind(SvsR,CRvsFR))
#
#
#
#=======  ###   Statistics   ### =======
#-------   ##     Q Zero     ## -------
#
Q0_ecosys_stat <- Rec15N %>%
  mutate(across(c("Plot", "MP"), as.character)) %>%
  mutate(across(c("Site", "MP", "Round"), as.factor))
#
# transform data
Q0_ecosys_stat <-  Q0_ecosys_stat %>%
  mutate(Recov = sysRec)
Q0_ecosys_stat <- Q0_ecosys_stat %>%
  mutate(logRecov = log(Recov+1), # Good for low percentage values.
         arcRecov = asin(sqrt(Recov/1000))) # General use is for this transformation. Recovery values are > 100, so transform by dividing by 1000 instead
#
# model:
lme0<-lme(arcRecov ~ Round*Site,
          random = ~1|Plot/Site,
          data = Q0_ecosys_stat, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme0), resid(lme0), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme0), main = "Normally distributed?")                 
qqline(resid(lme0), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme0)
par(mfrow = c(1,1))
#
# model output
Anova(lme0, type=2)
# sysRec, Arcsin transformation:
# Highly significant for Round (χ^2 = 50.0671, p = 5.95*10^-6). Trend for interaction (χ^2 = 23.1476, p = 0.05791)
#
#
# Per site
Q0_ecosys_stat_A <- Q0_ecosys_stat %>%
  filter(Site == "Abisko")
Q0_ecosys_stat_V <- Q0_ecosys_stat %>%
  filter(Site == "Vassijaure")
#
#
# Contrasts - Abisko
contrasts(Q0_ecosys_stat_A$Round) <- Contr_Abisko_MP
#
# Check if contrasts work, by using a two-way ANOVA
SysModel_alias <- aov(sysRec ~ Round, data = Q0_ecosys_stat_A)
Anova(SysModel_alias, type ="III")
# alias checks dependencies
alias(SysModel_alias)
#
# transform data
Q0_ecosys_stat_A <- Q0_ecosys_stat_A %>%
  mutate(Recov = sysRec)
Q0_ecosys_stat_A <- Q0_ecosys_stat_A %>%
  mutate(logRecov = log(Recov+1), # Good for low percentage values.
         arcRecov = asin(sqrt(Recov/1000))) # General use is for this transformation. Recovery values are > 100, so transform by dividing by 1000 instead
#
# model:
lme0_A<-lme(arcRecov ~ Round,
            random = ~1|Plot,
            data = Q0_ecosys_stat_A, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme0_A), resid(lme0_A), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme0_A), main = "Normally distributed?")                 
qqline(resid(lme0_A), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme0_A)
par(mfrow = c(1,1))
#
# model output
Anova(lme0_A, type=2)
summary(lme0_A)
# Highly significant for Snow Cold vs Warm (t = 2.60138, p = 0.0119)
# Significant for Summer vs Autumn (t = -2.05424, p = 0.0446)
#
Q0_season_A <- Q0_ecosys_stat_A %>%
  select(1:4, sysRec) %>%
  mutate(Snow = if_else(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 | MP == 10 | MP == 11 | MP == 12, "Snow","Clear"),
         SnowCW = case_when(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 ~ "Cold",
                              MP == 10 | MP == 11 | MP == 12 ~ "Warm",
                              TRUE ~ NA),
         FreeWC = case_when(MP == 1 | MP == 2 | MP == 13 | MP == 14 | MP == 15 ~ "Warm",
                              MP == 3 | MP == 4 ~ "Cold",
                              TRUE ~ NA))
#
#
Q0_season_A %>%
  filter(!is.na(SnowCW)) %>%
  summarise(sysRec = mean(sysRec), .by = c(SnowCW))
Q0_season_A %>%
  filter(!is.na(FreeWC)) %>%
  summarise(sysRec = mean(sysRec), .by = c(FreeWC))
#
#
Q0_season_A %>%
  summarise(sysRec_season = mean(sysRec), .by = c(Snow))
#
summarySE(Q0_season_A, measurevar = "sysRec", groupvars = c("Snow"))
summarySE(Q0_season_A, measurevar = "sysRec", groupvars = c("SnowCW"))
summarySE(Q0_season_A, measurevar = "sysRec", groupvars = c("FreeWC"))

Q0_season_A %>%
  ggplot(aes(x = Round, y = sysRec, fill = FreeWC)) +
  geom_boxplot() +
  facet_wrap(~Site)
#
#
# Contrasts - Vassijaure
contrasts(Q0_ecosys_stat_V$Round) <- Contr_Vassijaure_MP
#
# Check if contrasts work, by using a two-way ANOVA
SysModel_alias <- aov(sysRec ~ Round, data = Q0_ecosys_stat)
Anova(SysModel_alias, type ="III")
# alias checks dependencies
alias(SysModel_alias)
#
# transform data
Q0_ecosys_stat_V <-  Q0_ecosys_stat_V %>%
  mutate(Recov = sysRec)
Q0_ecosys_stat_V <- Q0_ecosys_stat_V %>%
  mutate(logRecov = log(Recov+1), # Good for low percentage values
         arcRecov = asin(sqrt(Recov/1000))) # Use is for this most transformations in percent. Recovery values are > 100, so transform by dividing by 1000 instead
#
# model:
lme0_V<-lme(arcRecov ~ Round,
            random = ~1|Plot,
            data = Q0_ecosys_stat_V, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme0_V), resid(lme0_V), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme0_V), main = "Normally distributed?")                 
qqline(resid(lme0_V), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme0_V)
par(mfrow = c(1,1))
#
# model output
Anova(lme0_V, type=2)
summary(lme0_V)
# Highly significant for Snow Cold vs Warm (t = 4.845361, p < 0.0001)
#
Q0_season_V <- Q0_ecosys_stat_V %>%
  select(1:4, sysRec) %>%
  mutate(Snow = if_else(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 | MP == 10 | MP == 11 | MP == 12 | MP == 13, "Snow","Clear"),
         SnowCW = case_when(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 ~ "Cold",
                            MP == 10 | MP == 11 | MP == 12 | MP == 13 ~ "Warm",
                            TRUE ~ NA),
         FreeWC = case_when(MP == 1 | MP == 2 | MP == 14 | MP == 15 ~ "Warm",
                            MP == 3 | MP == 4 ~ "Cold",
                            TRUE ~ NA))
#
#
Q0_season_V %>%
  filter(!is.na(SnowCW)) %>%
  summarise(sysRec = mean(sysRec), .by = c(SnowCW))
Q0_season_V %>%
  filter(!is.na(FreeWC)) %>%
  summarise(sysRec = mean(sysRec), .by = c(FreeWC))
#
#
Q0_season_V %>%
  summarise(sysRec_season = mean(sysRec), .by = c(Snow))
#
summarySE(Q0_season_V, measurevar = "sysRec", groupvars = c("Snow"))
summarySE(Q0_season_V, measurevar = "sysRec", groupvars = c("SnowCW"))
summarySE(Q0_season_V, measurevar = "sysRec", groupvars = c("FreeWC"))

Q0_season_V %>%
  ggplot(aes(x = Round, y = sysRec, fill = SnowCW)) +
  geom_boxplot() +
  facet_wrap(~Site)
#
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
Q1_veg_stat <- Rec15N %>%
  select(1:4, PlantRecovery, PlantR_frac) %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round"), as.factor))
#
# transform data
Q1_veg_stat <- Q1_veg_stat %>%
  mutate(Recov = PlantRecovery)
Q1_veg_stat <- Q1_veg_stat %>%
  mutate(logRecov = log(Recov+1), # Good for low percentage values. Use here as values are in the low end
         arcRecov = asin(sqrt(Recov/100))) # Use is for this most transformations in percent.
#
# model:
lme1<-lme(logRecov ~ Round*Site,
          random = ~1|Plot/Site,
          data = Q1_veg_stat, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1), resid(lme1), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1), main = "Normally distributed?")                 
qqline(resid(lme1), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1)
par(mfrow = c(1,1))
#
# model output
Anova(lme1, type=2)
# Log transformation
# Significant for Round (χ^2 = 27.671, p = 0.01573)
#
#
# Per site
Q1_veg_stat_A <- Q1_veg_stat %>%
  filter(Site == "Abisko")
Q1_veg_stat_V <- Q1_veg_stat %>%
  filter(Site == "Vassijaure")
#
#
# Contrasts Abisko
contrasts(Q1_veg_stat_A$Round) <- Contr_Abisko_MP
#
# transform data
Q1_veg_stat_A <- Q1_veg_stat_A %>%
  mutate(Recov = PlantRecovery)
Q1_veg_stat_A <- Q1_veg_stat_A %>%
  mutate(logRecov = log(Recov+1), # Good for low percentage values.
         arcRecov = asin(sqrt(Recov/100))) # General use is for this transformation.
#
# model:
lme1_A<-lme(logRecov ~ Round,
            random = ~1|Plot,
            data = Q1_veg_stat_A, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1_A), resid(lme1_A), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1_A), main = "Normally distributed?")                 
qqline(resid(lme1_A), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1_A)
par(mfrow = c(1,1))
#
# model output
Anova(lme1_A, type=2)
summary(lme1_A)
# Highly significant for Snow Cold vs Warm (t = 2.60138, p = 0.0119)
# Significant for Summer vs Autumn (t = -2.05424, p = 0.0446)
#
Q1_season_A <- Q1_veg_stat_A %>%
  select(1:4, PlantRecovery) %>%
  mutate(Snow = if_else(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 | MP == 10 | MP == 11 | MP == 12, "Snow","Clear"),
         SnowCW = case_when(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 ~ "Cold",
                            MP == 10 | MP == 11 | MP == 12 ~ "Warm",
                            TRUE ~ NA),
         FreeWC = case_when(MP == 1 | MP == 2 | MP == 13 | MP == 14 | MP == 15 ~ "Warm",
                            MP == 3 | MP == 4 ~ "Cold",
                            TRUE ~ NA))
#
Q1_season_A %>%
  summarise(PlantRecovery = mean(PlantRecovery), .by = c(Snow))
Q1_season_A %>%
  filter(!is.na(SnowCW)) %>%
  summarise(PlantRecovery = mean(PlantRecovery), .by = c(SnowCW))
Q1_season_A %>%
  filter(!is.na(FreeWC)) %>%
  summarise(PlantRecovery = mean(PlantRecovery), .by = c(FreeWC))
#
#
Q1_season_A %>%
  summarise(sysRec_season = mean(PlantRecovery), .by = c(Snow))
#
summarySE(Q1_season_A, measurevar = "PlantRecovery", groupvars = c("Snow"))
summarySE(Q1_season_A, measurevar = "PlantRecovery", groupvars = c("SnowCW"))
summarySE(Q1_season_A, measurevar = "PlantRecovery", groupvars = c("FreeWC"))
#
Q1_season_A %>%
  ggplot(aes(x = Round, y = PlantRecovery, fill = FreeWC)) +
  geom_boxplot() +
  facet_wrap(~Site)
#
#

# Contrasts Vassijaure
contrasts(Q1_veg_stat_V$Round) <- Contr_Vassijaure_MP
#
# transform data
Q1_veg_stat_V <- Q1_veg_stat_V %>%
  mutate(Recov = PlantRecovery)
Q1_veg_stat_V <- Q1_veg_stat_V %>%
  mutate(logRecov = log(Recov+1), # Good for low percentage values.
         arcRecov = asin(sqrt(Recov/100))) # General use is for this transformation.
#
# model:
lme1_V<-lme(logRecov ~ Round,
            random = ~1|Plot,
            data = Q1_veg_stat_V, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lme1_V), resid(lme1_V), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lme1_V), main = "Normally distributed?")                 
qqline(resid(lme1_V), main = "Homogeneity of Variances?", col = 2) #OK
plot(lme1_V)
par(mfrow = c(1,1))
#
# model output
Anova(lme1_V, type=2)
summary(lme1_V)
# Highly significant for Snow Cold vs Warm (t = 2.60138, p = 0.0119)
# Significant for Summer vs Autumn (t = -2.05424, p = 0.0446)
#
Q1_season_V <- Q1_veg_stat_V %>%
  select(1:4, PlantRecovery) %>%
  mutate(Snow = if_else(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 | MP == 10 | MP == 11 | MP == 12 | MP == 13, "Snow","Clear"),
         SnowCW = case_when(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 ~ "Cold",
                            MP == 10 | MP == 11 | MP == 12 | MP == 13 ~ "Warm",
                            TRUE ~ NA),
         FreeWC = case_when(MP == 1 | MP == 2 | MP == 14 | MP == 15 ~ "Warm",
                            MP == 3 | MP == 4 ~ "Cold",
                            TRUE ~ NA))
#
Q1_season_V %>%
  summarise(PlantRecovery = mean(PlantRecovery), .by = c(Snow))
Q1_season_V %>%
  filter(!is.na(SnowCW)) %>%
  summarise(PlantRecovery = mean(PlantRecovery), .by = c(SnowCW))
Q1_season_V %>%
  filter(!is.na(FreeWC)) %>%
  summarise(PlantRecovery = mean(PlantRecovery), .by = c(FreeWC))
#
#
Q1_season_V %>%
  summarise(sysRec_season = mean(PlantRecovery), .by = c(Snow))
#
summarySE(Q1_season_V, measurevar = "PlantRecovery", groupvars = c("Snow"))
summarySE(Q1_season_V, measurevar = "PlantRecovery", groupvars = c("SnowCW"))
summarySE(Q1_season_V, measurevar = "PlantRecovery", groupvars = c("FreeWC"))
#
Q1_season_V %>%
  ggplot(aes(x = Round, y = PlantRecovery, fill = FreeWC)) +
  geom_boxplot() +
  facet_wrap(~Site)
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
  left_join(vegroot15N_total_Plant, by = join_by(Site, Plot, MP, Round)) %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "Organ"), as.factor)) %>%
  select(1:7) %>%
  mutate(OrganRecovery_ofTotal = OrganRecovery, 
         OrganRecovery = OrganRecovery/PlantRecovery*100) %>%
  select(1:5,OrganRecovery)
#
# transform data
vegroot15N_Organ <- vegroot15N_Organ %>%
  mutate(logOrganRecovery = log(OrganRecovery+1), # Good for low percentage values.
         arcOrganRecovery = asin(sqrt((OrganRecovery)/100))) # Use is for this most transformations in percent.
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
#summary(lme1a)
# For organ recovery as fraction of whole plant recovery (Arcsin transformation):
# Highly significant for organ (χ^2 = 669.9609, p < 2.2*10^-16) 
# and all interactions (Round:Organ χ^2 = 190.0375, p < 2.2*10^-16 ; Site:Organ χ^2 = 16.2061, p = 0.0003026; Round:Site:Organ χ^2 = 64.6419, p = 0.0001006) except Round:Site
#
Q1a_season <- vegroot15N_Organ %>%
  select(1:4, Organ, OrganRecovery) %>%
  mutate(SeasonSW = case_when(MP == 1 | MP == 2 | MP == 13 | MP == 14 | MP == 15 ~ "Summer",
                              MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 ~ "Winter",
                              TRUE ~ NA),
         SeasonAS = case_when(MP == 3 | MP == 4 | MP == 5 ~ "Autumn",
                              MP == 10 | MP == 11 | MP == 12 ~ "Spring",
                              TRUE ~ NA),
         Snow = if_else(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 | MP == 10 | MP == 11 | MP == 12, "Snow", "Clear"))
#
# Significant for Summer vs Winter and Autumn vs Spring
Q1a_season %>%
  filter(!is.na(SeasonSW)) %>%
  summarise(PlantRecov_season = mean(OrganRecovery), .by = c(Organ, SeasonSW))
Q1a_season %>%
  filter(!is.na(SeasonAS)) %>%
  summarise(PlantRecov_season = mean(OrganRecovery), .by = c(SeasonAS))
# Summer 5.72 and winter 7.14. But only really visible for Abisko
# Autumn 6.56 and Spring 4.33, but larger for Abisko
#
# Significant for Snow covered season
Q1a_season %>%
  summarise(PlantRecov_season = mean(OrganRecovery), .by = c(Snow))
# Snow covered period: 6.09 and outside: 5.96
summarySE(Q1a_season, measurevar = "OrganRecovery", groupvars = c("Organ", "Snow"))
summarySE(Q1a_season, measurevar = "OrganRecovery", groupvars = c("Organ", "SeasonSW"))
summarySE(Q1a_season, measurevar = "OrganRecovery", groupvars = c("SeasonAS"))
summarySE(Q1a_season, measurevar = "OrganRecovery", groupvars = c("Organ", "Site"))
summarySE(Q1a_season, measurevar = "OrganRecovery", groupvars = c("Site", "MP", "Organ"))
#
Q1a_season %>%
  ggplot(aes(x = Round, y = OrganRecovery, fill = SeasonAS)) +
  geom_col() +
  facet_wrap(~Site)
#
Q1a_season %>%
  filter(Snow == "Clear") %>%
  summarise(Snow_sum = sum(OrganRecovery),
            Snow_count = length(OrganRecovery)) %>%
  mutate(Snow_avg = Snow_sum/Snow_count)
#
# Test correlation of biomass vs organ proportional recovery
vegroot15N_bio_organ <- vegroot15N %>%
  select(1:4,Species,Organ,Biomass_DW_g) %>%
  group_by(across(c("Site","Plot", "MP", "Round", "Organ"))) %>%
  summarise(Biomass = sum(Biomass_DW_g, na.rm = TRUE), .groups = "keep") %>%
  ungroup() %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "Organ"), as.factor))

vegroot15N_organBioCorr <- vegroot15N_Organ %>%
  select(1:4, Organ, OrganRecovery) %>%
  left_join(vegroot15N_bio_organ, by = join_by(Site, Plot, MP, Round, Organ))
vegroot15N_organBioCorr %>%
  ggplot(aes(x = log(Biomass), y = (OrganRecovery), color = Organ)) +
  geom_point() +
  geom_smooth(method = lm, span = 0.3, se = TRUE, alpha = 0.6) +
  #scale_color_viridis_d(labels = c("CR", "Plant", "TDN")) +
  facet_wrap( ~ Site)
#
CorrFunc <- function(xx) {
  return(data.frame(COR = cor(xx$Biomass, xx$OrganRecovery, method = "spearman")))
}
ddply(vegroot15N_organBioCorr, .(Organ, Site), CorrFunc)
# No correlation for most < 0.6
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
Q2_MBN_stat <- Rec15N %>%
  select("Site", "Plot", "MP", "Round", "R_MBN", "R_MBN_frac")
Q2_MBN_stat <- Q2_MBN_stat %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round"), as.factor))
#
# Transform data
Q2_MBN_stat <- Q2_MBN_stat %>%
  mutate(Recov = R_MBN_frac) 
Q2_MBN_stat <- Q2_MBN_stat %>%
  mutate(logRecov = log(Recov),
         arcRecov = asin(sqrt(Recov/max(100, max(Q2_MBN_stat$Recov)))))
#
#
#model:
lme2<-lme(arcRecov ~ Site*Round,
          random = ~1|Plot/Site,
          data = Q2_MBN_stat, na.action = na.exclude , method = "REML")
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
# Arcsin transformation
# Not significant for fraction
# For absolute recovery: Highly significant for Round (χ^2 = 39.8375, p = 0.0002704)
#
#
#
#=======  ###   Statistics   ### =======
#-------   ##       TDN      ##  -------
#
Q_TDN <- soil15N %>%
  left_join(Rec15N, by =  join_by(Site, Plot, MP)) %>%
  relocate(Round, .after = MP) %>%
  select(1:4, Nconc_microg_pr_gDW, NH4_microg_pr_gDW, NO3_microg_pr_gDW, R_TDN, R_TDN_frac) %>%
  mutate(across(c("Plot", "MP"), ~as.character(.x))) %>%
  mutate(across(c("Site", "MP", "Round"), ~as.factor(.x)))
#
# --- # TDN N concentration # ---
#
Q_TDN_Nconc <- Q_TDN %>%
  select(1:4, Nconc_microg_pr_gDW, NH4_microg_pr_gDW, NO3_microg_pr_gDW)
#
# transform data
Q_TDN_Nconc <- Q_TDN_Nconc %>%
  mutate(TDN = Nconc_microg_pr_gDW,
         NH4 = NH4_microg_pr_gDW,
         NO3 = NO3_microg_pr_gDW)
Q_TDN_Nconc <- Q_TDN_Nconc %>%
  mutate(logTDN = log(TDN),
         arcTDN = asin(sqrt(TDN/1000))) %>% # Max value is around 982 µg pr g
  mutate(logNH4 = log(NH4+1),
         arcNH4 = asin(sqrt((NH4+1)/1000))) %>%
  mutate(logNO3 = log(NO3+1),
         arcNO3 = asin(sqrt((NO3+1)/1000)))
#
# model:
lmeTDN_Nconc<-lme(logTDN ~ Round*Site,
            random = ~1|Plot/Site,
            data = Q_TDN_Nconc, na.action = na.exclude, method = "REML")
#
lmeTDN_NH4<-lme(logNH4 ~ Round*Site,
                  random = ~1|Plot/Site,
                  data = Q_TDN_Nconc, na.action = na.exclude, method = "REML")
#
lmeTDN_NO3<-lme(logNO3 ~ Round*Site,
                  random = ~1|Plot/Site,
                  data = Q_TDN_Nconc, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
# TDN
par(mfrow = c(1,2))
plot(fitted(lmeTDN_Nconc), resid(lmeTDN_Nconc), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lmeTDN_Nconc), main = "Normally distributed?")                 
qqline(resid(lmeTDN_Nconc), main = "Homogeneity of Variances?", col = 2) #OK
plot(lmeTDN_Nconc)
par(mfrow = c(1,1))
#
# model output
Anova(lmeTDN_Nconc, type=2)
# Not significant
#
# NH4
par(mfrow = c(1,2))
plot(fitted(lmeTDN_NH4), resid(lmeTDN_NH4), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lmeTDN_NH4), main = "Normally distributed?")                 
qqline(resid(lmeTDN_NH4), main = "Homogeneity of Variances?", col = 2) #OK
plot(lmeTDN_NH4)
par(mfrow = c(1,1))
#
# model output
Anova(lmeTDN_NH4, type=2)
#
# NO3
par(mfrow = c(1,2))
plot(fitted(lmeTDN_NO3), resid(lmeTDN_NO3), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lmeTDN_NO3), main = "Normally distributed?")                 
qqline(resid(lmeTDN_NO3), main = "Homogeneity of Variances?", col = 2) #OK
plot(lmeTDN_NO3)
par(mfrow = c(1,1))
#
# model output
Anova(lmeTDN_NO3, type=2)

Q_TDN_Nconc %>% ggplot(aes(x = Round, y = NO3)) + geom_boxplot()
#
#
# --- # TDN 15N recovery # ---
# 
# Not a part of the direct questions, but still related
# Model
# Response variable: Dissolved Nitrogen 15N pool of injected. TDN
# Factors: Time, Site, Time*Site
# Most important factor: Time
#
Q_TDN_stat <- Q_TDN %>%
  select(1:4, R_TDN, R_TDN_frac)
#
# transform data
Q_TDN_stat <- Q_TDN_stat %>%
  mutate(Recov = R_TDN_frac)
Q_TDN_stat <- Q_TDN_stat %>%
  mutate(logRecov = log(Recov),
         arcRecov = asin(sqrt(Recov/100)))
#
# model:
lmeTDN<-lme(logRecov ~ Round*Site,
          random = ~1|Plot/Site,
          data = Q_TDN_stat, na.action = na.exclude, method = "REML")
#
# Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lmeTDN), resid(lmeTDN), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lmeTDN), main = "Normally distributed?")                 
qqline(resid(lmeTDN), main = "Homogeneity of Variances?", col = 2) #OK
plot(lmeTDN)
par(mfrow = c(1,1))
#
# model output
Anova(lmeTDN, type=2)
# Highly significant for Round (χ^2 = 87.5184, p = 1.114e-12)
#
#
#
#=======  ###   Statistics   ### =======
#-------   ##  Competition   ##  -------
#
# Model
# Response variable: Microbial vs Plant vs TDN recovery of 15N (%)
# Factors: Time, Site, Part, Time*Site, Time*Part, Site*Part, Time*Site*Part
# Most important factor: Part and interaction
#
# Load data from excel instead of calculated combined
Compet_all <- Rec15N %>%
  select("Site", "Plot", "MP", "Round", "PlantR_frac", "R_MBN_frac", "R_TDN_frac") %>%
  pivot_longer(5:7, names_to = "sysPart", values_to = "Recov")
Compet_all <- Compet_all %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "sysPart"), as.factor))
#
# Contrasts - Part recovery
#              (Pl,Mic,TDN)
PlvsMic    <- c( 1, -1,  0)
PlvsMicTDN <- c( 2, -1, -1)
MicvsTDN   <- c( 0,  1, -1)
#
contrasts(Compet_all$sysPart)<-cbind(PlvsMic, PlvsMicTDN, MicvsTDN)
crossprod(cbind(PlvsMic, PlvsMicTDN, MicvsTDN))
#
# Month           (J, A, S, O, N, D, J, F, M, A, A, M, J, J, A) # Two times April
# MP              (1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15)
#
# Same contrasts for Site and Round
contrasts(Compet_all$Site)<-AvsV
# Same contrasts as for plant recovery:
# SummervsWinter, SpringvsAutumn, SnowvsNot
contrasts(Compet_all$Round)<-cbind(SummervsWinter,SpringvsAutumn,SnowvsNot, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh, cont11, cont12, cont13, cont14)
#contrasts(Compet_all$Round)<-contr.helmert
#
# Check if contrasts work
# Check contrasts are orthogonal
crossprod(cbind(SummervsWinter,SpringvsAutumn,SnowvsNot, JulvsJan, OctvsApr, Summervs2, SpringChA, SpringChV, AutumnCh, WinterCh))#, cont11, cont12, cont13, cont14))
#
# Alternative: Two-way ANOVA
# have time as a factor in a two-way ANOVA, combined with Site. As each sampling is destructive, the samples are technically independent of each other, although it does not account for the block design
CompetModel3 <- aov(Recov ~ Round*Site, data = Compet_all)
Anova(CompetModel3, type ="III")
#
#
# Transform data
Compet_all <- Compet_all %>%
  mutate(sqrtRecov = sqrt(Recov),
         logRecov = log(Recov),
         cubeRecov = (Recov)^(1/3),
         arcRecov = asin(sqrt(Recov/100)))
#
Compet_all %>%
  mutate(across("MP", ~ as.numeric(.x))) %>%
  ggplot(aes(x = MP, y = Recov, fill = sysPart)) +
  geom_point()
#
#
#model:
lmeCompet<-lme(Recov ~ Site*Round*sysPart,
          random = ~1|Plot/Site,
          data = Compet_all, na.action = na.exclude , method = "REML")
#
#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lmeCompet), resid(lmeCompet), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lmeCompet), main = "Normally distributed?")                 
qqline(resid(lmeCompet), main = "Homogeneity of Variances?", col = 2) #OK
plot(lmeCompet)
par(mfrow = c(1,1))
#
#model output
Anova(lmeCompet, type=2)
summary(lmeCompet)
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
  left_join(vegroot15N_total_Plant, by = join_by(Site, Plot, MP, Round)) %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "Organ"), as.factor)) %>%
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
Mic15N_sum <- summarySE(Rec15N, measurevar="R_MBN", groupvars=c("Site", "Round"), na.rm=TRUE)
TDN15N_sum <- summarySE(Rec15N, measurevar="R_TDN", groupvars=c("Site", "Round"), na.rm=TRUE)
sysRec_sum <- summarySE(Rec15N, measurevar="sysRec", groupvars=c("Site", "Round"), na.rm=TRUE)
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
# Total ecosystem recovery +/- 95% CI
sysRec_sum %>%  
  ggplot() + 
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP2$wstart, xmax=winterP2$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_errorbar(aes(x = Round, y = sysRec, ymin=sysRec-ci, ymax=sysRec+ci), position=position_dodge(.9)) +
  geom_col(aes(Round, sysRec),color = "black") +
  scale_x_discrete(labels = measuringPeriod) +
  facet_wrap( ~ Site, ncol = 2, scales = "free") + 
  labs(x = "Measuring period (MP)", y = expression("% of added "*{}^15*"N"), title = expression("Total ecosystem "*{}^15*"N tracer recovery")) + 
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(2, "lines"),axis.text.x=element_text(angle=60, hjust=1))
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
#
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
#
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
  rename(R_Plant_frac = "PlantR_frac") %>%
  pivot_longer(cols = c(R_Plant_frac, R_MBN_frac, R_TDN_frac), names_to = "Type", values_to = "Recovery") %>%
  select(1:4, Type, Recovery) %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x)))
#
Rec15N_3 <- Rec15N %>%
  rename(R_Plant = "PlantRecovery") %>%
  pivot_longer(cols = c(R_Plant, R_MBN, R_TDN), names_to = "Type", values_to = "Recovery") %>%
  select(1:4, Type, Recovery) %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x)))
#
# Plot with geom_smooth
Rec15N_2 %>%
  ggplot(aes(x = Day_of_harvest, y = Recovery, color = Type)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.4, fill = 'grey', inherit.aes = FALSE) +
  geom_point() +
  geom_smooth(span = 0.3, se = TRUE, alpha = 0.6) +
  scale_color_viridis_d(labels = c("Microbial", "Plant", "TDN")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  labs(x = "Time of harvest", y = expression("% of total system recovered "*{}^15*"N"), title = expression("Plant, microbial, and TDN "*{}^15*"N tracer recovery per part of the system")) +
  facet_wrap( ~ Site) +
  theme_classic(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))
#
Rec15N_3 %>%
  ggplot(aes(x = Day_of_harvest, y = Recovery, color = Type)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.4, fill = 'grey', inherit.aes = FALSE) +
  geom_point() +
  geom_smooth(span = 0.3, se = TRUE, alpha = 0.6) +
  scale_color_viridis_d(labels = c("Microbial", "Plant", "TDN")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  labs(x = "Time of harvest", y = expression("% of added "*{}^15*"N"), title = expression("Plant, microbial, and TDN "*{}^15*"N tracer recovery")) +
  facet_wrap( ~ Site) +
  theme_classic(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))
#
#
#
#-------   ## Environmental  ## -------
#
coreData_2 <- coreData %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x)))

coreData_snow <- summarySE(coreData, measurevar = "Snow_depth_plot_cm", groupvars = c("Site", "Round"))
coreData_snow <- coreData_snow %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x)))
#

# Soil mass
coreData %>%
  ggplot(aes(x = Round, y = Soil_RF_DW_g)) +
  geom_boxplot() +
  labs(x = "Time of harvest", y = "Root free soil g DW", title = "Soil mass for microbial") +
  facet_wrap(~Site, scales = "free") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))
#
# Core depth
# Boxplot
coreData %>%
  ggplot(aes(x = Round, y = -Soil_depth_cm)) +
  geom_boxplot() +
  coord_cartesian(y = c(-11,0)) + 
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10), labels = abs) +
  labs(x = "Time of harvest", y = "Depth of corer (cm)", title = "Corer depth as measured in the field at sampling") +
  facet_wrap(~Site, scales = "free") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))
#
# Point data
coreData %>%
  ggplot(aes(x = Round, y = -Soil_depth_cm)) +
  geom_point() +
  coord_cartesian(y = c(-11,0)) + 
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10), labels = abs) +
  labs(x = "Time of harvest", y = "Depth of corer (cm)", title = "Corer depth as measured in the field at sampling") +
  facet_wrap(~Site, scales = "free") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))
#
#
# Snow cover
coreData %>%
  ggplot(aes(x = Round, y = Snow_depth_plot_cm)) +
  geom_boxplot() +
  labs(x = "Time of harvest", y = "Snow cover (cm)", title = "Snow cover measured over the entire plot (around all 15 patches)") +
  facet_wrap(~Site, scales = "free") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))


coreData_snow %>%
  ggplot(aes(x = Day_of_harvest, y = Snow_depth_plot_cm, ymin = Snow_depth_plot_cm-ci, ymax = Snow_depth_plot_cm+ci, fill = Site)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_grey() +
  #scale_fill_viridis_d() +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  labs(x = "Time of harvest", y = "Snow cover (cm)", title = "Snow cover measured over the entire plot (around all 15 patches)") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))


#
# Select only snow-covered period: MP5-MP13 (Vassijaure end)
coreData_SnowP <- coreData %>%
  filter(MP == 5 | MP == 6 | MP == 7 | MP == 8 | MP == 9 | MP == 10 | MP == 11 | MP == 12)# | MP == 13)
lmeSnow<-lme(log(Snow_depth_plot_cm + 1) ~ Site*Round,
             random = ~1|Plot/Site,
             data = coreData_SnowP, na.action = na.exclude , method = "REML")
#
#Checking assumptions:
par(mfrow = c(1,2))
plot(fitted(lmeSnow), resid(lmeSnow), 
     xlab = "fitted", ylab = "residuals", main="Fitted vs. Residuals") 
qqnorm(resid(lmeSnow), main = "Normally distributed?")                 
qqline(resid(lmeSnow), main = "Homogeneity of Variances?", col = 2) #OK
plot(lmeSnow)
par(mfrow = c(1,1))
#
#model output
Anova(lmeSnow, type=2)
summary(lmeSnow)
#
#
# Soil moisture
coreData %>%
  ggplot(aes(x = Round, y = (SM_FW_g - SM_DW_g)/SM_DW_g*100)) +
  geom_boxplot() +
  labs(x = "Time of harvest", y = "Soil moisture (g pr g FW)", title = "Soil moisture") +
  facet_wrap(~Site, scales = "free") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))
#
#
#
# Testing ##----
#
vegroot15N_Organ_sum2 <- vegroot15N_Organ_sum %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, Organ, OrganRecovery, ci, Day_of_harvest)
#
vegroot15N_Organ_sum2 %>%
  ggplot(aes(x = Day_of_harvest, y = OrganRecovery, ymin = OrganRecovery-ci, ymax = OrganRecovery+ci, fill = Organ)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.4, fill = 'grey', inherit.aes = FALSE) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  #coord_cartesian(ylim = c(-110,75)) +
  scale_fill_viridis_d(labels = c("CR", "FR", "S")) +
  #scale_linetype(labels = c("Abisko", "Vassijaure")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  facet_wrap( ~ Site, ncol = 2, scales = "free") +
  labs(x = "Measuring period (MP)", y = expression("% of total plant recovered "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) +
  #guides(color = guide_legend(title = "Plant organ")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))

  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), labels = abs)
#
vegroot15N_Organ_sum2 %>%
  ggplot(aes(x = Day_of_harvest, y = OrganRecovery, ymin = OrganRecovery-ci, ymax = OrganRecovery+ci, fill=Organ, linetype=Site)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("CR", "FR", "S")) +
  scale_linetype(labels = c("Abisko", "Vassijaure")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  labs(x = "Time of harvest", y = expression("% of total plant recovered "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x=element_text(angle=60, hjust=1))
#
vegroot15N_Organ_sum2 %>%
  ggplot(aes(x = Day_of_harvest, y = OrganRecovery, ymin = OrganRecovery-ci, ymax = OrganRecovery+ci, fill=Organ, linetype=Organ)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("CR", "FR", "S")) +
  scale_linetype(labels = c("CR", "FR", "S")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))+
  labs(x = "Time of harvest", y = expression("% of total plant recovered "*{}^15*"N"), title = expression("Plant "*{}^15*"N tracer recovery")) +
  facet_wrap( ~ Site, ncol = 2) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))
#
# Biomass
test_biom <- vegroot15N %>%
  select(1:4, Species, Organ, Biomass_DW_g)
test_biom_sum <- summarySE(test_biom, measurevar = "Biomass_DW_g", groupvars = c("Site", "Round", "Organ"), na.rm = TRUE)
test_biom_sum <- test_biom_sum %>%
  left_join(DayOf, by = join_by(Site, Round)) %>%
  mutate(across(Day_of_harvest, ~ as.Date(.x))) %>%
  select(Site, Round, Organ, Biomass_DW_g, ci, Day_of_harvest)
#
test_biom_sum %>%
  ggplot(aes(x = Day_of_harvest, y = Biomass_DW_g, ymin = Biomass_DW_g-ci, ymax = Biomass_DW_g+ci, fill=Organ, linetype=Organ)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("CR", "FR", "S")) +
  scale_linetype(labels = c("CR", "FR", "S")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8))+
  labs(x = "Time of harvest", y = "Biomass (g)", title = "Plant biomass") +
  facet_wrap( ~ Site, ncol = 2) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))
#
# No CR
test_biom_sum %>%
  filter(Organ != "CR") %>%
  ggplot(aes(x = Day_of_harvest, y = Biomass_DW_g, ymin = Biomass_DW_g-ci, ymax = Biomass_DW_g+ci, fill=Organ, linetype=Organ)) +
  geom_rect(data=data.frame(variable=factor(1)), aes(xmin=winterP_date$wstart, xmax=winterP_date$wend, ymin=-Inf, ymax=Inf), alpha = 0.5, fill = 'grey', inherit.aes = FALSE) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_viridis_d(labels = c("FR", "S")) +
  scale_linetype(labels = c("FR", "S")) +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%Y-%b-%d") +
  scale_y_continuous(breaks = c(0, 1, 2, 3))+
  labs(x = "Time of harvest", y = "Biomass (g)", title = "Plant biomass") +
  facet_wrap( ~ Site, ncol = 2) +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))

#
# Trying out biomass vs recovery
test2 <- vegroot15N %>%
  select(1:4,Species,Organ,Biomass_DW_g) %>%
  group_by(across(c("Site","Plot", "MP", "Round", "Organ"))) %>%
  summarise(Biomass = sum(Biomass_DW_g, na.rm = TRUE), .groups = "keep") %>%
  ungroup() %>%
  mutate(across(c("Plot", "MP"), as.character))%>%
  mutate(across(c("Site", "MP", "Round", "Organ"), as.factor))

test <- vegroot15N_Organ %>%
  select(1:4, Organ, OrganRecovery) %>%
  left_join(test2, by = join_by(Site, Plot, MP, Round, Organ))
test %>%
  ggplot(aes(x = log(Biomass), y = (OrganRecovery), color = Organ)) +
  geom_point() +
  geom_smooth(method = lm, span = 0.3, se = TRUE, alpha = 0.6) +
  #scale_color_viridis_d(labels = c("CR", "Plant", "TDN")) +
  facet_wrap( ~ Site)
#
# TDN concentration changes
soil15N %>%
  left_join(coreData, by = join_by(Site,Plot,MP)) %>%
  select(1:11, Round) %>%
  ggplot(aes(x = Round, y = Nconc_microg_pr_gDW, fill = Extr_type)) +
  geom_boxplot() +
  labs(x = "Time of harvest", y = "TDN µg pr g DW", title = "TDN") +
  facet_wrap(~Site, scales = "free") +
  theme_classic(base_size = 20) +
  theme(panel.spacing = unit(1, "lines"), axis.text.x=element_text(angle=60, hjust=1))
#
#
#
#=======  ###  { The End }   ### =======