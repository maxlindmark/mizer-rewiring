#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. Analyisis of changes in yield with fishing warming
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(ggplot2)
library(devtools)
library(RColorBrewer)
library(RCurl)
library(magrittr)
library(viridis)
library(tidyr)
library(dplyr)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Install and reload local mizer package
devtools::load_all(".")

# Print package versions
# print(sessionInfo())
# other attached packages:
# mizer_1.1          testthat_2.3.0     patchwork_0.0.1    dplyr_0.8.3        tidyr_1.0.0        
# viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       RCurl_1.95-4.12   
# bitops_1.0-6       RColorBrewer_1.1-2 devtools_2.2.1     usethis_1.5.1      ggplot2_3.2.1  

# Load function for extracting size-at-age
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getSpectra.R", ssl.verifypeer = FALSE)
eval(parse(text = func))

#**** Read in parameters and data ==================================================
# Read in params object
params <- readRDS("baltic/params/mizer_param_calib.rds")

# Read in params object
ea <- read.csv("baltic/params/samples_activation_energy.csv")#[, 2:6]
ea <- ea %>% dplyr::rename("car" = "X.gro")

# Read in effort and temperature for projections
projectEffort <- read.csv("baltic/params/projectEffort.csv")[, 2:4]
projectTemp <- read.csv("baltic/params/projectTemp.csv")

projectEffort_m <- as.matrix(projectEffort)
rownames(projectEffort_m) <- 1:nrow(projectEffort)

# Define general parameters
dt <- 0.2
t_ref <- params@t_ref
kappa_ben <- params@kappa_ben
kappa <- params@kappa
w_bb_cutoff <- 20 # Not stored in mizerParams output
w_pp_cutoff <- 1 # Not stored in mizerParams outputs
r_pp <- 4 # Not stored in mizerParams output
r_bb <- 4 # Not stored in mizerParams output
t_max <- 100 # For FMSY plot!


#**** Update species params ========================================================
t <- params@species_params

t$ea_met <- mean(ea$met)
t$ea_int <- mean(ea$int)
t$ea_mor <- mean(ea$mor)

#t$ca_int <- -0.004 # Here we just use the fixed values
#t$ca_met <- 0.001 # Here we just use the fixed values

pars_no_res <- MizerParams(t, 
                           ea_gro = 0,
                           ea_car = 0, # -ea$gro[i] 
                           kappa_ben = kappa_ben,
                           kappa = kappa,
                           w_bb_cutoff = w_bb_cutoff,
                           w_pp_cutoff = w_pp_cutoff,
                           r_pp = r_pp,
                           r_bb = r_bb,
                           t_ref = t_ref)

pars_with_res <- MizerParams(t, 
                             ea_gro = mean(ea$gro),
                             ea_car = mean(ea$car), # -ea$gro[i] 
                             kappa_ben = kappa_ben,
                             kappa = kappa,
                             w_bb_cutoff = w_bb_cutoff,
                             w_pp_cutoff = w_pp_cutoff,
                             r_pp = r_pp,
                             r_bb = r_bb,
                             t_ref = t_ref)

pars_with_res_barnes <- MizerParams(t, 
                                    ea_gro = 0,
                                    ea_car = mean(ea$b_car),
                                    kappa_ben = kappa_ben,
                                    kappa = kappa,
                                    w_bb_cutoff = w_bb_cutoff,
                                    w_pp_cutoff = w_pp_cutoff,
                                    r_pp = r_pp,
                                    r_bb = r_bb,
                                    t_ref = t_ref)

# No physiological scaling scenarios:
t_no_ea <- t
t_no_ea$ea_int <- 0
t_no_ea$ea_met <- 0
t_no_ea$ea_mor <- 0

pars_with_res_barnes_np <- MizerParams(t_no_ea, 
                                       ea_gro = 0,
                                       ea_car = mean(ea$b_car),
                                       kappa_ben = kappa_ben,
                                       kappa = kappa,
                                       w_bb_cutoff = w_bb_cutoff,
                                       w_pp_cutoff = w_pp_cutoff,
                                       r_pp = r_pp,
                                       r_bb = r_bb,
                                       t_ref = t_ref)

pars_with_res_np <- MizerParams(t_no_ea, 
                                ea_gro = mean(ea$gro),
                                ea_car = mean(ea$car),
                                kappa_ben = kappa_ben,
                                kappa = kappa,
                                w_bb_cutoff = w_bb_cutoff,
                                w_pp_cutoff = w_pp_cutoff,
                                r_pp = r_pp,
                                r_bb = r_bb,
                                t_ref = t_ref)


# Define temperature-scenarios
consTemp <- projectTemp$temperature
start <- 1997-1914
consTemp[start:137] <- t_ref


# B. TEMP-DRIVEN CHANGE IN FISHERIES YIELD =========================================
# Here I will plot projected yield over time from time-varying effort and 
# temperature, i.e. not as in FMSY plot where we do a "bifurcation" plot.
# A more thorough exploration of the effect of temperature and fishing will
# be done using constant effort and temperature (not time-varying) - similar to how FMSY 
# is calculated. See scripts in baltic/R/analysis/constant_temperature_effort

# With temperature with no effects on the resource
m_no_resource_temp <- project(pars_no_res, 
                              dt = dt,
                              effort = projectEffort_m,
                              temperature = projectTemp$temperature,
                              diet_steps = 10) 

# With temperature with effects on the resource
m_with_resource_temp <- project(pars_with_res, 
                                dt = dt,
                                effort = projectEffort_m,
                                temperature = projectTemp$temperature,
                                diet_steps = 10) 

# With temperature with effects on the resource (BARNES)
m_with_resource_temp_b <- project(pars_with_res_barnes, 
                                 dt = dt,
                                 effort = projectEffort_m,
                                 temperature = projectTemp$temperature,
                                 diet_steps = 10) 

# No physiology: with temperature with effects on the resource
m_with_resource_temp_np <- project(pars_with_res_np, 
                                   dt = dt,
                                   effort = projectEffort_m,
                                   temperature = projectTemp$temperature,
                                   diet_steps = 10) 

# No physiology: With temperature with effects on the resource (BARNES)
m_with_resource_temp_b_np <- project(pars_with_res_barnes_np, 
                                     dt = dt,
                                     effort = projectEffort_m,
                                     temperature = projectTemp$temperature,
                                     diet_steps = 10) 

# Constant temperatures after calibration period (with res in pre calibration period)
m_cons_temp <- project(pars_with_res, 
                       dt = dt,
                       effort = projectEffort_m,
                       temperature = consTemp,
                       diet_steps = 10) 


# Predicted yield - no temp dep resource
no_resource_temp <- data.frame(getYield(m_no_resource_temp))
no_resource_temp$Year_ct <- as.numeric(rownames(getYield(m_no_resource_temp)))
no_resource_temp$Year <- no_resource_temp$Year_ct + (1914-1) # 1914 is so that 60 year burn-in leads to start at 1974
no_resource_temp$Scenario <- "Physio."

# Predicted yield - with temp on resource
with_resource_temp <- data.frame(getYield(m_with_resource_temp))
with_resource_temp$Year_ct <- as.numeric(rownames(getYield(m_with_resource_temp)))
with_resource_temp$Year <- with_resource_temp$Year_ct + (1914-1) # 1914 is so that 60 year burn-in leads to start at 1974
with_resource_temp$Scenario <- "Physio. + Resource (exp.)"

# Predicted yield - with temp on resource (Barnes)
with_resource_temp_b <- data.frame(getYield(m_with_resource_temp_b))
with_resource_temp_b$Year_ct <- as.numeric(rownames(getYield(m_with_resource_temp_b)))
with_resource_temp_b$Year <- with_resource_temp_b$Year_ct + (1914-1) # 1914 is so that 60 year burn-in leads to start at 1974
with_resource_temp_b$Scenario <- "Physio. + Resource (obs.)"

# No phys: Predicted yield - with temp on resource
with_resource_temp_np <- data.frame(getYield(m_with_resource_temp_np))
with_resource_temp_np$Year_ct <- as.numeric(rownames(getYield(m_with_resource_temp_np)))
with_resource_temp_np$Year <- with_resource_temp_np$Year_ct + (1914-1) # 1914 is so that 60 year burn-in leads to start at 1974
with_resource_temp_np$Scenario <- "Resource (exp.)"

# No phys: Predicted yield - with temp on resource (Barnes)
with_resource_temp_b_np <- data.frame(getYield(m_with_resource_temp_b_np))
with_resource_temp_b_np$Year_ct <- as.numeric(rownames(getYield(m_with_resource_temp_b_np)))
with_resource_temp_b_np$Year <- with_resource_temp_b_np$Year_ct + (1914-1) # 1914 is so that 60 year burn-in leads to start at 1974
with_resource_temp_b_np$Scenario <- "Resource (obs.)"

# Predicted yield - constant temperature at all
cons_temp <- data.frame(getYield(m_cons_temp))
cons_temp$Year_ct <- as.numeric(rownames(getYield(m_cons_temp)))
cons_temp$Year <- cons_temp$Year_ct + (1914-1) # 1914 is so that 60 year burn-in leads to start at 1974
cons_temp$Scenario <- "No warming"

# Combine
yield <- rbind(no_resource_temp, 
               with_resource_temp, 
               with_resource_temp_b,
               with_resource_temp_np, 
               with_resource_temp_b_np,
               cons_temp)

# Convert to long data frame (1 obs = 1 row)
# The first year with real effort is 1974. This is year 61 with centered time (1974-1914 +1),
# The last year before FMSY is 2012. In centered time this is 2012-1914 = 98
yield_l <- yield %>% 
  filter(Year > 1973) %>%
  gather(Species, yield_g.m2, 1:3)

# Scale up from g/m2 to 10^6 kg / Baltic
yield_l$Yield <- yield_l$yield_g.m2 * 2.49e+11 / (1e9) # See calibration document
yield_l <- yield_l %>% select(-yield_g.m2)

# Plot predicted and observed yield by species, normalize by max within species
pal <- rev(RColorBrewer::brewer.pal(n = 5, "Dark2"))
pal2 <- c("black", pal)

yield_l %>% 
  #filter(Year > 2002) %>% 
  #filter(Year > 2015) %>% 
  filter(Year > 2010) %>% 
  ggplot(., aes(Year, Yield, color = Scenario, linetype = Scenario, size = Scenario)) +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = pal2) +
  #scale_linetype_manual(values = c(2,1,1,1,1,1)) +
  scale_size_manual(values = c(0.5,1.3,1.3,1.3,1.3,1.3)) +
  theme(aspect.ratio = 1) +
  labs(y = "Yield (1000 tonnes/year)", x = "Year") +
  guides(linetype = FALSE) +
  theme_classic(base_size = 14) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 1/2) +
  NULL

#ggsave("baltic/figures/time_series_pred_yield.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")
