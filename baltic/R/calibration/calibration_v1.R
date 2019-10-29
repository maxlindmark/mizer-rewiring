#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.04.16: Max Lindmark
#
# Code for calibrating the default a Baltic Sea mizer model (1992-2002). This is a 
# standard calibration method, similar to Blanchard et al 2014 J.A.E. 
# 
# A. Load parameters and set up objects
#
# B. Load data
#
# C. Calibrate model to observed SSB, validate:
#    ssb, growth, recruitment, diet (see calibration protocol)
#
# D. Validate further by projecting and evaluate fit through time
#
#
#-Move below to analysis-script?
# Or where do I split?
# When I validate with time series I need temperature, so I should probably:
# 1. Add temperature of resource params (change code)
# 2. Read in temperature data in the Load data section
# 3. Do model validation with and without temperature, i.e. borrow code from below. 
#    Comparing correlations with and without temperature will be good enough as validatio
# 3. Analysis will be any alterations to the default setup. This code gets a default setup
#
# E. Project under climate & fishing scenarios
#
# F. Vary effort given warming scenarios
# 
# G. Test time (temperature)-varying kappa
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Questions in progress marked with ---

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# When doing a fresh start I need to check I'm in the right libpath to get the right mizer version
# .libPaths()
# .libPaths("C:/Program Files/R/R-3.5.0/library")

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
# [1] mizer_1.1 testthat_2.0.0 patchwork_0.0.1 dplyr_0.8.1 tidyr_0.8.3       
# [6] viridis_0.5.1 viridisLite_0.3.0 magrittr_1.5 RCurl_1.95-4.12 bitops_1.0-6      
# [11] RColorBrewer_1.1-2 usethis_1.4.0 devtools_2.0.2 ggplot2_3.1.1  


# B. READ DATA =====================================================================
#**** Species parameters ===========================================================
balticParams <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/params/species_params.csv"), sep = ";")

# Fix data format after reading
str(balticParams)
cols = c(2:4, 8) # Baltic area need special care...
balticParams[, cols] %<>% lapply(function(x) as.numeric(as.character(x)))
str(balticParams)

# Add area of Baltic (roughly Baltic proper) as a column
balticParams$sd25.29.32_m.2 <- 2.49e+11


#**** Empirical growth =============================================================
# Create a dataframe to hold empirical VBGE growth for comparison
age <- 0:20

empiri <- data.frame(
  Species = rep(c("Cod", "Sprat", "Herring"), each = 21),
  age = rep(age, length(unique(balticParams$species))),
  k_vb = rep(balticParams$k_vb, each = 21),
  w_inf = rep(balticParams$w_inf, each = 21),
  t_0 = rep(c(-0.936095, -2.97395, -3.739548), each = 21), # see VBGE script
  b = rep(c(3.07, 3.15, 3.14), each = 21),                 # see Appendix
  w_mat = rep(balticParams$w_mat, each = 21)) 

# Calculate weight at age
empiri$weight <- empiri$w_inf*(1-exp(-empiri$k_vb*(empiri$age - empiri$t_0)))^empiri$b

# ggplot(empiri, aes(age, weight)) + 
#   geom_line(size = 2) + 
#   facet_wrap(~Species, scales = "free_y", ncol = 3) +
#   geom_hline(data = empiri, aes(yintercept = w_mat), color = "gray30") +
#   theme_classic(base_size = 18)

# Read in growth data (currently only local since too big for github)
vbgedat <- read.csv("C:/R_STUDIO_PROJECTS/mizer-rewiring-baltic/baltic/data/BITS/clean_BITS.csv", sep = ",")
str(vbgedat)

vbgedat <- vbgedat %>% filter(Species %in% c("Cod", "Herring", "Sprat"))

vbgedat$Weight_g <- NA 

vbgedat$Weight_g <- ifelse(vbgedat$Species == "Herring",
                           0.0042*vbgedat$Length_cm^3.14, 
                           vbgedat$Weight_g)

vbgedat$Weight_g <- ifelse(vbgedat$Species == "Sprat",
                           0.0041*vbgedat$Length_cm^3.15, 
                           vbgedat$Weight_g)

vbgedat$Weight_g <- ifelse(vbgedat$Species == "Cod",
                           0.0078*vbgedat$Length_cm^3.07, 
                           vbgedat$Weight_g)

# Filter data to only predict within max age of species
maxcodage <- vbgedat %>% filter(Species == "Cod") %>% filter(AgeRings == max(AgeRings))
maxherage <- vbgedat %>% filter(Species == "Herring") %>% filter(AgeRings == max(AgeRings))
maxsprage <- vbgedat %>% filter(Species == "Sprat") %>% filter(AgeRings == max(AgeRings))

vbgedat$keep3 <- "N"

vbgedat$keep3 <- ifelse(vbgedat$Species == "Cod" & vbgedat$AgeRings < maxcodage$AgeRings[1],
                        "Y",
                        vbgedat$keep3)

vbgedat$keep3 <- ifelse(vbgedat$Species == "Herring" & vbgedat$AgeRings < maxherage$AgeRings[1],
                        "Y",
                        vbgedat$keep3)

vbgedat$keep3 <- ifelse(vbgedat$Species == "Sprat" & vbgedat$AgeRings < maxsprage$AgeRings[1],
                        "Y",
                        vbgedat$keep3)


#**** Stock assessment data ========================================================
ssb_f <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/data/SSB_F/SSB_F_data.csv"), sep = ";")

head(ssb_f)
str(ssb_f)

colnames(ssb_f)[5] <- "Fm" # Bad idea to call a column F in R

# Create dataframes with mean SSB and F for calibration period plotting
min_cal_yr <- 1992
max_cal_yr <- 2002

ref_time <- data.frame(Year = c(min_cal_yr, max_cal_yr),
                       TSB  = c(0, max(ssb_f$TSB, na.rm = T)),
                       SSB  = c(0, max(ssb_f$SSB, na.rm = T)),
                       Fm   = c(0, max(ssb_f$Fm,  na.rm = T)))

mean_ssb_F <- ssb_f %>% 
  dplyr::filter(Year >= min_cal_yr & Year <= max_cal_yr) %>% 
  dplyr::group_by(Species) %>% 
  dplyr::mutate(mean_SSB = mean(SSB),
                mean_F   = mean(Fm))

unique(mean_ssb_F$mean_F)
unique(mean_ssb_F$mean_SSB)


# C. CALIBRATE MODEL ================================================================
# This will help seeing the lines..
update_geom_defaults("line", list(size = 1.75))
col <- colorRampPalette(brewer.pal(5, "Dark2"))(5)

# Set defaults for regeneration and size ranges of background resources
r_pp <-  4
r_bb <-  4
w_bb_cutoff <- 20
w_pp_cutoff <- 1

# Set scalars on availability of background resources
balticParams$avail_PP[which(balticParams$species == "Cod")] <- 0.5
balticParams$avail_PP[which(balticParams$species == "Herring")] <- 0.5
balticParams$avail_PP[which(balticParams$species == "Sprat")] <- 1

balticParams$avail_BB[which(balticParams$species == "Cod")] <- 0.5
balticParams$avail_BB[which(balticParams$species == "Herring")] <- 0.5
balticParams$avail_BB[which(balticParams$species == "Sprat")] <- 0

# Refine predation parameters
# Using Reum et al 2019 Oikos parameters for forage fish, but specific value for cod (this study)
balticParams$beta[1] <- 426
balticParams$beta[2] <- 1000
balticParams$beta[3] <- 1000
balticParams$sigma[1] <- 5.6
balticParams$sigma[2] <- 1
balticParams$sigma[3] <- 1

# Set allometric erepro
balticParams$erepro <- 0.05 * balticParams$w_inf^(-0.5)

# Create effort vector
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])


#** 1. Find starting value for kappa ===============================================
# Given the default model, what should kappa be to get ssb in the same order of magnitude? This is just an iterative process, since well opimize r_max to minimize residual sum of squares (RSS) between predicted and observed SSB
dt <- 0.1

# These are the values I choose (lowest kappa with coexistence).
kappa_ben <- 2
kappa <- 2

# Create mizerParam object
params <- MizerParams(balticParams,
                      kappa_ben = kappa_ben,
                      kappa = kappa,
                      w_bb_cutoff = w_bb_cutoff,
                      w_pp_cutoff = w_pp_cutoff,
                      r_pp = r_pp,
                      r_bb = r_bb)

params@species_params$h

# Projectmodel
t_max <- 800

m1 <- project(params,
              temperature = rep(10, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

m1@params@linetype <- rep("solid", 6)

# Check dynamics and density
plotBiomass(m1) + theme_classic(base_size = 14)

plotSpectra(m1, algae = FALSE)


#**** Check growth =======================================================================
# Growth rates are extremely low
plotGrowthCurves(m1, max_age = 15) + 
  geom_point(data = subset(vbgedat, AgeRings < 16), 
             aes(AgeRings, Weight_g), size = 2, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = empiri, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(data = subset(empiri, age < 16), aes(age, weight), 
            color = col[2], size = 1.4, linetype = "twodash") +
  geom_line(aes(x = Age, y = value), 
            color = col[1], size = 1.4) +
  facet_wrap(~ Species, scales = "free", ncol = 3) +
  guides(color = FALSE, linetype = FALSE) +
  theme_classic(base_size = 12) + 
  theme(aspect.ratio = 3/4) +
  NULL


#** 2. Tune growth rate ============================================================
# We start by increasing the constant in the maximum consumption rate by a factor of 1.75
# If growth is more reasonable, tune down kappa (higher consumption allows for lower kappa to get coexistence).
params2 <- params

kappa_ben <- 1
kappa <- 1

# Increase maximum consumption rates by a factor 1.75
params2@species_params$species
params2@species_params$h[1] <- params@species_params$h[1] * 1.5 # cod
params2@species_params$h[2] <- params@species_params$h[2] * 2.2 # sprat 
params2@species_params$h[3] <- params@species_params$h[3] * 2.2 # herring 

# Remove gamma, because it needs to be recalculated using the new h. 
params2@species_params <- subset(params2@species_params,
                                 select = -gamma)

params2_upd <- MizerParams(params2@species_params,
                           kappa_ben = kappa_ben,
                           kappa = kappa,
                           w_bb_cutoff = w_bb_cutoff,
                           w_pp_cutoff = w_pp_cutoff,
                           r_pp = r_pp,
                           r_bb = r_bb)

# Project model
t_max <- 250

m2 <- project(params2_upd,
              dt = dt,
              temperature = rep(10, t_max),
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

m2@params@linetype <- rep("solid", 6)

# Check dynamics and density
plotBiomass(m2) + theme_classic(base_size = 14)


#**** Check growth =======================================================================
# Growth rates look slightly better, ok for now since they will change after calibrating r_max
plotGrowthCurves(m2, max_age = 15) + 
  geom_point(data = subset(vbgedat, AgeRings < 16), 
             aes(AgeRings, Weight_g), size = 2, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = empiri, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(aes(x = Age, y = value), 
            color = col[1], size = 1.4) +
  geom_line(data = subset(empiri, age < 16), aes(age, weight), 
            color = col[2], size = 1.4, linetype = "twodash") +
  facet_wrap(~ Species, scales = "free", ncol = 3) +
  guides(color = FALSE, linetype = FALSE) +
  theme_classic(base_size = 12) + 
  theme(aspect.ratio = 3/4) +
  NULL

# Feeding levels look ok
plotFeedingLevel(m2) + theme_classic(base_size = 14)

# After improving growth, reduce kappa from step one sequntially until coexistence is not possible and increase slightly. I landed on kappa = 1. 

# But now we need to fix so that the relative and absolute biomasses are realistic:
getSSB(m2)[t_max, ] * m2@params@species_params$sd25.29.32_m.2 / 1e9 
balticParams$AveSpawnBiomass

# The SSB's will be optimized in the next step.


#** 3. Optimize r_max ==============================================================
# Find r_max values that minimize the residual sum of squares between SSB from stock assessment and model output 
# Below I'm sourcing that functions below are written by by Asta Audzijonyte, Jon Reum and Julia Blanchard, that I've modified to better fit my model system (Asta, unpublished; Reum et al 2018 Oikos, Blanchard et al 2014, J. Applied. Ecol)

# Make sure MizerParams() in FunctionsForOptim.R uses the same parameter-setup that you've decided so far!


#**** Source functions =============================================================
# See the functions script for description of their arguments and such
script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/StartVector.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/FunctionsForOptim.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# With only r_max to be optimized, and three species, this function takes 7 minutes on my laptop
# Note also that I have NOT turned off warnings from creating mizerParams objects, so they will printed in the console.
# If this code takes to long to run, I will do it in a different script and store the .Rdata

system.time(optim_resultNC <- optim(
  par = start_vector,         # Vector of starting parameter values
  startpars = startpars,
  fn = calibratePar_Baltic,   # Function to feed parameter values to the mizer object, 
                              # Run the model and return the error. 
  modelParams = balticParams, # Pass the original parameter file 
                              # (Jon passed his sim object after the burn-in period, 
                              # But I don't run burn-in for now)
  lower = lower_bounds(
    startpars = startpars),   # Function to set lower bounds, 
                              # We pass a user constructed list 'startpars' that contains data 
  upper = upper_bounds(
    startpars = startpars),   # Same for upper bounds
  method = "L-BFGS-B",
  control = list(maxit = 10, 
                 REPORT = 1, 
                 trace = 6), 
  effort = effort)) 

# optimized r_max in case I don't want to run optim again:
# with allometric erepro (and h-scaling factor = 1.75):
# > exp(optim_resultNC$`par`)
# [1] 0.003937455 9.493599952 1.046058740

# with allometric erepro (and h-scaling factor = c(1.5, 2.2, 2.2)):
# > exp(optim_resultNC$`par`)
# [1] 0.01076359 5.02863794 0.51540528

# Compare with North Sea r_max (default)
balticParams$r_max
exp(optim_resultNC$`par`)

# Update Rmax
params3_upd <- params2_upd

# This is the final set of parameters I use after calibration has been done
params3_upd@species_params$r_max <- exp(optim_resultNC$`par`)

params3_upd@species_params

# Project model with optimized r_max
t_max <- 60

m3 <- project(params3_upd,
              dt = dt,
              temperature = rep(10, t_max),
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

m3@params@linetype <- rep("solid", 6)

# Check dynamics and density
plotBiomass(m3) + theme_classic(base_size = 14)

#plot(m3)

plotSpectra(m3, algae = FALSE) + theme_classic(base_size = 14)


#**** Check growth =======================================================================
# Check growth rate are reasonable
plotGrowthCurves(m3, max_age = 15) + 
  geom_point(data = subset(vbgedat, AgeRings < 16), 
             aes(AgeRings, Weight_g), size = 2, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = empiri, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(data = subset(empiri, age < 16), aes(age, weight), 
            color = col[2], size = 1.4, linetype = "twodash") +
  geom_line(aes(x = Age, y = value), 
            color = col[1], size = 1.4) +
  facet_wrap(~ Species, scales = "free", ncol = 3) +
  guides(color = FALSE, linetype = FALSE) +
  theme_classic(base_size = 12) + 
  theme(aspect.ratio = 3/4) +
  NULL

# They still look OK after calibration


#**** SSB pred/fit =================================================================
# Plot predicted vs fitted SSB. Take mean of last 20 time steps
obs <- balticParams$AveSpawnBiomass
pred <- colMeans(getSSB(m3)[c(I(dim(m3@n)[1]-20):dim(m3@n)[1]), ] ) * (balticParams$sd25.29.32_m.2) / (1e9)

ssb_eval <- data.frame(SSB = c(obs, pred),
                       Source = rep(c("Observed", "Predicted"), each = 3),
                       Species = rep(balticParams$species, 2))

p1 <- ggplot(ssb_eval, aes(Species, SSB, shape = Source, fill = Species)) + 
  geom_point(size = 6, alpha = 0.5) +
  scale_fill_manual(values = col) +
  scale_shape_manual(values = c(21, 24),
                     guide = guide_legend(override.aes = list(colour = "black", 
                                                              fill = "black",
                                                              size = 4))) + 
  labs(x = "", y = "Spawning stock biomass (millions kg)") +
  theme_classic(base_size = 14) +
  guides(fill = FALSE) +
  theme(legend.position = c(0.15, .85),
        legend.title = element_blank(),
        aspect.ratio = 1)

ssb_eval_l <- data.frame(obs = log10(obs), pred = log10(pred), Species = balticParams$species)

p2 <- ggplot(ssb_eval_l, aes(obs, pred, fill = Species, shape = Species)) +
  geom_point(size = 6, alpha = 0.5) +
  scale_fill_manual(values = col) +
  scale_shape_manual(values = c(21, 22, 24)) +
  labs(x = "Log10(Observed SSB)", y = "Log10(Predicted SSB)") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.15, .85),
        legend.title = element_blank(),
        aspect.ratio = 1)

p1 + p2  


#**** Recruitment & density dependence =============================================
# Evalute how much "density dependence" there is in the model from the stock-recruit relationship
# Calculate the density independent recruitment (total egg production) R_{p.i} before density dependence
rdi <- getRDI(m3@params,
              m3@n[t_max,,],
              m3@n_pp[t_max,],
              m3@n_bb[t_max,],
              m3@n_aa[t_max,],
              m3@intTempScalar[,,(t_max/dt)],
              m3@metTempScalar[,,(t_max/dt)])

# Calculate the flux entering the smallest size class of each species (recruitment - density dependence)
rdd <- getRDD(m3@params,
              m3@n[t_max,,],
              m3@n_pp[t_max,],
              m3@n_bb[t_max,],
              m3@n_aa[t_max,],
              sex_ratio = 0.5,
              m3@intTempScalar[,,(t_max/dt)],
              m3@metTempScalar[,,(t_max/dt)])

# rdi vs rdd
rdi / rdd
rdd / rdi

rec <- data.frame(Species = names(rdi),
                  "RDI/RDD" = rdi/rdd)

ggplot(rec, aes(Species, RDI.RDD, color = Species, shape = Species)) +
  scale_color_manual(values = col) +
  labs(x = "", y = "RDI / RDD") +
  theme_classic(base_size = 14) +
  geom_abline(intercept = 1, slope = 0, color = "gray30", 
              linetype = 2, size = 1) +
  geom_point(size = 5) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        aspect.ratio = 1)

# Get RDD to rmax ratio
rdd / params3_upd@species_params$r_max

# Get RDI to rmax ratio
rdi / params3_upd@species_params$r_max

# **** which is more important? rdd to rdi or rmax? rmax right? what's a good value for the ratio?

#**** Diet =========================================================================
# This is Jon's function 
# Looks OK, but clear that we need size-varying theta
plotDietComp(m3, prey = dimnames(m3@diet_comp)$prey[1:5]) + 
  theme_classic(base_size = 14) +
  scale_fill_manual(values = rev(col),
                    labels = c("Cod", "Sprat", "Herring", "Plankton", "Benthos")) +
  scale_x_continuous(name = "log10 predator mass (g)", expand = c(0,0)) +
  scale_y_continuous(name = "Proportion of diet by mass (g)", expand = c(0,0)) +
  NULL


# D. VALIDATE WITH TIME SERIES =====================================================
#**** Set up time varying effort ===================================================
# In this section I will try out and summarize different ways to get time series validation (burn in etc)

# This procedure is based on the Mizer vignette
# Read in historical F
f_history <- as.matrix(
  ssb_f %>% 
    filter(Year > 1973 & Year < 2013 & !Species == "Flounder") %>% 
    select(Year, Species, Fm) %>% 
    spread(Species, Fm) %>%
    select(Cod, Herring, Sprat)
)

rownames(f_history) <- 1974:2012  # If I wan the full time series, do: 1974:2012. Now I just to 10 years prior and after

# Before projecting forward, we want to remove the transient dynamics.
# Judging from the last projection, 40 yrs burn in seems ok.
# We also start the time series at 1974 (1992-2002 for calibration)
# That means we start the simulation to get a time series free from transients 
# by starting from 1974-burnin
plotBiomass(m3) + theme_classic(base_size = 14)

burnin <- 40

# Here I use the F at the first time step of the series, as in 'mizer'
initial_effort <- matrix(f_history[1, ], # colMeans(f_history[1:19, ])
                         byrow = TRUE,
                         nrow = burnin,
                         ncol = ncol(f_history), 
                         dimnames = list(1934:1973)) 

all_effort <- rbind(initial_effort, f_history)

dimnames(all_effort)

# Now project forward with the effort matrix
# When effort is provided as an array, the t_max argument is not used, but is taken from dimensions of effort.
# Temperature must match the time steps in length (temperature = rep(params@t_ref, times = t_max))
# from "project" code
time_effort <- as.numeric(dimnames(all_effort)[[1]])
t_max <- time_effort[length(time_effort)]

m4 <- project(params3_upd, 
              dt = 0.1,
              effort = all_effort,
              temperature = rep(10, 2012),
              diet_steps = 10,
              t_max = t_max) 

m4@params@linetype <- rep("solid", 6)

# Full time series including burn in
plotBiomass(m4) + theme_classic(base_size = 14)

# From 1974 and forward
plotBiomass(m4) + 
  theme_classic(base_size = 14) +
  xlim(1974, 2014) +
  ylim(1, 12)


#**** Compare with assessment data =================================================
pred_ssb <- data.frame(getSSB(m4))
pred_ssb$Year <- rownames(getSSB(m4))
pred_ssb$source <- "Predicted"

# Convert to long data frame (1 obs = 1 row)
pred_ssb_l <- pred_ssb %>% 
  filter(Year > 1973) %>% 
  gather(Species, ssb_g.m2, 1:3)

# Scale up from g/m2 to 10^6 kg / Baltic
pred_ssb_l$SSB <- pred_ssb_l$ssb_g.m2 * (balticParams$sd25.29.32_m.2) / (1e9)

pred_ssb_l <- pred_ssb_l %>% select(-ssb_g.m2)

str(pred_ssb_l)

# Now do observed ssb
obs_ssb_l <- ssb_f %>%
  filter(Year > 1973) %>% 
  select(Species, Year, SSB)

obs_ssb_l$source <- "Observed"

str(obs_ssb_l)

# And finally trawl survey index
# Need to add it on secondary axis in that case, maybe I won't add it afterall
# obs_tsb_l <- ssb_f %>%
#   filter(Year > 1973) %>% 
#   select(Species, Year, TSB)
# obs_tsb_l$source <- "Survey CPUE"
# colnames(obs_tsb_l)[3] <- "SSB" # Rename column for plotting purposes
# str(obs_tsb_l)

# Combine observed and predicted SSB
dat <- data.frame(rbind(obs_ssb_l, pred_ssb_l))
dat$Year <- as.integer(dat$Year)

# Normalize to maximum within species
dat <- dat %>% 
  drop_na() %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  mutate(test = max(SSB)) %>% 
  mutate(SSB_norm = SSB / max(SSB))

df <- data.frame(g = rep(c("a", "b"), each = 5),
                 y = runif(n=10))

df %>% 
  group_by(g) %>% 
  mutate(y_norm = y/max(y))

# Get mean data for calibration period
# For some reason this doesn't work with ggplot???
# mdat <- data.frame(filter(mean_ssb_F, !Species == "Flounder"))

# Plot predicted and observed ssb by species, normalize by max within species
ggplot(dat, aes(Year, SSB_norm, linetype = source, color = source)) +
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = 1),
            fill  = "gray50", alpha = 0.05) +
  geom_line(size = 1.7, alpha = 1) +
  # geom_point(data = mdat, aes(Year, mean_SSB), size = 1.2) + This is not really working atm
  scale_linetype_manual(values = c("solid", "twodash")) +
  scale_color_manual(values = col) +
  theme(aspect.ratio = 1) +
  ylab("SSB (1000 tonnes)") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  #coord_cartesian(xlim = c(1980, 2010)) + # This is the time frame used in Jacobsen et al 
  NULL


#**** Estimate FMSY from model - compare with assessment ===========================
# In lack of a better approach, I will just for-loop different F, extract Yield, plot
# over F. I will increase each species F separately, keeping the others at their mean

F_range <- seq(0.2, 1.4, 0.01)
t_max <- 75

#****** Cod ========================================================================
# Mean F in calibration time
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])

# FMSY from assessment
# effort = c(Cod = 0.55, 
#            Herring = 0.3, 
#            Sprat = 0.285)

# First do a test run
t <- project(params3_upd, 
             dt = 0.1,
             effort = effort, 
             temperature = rep(10, t_max),
             diet_steps = 10,
             t_max = t_max) 

plotSpectra(t)

length(t@params@linetype)

# Create empty data holder
codFmsy <- c()
Y <- c()
Fm <- c()
t <- c()

for(i in F_range) {
  
  effort[1] <- i
  
  t <- project(params3_upd, 
               dt = 0.2, # Use 0.2 here otherwise really slow
               effort = effort, 
               temperature = rep(10, t_max),
               diet_steps = 10,
               t_max = t_max) 
  
  Y <- mean(data.frame(getYield(t))$Cod[(t_max-20):t_max])
  Fm <- i
  t <- cbind(Y, Fm)
  
  codFmsy <- data.frame(rbind(t, codFmsy))
  
}

codFmsy$Species <- "Cod"

ggplot(codFmsy, aes(Fm, Y)) + geom_line()


#****** Sprat ======================================================================
# Mean F in calibration time
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])

# FMSY from assessment
# effort = c(Cod = 0.55, 
#            Herring = 0.3, 
#            Sprat = 0.285)

# Create empty data holder
sprFmsy <- c()
Y <- c()
Fm <- c()
t <- c()

for(i in F_range) {
  
  effort[3] <- i
  
  t <- project(params3_upd, 
               dt = 0.2,
               effort = effort, 
               temperature = rep(10, t_max),
               diet_steps = 10,
               t_max = t_max) 
  
  Y <- mean(data.frame(getYield(t))$Sprat[(t_max-20):t_max])
  Fm <- i
  t <- cbind(Y, Fm)
  
  sprFmsy <- data.frame(rbind(t, sprFmsy))
  
}

sprFmsy$Species <- "Sprat"

ggplot(sprFmsy, aes(Fm, Y)) + geom_line()


#****** Herring ====================================================================
# Mean F in calibration time
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])

# FMSY from assessment
# effort = c(Cod = 0.55, 
#            Herring = 0.3, 
#            Sprat = 0.285)

# Create empty data holder
herFmsy <- c()
Y <- c()
Fm <- c()
t <- c()

for(i in F_range) {
  
  effort[2] <- i
  
  t <- project(params3_upd, 
               dt = 0.2,
               effort = effort, 
               temperature = rep(10, t_max),
               diet_steps = 10,
               t_max = t_max) 
  
  Y <- mean(data.frame(getYield(t))$Herring[(t_max-20):t_max])
  Fm <- i
  t <- cbind(Y, Fm)
  
  herFmsy <- data.frame(rbind(t, herFmsy))
  
}

herFmsy$Species <- "Herring"

ggplot(herFmsy, aes(Fm, Y)) + geom_line()


#****** All together ===============================================================
Fmsy <- rbind(codFmsy, sprFmsy, herFmsy)

ggplot(Fmsy, aes(Fm, Y, color = Species)) + 
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  theme_classic(base_size = 14)

# Now plot FMSY from model and assessment
# Since FMSY is not available for the entire time period, I will take FMSY values
# from the assessment report from which I took historical data from (back-calculated)

# Cod (from Advice 2013)
codassesFMSY <- data.frame(Source = c("Single Species Assessment",
                                      "Multispecies Assessment",
                                      "Size Spectrum Model"),
                           FMSY = c(0.46, 0.55, codFmsy$Fm[codFmsy$Y == max(codFmsy$Y)]),
                           Species = rep("Cod", 3))

# Sprat (from Advice 2014)
sprassesFMSY <- data.frame(Source = c("Single Species Assessment",
                                      "Multispecies Assessment lwr",
                                      "Multispecies Assessment upr",
                                      "Size Spectrum Model"),
                           FMSY = c(0.29, 0.25, 0.32, sprFmsy$Fm[sprFmsy$Y == max(sprFmsy$Y)]),
                           Species = rep("Sprat", 4))


# Herring (from Advice 2014)
herassesFMSY <- data.frame(Source = c("Single Species Assessment",
                                      "Multispecies Assessment lwr",
                                      "Multispecies Assessment upr",
                                      "Size Spectrum Model"),
                           FMSY = c(0.26, 0.25, 0.35, herFmsy$Fm[herFmsy$Y == max(herFmsy$Y)]),
                           Species = rep("Herring", 4))

asses_mod_FMSY <- rbind(codassesFMSY, sprassesFMSY, herassesFMSY)


ggplot(asses_mod_FMSY, aes(Species, FMSY, shape = Source, fill = Source)) + 
  geom_point(size = 6, alpha = 0.7, color = "white") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = seq(21, 25)) +
  theme_classic(base_size = 16)



#**** Extra plots ==================================================================
# Plot predicted vs observed, all years
dat %>% 
  select(-SSB_norm) %>% 
  spread(source, SSB) %>% 
  mutate(calibration = ifelse(Year > (min_cal_yr-1) & Year < (max_cal_yr+1),
                              "within calibratio",
                              "outside calibration")) %>% 
  ggplot(., aes(Predicted, Observed, color = Year, shape = calibration)) +
  geom_abline(intercept = 0, slope = 1, 
              colour = "red", linetype = 1, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_viridis(option = "cividis")   +
  facet_wrap(~Species, scales = "free") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL

# Plot predicted vs observed, forecast only
dat %>%
  filter(Year > 1991) %>% 
  select(-SSB_norm) %>% 
  spread(source, SSB) %>% 
  mutate(calibration = ifelse(Year > (min_cal_yr-1) & Year < (max_cal_yr+1),
                              "within calibratio",
                              "outside calibration")) %>% 
  ggplot(., aes(Predicted, Observed, color = Year, shape = calibration)) +
  geom_abline(intercept = 0, slope = 1, 
              colour = "red", linetype = 1, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_viridis(option = "cividis") +
  facet_wrap(~Species, scales = "free") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL

# Plot SSB vs F
head(ssb_f)

ssb_f_2 <- ssb_f %>% 
  filter(Year > 1973 & Year < 2013) %>% 
  mutate(Fm_lag_1 = lag(Fm, 1),
         Fm_lag_2 = lag(Fm, 2),
         SSB_pred = filter(dat, source == "Predicted")$SSB) %>% 
  select(Species, Year, SSB, SSB_pred, Fm, Fm_lag_1, Fm_lag_2) %>% 
  gather(source, SSB, 3:4)

# No lag
ggplot(ssb_f_2, aes(Fm, SSB, color = source)) + 
  scale_color_manual(values = col[c(1,3)]) +
  geom_point(size = 4, alpha = 0.8) +
  stat_smooth(method = "gam") +
  facet_wrap(~Species, scales = "free") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL

# Lag 1
ggplot(ssb_f_2, aes(Fm_lag_1, SSB, color = source)) + 
  scale_color_manual(values = col[c(1,3)]) +
  geom_point(size = 4, alpha = 0.8) +
  stat_smooth(method = "gam") +
  facet_wrap(~Species, scales = "free") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL

# Lag 2
ggplot(ssb_f_2, aes(Fm_lag_2, SSB, color = source)) + 
  scale_color_manual(values = col[c(1,3)]) +
  geom_point(size = 4, alpha = 0.8) +
  stat_smooth(method = "gam") +
  facet_wrap(~Species, scales = "free") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL


#**** Brief summary ================================================================
# Overall I suppose this is a fairly OK calibration (process)
# growth rates became OK after increasing Cmax by a factor of 1.75

# These are the default h-values:
# > params@species_params$h
# [1] 20.733625  6.598427  6.875000

# These are the values increased by a factor of 1.75
# > params2@species_params$h
# [1] 31.10044 14.51654 15.12500

# These are the values with Gustav's new equation
# [1]  31.48966 148.51985  37.07342

# So the 1.7 are between the two equations, which seems OK

# ssb fits really well to the average F and SSB in the calibration period. This was done using r_max.
# recruitment looks ok to. There is some density dependence in the model (90% of r_max, r_phys = 10 * rec)
# rec:r_max could be lower, but the important thing is they respond to F
# They seem to do that actually more than in the data, so that's likely ok
# Yet, the predicted dynamics are changing less compared to observed, and in general it doesn't capture regime shifts (amplitude is off). However this is not strange given that the model essentially only uses fishing as input whereas environment likely casued the gadoid outburst seen in the data. 


#*** Below code should be move to an "analysis" script












# E. PROJECT UNDER CLIMATE AND FISHING SCENARIOS ===================================
#**** Set up time varying effort & constant temperature ============================
# In the temperature-invariant simulation I created a list holding historical
# effort values for all species up until 2012. Since I now want to project forward
# with temperature, I need to extend that to include FMSY-effort for all species
# up until 2050.
# Moreover, as t_max is taken from the effort-array and not the argument in 
# project, and temperature needs to have the same length as t_max, I here center
# effort to start at 1 and the actual year. I go back to normal scale when plotting
# First I need to test if that's correct:

# All effort until 2012
all_effort

# FMSY from assessment
asses_mod_FMSY

# Create effort for projection with temperature by appending FMSY to effort data
t_future <- 2050-2012

projectEffort_fwr_df <- data.frame(all_effort[1:t_future,])

projectEffort_fwr_df[, 1] <- 0.55  # Cod
projectEffort_fwr_df[, 2] <- 0.3   # Herring
projectEffort_fwr_df[, 3] <- 0.285 # Sprat

projectEffort_fwr <- as.matrix(projectEffort_fwr_df)
rownames(projectEffort_fwr) <- 2013:2050

projectEffort <- rbind(all_effort, projectEffort_fwr)

# Create vector of centered rownames
projectEffort_ct <- projectEffort
rownames(projectEffort_ct) <- 1:nrow(projectEffort)
projectEffort_ct

# Plot effort data
plotEffort <- data.frame(projectEffort)

plotEffort$Year <- as.numeric(as.character(rownames(projectEffort)))

plotEffort %>% 
  gather(Species, Effort, 1:3) %>% 
  ggplot(., aes(Year, Effort, color = Species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = 1.4),
            fill  = col[2], alpha = 0.1) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = 1974, 
                xmax = 2012,
                ymin = 0,
                ymax = 1.4),
            fill  = col[2], alpha = 0.05) +
  geom_line(size = 1.3) +
  theme_classic(base_size = 15) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = col[c(1,3,5)]) +
  theme(aspect.ratio = 3/4) +
  labs(x = "Year", y = "Fishing mortality (F)") +
  NULL

# First testing if I can reproduce constant-temperature scenarios with centered effort
m5 <- project(params3_upd, 
              dt = 0.1,
              effort = projectEffort_ct,
              temperature = rep(10, nrow(projectEffort_ct)),
              diet_steps = 10,
              t_max = t_max,
              t_ref = 10) 

m5@params@linetype <- rep("solid", 6)
m5@params@linecolour[1:3] <- col[c(1,5,3)]

# Full time series including burn in
plotBiomass(m5) + theme_classic(base_size = 14)

# From 1974 and forward (or after burn-in)
plotBiomass(m5) + 
  theme_classic(base_size = 14) +
  #xlim(1974, 2050) +
  xlim(41, 120) +
  #ylim(1, 12) +
  NULL


#**** Compare with assessment data =================================================
pred_ssb <- data.frame(getSSB(m5))
pred_ssb$Year_ct <- as.numeric(rownames(getSSB(m5)))
pred_ssb$Year <- pred_ssb$Year_ct + (1934-1) # 1934 is so that 40 year burn-in leads to start at 1974
pred_ssb$source <- "Predicted"
str(pred_ssb)

# Convert to long data frame (1 obs = 1 row)
# The first year with real effort is 1974. This is year 41 with centered time (1974-1934 +1),
# The last year before FMSY is 2012. In centered time this is 2012-1934 = 79
pred_ssb_l <- pred_ssb %>% 
  filter(Year > 1973 & Year < 2013) %>% # 1973 with non-centered year
  gather(Species, ssb_g.m2, 1:3)

# Note the projection goes to 2050 in raw time, but here I want to reproduce a specific figure

# Scale up from g/m2 to 10^6 kg / Baltic
pred_ssb_l$SSB <- pred_ssb_l$ssb_g.m2 * (balticParams$sd25.29.32_m.2) / (1e9)

pred_ssb_l <- pred_ssb_l %>% select(-ssb_g.m2)

str(pred_ssb_l)

# Now do observed ssb
str(ssb_f)
obs_ssb_l <- ssb_f %>%
  filter(Year > 1973) %>% 
  select(Species, Year, SSB)

obs_ssb_l$source <- "Observed"

str(obs_ssb_l)

# Centering year so that 1974 is Year_ct = 41
obs_ssb_l$Year_ct <- (obs_ssb_l$Year-(1974))+(41) # 1974 is first year, t_1 needs to be one. Then +40 to match predicted

# Note! SSB goes to 2014, wheres pred goes to 2012 in the current subset

# Combine observed and predicted SSB
dat <- data.frame(rbind(obs_ssb_l, pred_ssb_l))
dat$Year <- as.integer(dat$Year)

str(dat)
dat

# Normalize to maximum within species
dat <- dat %>% 
  drop_na() %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  mutate(test = max(SSB)) %>% 
  mutate(SSB_norm = SSB / max(SSB))

# Plot predicted and observed ssb by species, normalize by max within species
ggplot(dat, aes(Year, SSB_norm, linetype = source, color = source)) +
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = 1),
            fill  = col[2], alpha = 0.05) +
  geom_line(size = 1.7, alpha = 0.7) +
  # geom_point(data = mdat, aes(Year, mean_SSB), size = 1.2) + This is not really working atm
  scale_linetype_manual(values = c("solid", "twodash")) +
  scale_color_manual(values = col[c(1,3)]) +
  theme(aspect.ratio = 1) +
  labs(y = "SSB (1000 tonnes)", x = "Year") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  #coord_cartesian(xlim = c(1980, 2010)) + # This is the time frame used in Jacobsen et al
  NULL
# Looks like it produces the same plots as in the non-scaled scenario.


#**** Calculate and plot correlation coefficients ==================================
obs_df <- filter(dat, source == "Observed" & Year < 2013)
pred_df <- filter(dat, source == "Predicted")

cor_df <- data.frame(Year = obs_df$Year,
                     Obs = obs_df$SSB,
                     Pred = pred_df$SSB,
                     Species = obs_df$Species)

# Test
ggplot(cor_df, aes(Year, Obs, color = Species)) + geom_line()
ggplot(cor_df, aes(Year, Pred, color = Species)) + geom_line()

cor.test(subset(cor_df, Species == "Cod")$Pred, subset(cor_df, Species == "Cod")$Obs)
cor.test(subset(cor_df, Species == "Herring")$Pred, subset(cor_df, Species == "Herring")$Obs)
cor.test(subset(cor_df, Species == "Sprat")$Pred, subset(cor_df, Species == "Sprat")$Obs)

# All in one:
cors <- ddply(cor_df, c("Species"), summarise, cor = round(cor(Pred, Obs), 2))

# Plot correlation between predicted and observed
ggplot(cor_df, aes(Obs, Pred)) +
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
  geom_point(size = 3, alpha = 0.7, color = col[2]) +
  theme(aspect.ratio = 1) +
  labs(y = "Predicted", x = "Observed") +
  theme_classic(base_size = 14) +
  scale_y_continuous(expand = c(0, 0)) + 
  geom_text(data = cors, aes(label = paste("r = ", cor, sep = "")), 
            x = c(500, 1370, 1550), y = c(120, 595, 1020),
            fontface = "italic", size = 5) +
  theme(aspect.ratio = 3/4) +
  NULL


#**** Set up time varying temperature ==============================================
# Read data
temp_datRCP8.5 <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-baltic-sea/master/data/Climate/Test_RCP8.5_from_graph.csv"), sep = ";", stringsAsFactors = FALSE)

head(temp_datRCP8.5)
str(temp_datRCP8.5)

# Calculate average temperature by year
tempDat <- data.frame(temp_datRCP8.5 %>% 
                        dplyr::filter(year < 2050.5) %>% 
                        dplyr::mutate(Year.r = factor(round(year, digits = 0))) %>% 
                        dplyr::group_by(Year.r) %>% 
                        dplyr::summarize(mean_temp = mean(rel_temp_70.99Ave)) %>% 
                        dplyr::mutate(Year.num = as.numeric(as.character(Year.r))))

str(tempDat)
head(tempDat)
tail(tempDat)

# Plot temperature data
ggplot(tempDat, aes(Year.num, mean_temp)) +
  geom_hline(yintercept = 0, size = 1, alpha = 0.8, linetype = 2) +
  geom_line(size = 1.3, col = col[2]) +
  theme_classic(base_size = 15) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 3/4) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = -0.9,
                ymax = 1.8),
            fill  = col[2], alpha = 0.05) +
  labs(x = "Year", y = "Change in Baltic Sea SST (RCP8.5)") +
  NULL

# For calibration I use mean values, including temperature
# Which temperature do I use? Should be T_0, so that all other rates are constant
# T0 is 10 here, and T is T_dev + 10

tempDat$mean_temp_scaled <- tempDat$mean_temp + 10

# Now I need to match year for temperature and effort
projectEffort # Non-centered complete effort array
projectEffort_ct # Centered complete effort array

# I need a temperature-vector with same length as t_max.
# Repeat initial temp until from 1934-1974 (when effort starts)
projectTemp_df <- tempDat %>% 
  select(Year.r, mean_temp_scaled) %>% 
  rename(Year = Year.r,
         temperature = mean_temp_scaled)

projectTemp_burnin <- data.frame(Year = seq(from = 1934, to = 1969, by = 1),
                                 temperature = rep(projectTemp_df$temperature[1], length(seq(1934:1969))))

projectTemp <- rbind(projectTemp_burnin, projectTemp_df)
nrow(projectTemp)

# Set calibration-time to 10C
# ****** NO, I DON'T DO THIS WITH EFFORT. INSTEAD I SHOULD CALIBRATE TO MEAN TEMP, 
# AND SET MEAN TEMP = T_REF
# projectTemp$temperature <- ifelse(projectTemp$Year > 1991 & projectTemp$Year < 2003,
#                                   10,
#                                   projectTemp$temperature)

#**** Combine data and prepare for projection ======================================
params3_upd@species_params$ea_int <- 0.53
params3_upd@species_params$ea_met <- 0.57
params3_upd@species_params$ea_mor <- 0.57

params3_upd@species_params$ca_int <- -0.0025
params3_upd@species_params$ca_met <- -0.0019

# Project with temperatue and effort varying through time
m6 <- project(params3_upd, 
              dt = 0.1,
              effort = projectEffort_ct,
              temperature = projectTemp$temperature,
              diet_steps = 10,
              t_max = t_max,
              t_ref = 10) 

m6@params@linetype <- rep("solid", 6)
m6@params@linecolour[1:3] <- col[c(1,5,3)]

# Full time series including burn in
plotBiomass(m6) + theme_classic(base_size = 14)

# From 1974 and forward (or after burn-in)
plotBiomass(m6) + 
  theme_classic(base_size = 14) +
  xlim(41, 120) +
  NULL


#**** Compare with assessment data =================================================
pred_ssb <- data.frame(getSSB(m6))
pred_ssb$Year_ct <- as.numeric(rownames(getSSB(m6)))
pred_ssb$Year <- pred_ssb$Year_ct + (1934-1) # 1934 is so that 40 year burn-in leads to start at 1974
pred_ssb$source <- "Predicted"
str(pred_ssb)

# Convert to long data frame (1 obs = 1 row)
# The first year with real effort is 1974. This is year 41 with centered time (1974-1934 +1),
# The last year before FMSY is 2012. In centered time this is 2012-1934 = 79
pred_ssb_l <- pred_ssb %>% 
  filter(Year > 1973) %>%
  gather(Species, ssb_g.m2, 1:3)

# Scale up from g/m2 to 10^6 kg / Baltic
pred_ssb_l$SSB <- pred_ssb_l$ssb_g.m2 * (balticParams$sd25.29.32_m.2) / (1e9)

pred_ssb_l <- pred_ssb_l %>% select(-ssb_g.m2)

str(pred_ssb_l)

# Now do observed ssb
str(ssb_f)
obs_ssb_l <- ssb_f %>%
  filter(Year > 1973) %>% 
  select(Species, Year, SSB)

obs_ssb_l$source <- "Observed"

str(obs_ssb_l)

# Centering year so that 1974 is Year_ct = 41
obs_ssb_l$Year_ct <- (obs_ssb_l$Year-(1974))+(41) # 1974 is first year, t_1 needs to be one. Then +40 to match predicted

# Combine observed and predicted SSB
dat <- data.frame(rbind(obs_ssb_l, pred_ssb_l))
dat$Year <- as.integer(dat$Year)

str(dat)
dat

# Normalize to maximum within species
dat <- dat %>% 
  drop_na() %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  mutate(test = max(SSB)) %>% 
  mutate(SSB_norm = SSB / max(SSB))

# Plot predicted and observed ssb by species, normalize by max within species
dat %>% filter(Year < 2012) %>% 
  ggplot(., aes(Year, SSB_norm, linetype = source, color = source)) +
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = 1),
            fill  = col[2], alpha = 0.05) +
  geom_line(size = 1.7, alpha = 0.7) +
  # geom_point(data = mdat, aes(Year, mean_SSB), size = 1.2) + This is not really working atm
  scale_linetype_manual(values = c("solid", "twodash")) +
  scale_color_manual(values = col[c(1,3)]) +
  theme(aspect.ratio = 1) +
  labs(y = "SSB (1000 tonnes)", x = "Year") +
  theme_classic(base_size = 14) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 3/4) +
  #coord_cartesian(xlim = c(1980, 2010)) + # This is the time frame used in Jacobsen et al
  NULL
# Looks like it produces the same plots as in the non-scaled scenario.


#**** Calculate and plot correlation coefficients ==================================
obs_df <- filter(dat, source == "Observed" & Year < 2013)
pred_df <- filter(dat, source == "Predicted" & Year < 2013)

cor_df <- data.frame(Year = obs_df$Year,
                     Obs = obs_df$SSB,
                     Pred = pred_df$SSB,
                     Species = obs_df$Species)

# Test
ggplot(cor_df, aes(Year, Obs, color = Species)) + geom_line()
ggplot(cor_df, aes(Year, Pred, color = Species)) + geom_line()

cor.test(subset(cor_df, Species == "Cod")$Pred, subset(cor_df, Species == "Cod")$Obs)
cor.test(subset(cor_df, Species == "Herring")$Pred, subset(cor_df, Species == "Herring")$Obs)
cor.test(subset(cor_df, Species == "Sprat")$Pred, subset(cor_df, Species == "Sprat")$Obs)

# All in one:
cors <- ddply(cor_df, c("Species"), summarise, cor = round(cor(Pred, Obs), 2))

# Plot correlation between predicted and observed
ggplot(cor_df, aes(Obs, Pred)) +
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
  geom_point(size = 3, alpha = 0.7, color = col[2]) +
  theme(aspect.ratio = 1) +
  labs(y = "Predicted", x = "Observed") +
  theme_classic(base_size = 14) +
  scale_y_continuous(expand = c(0, 0)) + 
  geom_text(data = cors, aes(label = paste("r = ", cor, sep = "")), 
            x = c(500, 1370, 1550), y = c(120, 595, 1020),
            fontface = "italic", size = 5) +
  theme(aspect.ratio = 3/4) +
  NULL


#**** Plot growth rates ============================================================
# Here I need to extract the data that goes into the the plotGrowthCurves function.
# Copy and edit the function to only extract data. Do that for multiple models and plot

# 1. Changes in growth: calibration period relative to warming. 

# Define new function
getGrowth <- function(object, species,
                      max_age = 20, percentage = FALSE, print_it = TRUE) {
  if (is(object, "MizerSim")) {
    sim <- object
    if (missing(species)) {
      species <- dimnames(sim@n)$sp
    }
    # reorder list of species to coincide with order in sim
    idx <- which(dimnames(sim@n)$sp %in% species)
    species <- dimnames(sim@n)$sp[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list("Species" = species, "Age" = age))
    g <- getEGrowth(sim@params, sim@n[dim(sim@n)[1], , ], 
                    sim@n_pp[dim(sim@n)[1], ], sim@n_bb[dim(sim@n)[1], ], sim@n_aa[dim(sim@n)[1], ], sim@intTempScalar[,,1], sim@metTempScalar[,,1]) #AA
    for (j in 1:length(species)) {
      i <- idx[j]
      g_fn <- stats::approxfun(sim@params@w, g[i, ])
      myodefun <- function(t, state, parameters){
        return(list(g_fn(state)))
      }
      ws[j, ] <- deSolve::ode(y = sim@params@species_params$w_min[i],
                              times = age, func = myodefun)[, 2]
      if (percentage) {
        ws[j, ] <- ws[j, ] / sim@params@species_params$w_inf[i] * 100
      }
    }
    plot_dat <- reshape2::melt(ws)
    return(plot_dat) # Added this to extract growth data in data.frame
  }    
}

# Test if correct:
# growth_warm <- getGrowth(m6)
# plotGrowthCurves(m6, species = "Cod") + 
#   geom_line(data = subset(growth_warm, Species == "Cod"), aes(Age, value), color = "red", linetype = "dashed")

# Get growth data for all scenarios to be compared
# With projected temperature (until 2050)
growth_warm <- getGrowth(m6)
growth_warm$Scenario <- "Warm"

# With constant temperature (until 2050)
growth_con <- getGrowth(m5)
growth_con$Scenario <- "Constant"

# Only calibration scenario
# growth_cal <- getGrowth(m4)
# growth_cal$Scenario <- "Calibration"

# Combine
growth_df <- rbind(growth_con, growth_warm)

# Plot together
growth_df %>% filter(Age < 15) %>% 
  ggplot(., aes(Age, value, color = Scenario, linetype = Scenario)) + 
  geom_line(size = 1.5, alpha = 0.8) +
  labs(y = "Size (g)") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = c(brewer.pal(n = 8, name = "RdBu")[8], brewer.pal(n = 8, name = "RdBu")[1])) +
  scale_linetype_manual(values = c("solid", "twodash")) + 
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  NULL

# Test plotting % instead:
growth_warm <- getGrowth(m6)
growth_warm$Scenario <- "Warm"

# With constant temperature (until 2050)
growth_con <- getGrowth(m5)
growth_con$Scenario <- "Constant"

# Now with percentage instead
growth_df_wide <- data.frame(Age = growth_warm$Age,
                             value = growth_warm$value / growth_con$value,
                             Species = growth_warm$Species)

growth_df_wide %>% filter(Age < 15 & Age > 0) %>% 
  ggplot(., aes(Age, value)) + 
  geom_hline(yintercept = 1, size = 0.5, alpha = 0.8, color = "red") +
  geom_line(size = 1.5, alpha = 0.8, color = col[2]) +
  labs(y = "Relative Size (Warm/Constant)") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(limits = c(0.85, 1.05), expand = c(0, 0)) + 
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  NULL


#**** Plot size spectra ============================================================
# Here I need to extract the data that goes into the the plotSpectra function.
# Copy and edit the function to only extract data. Do that for multiple models and plot

# 2. Changes in size-spectra relative to calibration period. Fig. 4 from Blanchard et al 2012, but lines don’t represent countries but species. Extract numbers at size. Plot difference. Check the @n for spectra, don't know which one for yield

# Plot difference in spectra between warming and constant. 

# Create getSpectra function..
getSpectra <- function(object){
  time_range  <- max(as.numeric(dimnames(object@n)$time))
  time_elements <- get_time_elements(object, time_range)
  n <- apply(object@n[time_elements, , ,drop = FALSE], c(2, 3), mean)
  species <- balticParams$species
  n <- n[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
  power <- 1
  n <- sweep(n, 2, params@w^power, "*")
  specDat <- data.frame(w = rep(as.numeric(dimnames(n)$w), length(species)),
                        n = c(as.numeric(n[1, ]), 
                              as.numeric(n[2, ]),
                              as.numeric(n[3, ])),
                        species = rep(species, each = 100))
  return(specDat)
}

# Testing it's correct
# plotSpectra(m6, plankton = F, benthos = F, algae = F) + 
#   geom_line(data = subset(getSpectra(m6), n > 0), aes(w, n, group = species), linetype = 2, color = "red") +
#   NULL


# Compare spectra in m5 and m6
warmSpec <- getSpectra(m6)
warmSpec$Scenario <- "Warming"

conSpec <- getSpectra(m5)
conSpec$Scenario <- "Constant"

plotSpec <- rbind(conSpec, warmSpec)

# Plot separately
plotSpec %>% filter(n > 0) %>% 
  ggplot(., aes(w, n, linetype = Scenario, color = Scenario)) + 
  geom_line(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = col[c(2,4)]) +
  labs(y = "Warm/Constant") +
  facet_wrap(~species, scales = "free") +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  scale_y_log10() +
  scale_x_log10()  +
  NULL

# Plot ratio instead
spec_df_wide <- data.frame(w = conSpec$w,
                           n_diff = warmSpec$n / conSpec$n,
                           Species = conSpec$species)

spec_df_wide %>% filter(w < max(vbgedat$Weight_g)) %>% # Remove the large size-classes
  ggplot(., aes(w, n_diff, color = Species)) + 
  geom_hline(yintercept = 1, col = "red") +
  geom_line(size = 1.5, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE) +
  labs(y = "Warm/Constant", x = "W (g)") +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  scale_y_log10() +
  scale_x_log10()  +
  NULL


# F. VARY EFFORT GIVEN WARMING/NO WARMING ==========================================
# Here I should add a layer showing effort during calibration period instead - with and without warming
# I will do 5 sets of effort: 0.75*FMSY, 1.25FMSY, 1.5FMSY, 1.75*FSMY
# The control is no warming + FMSY

# FSMY effort
tail(projectEffort_ct, 1)[1]

# In this scenario, replace FMSY with the effort in the calibration

# 0.75FMSY
projectEffort_ct_075FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_075FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 0.75
projectEffort_ct_075FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 0.75
projectEffort_ct_075FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 0.75

# FMSY
projectEffort_ct_FMSY <- data.frame(projectEffort_ct)

# 1.25FMSY
projectEffort_ct_125FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_125FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 1.25
projectEffort_ct_125FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 1.25
projectEffort_ct_125FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 1.25

# 1.5FMSY
projectEffort_ct_15FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_15FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 1.5
projectEffort_ct_15FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 1.5
projectEffort_ct_15FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 1.5

# 1.75FMSY
projectEffort_ct_175FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_175FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 1.75
projectEffort_ct_175FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 1.75
projectEffort_ct_175FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 1.75

# Constant temperature models
m_cons_075 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_075FMSY),
                      temperature = rep(10, nrow(projectEffort_ct_075FMSY)),
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_cons_FMSY <- project(params3_upd, 
                       dt = 0.1,
                       effort = as.matrix(projectEffort_ct_FMSY),
                       temperature = rep(10, nrow(projectEffort_ct_FMSY)),
                       diet_steps = 10,
                       t_max = t_max,
                       t_ref = 10) 

m_cons_125 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_125FMSY),
                      temperature = rep(10, nrow(projectEffort_ct_125FMSY)),
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_cons_15 <- project(params3_upd, 
                     dt = 0.1,
                     effort = as.matrix(projectEffort_ct_15FMSY),
                     temperature = rep(10, nrow(projectEffort_ct_15FMSY)),
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = 10) 

m_cons_175 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_175FMSY),
                      temperature = rep(10, nrow(projectEffort_ct_175FMSY)),
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

# Empirical temperature models
m_warm_075 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_075FMSY),
                      temperature = projectTemp$temperature,
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_warm_FMSY <- project(params3_upd, 
                       dt = 0.1,
                       effort = as.matrix(projectEffort_ct_FMSY),
                       temperature = projectTemp$temperature,
                       diet_steps = 10,
                       t_max = t_max,
                       t_ref = 10) 

m_warm_125 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_125FMSY),
                      temperature = projectTemp$temperature,
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_warm_15 <- project(params3_upd, 
                     dt = 0.1,
                     effort = as.matrix(projectEffort_ct_15FMSY),
                     temperature = projectTemp$temperature,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = 10) 

m_warm_175 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_175FMSY),
                      temperature = projectTemp$temperature,
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 


# Extract size-spectra from all models
# Starting with reference:
ref_spec <- getSpectra(m_cons_FMSY)

# Constant temperature
spec_m_cons_075 <- getSpectra(m_cons_075)
spec_m_cons_075$Fm <- "0.75*FMSY"
spec_m_cons_075$Temperature <- "Constant"
spec_m_cons_075$n_rel <- spec_m_cons_075$n / ref_spec$n

spec_m_cons_125 <- getSpectra(m_cons_125)
spec_m_cons_125$Fm <- "1.25*FMSY"
spec_m_cons_125$Temperature <- "Constant"
spec_m_cons_125$n_rel <- spec_m_cons_125$n / ref_spec$n

spec_m_cons_15 <- getSpectra(m_cons_15)
spec_m_cons_15$Fm <- "1.50*FMSY"
spec_m_cons_15$Temperature <- "Constant"
spec_m_cons_15$n_rel <- spec_m_cons_15$n / ref_spec$n

spec_m_cons_175 <- getSpectra(m_cons_175)
spec_m_cons_175$Fm <- "1.75*FMSY"
spec_m_cons_175$Temperature <- "Constant"
spec_m_cons_175$n_rel <- spec_m_cons_175$n / ref_spec$n

# Varying temperature
spec_m_warm_075 <- getSpectra(m_warm_075)
spec_m_warm_075$Fm <- "0.75*FMSY"
spec_m_warm_075$Temperature <- "Warming"
spec_m_warm_075$n_rel <- spec_m_warm_075$n / ref_spec$n

spec_m_warm_FMSY <- getSpectra(m_warm_FMSY)
spec_m_warm_FMSY$Fm <- "1.00*FMSY"
spec_m_warm_FMSY$Temperature <- "Warming"
spec_m_warm_FMSY$n_rel <- spec_m_warm_FMSY$n / ref_spec$n

spec_m_warm_125 <- getSpectra(m_warm_125)
spec_m_warm_125$Fm <- "1.25*FMSY"
spec_m_warm_125$Temperature <- "Warming"
spec_m_warm_125$n_rel <- spec_m_warm_125$n / ref_spec$n

spec_m_warm_15 <- getSpectra(m_warm_15)
spec_m_warm_15$Fm <- "1.50*FMSY"
spec_m_warm_15$Temperature <- "Warming"
spec_m_warm_15$n_rel <- spec_m_warm_15$n / ref_spec$n

spec_m_warm_175 <- getSpectra(m_warm_175)
spec_m_warm_175$Fm <- "1.75*FMSY"
spec_m_warm_175$Temperature <- "Warming"
spec_m_warm_175$n_rel <- spec_m_warm_175$n / ref_spec$n

# Rbind all data
warmFishSpectra <- rbind(spec_m_warm_075,
                         spec_m_warm_FMSY,
                         spec_m_warm_125,
                         spec_m_warm_15,
                         spec_m_warm_175,
                         spec_m_cons_075,
                         spec_m_cons_125,
                         spec_m_cons_15,
                         spec_m_cons_175)

warmFishSpectra %>% filter(w < max(vbgedat$Weight_g) & n_rel > 0) %>% # Remove the large size-classes
  ggplot(., aes(w, n_rel, color = Fm)) + 
  facet_grid(Temperature~species) +
  geom_hline(yintercept = 1, col = "black", linetype = "dashed") +
  geom_line(size = 1.7, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE) +
  labs(y = "Warm/Constant", x = "W (g)") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  ylim(0.75, 1.35) +
  scale_x_log10() +
  NULL

# Different layout
warmFishSpectra %>% 
  filter(w < max(vbgedat$Weight_g) & n_rel > 0 & Fm %in% c("0.75*FMSY",
                                                           #"1.25*FMSY",
                                                           "1.50*FMSY")) %>% # Remove the large size-classes
  ggplot(., aes(w, n_rel, color = Fm, linetype = Temperature)) + 
  facet_wrap(~species, scales = "free_x") +
  geom_hline(yintercept = 1, col = "red", size = 0.5) +
  geom_line(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = col[c(1,4)]) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  labs(y = "# / FMSY & Constant Temp.", x = "W (g)") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  ylim(0.75, 1.35) +
  scale_x_log10() +
  NULL


# R color brewer colors:
warmFishSpectra %>% 
  filter(w < max(vbgedat$Weight_g) & n_rel > 0 & Fm %in% c("0.75*FMSY",
                                                           #"1.25*FMSY",
                                                           "1.50*FMSY")) %>% # Remove the large size-classes
  ggplot(., aes(w, n_rel, color = Temperature, linetype = Fm)) + 
  facet_wrap(~species, scales = "free_x") +
  geom_hline(yintercept = 1, color = "black", size = 0.7) +
  geom_line(size = 1.4, alpha = 0.8) +
  scale_color_manual(values = c(brewer.pal(n = 8, name = "RdBu")[8], brewer.pal(n = 8, name = "RdBu")[1])) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  labs(y = "# / FMSY & Constant Temp.", x = "W (g)") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  ylim(0.75, 1.35) +
  scale_x_log10() +
  NULL

# Interesting that fishing causes the accumulation of cod around 100g. Plot diet and compare 0.75*FMSY and 1.5*FMSY

#**** Plot Yield at 2050 ===========================================================
# 3. Yield. How about a simple value? Three columns (species), Three rows(F), two lines(warming or constant?)
# Create data frames of all models

# Warming
Yield_m_warm_075 <- data.frame(getYield(m_warm_075))[117, ]
Yield_m_warm_075$Temperature <- "Warming"
Yield_m_warm_075$Fm <- "0.75*FMSY"

Yield_m_warm_FMSY <- data.frame(getYield(m_warm_FMSY))[117, ]
Yield_m_warm_FMSY$Temperature <- "Warming"
Yield_m_warm_FMSY$Fm <- "1.00*FMSY"

Yield_m_warm_125 <- data.frame(getYield(m_warm_125))[117, ]
Yield_m_warm_125$Temperature <- "Warming"
Yield_m_warm_125$Fm <- "1.25*FMSY"

Yield_m_warm_15 <- data.frame(getYield(m_warm_15))[117, ]
Yield_m_warm_15$Temperature <- "Warming"
Yield_m_warm_15$Fm <- "1.50*FMSY"

Yield_m_warm_175 <- data.frame(getYield(m_warm_175))[117, ]
Yield_m_warm_175$Temperature <- "Warming"
Yield_m_warm_175$Fm <- "1.75*FMSY"

# Constant temperature
Yield_m_cons_075 <- data.frame(getYield(m_cons_075))[117, ]
Yield_m_cons_075$Temperature <- "Constant"
Yield_m_cons_075$Fm <- "0.75*FMSY"

Yield_m_cons_FMSY <- data.frame(getYield(m_cons_FMSY))[117, ]
Yield_m_cons_FMSY$Temperature <- "Constant"
Yield_m_cons_FMSY$Fm <- "1.00*FMSY"

Yield_m_cons_125 <- data.frame(getYield(m_cons_125))[117, ]
Yield_m_cons_125$Temperature <- "Constant"
Yield_m_cons_125$Fm <- "1.25*FMSY"

Yield_m_cons_15 <- data.frame(getYield(m_cons_15))[117, ]
Yield_m_cons_15$Temperature <- "Constant"
Yield_m_cons_15$Fm <- "1.50*FMSY"

Yield_m_cons_175 <- data.frame(getYield(m_cons_175))[117, ]
Yield_m_cons_175$Temperature <- "Constant"
Yield_m_cons_175$Fm <- "1.75*FMSY"

# Combine to single data frame
allYield <- rbind(Yield_m_cons_075,
                  Yield_m_cons_FMSY,
                  Yield_m_cons_125,
                  Yield_m_cons_15,
                  Yield_m_cons_175,
                  Yield_m_warm_075,
                  Yield_m_warm_FMSY,
                  Yield_m_warm_125,
                  Yield_m_warm_15,
                  Yield_m_warm_175)

# Make data frame long
allYield_long <- allYield %>% gather(Species, Yield, 1:3)

# Create relative yield (to Yield_m_cons_FMSY)
Yield_m_cons_FMSY

allYield_long$relYield <- 0

allYield_long$relYield <- ifelse(allYield_long$Species == "Cod",
                                 allYield_long$Yield / Yield_m_cons_FMSY$Cod,
                                 allYield_long$relYield)

allYield_long$relYield <- ifelse(allYield_long$Species == "Sprat",
                                 allYield_long$Yield / Yield_m_cons_FMSY$Sprat,
                                 allYield_long$relYield)

allYield_long$relYield <- ifelse(allYield_long$Species == "Herring",
                                 allYield_long$Yield / Yield_m_cons_FMSY$Herring,
                                 allYield_long$relYield)


ggplot(allYield_long, aes(Species, relYield, color = Temperature, shape = Fm)) +
  geom_jitter(size = 4, alpha = 0.6, width = 0.1, height = 0, stroke = 1.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~Species, scales = "free") +
  scale_color_manual(values = c(brewer.pal(n = 8, name = "RdBu")[8], brewer.pal(n = 8, name = "RdBu")[1])) +
  labs(y = "Relative Yield", x = "Species") +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 3/4) +
  NULL



# G. TEST EXTREME KAPPA AND LAMBA VALUES ===========================================
# Project with temperatue and effort varying through time
# This is equivalent to the m6-model with empirical temperature (centered to 10) and
# historical effort

params4_upd <- params3_upd

params4_upd@lambda <- 2.13
params4_upd@r_pp <- 200

# 1. How can we change r_pp when it's not in @params? Only hardwired changes?

# 2. Is r * (aw^b) the same as r*a*w^b ?

w <- seq(1, 10, 1)
a <- 0.1
b <- 0.7
r <- 4

r*a*w^b
r*(a*w^b)


params_test <- MizerParams(params4_upd@species_params,
                           kappa_ben = kappa_ben,
                           kappa = kappa,
                           w_bb_cutoff = w_bb_cutoff,
                           w_pp_cutoff = w_pp_cutoff,
                           r_pp = 4,
                           r_bb = r_bb)

params_test2 <- MizerParams(params4_upd@species_params,
                            kappa_ben = kappa_ben,
                            kappa = kappa,
                            w_bb_cutoff = w_bb_cutoff,
                            w_pp_cutoff = w_pp_cutoff,
                            r_pp = 999999,
                            r_bb = r_bb)

m7a <- project(params_test2, 
               dt = 0.1,
               effort = projectEffort_ct,
               temperature = projectTemp$temperature,
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10) 

m7b <- project(params_test, 
               dt = 0.1,
               effort = projectEffort_ct,
               temperature = projectTemp$temperature,
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10) 


plotBiomass(m7a)
plotBiomass(m7b)

tail(getYield(m7a))
tail(getYield(m7b))

str(m7@params)
str(m7@params@kappa)

m7@params@species_params
m7@params@kappa
m7@params@lambda
m7@params@rr_pp
m7@params@r_pp

# Here kappa is 6 and lambda is 2.13
# How much are kappa & lambda predicted to change according to Barnes?

# First delta-temp in time series
max(projectTemp$temperature) - min(projectTemp$temperature)

# How much does kappa and lambda change from default when delta T is 2.5?
# Lambda

# Lambda = -1.175 -0.002*T

df <- data.frame(temp  = c(10, 12),
                 m     = seq(1e-2, to = 1, length.out = 40),
                 log_a = c(9.8, 9.7),# 9.8-(2*0.045)
                 b     = c(-1.1, -1.104)) # -1.1-(0.002*2)


df$B <- df$log_a + log10(df$m)*df$b

ggplot(df, aes(m, B, color = factor(temp))) +
  geom_line() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  NULL


