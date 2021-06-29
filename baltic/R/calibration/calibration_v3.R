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
# D. Plots from calibration
#
# E. Validate with time series
#
# F. Save params object for analysis
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
library(plyr)
library(dplyr)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Install and reload local mizer package
#devtools::load_all(".") # THIS DOES NOT WORK ON NEW MAC; MAYBE NOT NEEDED ANYMORE SINCE
# I PUSHED CHANGES TO THE CODE ALREADY; NO NEED TO WORK IN LOCAL LIBRAY....

# Install the specific mizer version from github
# devtools::install_github("maxlindmark/mizer-rewiring", ref = "rewire-temp") 
library(mizer)

# Print package versions
# print(sessionInfo())
# other attached packages:
# mizer_1.1          testthat_2.3.0     patchwork_0.0.1    dplyr_0.8.3        tidyr_1.0.0        
#  viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       RCurl_1.95-4.12   
# bitops_1.0-6       RColorBrewer_1.1-2 devtools_2.2.1     usethis_1.5.1      ggplot2_3.2.1


# B. READ DATA =====================================================================
#**** Species parameters ===========================================================
balticParams <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/params/species_params.csv"), sep = ";")

# Fix data format after reading
str(balticParams)
cols = c(2:4, 8) # Baltic area need special care...
balticParams[, cols] %<>% lapply(function(x) as.numeric(as.character(x)))
str(balticParams)

# Add area of Baltic SD 25-29+32 (roughly Baltic proper) as a column
balticParams$sd25.29.32_m.2 <- 2.49e+11

# Read in activation energy data frame
ea <- read.csv("baltic/params/samples_activation_energy.csv")[, 2:6]
ea <- ea %>% dplyr::rename("car" = "X.gro")


#**** VBGE growth ==================================================================
# Read in vbge predictions
vbge_pred <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/data/vbge_pred.csv"))

# This gives age as 0+, 1+ etc etc. It affects t_0 if I had fitted the data on this age.
# This means that age in the data now corresponds to the number of birthdays, rather than
# the year they are on, which fits better with the modelled weight-at-age
vbge_pred$age <- vbge_pred$age + 1

# Read in growth data 
vbge_data <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/data/vbge_data.csv"))

# This gives age as 0+, 1+ etc etc.
vbge_data$age <- vbge_data$AgeRings + 1

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


#**** Temperature time series =======================================================
# This is a preliminary data set but from the correct model. Will clean this up later.
temp_datRCP8.5 <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/data/Climate/Test_RCP8.5_from_graph.csv"), sep = ";", stringsAsFactors = FALSE)

head(temp_datRCP8.5)
str(temp_datRCP8.5)

# Calculate average temperature by year
tempDat <- data.frame(temp_datRCP8.5 %>% 
                        dplyr::filter(year < 2050.5) %>% 
                        dplyr::mutate(Year.r = factor(round(year, digits = 0))) %>% 
                        dplyr::group_by(Year.r) %>% 
                        dplyr::summarize(mean_temp = mean(rel_temp_70.99Ave)) %>% 
                        dplyr::mutate(Year.num = as.numeric(as.character(Year.r))))

tempDat %>% 
  filter(Year.num < 2003 & Year.num > 1991) %>% 
  summarize(mean_temp = mean(mean_temp))

# The temperature time series shows the deviation from the average between 1970-1999.
# I want to add a certain constant to the temperature so that the time series is around 0
# at around 10C (arbitrary temperature)

plot(tempDat$Year.r, (tempDat$mean_temp + (10 + 0.1)))

plot(tempDat$Year.r, (tempDat$mean_temp + (9.57)))

tempDat %>% 
  dplyr::mutate(rel_T = mean_temp + 10.11562) %>% 
  dplyr::filter(Year.num > 1991 & Year.num < 2003) %>% 
  dplyr::summarize(mean_T = mean(rel_T))

# I want to rescale the relative temperature by adding a constant so that the mean 
# temperature in the calibration period is 10 C. Because then the temp scalar will
# be 1 in the calibration period, because t_ref in mizer is 10
t_ref <- 10.11562


# C. CALIBRATE MODEL ================================================================
# This will help seeing the lines..
update_geom_defaults("line", list(size = 1.75))
col <- colorRampPalette(brewer.pal(5, "Dark2"))(5)

# Set defaults for regeneration and size ranges of background resources
r_pp <-  4
r_bb <-  4
w_bb_cutoff <- 20
w_pp_cutoff <- 1

# Set erepro
balticParams$erepro <- 0.01

# Create effort vector
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])


#** 1. Find starting value for kappa ===============================================
# Given the default model, what should kappa be to get ssb in the same order of magnitude? 
# This is just an iterative process, since well opimize r_max to minimize residual sum of 
# squares (RSS) between predicted and observed SSB
dt <- 0.1

# These are the values I choose (lowest kappa with coexistence, found them iteratively).
kappa_ben <- 11
kappa <- 11

# Create mizerParam object
params <- MizerParams(balticParams,
                      kappa_ben = kappa_ben,
                      kappa = kappa,
                      w_bb_cutoff = w_bb_cutoff,
                      w_pp_cutoff = w_pp_cutoff,
                      r_pp = r_pp,
                      r_bb = r_bb,
                      t_ref = t_ref)

# Fix the bug with ks=0.2*h, should be ks=0.12*h
params@species_params$ks <- 0.12 * params@species_params$h

# Recalculate the params object with the new ks to get the correct metabolism
# using species_params instead of baltic params in MizerParams
params <- MizerParams(params@species_params,
                      kappa_ben = kappa_ben,
                      kappa = kappa,
                      w_bb_cutoff = w_bb_cutoff,
                      w_pp_cutoff = w_pp_cutoff,
                      r_pp = r_pp,
                      r_bb = r_bb,
                      t_ref = t_ref)

# Project model
t_max <- 300

m1 <- project(params,
              temperature = rep(t_ref, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

# Check dynamics and density
plotBiomass(m1) + theme_classic(base_size = 14)

# Check growth
# Growth rates are low
plotGrowthCurves(m1, max_age = 15) + 
  scale_color_manual(values = rep(col[1], 3)) +
  facet_wrap(~ Species, scales = "free", ncol = 3) +
  geom_point(data = subset(vbge_data, age < 16), 
             aes(age, Weight_g), size = 2, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = vbge_pred, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(aes(x = Age, y = value), 
            color = col[1], size = 1.4) +
  geom_line(data = subset(vbge_pred, age < 16), aes(age, weight), 
            color = col[2], size = 1.4, linetype = "twodash") +
  guides(color = FALSE, linetype = FALSE) +
  theme_classic(base_size = 12) + 
  theme(aspect.ratio = 3/4) +
  NULL


#** 2. Tune growth rate ============================================================
# We start by increasing the constant in the maximum consumption rate by a factor of 1.3
# If growth is more reasonable, you may tune down kappa (higher consumption allows for 
# lower kappa to get coexistence).
params2 <- params

# Increase maximum consumption rates by a factor 1.3
params2@species_params$h <- params@species_params$h * 1.3

# Remove gamma, because it needs to be recalculated using the new h. 
params2@species_params <- subset(params2@species_params,
                                 select = -gamma)

params2_upd <- MizerParams(params2@species_params,
                           kappa_ben = kappa_ben,
                           kappa = kappa,
                           w_bb_cutoff = w_bb_cutoff,
                           w_pp_cutoff = w_pp_cutoff,
                           r_pp = r_pp,
                           r_bb = r_bb,
                           t_ref = t_ref)

# Project model
t_max <- 250

m2 <- project(params2_upd,
              dt = dt,
              temperature = rep(t_ref, t_max),
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

# Check dynamics and density
plotBiomass(m2) + theme_classic(base_size = 14)

# Check growth 
# Growth rates look slightly better, ok for now since they will change after calibrating r_max
plotGrowthCurves(m2, max_age = 15) + 
  scale_color_manual(values = rep(col[1], 3)) +
  facet_wrap(~ Species, scales = "free", ncol = 3) +
  geom_point(data = subset(vbge_data, age < 16), 
             aes(age, Weight_g), size = 2, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = vbge_pred, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(aes(x = Age, y = value), 
            color = col[1], size = 1.4) +
  geom_line(data = subset(vbge_pred, age < 16), aes(age, weight), 
            color = col[2], size = 1.4, linetype = "twodash") +
  guides(color = FALSE, linetype = FALSE) +
  theme_classic(base_size = 12) + 
  theme(aspect.ratio = 3/4) +
  NULL

# Feeding levels look ok
plotFeedingLevel(m2) + theme_classic(base_size = 14)

# Now we need to check so that the relative and absolute biomasses are realistic:
# The SSB's will be optimized in the next step.
ssb_model <- getSSB(m2)[t_max, ] * m2@params@species_params$sd25.29.32_m.2 / 1e9 
ssb_data <- balticParams$AveSpawnBiomass
ssb_model/ssb_data

# We aim for being within an order of magnitude within the observed before using the 
# optim function to find the r_max'es. Need to retune kappa


#** 2b. Re-tune kappa ==============================================================
# Tune down kappa (higher consumption allows for lower kappa to get coexistence).
# This is to get calibrated biomasses closer to observed (1 order of magnitude)
params2b <- params2_upd

# These are the values I choose (lowest kappa with coexistence, found them iteratively).
kappa_ben2 <- 9
kappa2 <- 9

# Remove gamma, because it needs to be recalculated using the new kappa. 
params2b@species_params <- subset(params2b@species_params,
                                  select = -gamma)

params2b_upd <- MizerParams(params2b@species_params,
                            kappa_ben = kappa_ben2,
                            kappa = kappa2,
                            w_bb_cutoff = w_bb_cutoff,
                            w_pp_cutoff = w_pp_cutoff,
                            r_pp = r_pp,
                            r_bb = r_bb,
                            t_ref = t_ref)

# Just checking the relationship between h and ks still is true
params2b_upd@species_params$ks / ((params2b_upd@species_params$h/1.3)*0.12)

# Project model
t_max <- 250

m2b <- project(params2b_upd,
               dt = dt,
               temperature = rep(t_ref, t_max),
               effort = effort,
               diet_steps = 10,
               t_max = t_max) 

# Check dynamics and density
plotBiomass(m2b) + theme_classic(base_size = 14)

# Now we need to check so that the relative and absolute biomasses are realistic:
# The SSB's will be optimized in the next step.
ssb_model <- getSSB(m2b)[t_max, ] * m2b@params@species_params$sd25.29.32_m.2 / 1e9 
ssb_data <- balticParams$AveSpawnBiomass
ssb_model/ssb_data

# A bit better, almost 1 order of magnitude. Now do optimization to calibrate r_max


#** 3. Optimize r_max ==============================================================
# Find r_max values that minimize the residual sum of squares between SSB from stock 
# assessment and model output 
# Below I'm sourcing that functions below are written by by Asta Audzijonyte, Jon Reum 
# and Julia Blanchard, that I've modified to better fit my model system 
# (Asta, unpublished; Reum et al 2018 Oikos, Blanchard et al 2014, J. Applied. Ecol)

# Make sure MizerParams() in FunctionsForOptim.R uses the same parameter-setup that you've decided so far!

#**** Source functions =============================================================
# See the functions script for description of their arguments and such
script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/StartVector.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/FunctionsForOptim.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# With only r_max to be optimized, and three species, this function takes 3.8 minutes on my macbook pro (2019)
# Note also that I have NOT turned off warnings from creating mizerParams objects, so they will printed in the console.
# If this code takes to long to run, I will do it in a different script and store the .Rdata

# Add tuned erepro to BalticParams
#balticParams$erepro <- params2b_upd@species_params$erepro

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
# with allometric erepro (and h-scaling factor = 1.3):
exp(optim_resultNC$`par`)
# > exp(optim_resultNC$`par`)
# [1] 0.002801949 1.160167071 0.103787593

# Now update Rmax
params3 <- params2b_upd@species_params

# This is the final set of parameters I use after calibration has been done
params3$r_max <- exp(optim_resultNC$`par`)

params3_upd <- MizerParams(params3,
                           kappa_ben = kappa_ben2,
                           kappa = kappa2,
                           w_bb_cutoff = w_bb_cutoff,
                           w_pp_cutoff = w_pp_cutoff,
                           r_pp = r_pp,
                           r_bb = r_bb,
                           t_ref = t_ref)

#params3_upd@species_params$ks / ((params3_upd@species_params$h/1.3)*0.12)

# Project model with optimized r_max
t_max <- 65

m3 <- project(params3_upd,
              dt = dt,
              temperature = rep(t_ref, t_max),
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

# Check dynamics and density
plotBiomass(m3) + theme_classic(base_size = 14)

# Check growth
plotGrowthCurves(m3, max_age = 15) + 
  scale_color_manual(values = rep(col[1], 3)) +
  facet_wrap(~ Species, scales = "free", ncol = 1) +
  geom_point(data = subset(vbge_data, age < 16), 
             aes(age, Weight_g), size = 2, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = vbge_pred, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(aes(x = Age, y = value), 
            color = col[5], size = 1.3, alpha = 0.8) +
  geom_line(data = subset(vbge_pred, age < 16), aes(age, weight), 
            color = col[4], size = 1.3, linetype = "twodash", alpha = 0.8) +
  guides(color = FALSE, linetype = FALSE) +
  theme_classic(base_size = 14) + 
  theme(aspect.ratio = 1/2) +
  NULL

# They still look OK after calibration

# SSB pred/fit
ssb_model <- getSSB(m3)[t_max, ] * m3@params@species_params$sd25.29.32_m.2 / 1e9 
ssb_data <- balticParams$AveSpawnBiomass
ssb_model/ssb_data

> ssb_model/ssb_data
Cod    Sprat  Herring 
1.038210 1.044248 1.050059 

# Seems OK!

#** 4. Check recruitment & density dependence ======================================
# Evaluate how much density dependence there is in the model from the stock-recruit relationship
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
# > rdi / rdd
# Cod     Sprat   Herring 
# 554.29910  28.47088  79.35805 

# Get RDD to rmax ratio
rdd / params3_upd@species_params$r_max
# > rdd / params3_upd@species_params$r_max
# Cod     Sprat   Herring 
# 0.9981959 0.9648764 0.9873989 
# Very close to r_max...

# Get RDI to rmax ratio
rdi / params3_upd@species_params$r_max

# Tuning r_max given the default erepro results in all species being close to their 
# rmax, i.e. their recruitment is mostly controlled by the external parameter r_max
# This is evident in that the yield~fishing curves are very flat, and it is hard to 
# overfish the species. This is not something we want. Therefore, I will tune down 
# erepro from the calibratet model. When I do that, the fit to SSB will get worse,
# but the yield curves will improve. I will do it by hand through visual inspection.
# But it can also be done using the same optimization code I used for r_max, although
# with a different error function. I have for instance used a function that minimized
# difference between FMSY from model and data (while that gives good FMSY values, yield
# curves can still be flat). I recommend the user to think about which features of the
# model should be tuned given the questions and then define an appropriate error function.

# Yield curves look pretty flat (not shown here but code is below if you want to make the 
# plots for a different model setup you can run that code here). 
# See if I can make them look more realistic by tuning erepro down
params3b_upd <- params3_upd

params3b_upd@species_params

# I found these scalars after manually tuning erepro down and by looking that the R/R_max ratio
# as well as the yield curves and fit to SSB
params3b_upd@species_params$erepro <- params3_upd@species_params$erepro * c(0.005, 0.1, 0.08) # I want RDI / Rmax to be 1:10

t_max <- 500

m3b <- project(params3b_upd,
               dt = dt,
               temperature = rep(t_ref, t_max),
               effort = effort,
               diet_steps = 10,
               t_max = t_max) 

# Check dynamics and density
plotBiomass(m3b) + theme_classic(base_size = 14)

# SSB pred/fit
ssb_model <- getSSB(m3b)[t_max, ] * m3b@params@species_params$sd25.29.32_m.2 / 1e9 
ssb_data <- balticParams$AveSpawnBiomass
ssb_model/ssb_data
# > ssb_model/ssb_data
# Cod     Sprat   Herring 
# 0.6775174 0.8466963 1.0334025 
# A bit worse fit, but "better" density dependence, i.e. not too much...

# Re-evalute how much density dependence there is in the model from the stock-recruit relationship
# Calculate the density independent recruitment (total egg production) R_{p.i} before density dependence
rdi <- getRDI(m3b@params,
              m3b@n[t_max,,],
              m3b@n_pp[t_max,],
              m3b@n_bb[t_max,],
              m3b@n_aa[t_max,],
              m3b@intTempScalar[,,(t_max/dt)],
              m3b@metTempScalar[,,(t_max/dt)])

# Calculate the flux entering the smallest size class of each species (recruitment - density dependence)
rdd <- getRDD(m3b@params,
              m3b@n[t_max,,],
              m3b@n_pp[t_max,],
              m3b@n_bb[t_max,],
              m3b@n_aa[t_max,],
              sex_ratio = 0.5,
              m3b@intTempScalar[,,(t_max/dt)],
              m3b@metTempScalar[,,(t_max/dt)])

# rdi vs rdd
rdi / rdd
# > rdi / rdd
# Cod    Sprat  Herring 
# 2.780075 3.220745 7.160041

rdd / params3b_upd@species_params$r_max
# > rdd / params3b_upd@species_params$r_max
# Cod     Sprat   Herring 
# 0.6402974 0.6895128 0.8603360

rdi / params3b_upd@species_params$r_max
# > rdi / params3b_upd@species_params$r_max
# Cod    Sprat  Herring 
# 1.780075 2.220745 6.160041


#** 5. Estimate FMSY from model - compare with assessment ==========================
# In lack of a better approach, I will just for-loop different F, extract Yield, plot
# over F. I will increase each species F separately, keeping the others at their mean

F_range <- seq(0, 1.5, 0.02) # Can increase later, becomes too slow now
t_max <- 200

#**** Cod ==========================================================================
# Mean F in calibration time
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])

# Create empty data holder
codFmsy <- c()
Y <- c()
Fm <- c()
t <- c()

for(i in F_range) {
  
  effort[1] <- i
  
  t <- project(params3b_upd, 
               dt = 0.1,
               effort = effort, 
               temperature = rep(params3b_upd@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max) 
  
  Y <- mean(data.frame(getYield(t))$Cod[(t_max-20):t_max])
  Fm <- i
  t <- cbind(Y, Fm)
  
  codFmsy <- data.frame(rbind(t, codFmsy))
  
}

codFmsy$Species <- "Cod"

ggplot(codFmsy, aes(Fm, Y)) + geom_line()


#**** Sprat ========================================================================
# Mean F in calibration time
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])

# Create empty data holder
sprFmsy <- c()
Y <- c()
Fm <- c()
t <- c()

for(i in F_range) {
  
  effort[3] <- i
  
  t <- project(params3b_upd, 
               dt = 0.1,
               effort = effort, 
               temperature = rep(params3b_upd@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max) 
  
  Y <- mean(data.frame(getYield(t))$Sprat[(t_max-20):t_max])
  Fm <- i
  t <- cbind(Y, Fm)
  
  sprFmsy <- data.frame(rbind(t, sprFmsy))
  
}

sprFmsy$Species <- "Sprat"

ggplot(sprFmsy, aes(Fm, Y)) + geom_line()


#**** Herring ======================================================================
# Mean F in calibration time
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])

# Create empty data holder
herFmsy <- c()
Y <- c()
Fm <- c()
t <- c()

for(i in F_range) {
  
  effort[2] <- i
  
  t <- project(params3b_upd, 
               dt = 0.1,
               effort = effort, 
               temperature = rep(params3b_upd@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max) 
  
  Y <- mean(data.frame(getYield(t))$Herring[(t_max-20):t_max])
  Fm <- i
  t <- cbind(Y, Fm)
  
  herFmsy <- data.frame(rbind(t, herFmsy))
  
}

herFmsy$Species <- "Herring"

ggplot(herFmsy, aes(Fm, Y)) + geom_line()


# All together
Fmsy <- rbind(codFmsy, sprFmsy, herFmsy)


# D. PLOTS FROM CALIBRATION ========================================================
#** FMSY ===========================================================================
p1 <- Fmsy %>% filter(Y > 0.01) %>% ggplot(., aes(Fm, (Y*249), color = Species)) + 
  geom_line(alpha = 0.8) +
  scale_color_manual(values = rev(col)) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "Fishing mortality [1/year]", 
       y = "Yield [1000 tonnes/year]") +
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))

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


p2 <- ggplot(asses_mod_FMSY, aes(Species, FMSY, shape = Source, fill = Source)) + 
  geom_point(size = 6, alpha = 0.7, color = "white") +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  scale_shape_manual(values = seq(21, 25)) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))


pWord1 / pWord2
ggsave("baltic/figures/supp/FMSY.png", width = 6.5, height = 6.5, dpi = 600)


#** Growth =========================================================================
p3 <- plotGrowthCurves(m3b, max_age = 15) + 
  scale_color_manual(values = rep(col[1], 3)) +
  facet_wrap(~ Species, scales = "free", ncol = 1) +
  geom_point(data = subset(vbge_data, age < 16), 
             aes(age, Weight_g), size = 2, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = vbge_pred, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(aes(x = Age, y = value), 
            color = col[5], size = 1.3, alpha = 0.8) +
  geom_line(data = subset(vbge_pred, age < 16), aes(age, weight), 
            color = col[4], size = 1.3, linetype = "twodash", alpha = 0.8) +
  guides(color = FALSE, linetype = FALSE) +
  NULL

pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12),
                                       aspect.ratio = 1/2)
ggsave("baltic/figures/supp/growth_model_data.png", width = 6.5, height = 6.5, dpi = 600)


#** Spectra and feeding level ======================================================
p4 <- spect <- plotSpectra(m3b, algae = F) + 
  scale_color_manual(values = rev(col)) +
  NULL

pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12),
                                       aspect.ratio = 3/4)

# Check feeding level
p5 <- feedlev <- plotFeedingLevel(m3b) + 
  NULL

pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12),
                                       aspect.ratio = 3/4)

pWord4 / pWord5
ggsave("baltic/figures/supp/spect_feedlev.png", width = 6.5, height = 6.5, dpi = 600)


#** SSB fit ========================================================================
# Plot predicted vs fitted SSB. Take mean of last 20 time steps
obs <- balticParams$AveSpawnBiomass
pred <- colMeans(getSSB(m3b)[c(I(dim(m3b@n)[1]-20):dim(m3b@n)[1]), ] ) * (balticParams$sd25.29.32_m.2) / (1e9)

ssb_eval <- data.frame(SSB = c(obs, pred),
                       Source = rep(c("Observed", "Predicted"), each = 3),
                       Species = rep(balticParams$species, 2))

p6 <- ggplot(ssb_eval, aes(Species, SSB, shape = Source, fill = Species)) + 
  geom_point(size = 4, alpha = 0.8) +
  scale_fill_manual(values = rev(col)) +
  scale_shape_manual(values = c(24, 21),
                     guide = guide_legend(override.aes = list(colour = "black", 
                                                              fill = "black",
                                                              size = 4))) + 
  labs(x = "", y = "Spawning stock biomass\n(1000 tonnes)") +
  guides(fill = FALSE) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12),
                                       legend.position = c(.85, .2),
                                       legend.title = element_blank(),
                                       aspect.ratio = 1)

      
ssb_eval_l <- data.frame(obs = log10(obs), pred = log10(pred), Species = balticParams$species)

p7 <- ggplot(ssb_eval_l, aes(obs, pred, fill = Species, shape = Species)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_fill_manual(values = rev(col)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  labs(x = "Log10(Observed SSB)", y = "Log10(Predicted SSB)") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12),
                                       legend.position = c(.85, .2),
                                       legend.title = element_blank(),
                                       aspect.ratio = 1)

pWord6 + pWord7
ggsave("baltic/figures/supp/SSB_fit.png", width = 6.5, height = 6.5, dpi = 600)


#** Diet ===========================================================================
p8 <- plotDietComp(m3b, prey = dimnames(m3b@diet_comp)$prey[1:5]) + 
  scale_fill_manual(values = rev(col),
                    labels = c("Cod", "Sprat", "Herring", "Plankton", "Benthos")) +
  scale_x_continuous(name = "log10 predator mass (g)", expand = c(0,0)) +
  scale_y_continuous(name = "Proportion of diet by mass (g)", expand = c(0,0)) +
  NULL

pWord8 <- p8 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12),
                                       aspect.ratio = 1,
                                       legend.position = "bottom")

ggsave("baltic/figures/supp/diet.png", width = 6.5, height = 6.5, dpi = 600)


# E. VALIDATE WITH TIME SERIES =====================================================
#** Set up time varying effort =====================================================
# First I need to set-up the matricies holding time series of temperature and effort
# Note also that t_max is taken from the effort-array when used and not the argument 
# in project (that's just how it is in mizer), and temperature needs to have the same 
# length as t_max, I here center effort to start at 1 and the actual year. 
# I go back to normal scale when plotting.

# This procedure is based on the Mizer vignette
# Read in historical F
f_history <- as.matrix(
  ssb_f %>% 
    dplyr::filter(Year > 1973 & Year < 2013) %>% 
    dplyr::select(Year, Species, Fm) %>% 
    spread(Species, Fm) %>%
    dplyr::select(Cod, Herring, Sprat)
)

rownames(f_history) <- 1974:2012  # If I want the full time series, do: 1974:2012. Now I just to 10 years prior and after

# Before projecting forward, we want to remove the transient dynamics.
# Judging from the last projection, 100 yrs burn in seems ok.
# We also start the time series at 1974 (1992-2002 for calibration)
# That means we start the simulation to get a time series free from transients 
# by starting from 1974-burnin
plotBiomass(m3b)

burnin <- 100

# Here I use the F at the first time step of the series, as in 'mizer' vignette
initial_effort <- matrix(f_history[1, ], # colMeans(f_history[1:19, ])
                         byrow = TRUE,
                         nrow = burnin,
                         ncol = ncol(f_history),
                         dimnames = list(1874:1973)) # Must match with burn-in

all_effort <- rbind(initial_effort, f_history)

# All effort until 2012
all_effort

# Create effort for projection with temperature by appending FMSY to effort data
t_future <- 2050-2012

projectEffort_fwr_df <- data.frame(all_effort[1:t_future, ])

# Add in FMSY from the MIZER model.
# Note though that the model is calibrated to FMSY from stock assessment
# However, in order to not confound the analysis with FMSY differeing between models, 
# I will add FMSY from the mizer model

# FMSY from assessment
asses_mod_FMSY

size_spect_FMSY <- asses_mod_FMSY %>% filter(Source == "Size Spectrum Model")

projectEffort_fwr_df[, 1] <- filter(size_spect_FMSY, Species == "Cod")$FMSY
projectEffort_fwr_df[, 2] <- filter(size_spect_FMSY, Species == "Herring")$FMSY
projectEffort_fwr_df[, 3] <- filter(size_spect_FMSY, Species == "Sprat")$FMSY

projectEffort_fwr <- as.matrix(projectEffort_fwr_df)
rownames(projectEffort_fwr) <- 2013:2050

projectEffort <- rbind(all_effort, projectEffort_fwr)

# Create vector of centered rownames
projectEffort_ct <- projectEffort
rownames(projectEffort_ct) <- 1:nrow(projectEffort)
projectEffort_ct

# Create data frame for plotting later where Year is a column
plotEffort <- data.frame(projectEffort)

plotEffort$Year <- as.numeric(as.character(rownames(projectEffort)))


#** Set up time varying temperature ================================================
# The temperature data is relative to a mean (1970-1999). By adding t_ref here, the mean 
# temperature in the calibration time period becomes 10C, which is both reasonable and arbitrary
tempDat$mean_temp_scaled <- tempDat$mean_temp + m3b@params@t_ref

# Now I need to match year for temperature and effort
projectEffort_ct # Centered complete effort array
projectEffort    # Non-centered complete effort array

# I need a temperature-vector with same length as t_max.
# Repeat first element in temp-data from min(plotEffort$Year):1969, as temperature starts at 1970
projectTemp_df <- tempDat %>% 
  dplyr::select(Year.r, mean_temp_scaled) %>% 
  dplyr::rename(Year = Year.r,
                temperature = mean_temp_scaled)

projectTemp_burnin <- data.frame(Year = seq(from = min(plotEffort$Year), to = 1969, by = 1),
                                 temperature = rep(projectTemp_df$temperature[1], length(seq(min(plotEffort$Year):1969))))

projectTemp <- rbind(projectTemp_burnin, projectTemp_df)
projectTemp


#** Plot effort and temperature scenarios ==========================================
# Plot effort data
col <- RColorBrewer::brewer.pal("Dark2", n = 5)

p9 <- plotEffort %>% 
  gather(Species, Effort, 1:3) %>% 
  ggplot(., aes(Year, Effort, color = Species, linetype = Species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = 1974, 
                xmax = 2012,
                ymin = 0,
                ymax = 1.4),
            fill  = "gray90") +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = 1.4),
            fill  = "gray80") +
  geom_line(size = 1.2, alpha = 0.8) +
  coord_cartesian(expand = 0) +
  scale_color_manual(values = rev(col)) +
  labs(x = "Year", y = "Fishing mortality (F)") +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12), 
                                       aspect.ratio = 3/4, 
                                       legend.position = c(.2, .85),
                                       legend.title = element_blank())

# Plot temperature data as implemented in model (i.e. with t_ref until mid-point of calibration window)
projectEffort_m <- as.matrix(projectEffort)
consTemp <- projectTemp$temperature
start <- 1997 - min(plotEffort$Year)
consTemp[start:nrow(plotEffort)] <- t_ref

tempScen <- data.frame(Temperature = c(consTemp, projectTemp$temperature),
                       Scenario = rep(c("no warming", "warming"), each = length(consTemp)),
                       Year = 1:length(consTemp) + (min(plotEffort$Year) - 1))

col <- RColorBrewer::brewer.pal("Dark2", n = 5)
col <- RColorBrewer::brewer.pal("Set1", n = 3)[1:2]

p10 <- ggplot(tempScen, aes(Year, (Temperature), color = Scenario, linetype = Scenario)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 9.2,
                ymax = 11.8),
            fill  = "gray90") +
  geom_line(alpha = 0.8, size = 1.2) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = rev(col)) +
  coord_cartesian(expand = 0) +
  ylab(expression(paste("Relative temperature [", degree*C, "]"))) +
  NULL

pWord10 <- p10 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12), 
                                         aspect.ratio = 3/4, 
                                         legend.position = c(.25, .75),
                                         legend.title = element_blank())


# Plot effort and the two temperature-series temperature
pWord9 / pWord10

ggsave("baltic/figures/supp/effort_temp.png", width = 6.5, height = 6.5, dpi = 600)


#** Project with temperature and effort varying through time =======================
# Here I want to create two models
# m4_temp: time series with temperature + parameters for temperature
# m4_cons: constant temperature (t_ref, i.e. no temperature effect!)
# NOTE! This is not a good method when comparing the models as they will
# essentially have different starting values. It's ok now though because
# I want to see if they differ in their fit to assessment data, not to understand
# the effect of temperature per se
# In the next section I will evalute their fit to assessment data

# Update params to include non-default activation energies in the species params
temp_sp_param <- params3b_upd@species_params

temp_sp_param$ea_met <- mean(ea$met)
temp_sp_param$ea_mor <- mean(ea$mor)
temp_sp_param$ea_int <- mean(ea$int)

# No temperature effects on the resource
params_n_res <- MizerParams(temp_sp_param,
                            kappa_ben = kappa_ben2,
                            kappa = kappa2,
                            w_bb_cutoff = w_bb_cutoff,
                            w_pp_cutoff = w_pp_cutoff,
                            r_pp = r_pp,
                            r_bb = r_bb,
                            ea_gro = 0, # mean(ea$gro),
                            ea_car = 0, # mean(ea$car),
                            t_ref = t_ref)

m4_noRes <- project(params_n_res, 
                    dt = 0.1,
                    effort = projectEffort_ct,
                    temperature = projectTemp$temperature,
                    diet_steps = 10) 

# Full time series including burn in
plotBiomass(m4_noRes) + theme_classic(base_size = 14)

# From 1974 and forward 
plotBiomass(m4_noRes) + 
  theme_classic(base_size = 14) +
  xlim(41, 120) +
  NULL

# With temperature effects on the resource
params_w_res <- MizerParams(temp_sp_param,
                            kappa_ben = kappa_ben2,
                            kappa = kappa2,
                            w_bb_cutoff = w_bb_cutoff,
                            w_pp_cutoff = w_pp_cutoff,
                            r_pp = r_pp,
                            r_bb = r_bb,
                            ea_gro = mean(ea$gro),
                            ea_car = mean(ea$car),
                            t_ref = t_ref)

m4_wiRes <- project(params_w_res, 
                    dt = 0.1,
                    effort = projectEffort_ct,
                    temperature = projectTemp$temperature,
                    diet_steps = 10) 

# Full time series including burn in
plotBiomass(m4_wiRes) + theme_classic(base_size = 14)

# From 1974 and forward 
plotBiomass(m4_wiRes) + 
  theme_classic(base_size = 14) +
  xlim(41, 120) +
  NULL

# No temperature effect at all
m4_consTemp <- project(params_w_res, 
                       dt = 0.1,
                       effort = projectEffort_ct,
                       temperature = rep(t_ref, nrow(projectEffort_ct)),
                       diet_steps = 10) 

# Full time series including burn in
plotBiomass(m4_consTemp) + theme_classic(base_size = 14)

# From 1974 and forward 
plotBiomass(m4_consTemp) + 
  theme_classic(base_size = 14) +
  xlim(41, 120) +
  NULL


#**** Compare with assessment data =================================================
# Observed ssb
str(ssb_f)
obs_ssb_l <- ssb_f %>%
  filter(Year > 1973) %>% 
  select(Species, Year, SSB)

obs_ssb_l$Scenario <- "Stock assessment"

# Predicted ssb - no temp dep resource
# m4_noRes@params
# str(m4_noRes@carTempScalar)
# plot(m4_noRes@carTempScalar[1, ])

pred_ssb_noResT <- data.frame(getSSB(m4_noRes))
pred_ssb_noResT$Year_ct <- as.numeric(rownames(getSSB(m4_noRes)))
pred_ssb_noResT$Year <- pred_ssb_noResT$Year_ct + (min(plotEffort$Year)-1)
pred_ssb_noResT$Scenario <- "Physio."
str(pred_ssb_noResT)

# Predicted ssb - with temp on resource
# m4_wiRes@params
# str(m4_wiRes@carTempScalar)
# plot(m4_wiRes@carTempScalar[1, ])

pred_ssb_wiResT <- data.frame(getSSB(m4_wiRes))
pred_ssb_wiResT$Year_ct <- as.numeric(rownames(getSSB(m4_wiRes)))
pred_ssb_wiResT$Year <- pred_ssb_wiResT$Year_ct + (min(plotEffort$Year)-1) 
pred_ssb_wiResT$Scenario <- "Physio. + Resource (exp.)"
str(pred_ssb_wiResT)

# Predicted ssb - no temperature at all after calibration (t_ref)
pred_ssb_cons <- data.frame(getSSB(m4_consTemp))
pred_ssb_cons$Year_ct <- as.numeric(rownames(getSSB(m4_consTemp)))
pred_ssb_cons$Year <- pred_ssb_cons$Year_ct + (min(plotEffort$Year)-1) 
pred_ssb_cons$Scenario <- "Constant temp"
str(pred_ssb_cons)

# Combine
pred_ssb <- rbind(pred_ssb_noResT, pred_ssb_wiResT, pred_ssb_cons)

# Convert to long data frame (1 obs = 1 row)
pred_ssb_l <- pred_ssb %>% 
  dplyr::filter(Year > 1973) %>%
  gather(Species, ssb_g.m2, 1:3)

# Scale up from g/m2 to 10^6 kg / Baltic
pred_ssb_l$SSB <- pred_ssb_l$ssb_g.m2 * (balticParams$sd25.29.32_m.2) / (1e9)
pred_ssb_l <- pred_ssb_l %>% select(-c(ssb_g.m2, Year_ct))

# Combine observed and predicted SSB
dat <- data.frame(rbind(obs_ssb_l, pred_ssb_l))
dat$Year <- as.integer(dat$Year)

# Plot predicted and observed ssb by species, normalize by max within species
# Reorder factor levels
dat$Scenario <- factor(dat$Scenario, levels = c("Constant temp", "Physio.", "Physio. + Resource (exp.)", "Stock assessment"))

p11 <- dat %>% filter(Year < 2012 & Year > 1970) %>% 
  ggplot(., aes(Year, SSB, linetype = Scenario, color = Scenario, alpha = Scenario)) +
  facet_wrap(~ Species, ncol = 1, scales = "free") +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = 1),
            fill  = "gray90") +
  geom_line(size = 1.5) +
  scale_linetype_manual(values = c("twodash", "dashed", "dotted", "solid")) +
  scale_color_manual(values = c(rev(col)[1:3], "gray30")) +
  scale_alpha_manual(values = c(0.8, 0.8, 0.8, 0.5)) +
  labs(y = "Spawning stock biomass (1000 tonnes)", x = "Year") +
  scale_y_continuous(expand = c(0, 0)) +
  NULL

pWord11 <- p11 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12), 
                                         aspect.ratio = 1/2)

ggsave("baltic/figures/supp/time_series_pred_ssb.png", width = 6.5, height = 6.5, dpi = 600)


#**** Calculate and plot correlation coefficients ==================================
# Since the temperature-scenarios are so similar, I'm just calculating the correlations
obs_df <- filter(dat, Scenario == "Stock assessment" & Year < 2013)
pred_wTempR_df <- filter(dat, Scenario == "Physio. + Resource (exp.)" & Year < 2013)

# For the scenario with temperature-dependent resources
cor_df <- data.frame(Year = obs_df$Year,
                     Obs = obs_df$SSB,
                     pred_wTempR = pred_wTempR_df$SSB,
                     Species = obs_df$Species)

# Calculate correlations between predictions from _resource_temp and observations
cors_con <- ddply(cor_df, c("Species"), summarise, cor = round(cor(pred_wTempR, Obs), 2))

# Plot correlation between predicted and observed
p12 <- ggplot(cor_df, aes(Obs, pred_wTempR, color = Year)) +
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 0.7) +
  geom_point(size = 2) +
  labs(y = "Predicted", x = "Observed") +
  scale_y_continuous(expand = c(0, 0)) + 
  geom_text(data = cors_con, aes(label = paste("r = ", cor, sep = "")), 
            x = c(550, 1500, 1700), y = c(60, 600, 930),
            fontface = "italic", size = 4, inherit.aes = FALSE) +
  scale_color_viridis() +
  NULL

pWord12 <- p12 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12), 
                                         legend.text = element_text(size = 8),
                                         legend.position = "bottom",
                                         aspect.ratio = 1)

ggsave("baltic/figures/supp/obs_pred_corr.png", width = 6.5, height = 6.5, dpi = 600)


# F. SAVE OBJECTS FOR ANALYSIS =====================================================
#**** Mizer params =================================================================
mizer_param_calib <- params3b_upd
str(mizer_param_calib)

saveRDS(mizer_param_calib, file = "baltic/params/mizer_param_calib.rds") 

#**** Temperature and effort vectors ===============================================
write.csv(projectEffort_ct, file = "baltic/params/projectEffort.csv")
write.csv(projectTemp, file = "baltic/params/projectTemp.csv") 

params3b_upd@species_params


# G. TEST ==========================================================================
# TEST I can get the same biomasses when I scale R_max, Kappa and gamma as when I scale only biomasses afterwards...
# t <- params3_upd@species_params
# 
# t$r_max <- params3_upd@species_params$r_max * ((balticParams$sd25.29.32_m.2) / (1e9))
# t$gamma <- params3_upd@species_params$gamma / ((balticParams$sd25.29.32_m.2) / (1e9))
# 
# tt <- MizerParams(t,
#                   kappa_ben = kappa_ben2 * (balticParams$sd25.29.32_m.2)[1] / (1e9),
#                   kappa = kappa2 * (balticParams$sd25.29.32_m.2)[1] / (1e9),
#                   w_bb_cutoff = w_bb_cutoff,
#                   w_pp_cutoff = w_pp_cutoff,
#                   r_pp = r_pp,
#                   r_bb = r_bb,
#                   t_ref = t_ref)  
# 
# ttt <- project(tt,
#                dt = dt,
#                temperature = rep(t_ref, 65),
#                effort = effort,
#                diet_steps = 10,
#                t_max = 65) 
# 
# plot(ttt)
# colMeans(getSSB(ttt)[c(I(dim(ttt@n)[1]-20):dim(m3@n)[1]), ] )  * (balticParams$sd25.29.32_m.2) / (1e9)
# colMeans(getSSB(ttt)[c(I(dim(ttt@n)[1]-20):dim(m3@n)[1]), ] )  / (1e9)
# colMeans(getSSB(ttt)[c(I(dim(ttt@n)[1]-20):dim(m3@n)[1]), ] )
# These are the means in unit m^2 * scaling factor 249
# colMeans(getSSB(m3)[c(I(dim(m3@n)[1]-20):dim(m3@n)[1]), ] ) * (balticParams$sd25.29.32_m.2) / (1e9)

# i.e. the same as when I tune parameters directly.

# TEST if they can die from fishing
# effort = c(Cod = balticParams$AveEffort[1], 
#            Herring = balticParams$AveEffort[3], 
#            Sprat = balticParams$AveEffort[2])
# 
# teffort <- effort
# 
# teffort[1] <- 1.8
# 
# test <- project(params3b_upd,
#                 dt = 0.1,
#                 effort = teffort,
#                 temperature = rep(t_ref, t_max),
#                 diet_steps = 10,
#                 t_max = t_max)
# test@effort
# plotBiomass(test)
# plot(test)