#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.04.16: Max Lindmark
#
# Code for calibrating the model
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
#col <- viridis(n = 5)
col <- colorRampPalette(brewer.pal(5, "Dark2"))(5)
#col <- viridis(n = 5, option = "cividis")
#col <- viridis(n = 5, option = "magma")

# Set some standard parameters (some means overwriting existing parameters that were inherited from the NS model)
# 
r_pp <-  4
r_bb <-  4
w_bb_cutoff <- 20
w_pp_cutoff <- 1

balticParams$avail_PP[which(balticParams$species == "Cod")] <- 0.5
balticParams$avail_PP[which(balticParams$species == "Herring")] <- 0.5
balticParams$avail_PP[which(balticParams$species == "Sprat")] <- 1

balticParams$avail_BB[which(balticParams$species == "Cod")] <- 0.5
balticParams$avail_BB[which(balticParams$species == "Herring")] <- 0.5
balticParams$avail_BB[which(balticParams$species == "Sprat")] <- 0

# Set default predation window (leads to less reliance on background resource compared to NS parameters)
# balticParams$beta <- 100
# balticParams$sigma <- 1.3

# Using Reum et al 2019 Oikos parameters for forage fish, but specific value for cod (this study)
balticParams
balticParams$beta[1] <- 426
balticParams$beta[2] <- 1000
balticParams$beta[3] <- 1000
balticParams$sigma[1] <- 5.6
balticParams$sigma[2] <- 2
balticParams$sigma[3] <- 2

# Lower ereprop (default is one). 0.1 is not an uncommon value (Hartvig et al 2011)
#balticParams$erepro <- 0.1
balticParams$erepro <- 0.05 * balticParams$w_inf^(-0.5)

# Create effort vector
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])


#** 1. Find starting value for kappa ===============================================
# Given the default model, what should kappa be to get ssb in the same order of magnitude? This is just an iterative process, since well opimize r_max to minimize residual sum of squares (RSS) between predicted and observed SSB
dt <- 0.2

# These are the values I choose (lowest kappa with coexistence).
kappa_ben <- 50
kappa <- 50

# When using the new default for h, we need to use the exponent for the length-weight
# relationship, else it defaults to 3.
balticParams$b <- unique(empiri$b)
balticParams$t0 <- unique(empiri$t_0)

# Create mizerParam object
params <- MizerParams(balticParams,
                      kappa_ben = kappa_ben,
                      kappa = kappa,
                      w_bb_cutoff = w_bb_cutoff,
                      w_pp_cutoff = w_pp_cutoff,
                      r_pp = r_pp,
                      r_bb = r_bb)

params@species_params$h
params@species_params$alpha
params@n
params@p

# Project model
t_max <- 600

m1 <- project(params,
              temperature = rep(10, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

plot(m1)
plotSpectra(m1, algae = FALSE)
plotFeedingLevel(m1)
m1@params@linetype <- rep("solid", 6)

# Check dynamics and density
plotBiomass(m1) + theme_classic(base_size = 14)

#**** Check growth =======================================================================
# Growth rates are extremely low
plotGrowthCurves(m1, max_age = 15) + 
  geom_point(data = subset(vbgedat, AgeRings < 16), 
             aes(AgeRings, Weight_g), size = 3, fill = "gray50", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_hline(data = empiri, aes(yintercept = w_mat), 
             color = "black", size = 0.8, linetype = 2) +
  geom_line(data = subset(empiri, age < 16), aes(age, weight), 
            color = col[2], size = 2, linetype = "twodash") +
  geom_line(aes(x = Age, y = value), 
            color = col[1], size = 2.) +
  facet_wrap(~ Species, scales = "free", ncol = 3) +
  guides(color = FALSE, linetype = FALSE) +
  theme_classic(base_size = 12) + 
  theme(aspect.ratio = 3/4) +
  NULL
