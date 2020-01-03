# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. Analyisis of size-spectra & mortality with warming and different fishing scenario
#
# C. Plot
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

# Load function for extracting numbers-at-age
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getSpectra.R", ssl.verifypeer = FALSE)
eval(parse(text = func))

# Load function for extracting mortality-at-weight (predation)
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getMortality.R", ssl.verifypeer = FALSE)
eval(parse(text = func))


#**** Read in parameters and data ==================================================
# Read in params object
params <- readRDS("baltic/params/mizer_param_calib.rds")

# Read in params object
ea <- read.csv("baltic/params/samples_activation_energy.csv")[, 2:6]
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
                             #ea_car = mean(ea$car), # -ea$gro[i]
                             ea_car = mean(ea$car), # -ea$gro[i] 
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


# B. TEMP-DRIVEN CHANGE IN MORTALITY AND SIZE SPECTRA ==============================
# for-loop to take random samples for distributions representing activation energies
# Then compare that to a projection with a constant temperature

#**** Project reference scenario ===================================================
ref <- project(pars_no_res, 
               dt = dt,
               effort = projectEffort_m,
               temperature = consTemp,
               diet_steps = 10,
               t_max = t_max)

#**** Now do the same but with barnes temperatuer-dependence =======================

# size spectrum intercept correlatio with temperature: -0.002
# size spectrum slopes correlation with temperature: -0.045

# Assuming two degrees warming:


# Define general parameters
dt <- 0.2
t_ref <- params@t_ref
kappa_ben <- params@kappa_ben
kappa <- params@kappa
w_bb_cutoff <- 20 # Not stored in mizerParams output
w_pp_cutoff <- 1 # Not stored in mizerParams outputs
r_pp <- 4 # Not stored in mizerParams output
r_bb <- 4 # Not stored in mizerParams output

ref@params@lambda

barnes_sp_par <- pars_no_res@species_params

barnes_par <- pars_no_res <- MizerParams(barnes_sp_par, 
                                         ea_gro = 0,
                                         ea_car = 0, 
                                         kappa_ben = (kappa_ben - 0.045*2),
                                         kappa = (kappa - 0.045*2),
                                         # lambda = (ref@params@lambda - 0.002*2),
                                         # lambda_ben = (ref@params@lambda - 0.002*2),
                                         w_bb_cutoff = w_bb_cutoff,
                                         w_pp_cutoff = w_pp_cutoff,
                                         r_pp = r_pp,
                                         r_bb = r_bb,
                                         t_ref = t_ref)

barnes <- project(pars_no_res, 
                  dt = dt,
                  effort = projectEffort_m,
                  temperature = consTemp,
                  diet_steps = 10,
                  t_max = t_max)

tail(getSSB(barnes))
tail(getSSB(ref))

p1 <- plotGrowthCurves(barnes) + facet_wrap(~ Species, scales = "free", nrow = 3) + ggtitle("barnes")
p2 <- plotGrowthCurves(ref) + facet_wrap(~ Species, scales = "free", nrow = 3) + ggtitle("ref")

p1 + p2


## Potentially big difference... how to get in confidence intervals?

dat <- data.frame(readxl::read_excel("baltic/data/barnes_intercept.xlsx"))

head(dat)
dat$temp <- as.numeric(dat$temp)
dat$intercept <- as.numeric(dat$intercept)
str(dat)

plot(intercept ~ temp, data = dat)

summary(lm(intercept ~ temp, data = dat))

## Not very accurate


## For now, how to approximate the negative slope of intercept found in Barnes?
slope <- -0.045

temp <- seq(5, 25, 1)

intercept <- 10.271 + temp * slope

plot(intercept ~ temp)


script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/tempFun.R", ssl.verifypeer = FALSE)
eval(parse(text = script))


scal <- data.frame(scal = as.numeric(tempFun(temperature = temp, 
                                             t_ref = 10, 
                                             Ea = -0.41, 
                                             c_a = 0, 
                                             w = 10)),
                   temp = temp)

summary(lm(scal$scal ~ scal$temp))

plot(intercept ~ temp)
lines((scal$scal + 8.82) ~ scal$temp, col = "red")



