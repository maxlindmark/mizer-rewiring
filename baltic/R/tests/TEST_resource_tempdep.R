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
# D. Validate with time series
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


# C. CALIBRATE MODEL ================================================================
# This will help seeing the lines..
update_geom_defaults("line", list(size = 1.75))
col <- colorRampPalette(brewer.pal(5, "Dark2"))(5)

# Set defaults for regeneration and size ranges of background resources
r_pp <-  4
r_bb <-  4
w_bb_cutoff <- 20
w_pp_cutoff <- 1

# Set allometric erepro
balticParams$erepro <- 0.05 * balticParams$w_inf^(-0.5)

# Create effort vector
effort = c(Cod = balticParams$AveEffort[1], 
           Herring = balticParams$AveEffort[3], 
           Sprat = balticParams$AveEffort[2])


#** 1. Find starting value for kappa ===============================================
dt <- 0.1

# These are the values I choose (lowest kappa with coexistence, found them iteratively).
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

# Check parameter are there (scalars are stored in the sim, not param object)
#str(params) # parameters, yes

#--- Now try and project model with only defaults, i.e. activation energy 0 and the scalar thus = 1
t_max <- 27

m1 <- project(params,
              temperature = rep(10, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

# Check scalars are there
str(m1) # parameters, yes

tail(getSSB(m1), 2)
#> tail(getSSB(m1), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003235197 4.582229 0.1581226
# 27 0.0003894221 5.236273 0.1945952


#--- Now try changing the activation energy of the resources, but keeping the temperature at reference...
# Create mizerParam object
params2 <- MizerParams(balticParams,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb,
                       ea_gro = 0.63,
                       ea_car = -0.63)

m2 <- project(params2,
              temperature = rep(10, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

# The parameters where updated!
str(m2@params)

tail(getSSB(m2), 2)
#> tail(getSSB(m2), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003235197 4.582229 0.1581226
# 27 0.0003894221 5.236273 0.1945952

# Looks good because it's the same as m1... now try and change temperature as well



#--- Change temperature
# Create mizerParam object
params3 <- MizerParams(balticParams,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb,
                       ea_gro = 0.63,
                       ea_car = -0.63)

m3 <- project(params3,
              temperature = rep(12, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

tail(getSSB(m3), 2)
# > tail(getSSB(m3), 2)
# sp
# time          Cod    Sprat    Herring
# 26 0.0001177767 2.634839 0.09635425
# 27 0.0001401350 3.120203 0.11798553

# Something changed.. now we need to isolate the effect of temperature on resoruces and see if it's different from m2 (no temperature)


#--- Change temperature but set all activation energies to 0 except for resources
# Create mizerParam object
params4 <- MizerParams(balticParams,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb,
                       ea_gro = 0.63,
                       ea_car = -0.63)

params4@species_params$ea_met <- 0
params4@species_params$ea_int <- 0
params4@species_params$ea_mat <- 0
params4@species_params$ea_mor <- 0

params4@species_params

m4 <- project(params4,
              temperature = rep(12, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 


tail(getSSB(m4), 2)
# > tail(getSSB(m4), 2)
# sp
# time          Cod     Sprat    Herring
# 26 3.157357e-05 0.6280222 0.02527125
# 27 3.674400e-05 0.8162829 0.03101256

tail(getSSB(m1), 2)
#> tail(getSSB(m1), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003235197 4.582229 0.1581226
# 27 0.0003894221 5.236273 0.1945952

# Ok, biomasses are changing at 12C also when only the resource-parameters change, meaning they have an effect.
# If the resource activation energies are set to 0, I should get the same SSB as in m1

params5 <- MizerParams(balticParams,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb,
                       ea_gro = 0,
                       ea_car = 0)

params5@species_params$ea_met <- 0
params5@species_params$ea_int <- 0
params5@species_params$ea_mat <- 0
params5@species_params$ea_mor <- 0
params5@species_params

m5 <- project(params5,
              temperature = rep(12, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

tail(getSSB(m5), 2)
# > tail(getSSB(m5), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003235197 4.582229 0.1581226
# 27 0.0003894221 5.236273 0.1945952

tail(getSSB(m1), 2)
# > tail(getSSB(m1), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003235197 4.582229 0.1581226
# 27 0.0003894221 5.236273 0.1945952

# And they are in fact the same.


# Now turn off ind-rates activation energies. Increase temp. What is the new rr_pp and cc_pp (and rr_bb and cc_bb)? Do I get the same results as when I simply define kappa and lambda to those values?

# I will test with 15C. With an activation energy of 0.63 and -0.63, I should get the following scalar:.......
tempFun <- function(temperature, t_ref, Ea, c_a, w)
{
  temperature <- temperature + 273
  t_ref <- t_ref + 273
  temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature-(t_ref)))}) *exp((-Ea/8.617332e-5)*((1/temperature) - (1/(t_ref))))) 
  return(temperatureScalar)
}
tempFun(temperature = 15, t_ref = 10, Ea = 0.63, c_a = 0, w = 1)[1]

# .... i.e. 1.565956

# My current r_pp and r_bb are 4. 
r_pp
r_bb

# So I should get the same SSB's when I
# a) Set temperature to 15C and only having resource growth activation energies at 0.63 and the rest at 0
# b) Change the r_pp and r_bb and ALL other activation energies to 0

# (a)
params6 <- MizerParams(balticParams,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb,
                       ea_gro = 0.63,
                       ea_car = 0)

params6@species_params$ea_met <- 0
params6@species_params$ea_int <- 0
params6@species_params$ea_mat <- 0
params6@species_params$ea_mor <- 0
params6@species_params

m6 <- project(params6,
              temperature = rep(15, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

# And now fit the model r_pp and r_bb updated and not equal to 4 but 4*scalar, and also
# with ALL activation energies = 0
scal <- tempFun(temperature = 15, t_ref = 10, Ea = 0.63, c_a = 0, w = 1)[1]

params7 <- MizerParams(balticParams,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp * scal,
                       r_bb = r_bb * scal,
                       ea_gro = 0,
                       ea_car = 0)

params7@species_params$ea_met <- 0
params7@species_params$ea_int <- 0
params7@species_params$ea_mat <- 0
params7@species_params$ea_mor <- 0
params7@species_params

m7 <- project(params7,
              temperature = rep(15, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

tail(getSSB(m6), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003265373 4.711890 0.1608606
# 27 0.0003939234 5.399794 0.1986076

tail(getSSB(m7), 2)
# > tail(getSSB(m7), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003265373 4.711890 0.1608606
# 27 0.0003939234 5.399794 0.1986076

## Yey!! Exactly the same....


# TO DO: Test with time-varying temperature! First recap how that is constructed... then compare to a scalar that i calculate straight from the tempFun. In both cases, plot the scalar!

params8 <- MizerParams(balticParams,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb,
                       ea_gro = 0.63,
                       ea_car = -0.63)

params8@species_params$ea_met <- 0
params8@species_params$ea_int <- 0
params8@species_params$ea_mat <- 0
params8@species_params$ea_mor <- 0
params8@species_params

# Create a temperature-vector that is as long as t_max
temperature <- rnorm(mean = 10, sd = 2, t_max)

m8 <- project(params8,
              temperature = temperature,
              dt = 0.1,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

tail(getSSB(m8), 2)
# > tail(getSSB(m8), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0006575440 6.321939 0.2620813
# 27 0.0008102771 6.731523 0.3198012

# Now test if I can retrieve the scalar:
str(m8)

m8_groTempScalar <- m8@groTempScalar

str(m8_groTempScalar)
class(m8_groTempScalar)
head(m8_groTempScalar)
dim(m8_groTempScalar)

# Ok, so I want the columns now (270, for each iteration). There are 10 identical column names (dim names) because dt = 0.1.

# Now plot scalar as a function of temperature (large points):
plot(rep(temperature, each = 1/dt), m8_groTempScalar[1, ], cex = 2)

# And the equation:
fun_groTempScalar <- tempFun(temperature = temperature, 
                             t_ref = 10, Ea = 0.63, c_a = 0, 
                             w = 1)

str(fun_groTempScalar)
class(fun_groTempScalar)

# Add points to the plot
points(rep(temperature, each = 1/dt),
       rep(fun_groTempScalar[1,], each = 10), 
       col = "red", pch = 16)



# Lastly, test with time varying effort as well. Test that in the calibration script...












