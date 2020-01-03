#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.04.16: Max Lindmark
#
# Code for testing resource activation energies
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


#** 1.DEFINE PARAMS ================================================================
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


#** 2.PROJECT WITH DEFAULTS ========================================================
#--- Now try and project model with only defaults, i.e. activation energy 0 and the scalar thus = 1
t_max <- 27

params@species_params

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


#** 3.CHANGE EA ========================================================
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


#** 3.CHANGE EA AND TEMPERATURE ====================================================
#--- Change temperature and set all activation energies to 0 except for resources
# Create mizerParam object

params2@species_params

params_new <- balticParams

params_new$ea_met <- 0
params_new$ea_int <- 0
params_new$ea_mat <- 0
params_new$ea_mor <- 0

params_new$ca_met <- 0
params_new$ca_int <- 0
params_new$ca_mat <- 0
params_new$ca_mor <- 0

params4 <- MizerParams(params_new,
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb,
                       ea_gro = 0.63,
                       ea_car = -0.63)

str(params4)
params4@species_params

m4 <- project(params4,
              temperature = rep(15, t_max),
              dt = dt,
              effort = effort,
              diet_steps = 10,
              t_max = t_max) 

str(m4@carTempScalar)

tail(getSSB(m4), 2)
# > tail(getSSB(m4), 2)
# sp
# time            Cod       Sprat     Herring
# 26 0.000003021371 0.009925334 0.001357944
# 27 0.000003224152 0.011417240 0.001461222

# > tail(getSSB(m2), 2)
# sp
# time          Cod    Sprat   Herring
# 26 0.0003235197 4.582229 0.1581226
# 27 0.0003894221 5.236273 0.1945952

# > tail(getSSB(m1), 2)
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


# # TO DO: Test with time-varying temperature! First recap how that is constructed... then compare to a scalar that i calculate straight from the tempFun. In both cases, plot the scalar!
# params8 <- MizerParams(balticParams,
#                        kappa_ben = kappa_ben,
#                        kappa = kappa,
#                        w_bb_cutoff = w_bb_cutoff,
#                        w_pp_cutoff = w_pp_cutoff,
#                        r_pp = r_pp,
#                        r_bb = r_bb,
#                        ea_gro = 0.63,
#                        ea_car = -0.63)
# 
# params8@species_params$ea_met <- 0
# params8@species_params$ea_int <- 0
# params8@species_params$ea_mat <- 0
# params8@species_params$ea_mor <- 0
# params8@species_params
# 
# # Create a temperature-vector that is as long as t_max
# temperature <- rnorm(mean = 10, sd = 2, t_max)
# 
# m8 <- project(params8,
#               temperature = temperature,
#               dt = 0.1,
#               effort = effort,
#               diet_steps = 10,
#               t_max = t_max) 
# 
# tail(getSSB(m8), 2)
# # > tail(getSSB(m8), 2)
# # sp
# # time          Cod    Sprat   Herring
# # 26 0.0006575440 6.321939 0.2620813
# # 27 0.0008102771 6.731523 0.3198012
# 
# # Now test if I can retrieve the scalar:
# str(m8)
# 
# m8_groTempScalar <- m8@groTempScalar
# 
# str(m8_groTempScalar)
# class(m8_groTempScalar)
# head(m8_groTempScalar)
# dim(m8_groTempScalar)
# 
# # Ok, so I want the columns now (270, for each iteration). There are 10 identical column names (dim names) because dt = 0.1.
# 
# # Now plot scalar as a function of temperature (large points):
# plot(rep(temperature, each = 1/dt), m8_groTempScalar[1, ], cex = 2)
# 
# # And the equation:
# fun_groTempScalar <- tempFun(temperature = temperature, 
#                              t_ref = 10, Ea = 0.63, c_a = 0, 
#                              w = 1)
# 
# str(fun_groTempScalar)
# class(fun_groTempScalar)
# 
# # Add points to the plot
# points(rep(temperature, each = 1/dt),
#        rep(fun_groTempScalar[1,], each = 10), 
#        col = "red", pch = 16)


# A. TESTING WITH TIME-VARYING TEMPERATURE =========================================
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
func <- 
  getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getGrowth.R", 
         ssl.verifypeer = FALSE)
eval(parse(text = func))

# Load function for extracting mean weight by species
func <- 
  getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getSpeciesMeanWeight.R", 
         ssl.verifypeer = FALSE)
eval(parse(text = func))

# Load function for extracting raincloud plot
func <- 
  getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/raincloudPlot.R", 
         ssl.verifypeer = FALSE)
eval(parse(text = func))


#**** Read in parameters and data ==================================================
# Read in params object
params <- readRDS("baltic/params/mizer_param_calib.rds")

# Read in activation energy data frame
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

# Define temperature-scenarios
consTemp <- projectTemp$temperature

# The time series starts in 1914. From 1997 (mid point in calibration time), we want
# to fix the temperature at the mean of the calibration, i.e. t_ref. This insures
# we get comparable starting values for the models so that temperature is the only "treatment"
start <- 1997-1914
consTemp[start:137] <- t_ref

# Plot
col <- RColorBrewer::brewer.pal("Dark2", n = 5)
col <- RColorBrewer::brewer.pal("Set1", n = 3)[1:2]

tempScen <- data.frame(Temperature = c(consTemp, projectTemp$temperature),
                       Scenario = rep(c("no warming", "warming"), each = length(consTemp)),
                       Year = 1:length(consTemp) + 1913)

ggplot(tempScen, aes(Year, Temperature, color = Scenario, linetype = Scenario)) +
  geom_line(alpha = 0.8, size = 1.4) +
  theme_classic(base_size = 25) +
  scale_color_manual(values = rev(col)) +
  theme(legend.position=c(.2,.75),
        aspect.ratio = 3/4) +
  NULL

#ggsave("baltic/figures/supp/temperature_scenarios.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


# B. SIMULATE TEMP-DRIVEN CHANGE IN SIZE-AT-AGE ====================================
# for-loop to take random samples for distributions representing activation energies
# Then compare that to a projection with a constant temperature

#**** Project without temp (reference) =============================================
ref <- project(params, 
               dt = dt,
               effort = projectEffort_m,
               temperature = consTemp,
               diet_steps = 10,
               t_max = t_max)   


refGrowth <- getGrowth(ref)
refMeanWeight <- getMeanWeight(ref)
refSpeciesMeanWeight <- getSpeciesMeanWeight(ref)[nrow(projectEffort_m), ] 


#**** Barnes - with resource - no physiological scaling ============================
# NO LOOP NEEDED YET
# sim <- 1:200

t <- c()
tt <- c()
groj <- c()
growth <- c()
# data_list_with_res_no_phys_barn <- list()
# mean_weight_list_with_res_no_phys_barn <- list()

# for (i in sim) {

t <- params@species_params

t$ea_met <- 0
t$ea_int <- 0
t$ea_mor <- 0

tt <- MizerParams(t, 
                  ea_gro = 0, #0.41, #0.43, # ea$gro[i],
                  ea_car = -0.4, #-0.41, # 0,  # -0.41, # ea$car[i], # -ea$gro[i] 
                  kappa_ben = kappa_ben,
                  kappa = kappa,
                  w_bb_cutoff = w_bb_cutoff,
                  w_pp_cutoff = w_pp_cutoff,
                  r_pp = r_pp,
                  r_bb = r_bb,
                  t_ref = t_ref)

test_temp <- projectTemp$temperature
# test_temp[length(test_temp)- 2] <- 
# test_temp[length(test_temp)- 1] <- 20000
# test_temp[length(test_temp)] <- 20000
# test_temp[(length(test_temp)/2)] <- 20000
# test_temp[20:60] <- 20000
test_temp[(length(test_temp)/2):length(test_temp)] <- 20000000000

tt@ea_gro
tt@ea_car
tt@species_params$ea_mor

proj <- project(tt, 
                dt = dt,
                effort = projectEffort_m,
                temperature = test_temp, # projectTemp$temperature,
                #temperature = projectTemp$temperature,
                diet_steps = 10) # ,
#t_max = t_max)   

plotBiomass(proj)













