#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# TESTING EFFECTS OF RANDOM ACTIVATION ENERGY
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
# [1] mizer_1.1 testthat_2.0.0 patchwork_0.0.1 dplyr_0.8.1 tidyr_0.8.3       
# [6] viridis_0.5.1 viridisLite_0.3.0 magrittr_1.5 RCurl_1.95-4.12 bitops_1.0-6      
# [11] RColorBrewer_1.1-2 usethis_1.4.0 devtools_2.0.2 ggplot2_3.1.1  

# Load function for extracting size-at-age
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getGrowth.R", ssl.verifypeer = FALSE)
eval(parse(text = func))


#**** Read in parameters and data ==================================================
# Read in params object
params <- readRDS("baltic/params/mizer_param_calib.rds")

# Read in params object
ea <- read.csv("baltic/params/samples_activation_energy.csv")[, 2:6]

# Read in effort and temperature for projections
projectEffort <- read.csv("baltic/params/projectEffort.csv")[, 2:4]
projectTemp <- read.csv("baltic/params/projectTemp.csv")

projectEffort_m <- as.matrix(projectEffort)
rownames(projectEffort_m) <- 1:nrow(projectEffort)

# Define general parameters
t_max <- 2000
dt <- 0.2
t_ref <- 10
kappa_ben = 1
kappa = 1
w_bb_cutoff = 20
w_pp_cutoff = 1
r_pp = 4
r_bb = 4


# B. TEMP-DRIVEN CHANGE IN SIZE-AT-AGE =============================================
# for-loop to take random samples for distributions representing activation energies
# Then compare that to a projection with a constant temperature

#**** Project without temp (reference) =============================================
consTemp <- projectTemp$temperature
consTemp[] <- 10

ref <- project(params, 
               dt = dt,
               effort = projectEffort_m,
               temperature = consTemp,
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10)   

refGrowth <- getGrowth(ref)
refGrowth$scen <- "Ref"

#**** Project with temp and bad activation energies ================================
bad_par_sp <- params@species_params

bad_par_sp$ea_met <- max(ea$met)
bad_par_sp$ea_int <- min(ea$int)
bad_par_sp$ea_mor <- max(ea$mor)

bad_par <- MizerParams(bad_par_sp,
                       ea_gro = min(ea$gro),
                       ea_car = min(ea$car),
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb)

#bad_par_sp$ca_int <- -0.004 # Here we just use the fixed values
#bad_par_sp$ca_met <- 0.001 # Here we just use the fixed values

bad <- project(bad_par, 
               dt = dt,
               effort = projectEffort_m,
               temperature = projectTemp$temperature,
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10)   

badGrowth <- getGrowth(bad)
badGrowth$scen <- "Bad"

#**** Project with temp and good activation energies ===============================
good_par_sp <- params@species_params

good_par_sp$ea_met <- min(ea$met)
good_par_sp$ea_int <- max(ea$int)
good_par_sp$ea_mor <- min(ea$mor)

good_par <- MizerParams(good_par_sp,
                        ea_gro = max(ea$gro),
                        ea_car = max(ea$car),
                        kappa_ben = kappa_ben,
                        kappa = kappa,
                        w_bb_cutoff = w_bb_cutoff,
                        w_pp_cutoff = w_pp_cutoff,
                        r_pp = r_pp,
                        r_bb = r_bb)

good_par@species_params$ea_int[1]
good_par@ea_gro 
#bad_par_sp$ca_int <- -0.004 # Here we just use the fixed values
#bad_par_sp$ca_met <- 0.001 # Here we just use the fixed values

good <- project(good_par, 
                dt = dt,
                effort = projectEffort_m,
                temperature = projectTemp$temperature,
                diet_steps = 10,
                t_max = t_max,
                t_ref = 10)   

goodGrowth <- getGrowth(good)
goodGrowth$scen <- "Good"


#**** Plot all together ============================================================
df_all <- rbind(goodGrowth, badGrowth, refGrowth)

p1 <- df_all %>% filter(scen %in% c("Bad", "Ref")) %>% 
  ggplot(., aes(Age, value, color = scen, linetype = scen)) +
  geom_line(alpha = 0.7, size = 1.5) +
  facet_wrap(~Species, scales = "free") +
  theme_classic() +
  annotate("text", -Inf, Inf, 
           label = paste("met=", round(bad@params@species_params$ea_met[1], digits = 3),
                         "mor=", round(bad@params@species_params$ea_mor[1], digits = 3),
                         "int=", round(bad@params@species_params$ea_int[1], digits = 3),
                         "gro=", round(bad@params@ea_gro, digits = 3),
                         "car=", round(bad@params@ea_car, digits = 3),
                         sep = "\n"), 
           size = 4, fontface = "bold", hjust = -0.5, vjust = 1.3) +
  ggtitle("Bad") +
  NULL

p1

p2 <- df_all %>% filter(scen %in% c("Good", "Ref")) %>% 
  ggplot(., aes(Age, value, color = scen, linetype = scen)) +
  geom_line(alpha = 0.7, size = 1.5) +
  facet_wrap(~Species, scales = "free") +
  theme_classic() +
  annotate("text", -Inf, Inf, 
           label = paste("met=", round(good@params@species_params$ea_met[1], digits = 3),
                         "mor=", round(good@params@species_params$ea_mor[1], digits = 3),
                         "int=", round(good@params@species_params$ea_int[1], digits = 3),
                         "gro=", round(good@params@ea_gro, digits = 3),
                         "car=", round(good@params@ea_car, digits = 3),
                         sep = "\n"), 
           size = 4, fontface = "bold", hjust = -0.5, vjust = 1.3) +
  ggtitle("Good") +
  NULL

p2

p1/p2

# This doesn't look at all like we would predict given the temp-scalars. 
# Below Im plotting them:

#**** Plot temp scalars ============================================================
# Good
metScal <- as.matrix(good@metTempScalar[1,1,])
metScal_df <- data.frame(scal = metScal[, 1],
                         temp = as.numeric(as.character(rownames(metScal))),
                         param = "met",
                         scen = "good")

morScal <- as.matrix(good@morTempScalar[1,1,])
morScal_df <- data.frame(scal = morScal[, 1],
                         temp = as.numeric(as.character(rownames(morScal))),
                         param = "mor",
                         scen = "good")

intScal <- as.matrix(good@intTempScalar[1,1,])
intScal_df <- data.frame(scal = intScal[, 1],
                         temp = as.numeric(as.character(rownames(intScal))),
                         param = "int",
                         scen = "good")

carScal <- as.matrix(good@carTempScalar[2,])
carScal_df <- data.frame(scal = carScal[, 1],
                         temp = as.numeric(as.character(rownames(carScal))),
                         param = "car",
                         scen = "good")

groScal <- as.matrix(good@groTempScalar[2,])
groScal_df <- data.frame(scal = groScal[, 1],
                         temp = as.numeric(as.character(rownames(groScal))),
                         param = "gro",
                         scen = "good")

good_scalars <- rbind(metScal_df, morScal_df, intScal_df, carScal_df, groScal_df)

# Bad
metScal <- as.matrix(bad@metTempScalar[1,1,])
metScal_df <- data.frame(scal = metScal[, 1],
                         temp = as.numeric(as.character(rownames(metScal))),
                         param = "met",
                         scen = "bad")

morScal <- as.matrix(bad@morTempScalar[1,1,])
morScal_df <- data.frame(scal = morScal[, 1],
                         temp = as.numeric(as.character(rownames(morScal))),
                         param = "mor",
                         scen = "bad")

intScal <- as.matrix(bad@intTempScalar[1,1,])
intScal_df <- data.frame(scal = intScal[, 1],
                         temp = as.numeric(as.character(rownames(intScal))),
                         param = "int",
                         scen = "bad")

carScal <- as.matrix(bad@carTempScalar[2,])
carScal_df <- data.frame(scal = carScal[, 1],
                         temp = as.numeric(as.character(rownames(carScal))),
                         param = "car",
                         scen = "bad")

groScal <- as.matrix(bad@groTempScalar[2,])
groScal_df <- data.frame(scal = groScal[, 1],
                         temp = as.numeric(as.character(rownames(groScal))),
                         param = "gro",
                         scen = "bad")

bad_scalars <- rbind(metScal_df, morScal_df, intScal_df, carScal_df, groScal_df)

## All together
all_scalars <- rbind(good_scalars, bad_scalars)

ggplot(all_scalars, aes(temp, scal, color = scen)) + 
  facet_wrap(~ param) +
  geom_line(size = 1.5, alpha = 0.8) +
  theme_classic()


#**** Here I'm contrasting constant warm vs increasing temp, varying on param at the time. ====
# Good growth
test_par_sp2 <- bad_par_sp

test_par_sp2$ea_met <- 0.63
test_par_sp2$ea_mor <- 0.63
test_par_sp2$ea_int <- 0.9

test_par_sp2

test_par2 <- MizerParams(test_par_sp2,
                         ea_gro = 0,
                         ea_car = 0,
                         kappa_ben = kappa_ben,
                         kappa = kappa,
                         w_bb_cutoff = w_bb_cutoff,
                         w_pp_cutoff = w_pp_cutoff,
                         r_pp = r_pp,
                         r_bb = r_bb)

test2 <- project(test_par2, 
                 dt = dt,
                 effort = projectEffort_m,
                 temperature = rep(13, nrow(projectEffort_m)),
                 diet_steps = 10,
                 t_max = length(projectEffort_m),
                 t_ref = 10)   
 
testGrowth2 <- getGrowth(test2)
testGrowth2$scen <- "good"

# Bad growth
test_par_sp3 <- bad_par_sp

test_par_sp3$ea_met <- 0.63
test_par_sp3$ea_mor <- 0.63
test_par_sp3$ea_int <- 0.2

test_par_sp3

test_par3 <- MizerParams(test_par_sp3,
                         ea_gro = 0,
                         ea_car = 0,
                         kappa_ben = kappa_ben,
                         kappa = kappa,
                         w_bb_cutoff = w_bb_cutoff,
                         w_pp_cutoff = w_pp_cutoff,
                         r_pp = r_pp,
                         r_bb = r_bb)

test3 <- project(test_par3, 
                 dt = dt,
                 effort = projectEffort_m,
                 temperature = rep(13, nrow(projectEffort_m)),
                 diet_steps = 10,
                 t_max = length(projectEffort_m),
                 t_ref = 10)   

testGrowth3 <- getGrowth(test3)
testGrowth3$scen <- "bad"

df_all3 <- rbind(testGrowth2, testGrowth3, refGrowth)

ggplot(df_all3, aes(Age, value, color = scen)) + 
  facet_wrap(~ Species, scales = "free") +
  geom_line(size = 1.5, alpha = 0.8) +
  theme_classic()

# Ok, this looks correct! How about the true parameters and increasing temperature? 
#**** Project with temp and bad activation energies increasing warm =====================
bad5 <- project(bad_par, 
                dt = dt,
                effort = projectEffort_m,
                temperature = seq(from = 9, to = 11, length.out = nrow(projectEffort_m)),
                diet_steps = 10,
                t_max = t_max,
                t_ref = 10)   

bad5Growth <- getGrowth(bad5)
bad5Growth$scen <- "Bad"

#**** Project with temp and good activation energies incresing warm ====================
good5 <- project(good_par, 
                 dt = dt,
                 effort = projectEffort_m,
                 temperature = seq(from = 9, to = 11, length.out = nrow(projectEffort_m)),
                 diet_steps = 10,
                 t_max = t_max,
                 t_ref = 10)   

good5Growth <- getGrowth(good5)
good5Growth$scen <- "Good"

# Plot all
df_all5 <- rbind(bad5Growth, good5Growth, refGrowth)

p5 <- ggplot(df_all5, aes(Age, value, color = scen)) + 
  facet_wrap(~ Species, scales = "free") +
  geom_line(size = 1.5, alpha = 0.8) +
  theme_classic() +
  ggtitle("Increasing temp")


#**** Go back and plot the original runs again with empirical temperature ====================
df_all6 <- rbind(badGrowth, goodGrowth, refGrowth)

p6 <- ggplot(df_all6, aes(Age, value, color = scen)) + 
  facet_wrap(~ Species, scales = "free") +
  geom_line(size = 1.5, alpha = 0.8) +
  theme_classic() +
  ggtitle("Empirical temp")

p4/p5/p6


#**** How much extra years must I add to get steady state growth? ====================
str(projectEffort_m)
projectEffort_m2 <- data.frame(projectEffort_m)
app <- tail(projectEffort_m2, 1)
#projectEffort_m2 <- rbind(projectEffort_m2, rep(app, 500))
#projectEffort_m2 <- data.frame(rbind(as.matrix(projectEffort_m2), as.matrix(rep(app, 500))))
projectEffort_m2 <- projectEffort_m2 %>% add_row(Cod = rep(projectEffort_m2[137, 1], 500), 
                             Herring = rep(projectEffort_m2[137, 1], 500),
                             Sprat = rep(projectEffort_m2[137, 1], 500))

str(projectEffort_m2)
rownames(projectEffort_m2) <- 1:(nrow(projectEffort_m2))
projectEffort_m2 <- as.matrix(projectEffort_m2)
str(projectEffort_m2)
str(projectEffort_m)
projectEffort_m2
projectEffort_m
time_effort <- as.numeric(dimnames(projectEffort_m2)[[1]])
t_max <- time_effort[length(time_effort)]
t_max

projectTemp$temperature
#temperature2 <- c(projectTemp$temperature, rep(projectTemp$temperature[137], 50))
#temperature2 <- c(projectTemp$temperature, rep(11, 50))
#temperature2 <- rep(11, t_max)
temperature2 <- c(projectTemp$temperature, rep(projectTemp$temperature[137], 500))
str(temperature2)
temperature2
projectTemp$temperature

length(temperature2)
t_max

#**** Project with temp and bad activation energies increasing warm ================
bad7 <- project(bad_par, 
                dt = dt,
                effort = projectEffort_m2,
                temperature = temperature2,
                diet_steps = 10,
                t_max = t_max,
                t_ref = 10)   

bad7Growth <- getGrowth(bad7)
bad7Growth$scen <- "Bad"

#**** Project with temp and good activation energies incresing warm ====================
good7 <- project(good_par, 
                 dt = dt,
                 effort = projectEffort_m2,
                 temperature = temperature2,
                 diet_steps = 10,
                 t_max = t_max,
                 t_ref = 10)   

good7Growth <- getGrowth(good7)
good7Growth$scen <- "Good"

# Plot all
df_all7 <- rbind(bad7Growth, good7Growth, refGrowth)

ggplot(df_all7, aes(Age, value, color = scen)) + 
  facet_wrap(~ Species, scales = "free") +
  geom_line(size = 1.5, alpha = 0.8) +
  theme_classic() +
  ggtitle("Increasing temp")

# Still doesnt work with 500 additional years...


#**** How much can I perturb the constant temperature and still get predicted results? ====
projectTemp$temperature
temperature2 <- projectTemp$temperature
#temperature2[60:137] <- 11
#temperature2[20:137] <- 11
#temperature2[1:137] <- 11
#temperature2[2:137] <- 11
#temperature2[1] <- 11
#temperature2[2] <- 99999999999999
temperature2

#**** Project with temp and bad activation energies increasing warm ================
bad8 <- project(bad_par, 
                dt = dt,
                effort = projectEffort_m,
                temperature = temperature2,
                diet_steps = 10,
                t_max = t_max,
                t_ref = 10)   

bad8Growth <- getGrowth(bad8)
bad8Growth$scen <- "Bad"

#**** Project with temp and good activation energies incresing warm ====================
good8 <- project(good_par, 
                 dt = dt,
                 effort = projectEffort_m,
                 temperature = temperature2,
                 diet_steps = 10,
                 t_max = t_max,
                 t_ref = 10)   

good8Growth <- getGrowth(good8)
good8Growth$scen <- "Good"

# Plot all
df_all8 <- rbind(bad8Growth, good8Growth, refGrowth)

ggplot(df_all8, aes(Age, value, color = scen)) + 
  facet_wrap(~ Species, scales = "free") +
  geom_line(size = 1.5, alpha = 0.8) +
  theme_classic() +
  ggtitle("Increasing temp")

## THIS IS EXTREMELY STRANGE!!! IT SEEMS ONLY THE FIRST TEMPERATURE MATTERS
## So, since that is lower than t_ref, then the results make sense
## When the first value of the temperature vector is higher than t-ref,
## the results make sense even with empirical temperatures
## The temperature is though implemented correctly
## So it must be the growth function that is bugged
## First try to plot the built in growth function with empirical temperature
## to see if my getGrowth is bugged or working as intended but not the way
## I want it to.

p7 <- plotGrowthCurves(good8) + facet_wrap(~ Species, scales = "free")
p8 <- plotGrowthCurves(bad8) + facet_wrap(~ Species, scales = "free")
p7/p8 # here also bad is better than good, so theres something "wrong" with the growth func
plotBiomass(good8) # try this with the extrme temperature


# Is the problem that getGrowth and plotGrowthCurves do not take the length at age
# from the last iteration?
plotGrowthCurves

# g is created by the getEGrowth function, which uses:
# sim@metTempScalar[,,1]. What is the 1 here?
str(good8@metTempScalar)
good8@metTempScalar[,,nrow(projectEffort_m)]
good8@metTempScalar[,,length(time_temperature_dt)]

# it's 5*137, because (from project.R):
time_effort <- as.numeric(dimnames(projectEffort_m)[[1]])
t_max <- time_effort[length(time_effort)]
time_temperature_dt <- rep(projectTemp$temperature, length = t_max/dt, each = 1/dt) # works if t_max = length(temperature)
x_axis <- seq(length.out=(t_max/dt),from =1) 
temperature_dt <- matrix(time_temperature_dt, dimnames = list(x_axis, "temperature")) # without smoothing
temperature_dt
TESTmetTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
str(TESTmetTempScalar)

# Ok, now replacing 1 in sim@metTempScalar[,,1] with length(time_temperature_dt)
# in the getGrowth function... 

## OK, now I re-source the growth function...


## Does this also affect plotSpectra?
plotSpectra
# Does not seem like it!

#**** Random notes =================================================================
# Temperature needs to be length as t_max
# if (length(temperature) != t_max) {
#   stop("your temperature input vector is not the same length as t_max")

projectEffort_m <- as.matrix(projectEffort)
effort = projectEffort_m
temperature = projectTemp$temperature

tes <- data.frame(temp = projectTemp$temperature,
                  eff = effort[, 1],
                  yr = 1:length(temperature))


plot(tes$yr, tes$temp)  
plot(tes$yr, tes$eff)  

length(projectEffort_m)
length(temperature)

# here the lengths differ... how does that affect t_max

# line 69 in "Project.R"
#' If effort is specified as an array then the smallest time in the array is 
#' used as the initial time for the simulation. Otherwise the initial time is
#' set to 0. Also, if the effort is an array then the \code{t_max} argument is 
#' ignored and the maximum simulation time is the largest time of the effort
#' array.

# Ok, so what is t_max now that we have an effort array?
time_temperature_dt <- rep(temperature, length = t_max/dt, each = 1/dt) # works if t_max = length(temperature)
time_effort <- as.numeric(dimnames(projectEffort)[[1]])
t_max <- time_effort[length(time_effort)]

# No, it doesn't, because time_effort is the length of the first object in the list... so it matches temperature
# Can I plot temperature and effort in the data and in the sim-output?

# Here temperature is a list with t_max*number of temperature series (rates). This can be seen
# in MizerSimClass
str(good4@temperature)[1:137]
tempz <- good4@temperature[1:137]
length(tempz)
plot(tempz)

# Now try with original object. Here I rerun the first analysis to try and understand how 
# temperature is stored
bad_par_sp <- params@species_params

bad_par_sp$ea_met <- max(ea$met)
bad_par_sp$ea_int <- min(ea$int)
bad_par_sp$ea_mor <- max(ea$mor)

bad_par <- MizerParams(bad_par_sp,
                       ea_gro = min(ea$gro),
                       ea_car = min(ea$car),
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb)

#bad_par_sp$ca_int <- -0.004 # Here we just use the fixed values
#bad_par_sp$ca_met <- 0.001 # Here we just use the fixed values

bad <- project(bad_par, 
               dt = dt,
               effort = projectEffort_m,
               temperature = projectTemp$temperature,
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10)   

badGrowth <- getGrowth(bad)
badGrowth$scen <- "Bad"

str(bad@temperature)[1:137]
tempz <- bad@temperature[1:137]
tempz <- bad@temperature
tempz
length(tempz)
plot(tempz)
