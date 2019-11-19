#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. FMSY at two different temperatures
# 
# C. (HASHTAGGED) Analyisis of changes in yield with fishing warming
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
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getSpectra.R", ssl.verifypeer = FALSE)
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
t_ref <- (10 + 0.1156161)
kappa_ben <- 1
kappa <- 1
w_bb_cutoff <- 20
w_pp_cutoff <- 1
r_pp <- 4
r_bb <- 4

#**** Update species params ========================================================
t <- params@species_params

t$ea_met <- mean(ea$met)
t$ea_int <- mean(ea$int)
t$ea_mor <- mean(ea$mor)

t$ca_int <- -0.004 # Here we just use the fixed values
t$ca_met <- 0.001 # Here we just use the fixed values

pars_no_res <- MizerParams(t, 
                           ea_gro = 0,
                           ea_car = 0, # -ea$gro[i] 
                           kappa_ben = kappa_ben,
                           kappa = kappa,
                           w_bb_cutoff = w_bb_cutoff,
                           w_pp_cutoff = w_pp_cutoff,
                           r_pp = r_pp,
                           r_bb = r_bb)

pars_with_res <- MizerParams(t, 
                             ea_gro = mean(ea$gro),
                             ea_car = mean(ea$car), # -ea$gro[i] 
                             kappa_ben = kappa_ben,
                             kappa = kappa,
                             w_bb_cutoff = w_bb_cutoff,
                             w_pp_cutoff = w_pp_cutoff,
                             r_pp = r_pp,
                             r_bb = r_bb)


# B. EXTRACT FMSY FROM MODEL AT DIFFERENT TEMPERATURES =============================
# In lack of a better approach, I will just for-loop different F, extract Yield, plot
# over F. I will increase each species F separately, keeping the others at their mean

F_range <- seq(0.2, 1.4, 0.2) # Can decrease step later, becomes too slow now
t_max <- 100
index <- 1:length(F_range)

#**** Cod ==========================================================================
# Create effort vector
effort = c(Cod = params@species_params$AveEffort[1], 
           Herring = params@species_params$AveEffort[3], 
           Sprat = params@species_params$AveEffort[2])

# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  effort[1] <- F_range[i]
  
  r <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep(t_ref, t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)
  
  w <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep((t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)
  
  Y_ref <- mean(data.frame(getYield(r))$Cod[(t_max-20):t_max])
  Y_warm <- mean(data.frame(getYield(w))$Cod[(t_max-20):t_max])
  
  ssb_ref <- mean(data.frame(getSSB(r))$Cod[(t_max-20):t_max])
  ssb_warm <- mean(data.frame(getSSB(w))$Cod[(t_max-20):t_max])
  
  Fm <- F_range[i]
  
  td <- data.frame(biomass = rbind(Y_ref, Y_warm, ssb_ref, ssb_warm), 
                   type = rep(c("yield", "ssb"), each = 2), 
                   Fm = rep(Fm, 4), 
                   scen = rep(c("cold", "warm"), times = 2))
  
  data_list[[i]] <- td
  
}

codFmsy <- dplyr::bind_rows(data_list)

codFmsy$species <- "Cod"

ggplot(codFmsy, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


#**** Herring ======================================================================
# Create effort vector
effort = c(Cod = params@species_params$AveEffort[1], 
           Herring = params@species_params$AveEffort[3], 
           Sprat = params@species_params$AveEffort[2])

# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  effort[2] <- F_range[i]
  
  r <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep(t_ref, t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)
  
  w <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep((t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)
  
  Y_ref <- mean(data.frame(getYield(r))$Herring[(t_max-20):t_max])
  Y_warm <- mean(data.frame(getYield(w))$Herring[(t_max-20):t_max])
  
  ssb_ref <- mean(data.frame(getSSB(r))$Herring[(t_max-20):t_max])
  ssb_warm <- mean(data.frame(getSSB(w))$Herring[(t_max-20):t_max])
  
  Fm <- F_range[i]
  
  td <- data.frame(biomass = rbind(Y_ref, Y_warm, ssb_ref, ssb_warm), 
                   type = rep(c("yield", "ssb"), each = 2), 
                   Fm = rep(Fm, 4), 
                   scen = rep(c("cold", "warm"), times = 2))  
  
  data_list[[i]] <- td
  
}

herFmsy <- dplyr::bind_rows(data_list)

herFmsy$species <- "Herring"

ggplot(herFmsy, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


#**** Sprat ========================================================================
# Create effort vector
effort = c(Cod = params@species_params$AveEffort[1], 
           Herring = params@species_params$AveEffort[3], 
           Sprat = params@species_params$AveEffort[2])

# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  effort[3] <- F_range[i]
  
  r <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep(t_ref, t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)
  
  w <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep((t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)
  
  Y_ref <- mean(data.frame(getYield(r))$Sprat[(t_max-20):t_max])
  Y_warm <- mean(data.frame(getYield(w))$Sprat[(t_max-20):t_max])
  
  ssb_ref <- mean(data.frame(getSSB(r))$Sprat[(t_max-20):t_max])
  ssb_warm <- mean(data.frame(getSSB(w))$Sprat[(t_max-20):t_max])
  
  Fm <- F_range[i]
  
  td <- data.frame(biomass = rbind(Y_ref, Y_warm, ssb_ref, ssb_warm), 
                   type = rep(c("yield", "ssb"), each = 2), 
                   Fm = rep(Fm, 4), 
                   scen = rep(c("cold", "warm"), times = 2))  
  
  data_list[[i]] <- td
  
}

sprFmsy <- dplyr::bind_rows(data_list)

sprFmsy$species <- "sprat"

ggplot(sprFmsy, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


#**** All together =================================================================
col <- RColorBrewer::brewer.pal(n = 5, "Dark2")

Fmsy <- rbind(codFmsy, sprFmsy, herFmsy)

Fmsy %>% filter(biomass > 0.001) %>% 
  ggplot(., aes(Fm, biomass, linetype = type, color = scen)) + 
  geom_line(alpha = 0.6, size = 1.2) +
  facet_wrap(~ species, scales = "free") +
  scale_color_manual(values = rev(col)) +
  theme_classic(base_size = 10) +
  labs(x = "Fishing mortality [1/year]", y = "Biomass",
       color = "Scenario",
       linetype = "Metric") +
  theme(aspect.ratio = 3/4,
        legend.position = "bottom") +
  NULL

#ggsave("baltic/figures/FMSY_warm_cold.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


#**** Testing if getSSB results depebd on the fishing mortality ===================================================
test <- filter(Fmsy, Fm == 1)

ggplot(test, aes(Fm, biomass, color = scen)) +
  facet_wrap(species ~ type, scales = "free") +
  geom_point()

# Ok, here cold has higher ssb and yield, but in size_spectra.R, it's the opposite.
# The dfífferences are that here temperature is held constant and effort is average, i.e. not appended to a time series
# should that really affect warm/cold prediction?
projectEffort <- read.csv("baltic/params/projectEffort.csv")[, 2:4]
projectTemp <- read.csv("baltic/params/projectTemp.csv")

projectEffort_m <- as.matrix(projectEffort)
rownames(projectEffort_m) <- 1:nrow(projectEffort)

# Here params are the same but effort is not average effort but fmsy, and temperature is empirical
test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = projectTemp$temperature,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

consTemp <- projectTemp$temperature
consTemp[] <- t_ref
test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# Ok, so here warming leads to more ssb... why is not this captured in the FMSY code based on not time-varying effort and temp?
test_cold@params@species_params$AveEffort
tail(projectEffort_m, 1)

# Test to use the above script but with average effort appended instead of median FMSY from assessment.
# Then I can say if it's the effort itself, but I doubt it.
projectEffort <- read.csv("baltic/params/projectEffort.csv")[, 2:4]
projectTemp <- read.csv("baltic/params/projectTemp.csv")

projectEffort_m <- as.matrix(projectEffort)
rownames(projectEffort_m) <- 1:nrow(projectEffort)

projectEffort_new <- projectEffort_m
projectEffort_new[100:nrow(projectEffort_m), ][, 1] <- test_cold@params@species_params$AveEffort[1]
projectEffort_new[100:nrow(projectEffort_m), ][, 2] <- test_cold@params@species_params$AveEffort[3]
projectEffort_new[100:nrow(projectEffort_m), ][, 3] <- test_cold@params@species_params$AveEffort[2]

projectEffort_new

# Here I use the average effort appended, to see if it's the effort it self and not the time-varying part that leads to different results
test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_new,
                     temperature = projectTemp$temperature,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

consTemp <- projectTemp$temperature
consTemp[] <- t_ref
test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_new,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# Again here warming has a higher yield.
# Then it's only temperature that varies. Next try constant temperature and FMSY effort (as in spectra) to see if I can reproduce the current yield plot.
test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = (consTemp+2),
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

consTemp <- projectTemp$temperature
consTemp[] <- t_ref
test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# Aha! Here I get lower yields when temperature is constant but warmer (and I do get that with the average effort as well), so it seems to be an effect of that. 
# Why is that? Below I'm testing the magnitude of warming by doing two constant temperatures with a difference equal to empirical - t_ref.
delta <- projectTemp$temperature[137] - t_ref

test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = (consTemp+delta),
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

consTemp <- projectTemp$temperature
consTemp[] <- t_ref
test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# Here I get lower yield in the warming scenario still, as I always get when temperature is constant.
# Next I'm trying to extend to continuation by repeating the last temperature 50 times.
temperature2 <- projectTemp$temperature
temperature2[100:137] <- projectTemp$temperature[137] 
temperature2

test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = temperature2,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

consTemp <- projectTemp$temperature
consTemp[] <- t_ref
test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# Here yield is higher when it's warmer. NOT like in the constant scenario and NOT what we would predict from growth rates and size spectra BUT
# in line with the time-varying temperature-scenarios. 

# I'll now instead append values rather than change the last
t <- as.data.frame(projectEffort_m)[100:137, ]
projectEffort_d2 <- rbind(as.data.frame(projectEffort_m), t, t, t)
projectEffort_m2 <- as.matrix(projectEffort_d2)
str(projectEffort_m2)
rownames(projectEffort_m2) <- 1:nrow(projectEffort_d2)
projectEffort_m2

# Now do temperature
temp3 <- projectTemp$temperature
temp4 <- c(temp3, rep(projectTemp$temperature[137], nrow(projectEffort_m2)-137))
str(temp4)
temp4

test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m2,
                     temperature = temp4,
                     diet_steps = 10,
                     #t_max = t_max,
                     t_ref = t_ref)  

consTemp <- temp4
consTemp[] <- t_ref
test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m2,
                     temperature = consTemp,
                     diet_steps = 10,
                     #t_max = t_max,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# Here yield is higher with warming... Why can't I reproduce the effects of constant
# temperatures with time varying temperatures that match in magnitude?

# Basically, why are these different:
# Here I use the same time-varying effort, and two temperature vecors, but the time varying one approaches the constant temp

consTemp <- temp4
consTemp[] <- temp4[175]
plot(temp4)
lines(consTemp)

test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m2,
                     temperature = temp4,
                     diet_steps = 10,
                     t_ref = t_ref)  

test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m2,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 5)
tail(getSSB(test_cold), 5)

# Aha! Both are in equilbirum and both have EXACTLY the same parameters in the end.
plot(temp4)
lines(consTemp)

# That means that the initial values lead to different steady states and therefore 
# different predictions about the effect of temperature. 
# The scenario with the initially time-varying temperature has a lower temp, meaning
# higer carrying capacity and better scaling, which could affect the later iterations
# as well.

# What do we do now? Comparing two models with different starting values does not sound
# good. And the effect of temperature must match that of a non-time varying temperature
# treatment. 

# So I will now try and keep the time-varying temperature until maybe 1997, then let the 
# constant be the average of the calibration period and the warming be the increasing.

# What I expect to see now is that the ssb is declining with temperature, as that is what
# the constant temperature scenarios show. That is not what the other results (size spectra,
# growth) indicate, but that may be because I need to change those as well accordingly.

# Now do temperature
consTemp <- projectTemp$temperature
#consTemp[100:137] <- mean(consTemp[78:88]) # This corresponds to 1992-2002 in years rather than rownumbers
consTemp[100:137] <- mean(consTemp[100]) # This is 2014... to not get decrease in temp after calibration

plot(projectTemp$temperature)
points(consTemp, col = "red", cex = 0.5) 

test_warm <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = projectTemp$temperature,
                     diet_steps = 10,
                     t_ref = t_ref)  

test_cold <- project(pars_with_res, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_ref = t_ref)  

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# Still slightly more SLIGHTLY yield with warming. Maybe that's an effect of resource?
test_warm2 <- project(pars_no_res, 
                      dt = dt,
                      effort = projectEffort_m,
                      temperature = projectTemp$temperature,
                      diet_steps = 10,
                      t_ref = t_ref)  

test_cold2 <- project(pars_no_res, 
                      dt = dt,
                      effort = projectEffort_m,
                      temperature = consTemp,
                      diet_steps = 10,
                      t_ref = t_ref)  

tail(getSSB(test_warm2), 2)
tail(getSSB(test_cold2), 2)

# Still slightly more SLIGHTLY yield with warming. Maybe that's an effect of resource?
tail(getSSB(test_warm), 1) / tail(getSSB(test_cold), 1)
tail(getSSB(test_warm2), 1) / tail(getSSB(test_cold2), 1)

# The relative difference is the same because temperature is the same. if I increase temperature more, 
# even with the same factor, then the ration of warm/cold starts to decrease again 
temp5 <- projectTemp$temperature
temp5[100:137] <- projectTemp$temperature[100:137] * 1.05

plot(projectTemp$temperature, ylim = c(9, 13))
lines(consTemp, col = "blue")
lines(temp5, col = "red")

test_warm3 <- project(pars_no_res, 
                      dt = dt,
                      effort = projectEffort_m,
                      temperature = temp5,
                      diet_steps = 10,
                      t_ref = t_ref)  

test_cold3 <- project(pars_no_res, 
                      dt = dt,
                      effort = projectEffort_m,
                      temperature = consTemp,
                      diet_steps = 10,
                      t_ref = t_ref)  

tail(getSSB(test_warm3), 2)
tail(getSSB(test_cold3), 2)

# 
tail(getSSB(test_warm), 1) / tail(getSSB(test_cold), 1)
tail(getSSB(test_warm2), 1) / tail(getSSB(test_cold2), 1)
tail(getSSB(test_warm3), 1) / tail(getSSB(test_cold3), 1)

# The third scenario has a higher temperature increase in the time after the divergence.
# the first and second row is the same because the parameters change in the same way.
# Remember because the temperature is no longer constant in the "constant" scenario,
# the parameters matter.
# This leads to a higher ratio, so the effect is in fact linear.

# I think the STILL positive effect even after controlling for initial values is because
# we never consider no increase in growth but decrease in carrying capacity.
# TEST that using both constant and time varying effort and temperature.
t <- params@species_params
t$ea_met <- 0.8 #mean(ea$met)
t$ea_int <- 0.4 #mean(ea$int)
t$ea_mor <- 0.8 #mean(ea$mor)
t$ca_int <- -0.004 # Here we just use the fixed values
t$ca_met <- 0.001 # Here we just use the fixed values

pars_with_res2 <- MizerParams(t, 
                              ea_gro = 0, 
                              ea_car = mean(ea$car), 
                              kappa_ben = kappa_ben,
                              kappa = kappa,
                              w_bb_cutoff = w_bb_cutoff,
                              w_pp_cutoff = w_pp_cutoff,
                              r_pp = r_pp,
                              r_bb = r_bb,
                              t_ref = consTemp[137])

test_warm <- project(pars_with_res2, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = projectTemp$temperature,
                     diet_steps = 10)  

test_cold <- project(pars_with_res2, 
                     dt = dt,
                     effort = projectEffort_m,
                     temperature = consTemp,
                     diet_steps = 10,
                     t_ref = t_ref)  

plot(test_warm@carTempScalar[1, ])
points(test_cold@carTempScalar[1, ], col = "red", cex = 0.2)

tail(getSSB(test_warm), 2)
tail(getSSB(test_cold), 2)

# > tail(getSSB(test_warm), 2)
# sp
# time       Cod   Sprat  Herring
# 136 1.120812 4.55859 1.627390
# 137 1.120439 4.54112 1.627152
# > tail(getSSB(test_cold), 2)
# sp
# time       Cod    Sprat  Herring
# 136 1.102057 4.466224 1.596122
# 137 1.102108 4.468228 1.594283



# C. TEMP-DRIVEN CHANGE IN FISHERIES YIELD =========================================
# for-loop to take random samples for distributions representing activation energies
# Then compare that to a projection with a constant temperature

#**** Project reference scenario ===================================================
# consTemp <- projectTemp$temperature
# t_ref <- (10 + 0.1156161)
# consTemp[] <- t_ref
# 
# ref <- project(pars_with_res, 
#                dt = dt,
#                effort = projectEffort_m,
#                temperature = consTemp,
#                diet_steps = 10,
#                t_max = t_max,
#                t_ref = (10 + 0.1156161))   
# 
# refYield <- data.frame(Yield = c(getYield(ref)[dim(ref@effort)[1], 1],
#                                  getYield(ref)[dim(ref@effort)[1], 2],
#                                  getYield(ref)[dim(ref@effort)[1], 3]),
#                        Species = ref@params@species_params$species)
# 

# #**** for loop through different fishin effort (with temp dep resource) ============
# eff <- seq(0.8, 1.2, 0.02) # Factor for scaling fishing mortality
# temp <- seq(0.8, 1.2, 0.02)
# temp_eff <- expand.grid(data.frame(eff = eff, temp = temp))
# iter <- seq(from = 1, to = nrow(temp_eff))
# 
# projectEffort_new <- projectEffort_m
# projectTemp_new <- projectTemp$temperature
# tt <- c()
# groj <- c()
# yield <- c()
# data_list <- list()
# 
# # The projected fishing mortality starts at row 100
# 
# for (i in iter) {
#   
#   projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * temp_eff$eff[i]
#   
#   projectTemp_new[100:length(projectTemp$temperature)] <- projectTemp$temperature[100:length(projectTemp$temperature)] * temp_eff$temp[i]
#   
#   proj <- project(pars_with_res, 
#                   dt = dt,
#                   effort = projectEffort_new,
#                   temperature = projectTemp_new,
#                   diet_steps = 10,
#                   t_max = t_max,
#                   t_ref = t_ref)   
#   
#   # Extract yield at last iteration
#   proYield <- data.frame(Yield = c(getYield(proj)[dim(proj@effort)[1], 1],
#                                    getYield(proj)[dim(proj@effort)[1], 2],
#                                    getYield(proj)[dim(proj@effort)[1], 3]),
#                          
#                          Species = proj@params@species_params$species,
#                          Fm = proj@effort[dim(proj@effort)[1], ],
#                          temp = proj@temperature[dim(proj@temperature)[1], 1],
#                          Fm_scal = temp_eff$eff[i],
#                          temp_scal = temp_eff$temp[i])
#   
#   proYield$Yield_rel <- proYield$Yield / refYield$Yield
#   
#   data_list[[i]] <- proYield
#   
# }
# 
# str(data_list)
# 
# big_yield_data <- dplyr::bind_rows(data_list)
# 
# big_yield_data$Species <- factor(big_yield_data$Species, levels = c("Sprat", "Herring", "Cod"))
# 
# ggplot(big_yield_data, aes(temp_scal, Fm_scal, fill = Yield_rel)) + 
#   geom_tile(color = NA) +
#   facet_grid(~ Species, scales = "free") +
#   theme_classic(base_size = 10) +
#   scale_fill_viridis(option = "cividis") +
#   labs(x = "Temperature relative to\nRCP8.5 projection", 
#        y = "Fishing mortality relative\nto FMSY (assessment)",
#        fill = "Yield relative to\nFMSY +\nconstant temp.") +
#   coord_cartesian(expand = 0) +
#   theme(aspect.ratio = 3/4) +
#   NULL

#ggsave("baltic/figures/yield_heat.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")
