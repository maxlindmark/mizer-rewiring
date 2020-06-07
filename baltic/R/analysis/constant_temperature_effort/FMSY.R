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
# C. Analyisis of changes in yield with fishing warming (Heatmap)
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
#devtools::load_all(".") # THIS DOES NOT WORK ON NEW MAC; MAYBE NOT NEEDED ANYMORE SINCE
# I PUSHED CHANGES TO THE CODE ALREADY; NO NEED TO WORK IN LOCAL LIBRAY....

# Install the specific mizer version from github
# devtools::install_github("maxlindmark/mizer-rewiring", ref = "rewire-temp") 
library(mizer)

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
## TESTING TRAIT MODEL ## 
# params_trait <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5)
# simTrat <- project(params_trait, t_max=75, effort = 1)
# plot(simTrat)
# tail(getSSB(simTrat)) / tail(getYield(simTrat))
## END TEST, YIELD CAN BE LARGER THAN SSB... ## 

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
t_max <- 100 # For FMSY plot!


#**** Update species params ========================================================
t <- params@species_params

t$ea_met <- mean(ea$met)
t$ea_int <- mean(ea$int)
t$ea_mor <- mean(ea$mor)

#t$ca_int <- -0.004 # Here we just use the fixed values
#t$ca_met <- 0.001 # Here we just use the fixed values

pars_phys <- MizerParams(t, 
                         ea_gro = 0,
                         ea_car = 0, # -ea$gro[i] 
                         kappa_ben = kappa_ben,
                         kappa = kappa,
                         w_bb_cutoff = w_bb_cutoff,
                         w_pp_cutoff = w_pp_cutoff,
                         r_pp = r_pp,
                         r_bb = r_bb,
                         t_ref = t_ref)

pars_res_phys <- MizerParams(t, 
                             ea_gro = mean(ea$gro),
                             ea_car = mean(ea$car), # -ea$gro[i] 
                             kappa_ben = kappa_ben,
                             kappa = kappa,
                             w_bb_cutoff = w_bb_cutoff,
                             w_pp_cutoff = w_pp_cutoff,
                             r_pp = r_pp,
                             r_bb = r_bb,
                             t_ref = t_ref)

t_np <- params@species_params
t_np$ea_int <- 0
t_np$ea_met <- 0
t_np$ea_mat <- 0
t_np$ea_mor <- 0

pars_res <- MizerParams(t_np, 
                        ea_gro = mean(ea$gro),
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
start <- 1997-1874
consTemp[start:177] <- t_ref


# B. EXTRACT FMSY FROM MODEL AT DIFFERENT TEMPERATURES =============================
# Here I just for-loop different F, extract Yield, plotover F. 
# I will increase each species F separately, keeping the others at their mean, as
# in the calibration method.
# NOTE: Here we use constant temperatures and therefore the results may not be 
# exactly translateable to the time-varying effort and temperature projections.

F_range <- seq(0, 1.1, 0.01) # Can decrease step later, becomes too slow now
t_max <- 100
index <- 1:length(F_range)


#** Physio + Resource (exp) ========================================================
#**** Cod ==========================================================================
# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  effort = c(Cod = params@species_params$AveEffort[1], 
             Herring = params@species_params$AveEffort[3], 
             Sprat = params@species_params$AveEffort[2])
  
  effort[1] <- F_range[i]
  
  r <- project(pars_res_phys,
               dt = dt,
               effort = effort,
               temperature = rep(pars_res_phys@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_res_phys,
               dt = dt,
               effort = effort,
               temperature = rep((pars_res_phys@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

codFmsy_res_phys <- dplyr::bind_rows(data_list)

codFmsy_res_phys$species <- "Cod"

codFmsy_res_phys$scen2 <- "Physio. + Resource (exp.)"

ggplot(codFmsy_res_phys, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


#**** Herring ======================================================================
# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  # Create effort vector
  effort = c(Cod = params@species_params$AveEffort[1], 
             Herring = params@species_params$AveEffort[3], 
             Sprat = params@species_params$AveEffort[2])
  
  effort[2] <- F_range[i]
  
  r <- project(pars_res_phys,
               dt = dt,
               effort = effort,
               temperature = rep(pars_res_phys@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_res_phys,
               dt = dt,
               effort = effort,
               temperature = rep((pars_res_phys@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

herFmsy_res_phys <- dplyr::bind_rows(data_list)

herFmsy_res_phys$species <- "Herring"

herFmsy_res_phys$scen2 <- "Physio. + Resource (exp.)"

ggplot(herFmsy_res_phys, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


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
  
  r <- project(pars_res_phys,
               dt = dt,
               effort = effort,
               temperature = rep(pars_res_phys@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_res_phys,
               dt = dt,
               effort = effort,
               temperature = rep((pars_res_phys@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

sprFmsy_res_phys <- dplyr::bind_rows(data_list)

sprFmsy_res_phys$species <- "Sprat"

sprFmsy_res_phys$scen2 <- "Physio. + Resource (exp.)"

ggplot(sprFmsy_res_phys, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 



#** Resource (exp) ========================================================
#**** Cod ==========================================================================
# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  effort = c(Cod = params@species_params$AveEffort[1], 
             Herring = params@species_params$AveEffort[3], 
             Sprat = params@species_params$AveEffort[2])
  
  effort[1] <- F_range[i]
  
  r <- project(pars_res,
               dt = dt,
               effort = effort,
               temperature = rep(pars_res@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_res,
               dt = dt,
               effort = effort,
               temperature = rep((pars_res@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

codFmsy_res <- dplyr::bind_rows(data_list)

codFmsy_res$species <- "Cod"

codFmsy_res$scen2 <- "Resource (exp.)"

ggplot(codFmsy_res, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


#**** Herring ======================================================================
# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  # Create effort vector
  effort = c(Cod = params@species_params$AveEffort[1], 
             Herring = params@species_params$AveEffort[3], 
             Sprat = params@species_params$AveEffort[2])
  
  effort[2] <- F_range[i]
  
  r <- project(pars_res,
               dt = dt,
               effort = effort,
               temperature = rep(pars_res@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_res,
               dt = dt,
               effort = effort,
               temperature = rep((pars_res@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

herFmsy_res <- dplyr::bind_rows(data_list)

herFmsy_res$species <- "Herring"

herFmsy_res$scen2 <- "Resource (exp.)"

ggplot(herFmsy_res, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


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
  
  r <- project(pars_res,
               dt = dt,
               effort = effort,
               temperature = rep(pars_res@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_res,
               dt = dt,
               effort = effort,
               temperature = rep((pars_res@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

sprFmsy_res <- dplyr::bind_rows(data_list)

sprFmsy_res$species <- "Sprat"

sprFmsy_res$scen2 <- "Resource (exp.)"

ggplot(sprFmsy_res, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 



#** Physio =========================================================================
#**** Cod ==========================================================================
# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  effort = c(Cod = params@species_params$AveEffort[1], 
             Herring = params@species_params$AveEffort[3], 
             Sprat = params@species_params$AveEffort[2])
  
  effort[1] <- F_range[i]
  
  r <- project(pars_phys,
               dt = dt,
               effort = effort,
               temperature = rep(pars_phys@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_phys,
               dt = dt,
               effort = effort,
               temperature = rep((pars_phys@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

codFmsy_phys <- dplyr::bind_rows(data_list)

codFmsy_phys$species <- "Cod"

codFmsy_phys$scen2 <- "Physio."

ggplot(codFmsy_phys, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


#**** Herring ======================================================================
# Create empty data holder
Y <- c()
Fm <- c()
r <- c()
w <- c()
td <- c()
data_list <- list()

for(i in index) {
  
  # Create effort vector
  effort = c(Cod = params@species_params$AveEffort[1], 
             Herring = params@species_params$AveEffort[3], 
             Sprat = params@species_params$AveEffort[2])
  
  effort[2] <- F_range[i]
  
  r <- project(pars_phys,
               dt = dt,
               effort = effort,
               temperature = rep(pars_phys@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_phys,
               dt = dt,
               effort = effort,
               temperature = rep((pars_phys@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

herFmsy_phys <- dplyr::bind_rows(data_list)

herFmsy_phys$species <- "Herring"

herFmsy_phys$scen2 <- "Physio."

ggplot(herFmsy_phys, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


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
  
  r <- project(pars_phys,
               dt = dt,
               effort = effort,
               temperature = rep(pars_phys@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_phys,
               dt = dt,
               effort = effort,
               temperature = rep((pars_phys@t_ref + 2), t_max),
               diet_steps = 10,
               t_max = t_max)
  
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

sprFmsy_phys <- dplyr::bind_rows(data_list)

sprFmsy_phys$species <- "Sprat"

sprFmsy_phys$scen2 <- "Physio."

ggplot(sprFmsy_phys, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 



#**** All together =================================================================
#col <- RColorBrewer::brewer.pal(n = 5, "Dark2")
col <- RColorBrewer::brewer.pal(n = 3, "Set1")[1:2]

Fmsy <- rbind(codFmsy_phys, sprFmsy_phys, herFmsy_phys,
              codFmsy_res_phys, sprFmsy_res_phys, herFmsy_res_phys,
              codFmsy_res, sprFmsy_res, herFmsy_res)

Fmsy$species <- factor(Fmsy$species, levels = c("Sprat", "Herring", "Cod"))

# Get size-spectrum FMSY for all species:
Fmsy_sum <- Fmsy %>% 
  filter(type == "yield") %>% 
  group_by(species, scen, scen2) %>% 
  filter(biomass == max(biomass))

p1 <- Fmsy %>% 
  filter(type == "yield") %>% 
  filter(biomass > 0.001) %>% 
  ggplot(., aes(Fm, (biomass*240.342), linetype = scen2, color = scen)) + 
  geom_line(alpha = 0.6, size = 1.2) +
  facet_wrap(~ species, scales = "free") +
  scale_color_manual(values = c(col[2], col[1], col[3]), 
                     labels = c(expression("T"[ref]), 
                                expression(paste("T"[ref], "+2", degree*C)))) +
  labs(x = "Fishing mortality [1/year]", 
       y = "Yield [1000 tonnes/year]",
       color = "Scenario",
       linetype = "Metric") +
  geom_segment(data = filter(Fmsy_sum, scen == "warm" & scen2 == "Physio. + Resource (exp.)"), linetype = 3, 
               aes(x = Fm, xend = Fm, y = c(0, 0, 0), yend = c(120, 110, 70)), arrow = arrow(length = unit(0.2, "cm")),
               col = col[1], alpha  = 0.7) +
  geom_segment(data = filter(Fmsy_sum, scen == "warm" & scen2 == "Resource (exp.)"), linetype = 2, 
               aes(x = Fm, xend = Fm, y = c(0, 0, 0), yend = c(120, 110, 70)), arrow = arrow(length = unit(0.2, "cm")),
               col = col[1], alpha  = 0.7) +
  geom_segment(data = filter(Fmsy_sum, scen == "warm" & scen2 == "Physio."), linetype = 1, 
               aes(x = Fm, xend = Fm, y = c(0, 0, 0), yend = c(120, 110, 70)), arrow = arrow(length = unit(0.2, "cm")), 
               col = col[1], alpha  = 0.7) +
  geom_segment(data = filter(Fmsy_sum, scen == "cold" & scen2 == "Physio. + Resource (exp.)"), linetype = 1, 
               aes(x = Fm, xend = Fm, y = c(0, 0, 0), yend = c(120, 110, 70)), arrow = arrow(length = unit(0.2, "cm")), 
               col = col[2], alpha  = 0.7) +
  coord_cartesian(expand = 0) +
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text = element_text(size = 10),
                                      legend.position = "bottom",
                                      aspect.ratio = 3/4,
                                      legend.text = element_text(size = 6),
                                      legend.title = element_text(size = 10))

ggsave("baltic/figures/FMSY_warm_cold.png", width = 6.5, height = 6.5, dpi = 600)


p2 <- Fmsy %>% 
  filter(type == "ssb") %>% 
  filter(biomass > 0.001) %>% 
  ggplot(., aes(Fm, biomass*240.342, linetype = scen2, color = scen)) + 
  geom_line(alpha = 0.8, size = 1.2) +
  facet_wrap(~ species, scales = "free") +
  scale_color_manual(values = c(col[2], col[1], col[3]), 
                     labels = c(expression("T"[ref]), 
                                expression(paste("T"[ref], "+2", degree*C)))) +
  labs(x = "Fishing mortality [1/year]", 
       y = "SSB [1000 tonnes]",
       color = "Scenario",
       linetype = "Metric") +
  coord_cartesian(expand = 0) +
  NULL

pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text = element_text(size = 10),
                                      legend.position = "bottom",
                                      aspect.ratio = 3/4,
                                      legend.text = element_text(size = 6),
                                      legend.title = element_text(size = 10))

ggsave("baltic/figures/supp/SSB_warm_cold.png", width = 6.5, height = 6.5, dpi = 600)



# C. HEATMAP EFFORT~TEMPERATURE FOR LOOP ===========================================
simFMSY <- projectEffort_m[177, ]
baseEffort <- simFMSY

baseTemp <- t_ref
t_max <- 100

#**** Project reference scenario ===================================================
ref <- project(pars_res_phys,
               dt = dt,
               effort = baseEffort,
               temperature = rep(baseTemp, t_max),
               diet_steps = 10,
               t_max = t_max)

refYield <- data.frame(Yield = c(getYield(ref)[dim(ref@effort)[1], 1],
                                 getYield(ref)[dim(ref@effort)[1], 2],
                                 getYield(ref)[dim(ref@effort)[1], 3]),
                       Species = ref@params@species_params$species)


#**** for loop through different fishing effort (with temp. dep. resource) ============
data_list <- list()
temp <- seq(0.75, 1.25, 0.02)
eff <- seq(0.25, 3, 0.1) # Factor for scaling fishing mortality

temp_eff <- data.frame(expand.grid(eff = eff, temp = temp))
iter <- seq(from = 1, to = nrow(temp_eff))


## Cod
baseEffort_var <- baseEffort

for (i in iter) {

  baseEffort_var[1] <- baseEffort[1] * temp_eff$eff[i] #
  baseTemp_var <- baseTemp * temp_eff$temp[i]

  proj <- project(pars_res_phys,
                  dt = dt,
                  effort = baseEffort_var,
                  temperature = rep(baseTemp_var, t_max),
                  diet_steps = 10,
                  t_max = t_max)

  # Extract yield at last iteration
  proYield <- data.frame(Yield = getYield(proj)[dim(proj@effort)[1], 1],
                         Species = "Cod",
                         Fm = proj@effort[dim(proj@effort)[1], 1],
                         temp = proj@temperature[dim(proj@temperature)[1], 1],
                         Fm_scal = temp_eff$eff[i],
                         temp_scal = temp_eff$temp[i])

  proYield$Yield_rel <- proYield$Yield / refYield$Yield[1]

  data_list[[i]] <- proYield

}

# Add data
big_yield_data_cod <- dplyr::bind_rows(data_list)


## Herring
baseEffort_var <- baseEffort

for (i in iter) {

  baseEffort_var[2] <- baseEffort[2] * temp_eff$eff[i] #
  baseTemp_var <- baseTemp * temp_eff$temp[i]

  proj <- project(pars_res_phys,
                  dt = dt,
                  effort = baseEffort_var,
                  temperature = rep(baseTemp_var, t_max),
                  diet_steps = 10,
                  t_max = t_max)

  # Extract yield at last iteration
  proYield <- data.frame(Yield = getYield(proj)[dim(proj@effort)[1], 3],
                         Species = "Herring",
                         Fm = proj@effort[dim(proj@effort)[1], 3],
                         temp = proj@temperature[dim(proj@temperature)[1], 1],
                         Fm_scal = temp_eff$eff[i],
                         temp_scal = temp_eff$temp[i])

  proYield$Yield_rel <- proYield$Yield / refYield$Yield[3]

  data_list[[i]] <- proYield

}

# Add data
big_yield_data_herring <- dplyr::bind_rows(data_list)


## Sprat
baseEffort_var <- baseEffort

for (i in iter) {

  baseEffort_var[2] <- baseEffort[2] * temp_eff$eff[i] #
  baseTemp_var <- baseTemp * temp_eff$temp[i]

  proj <- project(pars_res_phys,
                  dt = dt,
                  effort = baseEffort_var,
                  temperature = rep(baseTemp_var, t_max),
                  diet_steps = 10,
                  t_max = t_max)

  # Extract yield at last iteration
  proYield <- data.frame(Yield = getYield(proj)[dim(proj@effort)[1], 2],
                         Species = "Sprat",
                         Fm = proj@effort[dim(proj@effort)[1], 2],
                         temp = proj@temperature[dim(proj@temperature)[1], 1],
                         Fm_scal = temp_eff$eff[i],
                         temp_scal = temp_eff$temp[i])

  proYield$Yield_rel <- proYield$Yield / refYield$Yield[2]

  data_list[[i]] <- proYield

}

# Add data
big_yield_data_sprat <- dplyr::bind_rows(data_list)

# Merge all data

big_yield_data <- rbind(big_yield_data_cod, big_yield_data_herring, big_yield_data_sprat)

big_yield_data %>% filter(temp_scal == 1)


# All together
p3 <- ggplot(big_yield_data, aes(temp_scal, Fm_scal, fill = Yield_rel)) +
  geom_tile(color = NA) +
  facet_grid(~ Species, scales = "free") +
  scale_fill_viridis() +
  labs(x = c(expression("Temperature factor to T"[ref])),
       y = "Fishing mortality relative\n to MSSM FMSY",
       fill = "Yield relative to FMSY +\nconstant temp.") +
  coord_cartesian(expand = 0) +
  NULL

pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text = element_text(size = 10),
                                      aspect.ratio = 3/4,
                                      legend.position = "bottom")

ggsave("baltic/figures/supp/yield_heat_v1.png", width = 6.5, height = 6.5, dpi = 600)


# All separate
p4 <- big_yield_data %>% filter(Species == "Cod") %>%
ggplot(., aes(temp_scal, Fm_scal, fill = Yield_rel)) +
  geom_tile(color = NA) +
  scale_fill_viridis() +
  labs(x = c(expression("Temperature factor to T"[ref])),
       y = "Fishing mortality relative\n to MSSM FMSY",
       fill = "Yield relative to FMSY +\nconstant temp.") +
  coord_cartesian(expand = 0) +
  ggtitle("Cod") +
  NULL

pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 10),
                                      axis.text = element_text(size = 8),
                                      legend.text = element_text(size = 6),
                                      aspect.ratio = 3/4,
                                      #legend.position = "bottom",
                                      legend.key.height = unit(0.75, "line"),
                                      legend.key.width = unit(0.5, "line"))


p5 <- big_yield_data %>% filter(Species == "Herring") %>%
ggplot(., aes(temp_scal, Fm_scal, fill = Yield_rel)) +
  geom_tile(color = NA) +
  scale_fill_viridis() +
  labs(x = c(expression("Temperature factor to T"[ref])),
       y = "Fishing mortality relative\n to MSSM FMSY",
       fill = "Yield relative to FMSY +\nconstant temp.") +
  coord_cartesian(expand = 0) +
  ggtitle("Herring") +
  NULL

pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 8),
                                       legend.text = element_text(size = 6),
                                       aspect.ratio = 3/4,
                                       #legend.position = "bottom",
                                       legend.key.height = unit(0.75, "line"),
                                       legend.key.width = unit(0.5, "line"))


p6 <- big_yield_data %>% filter(Species == "Sprat") %>%
  ggplot(., aes(temp_scal, Fm_scal, fill = Yield_rel)) +
  geom_tile(color = NA) +
  scale_fill_viridis() +
  labs(x = c(expression("Temperature factor to T"[ref])),
       y = "Fishing mortality relative\n to MSSM FMSY",
       fill = "Yield relative to FMSY +\nconstant temp.") +
  coord_cartesian(expand = 0) +
  ggtitle("Sprat") +
  NULL

pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 8),
                                       legend.text = element_text(size = 6),
                                       aspect.ratio = 3/4,
                                       #legend.position = "bottom",
                                       legend.key.height = unit(0.75, "line"),
                                       legend.key.width = unit(0.5, "line"))

pWord4 / pWord5 / pWord6

ggsave("baltic/figures/supp/yield_heat_v2.png", width = 6.5, height = 6.5, dpi = 600)


# The below code changes all fishing mortalities at the same time. I instead want to
# hold all other species at their FMSY from the size spectrum model (above code)

# data_list <- list()
# temp <- seq(0.75, 1.25, 0.05)
# eff <- seq(0.1, 3, length.out = length(temp)) # Factor for scaling fishing mortality
# temp_eff <- expand.grid(data.frame(eff = eff, temp = temp))
# iter <- seq(from = 1, to = nrow(temp_eff))
# 
# for (i in iter) {
#   
#   baseEffort_var <- baseEffort
#   baseEffort_var <- baseEffort * temp_eff$eff[i] #
#   #baseEffort_var <- baseEffort * temp_eff$eff[i] #
#   baseTemp_var <- baseTemp * temp_eff$temp[i]
#   
#   proj <- project(pars_res_phys,
#                   dt = dt,
#                   effort = baseEffort_var,
#                   temperature = rep(baseTemp_var, t_max),
#                   diet_steps = 10,
#                   t_max = t_max)
#   
#   # Extract yield at last iteration
#   proYield <- data.frame(Yield = c(getYield(proj)[dim(proj@effort)[1], 1],
#                                    getYield(proj)[dim(proj@effort)[1], 2],
#                                    getYield(proj)[dim(proj@effort)[1], 3]),
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
# # Add data 
# big_yield_data <- dplyr::bind_rows(data_list)
# 
# big_yield_data$Species <- factor(big_yield_data$Species, levels = c("Sprat", "Herring", "Cod"))
# 
# filter(big_yield_data, Species == "Cod")
# filter(big_yield_data, Species == "Sprat")
# 
# ggplot(big_yield_data, aes(temp_scal, Fm_scal, fill = Yield_rel)) +
#   geom_tile(color = NA) +
#   facet_grid(~ Species, scales = "free") +
#   theme_classic(base_size = 12) +
#   scale_fill_viridis() +
#   labs(x = c(expression("Temperature factor to T"[ref])),
#        y = "Fishing mortality relative\n to MSSM FMSY",
#        fill = "Yield relative to\naverage FMSY +\nconstant temp.") +
#   coord_cartesian(expand = 0) +
#   theme(aspect.ratio = 3/4,
#         legend.position = "bottom") +
#   NULL

#ggsave("baltic/figures/supp/yield_heat_ssmFMSY.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


#** SOME TESTS =====================================================================
plotYield(w) + xlim(50, 75)
plotBiomass(w) + xlim(50, 75)

tail(getBiomass(w), 10)
tail(getSSB(w), 10)
tail(getYield(w), 10)

test_eff <- effort
test_eff[1:3] <- 1

test <- project(pars_with_res,
                dt = dt,
                effort = test_eff,
                temperature = rep(pars_with_res@t_ref, t_max),
                diet_steps = 10,
                t_max = t_max)

plot(test)
tail(getSSB(w), 10)
tail(getYield(w), 10)

