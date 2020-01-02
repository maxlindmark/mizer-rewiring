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

# Define temperature-scenarios
consTemp <- projectTemp$temperature
start <- 1997-1914
consTemp[start:137] <- t_ref


# B. EXTRACT FMSY FROM MODEL AT DIFFERENT TEMPERATURES =============================
# Here I just for-loop different F, extract Yield, plotover F. 
# I will increase each species F separately, keeping the others at their mean, as
# in the calibration method.
# NOTE: Here we use constant temperatures and therefore the results may not be 
# exactly translateable to the time-varying effort and temperature projections.

F_range <- seq(0.2, 1.4, 0.01) # Can decrease step later, becomes too slow now
t_max <- 100
index <- 1:length(F_range)

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
  
  r <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep(pars_with_res@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep((pars_with_res@t_ref + 2), t_max),
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

codFmsy <- dplyr::bind_rows(data_list)

codFmsy$species <- "Cod"

ggplot(codFmsy, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


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
  
  r <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep(pars_with_res@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep((pars_with_res@t_ref + 2), t_max),
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
               temperature = rep(pars_with_res@t_ref, t_max),
               diet_steps = 10,
               t_max = t_max)
  
  w <- project(pars_with_res,
               dt = dt,
               effort = effort,
               temperature = rep((pars_with_res@t_ref + 2), t_max),
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

sprFmsy <- dplyr::bind_rows(data_list)

sprFmsy$species <- "Sprat"

ggplot(sprFmsy, aes(Fm, biomass, linetype = type, color = scen)) + geom_line() 


#**** All together =================================================================
col <- RColorBrewer::brewer.pal(n = 5, "Dark2")
col <- RColorBrewer::brewer.pal(n = 3, "Set1")[1:2]

Fmsy <- rbind(codFmsy, sprFmsy, herFmsy)

Fmsy$species <- factor(Fmsy$species, levels = c("Sprat", "Herring", "Cod"))

# Get size-spectrum FMSY for all species:
Fmsy_sum <- Fmsy %>% 
  filter(type == "yield") %>% 
  group_by(species, scen) %>% 
  filter(biomass == max(biomass))

Fmsy %>% 
  filter(biomass > 0.001) %>% 
  ggplot(., aes(Fm, biomass, linetype = type, color = scen)) + 
  geom_line(alpha = 0.6, size = 1.2) +
  facet_wrap(~ species, scales = "free") +
  scale_color_manual(values = rev(col), labels = c("T_ref", "T_ref + 2C")) +
  theme_classic(base_size = 12) +
  labs(x = "Fishing mortality [1/year]", y = "Biomass",
       color = "Scenario",
       linetype = "Metric") +
  theme(aspect.ratio = 3/4,
        legend.position = "bottom") +
  geom_segment(data = filter(Fmsy_sum, scen == "warm"), linetype = 1, 
               aes(x = Fm, xend = Fm, y = 0, yend = c(0.45, 0.8, 0.45)), arrow = arrow(length = unit(0.3, "cm")), col = col[1]) +
  geom_segment(data = filter(Fmsy_sum, scen == "cold"), linetype = 1, 
               aes(x = Fm, xend = Fm, y = 0, yend = c(0.45, 0.8, 0.45)), arrow = arrow(length = unit(0.3, "cm")), col = col[2]) +
  NULL

#ggsave("baltic/figures/FMSY_warm_cold.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


## IN THE GROWTH FIGURE, SHOULD I ALSO DO A 0 ACTIVATION ENERGY AS A CONTRAST? THATS THE NULL ACTUALLY!!!

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


# C. HEATMAP EFFORT~TEMPERATURE FOR LOOP ===========================================
baseEffort <- projectEffort_m[137, ] # This is the average of assessed and ssm FMSY

simFMSY <- projectEffort_m[137, ] # This is the average of assessed and ssm FMSY
simFMSY[] <- c(0.4, 0.52, 0.92)
baseEffort <- simFMSY

baseTemp <- t_ref
t_max <- 100

#**** Project reference scenario ===================================================
ref <- project(pars_with_res, 
               dt = dt,
               effort = baseEffort,
               temperature = rep(baseTemp, t_max),
               diet_steps = 10,
               t_max = t_max)

refYield <- data.frame(Yield = c(getYield(ref)[dim(ref@effort)[1], 1],
                                 getYield(ref)[dim(ref@effort)[1], 2],
                                 getYield(ref)[dim(ref@effort)[1], 3]),
                       Species = ref@params@species_params$species)



# #**** for loop through different fishin effort (with temp dep resource) ============
data_list <- list()
eff <- seq(0.8, 1.2, 0.05) # Factor for scaling fishing mortality
temp <- seq(0.8, 1.2, 0.05)
temp_eff <- expand.grid(data.frame(eff = eff, temp = temp))
iter <- seq(from = 1, to = nrow(temp_eff))

# projectEffort_new <- projectEffort_m
# projectTemp_new <- projectTemp$temperature
# tt <- c()
# groj <- c()
# yield <- c()
# data_list <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  baseEffort_var <- baseEffort * temp_eff$eff[i]
  baseTemp_var <- baseTemp * temp_eff$temp[i]
  
  proj <- project(pars_with_res,
                  dt = dt,
                  effort = baseEffort_var,
                  temperature = rep(baseTemp_var, t_max),
                  diet_steps = 10,
                  t_max = t_max)
  
  # Extract yield at last iteration
  proYield <- data.frame(Yield = c(getYield(proj)[dim(proj@effort)[1], 1],
                                   getYield(proj)[dim(proj@effort)[1], 2],
                                   getYield(proj)[dim(proj@effort)[1], 3]),
                         Species = proj@params@species_params$species,
                         Fm = proj@effort[dim(proj@effort)[1], ],
                         temp = proj@temperature[dim(proj@temperature)[1], 1],
                         Fm_scal = temp_eff$eff[i],
                         temp_scal = temp_eff$temp[i])
  
  proYield$Yield_rel <- proYield$Yield / refYield$Yield

  data_list[[i]] <- proYield

}

str(data_list)

big_yield_data <- dplyr::bind_rows(data_list)

big_yield_data$Species <- factor(big_yield_data$Species, levels = c("Sprat", "Herring", "Cod"))

ggplot(big_yield_data, aes(temp_scal, Fm_scal, fill = Yield_rel)) +
  geom_tile(color = NA) +
  facet_grid(~ Species, scales = "free") +
  theme_classic(base_size = 12) +
  scale_fill_viridis() +
  labs(x = "Temperature relative to\nRCP8.5 projection",
       #y = "Fishing mortality relative\nto average FMSY",
       y = "Fishing mortality relative\nto MSSM FMSY",
       fill = "Yield relative to\naverageFMSY +\nconstant temp.") +
  coord_cartesian(expand = 0) +
  theme(aspect.ratio = 3/4,
        legend.position = "bottom") +
  NULL

#ggsave("baltic/figures/supp/yield_heat.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")
#ggsave("baltic/figures/supp/yield_heat_ssmFMSY.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")

 