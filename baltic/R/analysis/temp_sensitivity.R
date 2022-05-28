#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2022.04.26: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v3 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. FMSY at two different temperatures
# 
# C. Analysis of changes in yield with fishing warming (Heatmap)
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
#devtools::load_all(".")

# Install the specific mizer version from github
# devtools::install_github("maxlindmark/mizer-rewiring", ref = "rewire-temp") 
library(mizer)

# Print package versions
# print(sessionInfo())
# other attached packages:
# mizer_1.1          testthat_2.3.0     patchwork_0.0.1    dplyr_0.8.3        tidyr_1.0.0        
# viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       RCurl_1.95-4.12   
# bitops_1.0-6       RColorBrewer_1.1-2 devtools_2.2.1     usethis_1.5.1      ggplot2_3.2.1  

# Load function for extracting abundance-at-size
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getSpectra.R", ssl.verifypeer = FALSE)
eval(parse(text = func))

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
t1 <- params@species_params
t2 <- params@species_params
t3 <- params@species_params

# t0: only metabolism
# t1: only intake
# t2: only mortality

t1$ea_met <- mean(ea$met)
t1$ea_int <- 0
t1$ea_mor <- 0

t2$ea_met <- 0
t2$ea_int <- mean(ea$int)
t2$ea_mor <- 0

t3$ea_met <- 0
t3$ea_int <- 0
t3$ea_mor <- mean(ea$mor)

# Without resource dependence
pars_nores_t1 <- MizerParams(t1, 
                             ea_gro = 0,
                             ea_car = 0, # -ea$gro[i] 
                             kappa_ben = kappa_ben,
                             kappa = kappa,
                             w_bb_cutoff = w_bb_cutoff,
                             w_pp_cutoff = w_pp_cutoff,
                             r_pp = r_pp,
                             r_bb = r_bb,
                             t_ref = t_ref)

pars_nores_t2 <- MizerParams(t2, 
                             ea_gro = 0,
                             ea_car = 0, # -ea$gro[i] 
                             kappa_ben = kappa_ben,
                             kappa = kappa,
                             w_bb_cutoff = w_bb_cutoff,
                             w_pp_cutoff = w_pp_cutoff,
                             r_pp = r_pp,
                             r_bb = r_bb,
                             t_ref = t_ref)

pars_nores_t3 <- MizerParams(t3, 
                             ea_gro = 0,
                             ea_car = 0, # -ea$gro[i] 
                             kappa_ben = kappa_ben,
                             kappa = kappa,
                             w_bb_cutoff = w_bb_cutoff,
                             w_pp_cutoff = w_pp_cutoff,
                             r_pp = r_pp,
                             r_bb = r_bb,
                             t_ref = t_ref)

# With resource dependence
pars_res_phys_t1 <- MizerParams(t1, 
                                ea_gro = mean(ea$gro),
                                ea_car = mean(ea$car), # -ea$gro[i] 
                                kappa_ben = kappa_ben,
                                kappa = kappa,
                                w_bb_cutoff = w_bb_cutoff,
                                w_pp_cutoff = w_pp_cutoff,
                                r_pp = r_pp,
                                r_bb = r_bb,
                                t_ref = t_ref)

pars_res_phys_t2 <- MizerParams(t2, 
                                ea_gro = mean(ea$gro),
                                ea_car = mean(ea$car), # -ea$gro[i] 
                                kappa_ben = kappa_ben,
                                kappa = kappa,
                                w_bb_cutoff = w_bb_cutoff,
                                w_pp_cutoff = w_pp_cutoff,
                                r_pp = r_pp,
                                r_bb = r_bb,
                                t_ref = t_ref)

pars_res_phys_t3 <- MizerParams(t3, 
                                ea_gro = mean(ea$gro),
                                ea_car = mean(ea$car), # -ea$gro[i] 
                                kappa_ben = kappa_ben,
                                kappa = kappa,
                                w_bb_cutoff = w_bb_cutoff,
                                w_pp_cutoff = w_pp_cutoff,
                                r_pp = r_pp,
                                r_bb = r_bb,
                                t_ref = t_ref)

# Define temperature & effort
baseEffort <- projectEffort_m[177, ]

baseTemp <- t_ref
t_max <- 200


#**** Project scenarios ============================================================
# Reference scenario
proj_ref <- project(pars_nores_t1,
                    dt = dt,
                    effort = baseEffort,
                    temperature = rep(baseTemp, t_max),
                    diet_steps = 10,
                    t_max = t_max)

# Without resource dependence
proj_nores_t1 <- project(pars_nores_t1,
                         dt = dt,
                         effort = baseEffort,
                         temperature = rep(baseTemp + 2, t_max),
                         diet_steps = 10,
                         t_max = t_max)

proj_nores_t2 <- project(pars_nores_t2,
                         dt = dt,
                         effort = baseEffort,
                         temperature = rep(baseTemp + 2, t_max),
                         diet_steps = 10,
                         t_max = t_max)

proj_nores_t3 <- project(pars_nores_t3,
                         dt = dt,
                         effort = baseEffort,
                         temperature = rep(baseTemp + 2, t_max),
                         diet_steps = 10,
                         t_max = t_max)

# With resource dependence
proj_res_t1 <- project(pars_res_phys_t1,
                       dt = dt,
                       effort = baseEffort,
                       temperature = rep(baseTemp + 2, t_max),
                       diet_steps = 10,
                       t_max = t_max)

proj_res_t2 <- project(pars_res_phys_t2,
                       dt = dt,
                       effort = baseEffort,
                       temperature = rep(baseTemp + 2, t_max),
                       diet_steps = 10,
                       t_max = t_max)

proj_res_t3 <- project(pars_res_phys_t3,
                       dt = dt,
                       effort = baseEffort,
                       temperature = rep(baseTemp + 2, t_max),
                       diet_steps = 10,
                       t_max = t_max)

#**** Calculate growth relative to the reference scenario =========================
growth_ref <- getGrowth(proj_ref)

growth_nores_t1 <- getGrowth(proj_nores_t1) %>% mutate(scen = "met", resource = "no resource temp dep")
growth_nores_t2 <- getGrowth(proj_nores_t2) %>% mutate(scen = "int", resource = "no resource temp dep")
growth_nores_t3 <- getGrowth(proj_nores_t3) %>% mutate(scen = "mor", resource = "no resource temp dep")

growth_res_t1 <- getGrowth(proj_res_t1) %>% mutate(scen = "met", resource = "resource temp dep")
growth_res_t2 <- getGrowth(proj_res_t2) %>% mutate(scen = "int", resource = "resource temp dep")
growth_res_t3 <- getGrowth(proj_res_t3) %>% mutate(scen = "mor", resource = "resource temp dep")

growth_nores_t1$re_growth <- growth_nores_t1$value / growth_ref$value
growth_nores_t2$re_growth <- growth_nores_t2$value / growth_ref$value
growth_nores_t3$re_growth <- growth_nores_t3$value / growth_ref$value

growth_res_t1$re_growth <- growth_res_t1$value / growth_ref$value
growth_res_t2$re_growth <- growth_res_t2$value / growth_ref$value
growth_res_t3$re_growth <- growth_res_t3$value / growth_ref$value

# Merge into big data
growth_big_dat <- rbind(growth_nores_t1, growth_nores_t2, growth_nores_t3,
                        growth_res_t1, growth_res_t2, growth_res_t3)

growth_big_dat$scen <- as.factor(growth_big_dat$scen)


#**** Plot absolute and relative size-at-age =======================================
# Reorder factor levels
growth_big_dat$Species <- factor(growth_big_dat$Species, levels = c("Sprat", "Herring", "Cod"))

# Plot growth curves (mean and ribbons)
p1 <- ggplot(growth_big_dat, aes(x = Age, y = value, color = factor(scen), linetype = resource)) +
  geom_line() +
  labs(y = "Body mass (g)") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_line(data = growth_ref, aes(Age, value), color = "black", 
            inherit.aes = FALSE, size = 0.5, alpha = 0.7, linetype = "dashed") + 
  guides(linetype = "none", color = "none") +
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 10),
                                       legend.position = "bottom",
                                       legend.title = element_blank())

pWord1


p2 <- growth_big_dat %>% filter(Age > 0 & Age < 16) %>%
  ggplot(., aes(x = Age, y = re_growth, color = factor(scen), linetype = resource)) +
  facet_wrap(~Species, scales = "free_y") +
  geom_line() +
  labs(y = "Body mass relative to\nconstant temperature") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_viridis(discrete = TRUE, name = "Scenario") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linetype = 1)),
         linetype = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 1, size = 0.3, linetype = "dashed", color = "black") +
  NULL

pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 10),
                                       legend.position = "bottom",
                                       legend.title = element_blank())

pWord2

# Plot together
pWord1 / pWord2

ggsave("baltic/figures/supp/growth_project_temp_sens.png", width = 6.5, height = 6.5, dpi = 600)



