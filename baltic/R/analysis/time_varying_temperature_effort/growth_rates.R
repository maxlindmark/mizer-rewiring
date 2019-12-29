#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. Individual growth trajectories in default fishing scenario and random temp-effects
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
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getGrowth.R", ssl.verifypeer = FALSE)
eval(parse(text = func))


#**** Read in parameters and data ==================================================
# Read in params object
params <- readRDS("baltic/params/mizer_param_calib.rds")

# Read in activation energy data frame
ea <- read.csv("baltic/params/samples_activation_energy.csv")[, 2:6]

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
# we get comparable starting values for the models so that temperature is the only
# "treatment"
start <- 1997-1914
consTemp[start:137] <- t_ref

# Plot
col <- RColorBrewer::brewer.pal("Dark2", n = 5)
tempScen <- data.frame(Temperature = c(consTemp, projectTemp$temperature),
                       Scenario = rep(c("warming", "no warming"), each = length(consTemp)),
                       Year = 1:length(consTemp) + 1913)

ggplot(tempScen, aes(Year, Temperature, color = Scenario, linetype = Scenario)) +
  geom_line(alpha = 0.8, size = 1.4) +
  theme_classic(base_size = 25) +
  scale_color_manual(values = rev(col)) +
  theme(legend.position=c(.2,.75),
        aspect.ratio = 3/4) +
  NULL

#ggsave("baltic/figures/supp/temperature_scenarios.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


# B. TEMP-DRIVEN CHANGE IN SIZE-AT-AGE =============================================
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


#**** for loop (with resource) =====================================================
sim <- 1:200

t <- c()
tt <- c()
groj <- c()
growth <- c()
data_list_with_res <- list()

for (i in sim) {
  
  t <- params@species_params
  
  t$ea_met <- ea$met[i]
  t$ea_int <- ea$int[i]
  t$ea_mor <- ea$mor[i]
  
  #t$ca_int <- -0.004 # Here we just use the fixed values
  #t$ca_met <- 0.001 # Here we just use the fixed values
  
  tt <- MizerParams(t, 
                    ea_gro = ea$gro[i],
                    ea_car = ea$car[i], # -ea$gro[i] 
                    kappa_ben = kappa_ben,
                    kappa = kappa,
                    w_bb_cutoff = w_bb_cutoff,
                    w_pp_cutoff = w_pp_cutoff,
                    r_pp = r_pp,
                    r_bb = r_bb,
                    t_ref = t_ref)
  
  proj <- project(tt, 
                  dt = dt,
                  effort = projectEffort_m,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max)   
  
  growth <- getGrowth(proj)
  
  growth$ea_met <- proj@params@species_params$ea_met[1]
  growth$ea_mor <- proj@params@species_params$ea_mor[1]
  growth$ea_int <- proj@params@species_params$ea_int[1]
  
  growth$ea_gro <- proj@params@ea_gro
  growth$ea_car <- proj@params@ea_car
  
  growth$sim <- i
  
  growth$re_growth <- growth$value / refGrowth$value
  
  data_list_with_res[[i]] <- growth
  
}

big_growth_data_w_r <- dplyr::bind_rows(data_list_with_res)


#**** for loop (no resource) =======================================================
sim <- 1:200

t <- c()
tt <- c()
groj <- c()
growth <- c()
data_list_no_res <- list()

for (i in sim) {
  
  t <- params@species_params
  
  t$ea_met <- ea$met[i]
  t$ea_int <- ea$int[i]
  t$ea_mor <- ea$mor[i]
  
  #t$ca_int <- -0.004 # Here we just use the fixed values
  #t$ca_met <- 0.001 # Here we just use the fixed values
  
  tt <- MizerParams(t, 
                    ea_gro = 0,
                    ea_car = 0,
                    kappa_ben = kappa_ben,
                    kappa = kappa,
                    w_bb_cutoff = w_bb_cutoff,
                    w_pp_cutoff = w_pp_cutoff,
                    r_pp = r_pp,
                    r_bb = r_bb,
                    t_ref = t_ref)
  
  proj <- project(tt, 
                  dt = dt,
                  effort = projectEffort_m,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max)   
  
  growth <- getGrowth(proj)
  
  growth$ea_met <- proj@params@species_params$ea_met[1]
  growth$ea_mor <- proj@params@species_params$ea_mor[1]
  growth$ea_int <- proj@params@species_params$ea_int[1]
  
  growth$ea_gro <- proj@params@ea_gro
  growth$ea_car <- proj@params@ea_car
  
  growth$sim <- i
  
  growth$re_growth <- growth$value / refGrowth$value
  
  data_list_no_res[[i]] <- growth
  
}

big_growth_data_no_r <- dplyr::bind_rows(data_list_no_res)


#**** Plot growth rates ============================================================
big_growth_data_w_r$scen <- "With resource temp. dep."
big_growth_data_no_r$scen <- "No resource temp. dep."

big_growth_data_w_r$sim <- paste("wr", big_growth_data_w_r$sim, sep = "")
big_growth_data_no_r$sim <- paste("nr", big_growth_data_no_r$sim, sep = "")

big_growth_data <- rbind(big_growth_data_w_r, big_growth_data_no_r)

big_growth_data$scen <- as.factor(big_growth_data$scen)

# Calculate means, max and min for plotting
mean_dat <- big_growth_data %>%
  dplyr::group_by(Species, Age, scen) %>% 
  dplyr::summarise(mean_val = mean(value),
                   max_val = max(value),
                   min_val = min(value))

col <- RColorBrewer::brewer.pal("Dark2", n = 5)

# Reorder factor levels
mean_dat$Species <- factor(mean_dat$Species, levels = c("Sprat", "Herring", "Cod"))

# Plot growth curves
p1 <- ggplot(mean_dat, aes(x = Age, ymin = min_val, ymax = max_val, fill = factor(scen))) +
  geom_line(data = mean_dat, aes(Age, mean_val, color = factor(scen)),
            inherit.aes = FALSE, size = 0.5) +
  geom_ribbon(alpha = 0.15, color = NA) +  
  labs(y = "Body mass (g)") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = rev(col)) +
  scale_fill_manual(values = rev(col)) +
  theme_classic(base_size = 14) +
  guides(color = FALSE, fill = FALSE) +
  geom_line(data = refGrowth, aes(Age, value), color = "black", 
           inherit.aes = FALSE, size = 0.5, alpha = 0.7, linetype = "dashed") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  NULL

p1 <- big_growth_data %>% filter(Age > 0 & Age < 16) %>% 
  ggplot(., aes(x = Age, y = value, color = factor(scen), group = sim)) +
  geom_line(size = 0.5, alpha = 0.03) +
  labs(y = "Body mass (g)") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = rev(col)) +
  theme_classic(base_size = 14) +
  guides(color = FALSE, fill = FALSE) +
  geom_line(data = filter(refGrowth, Age < 16), aes(Age, value), color = "black", 
            inherit.aes = FALSE, size = 0.5, alpha = 0.7, linetype = "dashed") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  NULL

p1

# Testing relative growth rates are correct (seems like larger difference there than in actual growth rates)
# big_growth_data %>% 
#   filter(Age > 0 & Age < 16 & Species == "Cod") %>% 
#   ggplot(., aes(x = Age, y = value, color = factor(scen), group = sim)) +
#   geom_line(size = 0.5, alpha = 0.03) +
#   labs(y = "Body mass (g)") +
#   scale_y_continuous(expand = c(0, 0)) + 
#   scale_color_manual(values = rev(col)) +
#   theme_classic(base_size = 14) +
#   guides(color = FALSE, fill = FALSE) +
#   geom_line(data = filter(refGrowth, Age < 16 & Species == "Cod"), aes(Age, value), color = "black", 
#             inherit.aes = FALSE, size = 0.5, alpha = 0.7, linetype = "dashed") +
#   theme(legend.position = "bottom",
#         legend.title = element_blank()) +
#   ylim(400, 1100) +
#   xlim(3, 4) +
#   NULL
  

# Plot relative growth curves
rel_dat <- big_growth_data %>%
  dplyr::group_by(Species, Age, scen) %>% 
  dplyr::summarise(mean_val = mean(re_growth),
                   max_val = max(re_growth),
                   min_val = min(re_growth))

rel_dat$Species <- factor(rel_dat$Species, levels = c("Sprat", "Herring", "Cod"))

p2 <- rel_dat %>% filter(Age > 0) %>% 
ggplot(., aes(x = Age, ymin = min_val, ymax = max_val, fill = factor(scen))) +
  geom_line(data = filter(rel_dat, Age > 0), aes(Age, mean_val, color = factor(scen)),
            inherit.aes = FALSE, size = 0.5) +
  geom_ribbon(alpha = 0.15, color = NA) +  
  labs(y = "Body mass relative to\nconstant temperature") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0), limits = c(0.95, 2.2)) + 
  scale_color_manual(values = rev(col)) +
  scale_fill_manual(values = rev(col)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  geom_hline(yintercept = 1, size = 0.3, linetype = "dashed", color = "black") +
  NULL

p2 <- big_growth_data %>% filter(Age > 0 & Age < 16) %>% 
  ggplot(., aes(x = Age, y = re_growth, color = factor(scen), group = sim)) +
  geom_line(alpha = 0.05, size = 0.5) +
  geom_line(data = filter(rel_dat, Age > 0 & Age < 16), aes(Age, mean_val, color = factor(scen)),
            inherit.aes = FALSE, size = 1, linetype = 2) +
  labs(y = "Body mass relative to\nconstant temperature") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0), limits = c(0.95, 2.2)) + 
  scale_color_manual(values = rev(col)) +
  scale_fill_manual(values = rev(col)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  geom_hline(yintercept = 1, size = 0.3, linetype = "dashed", color = "black") +
  NULL

p2

# Plot together
p1/p2

#ggsave("baltic/figures/growth_project.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")

# What are the activation energies in the MAX scenarios?
maxg <- big_growth_data %>%
  dplyr::filter(Age > 1) %>% # Need this since all length at age 0 are the same
  dplyr::group_by(Species, Age, scen) %>% 
  dplyr::filter(value == max(value)) %>% 
  droplevels()

nrow(maxg)
nrow(big_growth_data)
length(unique(maxg$sim))
length(unique(big_growth_data$sim))

ggplot(maxg, aes(Age, value, color = sim)) + 
  geom_line() + 
  facet_wrap(~ Species, scales = "free")

# Ok, so here the parameters are from iteration 137 and 173 for no resource and with 
# resource temperature-dependence
ea[136, ] # no resource temp dep
ea[62, ] # with resource temp dep

# What are the activation energies in the MIN scenarios?
maxg <- big_growth_data %>%
  dplyr::filter(Age > 1) %>% # Need this since all length at age 0 are the same
  dplyr::group_by(Species, Age, scen) %>% 
  dplyr::filter(value == min(value)) %>% 
  droplevels()

ggplot(maxg, aes(Age, value, color = sim)) + 
  geom_line() + 
  facet_wrap(~ Species, scales = "free")

# Ok, so here the parameters are from iteration 197 and 189 for no resource and with 
# resource temperature-dependence
ea[197, ] # no resource
ea[197, ] # with resource


#**** Testing I can reproduce a really bad example =================================
# (even without negative effects on carrying capacity)
t <- params@species_params

t$ea_met <- 0.8
t$ea_int <- 0.3
t$ea_mor <- 0.8

t$ca_int <- -0.004 # Here we just use the fixed values
t$ca_met <- 0.001 # Here we just use the fixed values

tt <- MizerParams(t, 
                  ea_gro = 0,
                  ea_car = 0,
                  kappa_ben = kappa_ben,
                  kappa = kappa,
                  w_bb_cutoff = w_bb_cutoff,
                  w_pp_cutoff = w_pp_cutoff,
                  r_pp = r_pp,
                  r_bb = r_bb,
                  t_ref = t_ref)

proj <- project(tt, 
                dt = dt,
                effort = projectEffort_m,
                temperature = projectTemp$temperature,
                diet_steps = 10,
                t_max = t_max)   

growth <- getGrowth(proj)

refGrowth2 <- refGrowth
refGrowth2$scen <- "ref"
growth$scen <- "bad"

kek <- rbind(refGrowth2, growth)

ggplot(kek, aes(Age, value, color = scen)) +
  geom_line() +
  facet_wrap(~ Species, scales = "free")

