#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. Simulate and extract individual growth trajectories in default fishing scenario 
# and random temp-effects
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
# func <- 
#   getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/raincloudPlot.R", 
#          ssl.verifypeer = FALSE)
# eval(parse(text = func))


#**** Read in parameters and data ==================================================
# Read in params object
params <- readRDS("baltic/params/mizer_param_calib.rds")

# Read in activation energy data frame
ea <- read.csv("baltic/params/samples_activation_energy.csv")#[, 2:6]
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

# The time series starts in 1874. From 1997 (mid point in calibration time), we want
# to fix the temperature at the mean of the calibration, i.e. t_ref. This ensures
# we get comparable starting values for the models so that temperature is the only "treatment"
start <- 1997-1874
consTemp[start:177] <- t_ref

# Plot
col <- RColorBrewer::brewer.pal("Dark2", n = 5)
col <- RColorBrewer::brewer.pal("Set1", n = 3)[1:2]

tempScen <- data.frame(Temperature = c(consTemp, projectTemp$temperature),
                       Scenario = rep(c("no warming", "warming"), each = length(consTemp)),
                       Year = 1:length(consTemp) + 1873)

tempScen %>% filter(Year < 2003 & Year > 1991) %>% summarize(mean_temp = mean(Temperature))

ggplot(tempScen, aes(Year, Temperature, color = Scenario, linetype = Scenario)) +
  geom_line(alpha = 0.8, size = 1.4) +
  scale_color_manual(values = rev(col)) +
  ylab("Temperature") +
  NULL


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
# sim <- 1:200
# 
# t <- c()
# tt <- c()
# groj <- c()
# growth <- c()
# data_list_with_res_no_phys_barn <- list()
# mean_weight_list_with_res_no_phys_barn <- list()
# 
# for (i in sim) {
# 
# t <- params@species_params
# 
# t$ea_met <- 0
# t$ea_int <- 0
# t$ea_mor <- 0
# 
# tt <- MizerParams(t, 
#                   ea_gro = 0, #0.41, #0.43, # ea$gro[i],
#                   ea_car = ea$b_car[i], #-0.41, # 0,  # -0.41, # ea$car[i], # -ea$gro[i] 
#                   kappa_ben = kappa_ben,
#                   kappa = kappa,
#                   w_bb_cutoff = w_bb_cutoff,
#                   w_pp_cutoff = w_pp_cutoff,
#                   r_pp = r_pp,
#                   r_bb = r_bb,
#                   t_ref = t_ref)
# 
# proj <- project(tt, 
#                 dt = dt,
#                 effort = projectEffort_m,
#                 temperature = projectTemp$temperature,
#                 diet_steps = 10) 
# 
# growth <- getGrowth(proj)
# 
# growth$ea_met <- proj@params@species_params$ea_met[1]
# growth$ea_mor <- proj@params@species_params$ea_mor[1]
# growth$ea_int <- proj@params@species_params$ea_int[1]
# 
# growth$ea_gro <- proj@params@ea_gro
# growth$ea_car <- proj@params@ea_car
# 
# growth$sim <- 1 # i
# 
# growth$re_growth <- growth$value / refGrowth$value
# 
# data_list_with_res_no_phys_barn[[i]] <- growth
#  
# mean_weight_list_with_res_no_phys_barn[[i]] <- data.frame(getSpeciesMeanWeight(proj)[nrow(projectEffort_m), ] )
# 
# }
# 
# big_growth_data_w_r_n_p_barn <- dplyr::bind_rows(data_list_with_res_no_phys_barn)
# big_mean_weight_data_w_r_n_p_barn <- dplyr::bind_rows(mean_weight_list_with_res_no_phys_barn)
# big_mean_weight_data_w_r_n_p_barn <- big_mean_weight_data_w_r_n_p_barn %>% 
#   dplyr::rename("mean_weight" = "getSpeciesMeanWeight.proj..nrow.projectEffort_m....")
# big_mean_weight_data_w_r_n_p_barn$species <- rep(ref@params@species_params$species, 200)


#**** for loop (Barnes - with resource) =====================================================
# sim <- 1:200
# 
# t <- c()
# tt <- c()
# groj <- c()
# growth <- c()
# data_list_with_res_barn <- list()
# mean_weight_list_with_res_barn <- list()
# 
# for (i in sim) {
#   
#   t <- params@species_params
#   
#   t$ea_met <- ea$met[i]
#   t$ea_int <- ea$int[i]
#   t$ea_mor <- ea$mor[i]
#   
#   tt <- MizerParams(t, 
#                     ea_gro = 0, #ea$gro[i],
#                     ea_car = ea$b_car[i],
#                     kappa_ben = kappa_ben,
#                     kappa = kappa,
#                     w_bb_cutoff = w_bb_cutoff,
#                     w_pp_cutoff = w_pp_cutoff,
#                     r_pp = r_pp,
#                     r_bb = r_bb,
#                     t_ref = t_ref)
#   
#   proj <- project(tt, 
#                   dt = dt,
#                   effort = projectEffort_m,
#                   temperature = projectTemp$temperature,
#                   diet_steps = 10,
#                   t_max = t_max)   
#   
#   growth <- getGrowth(proj)
#   
#   growth$ea_met <- proj@params@species_params$ea_met[1]
#   growth$ea_mor <- proj@params@species_params$ea_mor[1]
#   growth$ea_int <- proj@params@species_params$ea_int[1]
#   
#   growth$ea_gro <- proj@params@ea_gro
#   growth$ea_car <- proj@params@ea_car
#   
#   growth$sim <- i
#   
#   growth$re_growth <- growth$value / refGrowth$value
#   
#   data_list_with_res_barn[[i]] <- growth
#   
#   mean_weight_list_with_res_barn[[i]] <- data.frame(getSpeciesMeanWeight(proj)[nrow(projectEffort_m), ] )
#   
# }
# 
# big_growth_data_w_r_barn <- dplyr::bind_rows(data_list_with_res_barn)
# big_mean_weight_data_w_r_barn <- dplyr::bind_rows(mean_weight_list_with_res_barn)
# big_mean_weight_data_w_r_barn <- big_mean_weight_data_w_r_barn %>% 
#   dplyr::rename("mean_weight" = "getSpeciesMeanWeight.proj..nrow.projectEffort_m....")
# big_mean_weight_data_w_r_barn$species <- rep(ref@params@species_params$species, 200)



#**** for loop (with resource, MTE) =====================================================
sim <- 1:200

t <- c()
tt <- c()
groj <- c()
growth <- c()
data_list_with_res <- list()
mean_weight_list_with_res <- list()

for (i in sim) {
  
  t <- params@species_params
  
  t$ea_met <- ea$met[i]
  t$ea_int <- ea$int[i]
  t$ea_mor <- ea$mor[i]
  
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
  
  mean_weight_list_with_res[[i]] <- data.frame(getSpeciesMeanWeight(proj)[nrow(projectEffort_m), ] )
  
}

big_growth_data_w_r <- dplyr::bind_rows(data_list_with_res)
big_mean_weight_data_w_r <- dplyr::bind_rows(mean_weight_list_with_res)
big_mean_weight_data_w_r <- big_mean_weight_data_w_r %>% 
  dplyr::rename("mean_weight" = "getSpeciesMeanWeight.proj..nrow.projectEffort_m....")
big_mean_weight_data_w_r$species <- rep(ref@params@species_params$species, 200)


#**** for loop (with resource - no physiological scaling) ==========================
sim <- 1:200

t <- c()
tt <- c()
groj <- c()
growth <- c()
data_list_with_res_no_phys <- list()
mean_weight_list_with_res_no_phys <- list()

for (i in sim) {
  
  t <- params@species_params
  
  t$ea_met <- 0
  t$ea_int <- 0
  t$ea_mor <- 0

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
  
  data_list_with_res_no_phys[[i]] <- growth
  
  mean_weight_list_with_res_no_phys[[i]] <- data.frame(getSpeciesMeanWeight(proj)[nrow(projectEffort_m), ] )
  
}

big_growth_data_w_r_n_p <- dplyr::bind_rows(data_list_with_res_no_phys)
big_mean_weight_data_w_r_n_p <- dplyr::bind_rows(mean_weight_list_with_res_no_phys)
big_mean_weight_data_w_r_n_p <- big_mean_weight_data_w_r_n_p %>% 
  dplyr::rename("mean_weight" = "getSpeciesMeanWeight.proj..nrow.projectEffort_m....")
big_mean_weight_data_w_r_n_p$species <- rep(ref@params@species_params$species, 200)


#**** for loop (no resource) =======================================================
sim <- 1:200

t <- c()
tt <- c()
groj <- c()
growth <- c()
data_list_no_res <- list()
mean_weight_list_no_res <- list()

for (i in sim) {
  
  t <- params@species_params
  
  t$ea_met <- ea$met[i]
  t$ea_int <- ea$int[i]
  t$ea_mor <- ea$mor[i]
  
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
  
  mean_weight_list_no_res[[i]] <- data.frame(getSpeciesMeanWeight(proj)[nrow(projectEffort_m), ] )
  
}

big_growth_data_no_r <- dplyr::bind_rows(data_list_no_res)
big_mean_weight_data_no_r <- dplyr::bind_rows(mean_weight_list_no_res)
big_mean_weight_data_no_r <- big_mean_weight_data_no_r %>% 
  dplyr::rename("mean_weight" = "getSpeciesMeanWeight.proj..nrow.projectEffort_m....")
big_mean_weight_data_no_r$species <- rep(ref@params@species_params$species, 200)


# C. PLOT ==========================================================================
#** Plot growth rates ==============================================================
big_growth_data_w_r$scen <- "Physio. + Resource" # "With resource temp. dep."
big_growth_data_no_r$scen <- "Physio." # "No resource temp. dep."
big_growth_data_w_r_n_p$scen <- "Resource" # "With resource temp. dep. no. phys"
#big_growth_data_w_r_barn$scen <- "Physio. + Resource (obs.)" # "With resource temp. dep. barn"
#big_growth_data_w_r_n_p_barn$scen <- "Resource (obs.)" # "With resource temp. dep. barn. no. phys. "

big_growth_data_w_r$sim <- paste("wr", big_growth_data_w_r$sim, sep = "")
big_growth_data_no_r$sim <- paste("nr", big_growth_data_no_r$sim, sep = "")
big_growth_data_w_r_n_p$sim <- paste("wr_np", big_growth_data_w_r_n_p$sim, sep = "")
#big_growth_data_w_r_barn$sim <- paste("wr_bar", big_growth_data_w_r_barn$sim, sep = "")
#big_growth_data_w_r_n_p_barn$sim <- paste("wr_np_bar", big_growth_data_w_r_n_p_barn$sim, sep = "")

big_growth_data <- rbind(big_growth_data_w_r, 
                         big_growth_data_no_r, 
                         big_growth_data_w_r_n_p#,
                         #big_growth_data_w_r_barn,
                         #big_growth_data_w_r_n_p_barn
                         )

head(big_growth_data)
str(big_growth_data)

big_growth_data$scen <- as.factor(big_growth_data$scen)

# Calculate means, max and min for plotting
mean_dat <- data.frame(big_growth_data %>%
  dplyr::group_by(Species, Age, scen) %>% 
  dplyr::summarise(mean_val = mean(value),
                   max_val = max(value),
                   min_val = min(value),
                   val025 = quantile(value, probs = 0.025),
                   val975 = quantile(value, probs = 0.975)))

head(mean_dat)
str(mean_dat)


col <- RColorBrewer::brewer.pal("Dark2", n = 5)

# Reorder factor levels
mean_dat$Species <- factor(mean_dat$Species, levels = c("Sprat", "Herring", "Cod"))
big_growth_data$Species <- factor(big_growth_data$Species, levels = c("Sprat", "Herring", "Cod"))


test <- seq(1:10)

quantile(test, probs = 0.5)

quantile(test, probs = c(0.025, 0.975))

0.975 + 0.025


# Plot growth curves (mean and ribbons)
p1 <- ggplot(mean_dat, 
             # aes(x = Age, ymin = min_val, ymax = max_val, fill = factor(scen)) # min max version
             aes(x = Age, ymin = val025, ymax = val975, fill = factor(scen))  # percentile version
             ) +
  geom_line(data = mean_dat, aes(Age, mean_val, color = factor(scen)),
            inherit.aes = FALSE, size = 0.5) +
  geom_ribbon(alpha = 0.175, color = NA) +  
  labs(y = "Body mass [g]") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = rev(col)) +
  scale_fill_manual(values = rev(col)) +
  guides(color = FALSE, fill = FALSE) +
  geom_line(data = refGrowth, aes(Age, value), color = "black", 
           inherit.aes = FALSE, size = 0.5, alpha = 0.7, linetype = "dashed") +
  NULL

p1

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 10),
                                       legend.position = "bottom",
                                       legend.title = element_blank())


# Plot growth curves (all curves)
p1b <- big_growth_data %>% filter(Age > 0 & Age < 16) %>%
  ggplot(., aes(x = Age, y = value, color = factor(scen), group = sim)) +
  geom_line(size = 0.3, alpha = 0.05) +
  labs(y = "Body mass [g]") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = rev(col)) +
  theme_classic(base_size = 14) +
  guides(color = FALSE, fill = FALSE) +
  geom_line(data = filter(refGrowth, Age < 16), aes(Age, value), color = "black",
            inherit.aes = FALSE, size = 0.3, alpha = 1, linetype = "dashed") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  NULL


# Plot relative growth curves (ribbons)
rel_dat <- big_growth_data %>%
  dplyr::group_by(Species, Age, scen) %>% 
  dplyr::summarise(mean_val = mean(re_growth),
                   max_val = max(re_growth),
                   min_val = min(re_growth),
                   val025 = quantile(re_growth, probs = 0.025),
                   val975 = quantile(re_growth, probs = 0.975))

rel_dat$Species <- factor(rel_dat$Species, levels = c("Sprat", "Herring", "Cod"))

sort(unique(big_growth_data$scen))

p2 <- rel_dat %>% filter(Age > 0) %>% 
ggplot(., 
       #aes(x = Age, ymin = min_val, ymax = max_val, fill = factor(scen))
       aes(x = Age, ymin = val025, ymax = val975, fill = factor(scen))
       ) +
  geom_line(data = filter(rel_dat, Age > 0), aes(Age, mean_val, color = factor(scen)),
            inherit.aes = FALSE, size = 0.5, alpha = 1) +
  geom_ribbon(alpha = 0.175, color = NA) +  
  labs(y = "Body mass relative to\nconstant temperature") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0),
                     #limits = c(0.95, 2.2)
                     ) + 
  scale_color_manual(values = rev(col),
                     name = "Scenario") +
  scale_fill_manual(values = rev(col)) +
  guides(fill = FALSE,
         colour = guide_legend(nrow = 1,
                               override.aes = list(alpha = 1,
                                                   linetype = 1,
                                                   size = 2))) +
  geom_hline(yintercept = 1, size = 0.3, linetype = "dashed", color = "black") +
  NULL

pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 10),
                                       legend.position = "bottom")

# Plot relative growth curves (all curves)
p2b <- big_growth_data %>% filter(Age > 0 & Age < 16) %>%
  ggplot(., aes(x = Age, y = re_growth, color = factor(scen), group = sim)) +
  geom_line(alpha = 0.05, size = 1) +
  geom_line(data = filter(rel_dat, Age > 0 & Age < 16), aes(Age, mean_val, color = factor(scen)),
            inherit.aes = FALSE, size = 1.75, linetype = 2, alpha = 1) +
  labs(y = "Body mass relative to\nconstant temperature") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)
                     #, limits = c(0.95, 2.2)
                     ) +
  scale_color_manual(values = rev(col),
                     name = "Scenario") +
  guides(colour = guide_legend(nrow = 3,
                               override.aes = list(alpha = 1,
                                                   linetype = 1))) +
  guides(linetype = FALSE) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom"#,
        #legend.direction = "vertical"
        ) +
  geom_hline(yintercept = 1, size = 0.3, linetype = "dashed", color = "black") +
  NULL

p2b

# Plot together
pWord1 / pWord2

ggsave("baltic/figures/growth_project.png", width = 6.5, height = 6.5, dpi = 600)


#** Plot mean weights by species ===================================================
big_mean_weight_data_w_r$scen <- "Physio. + Resource" # "With resource temp. dep."
big_mean_weight_data_no_r$scen <- "Physio." # "No resource temp. dep."
big_mean_weight_data_w_r_n_p$scen <- "Resource" # "With resource temp. dep. no. phys"
#big_mean_weight_data_w_r_barn$scen <- "Physio. + Resource (obs.)" # "With resource temp. dep. barn"
#big_mean_weight_data_w_r_n_p_barn$scen <- "Resource (obs.)" # "With resource temp. dep. barn. no. phys. "

big_mean_weight_data <- rbind(big_mean_weight_data_w_r, 
                              big_mean_weight_data_no_r,
                              big_mean_weight_data_w_r_n_p#,
                              #big_mean_weight_data_w_r_barn,
                              #big_mean_weight_data_w_r_n_p_barn
                              )

big_mean_weight_data$scen <- as.factor(big_mean_weight_data$scen)

# Reorder factor levels
big_mean_weight_data$species <- factor(big_mean_weight_data$species, levels = c("Sprat", "Herring", "Cod"))

# Define palette
col <- rev(RColorBrewer::brewer.pal("Dark2", n = 5))

# Get reference mean weight
ref_w <- data.frame(species = c("Cod", "Sprat", "Herring"),
                    mean_weight = refSpeciesMeanWeight)

# Plot
p3 <- ggplot(big_mean_weight_data, aes(x = scen, y = mean_weight, fill = scen, colour = scen)) +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  #coord_flip() +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col, 
                    name = "Scenario") +
  # this doesn't work with facet_wrap's free scales...
  # geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = FALSE, alpha = 0.7) +
  geom_point(position = position_jitter(width = .15), size = 1.1, alpha = 0.9, shape = 21, color = "white") +
  geom_boxplot(aes(x = scen, y = mean_weight),
               outlier.shape = NA, alpha = 0.1, width = .2, color = "black", size = 0.5) +
  guides(fill = FALSE) +
  labs(x = "", y = "Mean weight [g]") +
  geom_hline(data = ref_w, aes(yintercept = mean_weight), linetype = 2) +
  guides(fill = guide_legend(#nrow = 3,
                             override.aes = list(alpha = 0.8))) +
  NULL

pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 10),
                                       axis.text.x = element_blank(),
                                       legend.direction = "vertical",
                                       legend.position = "bottom",
                                       legend.justification = c(1, 0),
                                       aspect.ratio = 3/4)

pWord3

ggsave("baltic/figures/mean_weight.png", width = 6.5, height = 6.5, dpi = 600)

# ggplot(big_mean_weight_data, aes(x = scen, y = mean_weight, fill = scen, colour = scen)) +
#   facet_wrap(~ species, scales = "free", nrow = 3) +
#   coord_flip() +
#   scale_color_manual(values = col) +
#   scale_fill_manual(values = col, name = "Scenario") +
#   geom_point(position = position_jitter(width = .15), size = 1.1, alpha = 0.7, shape = 21, color = "white") +
#   geom_boxplot(aes(x = scen, y = mean_weight),
#                outlier.shape = NA, alpha = 0.2, width = .2, color = "black", size = 0.5) +
#   guides(colour = FALSE) +
#   theme_classic(base_size = 14) +
#   theme(aspect.ratio = 3/4) +
#   labs(y = "", x = "Mean weight [g") +
#   geom_hline(data = ref_w, aes(yintercept = mean_weight), linetype = 2) +
#   guides(fill = guide_legend(override.aes = list(alpha = 0.8))) +
#   #theme(legend.position = "bottom") +
#   NULL


# D. Some tests ====================================================================
# Testing relative growth rates are correct (seems like larger difference there than in actual growth rates)
big_growth_data %>%
  filter(Age > 0 & Age < 16 & Species == "Cod") %>%
  ggplot(., aes(x = Age, y = value, color = factor(scen), group = sim)) +
  geom_line(size = 0.5, alpha = 0.03) +
  labs(y = "Body mass (g)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = rev(col)) +
  theme_classic(base_size = 14) +
  guides(color = FALSE, fill = FALSE) +
  geom_line(data = filter(refGrowth, Age < 16 & Species == "Cod"), aes(Age, value), color = "black",
            inherit.aes = FALSE, size = 0.5, alpha = 0.7, linetype = "dashed") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  ylim(400, 1100) +
  xlim(3, 4) +
  NULL


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

str(tt)

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

# The activation energy must be much lower for Cmax for it to be negative...

# Testing I can reproduce that by changing h directly, and not via temperature
t <- params@species_params

exp(0.8 * (((273.15+11) - (273.15+10)) / ((8.617332e-05) * (273.15+11) * (273.15+10))))
exp(0.3 * (((273.15+11) - (273.15+10)) / ((8.617332e-05) * (273.15+11) * (273.15+10))))

t$ks <- params@species_params$ks * 1.122307 # Or check what the scalar should be...
t$h <- params@species_params$h * 1.04422 # Or check what the scalar should be...

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
                effort = projectEffort_m[137, ],
                temperature = rep(t_ref, 200),
                diet_steps = 10,
                t_max = 200)   

growth <- getGrowth(proj)

# Now do the no-scaling scenario
t <- params@species_params

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
                effort = projectEffort_m[137, ],
                temperature = rep(t_ref+1.5, 200),
                diet_steps = 10,
                t_max = 200)    

growth2 <- getGrowth(proj)

growth$scen <- "bad"
growth2$scen <- "good"

kek <- rbind(refGrowth2, growth)

ggplot(kek, aes(Age, value, color = scen)) +
  geom_line() +
  facet_wrap(~ Species, scales = "free")
