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


#**** for loop ====================================================================
sim <- 1:100

t <- c()
tt <- c()
groj <- c()
growth <- c()
data_list_g <- list()

for (i in sim) {
  
  t <- params@species_params
  
  t$ea_met <- ea$met[i]
  t$ea_int <- ea$int[i]
  t$ea_mor <- ea$mor[i]
  
  t$ca_int <- -0.004 # Here we just use the fixed values
  t$ca_met <- 0.001 # Here we just use the fixed values
  
  tt <- MizerParams(t, 
                    ea_gro = ea$gro[i],
                    ea_car = ea$car[i],
                    kappa_ben = kappa_ben,
                    kappa = kappa,
                    w_bb_cutoff = w_bb_cutoff,
                    w_pp_cutoff = w_pp_cutoff,
                    r_pp = r_pp,
                    r_bb = r_bb)
  
  proj <- project(tt, 
                  dt = dt,
                  effort = projectEffort_m,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max,
                  t_ref = 10)   
  
  growth <- getGrowth(proj)
  
  growth$ea_met <- proj@params@species_params$ea_met[1]
  growth$ea_mor <- proj@params@species_params$ea_mor[1]
  growth$ea_int <- proj@params@species_params$ea_int[1]
  
  growth$ea_gro <- proj@params@ea_gro
  growth$ea_car <- proj@params@ea_car
  
  growth$sim <- i
  
  growth$re_growth <- growth$value - refGrowth$value
  
  data_list_g[[i]] <- growth
  
}

str(data_list_g)

big_growth_data <- dplyr::bind_rows(data_list_g)


#**** Plot growth rates ============================================================
blues <- RColorBrewer::brewer.pal(n = 5, "Blues")
reds <- RColorBrewer::brewer.pal(n = 5, "Reds")

mean_dat <- big_growth_data %>%
  dplyr::group_by(Species, Age) %>% 
  dplyr::summarise(mean_val = mean(value))

mean_dat

ggplot(big_growth_data, aes(Age, value, group = factor(sim))) +
  geom_line(size = 1, alpha = 0.05, color = "grey20") +
  labs(y = "Size (g)") +
  facet_wrap(~Species, scales = "free_y", ncol = 3) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_linetype_manual(values = c("solid", "twodash")) + 
  guides(color = FALSE) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  geom_line(data = refGrowth, aes(Age, value), color = "red", inherit.aes = FALSE,
            size = 0.2) +
  #geom_line(data = mean_dat, aes(Age, mean_val), color = "white", inherit.aes = FALSE,
  #          size = 0.2) +
  NULL



