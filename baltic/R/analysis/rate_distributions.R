#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.04.16: Max Lindmark
#
# Code for sampling from and plotting the distributions that we sample activation 
# energies from (1 for each rate that is temperature-dependent). The means and stds
# come from Lindmark et al (in prep), literature/metabolic theory.
#  
# A. Load parameters and set up objects
#
# B. Sample and save as .csv for analysis...
#
# C. Plot distributions
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(ggplot2)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Print package versions
# print(sessionInfo())
# other attached packages:
# [1] mizer_1.1 testthat_2.0.0 patchwork_0.0.1 dplyr_0.8.1 tidyr_0.8.3       
# [6] viridis_0.5.1 viridisLite_0.3.0 magrittr_1.5 RCurl_1.95-4.12 bitops_1.0-6      
# [11] RColorBrewer_1.1-2 usethis_1.4.0 devtools_2.0.2 ggplot2_3.1.1  

# Define function for shading densities
funcShaded <- function(x, lower_bound) {
  y = dnorm(x, mean = m, sd = std)
  y[x < lower_bound] <- NA
  return(y)
}


# B. SAMPLE DISTRIBUTIONS ==========================================================
set.seed(4361)

# We'll take 100 samples for now.
n <- 200

# Metabolic rate (Lindmark et al (in prep))
met_u <- 0.57
met_sd <- 0.018
met_q <- qnorm(c(0.025, 0.975), met_u, met_sd)
met <- rnorm(n = n, met_u, met_sd)

# Background Mortality (assumed to be = metabolic rate)
mor_u <- 0.57
mor_sd <- 0.018
mor <- rnorm(n = n, mor_u, mor_sd)

# Maximum consumption rate (Lindmark et al (in prep))
int_u <- 0.52
int_sd <- 0.03
int_q <- qnorm(c(0.025, 0.975), int_u, int_sd)
int <- rnorm(n = n, int_u, int_sd)

# Resource turnover rate (# van Savage et al (2004) AmNat. 95% of area within 0.56-1.14)
gro_u <- 0.8
gro_sd <- 0.13
gro_q <- qnorm(c(0.025, 0.975), gro_u, gro_sd)
gro <- rnorm(n = n, gro_u, gro_sd)

# Carrying capacity
car <- runif(n = n, -0.83, 0)

# Put them all together
ea <- data.frame(met, mor, int, gro, car)
ea_dat <- pivot_longer(ea, 1:5, names_to = "rate", values_to = "activation_energy")

# Plot samples
ggplot(ea_dat, aes(activation_energy)) + 
  facet_wrap(~rate, scales = "free") +
  geom_histogram() +
  theme_classic()

write.csv(ea, "baltic/params/samples_activation_energy.csv")

head(ea_dat)


# C. PLOT DISTRIBUTIONS ============================================================
#**** Metabolism ==================================================================
m <- met_u
std <- met_sd

meta <- ggplot(data.frame(x = c(0.5, 0.64)), aes(x)) + 
  stat_function(fun = dnorm, geom = "area",  fill = "grey30", 
                args = list(mean = m, sd = std)) +
  stat_function(fun = funcShaded, args = list(lower_bound = met_q[1]), 
                geom = "area", fill = "grey60") +
  stat_function(fun = funcShaded, args = list(lower_bound = met_q[2]), 
                geom = "area", fill = "grey30") +
  coord_cartesian(expand = 0) +
  ggtitle("Metabolic Rate") +
  labs(x = "Activation Energy") +
  theme_classic() +
  NULL

meta


#**** Mortality ==================================================================
# We assume this follows metabolism
m <- mor_u
std <- mor_sd

mort <- ggplot(data.frame(x = c(0.2, 1)), aes(x)) + 
  stat_function(fun = dnorm, geom = "area",  fill = "grey30", 
                args = list(mean = m, sd = std)) +
  stat_function(fun = funcShaded, args = list(lower_bound = mor_q[1]), 
                geom = "area", fill = "grey60") +
  stat_function(fun = funcShaded, args = list(lower_bound = mor_q[2]), 
                geom = "area", fill = "grey30") +
  coord_cartesian(expand = 0) +
  ggtitle("Metabolic Rate") +
  labs(x = "Activation Energy") +
  theme_classic() +
  NULL

mort


#**** Maximum consumption rate =====================================================
# Lindmark et al (in prep)
m <- int_u
std <- int_sd

inta <- ggplot(data.frame(x = c(0.2, 1)), aes(x)) + 
  stat_function(fun = dnorm, geom = "area",  fill = "grey30", 
                args = list(mean = m, sd = std)) +
  stat_function(fun = funcShaded, args = list(lower_bound = int_q[1]), 
                geom = "area", fill = "grey60") +
  stat_function(fun = funcShaded, args = list(lower_bound = int_q[2]), 
                geom = "area", fill = "grey30") +
  coord_cartesian(expand = 0) +
  ggtitle("Metabolic Rate") +
  labs(x = "Activation Energy") +
  theme_classic() +
  NULL

inta


#**** Resource turnover rate =====================================================
m <- gro_u
std <- gro_sd

grow <- ggplot(data.frame(x = c(0.2, 1)), aes(x)) + 
  stat_function(fun = dnorm, geom = "area",  fill = "grey30", 
                args = list(mean = m, sd = std)) +
  stat_function(fun = funcShaded, args = list(lower_bound = gro_q[1]), 
                geom = "area", fill = "grey60") +
  stat_function(fun = funcShaded, args = list(lower_bound = gro_q[2]), 
                geom = "area", fill = "grey30") +
  coord_cartesian(expand = 0) +
  ggtitle("Metabolic Rate") +
  labs(x = "Activation Energy") +
  theme_classic() +
  NULL

grow


