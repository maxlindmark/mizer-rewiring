#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.04.16: Max Lindmark
#
# Code for sampling from and plotting the distributions that we sample activation 
# energies from (1 for each rate that is temperature-dependent). The means and stds
# come from Lindmark et al (in prep), literature/metabolic theory.
#  
# A. Load parameters and set up objects
#
# B. Sample, plot and save as .csv for analysis...
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(tidyverse)
library(RCurl)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Print package versions
# print(sessionInfo())
# other attached packages:
# patchwork_0.0.1 RCurl_1.95-4.12 bitops_1.0-6    forcats_0.4.0   stringr_1.4.0   
# dplyr_0.8.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.0     tibble_2.1.3   
# ggplot2_3.2.1   tidyverse_1.2.1

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
met_u <- 0.62
met_sd <- 0.025
met_q <- qnorm(c(0.025, 0.975), met_u, met_sd)
met <- rnorm(n = n, met_u, met_sd)

# Background Mortality (assumed to be = metabolic rate)
mor_u <- 0.62
mor_sd <- 0.025
mor <- rnorm(n = n, mor_u, mor_sd)

# Maximum consumption rate (Lindmark et al (in prep))
int_u <- 0.69
int_sd <- 0.078
int_q <- qnorm(c(0.025, 0.975), int_u, int_sd)
int <- rnorm(n = n, int_u, int_sd)

# Resource turnover rate (# van Savage et al (2004) AmNat. 95% of area within 0.56-1.14)
gro_u <- 0.73
gro_sd <- 0.1
gro_q <- qnorm(c(0.025, 0.975), gro_u, gro_sd)
gro <- rnorm(n = n, gro_u, gro_sd)

# Resource carrying capacity (# Barnes et al  (201) JOURNAL OF PLANKTON RESEARCH. 95% of area within -(-0.99) - 0.59)
# b_car_u <- -0.7953
# b_car_sd <- 0.1
# b_car_q <- qnorm(c(0.025, 0.975), b_car_u, b_car_sd)
# b_car <- rnorm(n = n, b_car_u, b_car_sd)

# Carrying capacity (theoretical range, based on Gilbert et al (2014), ELE)
# NOT USED IN PROJECTION
car <- runif(n = n, -0.8, 0)

# Put them all together & rename factor levels for plotting
#ea <- data.frame(met, mor, int, gro, car) # This is with uniform carrying capacity
ea <- data.frame(met, mor, int, gro, -gro) #, b_car) # This is carrying = - growth

ea_dat <- ea %>% 
  pivot_longer(1:5, names_to = "rate", values_to = "activation_energy") %>% 
  mutate(rate = fct_recode(rate, 
                           #"Resource\ncarrying capacity" = "car",
                           "Resource\ncarrying capacity" = "X.gro",
                           #"Resource\ncarrying capacity (obs)" = "b_car",
                           "Resource growth rate" = "gro",
                           "Maximum\nconsumption rate" = "int",
                           "Metabolic rate" = "met",
                           "Background\nmortality rate" = "mor"))

# Plot samples
# Reorder levels
ea_dat$rate2 <- factor(ea_dat$rate, levels = c("Maximum\nconsumption rate", "Metabolic rate", "Background\nmortality rate", 
                                               "Resource growth rate", "Resource\ncarrying capacity"
                                               #, "Resource\ncarrying capacity (obs)"
                                               ))

p1 <- ggplot(ea_dat, aes(activation_energy)) + 
  facet_wrap(~rate2, scales = "free") +
  geom_histogram() +
  coord_cartesian(expand = 0) +
  labs(y = "Count", x = "Activation energy") +
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12),
                                       aspect.ratio = 1)

pWord1

ggsave("baltic/figures/supp/random_activation_energies.png", width = 6.5, height = 6.5, dpi = 600)

#write.csv(ea, "baltic/params/samples_activation_energy.csv")

