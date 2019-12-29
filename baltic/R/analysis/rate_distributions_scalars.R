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
# C. Plot temperature scalars
#
# D. Plot distributions
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
met_sd <- 0.03
met_q <- qnorm(c(0.025, 0.975), met_u, met_sd)
met <- rnorm(n = n, met_u, met_sd)

# Background Mortality (assumed to be = metabolic rate)
mor_u <- 0.62
mor_sd <- 0.03
mor <- rnorm(n = n, mor_u, mor_sd)

# Maximum consumption rate (Lindmark et al (in prep))
int_u <- 0.67
int_sd <- 0.08
int_q <- qnorm(c(0.025, 0.975), int_u, int_sd)
int <- rnorm(n = n, int_u, int_sd)

# Resource turnover rate (# van Savage et al (2004) AmNat. 95% of area within 0.56-1.14)
gro_u <- 0.8
gro_sd <- 0.13
gro_q <- qnorm(c(0.025, 0.975), gro_u, gro_sd)
gro <- rnorm(n = n, gro_u, gro_sd)

# Carrying capacity (theoretical range, based on Gilbert et al (2014), ELE)
car <- runif(n = n, -0.8, 0)

# Put them all together & rename factor levels for plotting
ea <- data.frame(met, mor, int, gro, car)
ea_dat <- ea %>% 
  pivot_longer(1:5, names_to = "rate", values_to = "activation_energy") %>% 
  mutate(rate = fct_recode(rate, 
                           "Resource\ncarrying capacity" = "car",
                           "Resource growth rate" = "gro",
                           "Maximum\nconsumption rate" = "int",
                           "Metabolic rate" = "met",
                           "Background\nmortality rate" = "mor"))

# Plot samples
# Reorder levels
ea_dat$rate2 <- factor(ea_dat$rate, levels = c("Maximum\nconsumption rate", "Metabolic rate", "Background\nmortality rate", 
                                               "Resource growth rate", "Resource\ncarrying capacity"))

ggplot(ea_dat, aes(activation_energy)) + 
  facet_wrap(~rate2, scales = "free") +
  geom_histogram() +
  theme_classic(base_size = 14) +
  coord_cartesian(expand = 0) +
  labs(y = "Count", x = "Activation energy") +
  NULL

#ggsave("baltic/figures/supp/random_activation_energies.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")

#write.csv(ea, "baltic/params/samples_activation_energy.csv")

head(ea_dat)


# C. PLOT TEMP SCALARS =============================================================
script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/tempFun.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Read in temperature data
temp_datRCP8.5 <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/data/Climate/Test_RCP8.5_from_graph.csv"), sep = ";", stringsAsFactors = FALSE)

# Calculate average temperature by year
tempDat <- data.frame(temp_datRCP8.5 %>% 
                        dplyr::filter(year < 2050.5) %>% 
                        dplyr::mutate(Year.r = factor(round(year, digits = 0))) %>% 
                        dplyr::group_by(Year.r) %>% 
                        dplyr::summarize(mean_temp = mean(rel_temp_70.99Ave)) %>% 
                        dplyr::mutate(Year.num = as.numeric(as.character(Year.r))))

t_ref <- 9.57

tempDat$mean_temp_scaled <- tempDat$mean_temp + t_ref

# For loop through scalars
data_list <- c()
scalDat <- c()
#temp <- seq(10, 12, 0.1)
temp <- tempDat$mean_temp_scaled
index <- seq(1:length(temp))

ea_loop <- ea_dat %>% 
  dplyr::group_by(rate) %>% 
  dplyr::summarize(max_ea = max(activation_energy),
            min_ea = min(activation_energy),
            mean_ea = mean(activation_energy))

df <- data.frame(max_ea = rep(ea_loop$max_ea, each = length(temp)),
                 min_ea = rep(ea_loop$min_ea, each = length(temp)),
                 mean_ea = rep(ea_loop$mean_ea, each = length(temp)),
                 rate = rep(ea_loop$rate, each = length(temp)),
                 temp = temp)

scalMin <- data.frame(scal = as.numeric(tempFun(temperature = df$temp, 
                                                t_ref = 10, 
                                                Ea = df$min_ea, 
                                                c_a = 0, 
                                                w = 10)),
                      temp = df$temp,
                      rate = df$rate,
                      val = "min")

scalMean <- data.frame(scal = as.numeric(tempFun(temperature = df$temp, 
                                                 t_ref = 10, 
                                                 Ea = df$mean_ea, 
                                                 c_a = 0, 
                                                 w = 10)),
                       temp = df$temp,
                       rate = df$rate,
                       val = "mean")

scalMax <- data.frame(scal = as.numeric(tempFun(temperature = df$temp, 
                                                t_ref = 10, 
                                                Ea = df$max_ea, 
                                                c_a = 0, 
                                                w = 10)),
                      temp = df$temp,
                      rate = df$rate,
                      val = "max")

dat <- rbind(scalMin, scalMean, scalMax)

dat2 <- pivot_wider(data = dat, names_from = val, values_from = scal)

# Plot mean min and max
col <- rev(RColorBrewer::brewer.pal("Dark2", n = 5))
dat2$rate <- factor(dat2$rate, levels = c("Maximum\nconsumption rate", "Metabolic rate", "Background\nmortality rate", 
                                          "Resource growth rate", "Resource\ncarrying capacity"))

ggplot(dat2, aes(x = temp, ymin = min, ymax = max, fill = rate)) + 
  geom_ribbon(alpha = 0.8, fill = "gray75") +
  geom_line(data = dat2, aes(temp, mean), color = "black", linetype = "dashed", size = 1) +
  facet_wrap(~rate, scales = "free") +
  theme_classic(base_size = 16) +
  guides(fill = FALSE) +
  labs(x = expression(paste("Temperature [", degree*C, "]")),
       y = "Rate scalars") +
  NULL

#ggsave("baltic/figures/supp/random_rate_scalar.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


# D. PLOT DISTRIBUTIONS ============================================================
#**** Metabolism ===================================================================
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


