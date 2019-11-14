#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. Analyisis of size-spectra with warming and different fishing scenario
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


#**** Project reference scenario ===================================================
consTemp <- projectTemp$temperature
t_ref <- (10 + 0.1156161)
consTemp[] <- t_ref

ref <- project(pars_with_res, 
               dt = dt,
               effort = projectEffort_m,
               temperature = consTemp,
               diet_steps = 10,
               t_max = t_max,
               t_ref = (10 + 0.1156161))   

refSpect <- getSpectra(ref)


#**** for loop through different fishin effort (with temp dep resource) ============
sim <- seq(0.6, 1.4, 0.2) # Factor for scaling fishing mortality
iter <- seq(from = 1, to = length(sim))

projectEffort_new <- projectEffort_m
tt <- c()
groj <- c()
spect <- c()
data_list_with_res <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_with_res, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max,
                  t_ref = (10 + 0.1156161))   
  
  spect <- getSpectra(proj)

  # Add in iteration
  spect$sim <- i
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  spect$re_spec <- spect$n / refSpect$n
  
  data_list_with_res[[i]] <- spect
  
}

big_spect_data_w_r <- dplyr::bind_rows(data_list_with_res)

### TEST
unique(big_spect_data_w_r$Fm)

big_spect_data_w_r$species <- factor(big_spect_data_w_r$species, levels = c("Sprat", "Herring", "Cod"))
options(scipen = 10000) # Set higher level before using scientific notation over normal

# Plotting relative spectra

pal <- viridis(n = 5, option = "cividis")

d2 <- RColorBrewer::brewer.pal("Dark2", n = 5)

#pal[3] <- "red"
pal[3] <- d2[4]

big_spect_data_w_r %>% 
  #mutate(Fvar = ifelse(Fm == 1, "n", "y")) %>% 
  #filter(n > 0 & Fvar == "y" & re_spec > 0.75 & re_spec < 2 & w > 0.01) %>% 
  filter(n > 0 & re_spec > 0.75 & re_spec < 2 & w > 0.01) %>% 
  ggplot(., aes(w, re_spec, color = factor(Fm), group = sim, 
                linetype = factor(Fm), alpha = factor(Fm))) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 1, alpha = 0.6) +
  geom_line(size = 1) + 
  #scale_colour_viridis(option = "cividis", discrete = TRUE) +
  scale_colour_manual(values = pal) +
  scale_alpha_manual(values = c(0.7, 0.7, 1, 0.7, 0.7)) +
  scale_linetype_manual(values = c("solid", "solid", "longdash", "solid", "solid")) +
  facet_wrap(~ species, scales = "free", ncol = 3) + 
  theme_classic() +
  scale_y_log10(breaks = sim) +
  scale_x_log10() +
  NULL

# Plotting spectra
# big_spect_data_w_r %>% 
#   filter(n > 0 & Fm %in% c(0.5, 0.75, 1, 1.25, 1.5) & re_spec > 0.75 & re_spec < 10 & w > 0.01) %>% 
#   filter(n > 0.0001 & Fm %in% c(0.5, 0.75, 1, 1.25, 1.5) & re_spec > 0.75 & re_spec < 10 & w > 0.01) %>% 
#   ggplot(., aes(w, n, color = factor(Fm), group = sim)) + 
#   geom_line(alpha = 0.8, size = 1) + 
#   scale_colour_viridis(option = "cividis", discrete = TRUE) +
#   facet_wrap(~ species, scales = "free", ncol = 3) + 
#   theme_classic() +
#   scale_y_log10() +
#   scale_x_log10() +
#   NULL

plotSpectra(proj)


#**** for loop through different fishin effort (with temp dep resource) ============
projectEffort_new <- projectEffort_m
tt <- c()
groj <- c()
spect <- c()
data_list_no_res <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_no_res, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max,
                  t_ref = (10 + 0.1156161))   
  
  spect <- getSpectra(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  spect$re_spec <- spect$n / refSpect$n
  
  data_list_no_res[[i]] <- spect
  
}

big_spect_data_no_r <- dplyr::bind_rows(data_list_no_res)


#### PLOT BELOW WITH FACET GRID


#**** Plot growth rates ============================================================
big_spect_data_w_r$scen <- "With resource temp. dep."
big_spect_data_no_r$scen <- "No resource temp. dep."

big_spect_data_w_r$sim <- paste("wr", big_spect_data_w_r$sim, sep = "")
big_spect_data_no_r$sim <- paste("nr", big_spect_data_no_r$sim, sep = "")

big_spect_data <- rbind(big_spect_data_w_r, big_spect_data_no_r)

big_spect_data$scen <- as.factor(big_spect_data$scen)

# Reorder factor levels
big_spect_data$species <- factor(big_spect_data$species, levels = c("Sprat", "Herring", "Cod"))

# Plot all together
big_spect_data %>% 
  # Note we are doing some filtering here to be able to plot all together. Look at how large 
  # cod abundance increases rapibdly in the low fishing scenarios.
  filter(n > 0 & re_spec > 0.75 & re_spec < 1.3 & w > 0.01) %>% 
  ggplot(., aes(w, re_spec, color = factor(Fm), group = sim, 
                linetype = factor(Fm), alpha = factor(Fm))) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_line(size = 1) + 
  scale_colour_manual(values = pal) +
  scale_alpha_manual(values = c(0.7, 0.7, 1, 0.7, 0.7)) +
  scale_linetype_manual(values = c("solid", "solid", "longdash", "solid", "solid")) +
  facet_grid(scen ~ species, scales = "free") + 
  theme_classic(base_size = 14) +
  scale_y_log10(breaks = sim) +
  scale_x_log10() +
  labs(x ="Weight (g)",
       y = "Relative abundance (warming/no warming)",
       color = "Fishing\nmortality", 
       alpha = "Fishing\nmortality", 
       linetype = "Fishing\nmortality") +
  NULL

#ggsave("baltic/figures/spectra_project.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")
