#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. Analyisis of size-spectra & mortality with warming and different fishing scenario
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

# Load function for extracting numbers-at-age
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getSpectra.R", ssl.verifypeer = FALSE)
eval(parse(text = func))

# Load function for extracting mortality-at-weight (predation)
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getMortality.R", ssl.verifypeer = FALSE)
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
t_ref <- (10 + 0.1156161)
kappa_ben <- 1
kappa <- 1
w_bb_cutoff <- 20
w_pp_cutoff <- 1
r_pp <- 4
r_bb <- 4


# B. TEMP-DRIVEN CHANGE IN ABUNDANCE-AT-WEIGHT =====================================
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
consTemp[] <- t_ref

ref <- project(pars_with_res, 
               dt = dt,
               effort = projectEffort_m,
               temperature = consTemp,
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)

refSpect <- getSpectra(ref)


#**** for loop through different fishing effort (with temp dep resource) ============
sim <- seq(0.8, 1.2, 0.1) # Factor for scaling fishing mortality
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
                  t_ref = t_ref)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Extract mortality (predation) - modified from plotM2
  #mort <- getMortality(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n #spect$n / refSpect$n
  
  data_list_with_res[[i]] <- spect
  
}

big_spect_data_w_r <- dplyr::bind_rows(data_list_with_res)


#**** for loop through different fishin effort (with temp dep resource) ============
projectEffort_new <- projectEffort_m
tt <- c()
groj <- c()
spect <- c()
data_list_no_res <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  #plot(y = as.numeric(projectEffort_m[, 1]), x = 1:137)
  #lines(y = as.numeric(projectEffort_new[, 1]), x = 1:137, col = "blue")
  
  proj <- project(pars_no_res, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max,
                  t_ref = t_ref)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Extract mortality (predation) - modified from plotM2
  # mort <- getMortality(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n
  
  # ggplot(spect, aes(log10(w), re_spec, color = factor(Fm))) +
  #   facet_wrap(~species, scales = "free", ncol = 1) + 
  #   geom_hline(yintercept = 1) +
  #   geom_line() +
  #   NULL
  # 
  # spect %>% filter(re_spec > 1 & species == "Cod")
  # 
  # test <- spect %>% filter(species == "Cod" & Fm == 1)
  # tt <- filter(getSpectra(ref), species == "Cod")
  # test2 <- data.frame(n = c(test$n, tt$n),
  #                     w = rep(test$w), 2,
  #                     scen = rep(c("warm", "ref"), each = 100))
  # 
  # ggplot(test2, aes(log10(w), n, color = factor(scen))) +
  #   geom_line(size = 1) +
  #   NULL
  # 
  # test2 %>% filter(log10(w) > 4)
  
  data_list_no_res[[i]] <- spect
  
}

big_spect_data_no_r <- dplyr::bind_rows(data_list_no_res)


#**** for loop through different fishing effort (constant temperature) =============
projectEffort_new <- projectEffort_m
tt <- c()
groj <- c()
spect <- c()
data_list_con_temp <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_with_res, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = consTemp,
                  diet_steps = 10,
                  t_ref = t_ref)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Extract mortality (predation) - modified from plotM2
  # mort <- getMortality(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n
  
  # test <- getSpectra(proj)$w / getSpectra(ref)$w
  # test
  
  data_list_con_temp[[i]] <- spect
  
}

big_spect_data_con_temp <- dplyr::bind_rows(data_list_con_temp)



#**** Plot size-spectra ============================================================
options(scipen = 10000) # Set higher level before using scientific notation over normal

big_spect_data_w_r$scen <- "With resource temp. dep."
big_spect_data_no_r$scen <- "No resource temp. dep."
big_spect_data_con_temp$scen <- "Constant temp."

big_spect_data_w_r$sim <- paste("wr", big_spect_data_w_r$sim, sep = "")
big_spect_data_no_r$sim <- paste("nr", big_spect_data_no_r$sim, sep = "")
big_spect_data_con_temp$sim <- paste("ct", big_spect_data_con_temp$sim, sep = "")

big_spect_data <- rbind(big_spect_data_w_r, big_spect_data_no_r, big_spect_data_con_temp)

big_spect_data$scen <- as.factor(big_spect_data$scen)

# Reorder factor levels
big_spect_data$species <- factor(big_spect_data$species, levels = c("Sprat", "Herring", "Cod"))

# Define colors
pal <- viridis(n = 5, option = "cividis")
pal[3] <- RColorBrewer::brewer.pal(n = 5, "Dark2")[4]

# Plot all together

plotdf <- select(proj@params@species_params, species, w_mat)

# Plot relative spectra
big_spect_data %>% 
  # Note we are doing some filtering here to be able to plot all together. Look at how large 
  # cod abundance increases rapibdly in the low fishing scenarios.
  #filter(n > 0 & re_spec > 0.85 & re_spec < 1.15 & w > 0.001) %>%
  filter(n > 0 & re_spec > 0.8 & re_spec < 1.2 & w > 0.001) %>%
  ggplot(., aes(w, re_spec, color = factor(Fm), group = sim, 
                linetype = factor(Fm), alpha = factor(Fm))) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_line(size = 1) + 
  scale_colour_manual(values = pal) +
  scale_alpha_manual(values = c(0.7, 0.7, 1, 0.7, 0.7)) +
  scale_linetype_manual(values = c("solid", "solid", "longdash", "solid", "solid")) +
  facet_grid(scen ~ species, scales = "free") +
  theme_classic(base_size = 13) +
  #scale_y_log10(breaks = sim) +
  scale_x_log10() +
  geom_vline(data = plotdf, aes(xintercept = w_mat), color = "red", linetype = "dotted") +
  labs(x ="Body mass (g)",
       y = "Relative abundance (warming/no warming)",
       color = "FMSY\nscaling\nfactor", 
       alpha = "FMSY\nscaling\nfactor", 
       linetype = "FMSY\nscaling\nfactor") +
  NULL

#ggsave("baltic/figures/spectra_project.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


# Plot spectra
col <- RColorBrewer::brewer.pal("Dark2", n = 5)
big_spect_data %>% 
  # Note we are doing some filtering here to be able to plot all together. Look at how large 
  # cod abundance increases rapibdly in the low fishing scenarios.
  #filter(n > 0 & re_spec > 0.85 & re_spec < 1.15 & w > 0.001) %>%
  filter(n > 0 & w > 0.001 & Fm == 1) %>%
  ggplot(., aes(w, n, color = factor(scen), group = sim, 
                linetype = factor(scen), alpha = factor(scen))) + 
  geom_line(size = 0.8) + 
  scale_colour_manual(values = rev(col)) +
  scale_alpha_manual(values = c(0.7, 0.7, 1, 0.7, 0.7)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted")) +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  theme_classic(base_size = 13) +
  scale_x_log10() +
  geom_vline(data = plotdf, aes(xintercept = w_mat), color = "grey10", linetype = "dotted") +
  labs(x ="Body mass (g)",
       y = "Abundance",
       color = "FMSY\nscaling\nfactor", 
       alpha = "FMSY\nscaling\nfactor", 
       linetype = "FMSY\nscaling\nfactor") +
  NULL


# More tests..
# What is the SSB here in these runs?
tt <- project(pars_with_res, 
              dt = dt,
              effort = projectEffort_m,
              temperature = projectTemp$temperature,
              diet_steps = 10,
              t_max = t_max,
              t_ref = t_ref)   

zz <- project(pars_no_res, 
              dt = dt,
              effort = projectEffort_m,
              temperature = projectTemp$temperature,
              diet_steps = 10,
              t_max = t_max,
              t_ref = t_ref)   

str(ref@effort)

# no resource - empirical temp - FMSy fishing
tail(getSSB(zz), 2)
# with resource - empirical temp - FMSy fishing
tail(getSSB(tt), 2)
# ref - constant temp - FMSy fishing
tail(getSSB(ref), 2)

# Ok, this fits with the spectra plot... SSB increases. So why does it not do that in the yield-script?

big_spect_data %>% filter(Fm == 1 & n > 0.0001) %>% 
  ggplot(., aes(log10(w), re_spec, color = scen)) + 
  geom_line(size = 0.8) +
  theme_classic() +
  scale_x_log10() +
  facet_wrap(~species, scales = "free", nrow = 3) + 
  NULL

big_spect_data %>% filter(re_spec > 1 & species == "Herring" & scen == "With resource temp. dep." & Fm == 1)


# #**** Plot predation mortality =====================================================
# Compare the predation mortalities from the two extremes:
# 
# projectEffort_new <- projectEffort_m
# projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * 1.2
# 
# proj <- project(pars_with_res, 
#                 dt = dt,
#                 effort = projectEffort_new,
#                 temperature = projectTemp$temperature,
#                 diet_steps = 10,
#                 t_ref = t_ref)   
# 
# warmMort <- getMortality(proj)
# warmMort$scen <- "warm"
# 
# refs <- project(pars_with_res, 
#                 dt = dt,
#                 effort = projectEffort_m,
#                 temperature = consTemp,
#                 diet_steps = 10,
#                 t_ref = t_ref)   
# 
# consMort <- getMortality(refs)
# consMort$scen <- "constant"
# 
# mortdat <- rbind(consMort, warmMort)
# 
# # More of less identical predation mortality
# ggplot(mortdat, aes(w, value, color = scen, linetype = scen)) + 
#   geom_line() +
#   scale_x_log10() +
#   NULL

# plotM2(ref)
# plotM2(proj)

