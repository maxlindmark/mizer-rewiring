# 2019.11.16: Max Lindmark
#
# Code for analyzing the Baltic Sea mizer model. The params-object is saved in the
# calibration_v1 code. 
# 
# A. Load libraries and read in data and parameters
#
# B. Analyisis of size-spectra & mortality with warming and different fishing scenario
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
devtools::load_all(".")

# Print package versions
# print(sessionInfo())
# other attached packages:
# mizer_1.1          testthat_2.3.0     patchwork_0.0.1    dplyr_0.8.3        tidyr_1.0.0        
# viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       RCurl_1.95-4.12   
# bitops_1.0-6       RColorBrewer_1.1-2 devtools_2.2.1     usethis_1.5.1      ggplot2_3.2.1 

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


#**** Update species params ========================================================
t <- params@species_params

t$ea_met <- mean(ea$met)
t$ea_int <- mean(ea$int)
t$ea_mor <- mean(ea$mor)

#t$ca_int <- -0.004 # Here we just use the fixed values
#t$ca_met <- 0.001 # Here we just use the fixed values

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

pars_with_res_barnes <- MizerParams(t, 
                                    ea_gro = 0,
                                    ea_car = mean(ea$b_car),
                                    kappa_ben = kappa_ben,
                                    kappa = kappa,
                                    w_bb_cutoff = w_bb_cutoff,
                                    w_pp_cutoff = w_pp_cutoff,
                                    r_pp = r_pp,
                                    r_bb = r_bb,
                                    t_ref = t_ref)


# No physiological scaling scenarios:
t_no_ea <- t
t_no_ea$ea_int <- 0
t_no_ea$ea_met <- 0
t_no_ea$ea_mor <- 0

pars_with_res_barnes_np <- MizerParams(t_no_ea, 
                                       ea_gro = 0,
                                       ea_car = mean(ea$b_car),
                                       kappa_ben = kappa_ben,
                                       kappa = kappa,
                                       w_bb_cutoff = w_bb_cutoff,
                                       w_pp_cutoff = w_pp_cutoff,
                                       r_pp = r_pp,
                                       r_bb = r_bb,
                                       t_ref = t_ref)

pars_with_res_np <- MizerParams(t_no_ea, 
                                ea_gro = mean(ea$gro),
                                ea_car = mean(ea$car),
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


# B. TEMP-DRIVEN CHANGE IN MORTALITY AND SIZE SPECTRA ==============================
# for-loop to take random samples for distributions representing activation energies
# Then compare that to a projection with a constant temperature

#**** Project reference scenario ===================================================
ref <- project(pars_no_res, 
               dt = dt,
               effort = projectEffort_m,
               temperature = consTemp,
               diet_steps = 10,
               t_max = t_max)

refSpect <- getSpectra(ref)
refMort <- getMortality(ref)

refSpect

# Testing mortality is correct
plot(ref)
ggplot(refMort, aes(w, value, linetype = Species)) + 
  geom_line(size = 1) + 
  theme_bw(base_size = 13) +
  scale_x_log10() +
  theme(legend.position = "bottom",
        aspect.ratio = 1/2) +
  NULL

# Create a little function to get the same structure as mortality... Clean up later
FL_df <- function(object){
  
  vec <- getFeedingLevel(object)
  str(vec)
  
  FL <- data.frame(w = rep(as.numeric(names(vec[137, 1, 1:100])), 3), 
                   value = c(as.numeric(vec[137, 1, 1:100]),
                             as.numeric(vec[137, 2, 1:100]),
                             as.numeric(vec[137, 3, 1:100])),
                   Species = rep(names(vec[137, 1:3, 1]), each = 100))
  
  FL <- arrange(FL, Species, w)
}

refFL <- FL_df(ref)

# plot(refFL$w, refMort$w); abline(0, 1, col = "red")
# unique(refFL$Species) 
# unique(refMort$Species)
# head(refFL)
# head(refMort)


#**** for loop through different fishing effort (with temp dep resource) ============
sim <- seq(0.8, 1.2, 0.1) # Factor for scaling fishing mortality
iter <- seq(from = 1, to = length(sim))

tt <- c()
groj <- c()
spect <- c()
data_list_with_res <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {

  projectEffort_new <- projectEffort_m
    
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_with_res, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Relative spectra
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n #spect$n / refSpect$n
  
  # Extract mortality (predation) - modified from plotM2
  mort <- getMortality(proj)
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in natural mortality
  spect$mort <- mort$value
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  # Relative mortality
  spect$re_mort <- mort$value / refMort$value #spect$n / refSpect$n
  
  # Add feeding level
  fl <- FL_df(proj)
  spect$feedingLevel <- fl$value
  
  # Relative feeding level
  spect$re_feedingLevel <- spect$feedingLevel / refFL$value 
  
  data_list_with_res[[i]] <- spect
  
}

big_spect_data_w_r <- dplyr::bind_rows(data_list_with_res)


#**** for loop through different fishing effort (with temp dep resource BARNES) ============
sim <- seq(0.8, 1.2, 0.1) # Factor for scaling fishing mortality
iter <- seq(from = 1, to = length(sim))

tt <- c()
groj <- c()
spect <- c()
data_list_with_res_barnes <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  projectEffort_new <- projectEffort_m
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_with_res_barnes, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Relative spectra
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n #spect$n / refSpect$n
  
  # Extract mortality (predation) - modified from plotM2
  mort <- getMortality(proj)
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in natural mortality
  spect$mort <- mort$value
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  # Relative mortality
  spect$re_mort <- mort$value / refMort$value #spect$n / refSpect$n
  
  # Add feeding level
  fl <- FL_df(proj)
  spect$feedingLevel <- fl$value
  
  # Relative feeding level
  spect$re_feedingLevel <- spect$feedingLevel / refFL$value 
  
  data_list_with_res_barnes[[i]] <- spect
  
}

big_spect_data_w_r_b <- dplyr::bind_rows(data_list_with_res_barnes)


#**** for loop through different fishing effort (with temp dep resource - NO PHYS) =
sim <- seq(0.8, 1.2, 0.1) # Factor for scaling fishing mortality
iter <- seq(from = 1, to = length(sim))

tt <- c()
groj <- c()
spect <- c()
data_list_with_res_np <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  projectEffort_new <- projectEffort_m
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_with_res_np, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Relative spectra
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n #spect$n / refSpect$n
  
  # Extract mortality (predation) - modified from plotM2
  mort <- getMortality(proj)
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in natural mortality
  spect$mort <- mort$value
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  # Relative mortality
  spect$re_mort <- mort$value / refMort$value #spect$n / refSpect$n
  
  # Add feeding level
  fl <- FL_df(proj)
  spect$feedingLevel <- fl$value
  
  # Relative feeding level
  spect$re_feedingLevel <- spect$feedingLevel / refFL$value 
  
  data_list_with_res_np[[i]] <- spect
  
}

big_spect_data_w_r_np <- dplyr::bind_rows(data_list_with_res_np)


#**** for loop through different fishing effort (with temp dep resource BARNES - NO PHYS) 
sim <- seq(0.8, 1.2, 0.1) # Factor for scaling fishing mortality
iter <- seq(from = 1, to = length(sim))

tt <- c()
groj <- c()
spect <- c()
data_list_with_res_barnes_np <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {
  
  projectEffort_new <- projectEffort_m
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_with_res_barnes_np, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Relative spectra
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n #spect$n / refSpect$n
  
  # Extract mortality (predation) - modified from plotM2
  mort <- getMortality(proj)
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in natural mortality
  spect$mort <- mort$value
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  # Relative mortality
  spect$re_mort <- mort$value / refMort$value #spect$n / refSpect$n
  
  # Add feeding level
  fl <- FL_df(proj)
  spect$feedingLevel <- fl$value
  
  # Relative feeding level
  spect$re_feedingLevel <- spect$feedingLevel / refFL$value 
  
  data_list_with_res_barnes_np[[i]] <- spect
  
}

big_spect_data_w_r_b_np <- dplyr::bind_rows(data_list_with_res_barnes_np)


#**** for loop through different fishing effort (no temp dep resource) ============
tt <- c()
groj <- c()
spect <- c()
data_list_no_res <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {

  projectEffort_new <- projectEffort_m
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  #plot(y = as.numeric(projectEffort_m[, 1]), x = 1:137)
  #lines(y = as.numeric(projectEffort_new[, 1]), x = 1:137, col = "blue")
  
  proj <- project(pars_no_res, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Relative spectra
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n #spect$n / refSpect$n
  
  # Extract mortality (predation) - modified from plotM2
  mort <- getMortality(proj)
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in natural mortality
  spect$mort <- mort$value
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  # Relative mortality
  spect$re_mort <- mort$value / refMort$value #spect$n / refSpect$n
  
  # Add feeding level
  fl <- FL_df(proj)
  spect$feedingLevel <- fl$value
  
  # Relative feeding level
  spect$re_feedingLevel <- spect$feedingLevel / refFL$value 

  data_list_no_res[[i]] <- spect
  
}

big_spect_data_no_r <- dplyr::bind_rows(data_list_no_res)


#**** for loop through different fishing effort (constant temperature) =============
tt <- c()
groj <- c()
spect <- c()
data_list_con_temp <- list()

# The projected fishing mortality starts at row 100

for (i in iter) {

  projectEffort_new <- projectEffort_m
  
  projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * sim[i]
  
  proj <- project(pars_with_res, 
                  dt = dt,
                  effort = projectEffort_new,
                  temperature = consTemp,
                  diet_steps = 10)   
  
  # Apply getSpectra function to get abundance at size
  spect <- getSpectra(proj)
  
  # Add in iteration
  spect$sim <- i
  
  # Relative spectra
  spect$re_spec <- getSpectra(proj)$n / getSpectra(ref)$n #spect$n / refSpect$n
  
  # Extract mortality (predation) - modified from plotM2
  mort <- getMortality(proj)
  
  # Add in fishing mortality
  spect$Fm <- rep(as.numeric(proj@effort[dim(proj@effort)[1], ]), each = nrow(spect)/length(unique(spect$species)))
  
  # Add in natural mortality
  spect$mort <- mort$value
  
  # Add in fishing mortality scalar
  spect$Fm <- sim[i]
  
  # Relative mortality
  spect$re_mort <- mort$value / refMort$value #spect$n / refSpect$n
  
  # Add feeding level
  fl <- FL_df(proj)
  spect$feedingLevel <- fl$value
  
  # Relative feeding level
  spect$re_feedingLevel <- spect$feedingLevel / refFL$value 

  data_list_con_temp[[i]] <- spect
  
}

big_spect_data_con_temp <- dplyr::bind_rows(data_list_con_temp)



# C. PLOT ==========================================================================
#**** Arrange data  ================================================================
options(scipen = 10000) # Set higher level before using scientific notation over normal

big_spect_data_w_r$scen <- "Physio. + Resource (exp.)"
big_spect_data_w_r_b$scen <- "Physio. + Resource (obs.)"
big_spect_data_w_r_np$scen <- "Resource (exp.)"
big_spect_data_w_r_b_np$scen <- "Resource (obs.)"
big_spect_data_no_r$scen <- "Physio."
big_spect_data_con_temp$scen <- "No warming"

big_spect_data_w_r$sim <- paste("wr", big_spect_data_w_r$sim, sep = "")
big_spect_data_w_r_b$sim <- paste("wr_b", big_spect_data_w_r_b$sim, sep = "")
big_spect_data_w_r_np$sim <- paste("wr_np", big_spect_data_w_r_np$sim, sep = "")
big_spect_data_w_r_b_np$sim <- paste("wr_b_np", big_spect_data_w_r_b_np$sim, sep = "")
big_spect_data_no_r$sim <- paste("nr", big_spect_data_no_r$sim, sep = "")
big_spect_data_con_temp$sim <- paste("ct", big_spect_data_con_temp$sim, sep = "")

big_spect_data <- rbind(big_spect_data_w_r, 
                        big_spect_data_w_r_b,
                        big_spect_data_w_r_np, 
                        big_spect_data_w_r_b_np,
                        big_spect_data_no_r, 
                        big_spect_data_con_temp)

big_spect_data$scen <- as.factor(big_spect_data$scen)

# Reorder factor levels
big_spect_data$species <- factor(big_spect_data$species, levels = c("Sprat", "Herring", "Cod"))


#**** Plot size-spectra ============================================================
# Define colors
#pal <- viridis(n = 5, option = "cividis")
#pal[3] <- RColorBrewer::brewer.pal(n = 5, "Dark2")[4]
pal <- RColorBrewer::brewer.pal(n = 5, "Dark2")

# Plot all together
plotdf <- select(proj@params@species_params, species, w_mat)

# Plot relative spectra
big_spect_data %>% 
  # Note we are doing some filtering here to be able to plot all together. Look at how large 
  # cod abundance increases rapibdly in the low fishing scenarios.
  #filter(n > 0 & re_spec > 0.85 & re_spec < 1.15 & w > 0.001) %>%
  #filter(n > 0 & re_spec > 0.8 & re_spec < 1.2 & w > 0.001) %>%
  filter(n > 0 & Fm %in% c(0.9, 1, 1.1) & w > 0.001 & re_spec > 0.8 & re_spec < 1.3 & w > 0.1) %>%
  ggplot(., aes(w, re_spec, color = factor(Fm), group = sim)) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_line(size = 1) + 
  scale_colour_manual(values = pal) +
  #scale_alpha_manual(values = c(0.7, 0.7, 1, 0.7, 0.7)) +
  #scale_linetype_manual(values = c("solid", "solid", "longdash", "solid", "solid")) +
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
  theme(strip.text.y = element_text(size = 5)) +
  NULL

#ggsave("baltic/figures/supp/spectra_FM_project.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


# No fishing and absolute spectra
pal <- rev(RColorBrewer::brewer.pal(n = 5, "Dark2"))
pal2 <- c("black", pal)

unique(big_spect_data$scen)

p1 <- big_spect_data %>% 
  filter(n > 0 & Fm == 1 & w > 0.1) %>%
  ggplot(., aes(w, (n*240.342), color = factor(scen), linetype = scen)) + 
  geom_line(size = 1) + 
  scale_colour_manual(values = pal2,
                      name = "Scenario") +
  theme_classic(base_size = 13) +
  scale_x_log10() +
  guides(linetype = FALSE) +
  scale_y_log10() +
  theme(legend.text = element_text(size = 8)) +
  facet_wrap(~species, scales = "free", nrow = 3) +
  labs(x ="Body mass (g)",
       y = "Abundance") +
  NULL

p1


# No fishing and relative spectra
pal <- RColorBrewer::brewer.pal(n = 5, "Dark2")

unique(big_spect_data$scen)

p2 <- big_spect_data %>% 
  filter(n > 0 & Fm == 1 & w > 0.1 & scen %in% c("Physio. + Resource (exp.)",
                                                 "Resource (exp.)",
                                                 "Resource (obs.)",
                                                 "Physio. + Resource (obs.)",
                                                 "Physio.")) %>%
  ggplot(., aes(w, re_spec, color = factor(scen), group = sim)) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_vline(data = plotdf, aes(xintercept = w_mat), color = "red", linetype = "dotted") +
  geom_line(size = 1) + 
  #scale_colour_viridis(option = "cividis", discrete = T) +
  scale_colour_manual(values = rev(pal[]),
                      name = "") +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  theme_classic(base_size = 13) +
  #scale_y_log10(breaks = sim) +
  scale_x_log10() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10)) +
  guides(color = FALSE) +
  labs(x ="Body mass (g)",
       y = "Relative abundance (warming/no warming)",
       color = "Scenario") +
  NULL

p2

# Plot relative and absolute in the same plot!
p2+p1 

#ggsave("baltic/figures/spectra_project.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")



#**** Plot mortality ===============================================================
# Absolute mortality
big_spect_data %>% 
  filter(n > 0 & Fm == 1 & scen %in% c("Physio. + Resource (exp.)",
                                       "Physio. + Resource (obs.)",
                                       "Resource (exp.)",
                                       "Resource (obs.)",
                                       "Physio.")) %>%
  ggplot(., aes(w, mort, color = factor(scen), linetype = scen, group = sim)) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_vline(data = plotdf, aes(xintercept = w_mat), color = "red", linetype = "dotted") +
  geom_line(size = 1) + 
  scale_colour_manual(values = rev(pal)) +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  theme_classic(base_size = 13) +
  scale_x_log10() +
  guides(linetype = FALSE) +
  theme(legend.position = "bottom",
        aspect.ratio = 1/2) +
  labs(x ="Body mass (g)",
       y = "Predation mortality",
       color = "Scenario") +
  NULL

# Relative mortality (filter really low values!)
big_spect_data %>% 
  filter(Fm == 1 & mort > 0.075 & scen %in% c("Physio. + Resource (exp.)",
                                              "Physio. + Resource (obs.)",
                                              "Resource (exp.)",
                                              "Resource (obs.)",
                                              "Physio.")) %>%
  ggplot(., aes(w, re_mort, color = factor(scen), linetype = scen, group = sim)) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_vline(data = plotdf, aes(xintercept = w_mat), color = "red", linetype = "dotted") +
  geom_line(size = 1) + 
  scale_colour_manual(values = rev(pal)) +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  theme_classic(base_size = 13) +
  scale_x_log10() +
  guides(linetype = FALSE) +
  theme(#legend.position = "bottom",
        aspect.ratio = 1/2) +
  labs(x ="Body mass (g)",
       y = "Relative predation mortality (warming/no warming)",
       color = "Scenario") +
  NULL
#ggsave("baltic/figures/supp/mort_project.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


#**** Plot feeding level ===========================================================
# Absolute feeding level
big_spect_data %>% 
  filter(Fm == 1 & scen %in% c("Physio. + Resource (exp.)",
                               "Physio. + Resource (obs.)",
                               "Resource (exp.)",
                               "Resource (obs.)",
                               "Physio.")) %>%
  ggplot(., aes(w, feedingLevel, color = factor(scen), group = sim)) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_line(size = 1) + 
  scale_colour_manual(values = rev(pal), 
                      labels = c("Physiology", "Physiology + Resource")) +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  theme_classic(base_size = 13) +
  scale_x_log10() +
  ylim(0, 1) +
  theme(legend.position = "bottom",
        aspect.ratio = 1/2) +
  labs(x ="Body mass (g)",
       y = "Feeding level",
       color = "Scenario") +
  NULL

# plotFeedingLevel(ref) + facet_wrap(~Species, nrow = 3)
# refFL %>% 
#   ggplot(., aes(w, value)) + 
#   geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
#   geom_line(size = 1) + 
#   scale_colour_manual(values = rev(pal)) +
#   facet_wrap(~ Species, scales = "free", nrow = 3) +
#   theme_classic(base_size = 13) +
#   scale_x_log10() +
#   ylim(0, 1) +
#   theme(legend.position = "bottom",
#         aspect.ratio = 1/2) +
#   labs(x ="Body mass (g)",
#        y = "Feeding level",
#        color = "Scenario") +
#   NULL

# Relative feeding level
big_spect_data %>% 
  filter(Fm == 1 & scen %in% c("Physio. + Resource (exp.)",
                               "Physio. + Resource (obs.)",
                               "Resource (exp.)",
                               "Resource (obs.)",
                               "Physio.")) %>%
  ggplot(., aes(w, re_feedingLevel, color = factor(scen), linetype = scen, group = sim)) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_line(size = 1) + 
  scale_colour_manual(values = rev(pal)) +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  theme_classic(base_size = 13) +
  scale_x_log10() +
  guides(linetype = FALSE) +
  theme(aspect.ratio = 1/2) +
  labs(x ="Body mass (g)",
       y = "Relative Feeding level (warming/no warming)",
       color = "Scenario") +
  NULL

#ggsave("baltic/figures/supp/feedingLevel_project.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


#** More tests =====================================================================
# What is the SSB here in these runs?
tt <- project(pars_with_res, 
              dt = dt,
              effort = projectEffort_m,
              temperature = projectTemp$temperature,
              diet_steps = 10,
              t_max = t_max)   

zz <- project(pars_no_res, 
              dt = dt,
              effort = projectEffort_m,
              temperature = projectTemp$temperature,
              diet_steps = 10,
              t_max = t_max)   

str(ref@effort)

# no resource - empirical temp - FMSy fishing
tail(getSSB(zz), 2)
# with resource - empirical temp - FMSy fishing
tail(getSSB(tt), 2)
# ref - constant temp - FMSy fishing
tail(getSSB(ref), 2)


# #**** Plot predation mortality =====================================================
# Compare the predation mortalities from the two extremes:

projectEffort_new <- projectEffort_m
projectEffort_new[100:nrow(projectEffort_m), ] <- projectEffort_m[100:nrow(projectEffort_m), ] * 1.2

proj <- project(pars_with_res,
                dt = dt,
                effort = projectEffort_new,
                temperature = projectTemp$temperature,
                diet_steps = 10,
                t_ref = t_ref)

warmMort <- getMortality(proj)
warmMort$scen <- "warm"

refs <- project(pars_with_res,
                dt = dt,
                effort = projectEffort_m,
                temperature = consTemp,
                diet_steps = 10,
                t_ref = t_ref)

consMort <- getMortality(refs)
consMort$scen <- "constant"

mortdat <- rbind(consMort, warmMort)

# More or less identical predation mortality
ggplot(mortdat, aes(w, value, color = scen, linetype = scen)) +
  geom_line() +
  scale_x_log10() +
  NULL



# plotM2(ref)
# plotM2(proj)

## Testing if I get the same change in spectra when scaling h directly, not via temperature
ref2 <- project(pars_no_res, 
                dt = dt,
                effort = projectEffort_m,
                temperature = consTemp,
                diet_steps = 10,
                t_max = t_max)

ref2Spect <- getSpectra(ref2)

# Create new mizer params object
test_pars <- pars_no_res@species_params
test_pars$h <- pars_no_res@species_params$h * 1.1

test_pars2 <- MizerParams(test_pars, 
                          ea_gro = 0,
                          ea_car = 0, 
                          kappa_ben = kappa_ben,
                          kappa = kappa,
                          w_bb_cutoff = w_bb_cutoff,
                          w_pp_cutoff = w_pp_cutoff,
                          r_pp = r_pp,
                          t_ref = t_ref)

test <- project(pars_with_res, 
                dt = dt,
                effort = projectEffort_m,
                temperature = consTemp,
                diet_steps = 10,
                t_max = t_max)   

testSpect <- getSpectra(test)

# Calculate relative spectra
testSpect$re_spec <- testSpect$n / ref2Spect$n

head(testSpect)
head(ref2Spect)

# Plot spectra
testSpect %>% 
  filter(n > 0 ) %>%
  ggplot(., aes(w, re_spec)) + 
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", size = 0.7, alpha = 0.6) +
  geom_vline(data = plotdf, aes(xintercept = w_mat), color = "red", linetype = "dotted") +
  geom_line(size = 1) + 
  scale_colour_manual(values = rev(pal)) +
  facet_wrap(~ species, scales = "free", nrow = 3) +
  theme_classic(base_size = 13) +
  scale_x_log10() +
  theme(legend.position = "bottom",
        aspect.ratio = 1/2,) +
  labs(x ="Body mass (g)",
       y = "Relative abundance (warming/no warming)",
       color = "Scenario") +
  NULL

# Yes! It happens also when i ONLY change h... but how this affects SSB and such depends on the
# large fish abundance in the population

tail(getSSB(test))
tail(getSSB(ref2))


