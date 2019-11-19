#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.11.16: Max Lindmark
#
# TESTING EFFECTS OF WARMING ON FOOD WEB METRICS
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
t_max <- 200
dt <- 0.2
t_ref <- 10
kappa_ben = 1
kappa = 1
w_bb_cutoff = 20
w_pp_cutoff = 1
r_pp = 4
r_bb = 4


# B. CONTRAST WARM AND COLD ========================================================
# Then compare that to a projection with a constant temperature

# Create effort vector
effort <- c(Cod = params@species_params$AveEffort[1], 
           Herring = params@species_params$AveEffort[3], 
           Sprat = params@species_params$AveEffort[2])

#**** Project without temp (reference) =============================================
ref <- project(params, 
               dt = dt,
               effort = effort,
               temperature = rep(t_ref, t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = t_ref)   


#**** Project with temp and bad activation energies ================================
bad_par_sp <- params@species_params

bad_par_sp$ea_met <- mean(ea$met)
bad_par_sp$ea_int <- mean(ea$int)
bad_par_sp$ea_mor <- mean(ea$mor)

bad_par <- MizerParams(bad_par_sp,
                       ea_gro = mean(ea$gro),
                       ea_car = mean(ea$car),
                       kappa_ben = kappa_ben,
                       kappa = kappa,
                       w_bb_cutoff = w_bb_cutoff,
                       w_pp_cutoff = w_pp_cutoff,
                       r_pp = r_pp,
                       r_bb = r_bb)

bad <- project(bad_par, 
               dt = dt,
               effort = effort,
               temperature = rep(12, t_max),
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10)   


#**** Compare metrics ==============================================================
#****** Growth =====================================================================
a <- plotGrowthCurves(ref) + facet_wrap(~ Species, scales = "free") + ggtitle("ref")
b <- plotGrowthCurves(bad) + facet_wrap(~ Species, scales = "free") + ggtitle("warm")

a/b # Growth is higher in warm than cold

#****** Mortality ==================================================================
a <- plotM2(ref) + ggtitle("ref")
b <- plotM2(bad) + ggtitle("warm")

a/b # Mortality is higher in warm than cold

#****** Size-spectra ==================================================================
# Load function for extracting numbers-at-age
func <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/getSpectra.R", ssl.verifypeer = FALSE)
eval(parse(text = func))

plotSpectra(ref)

a <- plotSpectra(ref, benthos = FALSE, algae = FALSE, plankton = FALSE) + facet_wrap(~ Species, scales = "free") + ggtitle("ref")
b <- plotSpectra(bad, benthos = FALSE, algae = FALSE, plankton = FALSE) + facet_wrap(~ Species, scales = "free") + ggtitle("warm")

aa <- getSpectra(ref)
bb <- getSpectra(bad)

aa$scen <- "cold"
bb$scen <- "warm"


df <- data.frame(rel = bb$n / aa$n,
                 w   = aa$w,
                 species = aa$species)

ggplot(df, aes(log10(w), rel)) + 
  facet_wrap(~ species, scales = "free") + 
  geom_line(size = 1) +
  geom_hline(yintercept = 1, color = "red") +
  NULL

## Plot together
cc <- rbind(aa, bb)

cc$t_hold <- "N"
cc$t_hold <- ifelse(cc$species == "Cod" & log10(cc$n) > -5,
                    "Y",
                    cc$t_hold)

cc$t_hold <- ifelse(cc$species == "Herring" & log10(cc$n) > -1.5,
                    "Y",
                    cc$t_hold)

cc$t_hold <- ifelse(cc$species == "Sprat" & log10(cc$n) > -0.6,
                    "Y",
                    cc$t_hold)

cc %>% filter(t_hold == "Y") %>% 
  ggplot(., aes(log10(w), log10(n), color = scen)) + 
  facet_wrap(~species, scales = "free", ncol = 1) + 
  geom_line(size = 0.5) +
  theme_classic() +
  NULL
  
a/b # Mortality is higher in warm than cold




