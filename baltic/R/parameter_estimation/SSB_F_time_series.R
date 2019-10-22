#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.03.06: Max Lindmark
#
# - The following code plots time series of F and (relative) abundance for 
# mizer-implementation of the Baltic Sea (Lindmark et al., in prep)
# 
# A. Load libraries
#
# B. Plot F and (relative) abundance
#      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Questions in progress marked with ***

#======== A. LOAD LIBRARIES =========================================================
rm(list = ls())

# Provide package names
pkgs <- c("ggplot2", "dplyr", "tidyr", "viridis", "devtools")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-baltic-sea/master/R/Functions/package_info.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
pkg_info(pkgs)

# Package versions
# 1 devtools  1.13.6
# 2    dplyr   0.7.5
# 3  ggplot2   3.0.0
# 4    tidyr   0.8.1
# 5  viridis   0.5.1

# devtools::install_github("thomasp85/patchwork")
library(patchwork)
# packageVersion("patchwork") # v.0.0.1


#======== B. PLOT F AND (RELATIVE) ABUNDANCE ========================================
# Read stock assessment and survey index data. See Appendix S1.
dat <- read.csv("Data/SSB_F/SSB_F_data.csv", sep = ";")

head(dat)

dat <- dat %>% rename(Fm = F) # Bad idea to call a column F in R

# Create dataframes with mean SSB and F for calibration period plotting
min_cal_yr <- 1992
max_cal_yr <- 2002

ref_time <- data.frame(Year = c(min_cal_yr, max_cal_yr),
                       TSB  = c(0, max(dat$TSB, na.rm = T)),
                       SSB  = c(0, max(dat$SSB, na.rm = T)),
                       Fm    = c(0, max(dat$Fm, na.rm = T)))

mean_ssb_F <- dat %>% 
  filter(Year >= min_cal_yr & Year <= max_cal_yr) %>% 
  group_by(Species) %>% 
  mutate(mean_SSB = mean(SSB),
         mean_F   = mean(Fm))

# Plot trawl survey index
tsb <- ggplot(dat, aes(Year, TSB, color = Species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = min(TSB),
                ymax = max(TSB)),
            fill  = "black", alpha = 0.1) +
  geom_line(size = 2) +
  #guides(color = FALSE) +
  theme(aspect.ratio = 1) +
  scale_color_viridis(discrete = TRUE) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL
  
# Plot SSB
ssb <- ggplot(dat, aes(Year, SSB, color = Species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = max(SSB)),
            fill  = "black", alpha = 0.1) +
  geom_line(size = 2) +
  theme(aspect.ratio = 1) +
  ylab("SSB (10^3 tonnes)") +
  scale_color_viridis(discrete = TRUE) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL

tsb + ssb


# Plot SSB and F, overlay mean for each period
ssb2 <- dat %>% 
  filter(Year > 1979 & Year < 2010 & Species %in% c("Cod", "Sprat", "Herring")) %>%
  ggplot(., aes(Year, SSB, color = Species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = max(SSB)),
            fill  = "black", alpha = 0.1) +
  geom_line(size = 1.2) +
  geom_point(data = filter(mean_ssb_F, Species %in% c("Cod", "Sprat", "Herring")), 
             aes(Year, mean_SSB), size = 1.2) +
  theme(aspect.ratio = 1) +
  ylab("SSB (1000 tonnes)") +
  scale_color_viridis(discrete = TRUE) +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(color = FALSE) +
  theme(aspect.ratio = 1) +
  NULL

f <- dat %>% 
  filter(Year > 1979 & Year < 2010 & Species %in% c("Cod", "Sprat", "Herring")) %>%
  ggplot(., aes(Year, Fm, color = Species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = 0,
                ymax = max(Fm)),
            fill  = "black", alpha = 0.1) +
  geom_line(size = 1.2) +
  geom_point(data = filter(mean_ssb_F, Species %in% c("Cod", "Sprat", "Herring")), 
             aes(Year, mean_F), size = 1.2) +
  theme(aspect.ratio = 1) +
  ylab("F") +
  scale_color_viridis(discrete = TRUE) +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 1) +
  NULL

p <- ssb2 + f

ggsave("Figures/SSB_F.tiff", plot = p, dpi = 300, width = 19, height = 19, units = "cm")


# Get mean SSB and F for Appendix and calibration data
mean_ssb_F %>% 
  group_by(Species) %>% 
  summarize(mean(mean_SSB),
            mean(mean_F))

# Here we have some options: 
# 1) Use SSB and ignore Flounder when calibrating
# 2) Use SSB and assume same ratio of Flounder:Cod in SSB as in survey
# I will try option 1