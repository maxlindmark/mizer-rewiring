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

#======== A. LOAD LIBRARIES =========================================================
# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(ggplot2)
library(devtools)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(viridis)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# UPDATE
# Print package versions
# print(sessionInfo())
# other attached packages:
# [1] mizer_1.1 nlstools_1.0-2 FSAdata_0.3.6 FSA_0.8.22 testthat_2.0.0    
# [6] patchwork_0.0.1 dplyr_0.8.3 viridis_0.5.1 viridisLite_0.3.0 tidyr_0.8.1       
# [11] magrittr_1.5 RCurl_1.95-4.10 bitops_1.0-6 RColorBrewer_1.1-2 usethis_1.4.0     
# [16] devtools_2.0.1 ggplot2_3.1.1      


#======== B. PLOT F AND (RELATIVE) ABUNDANCE ========================================
# Read stock assessment and survey index data. See Appendix S1.
dat <- read.csv("baltic/data/SSB_F/SSB_F_data.csv", sep = ";")

head(dat)

# Create dataframes with mean SSB and F for calibration period plotting
min_cal_yr <- 1992
max_cal_yr <- 2002

ref_time <- data.frame(Year = c(min_cal_yr, max_cal_yr),
                       TSB  = c(0, max(dat$TSB, na.rm = T)),
                       SSB  = c(0, max(dat$SSB, na.rm = T)),
                       Fm    = c(0, max(dat$Fm, na.rm = T)))

mean_ssb_F <- dat %>% 
  filter(Year >= min_cal_yr & Year <= max_cal_yr) %>% 
  dplyr::group_by(Species) %>% 
  dplyr::mutate(mean_SSB = mean(SSB),
         mean_F   = mean(Fm))

# Define palette
col <- rev(colorRampPalette(brewer.pal(5, "Dark2"))(5))

# Plot trawl survey index
tsb <- ggplot(dat, aes(Year, TSB, color = Species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE, 
            aes(xmin = min(Year), 
                xmax = max(Year),
                ymin = min(TSB),
                ymax = max(TSB)),
            fill = "black", alpha = 0.1) +
  geom_line(size = 2, alpha = 0.8) +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = col) +
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
  geom_line(size = 2, alpha = 0.8) +
  theme(aspect.ratio = 1) +
  ylab("SSB (10^3 tonnes)") +
  scale_color_manual(values = col) +
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
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(data = filter(mean_ssb_F, Species %in% c("Cod", "Sprat", "Herring")), 
             aes(Year, mean_SSB), size = 1.2) +
  theme(aspect.ratio = 1) +
  ylab("Spawning stock biomass\n(1000 tonnes)") +
  scale_color_manual(values = col) +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(color = FALSE) +
  theme(aspect.ratio = 1) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
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
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(data = filter(mean_ssb_F, Species %in% c("Cod", "Sprat", "Herring")), 
             aes(Year, mean_F), size = 1.2) +
  theme(aspect.ratio = 1) +
  ylab("F") +
  scale_color_manual(values = col) +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 1) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

ssb2 + f

#ggsave("baltic/figures/supp/SSB_F.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")


# Get mean SSB and F for Appendix and calibration data
mean_ssb_F %>% 
  group_by(Species) %>% 
  summarize(mean(mean_SSB),
            mean(mean_F))

# Here we have some options: 
# 1) Use SSB and ignore Flounder when calibrating
# 2) Use SSB and assume same ratio of Flounder:Cod in SSB as in survey
# I will try option 1