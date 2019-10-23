#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2018.11.20: Max Lindmark
#
# - The following code fits VBGE curves to length-at-age data for species used in 
# mizer-implementation of the Baltic Sea (Lindmark et al., in prep)
# 
# A. Load up libraries
#
# B. Cleanup BITS data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# When doing a fresh start I need to check I'm in the right libpath to get the right mizer version
# .libPaths()
# .libPaths("C:/Program Files/R/R-3.5.0/library")

# Load libraries, install if needed
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(FSA)
library(FSAdata)
library(nlstools)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Print package versions
# print(sessionInfo())
# other attached packages:
# [1] mizer_1.1 nlstools_1.0-2 FSAdata_0.3.6 FSA_0.8.22 testthat_2.0.0    
# [6] patchwork_0.0.1 dplyr_0.8.3 viridis_0.5.1 viridisLite_0.3.0 tidyr_0.8.1       
# [11] magrittr_1.5 RCurl_1.95-4.10 bitops_1.0-6 RColorBrewer_1.1-2 usethis_1.4.0     
# [16] devtools_2.0.1 ggplot2_3.1.1      


# B. CLEANUP BITS DATA =============================================================
dat <- read.csv("baltic/data/BITS/Exchange Data_2019-04-04 14_32_48_no_flounder.csv", sep = ";")

str(dat)
dat$AgeRings <- as.numeric(dat$AgeRings)

dat <- dat %>%
  filter(SpecCode %in% c("126436", "164712",
                         "126417", "161722",
                         "126425", "161789")
         & AgeRings > -1)

dat$Species <- "Cod"

dat$Species <- ifelse(dat$SpecCode %in% c("126417", "161722"),
                      "Herring",
                      dat$Species)

dat$Species <- ifelse(dat$SpecCode %in% c("126425", "161789"),
                      "Sprat",
                      dat$Species)

unique(dat$Species)

# Get length in cm for all individuals
# https://vocab.ices.dk/?ref=18
unique(dat$LngtCode)

ggplot(dat, aes(LngtCode)) +
  facet_wrap(~Species, ncol = 2, scales = "free_y") +
  geom_bar() +
  theme_bw(base_size = 14) +
  NULL

dat$Length_cm <- ifelse(dat$LngtCode %in% c(".", "0"), 
                        dat$LngtClass / 10,
                        dat$LngtClass)

ggplot(dat, aes(AgeRings, LngtClass)) +
  facet_wrap(~Species, ncol = 2, scales = "free_y") +
  geom_point() +
  theme_bw(base_size = 14) +
  NULL

# Check again with cm instead of LngtClass
ggplot(dat, aes(AgeRings, Length_cm)) +
  facet_wrap(~Species, ncol = 2, scales = "free_y") +
  geom_point() +
  theme_bw(base_size = 14)

# Some outliers, nothing evidently systematic though, so will just remove them
dat$keep <- "Y"

dat$keep <- ifelse(dat$Species == "Sprat" & dat$Length_cm > 20,
                   "N", 
                   dat$keep)

dat$keep <- ifelse(dat$Species == "Herring" & dat$Length_cm < 5,
                   "N", 
                   dat$keep)

dat <- dat %>% filter(Length_cm < 150 & keep == "Y")

# Check again
ggplot(dat, aes(AgeRings, Length_cm)) +
  facet_wrap(~Species, ncol = 2, scales = "free_y") +
  geom_point() +
  theme_bw(base_size = 14)

# Subset areas
#https://www.google.com/search?rlz=1C1CHBD_svSE722SE722&biw=1473&bih=838&tbm=isch&sa=1&ei=37ytXK_eFu2orgTws6XQBg&q=iceas+subdivision+rectangle+baltic&oq=iceas+subdivision+rectangle+baltic&gs_l=img.3...2184.2363..2597...0.0..0.134.266.0j2......1....1..gws-wiz-img.PqfHGtjaEaw#imgrc=aHrvX3FTrNJbpM:

dat$AreaCode2 <- dat$AreaCode 

dat <- dat %>%
  separate(AreaCode2, c("Area_s", "Code_s"), sep = 2)

head(dat)

ggplot(dat, aes(Code_s, fill = Area_s)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dat <- dat %>%
  filter(., Code_s %in% c("G5", "G6", "G7", "G8", "G9", "H0", "H1", "H2", "H3", "G4", "H5", "H6", "H7", "H8", "H9"))

ggplot(dat, aes(Code_s, fill = Area_s)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Filter out G4 38
dat$keep2 <- ifelse(dat$Area_s > 36 & dat$Area_s < 39 & dat$Code_s == "G4",
                    "N",
                    "Y")

dat <- dat %>% filter(keep2 == "Y")

ggplot(dat, aes(Code_s, fill = Area_s)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~Species, scale = "free_y") +
  theme_bw(base_size = 14) +
  NULL

str(dat)

#write.csv(dat, "data/BITS/clean_BITS.csv")


