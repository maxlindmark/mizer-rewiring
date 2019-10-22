#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2018.11.20: Max Lindmark
#
# - The following code fits VBGE curves to length-at-age data for species used in 
# mizer-implementation of the Baltic Sea (Lindmark et al., in prep)
# 
# A. Load up libraries
#
# B. Cleanup BITS data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES =========================================================
rm(list = ls())

# Provide package names
pkgs <- c("ggplot2", "dplyr", "tidyr", "FSA", "FSAdata", "nlstools")

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
# 1    dplyr   0.7.5
# 2      FSA  0.8.22
# 3  FSAdata   0.3.7
# 4  ggplot2   3.0.0
# 5 nlstools   1.0-2
# 6    tidyr   0.8.1

# See Appendix S1 and Readme.md for details about data sources.
# For individual-level data (age, length, weight) I use data from the BITS
#  survey (downloaded manually from DATRAS):
#  http://www.ices.dk/marine-data/data-portals/Pages/DATRAS.aspx

#    Selection: 
#      Exchange
#      CA
#      BITS
#      All quarters
#      1992-2002
#      All ships
#      Downloaded: 2019.04.04
#      Species codes:
#        https://datras.ices.dk/Data_products/qryspec.aspx
#        cod:      "126436", "164712"
#        flounder: "127141", "172894"
#        herring:  "126417", "161722"
#        sprat:    "126425", "161789"


#======== B. CLEANUP BITS DATA ======================================================
dat <- read.csv("Data/BITS/Exchange Data_2019-04-04 14_32_48.csv", sep = ",")

str(dat)
dat$AgeRings <- as.numeric(dat$AgeRings)

dat <- dat %>%
  filter(SpecCode %in% c("126436", "164712",
                         "127141", "172894",
                         "126417", "161722",
                         "126425", "161789")
         & AgeRings > -1)

dat$Species <- "Cod"

dat$Species <- ifelse(dat$SpecCode %in% c("127141", "172894"),
                      "Flounder",
                      dat$Species)

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

ggplot(dat, aes(AgeRings, Length_cm)) +
  facet_wrap(~Species, ncol = 2, scales = "free_y") +
  geom_point() +
  theme_bw(base_size = 14)

# Some errors I think, let's look at flounder
dat %>% 
  filter(Species == "Flounder") %>% 
  ggplot(., aes(AgeRings, Length_cm, color = LngtCode)) +
  facet_wrap(~Species, ncol = 2, scales = "free_y") +
  geom_point() +
  theme_bw(base_size = 14)

# Obviously a typo, will scale them down
dat$Length_cm <- ifelse(dat$Species == "Flounder" & dat$Length_cm > 100,
                        dat$Length_cm / 10,
                        dat$Length_cm)

# Check again
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


