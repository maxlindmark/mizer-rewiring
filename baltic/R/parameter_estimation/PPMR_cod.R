#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.03.19: Max Lindmark
#
# - The following code estimates predator-prey body mass ratio (PPMR) for Baltic Cod
# following the approach in Reum et al (2018) Ecology & Evolution based on stomach 
# data available at http://ecosystemdata.ices.dk/stomachdata/download.aspx 
# 
# A. Load up libraries and data
#
# B. Clean & Explore data
#
# C. Calculate PPMR
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES AND READ DATA ========
rm(list = ls())

# Load packages

# packageVersion("patchwork")
library(ggplot2)   # v3.0.0
library(dplyr)     # v0.7.5
library(tidyr)     # v0.8.1
library(patchwork) # v0.0.1

# Stomach data: http://ecosystemdata.ices.dk/stomachdata/download.aspx 
#  (filter only by cod, Gadus morhua)
dat <- read.csv("baltic/data/PPMR/cod/StomachDataSet201932090580.csv", sep = ",")

head(dat)
str(dat)

# Filter data that covers the calibration period
min_cal_yr <- 1992
max_cal_yr <- 2002

dat$calib_yr <- ifelse(dat$year > 1991 & dat$year < 2003, "Y", "N")

# Check sample size/year
ggplot(dat, aes(factor(year), fill = calib_yr)) + 
  geom_bar() +
  theme_bw(base_size = 18) +
  NULL

# Check digestion stage in calibration time period
# http://www.ices.dk/marine-data/data-portals/Documents/StomachData-Baltic.pdf
dat %>% 
  filter(calib_yr == "Y") %>% 
  ggplot(., aes(factor(year), fill = Prey_DigestionStage)) + 
  geom_bar() +
  theme_bw(base_size = 18) +
  NULL

# Check species in calibration time period and digestion == 0
dat %>% 
  filter(calib_yr == "Y" & Prey_DigestionStage == "0") %>% 
  ggplot(., aes(factor(year), fill = Prey_LatinName)) + 
  geom_bar() +
  theme_bw(base_size = 18) +
  NULL

# Check dataset-source in calibration time period and digestion == 0
dat %>% 
  filter(calib_yr == "Y" & Prey_DigestionStage == "0") %>% 
  ggplot(., aes(factor(year), fill = Dataset)) + 
  geom_bar() +
  theme_bw(base_size = 18) +
  NULL

# Check country in calibration time period and digestion == 0
dat %>% 
  filter(calib_yr == "Y" & Prey_DigestionStage == "0") %>% 
  ggplot(., aes(factor(year), fill = Country)) + 
  geom_bar() +
  theme_bw(base_size = 18) +
  NULL

# Ok, this means that for year prior to 2002, prey weight is not individual 
# weight but per prey species. So when there is n prey, I need to divide
# the weight by n prey

# Filter spatial distribution
# Split AreaCode-factor to more easily split by ICES sub division
dat$AreaCode2 <- dat$ICES_StatRec 

dat <- dat %>% separate(AreaCode2, c("Area_s", "Code_s"), sep = 2)
dat$Area_s <- as.numeric(dat$Area_s)

# Plot to see how many rows are outside the area
ggplot(dat, aes(Code_s, fill = factor(Area_s))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Filter to only get Sub Dvision 25-29 + 32
dat <- dat %>%
  filter(., Code_s %in% c("G4", "G5", "G6", "G7", "G8", "G9", "H0", "H1", 
                          "H2", "H3", "G4", "H5", "H6", "H7", "H8", "H9"))

# Check area codes again after filtering
ggplot(dat, aes(Code_s, fill = factor(Area_s))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Need to remove some rectangles within G4
dat$keep <- ifelse(dat$ICES_StatRec %in% c("37G4", "38G4"),
                   "N", "Y")

dat <- subset(dat, keep == "Y")

# Check area codes again after filtering
ggplot(dat, aes(Code_s, fill = factor(Area_s))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dat <- dat %>% filter(calib_yr == "Y" & 
                      Prey_DigestionStage == "0" & 
                      Stomach_Empty == 0)

# First need to filter stomachs with 0 "Prey_TotalNo" but with a weight, because 
# that must be an error. Because if it's also before 2000, I need to now the # of
# prey to get the mean weigth of individual prey

# hmm surprising number of 0's in total number of prey, even though I've filtered
# by non-empty stomachs AND there is a weight... However, since I'm anyway going 
# to use the average weights of the prey within each stomach, and that information 
# is the weight information after 2000, I only need to filter away 0 "Prey_TotalNo"
# before 2000...
ggplot(dat, aes(year, fill = factor(Prey_TotalNo))) +
  geom_bar()

dat$prey_NO_info <- ifelse(dat$Prey_TotalNo == 0 & dat$year < 2000,
                           "N", "Y")

dat <- dat %>% filter(prey_NO_info == "Y")

ggplot(dat, aes(year, fill = factor(Prey_TotalNo))) +
  geom_bar()

# Apparently SampleNo(FishID) is not the indivual ID, because that gives
# multiple rows of a single species in data before 2000 (when it was averaged by prey)
# Instead it's ICES_ItemID

# Calculate n prey by predator, plot summary
dat_sum <- dat %>% 
  group_by(ICES_ItemID) %>% 
  dplyr::mutate(n = n())

ggplot(dat_sum, aes(factor(n))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# In years after 2002, how can there be more Prey_TotalNo than number of rows within
# an individual? Each row should be an individual prey, so what does that column even mean?
ggplot(dat_sum, aes(factor(Prey_TotalNo), fill = factor(year))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Here I'm plotting by RecordType, (SS is code for single stomach)
# However I'm missing that information in this data-set...
ggplot(dat_sum, aes(factor(RecordType), fill = factor(year))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Divide prey weight by Prey_TotalNo if older than 2000
dat_sum$ind_weight <- ifelse(dat_sum$year < 2000, 
                             dat_sum$Prey_Weight / dat_sum$Prey_TotalNo,
                             dat_sum$Prey_Weight)

# Filter 0 weights
dat_sum <- dat_sum %>% filter(Prey_Weight > 0)

glimpse(dat_sum)

# It seems most records have 2 species but no one has more...
ggplot(dat_sum, aes(factor(n))) +
  geom_histogram(stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

head(data.frame(dat_sum), 25)

# Calculate average weight of prey within each individual (to avoid pseudoreplication)
dat_agg <- data.frame(dat_sum %>% 
  dplyr::group_by(ICES_ItemID) %>% 
  dplyr::mutate(mean_prey_w = mean(ind_weight)))

# Ok, looks like there's no weight in the data, need to convert from length!
dat_agg$pred_weight <- 0.01*dat_agg$Predator_Lengh.mean.^3

# Calculate PPMR
dat_agg$PPMR <- dat_agg$pred_weight / dat_agg$mean_prey_w

# Calculate log mean and standard deviation
dat_agg %>%
  distinct(ICES_ItemID, .keep_all = TRUE) %>% 
  summarize(mean   = mean(log10(PPMR)),
            sd     = sd(log10(PPMR)),
            mean_n = 10^(mean(log10(PPMR))),
            sd_n   = 10^(sd(log10(PPMR))))

# Plot log10(PPMR)
p1 <- dat_agg %>% 
  distinct(ICES_ItemID, .keep_all = TRUE) %>% 
  ggplot(., aes(log10(PPMR))) +
  geom_histogram(size = 3) +
  geom_vline(xintercept = 2.63, color = "red", linetype = 2, size = 1) +
  coord_cartesian(expand = 0) +
  annotate("text", label = "mean PPMR=426\nsd PPMR=5.63", 
           x = Inf, y = Inf, size = 4.5, col = "blue",
           vjust = 1, hjust = 1) +
  labs(y = "Count") +
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))
  
ggsave("baltic/figures/supp/PPMR.png", width = 3.5, height = 3.5, dpi = 600)

