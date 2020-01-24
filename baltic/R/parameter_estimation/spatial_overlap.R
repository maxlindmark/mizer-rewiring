#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2018.11.20: Max Lindmark
#
# - The following code estimates spatial overlap according to Schoener's method 
# (described in Blanchard et a (2014) Journal of Applied Ecology, based on trawl
# survey data (BITS, BIAS). This is used to scale the encounter rates in predatory
# interactions within and between species. 
# 
# A. Load up libraries and data
#
# B. Clean data
#
# C. Calculate spatial overlap index (pairwise) for interactions matrix
#
# D. Plot distribution maps for each species
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES AND READ DATA ===========================================
rm(list = ls())

# Provide package names
pkgs <- c("ggplot2", "tidyr", "dplyr", "maps", "viridis", "RCurl")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
# script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-baltic-sea/master/R/Functions/package_info.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# pkg_info(pkgs)

# Package versions
# 1   dplyr   0.7.5
# 2 ggplot2   3.0.0
# 3    maps   3.3.0
# 4   tidyr   0.8.1
# 5 viridis   0.5.1


#======== B. READ DATA ==============================================================
#====** BITS (Trawl Survey) ========
# See cleanup_BITS_data for post download processing of data
dat_bent <- read.csv("baltic/data/BITS/clean_BITS.csv", sep = ",")

head(dat_bent)

#====** Acoustic Survey ========
# BIAS (sprat, herring)
dat_pela <- read.csv("baltic/data/BIAS/sprat_herring_acoustic_abundance_1992_2002.csv", 
                     sep = ";", dec = ",")

head(dat_pela)
str(dat_pela)

# Standardize column names & filter year
dat_pela <- dplyr::rename(dat_pela, 
                          Year     = ANNUS,
                          Species  = SPECIES,
                          AreaCode = RECT)

dat_pela <- dat_pela %>% filter(Year >= 1991 & #max(dat_pela$Year)-2 &
                                Year <= 2002 & #max(dat_pela$Year) &
                                Sub_Div %in% c("25", "26", "27", "28_2", "29", "32"))

head(dat_pela)

# Before calculating the interaction matrix, I will prep the data so that it contains only abundance (pelagics) or CPUE (demersal) of mature and immature fish, AreaCode (ICES rectangle) and abundance of mature and immature fish. This process is different for each data sources (in pelagics abundance is a column whereas in demersal each observation is a unique row), as well as species (different size at maturity). I will therefore clean data by species. Then I apply a function to calculate the species overlap index (see section C.). 


#====** Cod ========
cod <- filter(dat_bent, Species == "Cod")

# Categorize immature/mature
cod$mature <- ifelse(cod$Length_cm > 30, "Adu", "Juv")

# nrow(filter(cod, AreaCode == "37G9" & mature == "Juv")), for testing below summary

# Summarize number of fish in each areacode
cod <- cod %>% 
  dplyr::group_by(AreaCode, mature) %>% 
  dplyr::summarize(n = n())%>% 
  dplyr::ungroup() 

# Go from long to wide to get n on the same row for each AreaCode
cod <- tidyr::spread(cod, mature, n)

# NA's here are really 0, so need to replace
cod[is.na(cod)] <- 0

cod$Species <- "Cod"


#====** Flounder ========
flo <- filter(dat_bent, Species == "Flounder")

# Categorize immature/mature
flo$mature <- ifelse(flo$Length_cm > 24.4, "Adu", "Juv")

# Summarize number of fish in each areacode
flo <- flo %>% 
  dplyr::group_by(AreaCode, mature) %>% 
  dplyr::summarize(n = n())%>% 
  dplyr::ungroup()

# Go from long to wide to get n on the same row for each AreaCode
flo <- spread(flo, mature, n)

# NA's here are really 0, so need to replace
flo[is.na(flo)] <- 0

flo$Species <- "Flounder"

head(flo)


#====** Herring ========
her <- subset(dat_pela, Species == "Herring")

head(her)

# We have abundance by age. Summarize over mature and immature.
# This is based on VBGE plot of length-at-age, maturation size (Appendix S1) and http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2018/WGBFAS/01%20WGBFAS%20Report%202018.pdf, table 4.2.8 (90% of 3 yrs are mature)
her$Juv <- rowSums(her[, 6:8])
her$Adu <- rowSums(her[, 8:14])

head(her)

her <- subset(her, select = c("Species", "AreaCode", "Juv", "Adu"))

# Summarize number of fish in each AreaCode
her <- her %>% 
  dplyr::group_by(AreaCode) %>% 
  dplyr::summarize(Juv = sum(Juv),
            Adu = sum(Adu))%>% 
  dplyr::ungroup()

# NA's here are really 0, so need to replace
her[is.na(her)] <- 0

her$Species <- "Herring"


#====** Sprat ========
spr <- subset(dat_pela, Species == "Sprat")

# We have abundance by age. Summarize over mature and immature.
# This is based on VBGE plot of length-at-age, maturation size (Appendix S1) and http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2018/WGBFAS/01%20WGBFAS%20Report%202018.pdf, section 7.8 (90% of 2 yrs are mature)
spr$Juv <- rowSums(spr[, 6:7])
spr$Adu <- rowSums(spr[, 8:14])

spr <- subset(spr, select = c("Species", "AreaCode", "Juv", "Adu"))

# Summarize number of fish in each areacode
spr <- spr %>% 
  dplyr::group_by(AreaCode) %>% 
  dplyr::summarize(Juv = sum(Juv),
                   Adu = sum(Adu)) %>% 
  dplyr::ungroup()


spr[is.na(spr)] <- 0

spr$Species <- "Sprat"

head(spr)


#======== C. ESTIMATE PAIRWISE SPATIAL OVERLAP ======================================
# This dataframe will be good to have further down... (species, juv and adu as columns)
df_intra <- bind_rows(cod, flo, spr, her)

ggplot(df_intra, aes(AreaCode, fill = Species)) +
  geom_bar()

# We want species + stage as columns and n as rows, one for each AreaCode
# Go from wide to long format to get "n" by species and lifestage (not Juv & Adu)
# Create new variable (spec + stage) because that is our grouping variable later
# Summarize n over lifestage (by Species and AreaCode)
# Make wide again so Species is a column
df_inter <- df_intra %>%
  gather(lifestage, n, 2:3) %>% 
  dplyr::mutate(Spec_str = paste(Species, lifestage, sep = "")) %>% 
  dplyr::group_by(AreaCode, Spec_str) %>% 
  dplyr::summarize(n = sum(n)) %>% 
  spread(Spec_str, n)

head(df_inter)

df_inter[is.na(df_inter)] <- 0

# Function to calculate proportion in each rectangle for species+stage combination (i and j), then find absolute values and calculate overlap index
overlap_inter <- function(spec_i, spec_j, data) {
  
  # Proportion of i & j in each rectangle
  spec_i_v <- select(data, "AreaCode", spec_i)[, 2] / 
    sum(select(data, c("AreaCode", spec_i))[, 2])
  
  spec_j_v <- select(data, c("AreaCode", spec_j))[, 2] / 
    sum(select(data, c("AreaCode", spec_j))[, 2])
  
  # Absolute difference in proportion for each rectangle
  diff_ab <- abs(spec_j_v - spec_i_v)
  
  data.frame(spec_i  = spec_i,
             spec_j  = spec_j,
             overlap = as.numeric((1 - 0.5*sum(diff_ab))))
  
}

# Testing that it becomes 1 if i and j are the same:
# overlap_inter(spec_i = "CodAdu", 
#               spec_j = "CodAdu",
#               data = df_inter)$overlap


# Apply the function pairwise (note one combination is enough since it's absolute values)

# Apologize in advance for the repetitive code.. can try and fix later!
# Create a matrix to hold the overlap-values
names <- unique(df_intra$Species)
inter <- matrix(nrow = 4, ncol = 4)
dimnames(inter) <- list(names, names)

# Fill in the diagonal (juv/adu overlap within species)
inter["Cod", "Cod"] <- overlap_inter(spec_i = "CodAdu", 
                                     spec_j = "CodJuv",
                                     data = df_inter)$overlap

inter["Sprat", "Sprat"] <- overlap_inter(spec_i = "SpratAdu",
                                         spec_j = "SpratJuv", 
                                         data = df_inter)$overlap

inter["Flounder", "Flounder"] <- overlap_inter(spec_i = "FlounderAdu", 
                                               spec_j = "FlounderJuv", 
                                               data = df_inter)$overlap

inter["Herring", "Herring"] <- overlap_inter(spec_i = "HerringAdu",
                                             spec_j = "HerringJuv", 
                                             data = df_inter)$overlap

# Fill in the non-diagonals (mean of overlap between across combinations of species and life stages (in total 4, hence I dive the sum of all theta values (elements in interactions matrix) by 4) 

# Cod / Flounder
inter["Cod", "Flounder"] <-
  (overlap_inter(spec_i = "CodJuv",
                 spec_j = "FlounderAdu", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodAdu",
                 spec_j = "FlounderJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodJuv",
                 spec_j = "FlounderJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodAdu",
                 spec_j = "FlounderAdu", 
                 data   = df_inter)$overlap) / 4

# The interactions between species are symmetrical
inter["Flounder", "Cod"] <- inter["Cod", "Flounder"]

# Cod / Sprat
inter["Cod", "Sprat"] <-
  (overlap_inter(spec_i = "CodJuv",
                 spec_j = "SpratAdu", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodAdu",
                 spec_j = "SpratJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodJuv",
                 spec_j = "SpratJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodAdu",
                 spec_j = "SpratAdu", 
                 data   = df_inter)$overlap) / 4

inter["Sprat", "Cod"] <- inter["Cod", "Sprat"]

# Cod / Herring
inter["Cod", "Herring"] <-
  (overlap_inter(spec_i = "CodJuv",
                 spec_j = "HerringAdu", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodAdu",
                 spec_j = "HerringJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodJuv",
                 spec_j = "HerringJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "CodAdu",
                 spec_j = "HerringAdu", 
                 data   = df_inter)$overlap) / 4

inter["Herring", "Cod"] <- inter["Cod", "Herring"]
  
# Sprat / Flounder
inter["Sprat", "Flounder"] <-
  (overlap_inter(spec_i = "FlounderJuv",
                 spec_j = "SpratAdu", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "FlounderAdu",
                 spec_j = "SpratJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "FlounderJuv",
                 spec_j = "SpratJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "FlounderAdu",
                 spec_j = "SpratAdu", 
                 data   = df_inter)$overlap) / 4

inter["Flounder", "Sprat"] <- inter["Sprat", "Flounder"]

# Sprat / Herring
inter["Sprat", "Herring"] <-
  (overlap_inter(spec_i = "HerringJuv",
                 spec_j = "SpratAdu", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "HerringAdu",
                 spec_j = "SpratJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "HerringJuv",
                 spec_j = "SpratJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "HerringAdu",
                 spec_j = "SpratAdu", 
                 data   = df_inter)$overlap) / 4

inter["Herring", "Sprat"] <- inter["Sprat", "Herring"]


# Flounder / Herring
inter["Flounder", "Herring"] <-
  (overlap_inter(spec_i = "HerringJuv",
                 spec_j = "FlounderAdu", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "HerringAdu",
                 spec_j = "FlounderJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "HerringJuv",
                 spec_j = "FlounderJuv", 
                 data   = df_inter)$overlap +
   overlap_inter(spec_i = "HerringAdu",
                 spec_j = "FlounderAdu", 
                 data   = df_inter)$overlap) / 4

inter["Herring", "Flounder"] <- inter["Flounder", "Herring"]

inter

# write.csv(data.frame(inter), "Params/inter.csv")


#======== D. PLOT DISTRIBUTION MAPS FOR EACH SPECIES ================================
# I will do this completely by species (life stages grouped). I restructure df_intra to get
# one column = one species and row is total abundance. Then I add lat and lon, and finally
# making data long again (one row = one observation) for easy ggplotting by species!
df_map <- df_intra %>% 
  gather(Life_stage, n, c("Adu", "Juv")) %>%
  dplyr::group_by(AreaCode, Species) %>% 
  dplyr::summarize(n = sum(n)) %>% 
  spread(Species, n)

df_map[is.na(df_map)] <- 0

head(df_map)

df_map$AreaCode2 <- df_map$AreaCode

df_map <- df_map %>% separate(AreaCode2, c("Area_s", "Code_s"), sep = 2)

df_map$Area_s <- as.numeric(df_map$Area_s)

head(df_map)

#====** Assign coordinates based on ICES rectangle ========
# Add latitude: calculate relation between latitude and rect area name
t <- data.frame(Area_s = head(sort(unique(df_map$Area_s))))
t$lat <- c(54.25, 54.75, 55.25, 55.75, 56.25, 56.75)
t$Area_cent <- t$Area_s - 36
summary(lm(lat ~ Area_cent, data = t))

df_map$lat <- 53.75 + (df_map$Area_s - 36)*0.5 # From summary

# Add longtitude (eh, don't think they codes actually mean anything so I'll just copy paste...)
df_map$lon <- 14.5
df_map$lon <- ifelse(df_map$Code_s == "G5", 15.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "G6", 16.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "G7", 17.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "G8", 18.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "G9", 19.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H0", 20.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H1", 21.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H2", 22.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H3", 23.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H4", 24.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H5", 25.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H6", 26.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H7", 27.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H8", 28.5, df_map$lon)
df_map$lon <- ifelse(df_map$Code_s == "H9", 29.5, df_map$lon)

#====** Plot maps ========
map1 <- fortify(map("world", 
                    resolution = 0, 
                    fill       = TRUE, 
                    bg         = "white", 
                    col        = "black", 
                    plot       = FALSE, 
                    xlim       = c(10,30.5), 
                    ylim       = c(53,70)))

# Calculate proportion in each rectangle, which we'll plot to standardize
df_map$Cod_p      <- df_map$Cod      / sum(df_map$Cod)
df_map$Flounder_p <- df_map$Flounder / sum(df_map$Flounder)
df_map$Herring_p  <- df_map$Herring  / sum(df_map$Herring)
df_map$Sprat_p    <- df_map$Sprat    / sum(df_map$Sprat)

# Reorganize data for plotting (long format, one row = one obs (rectangle))
p_dat <- df_map %>% 
  dplyr::select(AreaCode, lat, lon, Cod_p, Flounder_p, Herring_p, Sprat_p) %>% 
  dplyr::rename(Flounder = Flounder_p,
         Cod      = Cod_p,
         Sprat    = Sprat_p,
         Herring  = Herring_p) %>% 
  gather(key = Species, proportion, 4:7)

p_dat$proportion <- ifelse(p_dat$proportion == 0,
                           NA,
                           p_dat$proportion)

# Plot distribution maps, facet by species
p <- ggplot(data = p_dat, aes(x = lon, y = lat)) +
  geom_point(aes(size = proportion), color = viridis(n = 1, begin = 0.2)) + 
  geom_polygon(data = map1, aes(x = long, y = lat, group = group), 
               colour = "grey80", size = 0.9, fill = "grey80") +
  coord_cartesian(xlim = c(12, 30), ylim = c(54, 60.5)) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") + 
  scale_size_continuous(range = c(0.01, 4)) +
  #guides(size = FALSE) +
  facet_wrap(~ Species) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)+
  NULL

p

unique(p_dat$Species)

# Without flounder
p <- 
  p_dat %>% filter(Species %in% c("Cod", "Sprat", "Herring")) %>% 
  ggplot(., aes(x = lon, y = lat)) +
  geom_point(aes(size = proportion), color = viridis(n = 1, begin = 0.2)) + 
  geom_polygon(data = map1, aes(x = long, y = lat, group = group), 
               colour = "grey80", size = 0.9, fill = "grey80") +
  coord_cartesian(xlim = c(12, 30), ylim = c(54, 60.5)) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") + 
  scale_size_continuous(range = c(0.001, 3)) +
  #guides(size = FALSE) +
  facet_wrap(~ Species) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)+
  NULL

p

#ggsave("baltic/figures/supp/interaction_map.pdf", plot = p, dpi = 300, width = 15, height = 15, units = "cm")


# maybe manually add ICES rect? Through geom_segment?
