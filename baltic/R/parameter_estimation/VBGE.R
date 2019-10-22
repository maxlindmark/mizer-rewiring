#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2018.11.20: Max Lindmark
#
# - The following code fits VBGE curves to length-at-age data for species used in 
# mizer-implementation of the Baltic Sea (Lindmark et al., in prep)
# 
# A. Load libraries
#
# B. Fit VBGE curves by species to estimate K and L_inf for all species using Derek 
#    Ogles approach: http://derekogle.com/fishR/examples/oldFishRVignettes/VonBertalanffy.pdf
#      
# C. Plot fitted curves against empirical data 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Questions in progress marked with ***

#======== A. LOAD LIBRARIES =========================================================
rm(list = ls())

# Provide package names
pkgs <- c("ggplot2", 
          "dplyr", 
          "tidyr", 
          "FSA", 
          "FSAdata", 
          "nlstools", 
          "devtools", 
          "viridis")

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
# 3      FSA  0.8.21
# 4  FSAdata   0.3.6
# 5  ggplot2   3.0.0
# 6 nlstools   1.0-2
# 7    tidyr   0.8.1
# 8  viridis   0.5.1

#======== B. FIT VBGE BY SPECIES  ===================================================
#====** BITS (Trawl Survey) ========
# See cleanup_BITS_data for post download processing of data
dat <- read.csv("Data/BITS/clean_BITS.csv", sep = ",")

head(dat)

#====** Cod ========
cod <- dat %>% 
  filter(Species == "Cod") %>% 
  select("AgeRings", "IndWgt", "Species", "Length_cm", "Area_s", "Code_s")

head(cod)

ggplot(cod, aes(AgeRings, Length_cm)) +
  geom_point()

# Fit model
svTypical_c <- vbStarts(Length_cm ~ AgeRings, data = cod)

vbTypical <- Length_cm ~ Linf*(1-exp(-K*(AgeRings-t0)))

fitTypical_c <- nls(vbTypical, data = cod, start = svTypical_c)

fitPlot(fitTypical_c, xlab = "Age", ylab = "Total Length (cm)", main = "")

# Check model
overview(fitTypical_c)

# Check model assumptions
residPlot(fitTypical_c)

hist(residuals(fitTypical_c), main = "")

# Confidence intervals better acquired from bootstrapping in nls()
# bootTypical <- nlsBoot(fitTypical_c, niter = 200) 
# confint(bootTypical, plot = TRUE)

# Filter data for later plotting
cod$pred <- predict(fitTypical_c)

ggplot(cod, aes(AgeRings, Length_cm)) + 
  geom_point() +
  geom_line(aes(AgeRings, pred), col = "red", size = 2) +
  theme_classic(base_size = 20) +
  NULL

# Convert length to weight (see Appendix S1 & weight_length.R for paramters)
w_inf_cod <- (0.0078*summary(fitTypical_c)$coefficients[1]^3.07) # asymptotic weight
(0.0078*30^3.07) # weight at maturation


#====** Sprat ========
spr <- dat %>% 
  filter(Species == "Sprat") %>% 
  select("AgeRings", "IndWgt", "Species", "Length_cm", "Area_s", "Code_s")

head(spr)

ggplot(spr, aes(AgeRings, Length_cm)) +
  geom_point()

# Fit model
svTypical_s <- vbStarts(Length_cm ~ AgeRings, data = spr)

vbTypical <- Length_cm ~ Linf*(1-exp(-K*(AgeRings-t0)))

fitTypical_s <- nls(vbTypical, data = spr, start = svTypical_s)

fitPlot(fitTypical_s, xlab = "Age", ylab = "Total Length (cm)", main = "")

# Check model
overview(fitTypical_s)

# Check model assumptions
residPlot(fitTypical_s)

hist(residuals(fitTypical_s), main = "")

# Confidence intervals better acquired from bootstrapping in nls()
# bootTypical <- nlsBoot(fitTypical_s, niter = 200) 
# confint(bootTypical, plot = TRUE)

# Filter data for later plotting
spr$pred <- predict(fitTypical_s)

ggplot(spr, aes(AgeRings, Length_cm)) + 
  geom_point() +
  geom_line(aes(AgeRings, pred), col = "red", size = 2) +
  theme_classic(base_size = 20) +
  NULL

# Convert length to weight (see Appendix S1 & weight_length.R for paramters)
w_inf_spr <- (0.0041*summary(fitTypical_s)$coefficients[1]^3.15) # asymptotic weight
(0.0041*9^3.15) # weight at maturation


#====** Herring ========
her <- dat %>% 
  filter(Species == "Herring") %>% 
  select("AgeRings", "IndWgt", "Species", "Length_cm", "Area_s", "Code_s")

head(her)

ggplot(her, aes(AgeRings, Length_cm)) +
  geom_point()

# Fit model
svTypical_h <- vbStarts(Length_cm ~ AgeRings, data = her)

vbTypical <- Length_cm ~ Linf*(1-exp(-K*(AgeRings-t0)))

fitTypical_h <- nls(vbTypical, data = her, start = svTypical_h)

fitPlot(fitTypical_h, xlab = "Age", ylab = "Total Length (cm)", main = "")

# Check model
overview(fitTypical_h)

# Check model assumptions
residPlot(fitTypical_h)

hist(residuals(fitTypical_h), main = "")

# Confidence intervals better acquired from bootstrapping in nls()
# bootTypical <- nlsBoot(fitTypical_h, niter = 200) 
# confint(bootTypical, plot = TRUE)

# Filter data for later plotting
her$pred <- predict(fitTypical_h)

ggplot(her, aes(AgeRings, Length_cm)) + 
  geom_point() +
  geom_line(aes(AgeRings, pred), col = "red", size = 2) +
  theme_classic(base_size = 20) +
  NULL

# Convert length to weight (see Appendix S1 & weight_length.R for paramters)
w_inf_her <- (0.0042*summary(fitTypical_h)$coefficients[1]^3.14) # asymptotic weight
(0.0042*14^3.14) # weight at maturation


#======== C. PLOT FITTED CURVES AGAINST DATA ========================================
all_spec <- rbind(her, spr, cod)

#====** Length/weight ~ age ========
# Free y-axis
ggplot(all_spec, aes(AgeRings, Length_cm)) +
  facet_wrap(~ Species, ncol = 2, scales = "free_y") + 
  geom_point(size = 3, fill = "black", 
             color = "white", shape = 21, alpha = 0.05) + 
  geom_line(data = all_spec, aes(AgeRings, pred), size = 2, color = "red") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank()) +
  labs(x = "Age [yr]", y = "Length [cm]") +
  NULL

# Convert to weight-at-age
all_spec$Weight_g <- NA 

all_spec$Weight_g <- ifelse(all_spec$Species == "Herring",
                            0.0042*all_spec$Length_cm^3.14, 
                            all_spec$Weight_g)

all_spec$Weight_g <- ifelse(all_spec$Species == "Sprat",
                            0.0041*all_spec$Length_cm^3.15, 
                            all_spec$Weight_g)

all_spec$Weight_g <- ifelse(all_spec$Species == "Cod",
                            0.0078*all_spec$Length_cm^3.07, 
                            all_spec$Weight_g)

# Now convert predictions as well
all_spec$pred_g <- NA 

all_spec$pred_g <- ifelse(all_spec$Species == "Herring",
                          0.0042*all_spec$pred^3.14, 
                          all_spec$pred_g)

all_spec$pred_g <- ifelse(all_spec$Species == "Sprat",
                          0.0041*all_spec$pred^3.15, 
                          all_spec$pred_g)

all_spec$pred_g <- ifelse(all_spec$Species == "Cod",
                          0.0078*all_spec$pred^3.07, 
                          all_spec$pred_g)

ggplot(all_spec, aes(AgeRings, Weight_g)) +
  facet_wrap(~ Species, ncol = 2, scales = "free_y") + 
  geom_point(size = 3, fill = "black", 
             color = "white", shape = 21, alpha = 0.05) +
  geom_line(data = all_spec, aes(AgeRings, pred_g), size = 2, color = "red") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank()) +
  labs(x = "Age [yr]", y = "Weight [g]") +
  NULL

#====** Standardized length/weight ~ age ========
# Here I need to standardize against max (or asymptotic) weight so that I can use the same color density for all species
all_spec <- all_spec %>% 
  group_by(Species) %>% 
  mutate(weight_stand = Weight_g/max(Weight_g), 
         age_stand    = AgeRings/max(AgeRings),
         pred_stand   = pred_g/max(Weight_g))

# Plot standardized weight and age
p <- ggplot(all_spec, aes(age_stand, weight_stand)) +
  facet_wrap(~ Species, ncol = 2, scales = "free_y") + 
  geom_point(size = 3, fill = "black", 
             color = "white", shape = 21, alpha = 0.1) +
  geom_line(data = all_spec, aes(age_stand, pred_stand), 
            size = 2, color = viridis(n = 1, begin = 0.2)) +
  theme_bw(base_size = 12) +
  labs(x = "Age/max(Age)", y = "Weight/max(Weight)") +
  NULL

#ggsave("Figures/VBGE.tiff", plot = p, dpi = 300, width = 15, height = 15, units = "cm")

