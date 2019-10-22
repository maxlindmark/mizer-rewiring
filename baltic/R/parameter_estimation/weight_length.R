#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2018.11.20: Max Lindmark
#
# - The following code estimates length-weight relationships for the time period 
# used for calibration (1992-2002) using trawl survey data
#
# A. Load libraries
#
# B. Fit log-log model to get allometric length weight model for each species 
#      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Questions in progress marked with ***

# Largely following the approach in here: http://derekogle.com/fishR/examples/oldFishRVignettes/LengthWeight.pdf


#======== A. LOAD LIBRARIES =========================================================
rm(list = ls())

# Provide package names
pkgs <- c("ggplot2", "dplyr", "tidyr", "FSA", "devtools")

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
# 3      FSA  0.8.22
# 4  ggplot2   3.0.0
# 5    tidyr   0.8.1


#======== B. Fit allometric length-weight model =====================================
#====** BITS (Trawl Survey) ========
# See cleanup_BITS_data for post download processing of data
dat <- read.csv("Data/BITS/clean_BITS.csv", sep = ",")

head(dat)

#====** Cod ========
cod <- dat %>% 
  filter(Species == "Cod" & IndWgt > 0) %>% 
  select("AgeRings", "IndWgt", "Species", "Length_cm", "Area_s", "Code_s")

head(cod)
str(cod)

#====** Fit model ========
ggplot(cod, aes(Length_cm, IndWgt)) +
  geom_point(size = 4, alpha = 0.2) +
  theme_classic(base_size = 18) +
  NULL

cod$log_Length_cm <- log(cod$Length_cm)
cod$log_IndWgt <- log(cod$IndWgt)

c1 <- lm(log_IndWgt ~ log_Length_cm, data = cod)

fitPlot(c1, 
        xlab = "log Total Length (cm)", 
        ylab = "log Weight (g)")

summary(c1)

# Print coefficients 
coefficients(summary(c1))[2] # a
exp(coefficients(summary(c1))[1]) # b


#====** Flounder ========
flo <- dat %>% 
  filter(Species == "Flounder" & IndWgt > 0) %>% 
  select("AgeRings", "IndWgt", "Species", "Length_cm", "Area_s", "Code_s")

head(flo)
str(flo)

#====** Fit model ========
ggplot(flo, aes(Length_cm, IndWgt)) +
  geom_point(size = 4, alpha = 0.2) +
  theme_classic(base_size = 18) +
  NULL

flo$log_Length_cm <- log(flo$Length_cm)
flo$log_IndWgt <- log(flo$IndWgt)

f1 <- lm(log_IndWgt ~ log_Length_cm, data = flo)

fitPlot(f1, 
        xlab = "log Total Length (cm)", 
        ylab = "log Weight (g)")

summary(f1)

# Print coefficients 
coefficients(summary(f1))[2] # a
exp(coefficients(summary(f1))[1]) # b


#====** Sprat ========
spr <- dat %>% 
  filter(Species == "Sprat" & IndWgt > 0) %>% 
  select("AgeRings", "IndWgt", "Species", "Length_cm", "Area_s", "Code_s")

head(spr)
str(spr)

#====** Fit model ========
ggplot(spr, aes(Length_cm, IndWgt)) +
  geom_point(size = 4, alpha = 0.2) +
  theme_classic(base_size = 18) +
  NULL

spr$log_Length_cm <- log(spr$Length_cm)
spr$log_IndWgt <- log(spr$IndWgt)

s1 <- lm(log_IndWgt ~ log_Length_cm, data = spr)

fitPlot(s1, 
        xlab = "log Total Length (cm)", 
        ylab = "log Weight (g)")

summary(s1)

# Print coefficients 
coefficients(summary(s1))[2] # a
exp(coefficients(summary(s1))[1]) # b


#====** Herring ========
her <- dat %>% 
  filter(Species == "Herring" & IndWgt > 0) %>% 
  select("AgeRings", "IndWgt", "Species", "Length_cm", "Area_s", "Code_s")

head(her)
str(her)

#====** Fit model ========
ggplot(her, aes(Length_cm, IndWgt)) +
  geom_point(size = 4, alpha = 0.2) +
  theme_classic(base_size = 18) +
  NULL

her$log_Length_cm <- log(her$Length_cm)
her$log_IndWgt <- log(her$IndWgt)

h1 <- lm(log_IndWgt ~ log_Length_cm, data = her)

fitPlot(h1, 
        xlab = "log Total Length (cm)", 
        ylab = "log Weight (g)")

summary(h1)

# Print coefficients 
coefficients(summary(h1))[2] # a
exp(coefficients(summary(h1))[1]) # b
