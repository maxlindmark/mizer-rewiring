#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.01.02: Max Lindmark
#
# Code for fitting linear model of ln(r_max) ~ Arrhenius temp based on selected data
# from Savage et al 2004 AmNat and Barnes  to get activation energy
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())

dat <- read.csv("baltic/data/Savage_2004_r_max.csv")

head(dat)

m <- lm(ln_r ~ arrhenius_temp, data = dat)

summary(m)

# Call:
#   lm(formula = ln_r ~ arrhenius_temp, data = dat)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.20687 -0.31051  0.00112  0.32142  1.43704 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)     27.9822     4.0508   6.908 2.52e-08 ***
#   arrhenius_temp  -0.7346     0.1012  -7.260 8.18e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5759 on 40 degrees of freedom
# Multiple R-squared:  0.5685,	Adjusted R-squared:  0.5577 
# F-statistic:  52.7 on 1 and 40 DF,  p-value: 8.177e-09

## CI for slope:

# > 0.7346 - 2*0.1012 
# [1] 0.5322
# > 0.7346 + 2*0.1012 
# [1] 0.937

par(mfrow = c(2,2))
plot(m)


## Barnes
dat <- read.csv("baltic/data/barnes_intercept.csv")

head(dat)

colnames(dat)[1] <- "temp"
colnames(dat)[2] <- "intercept" # log10 of intercept

dat$normal_intercept <- 10^(dat$intercept)

dat$ln_intercept <- log(dat$normal_intercept)

dat$temp_arr <- 1/((dat$temp + 273.15) * 8.617332e-05)

m <- lm(intercept ~ temp_arr, data = dat) # log10 fit
summary(m)

m <- lm(ln_intercept ~ temp_arr, data = dat)
summary(m)

dat$pred <- predict(m)

# Call:
#   lm(formula = ln_intercept ~ temp_arr, data = dat)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.6051 -0.8654 -0.2320  0.6818  4.2818 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -10.0041     4.0980  -2.441   0.0154 *  
#   temp_arr      0.7953     0.1020   7.799 2.66e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.406 on 215 degrees of freedom
# Multiple R-squared:  0.2205,	Adjusted R-squared:  0.2169 
# F-statistic: 60.83 on 1 and 215 DF,  p-value: 2.665e-13

## CI for slope:
-0.7953 - 2*0.1020
#[1] -0.9993
-0.7953 + 2*0.1020
#[1] -0.5913

par(mfrow = c(2,2))
plot(m)
par(mfrow = c(1,1))

# What do I do with that??? Very small variation in the slope..

plot(dat$temp_arr, dat$ln_intercept)
lines(dat$temp_arr, dat$pred, col = "red")

library(ggplot2)
ggplot(dat, aes(temp_arr, ln_intercept)) +
  geom_point() +
  stat_smooth(method = "lm") +
  NULL



