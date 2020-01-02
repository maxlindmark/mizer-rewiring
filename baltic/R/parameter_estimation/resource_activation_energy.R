#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.01.02: Max Lindmark
#
# Code for fitting linear model of ln(r_max) ~ Arrhenius temp based on selected data
# from Savage et al 2004 AmNat to get activation energy
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

