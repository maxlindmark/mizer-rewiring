#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  # This code was provided by Asta Audzijonyte, here modified to my needs
#    To use this file:
#
# A. Read in your parameter data
#
# B. Comment out the lines in the starpars$ for those parameters that you don't want
#    to calibrate. *For now I only have r_max (removed the others from Asta)
#
# C. Decide on the multipliers for the parameter range to be explored: e.g. explore 
#    gamma values that are x2 times bigger and smaller than the original values. Here 
#    the multiplier would be 2. This sets the upper and lower bounds for the parameter 
#    values to explore
#
# D. Decide on the complete hard upper and lower bounds that would apply across all 
#    species, because using the multiplier for some parameters can produce unrealistic 
#    values, e.g. erepro of 0.7 multiplied by 2 would give 1.4. So you could set hard 
#    bounds for erepro to be below 1. 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. Read data
modelParams <- balticParams

# B. Create list to hold parameters to be optimized
startpars <- list()

# Which parameters do you want to optimise
startpars$names <- c("erepro")

# Log scale!
startpars$erepro <- log(modelParams$erepro) 

# C. Define r_max multiplier
startpars$multiplier <- 1000

# D. Define hard bounds (lower and upper)
startpars$upper_hard_bound <- 1

startpars$lower_hard_bound <- 1e-6

# Lastly, use the list to create a vector of parameters to be optimised. This is because optim takes a vector, not a list
start_vector <- c(startpars$erepro)
print(start_vector)
