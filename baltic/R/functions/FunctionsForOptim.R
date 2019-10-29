# This code was provided by Asta Audzijonyte, here modified to my needs
# The error function has is modified code from Reum et al (2018) Oikos

# General params ===================================================================
# Number of size groups 
no_size_groups = 100

# Timestep used in the integration 
dt = 0.1

# How many years to run the model for each optimiser round.
tmax = 200 

# Fishing mortality for calibrations. This is overall mean for calibration period, see abundance_F_time series
effort = c(Cod = 0.899, Herring = 0.340, Sprat = 0.306)


# Optimization params ==============================================================
# Initialize count of function evaluations
optimizer_count = 0

# How many years to calculate base the error on
meansteps.par <- 30 

# Upper and lower bound functions ==================================================
# This function will take the upper multiplier values and set upper bounds
upper_bounds <- function(startpars) {  
  
  # Set an empty vector for upper values
  upper_value <- list()
  
  if (length(startpars[['r_max']]) > 0) {
    startpars.t <- startpars$r_max
    t_ind <- which(startpars$names == "r_max")
    upper_value.t <- startpars.t + log(startpars$multiplier[t_ind])
    upper_value.t[c(which(upper_value.t > log(startpars$upper_hard_bound[t_ind])))] <- log(startpars$upper_hard_bound[t_ind])
    upper_value$r_max <- upper_value.t
  }

    return(unlist(upper_value))
}

# This function will take the lower multiplier values and set upper bounds
lower_bounds <- function(startpars) {  
  
  # Set an empty vector for upper values
  lower_value <- list()
  
  if(length(startpars[['r_max']]) > 0) {
    startpars.t <- startpars$r_max
    t_ind <- which(startpars$names == "r_max")
    lower_value.t <- startpars.t - log(startpars$multiplier[t_ind])
    lower_value.t[c(which(lower_value.t < log(startpars$lower_hard_bound[t_ind])))] <- log(startpars$lower_hard_bound[t_ind])
    lower_value$r_max <- lower_value.t
  }
  
  return(unlist(lower_value))
}


# Error functions ==================================================================
# Berror_SSB_TOT: calculates the residual sum of square (RSS) based on SSB for species and model output
errorSSB <- function(model_run, meansteps = meansteps.par){
  
  # Grab 1992-2002 average SSB form the species params df in mizerParam object
  bio_obs <- model_run@params@species_params$AveSpawnBiomass
  
  # bio_obs_n <- bio_obs/bio_obs[1] # use this if you only want to minimize error of ssb relative to cod ssb
  
  # Take SSB in meansteps - dim() : dim
  bset_ssb <- getSSB(model_run)[c(I(dim(model_run@n)[1]-meansteps):dim(model_run@n)[1]), ] 
  Bhat_ssb <- colMeans(bset_ssb)
  
  # Bhat_ssb_n <- Bhat_ssb/Bhat_ssb[1] # use this if you only want to minimize error of ssb relative to cod ssb
  
  # Scale up from g/m^2 to g/Baltic(SD25-29+32) (2.49342E+11 m^2). Then change unit to 10^6 kg/Baltic
  bio_pred <- Bhat_ssb * model_run@params@species_params$sd25.29.32_m.2 / (1e9)
  
  # Get difference between observed and predicted biomass
  bdif <- log(bio_pred) - log(bio_obs)
  
  b_error <- sum(na.omit(bdif^2))
  
  return(b_error)
}


# Function to run model ============================================================

run_model <- function(params, tmax, effort) {
  
  model_run <- project(params, t_max = tmax, effort = effort, dt = dt)
  
  return (model_run)
  
}


# Main optimization function =======================================================
# The optim function uses fn function to calculate errors
# This function will take parameter values from the optim function, will take my mizer object params to run the model (run_model function), and a few other thing

calibratePar_Baltic <- function(start_vector, 
                                modelParams, 
                                startpars, 
                                effort = effort) {
  
  # Function arguments: 
  # startpars       - vectors of initial parameters for optimising 
  # modelParams     - my csv parameter file
  # meansteps       - to pass to the error functions
  # plus various multipliers to give different weights to different errors
  
  # Keep runing count of function evaluations by optim()
  optimizer_count <- optimizer_count + 1  
  
  assign("optimizer_count", optimizer_count, pos = .GlobalEnv)
  
  # Put values into the right param slots
  no_sp <- length(modelParams$species)
  count = 0 

  if(length(startpars[['r_max']]) > 0) {
  modelParams$r_max <- exp(start_vector[c((count+1):(count + no_sp))])
  count = count + length(startpars[['r_max']])
  print("updated r_max")
  print(count)
  }

  #---------------#
  print("one step")
  
  params <- MizerParams(modelParams, 
                        no_w = no_size_groups,
                        store_kernel = F,
                        kappa_ben = 1,
                        kappa = 1,
                        w_bb_cutoff = 20,
                        w_pp_cutoff = 1,
                        r_pp = 4,
                        r_bb = 4)
  
  # Increase maximum consumption rates by a factor 1.75
  h <- params@species_params$h
  params@species_params$h <- h * 1.75
  
  # Remove gamma, because it needs to be recalculated using the new h. 
  params@species_params <- subset(params@species_params,
                                  select = -gamma)
  
  # Create another mizerParam object
  params_upd <- MizerParams(params@species_params,
                            no_w = no_size_groups,
                            store_kernel = F,
                            kappa_ben = 1,
                            kappa = 1,
                            w_bb_cutoff = 20,
                            w_pp_cutoff = 1,
                            r_pp = 4,
                            r_bb = 4)
  
  params_upd@species_params
  
  print("two step")
  
  # Run the model for tmax years with the updated parameters 
  model_run <- run_model(params_upd, tmax, effort)
  
  #-----------------#
  print("three step")
  
  # Get ssb error (see definition above)
  errorSSB <- errorSSB(model_run, tmax)
  print("errorSSB")
  print(errorSSB)
  
  print(paste(optimizer_count, "errorSSB:", errorSSB))
  
  # Print values of optimised parameters to see if they are moving anywhere 
  print(start_vector)
  return(errorSSB)
}


