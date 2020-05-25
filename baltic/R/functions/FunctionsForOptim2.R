# This code was provided by Asta Audzijonyte, here modified to my needs
# The error function has is modified code from Reum et al (2018) Oikos

# General params ===================================================================
# Number of size groups 
no_size_groups = 100

# Timestep used in the integration 
dt = 0.1

# How many years to run the model for each optimiser round.
t_max = 200 

# Reference temperature (rescaled projection, see calibration.v1)
t_ref <- 10.11562

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
  
  if (length(startpars[['erepro']]) > 0) {
    startpars.t <- startpars$erepro
    t_ind <- which(startpars$names == "erepro")
    upper_value.t <- startpars.t + log(startpars$multiplier[t_ind])
    upper_value.t[c(which(upper_value.t > log(startpars$upper_hard_bound[t_ind])))] <- log(startpars$upper_hard_bound[t_ind])
    upper_value$erepro <- upper_value.t
  }

    return(unlist(upper_value))
}

# This function will take the lower multiplier values and set upper bounds
lower_bounds <- function(startpars) {  
  
  # Set an empty vector for upper values
  lower_value <- list()
  
  if(length(startpars[['erepro']]) > 0) {
    startpars.t <- startpars$erepro
    t_ind <- which(startpars$names == "erepro")
    lower_value.t <- startpars.t - log(startpars$multiplier[t_ind])
    lower_value.t[c(which(lower_value.t < log(startpars$lower_hard_bound[t_ind])))] <- log(startpars$lower_hard_bound[t_ind])
    lower_value$erepro <- lower_value.t
  }
  
  return(unlist(lower_value))
}


# Error functions ==================================================================
# Berror_FMSY_TOT: calculates the residual sum of square (RSS) based on FMSY for species and model output
errorFMSY <- function(model_run, meansteps = meansteps.par){
  
  # Grab 1992-2002 average SSB form the species params df in mizerParam object
  FMSY_obs <- model_run@params@species_params$aveFMSY
  
  # bio_obs_n <- bio_obs/bio_obs[1] # use this if you only want to minimize error of ssb relative to cod ssb
  
  # Loop through fishing mortalities to find FMSY
              
              F_range <- seq(0, 1.2, 0.1) # Can increase later, becomes too slow now
              t_max_loop <- 200            
  
              # Cod
              # Mean F in calibration time
              effort = c(Cod = balticParams$AveEffort[1], Herring = balticParams$AveEffort[3], Sprat = balticParams$AveEffort[2])
              
              # Create empty data holder
              codFmsy <- c(); Y <- c(); Fm <- c(); t <- c()
              
              for(i in F_range) {
                effort[1] <- i
                t <- project(model_run, dt = 0.1, effort = effort, temperature = rep(model_run@t_ref, t_max_loop), diet_steps = 10, t_max = t_max_loop) 
                Y <- mean(data.frame(getYield(t))$Cod[(t_max-20):t_max])
                Fm <- i
                t <- cbind(Y, Fm)
                codFmsy <- data.frame(rbind(t, codFmsy))
              }
              
              codFmsy$Species <- "Cod"
  
              # Sprat
              # Mean F in calibration time
              effort = c(Cod = balticParams$AveEffort[1], Herring = balticParams$AveEffort[3], Sprat = balticParams$AveEffort[2])
              
              # Create empty data holder
              sprFmsy <- c(); Y <- c(); Fm <- c(); t <- c()
              
              for(i in F_range) {
                effort[3] <- i
                t <- project(model_run, dt = 0.1, effort = effort, temperature = rep(model_run@t_ref, t_max_loop), diet_steps = 10, t_max = t_max_loop) 
                Y <- mean(data.frame(getYield(t))$Sprat[(t_max-20):t_max])
                Fm <- i
                t <- cbind(Y, Fm)
                sprFmsy <- data.frame(rbind(t, sprFmsy))
              }

              sprFmsy$Species <- "Sprat"

              # Herring
              # Mean F in calibration time
              effort = c(Cod = balticParams$AveEffort[1], Herring = balticParams$AveEffort[3], Sprat = balticParams$AveEffort[2])
              
              # Create empty data holder
              herFmsy <- c(); Y <- c(); Fm <- c(); t <- c()
              
              for(i in F_range) {
                effort[2] <- i
                t <- project(model_run, dt = 0.1, effort = effort, temperature = rep(model_run@t_ref, t_max_loop), diet_steps = 10, t_max = t_max_loop) 
                Y <- mean(data.frame(getYield(t))$Herring[(t_max-20):t_max])
                Fm <- i
                t <- cbind(Y, Fm)
                herFmsy <- data.frame(rbind(t, herFmsy))
              }
              
              herFmsy$Species <- "Herring"
              
              # All together
              Fmsy <- rbind(codFmsy, sprFmsy, herFmsy)
              
              # Find FMSY
              fmsy <- Fmsy %>% dplyr::group_by(Species) %>% filter(Y == max(Y))
              

  # Get difference between observed and predicted FMSY
  bdif <- log(fmsy) - log(model_run@params@species_params$aveFMSY)
  
  b_error <- sum(na.omit(bdif^2))
  
  return(b_error)
}


# Function to run model ============================================================

run_model <- function(params, t_max, effort) {
  
  model_run <- project(params, t_max = t_max, effort = effort, dt = dt)
  
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

  if(length(startpars[['erepro']]) > 0) {
  modelParams$erepro <- exp(start_vector[c((count+1):(count + no_sp))])
  count = count + length(startpars[['erepro']])
  print("updated erepro")
  print(count)
  }

  #---------------#
  print("one step")
  
  params <- MizerParams(modelParams, 
                        no_w = no_size_groups,
                        store_kernel = F,
                        kappa_ben = 9,
                        kappa = 9,
                        w_bb_cutoff = 20,
                        w_pp_cutoff = 1,
                        r_pp = 4,
                        r_bb = 4,
                        t_ref = t_ref)
  
  # Recalculate ks so that it is 0.12*h and not 0.2*h (old hardwired bug...)
  params@species_params$ks <- params@species_params$h * 0.12
  
  # Increase maximum consumption rates by a factor (or species-specific factor)
  h <- params@species_params$h
  params@species_params$h <- h * 1.3
  
  # Remove gamma, because it needs to be recalculated using the new h. 
  params@species_params <- subset(params@species_params,
                                  select = -gamma)
  
  # Create another mizerParam object
  params_upd <- MizerParams(params@species_params,
                            no_w = no_size_groups,
                            store_kernel = F,
                            kappa_ben = 9,
                            kappa = 9,
                            w_bb_cutoff = 20,
                            w_pp_cutoff = 1,
                            r_pp = 4,
                            r_bb = 4,
                            t_ref = t_ref)
  
  params_upd@species_params
  
  print("two step")
  
  # Run the model for t_max years with the updated parameters 
  model_run <- run_model(params_upd, t_max, effort)
  
  #-----------------#
  print("three step")
  
  # Get ssb error (see definition above)
  errorSSB <- errorSSB(model_run, t_max)
  print("errorSSB")
  print(errorSSB)
  
  print(paste(optimizer_count, "errorSSB:", errorSSB))
  
  # Print values of optimised parameters to see if they are moving anywhere 
  print(start_vector)
  return(errorSSB)
}


