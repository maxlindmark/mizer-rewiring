# Define new function to get length at age
getGrowth <- function(object, species,
                      max_age = 20, percentage = FALSE, print_it = TRUE) {
  if (is(object, "MizerSim")) {
    sim <- object
    if (missing(species)) {
      species <- dimnames(sim@n)$sp
    }
    # reorder list of species to coincide with order in sim
    idx <- which(dimnames(sim@n)$sp %in% species)
    species <- dimnames(sim@n)$sp[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list("Species" = species, "Age" = age))
    g <- getEGrowth(sim@params, sim@n[dim(sim@n)[1], , ], 
                    sim@n_pp[dim(sim@n)[1], ], sim@n_bb[dim(sim@n)[1], ], sim@n_aa[dim(sim@n)[1], ], 
                    #sim@intTempScalar[,,1], sim@metTempScalar[,,1]) #AA
                    sim@intTempScalar[,,length(time_temperature_dt)], #ML
                    sim@metTempScalar[,, length(time_temperature_dt)]) #ML
    for (j in 1:length(species)) {
      i <- idx[j]
      g_fn <- stats::approxfun(sim@params@w, g[i, ])
      myodefun <- function(t, state, parameters){
        return(list(g_fn(state)))
      }
      ws[j, ] <- deSolve::ode(y = sim@params@species_params$w_min[i],
                              times = age, func = myodefun)[, 2]
      if (percentage) {
        ws[j, ] <- ws[j, ] / sim@params@species_params$w_inf[i] * 100
      }
    }
    plot_dat <- reshape2::melt(ws)
    return(plot_dat) # Added this to extract growth data in data.frame
  }    
}