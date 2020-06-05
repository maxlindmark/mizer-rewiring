# Read in get_time_elements in the global environment, can't get it to work from project_methods
get_time_elements <- function(sim, time_range, slot_name = "n"){
  if (!is(sim, "MizerSim"))
    stop("First argument to get_time_elements function must be of class MizerSim")
  time_range <- range(as.numeric(time_range))
  # Check that time range is even in object
  sim_times <- as.numeric(dimnames(sim@effort)$time)
  sim_time_range <- range(sim_times)
  if ( (time_range[1] < sim_time_range[1]) |
       (time_range[2] > sim_time_range[2]))
    stop("Time range is outside the time range of the model")
  time_elements <- (sim_times >= time_range[1]) & (sim_times <= time_range[2])
  names(time_elements) <- dimnames(sim@effort)$time
  return(time_elements)
}

# Create getSpectra function..
getSpectra <- function(object){
  time_range  <- max(as.numeric(dimnames(object@n)$time))
  time_elements <- get_time_elements(object, time_range)
  n <- apply(object@n[time_elements, , ,drop = FALSE], c(2, 3), mean)
  species <- object@params@species_params$species
  n <- n[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
  power <- 1
  n <- sweep(n, 2, params@w^power, "*")
  specDat <- data.frame(w = rep(as.numeric(dimnames(n)$w), length(species)),
                        n = c(as.numeric(n[1, ]), 
                              as.numeric(n[2, ]),
                              as.numeric(n[3, ])),
                        species = rep(species, each = 100))
  return(specDat)
}
