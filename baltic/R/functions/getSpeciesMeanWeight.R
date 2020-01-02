# Create getSpectra function..
getSpeciesMeanWeight <- function(sim, species = 1:nrow(sim@params@species_params), ...){
  check_species(sim, species)
  n_species <- getN(sim, ...)
  biomass_species <- getBiomass(sim, ...)
  #n_total <- apply(n_species[, species, drop = FALSE], 1, sum)
  #biomass_total <- apply(biomass_species[, species, drop = FALSE], 1, sum)
  #return(biomass_total / n_total)
  return(biomass_species / n_species)
}
