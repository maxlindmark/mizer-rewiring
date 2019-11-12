# Create getSpectra function..
getSpectra <- function(object){
  time_range  <- max(as.numeric(dimnames(object@n)$time))
  time_elements <- get_time_elements(object, time_range)
  n <- apply(object@n[time_elements, , ,drop = FALSE], c(2, 3), mean)
  species <- balticParams$species
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
