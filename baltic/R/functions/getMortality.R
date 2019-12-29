# Create getMortality function..
getMortality <- function(object){
m2_time <- getPredMort(object, time_range = dim(object@effort)[1], intakeScalar = object@intTempScalar[,,dim(object@effort)[1]], drop = FALSE)
m2 <- apply(m2_time, c(2, 3), mean)
m2 <- m2[as.character(dimnames(m2)[[1]]) %in% object@params@species_params$species, , 
         drop = FALSE]

mortDat <- data.frame(value = c(m2),
                      Species = dimnames(m2)[[1]],
                      w = rep(object@params@w, each = length(object@params@species_params$species)))

mortDat <- mortDat %>% arrange(Species)

return(mortDat)
}
