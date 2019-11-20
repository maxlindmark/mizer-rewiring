tempFun <- function(temperature, t_ref, Ea, c_a, w) # default are 0 for now as deactivation is buggy
{
  # tempFun returns a matrix with w (size) as columns and temperature as rows
  
  # equation
  # (w^(c_a*(temperature-t_ref)))  *exp((-Ea/8.617332e-5)*((1/temperature) - (1/t_ref)))
  # *(1/(w^(c_d*(temperature-t_d)))*exp((-Ed/8.617332e-5)*((1/temperature) - (1/t_d))))
  
  temperature <- temperature + 273 # converting to Kelvin from Celcius
  t_ref <- t_ref + 273
  
  temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature-(t_ref)))}) *exp((-Ea/8.617332e-5)*((1/temperature) - (1/(t_ref))))) 
  
  return(temperatureScalar)
}

# tempFun(temperature = c(5:15), t_ref = 10, Ea = 0, c_a = 0, w = c(1:10))
# tempFun(temperature = c(5:15), t_ref = 10, Ea = 0.63, c_a = 0, w = c(1:10))
# tempFun(temperature = c(5:15), t_ref = 10, Ea = 0.63, c_a = 0, w = c(1))
# tempFun(temperature = c(5:15), t_ref = 10, Ea = 0.63, c_a = -0.02, w = c(1:10))




