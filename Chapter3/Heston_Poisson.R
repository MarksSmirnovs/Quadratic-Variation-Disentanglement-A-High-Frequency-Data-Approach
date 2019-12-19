# Setup -------------------------------------------------------------------

library(dplyr)

source("Microstructure_noise_estimators/Heston_Sim.R")


# Compound Poisson --------------------------------------------------------

Compound_Pois = function(T_val, n , intensity, m , sigma_J){
  
  t = seq(0,T_val , length.out = n)
  
  #Number of jumps
  N = rpois(1 , intensity*T_val)
  
  #Distribution of jumps
  U = runif(n = N, min = 0 , max = T_val)
  
  #Simulation of jumps
  Jumps = rnorm(N , m , sigma_J)
  
  X = rep(0,n)
  
  for(i in 1:n) {
    
    Jump_Index = which(U < t[i])
    
    X[i] = sum(Jumps[Jump_Index])
  }
  
  X
  
}

#Testing:

#Parameters
# T = 250
# n = 1000
# x = 100
# 
# #Jumps
# intensity = 0.1
# m = 0.001
# delta = 0.004
# 
# #Simulation
# set.seed(1)
# Compound_Pois(T = T , n = n , intensity = intensity , m = m, delta = delta ) %>% plot(type = "l")


# Heston with jumps -------------------------------------------------------

Heston_Jumps_Sim = function(T_val, n , r , rho , alpha, lambda, sigma_v, S_0 , V_0, sigma_eps, intensity, m , sigma_J){

  #Simulating the Heston model plus noise:
  Heston_sim = Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)

  #Simulating the jumps:
  Jumps_sim = Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = sigma_J )

  Heston_sim$Jumps = Jumps_sim

  return(Heston_sim)

}

#Testing:

# Test = Heston_Jumps_Sim(T_val = 1/252,
#                  r = 0.05,
#                  alpha = 0.04*5,
#                  lambda = 5,
#                  sigma_v = 0.5,
#                  rho = -0.5,
#                  S_0 = 1,
#                  V_0 = 0.3,
#                  n = 23400,
#                  sigma_eps <-  0.0005,
#                  intensity = 0.014,
#                  m = 0,
#                  sigma_J = 0.0045)
# 
# par(mfrow=c(4,1))
# 
# plot(Test$Heston, type = "l")
# 
# plot(Test$Noise, type = "l")
# 
# plot(Test$Jumps, type = "l")
# 
# plot(Test$Heston + Test$Noise + Test$Jumps, type = "l")
# 
# par(mfrow=c(1,1))
