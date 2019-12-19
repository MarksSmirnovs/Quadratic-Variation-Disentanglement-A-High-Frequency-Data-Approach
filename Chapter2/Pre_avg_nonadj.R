PreAvg_Estimator <- function(X, t, delta, k_n){
  
  g_input <- function(x){min(x,1-x)}
  
  g <- function(i,k_n, g){g(i/k_n)}
  
  Delta_X <- function(i, X){X[i] - X[i-1]}
  
  bar_X <- function(i, k_n, g, g_select, Delta_X, X){
    
    summation <- 0
    for(idx in 2:(k_n-1)){
      summation <- summation + g(idx, k_n, g_select) * Delta_X(i + idx, X)   
    }
    summation
  }
  
  #Finite sample
  
  psi_1 <- 1
  psi_2 <- 1/12
  PHI_11 = 1/6
  PHI_12 =  1/96
  PHI_22 = 151/80640
  
  #PreAvg estimate
  


    
    bias_est <- function(k_n, g, g_input, Delta_X, X, t, delta, power){
      summation <- 0
      for(idx in 1:(round(t/delta) - k_n + 1)){
        summation <- summation + bar_X(idx, k_n, g, g_input, Delta_X, X)^power
      }
      summation
    }
    
    bias <- sum(diff(X)[1:(round(t/delta)-1)]^2)
    
    
    theta <- k_n * sqrt(delta)
    
   estimate <-  (sqrt(delta)/(theta*psi_2)) * bias_est(k_n, g, g_input, Delta_X, X, t, delta, power = 2) -
      ((psi_1*delta)/(2*theta^2*psi_2)) * bias
    
  
  
  bias_cross <- function(X, k_n, g, g_input, Delta_X){
    diff_X <- c(NA, diff(X))
    summation <- 0
    for(idx in 1:(round(t/delta - 2*k_n + 1))){
      
      new_term <- bar_X(idx, k_n, g, g_input, Delta_X, X)^2
      summation <- summation + new_term * sum(diff_X[(idx + k_n - 1):(idx + 2*k_n - 2)]^2)
    }
    summation
  }
  
  diff_x <- diff(X)
 
  variance <- ((4*PHI_22)/(3*theta*psi_2^4)) * bias_est(k_n, g, g_input, Delta_X, X, t, delta, power = 4) +
    ((4*delta)/(theta^3)) * ((PHI_12/(psi_2)^3) - (PHI_22*psi_1)/(psi_2^4)) * bias_cross(X, k_n, g, g_input, Delta_X) +
    (delta/(theta^3)) * ((PHI_11/(psi_2)^2) - 2 * ((PHI_12*psi_1)/(psi_2)^3) + ((PHI_22*psi_1^2)/(psi_2)^4)) *
    sum((diff_x[2:(round(t/delta)-3)]^2)*(diff_x[4:(round(t/delta) - 1)]^2))
  
  
  
  
  list(
    "Estimate" = estimate,
    "variance" = variance,
    "+-sigma" = sqrt(variance) * delta^(1/4),
    "True IV 95 % interval" = paste("[", round(estimate + qnorm(0.025)*sqrt(variance)*delta^(1/4),2),
                                    ",",  round(estimate + qnorm(0.975)*sqrt(variance)*delta^(1/4),2),  "]", sep = "")
  )
  
}

#   -----------------------------------------------------------------------

source("functions.R")
source("Microstructure_noise_estimators/Heston_Sim.R")


r = 0.05
alpha = 0.5
lambda = 0.5
sigma_v = 0.4
rho = -0.7
S_0 = 1
V_0 = 0.3
n = 100000
sigma_eps <-  0.0005
T_val = 1
delta = T_val/n

set.seed(123)
Sim_noise = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
                       lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)

### IV estimate
IV_fun <- function(vol_path, t, delta){
  eval_vals <- floor(t / delta)
  eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
  vol_path[eval_vals]
}

integrate(IV_fun, lower = 0, upper = T_val, vol_path = Sim_noise$Vol, delta = delta, subdivisions = 10000000)
###

diff(Sim_noise$Stock)^2 %>% sum
PreAvg_Estimator(X = Sim_noise$Stock, t = T_val, delta = delta, k_n = 51)


new_idx <- floor(seq(1,n, length.out = 30000))
new_x <- Sim_noise$Stock[new_idx]
new_delta <- T_val / length(new_x)

diff(new_x)^2 %>% sum
PreAvg_Estimator(X = new_x, t = T_val, delta = new_delta, k_n = 51)
