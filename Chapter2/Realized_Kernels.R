
# Setup -------------------------------------------------------------------

source("functions.R")

# Realized Kernel Estimation ----------------------------------------------

Realized_Kernel_Estimator = function(X, m = 2, K = 1200, alpha = 0.05, debug = F){
  if(debug){browser()}
  
  n = length(X)
  
  #Defineing help functions ----
  
  #Making a function for the Parzen kernel
  Parzen_kernel = function(x){
    
    if( 0 <= x & x< 0.5 ) k = 1 - 6*x^2 + 6*x^3
    if( 0.5 <= x & x <= 1) k = 2*(1-x)^3
    if( x > 1) k = 0
    
    return(k)
    
  }
  
  #Making the lagged return covariance function
  gamma_h = function(X,h){
    returns = diff(X)
   
    l = length(returns)
    
    gamma = sum( returns[(h+1):l]*returns[1:(l-h)] ) 
    
    return(gamma)
  }
  
  #End effects ----
  if(m > 0){
    
    X_0 = mean(X[1:m])
    
    X_T = mean(X[(n-m+1):n])
    
    X = c(X_0, X[-c(1:m,(n-m+1):n)] , X_T)
    
    n = n-2*(m-1)
    
  }
  
  #Finding the optimal bandwidth ----
  
  c_star = ((12^2)/0.269)^(1/5)
  
  IV_hat = IV_avg(X,K)
  
  sigma_eps = (sum(diff(X)^2))/(2*n)
  
  zhi = sqrt(sigma_eps) / sqrt(IV_hat)
  
  H = c_star*zhi^(0.8)*n^(0.6)
  
  #Computing the realized kernel estimate ----
  
  Index = which( c(1:n) / (H + 1) <= 1)
  
  RKs = lapply( c(0,Index) , function(h){
    
    RK = ifelse(h != 0,
                2*Parzen_kernel(h / (H+1 ) ) * gamma_h(X,h),
                Parzen_kernel(h / (H+1 ) ) * gamma_h(X,h) ) 
    
  })
  
  RK_IV = do.call(sum,RKs)
  
  #Estimating the variance of the error ----
  kappa = (12*(0.269^2))^(1/5)
  
  sigma_IQ = (sqrt(sigma_eps)*IV_hat^2)^(2/5)
  
  RK_IV_variance = 4*(kappa*sigma_IQ)^2
 
  
  #Making alpha level confidence intervals ----
  Lower = RK_IV + n^(-1/5)*kappa*sigma_IQ + n^(-1/5)*qnorm(alpha/2)*sqrt(RK_IV_variance)
  Upper = RK_IV + n^(-1/5)*kappa*sigma_IQ + n^(-1/5)*qnorm(1-alpha/2)*sqrt(RK_IV_variance)
  
  
  
  
  return(list(Estimate = RK_IV,
              ErrorVariance = RK_IV_variance,
              lower = Lower,
              upper = Upper))
}

# set.seed(42)
# Realized_Kernel_Estimator(rnorm(10000), debug = F)


# Testing -----------------------------------------------------------------

# source("Microstructure_noise_estimators/Heston_Sim.R")
# 
# r = 0.05
# alpha = 0.04*5
# lambda = 5
# sigma_v = 0.5
# rho = -0.5
# S_0 = 1
# V_0 = 0.3
# n = 1e6
# sigma_eps <-  0.0005
# T_val = 1/252
# delta = T_val/n
# # 
# # 
# # Count = 0
# # for (i in 1:100) {
# #   Sim_noise = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
# #                          lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)
# #   
# #   IV = Sim_noise$Vol %>% sum()*delta; IV
# #   
# #   #RV = diff(Sim_noise$Stock)^2 %>% sum(); RV
# #   
# #   Kernel = Realized_Kernel_Estimator(X = Sim_noise$Stock , debug = F); Kernel
# #   
# #   
# #   Count = ifelse( between(IV,Kernel$lower,Kernel$upper), Count +1 , Count)
# #   print(i)
# # }
# 
# IV_fun <- function(vol_path, t, delta){
#   eval_vals <- floor(t / delta)
#   eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
#   vol_path[eval_vals]
# }
# 
# {
#   cores <-  3
#   cluster <-  makeCluster(cores)
#   clusterEvalQ(cluster, c(library(dplyr)))
#   clusterExport(cluster, ls(), environment())
#   
#   
#   z1 = pbapply::pblapply(cl = cluster, X = 1:100, FUN = function(input){
#     
#     z1 <- Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
#                      lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)
#     z2 <- Realized_Kernel_Estimator(X = z1$Stock , debug = F);
#     
#     z3 <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = z1$Vol, delta = delta, subdivisions = 10000000)
#     
#     cat("\r", input)
#     
#     data.frame(
#       IV = z3$value,
#       RK = z2$Estimate,
#       upper = z2$upper,
#       lower = z2$lower,
#       conf = z2$upper >= z3$value & z3$value  >= z2$lower
#     )
#     
#   })
#   
#   stopCluster(cluster)
#   
#   z1 <- z1 %>% do.call(what = rbind)
#   
#   z1$conf %>% mean
# }
