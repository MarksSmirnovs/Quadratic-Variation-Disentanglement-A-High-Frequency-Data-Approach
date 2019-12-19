PreAvg_Estimator <- function(X, t, delta, k_n){

g_input <- function(x){min(x,1-x)}

g <- function(i,k_n, g){g(i/k_n)}

Delta_X <- function(i, X){X[i] - X[i-1]}

bar_X <- function(i, k_n, g, g_select, Delta_X, X){
 lapply((2:(k_n-1)), function(idx){g(idx, k_n, g_select) * Delta_X(i + idx, X)}) %>% do.call(what = cbind) %>% rowSums
}
  
#Finite sample

psi_1 <- k_n*lapply(1:k_n, function(idx){(g(idx + 1, k_n, g_input) - g(idx, k_n, g_input))^2}) %>% do.call(what = sum)


psi_2 <- (1/k_n)*lapply(1:(k_n - 1), function(idx){g(idx, k_n, g_input)^2}) %>% do.call(what = sum)


phi_1 <- function(i, k_n, g, g_input){
  
  summation <- 0
  for(idx in (i+1):(k_n)){
    summation <- summation + (g((i-1), k_n, g_input)  - g(i, k_n, g_input))*(g((i-idx-1), k_n, g_input)  - g((i-idx), k_n, g_input))
  }
  summation
}

phi_2 <- function(i, k_n, g, g_input){
  summation <- 0
  for(idx in (i+1):(k_n)){
    summation <- summation + g((i), k_n, g_input)*g((i-idx), k_n, g_input)
  }
  summation
}


PHI_11 <- function(k_n, g, g_input, phi_1){
  summation <- 0
  for(idx in 1:(k_n)){
    summation <- summation + (phi_1(idx, k_n, g, g_input)^2)
  }
  k_n*(summation - 1/2*(phi_1(0, k_n, g, g_input)^2))
}

PHI_12 <- function(k_n, g, g_input, phi_1, phi_2){
  summation <- 0
  for(idx in 1:(k_n)){
    summation <- summation + (phi_1(idx, k_n, g, g_input)*phi_2(idx, k_n, g, g_input))
  }
  (1/k_n)*(summation - 1/2*(phi_1(0, k_n, g, g_input)*phi_2(0, k_n, g, g_input)))
}

PHI_22 <- function(k_n, g, g_input, phi_1, phi_2){
  summation <- 0
  for(idx in 1:(k_n)){
    summation <- summation + (phi_2(idx, k_n, g, g_input)^2)
  }
  (1/k_n^3)*(summation - 1/2*(phi_2(0, k_n, g, g_input)^2))
}


#PreAvg estimate
theta <- k_n * sqrt(delta)

biased_est <- function(t, delta, k_n, Delta_X, X, bar_X){
  sum(bar_X((1:(floor(t/delta) - k_n + 1)), k_n, g, g_input, Delta_X, X)^2)
}

bias_correct <- function(t, delta, k_n, Delta_X, X, bar_X){
  sum(Delta_X((2:(floor(t/delta))), X)^2)
}

estimate <- (1 - (psi_1 * delta)/(2*theta^2*psi_2))^(-1) * 
  (((floor(t/delta) * sqrt(delta))/((floor(t/delta) - k_n + 2)*theta*psi_2))*biased_est(t, delta, k_n, Delta_X, X, bar_X) -
      ((psi_1*delta)/(2*theta^2*psi_2)) * bias_correct(t, delta, k_n, Delta_X, X, bar_X)
    )


#PreAvg variance
biased_est_bar <- function(t, delta, k_n, Delta_X, X, g, g_select){
  sum(bar_X(1:(floor(t/delta) - k_n + 1), k_n, g, g_select, Delta_X, X)^4)
}

biased_est_bar_double <- function(t, delta, k_n, Delta_X, X, g, g_select){
  diff_X <- c(NA,diff(X))
  sum(bar_X(1:(floor(t/delta) - 2*k_n + 1), k_n, g, g_select, Delta_X, X)^2 * 
  lapply(1:(floor(t/delta) - 2*k_n + 1), function(idx){sum(diff_X[(idx + k_n-1):(idx + 2*k_n - 1)]^2)}) %>%
    do.call(what = rbind)
  )
}

bias_cross_prod <- function(t, delta, k_n, Delta_X, X){
  sum(Delta_X((2:(floor(t/delta) - 2)), X)^2 * Delta_X((2:(floor(t/delta) - 2)) + 2, X)^2)
}
  
  
variance <- (1 - (psi_1 * delta)/(2*theta^2*psi_2))^(-2) * (
    ((4 * PHI_22(k_n, g, g_input, phi_1, phi_2) * floor(t/delta))/
       (3*theta*(psi_2)^4*(floor(t/delta) - k_n + 2))) *
      biased_est_bar(t, delta, k_n, Delta_X, X, g, g_input) + 
      ((4*delta*floor(t/delta))/(theta^3*(floor(t/delta) - k_n + 2))) * 
      ((PHI_12(k_n, g, g_input, phi_1, phi_2)/(psi_2)^3) -
         ((PHI_22(k_n, g, g_input, phi_1, phi_2)*psi_1 )/
            (psi_2)^4)) * biased_est_bar_double(t, delta, k_n, Delta_X, X, g, g_input) + 
      ((delta*floor(t/delta))/(theta^3*(floor(t/delta) - 2)))*
         ((PHI_11(k_n, g, g_input, phi_1)/(psi_2)^2) -
            ((2*PHI_12(k_n, g, g_input, phi_1, phi_2)*psi_1 )/
               (psi_2)^3) +
            ((PHI_22(k_n, g, g_input, phi_1, phi_2)*psi_1^2)/
               (psi_2)^4)) * bias_cross_prod(t, delta, k_n, Delta_X, X))

#Output
list(
  "Estimate" = estimate,
  "Estimate Variance" = variance,
  "+-sigma" = sqrt(variance)*delta^(1/4),
  "True IV 95 % interval" = paste("[", round(estimate + qnorm(0.025)*sqrt(variance)*delta^(1/4),2),
                                  ",",  round(estimate + qnorm(0.975)*sqrt(variance)*delta^(1/4),2),  "]", sep = ""),
  "lower" = (estimate + qnorm(0.025)*sqrt(variance)*delta^(1/4)),
  "upper" = (estimate + qnorm(0.975)*sqrt(variance)*delta^(1/4))
)
}

#   -----------------------------------------------------------------------

# source("functions.R")
# source("Microstructure_noise_estimators/Heston_Sim.R")
# 
# 
# r = 0.05
# alpha = 0.04*5
# lambda = 5
# sigma_v = 0.5
# rho = -0.5
# S_0 = 1
# V_0 = 0.3
# n = 23400
# sigma_eps <-  0.0005
# T_val = 1/252
# delta = T_val/n
# 
# # set.seed(123)
# # Sim_noise = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
# #                        lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)
# # 
# 
# ### IV estimate
# IV_fun <- function(vol_path, t, delta){
#   eval_vals <- floor(t / delta)
#   eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
#   vol_path[eval_vals]
# }
# #   
# # IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = Sim_noise$Vol, delta = delta, subdivisions = 10000000);IV
# # ###
# # 
# # diff(Sim_noise$Stock)^2 %>% sum
# # start <- now()
# # Pre_Avg <- PreAvg_Estimator(X = Sim_noise$Stock, t = T_val, delta = delta, k_n = 51);Pre_Avg
# # now()-start
# # 
# # abs(IV$value - Pre_Avg$Estimate)
# # 
# # new_idx <- floor(seq(1,n, length.out = 50000))
# # new_x <- Sim_noise$Stock[new_idx]
# # new_delta <- T_val / length(new_x)
# # 
# # diff(new_x)^2 %>% sum
# # PreAvg_Estimator(X = new_x, t = T_val, delta = new_delta, k_n = 51)
# #  
# 
# #   -----------------------------------------------------------------------
# 
# {
#   cores <-  3
#   cluster <-  makeCluster(cores)
#   clusterEvalQ(cluster, c(library(dplyr)))
#   clusterExport(cluster, ls(), environment())
#   
#   
#   z1 = pbapply::pblapply(cl = cluster, X = 1:1000, FUN = function(input){
#     
#     z1 <- Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
#                      lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)
#     z2 <- PreAvg_Estimator(X = z1$Stock, t = T_val, delta = delta, k_n = 51)
#     
#     z3 <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = z1$Vol, delta = delta, subdivisions = 10000000)
#     
#     cat("\r", input)
#     
#     data.frame(
#       IV = z3$value,
#       BMPV = z2$Estimate,
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
