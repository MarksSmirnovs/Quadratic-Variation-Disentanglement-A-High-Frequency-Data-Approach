# Source ------------------------------------------------------------------
source("functions.R")
source("Microstructure_noise_estimators/Heston_Sim.R")

# Helpfuns ----------------------------------------------------------------

BPV <- function(X, alpha = 0.05, delta){
  
  diff_X <- na.omit(diff(X))
  n <- length(diff_X)

  BPV <- (pi/2)*sum(abs(diff_X[2:n]) * abs(diff_X[1:(n-1)]))
  QPV <- delta^(-1) * ((sqrt(2)/sqrt(pi))^(-4))*sum(abs(diff_X[4:n]) * abs(diff_X[3:(n-1)]) *
                                                      abs(diff_X[2:(n-2)]) * abs(diff_X[1:(n-3)]))
  BPV_sigma <- ((pi^2)/4 + pi - 3)*QPV
  
  return_ls <- list()
  return_ls$BPV <- BPV
  return_ls$upper <- BPV + qnorm(1-alpha/2)*sqrt(BPV_sigma)*sqrt(delta)
  return_ls$lower <- BPV + qnorm(alpha/2)*sqrt(BPV_sigma)*sqrt(delta)
  
  return_ls
}

BPV_adj <- function(X, delta, BstrapM, BstrapSplit, alpha){
  BPV_adj <- BPV(X, delta = delta)$BPV - 2* (length(X)-1) * (sum(diff(X)^2)/(2*(length(X)-1)))
  
  n_Bstrap = function(X, BstrapSplit, delta){
    
    period_div <- split(X, ceiling(seq_along(X)/BstrapSplit))
    period_div <- Map(c, NA, period_div)
    Bstrapped_periods <- sample(length(period_div), replace = TRUE)
    Bstrapped_X <- period_div[Bstrapped_periods] %>% do.call(what = c)
    diff_X <- na.omit(diff(Bstrapped_X))
    n <- length(diff_X)
    (BPV(Bstrapped_X, delta = delta)$BPV - 2 * n *
      (sum(diff_X^2)/(2*n))) 
  }
  
  Bstrap_res <- lapply(1:BstrapM, FUN = function(inpt){n_Bstrap(X, BstrapSplit = BstrapSplit, delta = delta)}) %>% do.call(what = c)
  
  list(
    "BPV_adj" = BPV_adj,
    "upper_adj" = quantile(Bstrap_res, probs = 1-alpha/2),
    "lower_adj" = quantile(Bstrap_res, probs = alpha/2)
  )
}


#   -----------------------------------------------------------------------

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
# 
# 
# IV_fun <- function(vol_path, t, delta){
#   eval_vals <- floor(t / delta)
#   eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
#   vol_path[eval_vals]
# }
# 
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
#     z2 <- BPV(X = z1$Heston + z1$Noise, delta = delta)
#     z22 <- BPV_adj(X = z1$Heston + z1$Noise, delta = delta,
#                    BstrapM = 1000, BstrapSplit = 100, alpha = 0.05)
#     
#     z3 <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = z1$Vol,
#                     delta = delta, subdivisions = 10000000)
#     
#     data.frame(
#       IV = z3$value,
#       BPV = z2$BPV,
#       upper = z2$upper,
#       lower = z2$lower,
#       conf = z2$upper >= z3$value & z3$value  >= z2$lower,
#       BPV_adj = z22$BPV_adj,
#       upper_adj = z22$upper_adj,
#       lower_adj = z22$lower_adj,
#       conf_adj = z22$upper_adj >= z3$value & z3$value  >= z22$lower_adj
#     )
#     
#   })
#   
#   stopCluster(cluster)
#   
#   z1 <- z1 %>% do.call(what = rbind)
#   
#   print(z1$conf %>% mean)
#   print(z1$conf_adj %>% mean)
# }
