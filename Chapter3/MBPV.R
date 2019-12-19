# Source ------------------------------------------------------------------
source("functions.R")
source("Microstructure_noise_estimators/Heston_Sim.R")

# Helpfuns ----------------------------------------------------------------


MBPV <- function(Y, c_1, c_2, n, r, l, conf_alpha, T_val){

  #Parameters
  mu_1 <- sqrt(2)/sqrt(pi)   
  mu_2 <- 1           
  mu_r <- ifelse(r == 1, mu_1, ifelse(r == 2, 1, mean(abs(rnorm(100000000, 0 , 1))^r)))
  mu_l <- ifelse(l == 1, mu_1, ifelse(l == 2, 1, mean(abs(rnorm(100000000, 0 , 1))^l)))
  
  #Estimator
  MBV_bias <- function(K, M, Y, n, r,l){
    
    bar_Y = function(m, n, K, M, Y){
      
      (1/(floor(n/M) - K + 1)) * lapply(((m-1)*floor(n/M)):((m*floor(n/M)) - K), function(input){Y[input + K] -
          Y[input]}) %>% do.call(what = c) %>% sum
    }
    
    n^((r+l)/4 - 1/2) * lapply(2:M, function(input){(abs(bar_Y(m = input, n = n, K = K, M = M, Y = Y))^r)*
        (abs(bar_Y(m = input-1, n = n, K = K, M = M, Y = Y))^l)}) %>% do.call(what = c) %>% sum
  }
  
  if(missing(c_1) | missing(c_2)){
    c_1 <- 1
    c_2 <- 2
    K <- floor(c_1*n^(1/2))
    M <- floor(n^(1/2)/(c_1*c_2))
    IQ_pre <- MBV_bias(K = K, M = M, Y = Y, n = n, r = 4, l = 0)

    MBPV_constants <- function(param, mu_1, mu_2, IQ_pre){
      c_1 <- param[1]
      c_2 <- param[2]
    
      nu_1 <- (c_1*(3*c_2 - 4 + max((2-c_2)^3, 0)))/(3*(c_2-1)^2)

      loss <- ((c_1^2*c_2^2*(mu_2^2 + 2*mu_1^2*mu_2 - 3*mu_1^4))/(3*mu_1^4*nu_1^2))*IQ_pre

      if(loss < 0){1e1000}else{loss}
    }

    optim_pars <- optim(par = c(1,2), MBPV_constants,
                        mu_1 = mu_1, mu_2 = mu_2, IQ_pre = IQ_pre,
                        lower = c(0+1e-7, 1+1e-7), method = "L-BFGS-B")

    c_1 <- optim_pars$par[1]
    c_2 <- optim_pars$par[2]

    
  }
  
  K <- floor(c_1*n^(1/2))
  M <- floor(n^(1/2)/(c_1*c_2))
  nu_1 <- (c_1*(3*c_2 - 4 + max((2-c_2)^3, 0)))/(3*(c_2-1)^2)
  nu_2 <- (2*(min((c_2-1), 1)))/(c_1*(c_2-1)^2)
  
  
  #Noise estimate
  micro_noise <- (1/(2*n))*sum(diff(Y)^2)
  
  #MBPV estimate
  MBPV <-  (((c_1*c_2)/(mu_r*mu_l))*MBV_bias(K = K, M = M, Y = Y, n = n, r = r, l = l) - T_val*nu_2*micro_noise)/(nu_1)
  
  #Confidence intervals
  beta_sqrd <- ((c_1^2*c_2^2*(mu_2^2 + 2*mu_1^2*mu_2 - 3*mu_1^4))/(3*mu_1^4*nu_1^2)) * MBV_bias(K = K, M = M, Y = Y, n = n, r = 4, l = 0)
  
  list(
    "MBPV" = MBPV,
    "Beta" = sqrt(beta_sqrd),
    "upper" = MBPV + qnorm(1 - conf_alpha/2)*sqrt(beta_sqrd)*n^(-1/4),
    "lower" = MBPV + qnorm(conf_alpha/2)*sqrt(beta_sqrd)*n^(-1/4)
    #"MRV" = MRV
  )
}

#   -----------------------------------------------------------------------

r = 0.05
alpha = 0.04*5
lambda = 5
sigma_v = 0.5
rho = -0.5
S_0 = 1
V_0 = 0.3
n = 23400
sigma_eps <-  0.0005
T_val = 1/252
delta = T_val/n
# 
# # set.seed(123)
# Sim_noise = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
#                        lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
# 
# # 
# # ### IV estimate
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
# # 
# # start <- now()
# # MBPV_estimate <- MBPV(Y = Sim_noise$Heston + sigma_eps*rnorm(n + 1),
# #                       n = n, r = 1, l = 1, conf_alpha = 0.05, T_val = T_val, c_1 = 1.5, c_2 = 3);MBPV_estimate
# # now()-start
# # 
# # MBPV_estimate$MBPV
# # IV
# # 
# # #   -----------------------------------------------------------------------
# {
# cores <-  3
# cluster <-  makeCluster(cores)
# clusterEvalQ(cluster, c(library(dplyr)))
# clusterExport(cluster, ls(), environment())
# 
# 
# z1 = pbapply::pblapply(cl = cluster, X = 1:10, FUN = function(input){
# 
#   z1 <- Heston_Sim(T_val= T_val, n = n , r = r, rho = rho , alpha = alpha,
#                    lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
#   z2 <- MBPV_estimate <- MBPV(Y = z1$Heston + sigma_eps * rnorm(n + 1),
#                               n = n, r = 1, l = 1, conf_alpha = 0.05, T_val = T_val);MBPV_estimate
# 
#   z3 <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = z1$Vol, delta = delta, subdivisions = 10000000)
# 
#   cat("\r", input)
# 
#   data.frame(
#     IV = z3$value,
#     BMPV = z2$MBPV,
#     abs = abs(z2$MBPV - z3$value),
#     upper = z2$upper,
#     lower = z2$lower,
#     conf = z2$upper >= z3$value & z3$value  >= z2$lower
#   )
# 
# })
# 
# stopCluster(cluster)
# 
# z1 <- z1 %>% do.call(what = rbind)
# 
# }
# z1$conf %>% mean
# z1$abs %>% mean
