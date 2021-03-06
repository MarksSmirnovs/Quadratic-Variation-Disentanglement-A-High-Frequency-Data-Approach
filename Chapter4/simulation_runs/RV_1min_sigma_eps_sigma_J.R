
# Sourcing the estimators ---------------------------------------------------

source("Microstructure_noise_estimators/TSRV.R")
source("Microstructure_noise_estimators/Pre_average.R")
source("Microstructure_noise_estimators/Realized_Kernels.R")
source("Microstructure_noise_estimators/Heston_Sim.R")

source("Jumps/BPV_adj.R")
source("Jumps/MBPV.R")

IV_fun <- function(vol_path, t, delta){
  eval_vals <- floor(t / delta)
  eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
  vol_path[eval_vals]
}

# Parameters --------------------------------------------------------------

source("Jumps/Heston_Poisson.R")

r = 0.05
alpha = 0.04*5
lambda = 5
sigma_v = 0.5
rho = -0.5
S_0 = 1
V_0 = 0.3
n = 23400
T_val = 1/252
delta = T_val/n

m = 0
intensity =  50/T_val

#Make grid
sigma_eps_grid <- c(0, 0.00025, 0.0005, 0.001)
sigma_J_grid = c(0, 0.0016, 0.0025, 0.0032)

Paramater_sets = expand.grid(sigma_eps = sigma_eps_grid, sigma_J = sigma_J_grid)

# Monte Carlo --------------------------------------------------------------

M = 10000

cores <-  10
cluster <-  makeCluster(cores)
clusterEvalQ(cluster, c(library(dplyr)))
clusterExport(cl = cluster, ls(), environment())

# sim_result = pbapply::pblapply(cl = cluster, X = 1:nrow(Paramater_sets), function(params_idx){
sim_result = lapply(1:nrow(Paramater_sets), function(params_idx){
  
  # single_sim <- lapply(X = 1:M, FUN = function(params){
  single_sim <- pbapply::pblapply(cl = cluster, X = 1:M, FUN = function(params){
    
    sigma_eps <- Paramater_sets[params_idx,]$sigma_eps
    sigma_J <- Paramater_sets[params_idx,]$sigma_J
    
    hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                           lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
    
    micro_noise <- sigma_eps * rnorm(n + 1)
    
    Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = sigma_J)
    
    Path <- hest_sim$Heston + micro_noise + Jumps_sim
    
    n_seq <- floor(seq(1, length(Path), length.out = 390 + 1))
    Path <- Path[n_seq]
    n <- length(Path) - 1
    delta <-  T_val/n
    
    RV <- sum(diff(Path)^2)
    RV_conf_std <- sqrt(2/3 * sum(diff(Path)^4))
    RV_conf_upper <- RV + qnorm(0.975) * RV_conf_std
    RV_conf_lower <- RV + qnorm(0.025) * RV_conf_std
    
    IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
    
    data.frame(RV = abs(RV - IV)/IV, RV_hit = (IV >= RV_conf_lower & IV <= RV_conf_upper)
    )
  }) %>% do.call(what = rbind) %>% colMeans()
  
  return_single_sim <- data.frame(single_sim)
  colnames(return_single_sim) <- paste("sigma_eps = ", Paramater_sets[params_idx,]$sigma_eps, "|sigma_J = ", Paramater_sets[params_idx,]$sigma_J, sep = "")
  return_single_sim
}) %>% do.call(what = cbind)
stopCluster(cluster)

saveRDS(sim_result, "sim_study_RV_1min.rds")