
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
intensity = 50/T_val

#Make grid
n_grid = c(4680, 13000, 23400)
sigma_eps_grid = c(0, 0.00025, 0.0005, 0.001)
sigma_J <- 0.0016

Paramater_sets = expand.grid(n_sampl = n_grid, sigma_eps = sigma_eps_grid)

# Monte Carlo --------------------------------------------------------------

M = 10000

cores <-  25
cluster <-  makeCluster(cores)
clusterEvalQ(cluster, c(library(dplyr)))
clusterExport(cl = cluster, ls(), environment())

# sim_result = pbapply::pblapply(cl = cluster, X = 1:nrow(Paramater_sets), function(params_idx){
sim_result = lapply(1:nrow(Paramater_sets), function(params_idx){
  
  # single_sim <- lapply(X = 1:M, FUN = function(params){
  single_sim <- pbapply::pblapply(cl = cluster, X = 1:M, FUN = function(params){
    
    n_sampl <- Paramater_sets[params_idx,]$n_sampl
    sigma_eps <- Paramater_sets[params_idx,]$sigma_eps
    
    hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                           lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
    
    micro_noise <- sigma_eps * rnorm(n + 1)
    
    Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = sigma_J)
    
    Path <- hest_sim$Heston + micro_noise + Jumps_sim
    
    if(n_sampl != 13000){
      n_seq <- floor(seq(1, length(Path), length.out = n_sampl + 1))
      Path <- Path[n_seq]
      n <- length(Path) - 1
      delta <-  T_val/n
    }else{
      NA_places <- sample(x = 2:length(Path), size = 10400, replace =  FALSE)
      Path[NA_places] <- NA
      Path <- zoo::na.locf(Path, fromLast = FALSE)
      n <- length(Path) - 1
      delta <-  T_val/n
    }
    
    TSRV <- tryCatch(TSRV_estimator_optimal(X = Path, 5, 1, 3, T_val), error = function(cond) list(TSRV = NA, upper = NA, lower = NA))
    RV <- sum(diff(Path)^2)
    RV_conf_std <- sqrt(2/3 * sum(diff(Path)^4))
    RV_conf_upper <- RV + qnorm(0.975) * RV_conf_std
    RV_conf_lower <- RV + qnorm(0.025) * RV_conf_std
    
    IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
    
    data.frame(RV = abs(RV - IV)/IV, RV_hit = (IV >= RV_conf_lower & IV <= RV_conf_upper),
               TSRV = abs(TSRV$TSRV - IV)/IV, TSRV_hit = (IV >= TSRV$lower & IV <= TSRV$upper)
    )
  }) %>% do.call(what = rbind) %>% colMeans()
  
  return_single_sim <- data.frame(single_sim)
  colnames(return_single_sim) <- paste("n_sampl = ", Paramater_sets[params_idx,]$n_sampl, "|sigma_J = ", Paramater_sets[params_idx,]$sigma_J, sep = "")
  return_single_sim
}) %>% do.call(what = cbind)
stopCluster(cluster)

saveRDS(sim_result, "tsrv_corrected_coverage_sigma_eps_freq.rds")