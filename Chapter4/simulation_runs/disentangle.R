source("Microstructure_noise_estimators/TSRV.R")
source("Microstructure_noise_estimators/Pre_average.R")
source("Microstructure_noise_estimators/Realized_Kernels.R")
source("Microstructure_noise_estimators/Heston_Sim.R")

source("Jumps/BPV_adj.R")
source("Jumps/MBPV.R")
source("Jumps/Heston_Poisson.R")
source("functions.R")

IV_fun <- function(vol_path, t, delta){
  eval_vals <- floor(t / delta)
  eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
  vol_path[eval_vals]
}


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

intensity = 50/T_val
m = 0
sigma_eps_grid <- c(0, 0.00025, 0.0005, 0.001)
sigma_J_grid = c(0, 0.0016, 0.0025, 0.0032)
Paramater_sets = expand.grid(sigma_eps = sigma_eps_grid, sigma_J = sigma_J_grid)

M = 10000

cores <- 10
cluster <- makeCluster(cores)
clusterEvalQ(cluster, c(library(dplyr)))
clusterExport(cluster, ls(), environment())

Montes = lapply(1:nrow(Paramater_sets), function(params_idx){
  single_sim <- pbapply::pblapply(cl = cluster, X = 1:M, FUN = function(idx){
    
    sigma_eps <- Paramater_sets[params_idx,]$sigma_eps
    sigma_J <- Paramater_sets[params_idx,]$sigma_J
    
    hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                           lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
    noise <- sigma_eps * rnorm(n + 1)
    
    jumps <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = sigma_J)
    
    path <- hest_sim$Heston + noise + jumps
    
    IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
    
    true_jumps <- sum(diff(unique(jumps))^2)
    
    IV_MBPV <- MBPV(path, T_val = T_val, r = 1, l = 1, n = n, conf_alpha = 0.05)$MBPV
    QV_Kernel <- Realized_Kernel_Estimator(X = path)$Estimate
    QV_TSRV <- TSRV_estimator_optimal(X = path, 5, 1, 3, T_val)$TSRV
    QV_PreAvg <- PreAvg_Estimator(X = path, t = T_val, delta = delta, k_n = 51)$Estimate
    
    if(sigma_eps == 0.001){
      C <- 1/18
    }else if(sigma_eps == 0){
      C <- 0
    }else{
      C <- 1/19
    }
    M <- floor(C*n^(1/2))
    
    if(M != 0){
      preAvgPath <- zoo::rollmean(path, M-1, align = "left")
      
      preAvgPath <- preAvgPath[seq(1, length(preAvgPath), M)]
      delta <- T_val/length(preAvgPath)
    }else{
      preAvgPath = path
    }
    
    
    IV_BPV <- BPV(preAvgPath, delta = delta)$BPV
    
    if(sigma_eps == 0.001){
      substr <- IV_MBPV
    }else if(sigma_eps == 0.00025 & Sigma_J %in% c(0, 0.0016) ){
      substr <- IV_MBPV
    }else{
      substr <- IV_BPV
    }
    
    RV_preavg <- sum(diff(preAvgPath)^2)
    
    data.frame(
      "Pre_Average_Sample_MBPV" = abs((RV_preavg - IV_MBPV) - true_jumps)/true_jumps,
      "Pre_Average_Sample_BPV" = abs((RV_preavg - IV_BPV) - true_jumps)/true_jumps,
      "Kernel" = abs((QV_Kernel - IV_MBPV) - true_jumps)/true_jumps,
      "TSRV" = abs((QV_TSRV - IV_MBPV) - true_jumps)/true_jumps,
      "PreAvg" = abs((QV_PreAvg - IV_MBPV) - true_jumps)/true_jumps
    )
  }) %>% do.call(what = rbind) %>% colMeans()
  
  return_single_sim <- data.frame(single_sim)
  colnames(return_single_sim) <- paste("sigma_eps = ", Paramater_sets[params_idx,]$sigma_eps, "|sigma_J = ", Paramater_sets[params_idx,]$sigma_J, sep = "")
  return_single_sim
}) %>% do.call(what = cbind)

stopCluster(cluster)


saveRDS(Montes, "QV_sim.rds")