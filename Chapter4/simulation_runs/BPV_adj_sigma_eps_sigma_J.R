# Sourcing the Jump test and functions ---------------------------------------------------

source("Jumps/BPV_adj.R")

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

intensity = 50/T_val 
m = 0

#Make grid

Boots=c(5,10,25,50,100,500,1000)
sigma_eps_grid = c(0, 0.00025, 0.0005, 0.001)
jumps_std_grid =  c(0, 0.0016, 0.0025, 0.0032)

Paramater_sets = expand.grid(sigma_eps = sigma_eps_grid, jump_std = jumps_std_grid)


# Monte Carlo -------------------------------------------------------------

M = 10000
B_cap <- 1000
K = 100

cores <-  30
cluster <-  makeCluster(cores)
clusterEvalQ(cluster, c(library(dplyr),
                        library(Matrix)))

clusterExport(cl = cluster, ls(), environment())

sim_result = lapply( X = 1:nrow(Paramater_sets), function(params_idx){
  
  sigma_eps <- Paramater_sets[params_idx,]$sigma_eps
  jump_std <- Paramater_sets[params_idx,]$jump_std
  
  #Find the best number of bootstraps:
  # Coverage <- parLapply(cl = cluster,X = 1:K, fun = function(idx){
  
  
  clusterExport(cl = cluster, c("sigma_eps","jump_std"), environment())
  
  Coverage <- pbapply::pblapply(cl = cluster,X = 1:K, FUN = function(idx){
    # Coverage <- lapply(X = 1:K, FUN = function(idx){
    single_boots <- lapply(Boots, function(bts){
     
      hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                             lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
      
      micro_noise <- sigma_eps * rnorm(n + 1)
      
      Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = jump_std)
      
      Path <- hest_sim$Heston + micro_noise + Jumps_sim
      
      #BPV_preAvg_estimate =BPV_preAvg(X = Path, sigma_eps = sigma_eps,  n = n, BstrapM = B_cap, BstrapSplit = bts, alpha = 0.05)
      BPV_adj_estimate=BPV_adj(Path,delta,BstrapM=B_cap,BstrapSplit=bts,0.05)
      IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
      
      BPV_adj_hit = (IV >= BPV_adj_estimate$lower_adj & IV <= BPV_adj_estimate$upper_adj) %>% data.frame
      colnames(BPV_adj_hit) <- bts
      BPV_adj_hit
    }) %>% do.call(what = cbind)
    single_boots
  }) %>% do.call(what = rbind) %>% colMeans()
  
  
     #browser()
  Best_Boot = which.min(abs(Coverage - 0.95)) %>% names %>% as.numeric()
  
  #Simulate with the best boot
  
  
  single_sim <- pbapply::pblapply(cl = cluster,X = 1:M, FUN = function(params){
     #single_sim <- lapply(X = 1:M, FUN = function(params){
      # browser()
    sigma_eps <- Paramater_sets[params_idx,]$sigma_eps
    jump_std <- Paramater_sets[params_idx,]$jump_std
    
    hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                           lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
    
    micro_noise <- sigma_eps * rnorm(n + 1)
    
    Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = jump_std)
    
    Path <- hest_sim$Heston + micro_noise + Jumps_sim
    
    BPV_adj_estimate=BPV_adj(Path,delta,BstrapM=B_cap,BstrapSplit=Best_Boot,0.05)
    IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
    
    data.frame(BPV_adj_estimate = abs(IV-BPV_adj_estimate$BPV_adj)/IV,
               BPV_adj_hit = (IV >= BPV_adj_estimate$lower_adj & IV <= BPV_adj_estimate$upper_adj))
  }) %>% do.call(what = rbind) %>% colMeans()
  
  return_single_sim <- data.frame(single_sim[1],single_sim[2], Best_Boot)
  rownames(return_single_sim) <- paste(Paramater_sets[params_idx,], collapse = "|")
  colnames(return_single_sim) <-c("BPV_adj","BPV_adj_Hits","Best Boot")
  
  total <- try(readRDS("BPV_adj_Sigma_eps_Sigma_J_Boot.rds"))
  if(class(total) == "try-error"){
    saveRDS(return_single_sim, "BPV_adj_Sigma_eps_Sigma_J_Boot.rds")
  }else{
    saveRDS(rbind(total, return_single_sim), "BPV_adj_Sigma_eps_Sigma_J_Boot.rds")
  }
  return_single_sim
}) %>% do.call(what = rbind)
stopCluster(cluster)


saveRDS(sim_result, file = "BPV_adj_Sigma_eps_Sigma_J_Boot.rds")


