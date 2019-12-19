
# Sourcing the estimators ---------------------------------------------------

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


cores <- 4
cluster <- makeCluster(cores)
clusterEvalQ(cluster, c(library(dplyr)))
clusterExport(cl = cluster, ls(), environment())


# jump_std_grid <- c(0.00115, seq(0.00075, 0.0015, length.out = 5))
jump_std_grid <- c(0.0032)

result <- c()
for(i in 1:length(jump_std_grid)){
  jumps_std = jump_std_grid[i]
  clusterExport(cluster, "jumps_std", environment())
  
  sim_res <- pbapply::pblapply(X = 1:1000, cl = cluster, FUN = function(inpt){
    hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                           lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
    
    Jumps_sim_s <- Compound_Pois(T_val = T_val , n = n + 1,
                                 intensity = intensity , m = m, sigma_J = jumps_std)
    
    IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
    
    
    iv_ratio_s <- sum(diff(Jumps_sim_s)^2)/(IV + sum(diff(Jumps_sim_s)^2)) #%>% mean 
    
    
    
    data.frame(
      iv_ratio_s = iv_ratio_s
    )
  })
  
  single <- sim_res %>% do.call(what = rbind) %>%
    colMeans() %>% round(5) %>% data.frame
  colnames(single) <- paste(jumps_std)
  result <- c(result, single)
}


stopCluster(cluster)
print(result)
