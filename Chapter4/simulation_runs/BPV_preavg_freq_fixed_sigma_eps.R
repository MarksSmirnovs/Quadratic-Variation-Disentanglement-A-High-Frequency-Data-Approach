# Sourcing the Jump test and functions ---------------------------------------------------

source("Jumps/BPV_adj.R")

IV_fun <- function(vol_path, t, delta){
  eval_vals <- floor(t / delta)
  eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
  vol_path[eval_vals]
}




BPV_preAvg <- function(X, sigma_eps, n, BstrapM, BstrapSplit, alpha){
  
  n_Bstrap = function(X, BstrapSplit, delta){
    
    period_div <- split(X, ceiling(seq_along(X)/BstrapSplit))
    period_div <- Map(c, NA, period_div)
    Bstrapped_periods <- sample(length(period_div), replace = TRUE)
    Bstrapped_X <- period_div[Bstrapped_periods] %>% do.call(what = c)
    diff_X <- na.omit(diff(Bstrapped_X))
    n <- length(diff_X)
    (BPV(Bstrapped_X, delta = delta)$BPV) 
  }
  
  C <- switch(as.character(sigma_eps),
              "5e-04" = {1/19},
              "0.00025" = {1/19},
              "0.001" = {1/18},
              "0" = 0
  )
  
  M <- floor(C*n^(1/2))
  if(M != 0){ 
    X <- zoo::rollmean(X, M, align = "left")
    
    X <- X[seq(1, length(X), M)]
    delta <- T_val/length(X)
  }
  
  BPV_preavg <- BPV(X, delta = delta)$BPV
  
  Bstrap_res <- lapply(1:BstrapM, FUN = function(inpt){n_Bstrap(X, BstrapSplit = BstrapSplit, delta = delta)}) %>% do.call(what = c)
  
  list(
    "BPV_preavg" = BPV_preavg,
    "upper_adj" = quantile(Bstrap_res, probs = 1-alpha/2),
    "lower_adj" = quantile(Bstrap_res, probs = alpha/2)
  )
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
jump_std_grid = c(0, 0.0016, 0.0025, 0.0032)
n_grid = c(4680, 13000, 23400)
sigma_eps <- 0.0005

Paramater_sets = expand.grid(jump_std = jump_std_grid, n_freq = n_grid)


# Monte Carlo -------------------------------------------------------------

M = 10000
B_cap <- 1000
K = 100

# M = 10
# B_cap <- 1
# K = 1

cores <-  25
cluster <-  makeCluster(cores)
clusterEvalQ(cluster, c(library(dplyr),
                        library(Matrix)))

clusterExport(cl = cluster, ls(), environment())

sim_result = lapply( X = 1:nrow(Paramater_sets), function(params_idx){
  
  jump_std <- Paramater_sets[params_idx,]$jump_std
  n_sampl <- Paramater_sets[params_idx,]$n_freq
  
  #Find the best number of bootstraps:
  # Coverage <- parLapply(cl = cluster,X = 1:K, fun = function(idx){
  
  
  clusterExport(cl = cluster, c("jump_std", "n_sampl"), environment())
  
  Coverage <- pbapply::pblapply(cl = cluster,X = 1:K, FUN = function(idx){
    # Coverage <- lapply(X = 1:K, FUN = function(idx){
    single_boots <- lapply(Boots, function(bts){
      
      hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                             lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
      
      micro_noise <- sigma_eps * rnorm(n + 1)
      
      Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = jump_std)
      
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
      
      BPV_preAvg_estimate =BPV_preAvg(X = Path, sigma_eps = sigma_eps,  n = n, BstrapM = B_cap, BstrapSplit = bts, alpha = 0.05)
      IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
      
      BPV_preAvg_hit = (IV >= BPV_preAvg_estimate$lower_adj & IV <= BPV_preAvg_estimate$upper_adj) %>% data.frame
      colnames(BPV_preAvg_hit) <- bts
      BPV_preAvg_hit
      
    }) %>% do.call(what = cbind)
    single_boots
  }) %>% do.call(what = rbind) %>% colMeans()
  
  
  
  Best_Boot = which.min(abs(Coverage - 0.95)) %>% names %>% as.numeric()
  
  #Simulate with the best boot
  
  
  single_sim <- pbapply::pblapply(cl = cluster,X = 1:M, FUN = function(params){
    # single_sim <- lapply(X = 1:M, FUN = function(params){
    
    hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                           lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
    
    micro_noise <- sigma_eps * rnorm(n + 1)
    
    Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = jump_std)
    
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
    
    BPV_preAvg_estimate=BPV_preAvg(X = Path, n = n, sigma_eps = sigma_eps, BstrapM = B_cap, BstrapSplit = Best_Boot, alpha = 0.05)
    IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
    
    data.frame(BPV_preAvg = abs(IV-BPV_preAvg_estimate$BPV_preavg) /IV,
               BPV_preAvg_hit = (IV >= BPV_preAvg_estimate$lower_adj & IV <= BPV_preAvg_estimate$upper_adj))
  }) %>% do.call(what = rbind) %>% colMeans()
  
  return_single_sim <- data.frame(single_sim[1],single_sim[2], Best_Boot)
  rownames(return_single_sim) <- paste(Paramater_sets[params_idx,], collapse = "|")
  colnames(return_single_sim) <-c("BPV_preAvg","BPV_preAvg_Hits","Best Boot")
  
  total <- try(readRDS("BPV_preAvg_Sigma_J_n_freq.rds"))
  if(class(total) == "try-error"){
    saveRDS(return_single_sim, "BPV_preAvg_Sigma_J_n_freq.rds")
  }else{
    saveRDS(rbind(total, return_single_sim), "BPV_preAvg_Sigma_J_n_freq.rds")
  }
  return_single_sim
}) %>% do.call(what = rbind)
stopCluster(cluster)


saveRDS(sim_result, file = "BPV_preAvg_Sigma_J_n_freq_total.rds")