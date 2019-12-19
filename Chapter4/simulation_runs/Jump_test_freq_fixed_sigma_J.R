# Sourcing the Jump test and functions ---------------------------------------------------

source("Jumps/jump_test_v2.R")
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
T_val = 1/252
n=23400
delta = T_val/n

J_int = 12600 
sigma_eps = 0.0005 
sigma_J=0.0016
m = 0
C=1/19

#Make grid

#C values are chosen from Lee and Mykland table 5 page 402

sigma_eps_grid=c(0,0.00025,0.0005,0.001)
#sigma_J_grid = c(0, 0.0016,0.0025,0.0032)
n_grid=c(4680,13000,23400)

C_grid=data.frame(C=c(rep(1/19,9),rep(1/18,3)))

 Paramater_sets = expand.grid(n = n_grid,
                              sigma_eps = sigma_eps_grid
 )%>%cbind(C_grid)




# Monte Carlo --------------------------------------------------------------

M = 10000

cores <-  12
cluster <-  makeCluster(cores)
clusterEvalQ(cluster, c(library(dplyr),
                        library(Matrix)
)
)
clusterExport(cl = cluster, ls(), environment())

sim_result = pbapply::pblapply(cl = cluster, X = 1:nrow(Paramater_sets), function(params_idx){
  #sim_result = lapply(1:nrow(Paramater_sets), function(params_idx){
  #browser()
  single_sim <- lapply(X = 1:M, FUN = function(params){
    #browser()
    sigma_eps <- Paramater_sets[params_idx,]$sigma_eps
    C  <- Paramater_sets[params_idx,]$C
    #sigma_J=Paramater_sets[params_idx,]$sigma_J
    n_sampl=13000
    n_sampl=Paramater_sets[params_idx,]$n
    
    hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                           lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)
    
    micro_noise <- sigma_eps * rnorm(n + 1)
    
    Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = J_int , m = m, sigma_J = sigma_J)

    Path <- hest_sim$Heston + micro_noise + Jumps_sim
    
    if(n_sampl != 13000){
      n_seq <- floor(seq(1, length(Path), length.out = n_sampl + 1))
      Path <- Path[n_seq] 
      delta <-  T_val/(length(Path) - 1)
    }else{
      NA_places <- sample(x = 2:length(Path), size = 10400, replace =  FALSE)
      Path[NA_places] <- NA
      Path <- zoo::na.locf(Path, fromLast = FALSE)
      delta <-  T_val/(length(Path) - 1)
    }
    
    
    #BPV_adj_estimate=BPV_adj(Path,delta,1000,Boot,0.05)
    IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
    Hit_est=Jump_test(C,Path,0.05,delta)$Rejection
    Hit_real=Jump_test2(C,Path,0.05,delta,IV,sigma_eps)$Rejection
    # ,
    # BPV_adj = abs(BPV_adj_estimate$BPV_adj - IV)/IV,
    # RV_hit = (IV >= BPV_adj_estimate$lower_adj & IV <= BPV_adj_estimate$upper_adj)
    
    data.frame(Hit_est= Hit_est,
               Hit_real=Hit_real)
  }) %>% do.call(what = rbind) %>% colMeans()
  
  return_single_sim <- data.frame(single_sim[1],single_sim[2])
  rownames(return_single_sim) <- paste("sigma_eps = ", Paramater_sets[params_idx,]$sigma_eps,
                                       "|n_sampl= ",  Paramater_sets[params_idx,]$n,
                                       sep = "")
  colnames(return_single_sim) <-c("Rejections est","Rejections real")
  return_single_sim
}) %>% do.call(what = rbind)
stopCluster(cluster)

saveRDS(sim_result,"sim_result_jumptest_eps_n.rds")

# 
# for (i in 1:nrow(Paramater_sets)) {
#   
# clusterExport(cluster, ls(), environment())  
# 
# Simulation = pbapply::pblapply(cl = cluster, X = 1:M, FUN = function(x){
#   
#   #Simulate path
#   set.seed(x)
#   Sim = Heston_Jumps_Sim(Heston_Jumps_Sim(T_val = T_val,
#                                           r = r,
#                                           alpha = alpha,
#                                           lambda = lambda,
#                                           sigma_v = sigma_v,
#                                           rho = rho,
#                                           S_0 = S_0,
#                                           V_0 = V_0,
#                                           n = n,
#                                           sigma_eps <- Paramater_sets[i,1],
#                                           intensity = intensity,
#                                           m = m,
#                                           sigma_J = Paramater_sets[i,2]))
#   
#   
#   
#   Path = Sim$Heston + Sim$Noise + Sim$Jumps
#   
#   #Apply estimators ----
#   
#   #TSRV:
#   TSRV_estimator_optimal(X = Path, 5, 1, 3, T_val)
#     
#   
#    
# })
# 
# }

