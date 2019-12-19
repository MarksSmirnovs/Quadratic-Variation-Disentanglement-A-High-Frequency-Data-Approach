# Function ---------------------------------------------------------------

Heston_Sim = function(T_val, n , r , rho , alpha, lambda, sigma_v, S_0 , V_0){
 
  BM_Sim = function(T_val,n,d=1){
    
    BM_path = matrix(data = 0 , nrow = n , ncol = d)
    
    for (i in 2:n) {
      
      BM_path[i] = BM_path[i-1] + sqrt(T_val/n)*rnorm(d)
      
    }
    
    BM_path
    
  }
  
  BM_v = BM_Sim(T_val,n+1)
  BM = rho*BM_v + sqrt(1-rho^2)*BM_Sim(T_val,n+1)
  
  dBM_v = diff(BM_v)
  dBM = diff(BM)
  
  delta = T_val/n
  
  V = rep(V_0,n+1)
  S = rep(S_0,n+1)
  
  for (i in 1:n) {
    V[i+1] = V[i] + (alpha-lambda*V[i])*delta + sigma_v * sqrt(abs(V[i])) * dBM_v[i]
    
    S[i+1] = S[i] + r*delta + sqrt(abs(V[i+1])) * dBM[i] 
  }
  
  return(list(Heston = S,
              Vol = V))
  
}


# Example -----------------------------------------------------------------

# r = 0.05
# alpha = 0.5
# lambda = 0.5
# sigma_v = 0.4
# rho = -0.7
# S_0 = 1
# V_0 = 0.3
# n = 1000
# sigma_eps <-  0.01
# 
# Sim = Heston_Sim(T = T, n = n , r = r, rho = rho , alpha = alpha,
#            lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)
# 
# Sim$Stock %>% plot(type = "l")
# Sim$Vol %>% lines(col = "red")
