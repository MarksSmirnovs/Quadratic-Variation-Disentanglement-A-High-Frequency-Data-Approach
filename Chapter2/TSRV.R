#setwd("~/P9/Code_Library")
source("functions.R")
source("Microstructure_noise_estimators/Heston_Sim.R")


#X-- The log prices
#M,c_1,I--- Used for finding the optimal number of groups
#I should be 3,4. c_1 and M can be arbitrary as long as K_1,K_2<1/M*n


TSRV_estimator_optimal=function(X,M,c_1,I,T_val , alpha = 0.05){

  #---Calculate the IV^avg estimator
  IV_avg_estimator=function(X,K){ 
    n=length(X)
    
    #Create the groups
    G=list()
    for(i in 1:K){
      G[[i]]=seq(i,n,K) 
    }
    
    #Calculate the RV for the various groups
    IV_avg=c()
    for(i in 1:K){
      X_i=X[G[[i]]]
      IV_avg[i]=sum(diff(X_i)^2)
    }
    
    
    return(mean(IV_avg))
  }
  
  
  #---calculate the TSRV estimator based on a particular number of groups K
  TSRV_estimator=function(X,K){
    
    n=length(X)
    
    IV_avg=IV_avg_estimator(X,K)
    
    #Calculate the number of observations in each group. 
    n_bar=floor(n/K)
    
    #Estimate the variance of the noise term
    sigma_eps_hat=(n_bar/n)*sum(diff(X)^2)
    
    #substract the variamce of the noise term to get an unbiased estimate. 
    TSRV=IV_avg-sigma_eps_hat
    
    return(TSRV)
    
  }
  
  
  #---Calculate the optimal c
  TSRV_c_estimator=function(X=X,M=M,c_1=c_1,I=I,T_val=T_val,spread=FALSE){ 
    n=length(X)
    
    #Find the number of groups for the two groupings.
    K_1=floor(c_1*n^{2/3})
    
    K_2=floor(c_1*I*n^{2/3})
    
    TSRV_1=numeric()
    TSRV_2=numeric()
    
    #Calculate the two TSRV estimators needed for the estimate of s
    for(i in 1:M){
      
      T_m=(i/M)*n 
      T_m_previous=((i-1)/M)*n+1 
      
      X_m=X[1:T_m]
      
      TSRV_1[i]=TSRV_estimator(X_m,K_1) 
      TSRV_2[i]=TSRV_estimator(X_m,K_2) 
      
    }
    
    #Find the estimate of s
    s_0_hat=sum((diff(TSRV_1)-diff(TSRV_2))^2)*n^(1/3) 
    
    #Find estimate of sigms_eps
    sigma_eps_hat=sum(diff(X)^2)/(2*n)
    
    
    #Calculate the spread if needed 
    
    if(spread==TRUE){
      
      s_hat=(8*(c_1^{-2}-c_1^{-2}*(I^{-2}-I^{-1}+1)/(I^{1/2}-1)^2)*sigma_eps_hat^2+s_0_hat/(I^{1/2}-1)^2)%>%sqrt()
      return(s_hat)
    }
    
    #Find the estimate of eta
    eta_hat=(s_0_hat-8*sigma_eps_hat^2*c_1^{-2}*(1+I^{-2}-I^{-1}))/(c_1*(I^{1/2}-1)^2*T_val) 
    
    #Find the optimal c
    c=((16*sigma_eps_hat^2)/(T_val*eta_hat))^{1/3}
    c = ifelse(c < 1, 1, c)
    
    return(c)
  }
  
  #---Find the optimal number of groups and calculate the TSRV estimator
 
  c=TSRV_c_estimator(X,M,c_1,I,T_val)
  
  K=ceiling(c*n^{2/3})
  
  TSRV=TSRV_estimator(X,K)
  
  TSRV_sigma=TSRV_c_estimator(X,M,c,I,T_val,spread = TRUE)

  output=list()
  output$TSRV=TSRV
  output$TSRV_sigma=TSRV_sigma
  output$upper=round(TSRV + qnorm(1-alpha/2)*TSRV_sigma/(length(X)^{1/6}),10)
  output$lower=round(TSRV + qnorm(alpha)*TSRV_sigma/(length(X)^{1/6}),10)
  return(output)
  
}
# Example/test ------------------------------------------------------------


# ## Zhangs parametre.
# r = 0.05
# lambda = 5
# alpha = 0.4*lambda
# sigma_v = 0.5
# rho = -0.5
# S_0 = 1
# V_0 = 0.3
# n =6.5*60*60
# sigma_eps <-  0.0005
# T_val = 1/252
# sigma_eps= 0.0005
# delta = T_val/n
# 
# # test=numeric()
# 
# # for(i in 1:1000){
# #set.seed(123)
# # sim = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
#                  # lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0,sigma_eps = sigma_eps)
# # 
# # X=sim[[1]]
# # Vol=sim[[2]]
# # 
# # IV_est=TSRV_estimator_optimal(sim$Stock,5,1,3,T_val)
# # 
# # 
# # 
# # IV_fun <- function(vol_path, t, delta){
# #   eval_vals <- floor(t / delta)
# #   eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
# #   vol_path[eval_vals]
# # }
# # 
# # IV=integrate(IV_fun, lower = 0, upper = T_val, vol_path = Vol, delta = T_val/n, subdivisions = 10000000)
# # 
# # test[i]=ifelse(IV$value>IV_est$lower&IV$value<IV_est$upper,1,0)
# # 
# # }
# 
# 
# IV_fun <- function(vol_path, t, delta){
#   eval_vals <- floor(t / delta)
#   eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
#   vol_path[eval_vals]
# }
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
#     z2 <- TSRV_estimator_optimal(z1$Stock,5,1,3,T_val)
#     
#     if(class(z2) == "try-error"){
#       z2 <- NULL
#       z2$value = NA
#       z2$TSRV = NA
#       z2$upper = NA
#       z2$lower = NA
#     }
#     
#     z3 <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = z1$Vol, delta = delta, subdivisions = 10000000)
#     
#     
#     
#     data.frame(
#       IV = z3$value,
#       TSRV = z2$TSRV,
#       upper = z2$upper,
#       lower = z2$lower,
#       conf = z2$upper >= z3$value & z3$value  >= z2$lower
#     )
#     
#   })
#   
#   stopCluster(cluster)
#   
#   z1 <- z1 %>% do.call(what = rbind)
#   
#   z1$conf %>%  mean
# }

