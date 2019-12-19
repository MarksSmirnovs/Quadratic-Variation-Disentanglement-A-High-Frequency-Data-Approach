setwd("~/P9/Code_Library")
source("functions.R")
source("Microstructure_noise_estimators/Heston_Sim.R")

IV_avg_estimator=function(X,K){ #Put it inside into TSRV estimator, so we only need to call one function
  n=length(X)
  
  G=list()
  for(i in 1:K){
    G[[i]]=seq(i,n,K) #Subgrid overlaps fx. K = 4 n = 10, i = 1 & i = 5
  }

  
  IV_avg=c()
  for(i in 1:K){
    X_i=X[G[[i]]]
    IV_avg[i]=sum(diff(X_i)^2)
  }


  return(mean(IV_avg))
}


TSRV_estimator=function(X,K){
  IV_avg=IV_avg_estimator(X,K)
  
  n_bar=n/K%>%floor()
  
  sigma_eps_hat=(n_bar/n)*sum(diff(X)^2)
  
  TSRV=IV_avg-sigma_eps_hat
  return(TSRV)
}


TSRV_c_estimator=function(X,M,c_1,I,T){ #Put it into TSRV estimator for one function call (also avoids calculating sigma eps multiple times)
  n=length(X)
  K_1=floor(c_1*n^{2/3})

  K_2=c_1*I*n^{2/3}%>%floor()

  TSRV_1=numeric()
  TSRV_2=numeric()
  
  for(i in 1:M){

    T_m=(i/M)*n # Multiplikation has higher priority than division
    T_m_previous=((i-1)/M)*n+1 # Multiplikation has higher priority than division
    
    X_m=X[T_m_previous:T_m]
    
    TSRV_1[i]=TSRV_estimator(X_m,K_1) #This is the  estimator for K groups (original grid)
    TSRV_2[i]=TSRV_estimator(X_m,K_2) #This estimator should be for K * I groups 
    
  }
  
  s_hat=sum((diff(TSRV_1)-diff(TSRV_2))^2)*n^(1/3) #Should it be n^(1/3) * sum((diff(TSRV_1) - diff(TSRV_2))^2)
  sigma_eps_hat=sum(diff(X)^2)/(2*n)
  
  eta_hat=(s_hat-8*sigma_eps_hat^2*c_1^{-2}*(1+I^{-2}-I^{-1}))/(c_1*(I^{1/2}-1)^2*T) #Fix T to something else
  
  c=((16*sigma_eps_hat^2)/(T*eta_hat))^{1/3} #fix T
  
  return(c)
}


# Example/test ------------------------------------------------------------


# Zhangs parametre.
r = 0.05
lambda = 5
alpha = 0.4*lambda
sigma_v = 0.5
rho = -0.5
S_0 = 1
V_0 = 0.3
n = 10000000
sigma_eps <-  0.0005
T_val = 1/252
sigma_eps=0.0005






set.seed(123)
sim = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
                 lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0,sigma_eps = sigma_eps)

X=sim[[1]]
Vol=sim[[2]]


c=TSRV_c_estimator(X,5,1,3,T_val)

K=c*n^{2/3}%>%floor()

IV_est=TSRV_estimator(X,K)

IV_fun <- function(vol_path, t, delta){
  eval_vals <- floor(t / delta)
  eval_vals <- ifelse(eval_vals == 0, 1, eval_vals)
  vol_path[eval_vals]
}



IV=integrate(IV_fun, lower = 0, upper = T_val, vol_path = Vol, delta = T_val/n, subdivisions = 10000000)

IV_est-IV$value


