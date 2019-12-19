# Setup -------------------------------------------------------------------

source("functions.R")
source("Jumps/MBPV.R")
library(dplyr)
library(Matrix)
library(evd)

# Test ----------------------------------------------

Jump_test=function(C,X,alpha,delta){
  #browser()
  n=length(X)
  M= ifelse(C*sqrt(n)>1,floor(C*sqrt(n)),1)
  
  Leftover = n %% M 
  
  if(Leftover > 0){X = X[-c((n-Leftover + 1):n)]}
  
  #Compute variance
  IV=MBPV(Y = X, n = n-Leftover, r = 1, l = 1, conf_alpha = 0.05,T_val = n*delta)$MBPV
  returns=diff(X)
  q_2_sqrd=sum(returns^2)/(n-1)
  V=2/3*IV*C^2+q_2_sqrd
  #V=sigma_eps
  
  f1 <- function(x, y){
    rep(x,y)
  }
  
  n_means = (n-Leftover)/M
  
  j=1:(n-M-Leftover)
  i=sapply(1:(n_means - 1),f1,y=M)%>%array()
  x <- rep(1,(n-M-Leftover))
  
  zeros=matrix(0L, nrow = max(i), ncol = M)
  
  A <- sparseMatrix(i=i, j=j,x=x)
  A_1=cbind(A,zeros)
  A_2=cbind(zeros,A)
  
  P_1=A_1%*%X*1/M
  P_2=A_2%*%X*1/M
  
  Chi=(abs(P_2-P_1)*sqrt(M)/sqrt(V))
  B=1/sqrt(2*log(n_means))
  A=sqrt(2*log(n_means))-(log(pi)+log(log(n_means)))/(2*(2*log(n_means))^(1/2))
  Xi=(max(Chi)-A)/B
  pvalue=pgumbel(Xi,lower.tail = F)
  
  rejections=ifelse(pvalue<alpha,1,0)
  
  output=list()
  output$Test=Xi
  output$Pvalue=pvalue
  output$Rejection=rejections
  return(output)
  
  
}
# Simulation test ----------------------------------------------



#source("Microstructure_noise_estimators/Heston_Sim.R")

# r = 0.05
# alpha = 0.04*5
# lambda = 5
# sigma_v = 0.5
# rho = -0.5
# S_0 = 1
# V_0 = 0.3
# n = 10000
# sigma_eps <-  0.0005
# T_val = 1/252
# delta = T_val/n
# 
# 
# Sim_noise = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
#                        lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)


# starttime1=Sys.time()
# Jump_test(1/19,Sim_noise$Stock,alpha=0.05)
# endtime1=Sys.time()
# 
# starttime2=Sys.time()
# Jump_test2(1/19,Sim_noise$Stock,alpha=0.05)
# endtime2=Sys.time()
# 
# endtime1-starttime1
# endtime2-starttime2

# tests=c()
# for(i in 1:1000){
#   
#   
#   Sim_noise = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
#                          lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)
#   
#   
#   
#   
#   tests[i]=Jump_test2(1/19,Sim_noise$Heston+Sim_noise$Noise,alpha=0.05,delta)$Rejection
#   
# }
# tests
# sum(tests%>%na.omit())

#  test of sparse matrix function----------------------------------------------

# 
# source("Microstructure_noise_estimators/Heston_Sim.R")
# 
# r = 0.05
# alpha = 0.04*5
# lambda = 5
# sigma_v = 0.5
# rho = -0.5
# S_0 = 1
# V_0 = 0.3
# n = 1e6
# sigma_eps <-  0.0005
# T_val = 1
# delta = T_val/n
# 
# 
# r = 0.05
# alpha = 0.04*5
# lambda = 5
# sigma_v = 0.5
# rho = -0.5
# S_0 = 1
# V_0 = 0.3
# n = 1000
# sigma_eps <-  0.0005
# T_val = 1
# delta = T_val/n
# 
# 
# Sim_noise = Heston_Sim(T = T_val, n = n , r = r, rho = rho , alpha = alpha,
#                        lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0, sigma_eps = sigma_eps)
# length(Sim_noise$Stock[3])
# M=floor(1/8*sqrt(n))
# 
# 
# 
# f1 <- function(x, y){
#   rep(x,y)
# }
# 
# f2 <- function(x, y){
#   x:(x+y-1)
# }
# 
# 
# 
# # rows=rep(1,M)
# # cols=1:M
# # for(i in 2:n){
# #   rows=append(rows,rep(i,M))
# #   cols=append(cols,i:(M+i-1))
# # }
# 
# floor(n/M)
# 
# j=1:(n-M+1)
# i=sapply(1:(floor(n/M)),f1,y=M)[1:length(j)]%>%array()
# x <- rep(1,(n-M+1))
# length(i)
# 
# zeros=matrix(0L, nrow = max(i), ncol = M)
# 
# A <- sparseMatrix(i=i, j=j,x=x)
# 
# 
# A_1=cbind(A,zeros)%>%as.matrix()
# 
# A_2=cbind(zeros,A)%>%as.matrix()
# 
# #%>%as.matrix()
# 
# test1=A_1%*%Sim_noise$Stock
# test2=A_2%*%Sim_noise$Stock
# 
# test3=test2-test1%>%as.numeric()


