# Sourcing the Jump test and functions ---------------------------------------------------
source("functions.R")
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

intensity = 50/T_val
m = 0

#Make grid
sigma_eps = c(0.0005)
sigma_J = 0



#   -----------------------------------------------------------------------


hest_sim <- Heston_Sim(T_val = T_val, n = n , r = r, rho = rho , alpha = alpha, 
                       lambda = lambda, sigma_v = sigma_v, S_0 = S_0 , V_0 = V_0)

micro_noise <- sigma_eps * rnorm(n + 1)

Jumps_sim <- Compound_Pois(T_val = T_val , n = n + 1, intensity = intensity , m = m, sigma_J = sigma_J)

Path <- hest_sim$Heston + micro_noise + Jumps_sim

###
C <- 1/19
M <- floor(C*n^(1/2))

PreAvg_path <- zoo::rollmean(Path, M, align = "left")

PreAvg_path <- PreAvg_path[seq(1, length(PreAvg_path), M)]
delta <- T_val/length(PreAvg_path)

L_stat <- diff(PreAvg_path)


IV <- integrate(IV_fun, lower = 0, upper = T_val, vol_path = hest_sim$Vol, delta = delta, subdivisions = 10000000)$value
IV_MBPV <- MBPV(Y = Path, n = n, r = 1, l = 1, conf_alpha = 0.05, T_val = T_val)$MBPV

V_n <- 2/3*C^2*IV + 2*sigma_eps^2
X_stat <- (sqrt(M)/sqrt(V_n)) * L_stat


V_n_MBPV <- 2/3*C^2*IV_MBPV + 2*sqrt(sum(diff(Path)^2)/(2*n))^2
X_stat_MBPV <- (sqrt(M)/sqrt(V_n_MBPV)) * L_stat

empeos <- X_stat

theos <- rnorm(length(X_stat)) 

p1 <- QQ_plot_base(y = empeos, x = theos, xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles") + ggtitle(expression("a) True"* " IV and "*sigma[epsilon]))

empeos_MBPV <- X_stat_MBPV
p2 <- QQ_plot_base(y = empeos_MBPV, x = theos, xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles") + ggtitle(expression("b) Estimated"* " IV and "*sigma[epsilon]))

gridExtra::grid.arrange(p1, p2, layout_matrix = rbind(c(1,2)))

