#   -----------------------------------------------------------------------
source("functions.R")
source("Jumps/BPV_adj.R")


#   -----------------------------------------------------------------------
BPV_preAvg <- function(X, sigma_eps, n, BstrapM, BstrapSplit, alpha,T_val){
  
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
              "0.0005" = {1/19},
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