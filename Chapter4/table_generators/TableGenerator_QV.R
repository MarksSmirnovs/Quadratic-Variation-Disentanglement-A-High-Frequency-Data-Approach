library(dplyr)
library(tibble)
library(kableExtra)
library(magrittr)
library(stringr)

SimData <- readRDS("SimulationStudy/results/QV_sim.rds")

  #Transpose data
  SimData = t(SimData) %>% as.data.frame() %>% round(3) %>% rownames_to_column("Parameters")
  
  df = SimData
  
  sigma_eps <- c()
  jump_stds <- c()
  
  #Format parameters
  for(i in 1:nrow(SimData)){
    
    Params = strsplit(df[i,1], "[|]")
    
    sigma_eps %<>% c(paste0(gsub("sigma_eps","$\\sigma_\\varepsilon", Params[[1]][1] , fixed = TRUE),"$"))
    jump_stds %<>% c( paste0("\\rotatebox{90}{",gsub("sigma_J","\\tiny $\\sigma_J", Params[[1]][2], fixed = T),"$}"))
    
    df[i,1] = sigma_eps[i]
    
  }
  
  for(i in 1:(nrow(SimData))) {
    #Estimates
    Min = min(as.numeric(df[i,-1]))
    
    Red = which(as.numeric(df[i,-1]) == Min)
    
    df[i,Red+1] = paste0( "\\textcolor{red}{", df[i,Red+1] ,"}")
    
  }
  
  n_sigma = sigma_eps %>% unique() %>% length()
  jump_stds %<>% unique()
  n_jump_stds = jump_stds %>% length()
  
  
  #Add row for groups:
  Groups = c()
  
  for (i in 1:n_jump_stds) {
    
    Groups %<>% c(rep(jump_stds[i], n_sigma))
    
  }
  
  df = cbind(Groups,df)
  
  colnames(df)[1:2] <- ""
  colnames(df) <- gsub(pattern = "PreAvg", replacement = "$\\IV\\^{PRE}$", x = colnames(df))
  colnames(df) <- gsub(pattern = "Kernel", replacement = "$\\IV\\^{Kernel}$", x = colnames(df))
  colnames(df) <- gsub(pattern = "TSRV", replacement = "$\\IV\\^{TSRV}$", x = colnames(df))
  colnames(df) <- gsub(pattern = "Pre_Average_Sample_MBPV", replacement = "PreAvg MBPV", x = colnames(df))
  colnames(df) <- gsub(pattern = "Pre_Average_Sample_BPV", replacement = "PreAvg BPV", x = colnames(df))

  
  #Make kable
  Latex = kable(df, "latex", booktabs = T , escape = F) %>% kable_styling()
  
  #Insert rows
  for (i in 1:(n_jump_stds-1)){
    Latex %<>% row_spec(i*n_sigma, hline_after = T)
  }
  
  Latex %<>%
    row_spec(0, bold = T) %>% 
    collapse_rows(columns = 1:2, latex_hline = "none", valign = "middle") %>% 
    column_spec(2,border_right = T) 
  
  
cat(Latex)




