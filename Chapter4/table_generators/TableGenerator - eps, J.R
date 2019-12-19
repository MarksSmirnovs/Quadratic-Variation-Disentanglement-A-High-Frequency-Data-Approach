
# Setup -------------------------------------------------------------------

library(dplyr)
library(tibble)
library(kableExtra)
library(magrittr)
library(stringr)

Data = readRDS("SimulationStudy/results/sim_result_ej.rds")

TableGenerator = function(SimData){
  
  #Transpose data
  SimData = t(SimData) %>% as.data.frame() %>% round(3) %>% rownames_to_column("Parameters")
  
  #Remove MBPV Fixed
  SimData %<>% select(-MBPV_estimate_fix , - MBPV_estimate_fix_hit)  
  
  #Make df to insert into
  df = matrix(0, ncol = ncol(SimData)/2 + 1 , nrow = nrow(SimData)*2 ) %>% as.data.frame()
  
  names_seq = c(1, seq(2,12,2))
  
  names(df) = names(SimData)[names_seq]
  
  names(df)[5] = "Kernel"
  names(df)[6] = "MBPV"
  
  #Insert into df
  for(i in 1:(nrow(SimData))) {
    
      df[2*(i-1)+1,] = c(SimData[i,1], paste0(sprintf("%.3f",SimData[i,][names_seq][-1])) )
      
      df[2*i,] = c(SimData[i,1] , paste0("(", sprintf("%.3f",SimData[i,][seq(3,13,2)]),")" ) )
      
  }
  
  #Make best estimates and conf ints red
  for(i in 1:(nrow(SimData))) {
    #Estimates
    Min = min(as.numeric(df[2*(i-1)+1,-1]))
    
    Red = which(as.numeric(df[2*(i-1)+1,-1]) == Min)
    
    df[2*(i-1)+1,Red+1] = paste0( "\\textcolor{red}{", df[2*(i-1)+1,Red+1] ,"}")
    
    #Conf ints
    
    Red_int = abs(df[2*i,-1] %>% gsub("[(]","",.) %>% gsub("[)]","",.) %>% as.numeric()-0.95) %>% which.min()
    
    df[2*i,Red_int+1] = paste0( "\\textcolor{red}{",df[2*i,Red_int+1] ,"}")
    
  }

  sigma_eps =  c()
  jump_stds = c()
  
  #Format parameters
  for(i in 1:nrow(SimData)){
    
    Params = strsplit(df[2*(i-1)+1,1], "[|]")
    
    sigma_eps %<>% c(paste0(gsub("sigma_eps","$\\sigma_\\varepsilon", Params[[1]][1], fixed = TRUE),"$"))
    jump_stds %<>% c( paste0("\\rotatebox{90}{",gsub("jump_std","$\\sigma_J", Params[[1]][2], fixed = T),"$}"))
    
    df[2*(i-1)+1,1] = sigma_eps[i]
    
    df[2*i,1] = ""
    
  }
  
  n_sigma = sigma_eps %>% unique() %>% length()
  jump_stds %<>% unique()
  n_jump_stds = jump_stds %>% length()
 
  
  #Add row for groups:
  Groups = c()
  
  for (i in 1:n_jump_stds) {
    
    Groups %<>% c(rep(jump_stds[i], 2*n_sigma))
    
  }
  
  df = cbind(Groups,df)
  
  df = df[,c(1:6,8,7)]
  
  names(df)[c(1,2)] = ""
  
  #Make kable
  Latex = kable(df, "latex", booktabs = T , escape = F) %>% kable_styling()
  
  #Insert rows
  for (i in 1:(n_jump_stds-1)){
    Latex %<>% row_spec(i*2*n_sigma, hline_after = T)
  }
  
  Latex %<>%
    row_spec(0, bold = T) %>% 
    collapse_rows(columns = 1:2, latex_hline = "none", valign = "middle") %>% 
    column_spec(2,border_right = T) 

  
  return(Latex)
}

cat(TableGenerator(Data))
