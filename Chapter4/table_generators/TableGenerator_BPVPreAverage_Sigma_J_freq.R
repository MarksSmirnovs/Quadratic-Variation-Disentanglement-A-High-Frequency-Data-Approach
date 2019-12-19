
# Setup -------------------------------------------------------------------

library(dplyr)
library(tibble)
library(kableExtra)
library(magrittr)
library(stringr)


SimData <- readRDS("SimulationStudy/results/BPV_preAvg_Sigma_J_n_freq_total.rds")


#Transpose data
SimData = SimData %>% as.data.frame() %>% round(3) %>% rownames_to_column("Parameters")

SimData <- SimData %>% arrange(as.numeric(sub("\\|.*", "\\1", Parameters)))

df = SimData

sigma_eps <- c()
jump_stds <- c()



#Format parameters
for(i in 1:nrow(SimData)){
  
  Params = strsplit(df[i,1], "[|]")
  
  sigma_eps %<>% c(paste0("$n = ", sub("*\\|.*", "\\1",  Params[[1]][2]), "$"))
  jump_stds %<>% c(paste0("\\rotatebox{90}{\\tiny $\\sigma_J = ", sub(".*\\|", "\\1", Params[[1]][1]), "$}"))
  
  
  df[i,1] = sigma_eps[i]
  
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
colnames(df) <- gsub(pattern = "BPV_adj_Hits", replacement = "$Coverage$", x = colnames(df))
colnames(df) <- gsub(pattern = "BPV_adj", replacement = "$\\IV\\^{BPV,adj}$", x = colnames(df))
colnames(df) <- gsub(pattern = "Best Boot", replacement = "$B$", x = colnames(df))



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



#   -----------------------------------------------------------------------

rm(list = ls())



# Setup -------------------------------------------------------------------

library(dplyr)
library(tibble)
library(kableExtra)
library(magrittr)
library(stringr)


SimData <- readRDS("SimulationStudy/results/BPV_preAvg_Sigma_eps_n_freq_total.rds")


#Transpose data
SimData = SimData %>% as.data.frame() %>% round(3) %>% rownames_to_column("Parameters")

SimData <- SimData %>% arrange(as.numeric(sub("\\|.*", "\\1", Parameters)))

df = SimData

sigma_eps <- c()
jump_stds <- c()



#Format parameters
for(i in 1:nrow(SimData)){
  
  Params = strsplit(df[i,1], "[|]")
  
  sigma_eps %<>% c(paste0("$n = ", sub("*\\|.*", "\\1",  Params[[1]][2]), "$"))
  jump_stds %<>% c(paste0("\\rotatebox{90}{\\tiny $\\sigma_\\varepsilon = ", sub(".*\\|", "\\1", Params[[1]][1]), "$}"))
  
  
  df[i,1] = sigma_eps[i]
  
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
colnames(df) <- gsub(pattern = "BPV_adj_Hits", replacement = "$Coverage$", x = colnames(df))
colnames(df) <- gsub(pattern = "BPV_adj", replacement = "$\\IV\\^{BPV,adj}$", x = colnames(df))
colnames(df) <- gsub(pattern = "Best Boot", replacement = "$B$", x = colnames(df))



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

