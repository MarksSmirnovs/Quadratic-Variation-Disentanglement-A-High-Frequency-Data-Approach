
# Setup -------------------------------------------------------------------

library(dplyr)
library(tibble)
library(kableExtra)
library(magrittr)
library(stringr)

Data = readRDS("SimulationStudy/sim_result_ej.rds")
RV_df <- readRDS("SimulationStudy/results/sim_study_RV_5min.rds")

TableGenerator = function(SimData){
  
  #Transpose data
  SimData = t(SimData) %>% as.data.frame() %>% round(3) %>% rownames_to_column("Parameters")
  
  #Remove MBPV Fixed
  SimData %<>% select(-MBPV_estimate_fix , - MBPV_estimate_fix_hit)  
  
  SimData <- SimData[,!grepl("hit", colnames(SimData))] %>% data.frame
  SimData <- split(SimData, seq_len(nrow(SimData)))
  SimData <- lapply(SimData, function(idx) data.frame(Parameters = idx$Parameters, RAD =  min(idx[-1]))) %>% do.call(what = rbind) %>% 
   mutate(Parameters = as.character(Parameters),
          Parameters = gsub(pattern = "jump_std", replacement = "sigma_J", x = Parameters)) 
    

    
  RV_df <- t(RV_df) %>% data.frame(check.rows = FALSE)
  RV_df <- RV_df %>%  mutate(Parameters = rownames(RV_df)) %>% select(-RV_hit)


  SimData <- merge(RV_df, SimData, by = "Parameters", sort = FALSE) %>% mutate(RAD_perc = round((RV - RAD)/RAD,3) ) %>%
  select(Parameters, RAD_perc)
  
  rownames(SimData) <- SimData$Parameters 
  SimData <- SimData %>% select(RAD_perc)
  
  rw_names <- rownames(SimData)
  rnames <- sub(".*sigma_J =  *(.*?)", "\\1", rw_names) %>% unique
  cnames <- sub(".*sigma_eps = *(.*?) *\\|.*", "\\1", rw_names) %>% unique
  latex_mat <- (matrix(NA, nrow = length(rnames),
                       ncol = cnames  %>% length))
  
  colnames(latex_mat) <- cnames
  
 
  row.names(latex_mat) <- rnames
  
  
  
  for(idx in 1:nrow(SimData)){
    col_idx <- which(colnames(latex_mat) == paste(sub(".*sigma_eps = *(.*?) *\\|.*", "\\1", rw_names[idx])))
    r_idx <- which(rownames(latex_mat) == paste(sub(".*sigma_J =  *(.*?)", "\\1", rw_names[idx])))
    latex_mat[r_idx, col_idx] <- SimData[idx, 1]
  }
  
  
  colnames(latex_mat) <- paste0("$\\sigma_{\\varepsilon}  =", colnames(latex_mat), "$")
  
  rownames(latex_mat) <- ifelse(rownames(latex_mat) != " ", paste0("$\\sigma_J=", rownames(latex_mat), "$"), "")
  
  Latex <- kable(latex_mat, "latex", booktabs = T , escape = F) %>% kable_styling() %>%  column_spec(1,border_right = T) 
  
  
  return(Latex)
}

cat(TableGenerator(Data))
