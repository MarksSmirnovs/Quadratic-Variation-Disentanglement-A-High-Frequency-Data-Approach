library(dplyr)
library(tibble)
library(kableExtra)
library(magrittr)
library(stringr)

Data <- readRDS("SimulationStudy/results/sim_study_RV_1min.rds")

Data  %<>% round(3) %>% t()

rw_names <- rownames(Data)
rnames <- sub(".*sigma_J =  *(.*?)", "\\1", rw_names) %>% unique
cnames <- sub(".*sigma_eps = *(.*?) *\\|.*", "\\1", rw_names) %>% unique
latex_mat <- (matrix(NA, nrow = length(rnames)*2,
                     ncol = cnames  %>% length))

colnames(latex_mat) <- cnames

for(idx in seq(1, length(rnames)*2, 2)){
  rnames  %<>% append(values =  " ", after = idx)
}

row.names(latex_mat) <- rnames



for(idx in 1:nrow(Data)){
  col_idx <- which(colnames(latex_mat) == paste(sub(".*sigma_eps = *(.*?) *\\|.*", "\\1", rw_names[idx])))
  r_idx <- which(rownames(latex_mat) == paste(sub(".*sigma_J =  *(.*?)", "\\1", rw_names[idx])))
  latex_mat[r_idx, col_idx] <- Data[idx, 1]
  latex_mat[r_idx + 1, col_idx] <- paste("(",sprintf("%.3f", Data[idx, 2]), ")", sep = "")
 
}


colnames(latex_mat) <- paste0("$\\sigma_{\\varepsilon}  =", colnames(latex_mat), "$")

rownames(latex_mat) <- ifelse(rownames(latex_mat) != " ", paste0("$\\sigma_J=", rownames(latex_mat), "$"), "")
# 
# for(idx in 1:nrow(latex_mat)){
#   if(rownames(latex_mat)[idx] != ""){
#     latex_mat[idx, which.min(latex_mat[idx,])] <- paste("\\textcolor{red}{",
#                                                         latex_mat[idx, which.min(latex_mat[idx,])],
#                                                         "}", sep ="")
#   
#   }else{
#     latex_mat[idx, which.min(latex_mat[idx,])] <- paste("\\textcolor{red}{(",
#                                                                latex_mat[idx, which.min(abs(0.95 - as.numeric(latex_mat[idx,])))],
#                                                           ")}", sep ="")  
# 
#   }
# }



# 
# as.numeric(latex_mat[idx,])
# 
# sub("(", "", latex_mat[idx,])
# 
# latex_mat[idx] %>% sub("(", "")

kable(latex_mat, "latex", booktabs = T , escape = F) %>% kable_styling() %>%  column_spec(1,border_right = T) 

 