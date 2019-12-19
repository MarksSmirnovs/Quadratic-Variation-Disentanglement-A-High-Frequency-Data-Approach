
# Setup -------------------------------------------------------------------

library(dplyr)
library(tibble)
library(kableExtra)
library(magrittr)
library(stringr)

SimData <- readRDS("SimulationStudy/results/BPV_preAvg_sim1.rds")
SimData_2 <- readRDS("SimulationStudy/results/BPV_preAvg_sim2.rds")

SimData <- rbind(SimData, SimData_2)


Data <- SimData %>% round(3)

rw_names <- rownames(Data)
rnames <- sub(".*\\|", "", rw_names) %>% unique
cnames <- sub("\\|.*", "", rw_names) %>% unique
latex_mat <- (matrix(NA, nrow = length(rnames)*2,
                     ncol = cnames  %>% length))

for(idx in seq(1, length(rnames)*2, 2)){
  rnames  %<>% append(values =  " ", after = idx)
}

colnames(latex_mat) <- cnames

row.names(latex_mat) <- rnames


for(idx in 1:nrow(Data)){
  col_idx <- which(colnames(latex_mat) == paste(sub("\\|.*", "", rw_names[idx])))
  r_idx <- which(rownames(latex_mat) == paste(sub(".*\\|", "", rw_names[idx])))
  latex_mat[r_idx, col_idx] <- Data[idx, 2]
  latex_mat[r_idx + 1, col_idx] <- paste("(",sprintf("%.3f", Data[idx, 3]), ")", sep = "")
  
}


colnames(latex_mat) <- paste0("$\\sigma_{\\varepsilon}  =", colnames(latex_mat), "$")

rownames(latex_mat) <- ifelse(rownames(latex_mat) != " ", paste0("$\\sigma_J=", rownames(latex_mat), "$"), "")


kable(latex_mat, "latex", booktabs = T , escape = F) %>% kable_styling() %>%  column_spec(1,border_right = T) 


