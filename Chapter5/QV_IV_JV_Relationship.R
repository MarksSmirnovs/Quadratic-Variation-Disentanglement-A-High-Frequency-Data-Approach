# Source ----------------------------------------------------------------
source("functions.R")
source("Jumps/jump_test_v2.R")
source("Microstructure_noise_estimators/Pre_average.R")
source("Jumps/BPV_preAvg.R")

#   -----------------------------------------------------------------------

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = paste0("Data/Synced_data_2007.db"))
synced_data <- DBI::dbGetQuery(con, "select * from SPY_2007_1sec") %>% collect()
DBI::dbDisconnect(con)

synced_data <- synced_data[seq(1, nrow(synced_data), 5),]

T_val <- 1/(synced_data$Date %>% as.Date() %>% unique %>% length())
n <- 4680
delta <- T_val/n
synced_data <- synced_data %>% mutate(Price = log(Price)) 

daily_return <- synced_data %>% group_by(Date = as.Date(Date)) %>%
  summarise(return = mean(diff(log(Price))))

# Return variation  -----------------------------------------------------------------------

QV_IV_JV <- readRDS("Application/jump_variation_result_2007.RDS")
dates <- QV_IV_JV %>% filter(JV_proportion < 0) %>% .$Date

df <- merge(daily_return, QV_IV_JV, by = "Date")

ggplot(df %>% filter(!(Date  %in%  dates))) + 
  geom_line(aes(x = Date, y = return, col = "Return"),size = 0.5) +
  geom_line(aes(x = Date, y = IV/300, col = "IV"), alpha = 0.4, size = 1) +
  geom_line(aes(x = Date, y = JV/300, col = "JV"), alpha = 0.4, size = 1)  +
  annotate("rect", xmin = as.Date("2007-07-22"), xmax =as.Date("2007-09-20"), ymin = -10, ymax =10,
           alpha = .2) + 
  annotate("rect", xmin = as.Date("2007-10-15"), xmax =as.Date("2007-12-20"), ymin = -10, ymax =10,
           alpha = .2) + 
  coord_cartesian(ylim=c(min(df$return), max(df$return)*1.15)) +  
  scale_color_manual(labels = c("Return", "IV", "JV"),
                     breaks = c("Return", "IV", "JV"),
                     values = c("Red",  "Blue", "Black"),
                     name = "") +
  # theme(axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank()) + 
  xlab("Date") + ylab("Daily Average Log Return")

# Correlation structure ---------------------------------------------------

VIX <- read_csv("Application/vixcurrent.csv", skip = 1)
VIX <- VIX %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>% 
  filter(Date >= "2007-01-01", Date < "2008-01-01") %>% select(Date, 'VIX Close') %>%
  dplyr::rename('VIX' = 'VIX Close')


predictors <- readRDS("Application/jump_variation_result_2007.RDS")
predictors <- predictors %>% select(-JV_proportion) %>%
  merge(VIX, by = "Date") %>% merge(daily_return, by = "Date") %>% 
  mutate(
    IV_1 = lag(IV, 1),
    JV_1 = lag(JV, 1),
    return_1 = lag(return, 1),
    VIX_1 = lag(VIX, 1),
    
    IV_5 = rollmean(x = IV_1, k = 5, align = "right", fill = NA),
    JV_5 = rollmean(x = JV_1, k = 5, align = "right", fill = NA),
    return_5 = rollmean(x = return_1, k = 5, align = "right", fill = NA),
    VIX_5 = rollmean(x = VIX_1, k = 5, align = "right", fill = NA),

    IV_22 = rollmean(x = IV_1, k = 22, align = "right", fill = NA),
    JV_22 = rollmean(x = JV_1, k = 22, align = "right", fill = NA),
    return_22 = rollmean(x = return_1, k = 22, align = "right", fill = NA),
    VIX_22 = rollmean(x = VIX_1, k = 22, align = "right", fill = NA),
    
    return_1 = ifelse(return_1 < 0, return_1, 0),
    return_5 = ifelse(return_5 < 0, return_5, 0),
    return_22 = ifelse(return_22 < 0, return_22, 0)
  ) %>% select(-IV, -JV, -return, -VIX) %>% na.omit

saveRDS(predictors,"Application/predictors.RDS")

# Pearson correlation  -----------------------------------------------------------------------

Predictors_cor = predictors[,c(1,2,3,7,11,4,8,12,6,10,14,5,9,13) ]

names(Predictors_cor) = gsub("return","r^-",names(Predictors_cor) )

names(Predictors_cor) = gsub("_","_{",names(Predictors_cor) )

cor_mat = matrix(0 , nrow = 2, ncol = 12)

for (i in 1:12) {
  
  cor_mat[1,i] = paste0("$",names(Predictors_cor)[i+2],"}$")
  
  cor_mat[2,i] = cor(Predictors_cor$QV, Predictors_cor[,i+2])
  
}

cor_mat[2,] %<>% as.numeric() %>%  round(2)

Latex_cor = kable(cor_mat, "latex", booktabs = T , escape = F) %>% kable_styling() %>% cat()

# cormat <- round(cor(predictors %>% select(-Date)),2)
# cormat[lower.tri(cormat)]<- NA
# cormat <- get_upper_tri(cormat)
# melted_cormat <- melt(cormat, na.rm = TRUE) %>% filter(Var1 == "QV")
# 
# 
# ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="Pearson Correlation") +
#   theme_minimal()+ 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1))+
#   coord_fixed() + 
#   geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.justification = c(1, 0),
#     legend.position = c(0.6, 0.9),
#     legend.direction = "horizontal")+
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
#                                title.position = "top", title.hjust = 0.5))


#  Scatter plots ----------------------------------------------------------


p1 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = IV_1))

p2 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = IV_5))

p3 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = IV_22))

p4 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = JV_1))

p5 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = JV_5))

p6 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = JV_22))

p7 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = VIX_1))

p8 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = VIX_5))

p9 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = VIX_22))

p10 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = return_1))

p11 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = return_5))

p12 <- ggplot(predictors) +
  geom_point(aes(y = QV, x = return_22))


gridExtra::grid.arrange(
  p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
  layout_matrix = rbind(c(1,2,3),c(4,5,6), c(7,8,9), c(10,11,12))
  #layout_matrix = rbind(c(1,2),c(3,4), c(5,6), c(7,8), c(9,10), c(11, 12))
)

#   -----------------------------------------------------------------------

ggplot(predictors) +
  geom_point(aes(y = QV, x = IV_22)) +
  geom_smooth(aes(y = QV, x = IV_22), se = FALSE, col = "red")

ggplot(predictors) +
  geom_point(aes(y = QV, x = VIX_22)) +
  geom_smooth(aes(y = QV, x = VIX_22), se = FALSE, col = "red")


