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


# Jump test ---------------------------------------------------------------

# jump_tests <- synced_data %>%
#   group_by(Date = as.Date(Date)) %>% summarise(
#     Jump = Jump_test(C = 1/19, X = Price, alpha = 0.05, delta = delta)$Rejection
#   )
# 
# saveRDS(jump_tests, file = "Application/jump_tests_result_2007.RDS")

jump_tests <- readRDS("Application/jump_tests_result_2007.RDS")
jump_tests$Jump %>% mean


ggplot(jump_tests) + geom_point(aes(y = Jump, x = Date)) + ylab("Null Hypothesis Rejection") + xlab("Date")

# JV ----------------------------------------------------------------------

# C <- 1/19
# M <- floor(C*n^(1/2))
# 
# jump_variation <- synced_data   %>%
#   group_by(Date = as.Date(Date)) %>% mutate(preAvgPath = zoo::rollmean(Price, M, align = "left", fill = NA)) %>%
#   summarise(
#     QV = PreAvg_Estimator(X = Price, t = T_val, delta = delta, k_n = 51)$Estimate,
#     IV = BPV((preAvgPath %>% na.omit)[seq(1, n(), M)],
#              delta = T_val/length((preAvgPath %>% na.omit)[seq(1, n(), M)]))$BPV,
#     JV = QV - IV,
#     JV_proportion = JV/QV
#   )
# 
# saveRDS(jump_variation, file = "Application/jump_variation_result_2007.RDS")

jump_variation <- readRDS("Application/jump_variation_result_2007.RDS")

ggplot(jump_variation) + geom_line(aes(x = Date, y = JV_proportion))

dates <- jump_variation %>% filter(JV_proportion < 0) %>% .$Date

synced_data %>% filter(as.Date(Date)  %in% dates) %>%
  group_by(Date = as.Date(Date)) %>% summarise(
    length(unique(Price))
)
 
ggplot(jump_variation %>% filter(!(Date  %in%  dates))) + geom_line(aes(x = Date, y = JV_proportion)) + 
  xlab("Date") + ylab(expression(paste(frac(JV, QV)))) + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))


jump_variation %>% filter(!(Date  %in%  dates)) %>% summarise(mean(JV_proportion)) %>% round(3)


