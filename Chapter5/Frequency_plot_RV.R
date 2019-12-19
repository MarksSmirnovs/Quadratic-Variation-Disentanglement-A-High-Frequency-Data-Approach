# Source ------------------------------------------------------------------

source("functions.R")


# Doc ---------------------------------------------------------------------

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = paste0("Data/Synced_data.db"))
synced_data <- DBI::dbGetQuery(con, "select * from SPY_2008_1sec")
DBI::dbDisconnect(con)


synced_data <- synced_data %>% mutate(
  Date = as.POSIXct(Date, origin = "1970-01-01")
)

synced_data %>% filter(Date >= "2008-01-01", Date < "2008-10-01") %>%
  mutate(price = log(price))  %>% group_by(as.Date(Date)) %>% summarise(
    unique_count = length(unique(price)),
     RV = sqrt(sum(diff(price)^2)/(2*unique_count))
  ) %>% summarise(
    mean(RV)
  )

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = paste0("Application/SPY_2008_Februar_us.db"))
feb_sync <- DBI::dbGetQuery(con, "select * from '1 sec'")
DBI::dbDisconnect(con)

feb_sync %>% group_by(as.Date(Date)) %>%  summarise(
  RV = sum(diff(log(Price))^2)/(2*n()) %>% sqrt
) 


freq <- lapply(c(1, seq(10, 60, 5), seq(120, 600, 60)), function(idx){

  RV <-  feb_sync[seq(1, nrow(feb_sync), by = idx),] %>% group_by(as.Date(Date)) %>%  summarise(
    RV = sqrt(sum(diff(log(Price))^2)/(2*n()))
  ) %>% summarise(
    RV = mean(RV)
  )
  colnames(RV) <-  paste0("sec_", idx)
  RV
}) %>% do.call(what = cbind)


freq <- t(freq) %>% data.frame %>%  mutate(sec = colnames(freq))
colnames(freq)[1] <- "RV"

plot(freq$RV)




