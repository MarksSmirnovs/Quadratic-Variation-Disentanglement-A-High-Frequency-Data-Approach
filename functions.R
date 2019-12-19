# Packload ----------------------------------------------------------------
library(dplyr)
library(lubridate)
library(zoo)
library(parallel)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(dbplyr)
library(RSQLite)
library(DBI)
library(Matrix)
library(evd)
library(readr)
library(kableExtra)

# Utility functions -------------------------------------------------------

#Reads & creates a datetime column, the path must be the folder in which the data is contained.
#It reads every csv file, so seperate different time series into different folders
read_csv_folder <- function(path){
  
  file_names <- list.files(path = path, pattern="*.csv")
  
  data_return <- lapply(file_names, function(x) read.csv(paste(path, x, sep = ""), stringsAsFactors = FALSE) %>%
                          mutate(date = sub(".*_ *(.*?) *.csv.*", "\\1",  x)) %>%
                          mutate(date = as.POSIXct(paste(date, utcsec), format = "%Y%m%d %H:%M:%S")) %>% 
                          select(-utcsec))
  
  data_return <- do.call(rbind, data_return)
  return(data_return)
}


# High resolution PNG  ----------------------------------------------------
HR_Png <- function(name, figure){
  png(paste(name,".png", sep =""), units = "px", width = 1600, height = 800, res = 300)
  print(figure)
  dev.off()
}


# Q-Q plot ----------------------------------------------------------------

QQ_plot_base <- function(x, y, xlab, ylab){
  order_x <- x[order(x)]
  order_y <- y[order(y)]
  ggplot() + geom_point(aes(x = order_x, y = order_y)) +
    geom_abline(intercept = 0, slope = 1, color = "red", size = 1) + 
    xlab(paste(xlab)) + ylab(paste(ylab))
}



# Data sync function ------------------------------------------------------

sync_data <- function(data, tick_time, cores = 2){
  
  all_dates <- as.Date(data$datetime) %>% unique()
  all_hours <- hour(data$datetime) %>% unique()
  
  tick_grid <- data.frame("tick_grid" = seq(as.POSIXct(paste(all_dates[1], paste("0",min(all_hours),":00:00", sep = "")), tz = "UTC"),
                                            as.POSIXct(paste(all_dates[length(all_dates)],
                                                             paste(max(all_hours)+1,":00:00", sep = "")), tz = "UTC"), paste(tick_time)))
  tick_grid <- subset(tick_grid, format(tick_grid, format = "%H:%M:%S") >= "09:30:00" &
                        format(tick_grid, format = "%H:%M:%S") <= "16:00:00") %>% filter(!chron::is.weekend(tick_grid))
  
  first <-  data.frame(Date = tick_grid[1,1], price = data$price[1])
  
  
  
  cluster <- makeCluster(cores)
  clusterEvalQ(cluster, c(library(dplyr), library(lubridate)))
  clusterExport(cluster, ls(), environment())
  
  synced_data <- pbapply::pblapply(X = 2:nrow(tick_grid), cl = cluster , FUN = function(i){
    prices_bw_ticks <- subset(data, datetime >= tick_grid[i-1,1] & datetime <= tick_grid[i,1] )
    if(nrow(prices_bw_ticks) == 0){
      syn_tick <- data.frame(Date =  tick_grid[i,1], price = NA)
    }else{
      prices_bw_ticks <- prices_bw_ticks[order(prices_bw_ticks$datetime, decreasing = TRUE),]
      syn_tick <- data.frame(Date =  tick_grid[i,1], price = prices_bw_ticks$price[1])
    }
    syn_tick 
    
  }) %>% data.table::rbindlist()
  stopCluster(cluster)
  
  synced_data <- rbind(first, synced_data) %>% mutate(price = zoo::na.locf(price))
  synced_data 
  return(synced_data)
}

# Get from DB -------------------------------------------------------------

#Gets data from SPY data bases for all dates from start_date to end_date.
#sync chooses which data base to get the data from. Default is unsynchronized.

get_DB_data = function(start_date, end_date, sync = "us" ){
  
  Dates = seq.Date(from = as.Date(start_date), to = as.Date(end_date) , by = 1)
  
  Years = year(Dates) %>% unique
  
  Collected_all = lapply(Years, function(x){
    
    Dates_in_year = Dates[year(Dates) == x]
    
    con = dbConnect(RSQLite::SQLite(), dbname = paste0("./Data/SPY_",x,"_",sync,".db") )
    
    Dates_in_DB = dbListTables(con) %>% as.Date()
    
    Dates_to_collect = Dates_in_DB[which(Dates_in_DB %in% Dates_in_year)]
    
    #Collect data from DB and make large df
    
    Collected_year = lapply(Dates_to_collect, function(y){
      
      Collected = tbl(con, as.character(y)) %>% 
        collect() %>% 
        mutate(datetime = as.POSIXct(datetime, origin="1970-01-01" , tz = "GMT"))
      
    })
    
    dbDisconnect(con)
    
    return(do.call(rbind,Collected_year))
    
  })
  
  return(do.call(rbind,Collected_all))
  
}


# TSRV Estimator ----------------------------------------------------------

#IV average estimatoren

IV_avg=function(X,K){
  n=length(X)
  
  G=list()
  for(i in 1:K){
    G[[i]]=seq(i,n,K)
  }
  
  
  IV_avg=c()
  for(i in 1:K){
    X_i=X[G[[i]]]
    IV_avg[i]=sum(diff(X_i)^2)
  }
  
  
  return(mean(IV_avg))
}


#Example:
#get_DB_data(start_date = "2008-12-25", end_date = "2009-01-10")
