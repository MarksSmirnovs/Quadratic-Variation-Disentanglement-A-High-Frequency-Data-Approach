
# Setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pacman::p_load(tidyverse,
               lubridate,
               magrittr,
               dbplyr,
               RSQLite,
               DBI,
               zoo)

#Local path for SP500
SPY_dir = "C:/Users/Mathias/Desktop/High Frequency Stock Data/SPY"

Files = list.files(SPY_dir)[-1] 

#Select files from 2005

Years_to_remove = 1999:2004

for (i in seq_along(Years_to_remove)) {
  
  Index = grepl(Years_to_remove[i],Files)
  
  Files = Files[!Index]
  
}

rm(Index,i,Years_to_remove)


# Make functions for cleaning ---------------------------------------------

# Barndorff-Nielsen et al. (2009) Trade data cleaning
Clean_data = function(data,date){
  
  names(data)[1] = "datetime"
  
  Cleaned = data %>% 
    filter(hms(datetime) >= hms("09:30:00"), 
           hms(datetime) <= hms("16:00:00")  
  ) %>% 
    mutate(datetime = ymd_hms(paste(date, datetime, sep = " ")),
    cond = trimws(cond,"both")
  ) %>% 
    filter(price > 0,
           ex %in% c("N","T"),
           corr == 0,
           cond %in% c(NA,"","E","F")
  ) %>% 
    group_by(datetime) %>% 
    summarise(
      price = sum(price * volume)/sum(volume) 
  ) 
  
  # Remove outliers
  Index = c()

  n = nrow(Cleaned)
  
  for(i in 1:n){
    if(i <= 25){
      id = 1:(i+25)
      id = id[-i]
      med = median(Cleaned$price[id])
    }
    else if( (n-i) <= 25){
      id = (i-25):n
      id = id[-i]
      med = median(data$price[id])
    }
    else{
      id = (i-25):(i+25)
      id = id[-i]
      med = median(Cleaned$price[id])
    }
    
    mad = mean(abs(Cleaned$price[id]-med))
    
    if(Cleaned$price[i] > (med + 10*mad) | Cleaned$price[i] < (med - 10*mad)){
      Index = c(Index,i)
    }
    
  }
  
  if(!(length(Index)==0)){
    Cleaned = Cleaned[-Index,]
  }
  
  return(Cleaned)
}



# Write data to a data base for each year ---------------------------------

write_cleaned_to_DB = function(files){
  
  Years = c(2005:2009)
  
  lapply(Years, function(x){
  
  DB <- dbConnect(RSQLite::SQLite(), dbname = paste0("./SPY_",x,"_us.db"))
  
  files_year = files[grepl(x, substr(Files,5,12))]
  
  for (i in seq_along(files_year)) {
    
    #Load data
    Raw_data = read.csv( paste0(SPY_dir,"/", files_year[i]) )
    
    #Date
    Date = as.character(ymd(substr(files_year[i],5,12)))
    
    #Clean data
    Cleaned_Data = Clean_data(data = Raw_data, date = Date)
    
    #Insert as table in data base
    dbWriteTable(conn = DB, name = Date, value = Cleaned_Data, append = T, overwrite = F)
    
    #Keeping track of progress
    cat("File completed = ", files_year[i], "\n")
  }
  
  dbDisconnect(DB)
  
  })
  
}

write_cleaned_to_DB(Files)
