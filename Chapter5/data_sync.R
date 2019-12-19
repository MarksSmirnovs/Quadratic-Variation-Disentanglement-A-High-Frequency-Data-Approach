# Source ----------------------------------------------------------------
source("functions.R")

# Params ------------------------------------------------------------------
start_date <- "2008-01-01"
end_date <- "2008-12-31"

# Doc ---------------------------------------------------------------------

Data <-  get_DB_data(start_date = start_date, end_date = end_date)

synced_data <- sync_data(Data, tick_time = "1 sec", cores = 3)

DB <- dbConnect(RSQLite::SQLite(), dbname = paste0("Data/Synced_data.db"))
dbWriteTable(conn = DB, name = "SPY_2008_1sec", value = synced_data, append = T, overwrite = F)
dbDisconnect(DB)
