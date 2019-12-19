# Sources -----------------------------------------------------------------
source("Code_Library/functions.R")
load("Code_Library/sp500_2003.Rdata")

# Barndorff-Nielsen et al. (2009) Trade data cleaning -------------------------------------------------------------------
                                                                      
cond_names <- sp500_2003$cond %>% unique()

SP500_2003_cleaned <- sp500_2003 %>% filter(
                      ex %in% c("N","T"), #P3
                      corr == 0, #T1
                      cond %in% c(NA, cond_names[3], cond_names[5], cond_names[10]) #T2
                      ) %>% group_by(date) %>% summarise(    #T3
                        price = sum(price * volume)/sum(volume)  #weighted average price
                      ) %>% mutate(roll_median=rollmedian(price, 51,  align = "center",na.pad = TRUE
                      ))%>%mutate(MAD=zoo::rollmean(abs(price-roll_median), 51,align = 
                                                      "center",na.pad = TRUE)) %>%
                      filter(price <= 10*MAD+roll_median, 
                             price >= -10*MAD+roll_median)

# Bauwens et al. (2012) Synchronization -----------------------------------

SP500_2003_cleaned_n_synced <- sync_data(data = SP500_2003_cleaned, tick_time = "1 min")

trading_days <- as.Date(sp500_2003$date, tz = "CET") %>% unique()
SP500_2003_cleaned_n_synced <- SP500_2003_cleaned_n_synced %>% filter(as.Date(Date, tz = "CET") %in% trading_days)

save(SP500_2003_cleaned_n_synced, file = "Code_Library/sp500_2003_cleaned_n_sync.Rdata", compress = "xz")


# Comparison plot ---------------------------------------------------------
rm(list = ls())
load("Code_Library/sp500_2003.Rdata")
load("Code_Library/sp500_2003_cleaned_n_sync.Rdata")
source("Code_Library/libs.R")

min_price_clean <- min(SP500_2003_cleaned_n_synced$Price)
max_price_clean <- max(SP500_2003_cleaned_n_synced$Price)

original <- ggplot(data = sp500_2003, mapping = aes(x = date, y = price)) +
  geom_line() +
  coord_cartesian(xlim = NULL, ylim = c(min_price_clean, max_price_clean), expand = TRUE,
                  default = FALSE, clip = "on") +
  xlab("Date") + ylab("Price")

cleaned <- ggplot(data = SP500_2003_cleaned_n_synced, mapping = aes(x = Date, y = Price)) +
  geom_line()

png("sp500_2003_cleaned.png", units = "px", width = 1600, height = 800, res = 300)
cleaned
dev.off()
png("sp500_2003_origin.png", units = "px", width = 1600, height = 800, res = 300)
original
dev.off()


# 2008 --------------------------------------------------------------------
# Fetching and cleaning data in chunks due to memory limit.

# load("sp500_2008_1.Rdata")
# SP500_2008_1_cl <- sp500_2008_1 %>% filter(
#   ex %in% c("N","T"), #P3
#   corr == 0, #T1
#   cond %in% c(NA, "  ", " E", "") #T2
# ) %>% select(date, volume, price)
# 
# save(SP500_2008_1_cl, file = "SP500_2008_1_cl.Rdata", compress = "xz")
# rm("sp500_2008_1"); rm("SP500_2008_1_cl")
# 
# 
# load("sp500_2008_2.Rdata")
# SP500_2008_2_cl <- sp500_2008_2 %>% filter(
#   ex %in% c("N","T"), #P3
#   corr == 0, #T1
#   cond %in% c(NA, "  ", " E", "") #T2
# ) %>% select(date, volume, price)
# 
# save(SP500_2008_2_cl, file = "SP500_2008_2_cl.Rdata", compress = "xz")
# rm("sp500_2008_2"); rm("SP500_2008_2_cl")
# 
# 
# load("sp500_2008_3.Rdata")
# SP500_2008_3_cl <- sp500_2008_3 %>% filter(
#   ex %in% c("N","T"), #P3
#   corr == 0, #T1
#   cond %in% c(NA, "  ", " E", "") #T2
# ) %>% select(date, volume, price)
# 
# save(SP500_2008_3_cl, file = "SP500_2008_3_cl.Rdata", compress = "xz")
# rm("sp500_2008_3"); rm("SP500_2008_3_cl")
# 
# 
# load("sp500_2008_4.Rdata")
# SP500_2008_4_cl <- sp500_2008_4 %>% filter(
#   ex %in% c("N","T"), #P3
#   corr == 0, #T1
#   cond %in% c(NA, "  ", " E", "") #T2
# ) %>% select(date, volume, price)
# 
# save(SP500_2008_4_cl, file = "SP500_2008_4_cl.Rdata", compress = "xz")
# rm("sp500_2008_4"); rm("SP500_2008_4_cl")
# 
# load("sp500_2008_1_cl.Rdata")
# load("sp500_2008_2_cl.Rdata")
# load("sp500_2008_3_cl.Rdata")
# load("sp500_2008_4_cl.Rdata")
# 
# 
# SP500_2008 <- rbind(SP500_2008_1_cl, SP500_2008_2_cl, SP500_2008_3_cl, SP500_2008_4_cl)
# 
# save(SP500_2008, file = "SP500_2008_total_cl.Rdata", compress = "xz")
# 
# load(file = "SP500_2008_total_cl.Rdata")
# 
# SP500_2008_cleaned <- SP500_2008 %>% group_by(date) %>% summarise(    #T3
#   price = sum(price * volume)/sum(volume)  #weighted average price
# ) %>% mutate(roll_median=rollmedian(price, 51,  align = "center",na.pad = TRUE
# ))%>%mutate(MAD=zoo::rollmean(abs(price-roll_median), 51,align =
#                                 "center",na.pad = TRUE)) %>%
#   filter(price <= 10*MAD+roll_median,
#          price >= -10*MAD+roll_median)
# 
# save(SP500_2008_cleaned, file = "SP500_2008_cleaned.Rdata", compress = "xz")

# Bauwens et al. (2012) Synchronization -----------------------------------
load("SP500_2008_cleaned.Rdata")

SP500_2008_cleaned_n_synced <- sync_data(data = SP500_2008_cleaned, tick_time = "1 min")

trading_days <- as.Date(SP500_2008_cleaned$date, tz = "CET") %>% unique()
SP500_2008_cleaned_n_synced <- SP500_2008_cleaned_n_synced %>% filter(as.Date(Date, tz = "CET") %in% trading_days)

save(SP500_2008_cleaned_n_synced, file = "Code_Library/sp500_2008_cleaned_n_sync.Rdata", compress = "xz")
