
# Setup -------------------------------------------------------------------

source("functions.R")

# Only February ----------------------------------------------------------


#Get the synchronized data
#DB <- dbConnect(RSQLite::SQLite(), dbname = paste0("Application./SPY_2008_Februar_us.db"))
DB <- dbConnect(RSQLite::SQLite(), dbname = paste0("Application./SPY_2008_Marts.db"))


frequens=c(1,seq(5,60,5),seq(120,600,60))%>%lapply(function(i){
  paste(i, "sec")
})%>%do.call(what=rbind)

Data=lapply(1:22,function(i){
  tbl(DB,frequens[i])%>%collect()   
}
)



dbDisconnect(DB)

#set up vector witl all the dates.
dates=Data[1]%>%as.data.frame()%>%mutate("Date2"=as.Date(Date))%>%
  select(Date2)%>%unique()





#calculate the Rv at different frequencies.
RV=c()

for(i in 1:22){
  data=Data[i]%>%as.data.frame()
  data_new=data.frame(data,"Date2"=as.Date(data$Date))
  RV_day=c()
  
  for(j in 1:length(dates$Date2)){
    data_day=filter(data_new,Date2==dates$Date2[j])
    RV_day[j]=sum(diff(log(data_day$Price))^2)
  }
  
  RV[i]=mean(RV_day)
}

plot(RV)




# Whole 2007 --------------------------------------------------------------


#get the data
DB <- dbConnect(RSQLite::SQLite(), dbname = paste0("Data/Synced_data_2007.db"))
dbListTables(DB)

data=tbl(DB,"SPY_2007_1sec")%>%collect()

dbDisconnect(DB)

#Plot af det cleaned data.
index=seq(1,length(data_test$Price),300)

frequens=c(1,seq(5,60,5),seq(120,900,60))

#calculate the Rv at different frequencies.
RV=c()
for(i in 1:length(frequens)){
  index=seq(1,length(data$Price),frequens[i])
  data_new=data[index,]
  RV[i]=sum(diff(log(data_new$Price))^2)
}

plot(frequens,RV)

frequens=c(1,seq(5,60,5),seq(120,900,60))%>%lapply(function(i){
  paste(i, "sec")
})%>%do.call(what=rbind)

df=data.frame(x=1:length(frequens),y=RV)

ggplot(df,aes(x=x,y=y))+ geom_point()+scale_x_continuous(breaks = 1:length(frequens),label=frequens)+
  xlab("Frequency") + ylab("Realized Variance")+ theme(axis.text.x = element_text(angle = 50, hjust = 1))


#Estimate of sigma_epsilon using synched data

sigma=sqrt(sum(diff(log(data$Price))^2)/(2*length(diff(log(data$Price)))));sigma


#Using unsynched data
data_us = get_DB_data("2007-01-01", "2007-12-31")

sigma_us = sqrt(sum(diff(log(data_us$price))^2)/(2*length(diff(log(data_us$price)))));sigma_us




