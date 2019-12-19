source("functions.R")

#Get the data
DB <- dbConnect(RSQLite::SQLite(), dbname = paste0("Data/Synced_data_2007.db"))
dbListTables(DB)

data=tbl(DB,"SPY_2007_1sec")%>%collect()

dbDisconnect(DB)

dates=data$Date%>%as.Date()%>%unique()

#add date variable.
data_new=data%>%mutate(Date2=as.Date(Date))


#calculate average daily log price

r=c()

for( i in 1:length(dates)){
  data_day=filter(data_new,Date2==dates[i])
    r[i]=mean(log(data_day$Price))
}


IV_BPV_avg_5sec = readRDS("Application/IV_BPV_avg_5sec.RDS")
jump_variation_result_2007 = readRDS("Application/jump_variation_result_2007.RDS")

#Set up various series.
IV=IV_BPV_avg_5sec
IV_prop=IV_BPV_avg_5sec$IV/jump_variation_result_2007$QV
IV_diff=(abs(IV$IV-IV$upper)+abs(IV$IV-IV$lower))/(2*IV$IV)



#Remoeves extreme outliers.
IV[which(IV_prop>1),]=NA
r[which(IV_prop>1)]=NA
dates[which(IV_prop>1)]=NA
IV_diff[which(IV_prop>1)]=NA
IV_prop[which(IV_prop>1)]=NA

library("reshape2")
min(as.numeric(IV$sd))

#Plot the average deviation from the confidence intervals.
df=data.frame(x=dates,y=IV_diff)%>%na.omit()

ggplot(df,aes(x=x,y=y))+geom_line()+
  xlab("Date") + ylab("")

#Plot the daily log price.
df=data.frame(x=dates,y=r)%>%na.omit()

ggplot(df,aes(x=x,y=y))+ geom_line()+xlab("Date") + ylab("Daily Log Price")

#Plot the daily log returns.
df=data.frame(x=dates[-length(dates)],y=diff(r))%>%na.omit()

ggplot(df,aes(x=x,y=y))+ geom_line()+xlab("Date") + ylab("Daily Log Returns")

#plot the Iv standard deviation.
df=data.frame(x=dates,y=IV$sd)%>%na.omit()

ggplot(df,aes(x=x,y=y))+ geom_line()+xlab("Date") + ylab("IV standard deviation")

#plot the IV
df=data.frame(x=dates,y=IV$IV)%>%na.omit()

ggplot(df,aes(x=x,y=y))+geom_line()+
  xlab("Date") + ylab(expression(IV))

#Plot the IV along with is confidence intervals.
df=data.frame(date=dates,IV=IV$IV,lower=IV$lower,upper=IV$upper)%>%na.omit()%>%melt(id="date")

ggplot(df,aes(x=date,y=value,colour=variable))+geom_line()+
  scale_colour_discrete(name = "", labels = c("IV", "2.5%", "97.5%"))+
  xlab("Date") + ylab(expression(IV))


#Plot the proportion of QV due to IV.
df=data.frame(x=dates,y=IV_prop)%>%na.omit()

ggplot(df,aes(x=x,y=y))+ geom_line()+xlab("Date") + ylab(expression(over(IV, QV)))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))



