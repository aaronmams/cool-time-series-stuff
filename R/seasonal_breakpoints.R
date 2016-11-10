##################################################
##################################################
##################################################
##################################################
#Script to analyze structural breaks with seasonal data
# under different conditions

library(strucchange)
library(dplyr)
library(ggplot2)
library(lubridate)

##################################################
##################################################
##################################################
##################################################
##################################################

#first simulate a monthly time-series with a 
# constant level and structural break in year 
# 5 of a ten  year time series

monthly.df <- tbl_df(data.frame(year=rep(2000:2010,12),month=c(rep(1,11),rep(2,11),rep(3,11),
                                                        rep(4,11),rep(5,11),rep(6,11),
                                                        rep(7,11),rep(8,11),rep(9,11),
                                                        rep(10,11),rep(11,11),rep(12,11)))) %>%
              mutate(jan=ifelse(month==1,1,0),feb=ifelse(month==2,1,0),march=ifelse(month==3,1,0),
                     april=ifelse(month==4,1,0),may=ifelse(month==5,1,0),june=ifelse(month==6,1,0),
                     july=ifelse(month==7,1,0),aug=ifelse(month==8,1,0),sept=ifelse(month==9,1,0),
                     oct=ifelse(month==10,1,0),nov=ifelse(month==11,1,0),dec=ifelse(month==12,1,0)) %>%
              mutate(beta1=ifelse(year<2005,10,20),
                     beta2=ifelse(year<2005,20,40),
                     beta3=ifelse(year<2005,30,60),
                     beta4=ifelse(year<2005,40,80),
                     beta5=ifelse(year<2005,50,100),
                     beta6=ifelse(year<2005,60,120),
                     beta7=ifelse(year<2005,70,140),
                     beta8=ifelse(year<2005,80,160),
                     beta9=ifelse(year<2005,90,180),
                     beta10=ifelse(year<2005,100,200),
                     beta11=ifelse(year<2005,25,50),
                     beta12=ifelse(year<2005,15,30),e=rnorm(132,10,10)) %>%
              mutate(trips=(jan*beta1)+(feb*beta2)+(march*beta3)+(april*beta4)+(may*beta5)+
                       (june*beta6)+(july*beta7)+(aug*beta8)+(sept*beta9)+
                       (oct*beta10)+(nov*beta11)+(dec*beta12)+e) %>%
              mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d"))

ggplot(monthly.df,aes(x=date,y=trips)) + geom_line() + geom_point()

fstats <- Fstats(trips~factor(month),data=monthly.df)
bp <- breakpoints(fstats)


