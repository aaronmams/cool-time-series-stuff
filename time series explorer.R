#Currently the script covers material related to the following blog posts:

#http://thesamuelsoncondition.com/2016/01/16/time-series-part-i-some-popular-diagnostics/
#http://thesamuelsoncondition.com/2016/01/22/time-series-ii-autocorrelation-serial-correlation-error-structures/

#Necessary updates include:

#1. expirimenting with the different structural break and regime switching models for
#    seasonal time-series data: MSwM, bfast, strucchange packages, etc.s

#2.  improving the data generation function to be a bit more general.

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


#In this script we set up a hands on practical for illustrating some important
# time-series methods

#First we simulate two pretty simple series with some interesting time-series
# properties.

#Then we run through some popular time-series diagnostics for stationarity,
# seasonality, and autoregressive order

#Finally, we try to apply some formal models to the two series to detect
# structural change or regime change.

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#Preamble
require(MSwM)
require(ggplot2)
require(dplyr)
require(bfast)
require(data.table)
require(strucchange)
require(lmtest)
require(lawstat)
require(tseries)
require(matrixcalc)
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

#Generate Data
# write the data generation as a function which will allow us to 
# either 
#1. read from the existing data and reproduce the analysis
#    on the blog, or
#2.  generate some new data to play with

#generate 2 really simple monthly process
seasonal.data <- function(mu,sigma1,annual.total=12000,
                          shares.before=c(0.02,0.03,0.04,0.05,0.05,0.06,0.15,0.2,0.20,0.1,0.05,0.05),
                          shares.after=c(0.01,0.01,0.02,0.02,0.02,0.05,0.1,0.2,0.30,0.2,0.05,0.02),
                          sigma2=75){

#inputs
  # mu    the annual mean of the first series...the first series is generated according to a 
  #       process that spreads the annual total evenly among all the months before the break, then
  #       after the break spreads 25% of the annual total evenly among 9 months and spreads 75% of
  #       the annual total among the 3 summer months.  Note that this series is assumed stationary
  #       in the sense that the annual total doesn't change after the break.  So once the parameter mu
  #       is set, the monthly mean values for every month in the series can be recovered.
  
  # sigma1   a 3X1 vector of variance parameters.  
  #           sigma1[1]: noise around the monthly mean values in the first part of the first series.
  #           sigma[2]: noise around the monthly mean values after the break for the 9 months where
  #                     where 25% of the annual total is evenly distributed among 9 months.
  #           sigma[3]: noise around the monthly mean values after the break for the 3 months where 75% 
  #                     of the annual total is evenly distributed over 3 months.
  
  # annual.total  the annual total for the 2nd series that will be apportioned among the months
  
  # shares.before   a 1X12 vector defining the seasonality of the 2nd series in the period before the break.
  #                   values correspond to the % of the annual total in each month
  
  # shares.after    a 1X12 vector defining the seasonality of the 2nd series in the period after the break.
  #                   values correspond to the % of the annual total in each month
  
  # sigma2          a numeric constant defining the noise around the monthly means for the 2nd series.  In
  #                   generating the 2nd series the variance for each monthly observation is taken to be a 
  #                   constant fraction of the mean.  The value sigma2 defines this fraction.  Larger values will 
  #                   produce less noise in the series and smaller values will produce series with more noise.
  
  #6 years worth of approximately evenly distributed landings
  year <- c(rep(2005,12),rep(2006,12),rep(2007,12),rep(2008,12),rep(2009,12),rep(2010,12))
  month <- c(rep(1:12,6))
  y <- rnorm(72,mu/12,sigma1[1])
  df <- data.frame(year=year,month=month,lbs=y) 

mu2 <- mu*0.25
mu3 <- mu*0.75
  #4 years worth of landings that are evenly distributed throughout 9 months and
  # spike up considerably in the months June-August
  year2 <- c(rep(2011,12),rep(2012,12),rep(2013,12),rep(2014,12))
  month2 <- c(rep(1:12,4))
  
  #9 month X 4 years of data distributed N(333.333,0.666)
  month1 <- data.frame(lbs=c(rnorm(36,mu2/9,sigma1[2])),month=rep(c(1,2,3,4,5,9,10,11,12),4),
                       year=c(rep(2011,9),rep(2012,9),rep(2013,9),rep(2014,9)))
  #3 months X 4 years of data distributed N(3000,6)
  month2 <- data.frame(lbs=c(rnorm(12,mu3/3,sigma1[3])),month=rep(c(6,7,8),4),
                       year=c(rep(2011,3),rep(2012,3),rep(2013,3),rep(2014,3)))
  
  df2 <- data.frame(rbind(month1,month2))
  
  z <- tbl_df(rbind(df,df2)) 
  #====================================================================================
  #====================================================================================
  #====================================================================================
  #====================================================================================
  #next is a series that is seasonal in nature in both periods but 
  # there is an abrupt shift in the intensity of seasonality after the 
  # 6th year
  
  year <- c(rep(2005,12),rep(2006,12),rep(2007,12),rep(2008,12),rep(2009,12),rep(2010,12))
  month <- c(rep(1:12,6))
  
  y <- data.frame(year=year,month=month)
  
  total.lbs <- annual.total
  #mean.shares <- c(0.02,0.03,0.04,0.05,0.05,0.06,0.15,0.2,0.20,0.1,0.05,0.05)
  mean.shares <- shares.before
  mean.lbs <- mean.shares*total.lbs
  
  df <- rbind(data.frame(year=c(2005:2010),month=1,lbs=rnorm(6,mean.lbs[1],mean.lbs[1]/sigma2)),
              data.frame(year=c(2005:2010),month=2,lbs=rnorm(6,mean.lbs[2],mean.lbs[2]/sigma2)),
              data.frame(year=c(2005:2010),month=3,lbs=rnorm(6,mean.lbs[3],mean.lbs[3]/sigma2)),
              data.frame(year=c(2005:2010),month=4,lbs=rnorm(6,mean.lbs[4],mean.lbs[4]/sigma2)),
              data.frame(year=c(2005:2010),month=5,lbs=rnorm(6,mean.lbs[5],mean.lbs[5]/sigma2)),
              data.frame(year=c(2005:2010),month=6,lbs=rnorm(6,mean.lbs[6],mean.lbs[6]/sigma2)),
              data.frame(year=c(2005:2010),month=7,lbs=rnorm(6,mean.lbs[7],mean.lbs[7]/sigma2)),
              data.frame(year=c(2005:2010),month=8,lbs=rnorm(6,mean.lbs[8],mean.lbs[8]/sigma2)),
              data.frame(year=c(2005:2010),month=9,lbs=rnorm(6,mean.lbs[9],mean.lbs[9]/sigma2)),
              data.frame(year=c(2005:2010),month=10,lbs=rnorm(6,mean.lbs[10],mean.lbs[10]/sigma2)),
              data.frame(year=c(2005:2010),month=11,lbs=rnorm(6,mean.lbs[11],mean.lbs[11]/sigma2)),
              data.frame(year=c(2005:2010),month=12,lbs=rnorm(6,mean.lbs[12],mean.lbs[12]/sigma2))
  )
  
  y <- merge(y,df,by=c('year','month'))
  
  
  #Now change the process abruptly
  #mean.shares2 <- c(0.01,0.01,0.02,0.02,0.02,0.05,0.1,0.2,0.30,0.2,0.05,0.02)
  mean.shares2 <- shares.after
  mean.lbs2 <- mean.shares2*annual.total
  
  df <- rbind(data.frame(year=c(2011:2014),month=1,lbs=rnorm(4,mean.lbs2[1],mean.lbs2[1]/sigma2)),
              data.frame(year=c(2011:2014),month=2,lbs=rnorm(4,mean.lbs2[2],mean.lbs2[2]/sigma2)),
              data.frame(year=c(2011:2014),month=3,lbs=rnorm(4,mean.lbs2[3],mean.lbs2[3]/sigma2)),
              data.frame(year=c(2011:2014),month=4,lbs=rnorm(4,mean.lbs2[4],mean.lbs2[4]/sigma2)),
              data.frame(year=c(2011:2014),month=5,lbs=rnorm(4,mean.lbs2[5],mean.lbs2[5]/sigma2)),
              data.frame(year=c(2011:2014),month=6,lbs=rnorm(4,mean.lbs2[6],mean.lbs2[6]/sigma2)),
              data.frame(year=c(2011:2014),month=7,lbs=rnorm(4,mean.lbs2[7],mean.lbs2[7]/sigma2)),
              data.frame(year=c(2011:2014),month=8,lbs=rnorm(4,mean.lbs2[8],mean.lbs2[8]/sigma2)),
              data.frame(year=c(2011:2014),month=9,lbs=rnorm(4,mean.lbs2[9],mean.lbs2[9]/sigma2)),
              data.frame(year=c(2011:2014),month=10,lbs=rnorm(4,mean.lbs2[10],mean.lbs2[10]/sigma2)),
              data.frame(year=c(2011:2014),month=11,lbs=rnorm(4,mean.lbs2[11],mean.lbs2[11]/sigma2)),
              data.frame(year=c(2011:2014),month=12,lbs=rnorm(4,mean.lbs2[12],mean.lbs2[12]/sigma2))
  )
  
  y2 <- data.frame(month=c(rep(1:12,4)),year=c(rep(2011,12),rep(2012,12),rep(2013,12),rep(2014,12)))
  y2 <- merge(y2,df,by=c('year','month'))
  
  z2 <- tbl_df(rbind(y,y2))
  
  #===============================================================================
  #===============================================================================
  #===============================================================================
  #===============================================================================
  #===============================================================================
  
  #name each series
  z$series <- 1
  z2$series <- 2
  
  #combine in one data frame
  z <- rbind(z,z2)
  
  #add monthly dummy variables and make some derived columns
  z <- z %>% mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
    mutate(jan=ifelse(month==1,1,0),feb=ifelse(month==2,1,0),mar=ifelse(month==3,1,0),
           apr=ifelse(month==4,1,0),may=ifelse(month==5,1,0),jun=ifelse(month==6,1,0),
           jul=ifelse(month==7,1,0),aug=ifelse(month==8,1,0),sep=ifelse(month==9,1,0),
           oct=ifelse(month==10,1,0),nov=ifelse(month==11,1,0),dec=ifelse(month==12,1,0),
           regime=ifelse(year>2010,1,0),regime2=ifelse(year>2009,1,0),regime3=ifelse(year>2008,1,0)) %>%
    group_by(series) %>% arrange(year,month,series) %>%
    mutate(lag1 = lag(lbs,1),diff=lbs-lag1,log.diff=log(lbs)-log(lag1),
           lag12=lag(lbs,12),diff12=lbs-lag12,g12 = log(lbs)-log(lag12)) %>%
    group_by(year,series) %>% mutate(share=lbs/sum(lbs),annualmean=mean(lbs)) %>% arrange(year,month) %>%
    group_by(series) %>% mutate(meanlbs=mean(lbs))
  return(z)
}




#====================================================================================
#====================================================================================
#====================================================================================
#====================================================================================
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################



#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

#here you can choose to run the simulation function above to generate a new 
# set of data...or you can read from a saved data frame to work with the same
# data I used in the blog post
z <- seasonal.data(mu=12000,sigma1=c(10,10,30))

#save the data for replication
#saveRDS(z,file="test_data.RDA")


#UNCOMMENT NEXT LINE TO WORK WITH THE SAME DATA SET I USE IN THE BLOG POST
z <- readRDS(file="test_data.RDA")

#plot the two series
ggplot(z,aes(x=date,y=lbs)) + geom_line()  + facet_wrap(~series) + 
  geom_hline(aes(yintercept=meanlbs),data=z,color='red') + 
  ylab('generic y')

#plot differenced 2nd series
ggplot(subset(z,series==2 & year < 2011),aes(x=date,y=g12)) + geom_line()

#==================================================================
#==================================================================
#==================================================================
#==================================================================
#==================================================================

#==================================================================
#note that the mean annual landings statistically unchanged
avglbs <- z %>% group_by(year,series) %>% summarise(totallbs=sum(lbs))
t.test(avglbs$totallbs[which(avglbs$series==1)][1:6],avglbs$totallbs[which(avglbs$series==1)][7:10])

t.test(avglbs$totallbs[which(avglbs$series==2)][1:6],avglbs$totallbs[which(avglbs$series==2)][7:10])

#==================================================================

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
#first a breif digression on what stationary and non-station series look like

#stationary AR(1) process
set.seed(1)
x <- w <- rnorm(120)
for (t in 2:120) x[t] <- (x[t-1]*0.25)+w[t]

#a trend stationary process
set.seed(1)
x.trend <- w.trend <- rnorm(120)
for (t in 2:120)x.trend[t] <- (x.trend[t-1]*0.25)+w.trend[t] + 0.01*t

#a random walk
set.seed(1)
y <- eps <- rnorm(120)
for (t in 2:120)y[t] <- (y[t-1])+eps[t]

#a random walk with drift
set.seed(1)
rwd <- rw <- rnorm(120)
for (t in 2:120)rwd[t] <- (rwd[t-1])+rw[t]+0.05


example <- data.frame(rbind(data.frame(t=c(1:120),x=x,type="AR(1) stationary"),
                            data.frame(t=c(1:120),x=y,type="Random Walk"),
                            data.frame(t=c(1:120),x=x.trend,type="Trend Stationary"),
                            data.frame(t=c(1:120),x=rwd,type="Random Walk with Drift")))

ggplot(example,aes(x=t,y=x)) + geom_line() + facet_wrap(~type,scales="free")

#how about a stochastic stationary seasonal process?
require(stats)
#use the arima.sim function to simulate:
# 1). a seasonal AR(1) process with monthly frequency
seasonal.ar1.sim <- arima.sim(list(order = c(12,0,0), ar = c(rep(0,11),0.9)), n = 120)

# 2). a basic AR(1) process with the same ar parameter for comparison
ar1.sim <- arima.sim(list(order = c(1,0,0), ar = 0.9), n = 120)

#combine them for plotting
s.tmp <- data.frame(rbind(data.frame(t=c(1:120),x=c(seasonal.ar1.sim),type="Seasonal AR(1)"),
                          data.frame(t=c(1:120),x=c(ar1.sim),type="Non-Seasonal AR(1)")))

#plot the seasonal AR(1) and the regular AR(1)
ggplot(s.tmp,aes(x=t,y=x)) + geom_line() + facet_wrap(~type) 

#add a month column
s.tmp$month <- rep(c('jan','feb','mar','apr','may','june','july','aug','sept','oct','nov','dec')
                   ,10)

#run the dummy variable model with the seasonal AR(1) data
det.season <- lm(x~factor(month),
                 data=s.tmp[which(s.tmp$type=="Seasonal AR(1)"),])

#fit a seasonal AR(1) model to the simulated data
tsx <- s.tmp$x[which(s.tmp$type=="Seasonal AR(1)")]
tsx <- ts(tsx,start(c(2000,12)),frequency=12)
ar1.season <- arima(tsx,
                    order=c(0,0,0),
                    seasonal=list(order=c(1,0,0),period=12))


#Compare the RMSE of the dummy variable regression v. the seasonal AR(1) model
# for modeling the simulated seasonal AR(1) data
lm.df <- tbl_df(data.frame(y=s.tmp$x[which(s.tmp$type=="Seasonal AR(1)")],
                           yhat=det.season$fitted.values)) %>% 
  mutate(error.sq=(y-yhat)^2) %>%
  summarise(mean.errorsq=mean(error.sq)) %>% mutate(RMSE=sqrt(mean.errorsq))


#I don't want to mess with the math of simulating a non-stationary stochastic 
# seasonal process...but can we draw one?

# if the seasonal peak of a process were evolving over time (maybe due to climate change)
# we might end up with something 
x.ns <- w.ns <- rnorm(120)
for (t in 13:120) x.ns[t] <- (x.ns[t-12])+w.ns[t]

df.ns <- data.frame(t=c(1:120),x=x.ns,month=rep(c('jan','feb','mar','apr','may','jun','jul',
                                                  'aug','sep','oct','nov','dec'),10),
                    year=c(rep(1,12),rep(2,12),rep(3,12),rep(4,12),rep(5,12),rep(6,12),rep(7,12),rep(8,12),
                           rep(9,12),rep(10,12)))

ggplot(df.ns,aes(x=t,y=x)) + geom_line() + geom_vline(xintercept=c(12,24,36,48))

#plot the periodogram
spec <- spec.pgram(df.ns$x,log="yes")
ggplot(data.frame(freq=spec$freq,spec=spec$spec),aes(x=freq,y=spec)) + 
  geom_line() + theme_bw() + ggtitle("Periodogram of Monthly Fake Data")
#identify the period of the first spike (location and period)
1/spec$freq[which(spec$spec==max(spec$spec))]

#------------------------------------------------------------------------------------
#try this a different way, start with a deterministic seasonal process like series 2
# here we have some monthly landings distributed throughout the year with a peak around
# late summer early fall.  We might imagine a climate process that, because it is tending to
# get warmer earlier in the year, fishing activity is tending to move toward the late winter/spring months.
# the following code generates a process that transitions from a seasonal summer peak to a seasonal
# early spring peak
lbs <- c(200,250,250,300,300,325,325,350,325,325,300,300)
phi <- c(1.02,1.02,1.03,1.02,1.02,0.98,0.98,0.98,0.98,0.985,0.98,0.99)
data <- data.frame(year=rep(2000,12),month=c(1:12),lbs=lbs)
for(t in 2001:2009){
  eps <- rnorm(12,0,12)
  lbs <- (lbs*phi) + eps
  dnow <- data.frame(year=rep(t,12),month=c(1:12),lbs=lbs)
  data <- rbind(data,dnow)
}
data$date <- as.Date(paste(data$year,"-",data$month,"-","01",sep=""),format="%Y-%m-%d")
data$t <- c(1:120)
ggplot(data,aes(x=t,y=lbs))+geom_line() + geom_point() + geom_vline(xintercept=c(1,13,25,97,109),color="red") +
  annotate("text",x=13,y=375,
           label=as.character("Jan 01")) + 
  annotate("text",x=25,y=375,
           label="Jan 02") + 
  annotate("text",x=97,y=375,
           label="Jan 08") + 
  annotate("text",x=109,y=375,
           label="Jan 09")  
  
#-------------------------------------------------------------------------------------------


######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################



#Some time-series diagnostics.  Here we will plot the series in levels, differences,
# and logs

# We will use a roll-your-own type procedure to examine whether the mean
# of the series changes depending on the window of time the mean is
# calculated over

# to address seasonality we will plot the periodogram


#==================================================================
#Step 1: plot the series
#plot the series in 1st differences, and log differences
ggplot(z,aes(x=date,y=diff)) + geom_line() + facet_wrap(~series,scales='free') + 
  ggtitle("Some fake tourism data") + 
  theme(plot.title = element_text(lineheight=.8))  

#plot the log difference
ggplot(z,aes(x=date,y=log.diff)) + geom_line() + facet_wrap(~series,scales='free') + 
  ggtitle("Some fake tourism data") + 
  theme(plot.title = element_text(lineheight=.8))  

#plot the 12th difference 
ggplot(z,aes(x=date,y=diff12)) + geom_line() + 
  facet_wrap(~series) + ggtitle("Twelth order difference") + theme(plot.title=element_text(lineheight=0.8))

#===================================================================


#=====================================================================
# a function to test for equality of means calculated over different windows
# of the time-series

#we'll focuse here on series 1 for a minute
s1 <- z[which(z$series==1),]

welch.t <- function(i,n){
  x <- s1$lbs[i[1]:n[1]]
  y <- s1$lbs[i[2]:n[2]]
  welch.t <- (mean(x)-mean(y))/(sqrt((var(x)/length(x))+(var(y)/length(y))))
  return(data.frame(xbar=mean(x),ybar=mean(y),welch.t=welch.t,xstart=s1$date[i[1]],
                    xend=s1$date[n[1]],ystart=s1$date[i[2]],yend=s1$date[n[2]]))
}

yrs.even <- data.frame(start=c(1,12,72),end=c(120,48,96))
yrs <- data.frame(start=c(48,1,1,15,15,15,78,111),end=c(120,41,77,120,31,88,81,120))

meanstest <- list()
counter=0
for(k in 1:nrow(yrs)){
  for(l in 1:nrow(yrs)){
    counter= counter+1
    meanstest[[counter]]<-  welch.t(i=c(yrs$start[k],yrs$start[l]),
                                    n=c(yrs$end[k],yrs$end[l]))
  }
}

#plot the series along with some paired windows for which the mean is not statistically
# different along with some for which it is.

meanstest <- data.frame(rbindlist(meanstest))
plot.tmp <- meanstest[54,]

s1$t <- seq(1:nrow(s1))

#this is pretty contrived but since I know the 54th element of the
# "meanstest" object test for a significant difference between the
# mean of the series calculated from 2011/06/01-2011/08/01 and 
# from 2008/01/01-2014/01/01, I'm just going to pull that element out
# and illustrate the different means.
ggplot(s1,aes(x=as.numeric(t),y=lbs)) + geom_line() +
  xlab("t") + ylab("fake data") + 
  geom_segment(aes(x = 78, y = 2000, xend = 81, yend = 2000,
                   colour = "segment"),data = plot.tmp,
               arrow = arrow(length = unit(0.03, "npc"))) + 
  geom_segment(aes(x = 15, y = 998, xend = 88, yend = 998, 
                   colour = "segment",lineend='square'), data = plot.tmp,
               arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text",x=70,y=2000,
           label=as.character(expression(~mu==2329)),parse=T) + 
  annotate("text",x=42,y=960,
           label=as.character(expression(~mu==963)),parse=T)  +
  scale_color_discrete(guide=FALSE)



#formal tests
kpss.test(z$lbs[which(z$series==1)])
kpss.test(z$lbs[which(z$series==1)],null="Trend")
kpss.test(z$lbs[which(z$series==2)])
kpss.test(z$lbs[which(z$series==2)],null="Trend")

adf.test(z$lbs[which(z$series==1)],k=1)
adf.test(z$lbs[which(z$series==1)],k=12)
adf.test(z$lbs[which(z$series==1)])

adf.test(z$lbs[which(z$series==2)],k=1)
adf.test(z$lbs[which(z$series==2)],k=12)
adf.test(z$lbs[which(z$series==2)])

#=====================================================================
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#Autocorrelation Structure

#----------------------------------------------------------------------
# what does it mean to have autocorrelated errors?

#let revisit our AR(1) with trend process
set.seed(1)
x.trend <- w.trend <- rnorm(120)
for (t in 2:120)x.trend[t] <- (x.trend[t-1]*0.25)+w.trend[t] + 0.1*t
tmp <- data.frame(t=c(1:120),x=x.trend)
ggplot(tmp,aes(x=t,y=x)) + geom_line() + stat_smooth(method='lm')

tmp$pred <- lm(x~t,data=tmp)$fitted.values
tmp$error <- tmp$x-tmp$pred
ggplot(tmp,aes(x=t,y=error)) + geom_line()

#are the errors correlated
tmp <- tbl_df(tmp) %>% mutate(e1=lag(error))

cor.test(tmp$error[2:nrow(tmp)],tmp$e1[2:nrow(tmp)])

#ghetto runs test...are positive values followed by positive values?
tmp <- tmp %>% mutate(run=ifelse(error>0,1,ifelse(error==0,0,-1)))
ggplot(tmp,aes(x=t,y=run)) + geom_line()

#---------------------------------------------------------------------


#----------------------------------------------------------------------
#Graphical Examination: the autocorrelation and partial autocorrelation
# functions
acf(z$lbs[which(z$series==1)])
acf(z$lbs[which(z$series==2)])

#series 1, the non-seasonal portion
acf(z$lbs[which(z$series==1 & z$year < 2011)])
# series 1, the really seasonal part
acf(z$lbs[which(z$series==1 & z$year > 2010)])


#partial acf
pacf(z$lbs[which(z$series==1)])
pacf(z$lbs[which(z$series==2)])

#series 1, the non-seasonal portion
pacf(z$lbs[which(z$series==1 & z$year < 2011)])
# series 1, the really seasonal part
pacf(z$lbs[which(z$series==1 & z$year > 2010)])

#------------------------------------------------------------------------


# Formal Test: Brusch-Godfrey Test
# for BG test we estimate the auxiliary regression:

# y(t) = alpha + beta*DM + rho1*u(t-1) + rho2*u(t-1)

#  where u() are the residuals from the original regression of y on the 
# monthly dummy variables

#----------------------------------------------------------------
#start with the first series and only the time before the break
df.tmp <- seasonal.data(sigma1=10,sigma2=100,nmonth1=72)
df.tmp <- df.tmp[which(df.tmp$series==2 & df.tmp$year<2011),]
ggplot(df.tmp,aes(x=date,y=lbs)) + geom_line()


mod.pre <- lm(lbs~feb+mar+apr+may+jun+jul+aug+sep+oct+nov+dec,
              data=df.tmp)
df.tmp$resid <- mod.pre$residuals
df.tmp$fitted <- mod.pre$fitted.values

ggplot(df.tmp,aes(x=date,y=resid))+geom_point()
ggplot(df.tmp,aes(x=date,y=lbs)) + geom_line() + geom_line(data=df.tmp,aes(x=date,y=fitted),color="red")



#now create the lagged residuals for the BG test
df.tmp <- df.tmp %>% mutate(r1=lag(resid,1),r2=lag(resid,2),
                            r3=lag(resid,3),r4=lag(resid,4),
                            r5=lag(resid,5),r6=lag(resid,6),
                            r7=lag(resid,7),r8=lag(resid,8),
                            r9=lag(resid,9),r10=lag(resid,10),
                            r11=lag(resid,11),r12=lag(resid,12))


BG.test <- lm(resid~feb+mar+apr+may+jun+jul+aug+sep+oct+nov+dec+r1+
                r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12,data=df.tmp)

#the relevant test statistic is:
#   (T-p)*R^2
# which is distributed approximately Chi-Squared(p) and test the null
# hypothesis that all rho = 0 (not autocorrelation)

BG.teststat <- (nrow(df.tmp)-12)*summary(BG.test)$r.squared

#chi-square critical values for p=12 (number of regressors in the model) 
# degrees of freedom is 18.5 for alpha = 0.9

#run the BG test in R's lmtest package
require(lmtest)
bgtest(mod.pre,order=12)
bgtest(mod.pre,order=1)
bgtest(df.tmp$lbs~1,order=12)

#compare to DW test
dwtest(mod.pre)

#runs test on the residuals
require(lawstat)
lawstat::runs.test(df.tmp$resid)
#-------------------------------------------------------------------------------------


#----------------------------------------------------------------------
#let's look at the errors of a regression on the first part of the
# second series
s2 <- z[which(z$series==2 & z$year<2011),]
s2$resid <- lm(lbs~1,data=s2)$residuals
lawstat::runs.test(s2$resid,plot.it=TRUE)

#----------------------------------------------------------------------


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################




