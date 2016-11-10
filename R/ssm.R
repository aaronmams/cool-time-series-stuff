#This script provides some examples of state-space modeling using the 
# KFAS package.

#Outline:

#1. Using data from the West Coast Groundfish Fishery, we estimate a structural 
#   time-series model in state-space form.  The structural time-series model decomposes
#   number of fishing trips per month from 1994-2014 into a level component and a 
#   seasonal component.

#1A. We compare this model to a linear seasonal model where the seasonality is captured 
#    montly dummy variables and we define a structural break in the series.  This is a bit
#    of a straw-man but the state-space model is really good at picking up smooth
#    changes in the seasonal effect....for purely illustrative purposes we want to 
#    evaluate the state-space model against something that (conceptually anyway) should
#    do a good job of picking up abrupt changes in the seasonal effect.

#2.  A simulation exercise.  We simulate some data from a seasonal time-series process
#    where the seasonal influences change dramatically at a certian point in time.  This 
#    is done as an illustration of how the state-space model performs in the 
#    presence of a discrete structural change in seasonality.

library(Quandl)
library(KFAS)
library(dplyr)
library(ggplot2)
library(lubridate)
library("ggfortify")
library(xts)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#A GAUSSIAN STATE-SPACE SEASONAL MODEL FOR FISHING TRIPS:

#Read in the groundfish effort data------------------------------
df.monthly <- readRDS('data/gf_monthly.RDA')
df.monthly <- readRDS('/Users/aaronmamula/Downloads/gf_monthly.RDA')
#----------------------------------------------------------------

#create a univariate ts object---------------------------------
trips <- df.monthly$ntrips
#--------------------------------------------------------------

#do some other convenience operations---------------------------
df.monthly$date <- as.Date(paste(df.monthly$year,"-",df.monthly$month,"-","01",sep=""),format="%Y-%m-%d")
dates <- df.monthly$date
trips <- xts(trips,order.by=dates)

df.monthly <- df.monthly %>% mutate(july=ifelse(month==7,1,0))
#-----------------------------------------------------------------

#-------------------------------------------------------------------
# a few quick illustative plots
ggplot(df.monthly,aes(x=date,y=ntrips)) + geom_line() + 
  geom_point(aes(color=factor(july))) +  
  scale_color_manual(values=c('red','black')) +
  theme_bw() +
  theme(legend.position='none') 

ggplot(subset(df.monthly,area=='north'),aes(x=date,y=tripshare)) + geom_line() + 
  geom_point(aes(color=factor(july)))
#---------------------------------------------------------------------

#-----------------------------------------------------------------------
# A state-space model with a local time-varying level and time varying
# seasonal effects:
tripsmodel<-SSModel(trips ~ SSMtrend(degree = 1, Q=list(matrix(NA))) + 
                       SSMseasonal(period=12, sea.type="dummy", Q = NA), H = NA)
str(tripsmodel)
#-----------------------------------------------------------------------

#------------------------------------------------------------------------
#Estimate the model:
tripsFit<-fitSSM(tripsmodel,inits=c(0.1,0.05, 0.001),method='BFGS')$model

#Recover the smoothed state estimates
tripsSmooth <- KFS(tripsFit,smooth= c('state', 'mean','disturbance'))

#Get smoothed estimates for the level
tripsLevel <-signal(tripsSmooth, states = 'level')
tripsLevel$signal <- xts(tripsLevel$signal, order.by = dates)

#plot level
autoplot(cbind(trips, tripsLevel$signal),facets = FALSE)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#plot fitted values with raw data
pred <- data.frame(date=dates,kfas=fitted(tripsSmooth),y=as.numeric(trips))

ggplot(pred,aes(x=date,y=y)) + geom_line() + 
  geom_line(aes(x=date,y=kfas),data=pred,color='red') +theme_bw()
#------------------------------------------------------------------------------

#---------------------------------------------------------------------------
#plot the seasonal
trips.sea <- as.numeric(tripsSmooth$alphahat[,2])
plot.df <- data.frame(date=dates,seasonal=trips.sea)
#add a reference point to the seasonal data frame
plot.df <- plot.df %>% mutate(month=month(date),july=ifelse(month==7,1,0))

ggplot(plot.df,aes(x=date,y=seasonal)) + geom_line() + 
  geom_point(aes(color=factor(july))) + scale_color_manual(values=c('red','black'))+
  theme_bw() + 
  theme(legend.position="none")

#----------------------------------------------------------------------------


#plot level
level <- as.numeric(tripsSmooth$alphahat[,1])
plot.df <- data.frame(date=dates,level=level,z='smoothed')
                
ggplot(plot.df,aes(x=date,y=level)) + geom_line() + 
  geom_line(data=subset(df.monthly),aes(x=date,y=ntrips,color="red"))


#==========================================================================
#COMPARISON OF THE STATE-SPACE MODEL WITH A DUMMY VARIABLE REGRESSION WITH
# A STRUCTURAL BREAK
#==========================================================================
#get mean absolute error from the KFAS model and compare it to the
# dummy variable model with structural break in 2010.

#the dummy variable model...I know this is clunky but I 
# prefer the traditional dummy variable set up.
chow.df <- df.monthly  %>%
            mutate(jan=ifelse(month==1,1,0),
                   feb=ifelse(month==2,1,0),
                   march=ifelse(month==3,1,0),
                   april=ifelse(month==4,1,0),
                   may=ifelse(month==5,1,0),
                   june=ifelse(month==6,1,0),
                   july=ifelse(month==7,1,0),
                   aug=ifelse(month==8,1,0),
                   sept=ifelse(month==9,1,0),
                   oct=ifelse(month==10,1,0),
                   nov=ifelse(month==11,1,0),
                   dec=ifelse(month==12,1,0))

model.pre <- lm(ntrips~feb+march+april+may+june+july+aug+
                  sept+oct+nov+dec+factor(year),data=subset(chow.df,year<=2010))
model.post <- lm(ntrips~feb+march+april+may+june+july+aug+
                   sept+oct+nov+dec+factor(year),data=subset(chow.df,year>2010))

yhat <- rbind(
data.frame(date=chow.df$date[chow.df$year<=2010],
              that=predict(model.pre,newdata=chow.df[chow.df$year<=2010,])),
data.frame(date=chow.df$date[chow.df$year>2010],
           that=predict(model.post,newdata=chow.df[chow.df$year>2010,]))
)

yhat$y <- as.numeric(trips)
yhat$model <- 'Seasonal Dummy'

ggplot(yhat,aes(x=date,y=y)) + geom_line() + geom_line(aes(x=date,y=that,color='red'))

# the KFAS model
kfas.pred <- data.frame(date=dates,that=fitted(tripsSmooth),y=as.numeric(trips),model='KFAS')

#combine the two
model.comp <- rbind(yhat,kfas.pred)

model.comp$eps <- model.comp$y-model.comp$that

ggplot(model.comp,aes(x=eps)) + geom_histogram() + facet_wrap(~model)
ggplot(model.comp,aes(x=date,y=eps,color=model)) + geom_line() + theme_bw()


#get mean absolute in-sample prediction error:
tbl_df(model.comp) %>% mutate(abs.error=abs(eps)) %>% group_by(model) %>%
              summarise(mean(abs.error,na.rm=T))

#==========================================================================


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
#A SIMULATION EXPERIMENT:

#----------------------------------------------------------------
#simulate data that is seasonal but with seasonality that shifts at some point
means <- c(10,20,30,40,50,60,70,40,30,20,20,10)
means2 <- c(70,60,50,40,30,20,10,20,30,40,60,70)

df <- rbind(
  data.frame(trips.s=rnorm(10,means[1],5),year=c(1995:2004),month=1),
  data.frame(trips.s= rnorm(10,means[2],5),year=c(1995:2004),month=2),
  data.frame(trips.s=rnorm(10,means[3],5),year=c(1995:2004),month=3),
  data.frame(trips.s=rnorm(10,means[4],5),year=c(1995:2004),month=4),
  data.frame(trips.s=rnorm(10,means[5],5),year=c(1995:2004),month=5),
  data.frame(trips.s=rnorm(10,means[6],5),year=c(1995:2004),month=6),
  data.frame(trips.s=rnorm(10,means[7],5),year=c(1995:2004),month=7),
  data.frame(trips.s=rnorm(10,means[8],5),year=c(1995:2004),month=8),
  data.frame(trips.s=rnorm(10,means[9],5),year=c(1995:2004),month=9),
  data.frame(trips.s=rnorm(10,means[10],5),year=c(1995:2004),month=10),
  data.frame(trips.s=rnorm(10,means[11],5),year=c(1995:2004),month=11),
  data.frame(trips.s=rnorm(10,means[12],5),year=c(1995:2004),month=12)
)

df <- tbl_df(df) %>% mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month) %>% mutate(july=ifelse(month==7,1,0))

df.post <- rbind(
  data.frame(trips.s=rnorm(5,means2[1],5),year=c(2005:2009),month=1),
  data.frame(trips.s= rnorm(5,means2[2],5),year=c(2005:2009),month=2),
  data.frame(trips.s=rnorm(5,means2[3],5),year=c(2005:2009),month=3),
  data.frame(trips.s=rnorm(5,means2[4],5),year=c(2005:2009),month=4),
  data.frame(trips.s=rnorm(5,means2[5],5),year=c(2005:2009),month=5),
  data.frame(trips.s=rnorm(5,means2[6],5),year=c(2005:2009),month=6),
  data.frame(trips.s=rnorm(5,means2[7],5),year=c(2005:2009),month=7),
  data.frame(trips.s=rnorm(5,means2[8],5),year=c(2005:2009),month=8),
  data.frame(trips.s=rnorm(5,means2[9],5),year=c(2005:2009),month=9),
  data.frame(trips.s=rnorm(5,means2[10],5),year=c(2005:2009),month=10),
  data.frame(trips.s=rnorm(5,means2[11],5),year=c(2005:2009),month=11),
  data.frame(trips.s=rnorm(5,means2[12],5),year=c(2005:2009),month=12)
)

df.post <- tbl_df(df.post) %>% mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month) %>% mutate(july=ifelse(month==7,1,0))

#----------------------------------------------------------------
# slightly different simulation
error <- rnorm(120,0,5)
chow2.df <- tbl_df(data.frame(year=c(rep(1995,12),rep(1996,12),rep(1997,12),rep(1998,12),
                                     rep(1999,12),rep(2000,12),rep(2001,12),
                                     rep(2002,12),rep(2003,12),rep(2004,12)),
                              month=rep(1:12,10),
                              int=11.206,error=error)) %>%
  mutate(
    feb=ifelse(month==2,7.8,0),
    march=ifelse(month==3,17.47,0),
    april=ifelse(month==4,30.5,0),
    may=ifelse(month==5,38.2,0),
    june=ifelse(month==6,49.4,0),
    july=ifelse(month==7,59.3,0),
    aug=ifelse(month==8,26.2,0),
    sept=ifelse(month==9,19.2,0),
    oct=ifelse(month==10,11.8,0),
    nov=ifelse(month==11,7.4,0),
    dec=ifelse(month==12,0.8,0),
    trips.s=int+feb+march+april+may+june+july+aug+sept+oct+nov+dec+error,
    date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month)
error <- rnorm(72,0,5)              
chow2.df2 <- tbl_df(data.frame(year=c(rep(2005,12),rep(2006,12),rep(2007,12),rep(2008,12),
                                      rep(2009,12),rep(2010,12)),
                               month=rep(1:12,6),
                               int=69.6,error=error)) %>%
  mutate(
    feb=ifelse(month==2,-10.6,0),
    march=ifelse(month==3,-16.5,0),
    april=ifelse(month==4,-27.4,0),
    may=ifelse(month==5,-37.8,0),
    june=ifelse(month==6,-50.5,0),
    july=ifelse(month==7,-61.6,0),
    aug=ifelse(month==8,-50.35,0),
    sept=ifelse(month==9,-37.11,0),
    oct=ifelse(month==10,-26.5,0),
    nov=ifelse(month==11,-8.5,0),
    dec=ifelse(month==12,-3.1,0),
    trips.s=int+feb+march+april+may+june+july+aug+sept+oct+nov+dec+error,
    date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month)

chow2.df <- chow2.df %>% select(trips.s,year,month,date) %>% mutate(july=ifelse(month==7,1,0))
chow2.df2 <- chow2.df2 %>% select(trips.s,year,month,date) %>% mutate(july=ifelse(month==7,1,0))

#---------------------------------------------------------------



#combine the two data frames
ss.data <- rbind(df,df.post)
#ss.data <- rbind(chow2.df,chow2.df2)


#give the data a trend
trend <- 2
c <- c(200,rep(0,length(ss.data$trips.s)-1))
for(i in 2:length(c)){
  c[i]<- c[i-1] + trend + rnorm(1,0,1)
}

ss.data$c <- c
ss.data$trips <- ss.data$c + ss.data$trips.s
ss.data$t <- seq(1:nrow(ss.data))

ggplot(ss.data,aes(x=date,y=trips.s))+geom_line()+geom_point(aes(color=factor(july))) + 
  scale_color_manual(values=c('red','black')) + theme_bw() + 
  theme(legend.position="none")

ggplot(ss.data,aes(x=date,y=trips))+geom_line()+geom_point(aes(color=factor(july))) + 
  scale_color_manual(values=c('red','black')) + theme_bw() + 
  theme(legend.position="none")

#--------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#specify a state-space model with seasonals only
dates <- ss.data$date
trips <- xts(ss.data$trips.s,order.by=dates)
trips.trend <- xts(ss.data$trips, order.by=dates)

tripsModel<-SSModel(trips ~ SSMseasonal(period=12, sea.type="dummy", Q = NA), H = NA)

#state space model with local level and seasonals
tripsModel_trend<-SSModel(trips.trend ~ SSMtrend(degree = 1, Q=list(matrix(NA))) + 
                      SSMseasonal(period=12, sea.type="dummy", Q = NA), H = NA)

#---------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#Get the likelihood maximizing parameter values for the model
tripsFit<-fitSSM(tripsModel,inits=c(0.05, 0.001),method='BFGS')$model
tripsFit_trend<-fitSSM(tripsModel_trend,inits=c(0.1,0.05, 0.001),method='BFGS')$model

#Get the smoothed estimates of the state values
tripsSmooth <- KFS(tripsFit,smooth= c('state', 'mean','disturbance'))
tripsSmooth_trend <- KFS(tripsFit_trend,smooth= c('state', 'mean','disturbance'))

#Get smoothed estimates of the level
#tripsLevel <-signal(tripsSmooth, states = 'level')
#tripsLevel$signal <- xts(tripsLevel$signal, order.by = dates)
#autoplot(cbind(trips, tripsLevel$signal),facets = FALSE)
#-----------------------------------------------------------------------------

#---------------------------------------------------------------------------
#extract the seasonal
trips.sea <- as.numeric(tripsSmooth$alphahat[,2])
plot.df <- data.frame(date=dates,seasonal=trips.sea)

trips.sea_trend<- as.numeric(tripsSmooth_trend$alphahat[,2])

#add a reference point to the seasonal data frame
plot.df <- plot.df %>% mutate(month=month(date),july=ifelse(month==7,1,0))

ggplot(plot.df,aes(x=date,y=seasonal)) + geom_line() + geom_point(aes(color=factor(july))) +
  scale_color_manual(values=c('red','black')) + 
  theme_bw() + 
  theme(legend.position="none")

#-----------------------------------------------------------------------------

#===========================================================
#Compare our state-space model with a dummy variable regression
# with structural break:

chow.df <- ss.data %>%
  mutate(jan=ifelse(month==1,1,0),
         feb=ifelse(month==2,1,0),
         march=ifelse(month==3,1,0),
         april=ifelse(month==4,1,0),
         may=ifelse(month==5,1,0),
         june=ifelse(month==6,1,0),
         july=ifelse(month==7,1,0),
         aug=ifelse(month==8,1,0),
         sept=ifelse(month==9,1,0),
         oct=ifelse(month==10,1,0),
         nov=ifelse(month==11,1,0),
         dec=ifelse(month==12,1,0))


model.pre <- lm(trips.s~feb+march+april+may+june+july+aug+
                  sept+oct+nov+dec,data=subset(chow.df,year<=2004))
model.pre.trend <- lm(trips~feb+march+april+may+june+july+aug+
                  sept+oct+nov+dec+t,data=subset(chow.df,year<=2004))

model.post <- lm(trips.s~feb+march+april+may+june+july+aug+
                   sept+oct+nov+dec,data=subset(chow.df,year>2004))

model.post.trend <- lm(trips~feb+march+april+may+june+july+aug+
                   sept+oct+nov+dec+t,data=subset(chow.df,year>2004))


yhat <- rbind(
  data.frame(date=chow.df$date[chow.df$year<=2004],
             that=predict(model.pre,newdata=chow.df[chow.df$year<=2004,])),
  data.frame(date=chow.df$date[chow.df$year>2004],
             that=predict(model.post,newdata=chow.df[chow.df$year>2004,]))
)

#quick plot of the seasonal fit from the dummy variable model
plot.tmp <- tbl_df(yhat) %>% mutate(year=year(date),month=month(date)) %>%
              mutate(july=ifelse(month==7,1,0))
ggplot(plot.tmp,aes(x=date,y=that)) + geom_line() + geom_point(aes(color=factor(july)))


yhat_trend <- rbind(
  data.frame(date=chow.df$date[chow.df$year<=2004],
             that=predict(model.pre.trend,newdata=chow.df[chow.df$year<=2004,])),
  data.frame(date=chow.df$date[chow.df$year>2004],
             that=predict(model.post.trend,newdata=chow.df[chow.df$year>2004,]))
)


yhat$y <- as.numeric(trips)
yhat$model <- 'Seasonal Dummy'

yhat_trend$y <- as.numeric(trips.trend)
yhat_trend$model <- 'Seasonal Dummy w/trend'

# the KFAS model
kfas.pred <- data.frame(date=dates,that=fitted(tripsSmooth),y=as.numeric(trips),model='KFAS')
kfas.pred_trend <- data.frame(date=dates,that=fitted(tripsSmooth_trend),
                              y=as.numeric(trips.trend),model='KFAS trend')

#combine the two
model.comp <- rbind(yhat,yhat_trend,kfas.pred,kfas.pred_trend)

model.comp$eps <- model.comp$y-model.comp$that

tbl_df(model.comp) %>% mutate(abs.error=abs(eps)) %>% group_by(model) %>%
  summarise(mean(abs.error,na.rm=T))

#we can get some insight into the really good fit of the state-space model
# despite it's failure to estimate the seasonal effects correctly:
tbl_df(data.frame(date=dates,level=tripsSmooth$alphahat[,1],
                  state=tripsSmooth$alphahat[,2],
                  pred=fitted(tripsSmooth),
                  y=trips,seasonal=ss.data$trips.s,y.level=ss.data$c)) %>% 
  mutate(level_plus_state=level+state, month=month(date)) %>% filter(month==7)



#plot models with trend
ggplot(subset(model.comp,model=='KFAS trend'),aes(x=date,y=y)) + geom_line() + 
  geom_line(aes(x=date,y=that),color='red')


#would we have found the 'right' structural break if we didn't know
# where it was

break.test <- function(yr){
sb.df <- chow.df %>% mutate(sbreak=ifelse(year==yr,1,0))

  sb <- lm(trips.s~feb+march+april+may+june+july+aug+sept+oct+nov+dec+
           feb*sbreak+march*sbreak+
           april*sbreak+may*sbreak+june*sbreak+july*sbreak+
           aug*sbreak+sept*sbreak+sept*sbreak+oct*sbreak+
           nov*sbreak+dec*sbreak,data=sb.df)
return(extractAIC(sb))  
}
















