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
  geom_point(aes(color=factor(july)))

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
  geom_line(aes(x=date,y=kfas),data=pred,color='red') 
#------------------------------------------------------------------------------

#---------------------------------------------------------------------------
#plot the seasonal
trips.sea <- as.numeric(tripsSmooth$alphahat[,2])
plot.df <- data.frame(date=dates,seasonal=trips.sea)
#add a reference point to the seasonal data frame
plot.df <- plot.df %>% mutate(month=month(date),july=ifelse(month==7,1,0))

ggplot(plot.df,aes(x=date,y=seasonal)) + geom_line() + geom_point(aes(color=factor(july)))

#----------------------------------------------------------------------------


#plot level
level <- as.numeric(tripsSmooth$alphahat[,1])
plot.df <- data.frame(date=dates,level=level,z='smoothed')
                
ggplot(plot.df,aes(x=date,y=level)) + geom_line() + 
  geom_line(data=subset(df.monthly,area=='south'),aes(x=date,y=ntrips,color="red"))


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
ggplot(model.comp,aes(x=date,y=eps,color=model)) + geom_line()


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
  data.frame(trips=rnorm(10,means[1],5),year=c(1995:2004),month=1),
  data.frame(trips= rnorm(10,means[2],5),year=c(1995:2004),month=2),
  data.frame(trips=rnorm(10,means[3],5),year=c(1995:2004),month=3),
  data.frame(trips=rnorm(10,means[4],5),year=c(1995:2004),month=4),
  data.frame(trips=rnorm(10,means[5],5),year=c(1995:2004),month=5),
  data.frame(trips=rnorm(10,means[6],5),year=c(1995:2004),month=6),
  data.frame(trips=rnorm(10,means[7],5),year=c(1995:2004),month=7),
  data.frame(trips=rnorm(10,means[8],5),year=c(1995:2004),month=8),
  data.frame(trips=rnorm(10,means[9],5),year=c(1995:2004),month=9),
  data.frame(trips=rnorm(10,means[10],5),year=c(1995:2004),month=10),
  data.frame(trips=rnorm(10,means[11],5),year=c(1995:2004),month=11),
  data.frame(trips=rnorm(10,means[12],5),year=c(1995:2004),month=12)
)

df <- tbl_df(df) %>% mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month) %>% mutate(july=ifelse(month==7,1,0))

df.post <- rbind(
  data.frame(trips=rnorm(5,means2[1],5),year=c(2005:2009),month=1),
  data.frame(trips= rnorm(5,means2[2],5),year=c(2005:2009),month=2),
  data.frame(trips=rnorm(5,means2[3],5),year=c(2005:2009),month=3),
  data.frame(trips=rnorm(5,means2[4],5),year=c(2005:2009),month=4),
  data.frame(trips=rnorm(5,means2[5],5),year=c(2005:2009),month=5),
  data.frame(trips=rnorm(5,means2[6],5),year=c(2005:2009),month=6),
  data.frame(trips=rnorm(5,means2[7],5),year=c(2005:2009),month=7),
  data.frame(trips=rnorm(5,means2[8],5),year=c(2005:2009),month=8),
  data.frame(trips=rnorm(5,means2[9],5),year=c(2005:2009),month=9),
  data.frame(trips=rnorm(5,means2[10],5),year=c(2005:2009),month=10),
  data.frame(trips=rnorm(5,means2[11],5),year=c(2005:2009),month=11),
  data.frame(trips=rnorm(5,means2[12],5),year=c(2005:2009),month=12)
)

df.post <- tbl_df(df.post) %>% mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month) %>% mutate(july=ifelse(month==7,1,0))

#combine the two data frames
ss.data <- rbind(df,df.post)

ggplot(ss.data,aes(x=date,y=trips))+geom_line()+geom_point(aes(color=factor(july)))
#--------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#specify a state-space model with seasonals
dates <- ss.data$date
trips <- xts(ss.data$trips,order.by=dates)

tripsModel<-SSModel(trips ~ SSMtrend(degree = 1, Q=list(matrix(NA))) + 
                      SSMseasonal(period=12, sea.type="dummy", Q = NA), H = NA)
#---------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#Get the likelihood maximizing parameter values for the model
tripsFit<-fitSSM(tripsModel,inits=c(0.1,0.05, 0.001),method='BFGS')$model
#Get the smoothed estimates of the state values
tripsSmooth <- KFS(tripsFit,smooth= c('state', 'mean','disturbance'))
#Get smoothed estimates of the level
tripsLevel <-signal(tripsSmooth, states = 'level')
tripsLevel$signal <- xts(tripsLevel$signal, order.by = dates)
autoplot(cbind(trips, tripsLevel$signal),facets = FALSE)
#-----------------------------------------------------------------------------

#---------------------------------------------------------------------------
#extract the seasonal
trips.sea <- as.numeric(tripsSmooth$alphahat[,2])
plot.df <- data.frame(date=dates,seasonal=trips.sea)

#add a reference point to the seasonal data frame
plot.df <- plot.df %>% mutate(month=month(date),july=ifelse(month==7,1,0))

ggplot(plot.df,aes(x=date,y=seasonal)) + geom_line() + geom_point(aes(color=factor(july)))
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

model.pre <- lm(trips~feb+march+april+may+june+july+aug+
                  sept+oct+nov+dec,data=subset(chow.df,year<=2004))
model.post <- lm(trips~feb+march+april+may+june+july+aug+
                   sept+oct+nov+dec,data=subset(chow.df,year>2004))

yhat <- rbind(
  data.frame(date=chow.df$date[chow.df$year<=2004],
             that=predict(model.pre,newdata=chow.df[chow.df$year<=2004,])),
  data.frame(date=chow.df$date[chow.df$year>2004],
             that=predict(model.post,newdata=chow.df[chow.df$year>2004,]))
)

yhat$y <- as.numeric(trips)
yhat$model <- 'Seasonal Dummy'

# the KFAS model
kfas.pred <- data.frame(date=dates,that=fitted(tripsSmooth),y=as.numeric(trips),model='KFAS')

#combine the two
model.comp <- rbind(yhat,kfas.pred)

model.comp$eps <- model.comp$y-model.comp$that

tbl_df(model.comp) %>% mutate(abs.error=abs(eps)) %>% group_by(model) %>%
  summarise(mean(abs.error,na.rm=T))

#we can get some insight into the really good fit of the state-space model
# despite it's failure to estimate the seasonal effects correctly:
tbl_df(data.frame(date=dates,level=tripsSmooth$alphahat[,1],
                  state=tripsSmooth$alphahat[,2],
                  pred=fitted(tripsSmooth),
                  y=trips)) %>% 
  mutate(level_plus_state=level+state, month=month(date)) %>% filter(month==7)
























