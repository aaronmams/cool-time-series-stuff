######################################################################
######################################################################
#######################################################################
#Estimations

#this script contains the code we need to run the Hamilton Filter and 
# Hamilton Smoother for estimating parameters of a Markov Regime Switching
# Model.  It also contains some code to compare the Hamilton Fitler with
# result from the R package 'MSWM'


require(MSwM)
require(ggplot2)
require(dplyr)
######################################################################
######################################################################
#######################################################################
######################################################################
######################################################################
#######################################################################
######################################################################
######################################################################
#######################################################################
#First order of business is to test a Hamilton Filtering algorithm on a 
# simple problem.  I will follow the lead of Matt Brigada here, 
# https://github.com/Matt-Brigida/R-Finance-2015-MSCP
# and use a model of natural gas prices and oil prices with a 2-regime
# model.

#---------------------------------------------------------------------
#use this function to get the same data that Matt Brigida uses
require(XML)
require(xts)
library(EIAdata)

getMonEIA <- function(ID, key) {
  
  ID <- unlist(strsplit(ID, ";"))
  key <- unlist(strsplit(key, ";"))
  
  url <- paste("http://api.eia.gov/series?series_id=", ID, "&api_key=", key, 
               "&out=xml", sep = "")
  
  doc <- xmlParse(file = url, isURL = TRUE)
  
  df <- xmlToDataFrame(nodes = getNodeSet(doc, "//data/row"))
  
  df <- plyr::arrange(df, df$date)
  
  date <- as.Date(paste(as.character(levels(df[, 1]))[df[, 1]], "01", sep = ""), 
                  "%Y%m%d")
  values <- as.numeric(levels(df[, -1]))[df[, -1]]
  
  tmp <- data.frame(date=date,x=values)
  names(tmp) <- c('date',paste(ID))
  return(tmp)
  #xts_data <- xts(values, order.by = date)
  #names(xts_data) <- sapply(strsplit(ID, "-"), paste, collapse = ".")
  
  #assign(sapply(strsplit(ID, "-"), paste, collapse = "."), xts_data, envir = .GlobalEnv)
}
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#in order to pull the data we will need to provide a key.  The EIA will send you a key for
# accessing their open data API....details can be found here:
#http://complete-markets.com/2014/01/r-functions-to-interact-with-the-eias-application-programming-interface-api/

key <- '80DAEF6FDCF88254A3D4E2F5F4CBF99C'
#-----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
# use the function above to get monthly data from the EIA on natural gas prices and West
# Texas Crude Oil...then filter the series so they are the same length

lng <- getMonEIA(ID='NG.RNGWHHD.M',key=key)
names(lng) <- c('date','lng')
oil <- getMonEIA(ID='PET.RWTC.M',key=key)
names(oil) <- c('date','oil')

start <- max(c(min(lng$date),min(oil$date)))
end <- min(max(lng$date),max(oil$date))

lnoil <- log(as.vector(oil$oil[which(oil$date>=start & oil$date<=end)]))
lnng <- log(as.vector(lng$lng[which(lng$date>=start & lng$date<=end)]))
#-------------------------------------------------------------------------------------------

################################################################################
################################################################################
################################################################################
# The Hamilton Filter

#need to write a function to do the calculation then optimize that function

#NOTE: this is not all that general.  It is set up to deal with a simple univariate
# regression model where:

# y(t) = alpha1 + beta1*x(t) if we are in regime 1
# y(t) = alpha2 + beta2*x(t) if we are in regime 2

# we need to generalize this later to be amenable to different regression forms

mrs.est <- function(theta,x,y){
  alpha1 <- theta[1]
  alpha2 <- theta[2]
  alpha3 <- theta[3]
  alpha4 <- theta[4]
  alpha5 <- theta[5]
  alpha6 <- theta[6]
  
  p11 <- 1/(1+exp(-theta[7]))
  p22 <- 1/(1+exp(-theta[8]))
  
  
  #in order to make inference about what state we are in in period t we need the conditional
  # densities given the information set through t-1
  f1 <- (1/(alpha5*sqrt(2*pi)))*exp(-((y-alpha1-(alpha3*x))^2)/(2*(alpha5^2)))
  f2 <- (1/(alpha6*sqrt(2*pi)))*exp(-((y-alpha2-(alpha4*x))^2)/(2*(alpha6^2)))
  f <- matrix(cbind(f1,f2),nc=2)
  
  #S.forecast is the state value looking forward conditional on info up to time t
  #S.inf is the updated state value
  S.forecast <- rep(0,2*length(y))
  S.forecast <- matrix(S.forecast,nrow=(length(y)),ncol=2)
  
  S.inf <- S.forecast
  o.v <- c(1,1)
  
  P<-matrix(c(p11,(1-p11),(1-p22),p22),nr=2,nc=2)
  model.lik <- rep(0,length(y))
  
  S.inf[1,] <- (c(p11,p22)*f[1,])/(o.v %*% (c(p11,p22)*f[1,]))
#  Matt Brigada's starting point
  S.inf[1,] <- (c(p11,p22)*f[1,])/(o.v %*% (c(p11,p22)*f[1,]))
  
  for(t in 1:(length(y)-1)){
    #in time t we first make our forecast of the state in t+1 based on the 
    # data up to time t, then we update that forecast based on the data
    # available in t+1
    S.forecast[t+1,] <- P%*%S.inf[t,]    
    S.inf[t+1,] <- (S.forecast[t+1,]*f[t+1,])/(S.forecast[t+1,] %*% f[t+1,])
    model.lik[t+1] <- o.v%*%(S.forecast[t+1,]*f[t+1,])
  }
  
  logl <- sum(log(model.lik[2:length(model.lik)]))
  return(-logl)
}
################################################################################
################################################################################
################################################################################




################################################################################
################################################################################
################################################################################
#The Hamilton Smoother:  the smoothed probabilities are obtained by:

#1.  using the Hamilton Filter to get maximum likelihood parameter estimates
#2.  using the maximum likelihood parameter estimates to run the Filter again...
#2a  but in the smoothing step we work from T-1 back to t=1


#==============================================================
#==============================================================
#to get the smoothed probabilities we run the filter again using
# the maximum likelihood values but we need to save both the  
# forecasted state probabilities (St given info in t-1) and 
# the updated state probabilities (St updated to reflect info in time t)
#=================================================================
#=================================================================
#get the maximum likelihood estimates

ham.smooth<-function(theta,y,x){
  alpha1 <- theta[1]
  alpha2 <- theta[2]
  alpha3 <- theta[3]
  alpha4 <- theta[4]
  alpha5 <- theta[5]
  alpha6 <- theta[6]
  
  #  p11 <- 1/(1+exp(-theta[7]))
  #  p22 <- 1/(1+exp(-theta[8]))
  
  p11 <- theta[7]
  p22 <- theta[8]
  
  #in order to make inference about what state we are in in period t we need the conditional
  # densities given the information set through t-1
  f1 <- (1/(alpha5*sqrt(2*pi)))*exp(-((y-alpha1-(alpha3*x))^2)/(2*(alpha5^2)))
  f2 <- (1/(alpha6*sqrt(2*pi)))*exp(-((y-alpha2-(alpha4*x))^2)/(2*(alpha6^2)))
  f <- matrix(cbind(f1,f2),nc=2)
  
  #S.forecast is the state value looking forward conditional on info up to time t
  #S.inf is the updated state value
  S.forecast <- rep(0,2*length(y))
  S.forecast <- matrix(S.forecast,nrow=(length(y)),ncol=2)
  
  S.inf <- S.forecast
  o.v <- c(1,1)
  
  P<-matrix(c(p11,(1-p11),(1-p22),p22),nr=2,nc=2)
  model.lik <- rep(0,length(y))
  
  S.inf[1,] <- (c(p11,p22)*f[1,])/(o.v %*% (c(p11,p22)*f[1,]))
  
  for(t in 1:(length(y)-1)){
    #in time t we first make our forecast of the state in t+1 based on the 
    # data up to time t, then we update that forecast based on the data
    # available in t+1
    S.forecast[t+1,] <- P%*%S.inf[t,]    
    S.inf[t+1,] <- (S.forecast[t+1,]*f[t+1,])/(S.forecast[t+1,] %*% f[t+1,])
    model.lik[t+1] <- o.v%*%(S.forecast[t+1,]*f[t+1,])
  }
  
  
  #the smoother works kind of like running the filter in reverse
  # we start with the last value of S from the filter recursion...
  # this is the value of S in the last period, updated to reflect
  # information available in that period...we then work backwards from 
  # there
  T<- length(y)
  P.smooth <- data.frame(s1=rep(0,T),s2=rep(0,T))
  P.smooth[T,] <- S.inf[T,]
  for(is in (T-1):1){
    
    #for clarify we can think of this problem has having 4 states...
    #note we word these as current period and next period because
    # the smoother works from T back to 1
    
    #1. probability that we observe S(t)=1, S(t+1)=1 
    #2. probability that we observe S(t)=1, S(t+1)=2
    #3. probability that we observe S(t)=2, S(t+1)=1
    #4. probability that we observe S(t)=2, S(t+1)=2
    
    #for #1 P[S(t)=1,S(t+1)=1|I(t+1)] = {P[S(t+1)=1|I(t+1)]*P[S(t)=1|I(t)]*P[S(t+1)=1|S(t)=1]}/P[S(t+1)=1|I(t)]
    p1 <- (S.inf[is+1,1]*S.inf[is,1]*p11)/S.forecast[is+1,1]
    
    #for #2 we have
    p2 <- (S.inf[is+1,2]*S.inf[is,1]*(1-p11))/S.forecast[is+1,2]
    
    #for #3 we have
    p3 <- (S.inf[is+1,1]*S.inf[is,2]*(1-p22))/S.forecast[is+1,1]
    
    #for #4 we have
    p4 <- (S.inf[is+1,2]*S.inf[is,2]*p22)/S.forecast[is+1,2]
    
    P.smooth[is,1] <- p1 + p2
    P.smooth[is,2] <- p3 + p4
    
  }
return(P.smooth)  
}
#==========================================================

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#########################################################################################
#########################################################################################
#Use the Hamilton Filter to get maximum likelihood parameters of the following model:

# y(t) = alpha1 + beta1*x(t) if regime = 1
# y(t) = alpha2 + beta2*x(t) if regime = 2

#where:
# y=price of natural gas
# x=price of oil

theta.start <- c(0.1,-0.01,0.4,0.3,0.1,0.2,0.5,0.5)
opt.mams <- optim(theta.start,mrs.est,x=lnoil,y=lnng,hessian=T,control=list(maxit=20000))

theta <- opt.mams$par
#########################################################################################
#########################################################################################


#-------------------------------------------------------------------
#first I'm going to QA my smoother by comparing smoothed state
# probabilities from my approach versus the MSwM package.  To do
# this I need to run the filter using the same estimated parameters
# that the MSwM gets
mod1 <- lm(lnng~lnoil)
mod1.mswm <- msmFit(mod1,k=2,p=0,sw=rep(TRUE,3),control=list(parallel=F))
summary(mod1.mswm)
dev.off()
plotProb(mod1.mswm,which=2)

sigma.mswm <- mod1.mswm@std
theta.smooth <- c(mod1.mswm@Coef[1,1],mod1.mswm@Coef[2,1],mod1.mswm@Coef[1,2],mod1.mswm@Coef[2,2],sigma.mswm,
           mod1.mswm@transMat[1,1],mod1.mswm@transMat[2,2])

tmp <- ham.smooth(theta=theta.smooth,y=lnng,x=lnoil)
tmp$t <- seq(1:nrow(tmp))
tmp$model="Mams"

mswm.smooth <- data.frame(mod1.mswm@Fit@smoProb[2:nrow(mod1.mswm@Fit@smoProb),])
names(mswm.smooth) <- c('s1','s2')
mswm.smooth$t <- seq(1:nrow(mswm.smooth))
mswm.smooth$model <- 'MSwM'

tmp <- rbind(tmp,mswm.smooth)
  ggplot(tmp,aes(x=t,y=s2,color=model)) + geom_line() 
#---------------------------------------------------------------------
###########################################################################################
###########################################################################################
###########################################################################################

###########################################################################################
###########################################################################################
###########################################################################################
#Now let's reproduce the plots that the MSwM package makes by adding 'recession bars' to
# indicate which regime my Hamilton Filter/Hamilton Smoother says we are in...



###########################################################################################
###########################################################################################
###########################################################################################



###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
#Apply the MSWM Package to our fake seasonal data
#==================================================================================
#==================================================================================
#==================================================================================
#As a first pass let's see if the Markov Regime Switching Model in the 
# MSwM package can identify the break in our seasonal data?

#MRSS with no autoregressive term---------------------------------
mod1 <- lm(lbs~feb+mar+apr+may+jun+jul+aug+sep+oct+nov+dec,data=z[which(z$series==1),])
mod1.mswm <- msmFit(mod1,k=2,p=0,sw=rep(TRUE,13),control=list(parallel=F))
summary(mod1.mswm)
dev.off()
plotProb(mod1.mswm,which=2)
#-------------------------------------------------------------------

#MRSS with AR(1) term---------------------------------
mod2 <- lm(lbs~feb+mar+apr+may+jun+jul+aug+sep+oct+nov+dec,data=z[which(z$series==1),])
mod2.mswm <- msmFit(mod2,k=2,p=1,sw=rep(TRUE,14),control=list(parallel=F))
summary(mod2.mswm)
dev.off()
plotProb(mod2.mswm,which=2)
#-------------------------------------------------------------------

#----------------------------------------------------------------
#MRSS model estimated in first differences
mod.dff <- lm(diff~feb+mar+apr+may+jun+jul+aug+sep+oct+nov+dec,data=z[which(z$series==1),])
moddiff.mswm <- msmFit(mod.dff,k=2,p=1,sw=rep(TRUE,14),control=list(parallel=F))
summary(moddiff.mswm)
dev.off()
plotProb(moddiff.mswm,which=2)

#----------------------------------------------------------------


#MRSS with just an intercept and AR(1) term--------------------

mod3 <- lm(lbs~1,data=z)
mod3.mswm <- msmFit(mod3,k=2,p=1,sw=rep(TRUE,3),control=list(parallel=F))
dev.off()
plotProb(mod3.mswm,which=2)
summary(mod3.mswm)
#--------------------------------------------------------------------

