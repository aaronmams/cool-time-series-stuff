# cool-time-series-stuff

This project contains scripts supporting the analysis in these three blog posts:

* http://thesamuelsoncondition.com/2016/01/23/time-series-iii-deterministicstochastic-seasonality/
* http://thesamuelsoncondition.com/2016/01/22/time-series-ii-autocorrelation-serial-correlation-error-structures/
* http://thesamuelsoncondition.com/2016/01/16/time-series-part-i-some-popular-diagnostics/

The purpose of these exercises is to allow me to expiriment with/refresh my knowledge of important time-series methods.  The particular research question motivating this exploration is: how can we identify whether a pattern of seasonality has changed as a result of some policy?

Some related subquestions:

* how do we evaluate seasonality in a time-series?
* how do we determine if the seasonal pattern is stable through time or evolving?
* how do we test for structural breaks in data when the number and location of breaks is unknown?
* if we have several structural breaks in a time-series are they really breaks? or are they the result of a dynamic process?  i.e. do we have i) a regime switching process, ii) a state-dependent process, or iii) a process with multiple discrete breaks?
* what are the consequences of choosing the wrong model? i.e. if we model a process with multiple discrete breaks as a Markov regime switching process what are the consquences of this model mis-specification?


This is a work in progress.  Please feel free to suggest improvements to my code through a pull-request.  The script **time series explorer.R** should allow you to replicate all of the analysis in the aforementioned blog posts. 

## Getting Started

The backbone of this project is the function **seasonal.data()**.  Currently, this function is not terribly flexible.  It will generate a single data frame containing two different time-series, each with 10 years of monthly data.  Both series will have an abrupt change in the seasonal pattern following year 6.  

The first series generated will have 6 years of data wherein the annual total is approximately evenly distributed among the twelve months.  In the last 4 years of data the process will transitions such that 25% of the annual total is approximately evenly distributed among the months (jan, feb, march, april, may, sept, oct, nov, dec) and 75% of the annual total is approximately evenly distributed over the three summer months (june, july, aug).  

This series allows us to apply different time-series methods to a process which abruptly changes from having no seasonality to having notably seasonality.  We will be able to see how this type of shift affects different time-series diagnostics and models.

The second series generated will have 6 years of data wherein the annual total is distributed seasonally according to the parameters in the input *shares.before*.  In the last 4 years of the series the seasonality changes to a process where the annual total is distributed according to the parameters in the input *shares.after*. 

This series will allow us to see how a changing pattern of seasonality affect different time-series diagnostics and models.

The critical feature of both data series (for my purposes, at least) is that both series have an annual total that is approximately equal for every year in the series.  

##Using the seasonal.data() function

The function **seasonal.data()** accepts the following inputs:

* mu - this is the annual total for the first series.  this is the amount that will be divided up into monthly shares in each year of the simulated data.
* sigma1 - this is a 1X3 vector of variance parameters for the first time-series.  
  * sigma1[1] - this determines the noise that will be generated around the monthly means in the first 6 years of the first series
  * sigma1[2] - this determines the noise that will be generated around the monthly means in the second 4 years of the data, for the 9 months where 25% of the annual total is apportioned approximately equally.
  * sigma1[3] - this determines the noise that will be generated around the monthly means in the second 4 years of the data, for the 3 months where 75% of the annual total is apportioned approximately equally
* annual.total - this is the annual total for the second series
* shares.before - this is a 1X12 vector that determines the monthly shares of the total in the first 6 years of the series.  Example: shares.before=c(0.01,0.01,0.05,0.05,0.01,0.01,0.2,0.2,0.1,0.1,0.2,0.06) says that on average january gets 10% of *annual.total*, february gets 10%, march gets 5%, april gets 5%,....., and december gets 6% of *annual.total*.
* shares.after - this is a 1X12 vector that determines the monthly shares of the total in the last 4 years of the series.
* sigma2 - a numeric constant.  This is a number that determines the noise around the monthly means for the second series.  In the second series the variance for the monthly distributions is proportional to the mean.  The constant sigma2 defines what proportion of the mean the variance will be.  Larger values will produce less noise around the monthly mean values and smaller values will produce more noise.

##Estimations

The script **time series models.R** has code do a few different maximum likelihood-type estimation of Markov Regime Switching
Models.  I'm working on beefing this up but right now it just has a simple Hamilton Filter for a 2-state model of natural gas prices and oil prices that I wrote...plus a few snippets of code that will give you the flavor of how to use the **MSWM Package** to estimate a Markov Regime Switching Model.