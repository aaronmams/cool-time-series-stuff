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