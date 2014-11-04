## Continuous time model

######################################################################
## Couples cohort transmission simulator
###################################################################### 
## Steve Bellan, 2014, steve.bellan@gmail.com
if(Sys.info()['nodename']=='stevebemacbook3') setwd('~/Documents/R Repos/AcuteRetroSim/') else setwd('Rakai/AcuteRetroSim/')
rm(list=ls(all=T))                           # clear workspace
args=(commandArgs(TRUE))                # load arguments from R CMD BATCH 
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
      }  }else{ seed <- 1}
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
set.seed(seed)
SimTime <- 60*30 ## seconds to simulate for
maxN <- 10^4

parm <- simParmSamp(1)

rtime <- rexp(1)
rhaz_bmb <- rtime/parm$bmb



We really need Google vim/emacs.

JD

----------------------------------------------------------------------

We will choose time to events by:

* picking a unitless interval

* dividing by β to get a prevalence-time interval, p

* rescaling time using our spline functions to get the actual waiting time.

We will deal with the "inverse-approximation glitch" by calculating
time until the next event using:

finv(f(t)+p) - finv(f(t)),

where t is current time, f is the smoothed cumulative prevalence-time
function, and finv is its smoothed sort-of inverse. This ensures that
in the limit p→0, the waiting time goes to zero as well.

----------------------------------------------------------------------

The flow for each sub-couple is:

* For each individual:
** Choose a pre-couple infection time
** discard infection if it is after their coupling time, ELSE
*** discard whole couple if death time is before survey

* IF everybody is infected already, STOP

* For each individual not infected yet
** Choose an extra-couple infection time
** discard infection if it is after their coupling time
** Don't recalculate death times yet (one may change)

* IF nobody is infected yet, STOP

* For first individual infected
** pick death time if not picked yet; discard whole couple if death
time is before survey
** Otherwise, pick a within-couple infection time and update partner's
infection time if this one is earlier

* For other individual
** Pick death time and discard if necessary

----------------------------------------------------------------------

Still need:

* Latin hypercube notes

* Weibull non-picking
