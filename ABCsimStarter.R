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

parms <- simParmSamp(1000)
startTime <- Sys.time()

ii <- 1
while(as.numeric(Sys.time() - startTime) < SimTime) { ## 
    temp <- with(parms[ii,], psrun(maxN = 400, jobnum = seed, pars = parms[ii, hazs], 
                                   acute.sc = acute.sc, dur.ac = dur.ac, het.gen.sd = het.gen.sd, browse = F, nc = 1))
    ii <- ii+1
}

load(temp)                              # load output of simulations
names(output)                           # list objects summarized below
## jobnum=which job; evout=line list of each couple; tss=time series of pseudopopulation, rakacRR=
## retrospective cohort estimate of acute infectivity; pars=input parameters; tmar=marriage cohort
## dates; each=# couples per marriage cohort
## 
head(output$evout)                      # columns explained below
rm(list=ls(all=T))                           # clear workspace
