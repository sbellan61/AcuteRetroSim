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
SimMinutes <- .5 ## minutes to simulate for
maxN <- 10^4


startTime <- Sys.time()

simParmSamp(1)

abcSimSumStat(maxN=2000, browse=F)

ii <- 1
timeTaken <- as.numeric(difftime(Sys.time(), startTime, units='mins'))
while(timeTaken < SimMinutes) { ## 
    temp <- with(parms[ii,], psrun(maxN = 1000, jobnum = seed, pars = parms[ii, hazs], save.new = T, ret
                                   acute.sc = acute.sc, dur.ac = dur.ac, het.gen.sd = het.gen.sd, browse = F, nc = 12))
    ii <- ii+1
    timeTaken <- as.numeric(difftime(Sys.time(), startTime, units='mins'))
    print(timeTaken)
}

rm(list=ls(all=T)) ## clear workspace


