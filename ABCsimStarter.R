######################################################################
## Couples cohort transmission simulator
## Simulate from priors for ABC particle swarm
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
SimMinutes <- 10 ## minutes to simulate for
maxN <- 5000

startTime <- Sys.time()
ii <- 1
timeTaken <- as.numeric(difftime(Sys.time(), startTime, units='mins'))
rcohsList <- list()
while(timeTaken < SimMinutes) { ##
    temprcoh <- retroCohSim(parms = simParmSamp(1), maxN=maxN, browse=F)
    rcohsList[[ii]] <- temprcoh
    ii <- ii+1
    timeTaken <- as.numeric(difftime(Sys.time(), startTime, units='mins'))
    print(timeTaken)
}

save(rcohsList, file = paste0('results/testDir/rcohsList-',seed,'.Rdata'))

lapply(rcohsList, function(x) sbmod.to.wdat(x$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=F, giveProp=T))

testPars <- simParmSamp(1)
testPars
testPars2 <- testPars
testPars2[c('acute.sc','dur.ac','het.gen.sd')] <- c(5,3,2)
testPars2[c('bmp','bfp')] <- testPars2[c('bmp','bfp')]*10

sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
testPars2
rcohsim <- retroCohSim(parms = testPars2, maxN=500, browse=F)
sbmod.to.wdat(rcohsim$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=F, giveProp=T)
sbmod.to.wdat(rcohsim$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=T, giveProp=T)

wtab.rlp
lapply(wtab.rlp, function(x) x$p)
