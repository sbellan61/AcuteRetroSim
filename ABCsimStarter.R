######################################################################
## Couples cohort transmission simulator
## Simulate from priors for ABC particle swarm
###################################################################### 
## Steve Bellan, 2014, steve.bellan@gmail.com
rm(list=ls(all=T))                           # clear workspace
args=(commandArgs(TRUE))                # load arguments from R CMD BATCH 
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }  }else{ 
        seed <- 1; batch.dir <- file.path('results','testDir')
        batch <- 1
        ncores <- 12
        if(Sys.info()['nodename']=='stevebemacbook3') setwd('~/Documents/R Repos/AcuteRetroSim/') else setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
    }
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
set.seed(seed)
SimMinutes <- 12*60 ## minutes to simulate for
maxN <- 500#4785 ## Rakai couple population (from Porter et al. 2004), seems like this is size from which retrospective cohort was identified

if(batch==1) {
    simParms <- simParmSamp(SimMinutes*20) ## each sim takes a minute, so this is conservatively large number of samples
}else{
    load(file.path(out.dir, paste0('IntermedDistr',batch-1,'.Rdata'))) ## load last batch
    lastParticleSamp <- pmatChosen[sample(1:nrow(pmatChosen), SimMinutes*100, prob = pmatChosen$weight, replace=T),parnms] ## sample according to weights
    simParms <- perturbParticle(lastParticleSamp, sds = sds, browse=F) ## perturb particles
    if(nrow(simParms) < SimMinutes*20) {
        print('too few particles after restricting by prior')
    }else{
        simParms <- simParms[1:(SimMinutes*20),] ## reduce size to avoid having big stuff in memory
        rm(pmatChosen,sds,lastParticleSamp); gc()
    }
}

print(paste('seed is', seed))

startTime <- Sys.time()
ii <- 1
timeTaken <- as.numeric(difftime(Sys.time(), startTime, units='mins'))
rcohsList <- list()
while(timeTaken < SimMinutes) { ##
    temprcoh <- retroCohSim(parms = simParms[ii,], seed = seed, maxN=maxN, browse=F, nc = ncores)
    plausible <- with(temprcoh$rakll, sum(phase=='inc') > 0 & sum(inf[phase=='inc']) > 0 & 
                          sum(phase=='prev') > 0 & sum(inf[phase=='prev']) > 0) ## otherwise causes errors later on
    if(plausible) {
        rcohsList[[ii]] <- temprcoh
        ii <- ii+1
    }
    timeTaken <- as.numeric(difftime(Sys.time(), startTime, units='mins'))
    print(paste('on', ii,'taken', round(timeTaken,2),'mins'))
    if(ii%%30==0) save(rcohsList, file = file.path(batch.dir, paste0('rcohsList-',seed,'.Rdata'))) ## save every 30
}

save(rcohsList, file = file.path(batch.dir, paste0('rcohsList-',seed,'.Rdata')))

## head(temprcoh$rakll)
## temprcoh$rakll <- within(temprcoh$rakll, {inf[phase=='inc'] <- 0})
## wt <- sbmod.to.wdat(temprcoh$rakll, simpPois=T, browse=F, condRakai=T, giveLate=F)
## wt
## ct <- contTabFxn(wt)
## gS <- gSumStat(ct)

## duplicated(do.call(rbind.data.frame, lapply(rcohsList, '[[','pars')))
## lapply(rcohsList, function(x) sbmod.to.wdat(x$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=T, giveProp=T))
## sbmod.to.wdat(rcohsList[[1]]$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=F, giveProp=T)
## sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
## temprcoh <- retroCohSim(parms = simParmSamp(1), seed = seed, maxN=1200, browse=F)
## sbmod.to.wdat(temprcoh$rakll, excl.by.err = T, browse=F, simpPois=T, giveLate=F, condRakai=T, giveProp=T)
