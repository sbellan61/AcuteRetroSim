####################################################################################################
## Plot posterior results of ABC for manuscript
####################################################################################################
## Collect batch of ABC results & prepare for next batch & do some diagnostics.
rm(list=ls()); gc()
setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
ncores <- 12
finalbatch <- 4

pmc0 <- addEHM(simParmSamp(4*10^4)) ## sample prior distribution
pmc0 <- within(pmc0, {for(ll in logParms) assign(paste0('log',ll), log(get(ll))); rm(ll)}) ## add log variables to priors
for(bb in 1:finalbatch) { ## load all batches as pmc1, ..., pmcN
    load(file=file.path(out.dir, paste0('IntermedDistr',bb,'.Rdata'))) ## Load last distribution (already filtered)
    assign(paste0('pmc',bb), pmatChosen) ## pmc = parameter matrix chosen from ABC, number indicates ABC batch (intermediate distribution)
}


## Show approach from prior to posterior through intermediate distributions for EHMacute, het.gen.sd,
pdf(file.path(fig.dir, 'posterior approach.pdf'), w = 6.83, h=4)
par(mfrow=c(4,4))
plot(
graphics.off()





## Look at CIs nonlog
apply(pmatChosen[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
apply(priorParms[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
with(pmatChosen, cor.test(EHMacute, het.gen.sd))
