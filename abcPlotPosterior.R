####################################################################################################
## Plot posterior results of ABC for manuscript
####################################################################################################
## Collect batch of ABC results & prepare for next batch & do some diagnostics.
rm(list=ls()); gc()
setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
ncores <- 12
finalbatch <- 4
parsDo <- c('EHMacute','het.gen.sd','acute.sc','dur.ac','bp')
logDo <- c(F,F,T,T,T)
nmsDo <- expression(EHM[acute], sigma[lambda], RH[acute], d[acute], lambda)
npars <- length(parsDo)
cols <- rainbow(npars)
distrNms <- c('Prior', paste0('Intermediate ', 1:(finalbatch-1)), 'Posterior')
xlims <- list(EHMacute=c(0,100), het.gen.sd=c(0,3), acute.sc=c(.5,200), dur.ac=c(.5,10), 
              bp = exp(c(logGHazLoBound['bp'], logGHazHiBound['bp'])))

pmc0 <- addEHM(simParmSamp(4*10^4)) ## sample prior distribution
pmc0 <- within(pmc0, {for(ll in logParms) assign(paste0('log',ll), log(get(ll))); rm(ll)}) ## add log variables to priors
for(bb in 1:finalbatch) { ## load all batches as pmc1, ..., pmcN
    load(file=file.path(out.dir, paste0('IntermedDistr',bb,'.Rdata'))) ## Load last distribution (already filtered)
    assign(paste0('pmc',bb), pmatChosen) ## pmc = parameter matrix chosen from ABC, number indicates ABC batch (intermediate distribution)
}

rgs <- apply(pmc0[,parsDo], 2, range)

## Show approach from prior to posterior through intermediate distributions for EHMacute, het.gen.sd,
pdf(file.path(fig.dir, 'posterior approach.pdf'), w = 6.83, h=8)
par(mfrow=c(3,2), lwd=2, bty='n')
for(pp in 1:npars) for(bb in 0:finalbatch) {
    if(bb==0) {
        plot(get(paste0('Dens',parsDo[pp],bb)), xlab=nmsDo[pp], log=ifelse(logDo[pp],'x',''), ylab = 'density', 
             col=cols[bb+1], ylim = c(0,get(paste0('Max',parsDo[pp]))), las = 1, main = '', xlim=xlims[[pp]])
    }else{
        lines(get(paste0('Dens',parsDo[pp],bb)), col = cols[bb+1])
    }
    if(pp==1)legend('topright', distrNms, col = cols, lwd = 2, bty = 'n')
}
graphics.off()


for(bb in 0:finalbatch) for(ii in 1:npars)  {
    pp <- parsDo[ii]
    temp <- with(get(paste0('pmc',bb)), density(get(pp), from = rgs[1,pp], to = rgs[2,pp]))
    if(logDo[ii]) temphist <- with(get(paste0('pmc',bb)), hist(log(get(pp)), plot=F, breaks = 20))
    if(!logDo[ii]) temphist <- with(get(paste0('pmc',bb)), hist(get(pp), plot=F, breaks = 20))
    assign(paste0('Dens',pp,bb), temp)
    assign(paste0('Hist',pp,bb), temphist)
    if(bb==0) assign(paste0('Max',pp), max(temp$y)) else assign(paste0('Max',pp), max(temp$y, get(paste0('Max',pp))))
}

## Show approach from prior to posterior through intermediate distributions for EHMacute, het.gen.sd,
pdf(file.path(fig.dir, 'posterior approach hist.pdf'), w = 6.83, h=8)
par(mfrow=c(3,2), lwd=2, bty='n')
for(pp in 1:npars) {
    for(bb in 0:finalbatch) {
        temphist <- get(paste0('Hist',parsDo[pp],bb))
        if(bb==0) {
            if(logDo[pp]) xs <- exp(temphist$mids) else xs <- temphist$mids
            plot(xs, temphist$density, xlab=nmsDo[pp], log=ifelse(logDo[pp],'x',''), ylab = 'density', 
                 col=cols[bb+1], ylim = c(0,get(paste0('Max',parsDo[pp]))), las = 1, main = '', xlim=xlims[[pp]])
        }else{
            if(logDo[pp]) xs <- exp(temphist$mids) else xs <- temphist$mids
            lines(xs, temphist$density, col = cols[bb+1])
        }
        if(pp==1)legend('topright', distrNms, col = cols, lwd = 2, bty = 'n')
    }}
graphics.off()


## Show approach from prior to posterior through intermediate distributions for EHMacute, het.gen.sd,
pdf(file.path(fig.dir, 'posterior approach hist.pdf'), w = 6.83, h=8)
par(mfrow=c(3,2), lwd=2, bty='n')
for(pp in 1:npars) {#for(bb in 0:finalbatch) {
hist(
}
graphics.off()

pdf(file.path(fig.dir, 'posterior approach hist.pdf'), w = 6.83, h=8)
        with(Histhet.gen.sd0, plot(mids, density, xlab=nmsDo[pp], log=ifelse(logDo[pp],'x',''), ylab = 'density', 
                            col=cols[bb], ylim = c(0,get(paste0('Max',parsDo[pp]))), las = 1, main = '', xlim=xlims[[pp]]))
#hist(pmc0$het.gen.sd)
graphics.off()



## Look at CIs nonlog
apply(pmatChosen[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
apply(priorParms[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
with(pmatChosen, cor.test(EHMacute, het.gen.sd))
