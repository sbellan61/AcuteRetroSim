####################################################################################################
## Plot posterior results of ABC for manuscript
####################################################################################################
## Collect batch of ABC results & prepare for next batch & do some diagnostics.
rm(list=ls()); gc()
setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
ncores <- 12
finalbatch <- 4

transfEHM <- TRUE
parsDo <- c('acute.sc','dur.ac','EHMacute','bp','het.gen.sd')
logDo <- c(T,T,F,T,F) ## which variables to show on log scale
flatPriors <- c(T,T,F,T,T) ## which variables have flat priors
nmsDo <- expression( RH[acute], d[acute],EHM[acute],lambda, sigma[lambda] ) ## variable names
npars <- length(parsDo)
cols <- rainbow(npars)
distrNms <- c('Prior', paste0('Intermediate ', 1:(finalbatch-1)), 'Posterior') ## distribution names
xlims <- list(acute.sc=c(.5,200), dur.ac=c(.5,10), EHMacute=c(-5,100,500), 
              bp = exp(c(logGHazLoBound['bp'], logGHazHiBound['bp'])),het.gen.sd=c(0,3))
xticks <- list(acute.sc=c(.2,1,5,10,50,100,200), dur.ac=c(.5,1:10), EHMacute=c(-5,1,2,5,10,20,50,100,200,500), 
               bp = signif(exp(seq(logGHazLoBound['bp'], logGHazHiBound['bp'], l = 5)),2), het.gen.sd=seq(0,3,by=.5))
flatDens <- list(rhAcutePrior(parms = 5), durAcutePrior(parms=3),EHMacute=NA, 
                 bp = logGHazPrior(parms=logGHazPrior(1))[3], het.gen.sd = hetGenSDPrior(parms=2))
logxticks <- xticks
for(ii in which(logDo)) logxticks[[ii]] <- log(xticks[[ii]])

rgs <- list(acute.sc=c(.5,200), dur.ac=c(.5,8), EHMacute=c(-5,500), bp=exp(c(logGHazLoBound[3],  logGHazHiBound[3])), het.gen.sd=c(0,3))
for(ii in which(logDo)) rgs[[ii]] <- log(rgs[[ii]])
if(transfEHM) {
    rgs[['EHMacute']] <-  log(rgs[['EHMacute']]+10)
    ehmticks <- log(xticks[['EHMacute']]+10)
    xlims[['EHMacute']] <- log(xlims[['EHMacute']]+10)
}                

pmc0 <- addEHM(simParmSamp(4*10^4)) ## sample prior distribution
pmc0 <- within(pmc0, {for(ll in logParms) assign(paste0('log',ll), log(get(ll))); rm(ll)}) ## add log variables to priors
for(bb in 1:finalbatch) { ## load all batches as pmc1, ..., pmcN
    load(file=file.path(out.dir, paste0('IntermedDistr',bb,'.Rdata'))) ## Load last distribution (already filtered)
    assign(paste0('pmc',bb), pmatChosen) ## pmc = parameter matrix chosen from ABC, number indicates ABC batch (intermediate distribution)
}

for(bb in 0:finalbatch) for(ii in 1:npars)  {
    pp <- parsDo[ii]
    tempPmc <- get(paste0('pmc',bb))
    if(logDo[ii]) {
        temp <- with(tempPmc, density(log(get(pp)), from = rgs[[pp]][1], to = rgs[[pp]][2]))
    }else{
        temp <- with(tempPmc, density(get(pp), from = rgs[[pp]][1], to = rgs[[pp]][2]))
    }
    if(transfEHM & pp=='EHMacute') temp <- with(tempPmc, density(log(get(pp)+10), from = rgs[[pp]][1], to = rgs[[pp]][2]))
    assign(paste0('Dens',pp,bb), temp)
    if(bb==0) assign(paste0('Max',pp), max(temp$y)) else assign(paste0('Max',pp), max(temp$y, get(paste0('Max',pp))))
}

## Show approach from prior to posterior through intermediate distributions for EHMacute, het.gen.sd,
toDo <- 1:npars
#toDo <- (1:npars)[-3]
for(ff in 1:2) {
    if(ff==11) pdf(file.path(fig.dir, 'posterior approach EHM.pdf'), w = 6.83, h=5) else{
        png(file.path(fig.dir, 'posterior approach EHM.png'), w = 6.83, h=5, units='in',res=200)}
    par(mfrow=c(2,3), lwd=2, bty='n', mar = c(5,4,2,.5),'ps'=12)
    for(pp in toDo) for(bb in 0:finalbatch) {
        tempDens <- get(paste0('Dens',parsDo[pp],bb))
        if(bb==0) {
            plot(tempDens, xlab=nmsDo[pp], ylab = 'density', xlim = rgs[[pp]], type = ifelse(flatPriors[pp],'n','l'),
                 col=cols[bb+1], ylim = c(0,get(paste0('Max',parsDo[pp]))), las = 1, main = paste0('(',LETTERS[pp],')'), xaxt='n')
            if(flatPriors[pp]) segments(rgs[[pp]][1], flatDens[[pp]],rgs[[pp]][2], flatDens[[pp]], col = cols[bb+1])
            if(logDo[pp]) axis(1, at = logxticks[[pp]], label = xticks[[pp]]) else{
                if(transfEHM & parsDo[pp]=='EHMacute') axis(1, at = ehmticks, label=xticks[[pp]]) else axis(1, at = xticks[[pp]]) 
            }
        }else{
            lines(tempDens, col = cols[bb+1])
        }
    }
    par(mar=rep(2,4))
    plot(0,0, type = 'n', bty ='n', axes=F)
    legend('topright', distrNms, col = cols, lwd = 2, bty = 'n')
}
graphics.off()


## Without EHMacute
toDo <- (1:npars)[-3]
for(ff in 1:2) {
    if(ff==11) pdf(file.path(fig.dir, 'posterior approach.pdf'), w = 6.83, h=5) else{
        png(file.path(fig.dir, 'posterior approach.png'), w = 6.83, h=5, units='in',res=200)}
    par(mfrow=c(2,2), lwd=2, bty='n', mar = c(5,4,2,.5),'ps'=12)
    for(ipp in 1:length(toDo)) for(bb in 0:finalbatch) {
        pp <- toDo[ipp]
        tempDens <- get(paste0('Dens',parsDo[pp],bb))
        if(bb==0) {
            plot(tempDens, xlab=nmsDo[pp], ylab = 'density', xlim = rgs[[pp]], type = ifelse(flatPriors[pp],'n','l'),
                 col=cols[bb+1], ylim = c(0,get(paste0('Max',parsDo[pp]))), las = 1, main = paste0('(',LETTERS[ipp],')'), xaxt='n')
            if(flatPriors[pp]) segments(rgs[[pp]][1], flatDens[[pp]],rgs[[pp]][2], flatDens[[pp]], col = cols[bb+1])
            if(logDo[pp]) axis(1, at = logxticks[[pp]], label = xticks[[pp]]) else{
                if(transfEHM & parsDo[pp]=='EHMacute') axis(1, at = ehmticks, label=xticks[[pp]]) else axis(1, at = xticks[[pp]]) 
            }
        }else{
            lines(tempDens, col = cols[bb+1])
        }
    }
    if(ipp==4) legend('topleft', distrNms, col = cols, lwd = 2, bty = 'n', cex = .8)
}
graphics.off()

## CI convergence

for(bb in 0:finalbatch) {
    tempArr <- apply(get(paste0('pmc',bb))[,c('EHMacute','acute.sc','dur.ac','het.gen.sd','bp')],
                     2, function(x) quantile(x,c(.025,.5,.975)))
    if(bb==0) ciArr <- tempArr else ciArr <- abind(ciArr,tempArr, along=3)
}
ciArr[,'EHMacute',]

## Table of proportion of couples by interval from final distribution
load(file=file.path(out.dir, paste0('IntermedDistr',finalbatch,'.Rdata'))) ## Load last distribution (already filtered)
propFromWtab <- function(x)  x[[1]][,3,]  / apply(x[[1]][,3:4,],c(1,3),sum)
WpropsRl <- signif(with(wtab.rl, cbind(with(inct, i/n), with(prevt, i/n))),2)
Wprops <- lapply(gtabs, propFromWtab)
Wprops <- array(unlist(Wprops), dim=c(4,2,length(Wprops)))
WpropsMed <- signif(apply(Wprops, 1:2, function(x) quantile(x, .5, na.rm=T)),2)
WpropsLCI <- signif(apply(Wprops, 1:2, function(x) quantile(x, .025, na.rm=T)),2)
WpropsUCI <- signif(apply(Wprops, 1:2, function(x) quantile(x, .975, na.rm=T)),2)

cistrg <- function(med,lci,uci) paste0(med, ' (',lci,', ',uci,')')
WtabFit <- data.frame(WawInc=WpropsRl[,1], PostInc=cistrg(WpropsMed[,1],WpropsLCI[,1],WpropsUCI[,1]), 
                      WawInc=WpropsRl[,2], PostInc=cistrg(WpropsMed[,2],WpropsLCI[,2],WpropsUCI[,2]))
write.csv(WtabFit, file.path(fig.dir, 'Fit to proportions in Wawer table.csv'))

## Table of RHacute & dacute values for future modeling work

## Look at CIs nonlog
apply(pmatChosen[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
with(pmatChosen, cor.test(EHMacute, het.gen.sd))
