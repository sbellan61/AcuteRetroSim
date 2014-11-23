####################################################################################################
## Plot posterior results of ABC for manuscript
####################################################################################################
## Collect batch of ABC results & prepare for next batch & do some diagnostics.
rm(list=ls()); gc()
setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
ncores <- 12
finalbatch <- 5

transfEHM <- TRUE
parsDo <- c('acute.sc','dur.ac','EHMacute','bp','het.gen.sd')
logDo <- c(T,T,F,T,F) ## which variables to show on log scale
flatPriors <- c(T,T,F,T,T) ## which variables have flat priors
nmsDo <- expression( RH[acute], d[acute],EHM[acute],lambda, sigma[lambda] ) ## variable names
npars <- length(parsDo)
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

adjustBW <- c(rep(2,5))
for(bb in 0:finalbatch) for(ii in 1:npars)  {
    pp <- parsDo[ii]
    tempPmc <- get(paste0('pmc',bb))
    if(logDo[ii]) {
        temp <- with(tempPmc, density(log(get(pp)), from = rgs[[pp]][1], to = rgs[[pp]][2], adjust=adjustBW[ii]))
    }else{
        temp <- with(tempPmc, density(get(pp), from = rgs[[pp]][1], to = rgs[[pp]][2], adjust=adjustBW[ii]))
    }
    if(transfEHM & pp=='EHMacute') temp <- with(tempPmc, density(log(get(pp)+10), from = rgs[[pp]][1], to = rgs[[pp]][2], adjust=adjustBW[ii]))
    assign(paste0('Dens',pp,bb), temp)
    if(bb==0) assign(paste0('Max',pp), max(temp$y)) else assign(paste0('Max',pp), max(temp$y, get(paste0('Max',pp))))
}

## Show approach from prior to posterior through intermediate distributions for EHMacute, het.gen.sd,
toDo <- 1:npars
cols <- c('black',colorRampPalette(c('gray','yellow','red'))(finalbatch))
cols <- c('black',rainbow(finalbatch))
#toDo <- (1:npars)[-3]
for(ff in 1:2) {
    if(ff==1) pdf(file.path(fig.dir, 'posterior approach EHM.pdf'), w = 6.83, h=5) else{
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
    legend('topright', distrNms, col = cols, lwd = 2, bty = 'n', cex=.8)
}
graphics.off()


## Without EHMacute
toDo <- (1:npars)[-3]
for(ff in 1:2) {
    if(ff==1) pdf(file.path(fig.dir, 'posterior approach.pdf'), w = 6.83, h=5) else{
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
ciArr[,'EHMacute',] ## still shrinking in batch 3->4
ciArr ## still shrinking in batch 3->4
for(bb in 0:finalbatch) write.csv(ciArr[,,bb], file.path(fig.dir, paste0('ciArr',bb,'.csv')))
ciArr[,'bp',6]*12

## Get median transmission rates (i.e. mean on log scale) for paper
hazpost <- get(paste0('pmc',finalbatch))[,c('bp','het.gen.sd')]
hazpost <- within(hazpost, {
    medianbp <- bp / exp(het.gen.sd^2 / 2) ## mean logbp is median bp
    annualMedian_bp <- medianbp*12
    annualMean_bp <- bp*12
    bp025 <- exp(qnorm(.025,0,2))*annualMedian_bp
    bp975 <- exp(qnorm(.975,0,2))*annualMedian_bp
})
apply(hazpost, 2, function(x) quantile(x, c(.025,.5,.975)))

finalbatch <- 5
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
WtabFit

for(ff in 1:2) {
    if(ff==1) pdf(file.path(fig.dir, 'GoF.pdf'), w = 3.5, h=5.5) else{
              png(file.path(fig.dir, 'GoF.png'), w = 3.5, h=5.5,units='in',res=200) }
    par(mfcol=c(4,2), mar = c(2,2,2,.5), oma = c(3,3,0,0), 'ps'=10)
    for(jj in 1:2) for(ii in 1:4) {
        main <- ifelse(ii==1, c('incident couples','prevalent couples')[jj], '')
        hist(Wprops[ii,jj,], breaks = ifelse(jj==1&ii==2,5,10), xlab = '', col = 'black', main = main, yaxt='n', xlim = c(0,1))
        abline(v=WpropsRl[ii,jj], col = 'red')
        if(jj==1) mtext(paste0('interval',ii) ,2, 3)
    }
    mtext('proportion seroconverting',1, 1, outer=T)
    mtext('frequency',2, -1, outer=T)
    graphics.off()
}

## Table of RHacute & dacute values for future modeling work

## Look at CIs nonlog
apply(pmatChosen[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
with(pmatChosen, cor.test(EHMacute, het.gen.sd))

## How many simulations were run?
numSims <- rep(NA,finalbatch)
for(bb in 1:finalbatch) {
    load(file=file.path(out.dir, paste0('pmatLs',bb,'.Rdata'))) ## Load last distribution (already filtered)
    numSims[bb] <- sum(unlist(lapply(pmatLs, function(x) nrow(x$pmat))))
}
numSims
sum(numSims)*1.2/60/24
