####################################################################################################
## Collect batch of ABC results & prepare for next batch & do some diagnostics.
## Steve Bellan, 2014, steve.bellan@gmail.com
####################################################################################################
rm(list=ls()); gc()
setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
ncores <- 12
batch <- 5
in.dir <- paste0('results/abcBatch',batch) ## collect results from here (not in repo, too big)

fls <- list.files(in.dir, pattern='Rdata', full.names=T)
## tst <- collectSumStat(fls[1], browse=F, returnGtable=T, ncores=1) ## help with debugging, collect one result

## Collect results and put them in a matrix and a list of Wawer-type tables
notZeroSize <- file.info(fls)$size!=0
pmatLs <- mclapply(fls[notZeroSize], collectSumStat, returnGtable=T, mc.cores=12, ncores = 1)
save(pmatLs, file = file.path(out.dir,paste0('pmatLs',batch,'.Rdata')))
pmat <- do.call(rbind.data.frame, lapply(pmatLs, '[[', 'pmat'))
gtabs <- do.call(c, lapply(pmatLs, '[[', 'Gtable'))
ord <- order(pmat$gVals)
pmat <- pmat[ord,]
gtabs <- gtabs[ord]
pmat <- within(pmat, {for(ll in logParms) assign(paste0('log',ll), log(get(ll))); rm(ll)})
nrow(pmat)
sum(duplicated(pmat[,parnms])) ## make sure we have unique parm draws
save(gtabs,pmat, pmatLs, file = file.path(out.dir,paste0('pmatLs',batch,'.Rdata'))) ## not in repo, too big

## Look at Gtables for the best ones to make sure they seem to be good matchest
head(pmat[,c('acute.sc','dur.ac','het.gen.sd','enoughHet','inc','prev','gVals')],10)
gtabs[1]

## load(file = file.path(out.dir,paste0('pmatLs',batch,'.Rdata')))
excludeNA <- !is.na(pmat$gVals)
pmat <- pmat[excludeNA,]
gtabs <- gtabs[excludeNA]

## Histogram of g value statistics
pdf(file.path(fig.dir, 'hist gval', paste0('hist gval',batch,'.pdf')))
hist(pmat$gVals, breaks = 0:max(pmat$gVals+1), col = 'black', xlab = 'G stat', main='')
graphics.off()

## Look at hazard & number of couple parameters collinearity
bdir <- file.path(fig.dir,'bs')
sbpairs(pmat[,c('bb','inc')], file.path(bdir, paste0('bbInc',batch)), do.jpeg=T)
sbpairs(pmat[,c('be','inc')], file.path(bdir, paste0('beInc',batch)), do.jpeg=T)
sbpairs(pmat[,c('bp','inc')], file.path(bdir, paste0('bpInc',batch)), do.jpeg=T)
sbpairs(pmat[,c('bb','prev')], file.path(bdir, paste0('bbPrev',batch)), do.jpeg=T)
sbpairs(pmat[,c('be','prev')], file.path(bdir, paste0('bePrev',batch)), do.jpeg=T)
sbpairs(pmat[,c('bp','prev')], file.path(bdir, paste0('bpPrev',batch)), do.jpeg=T)

## Select intermediate distribution based on threshold criteria
qchisq(.95,8) ## for 8 df, < 15 suggests simulation & data are not significantly different. We'll be slightly more conservative.
cutfs <- c(15,15,10,10,10) ## cutoffs for each batch
cutf <- cutfs[batch]
## Match G stats
if(batch>4) { ## weight G table for preferential fitting of first interval observations to
    ## counterbalance the overweighting of later ones caused by the other criteria &
    ## loss-to-follow-up mismatch
    gtabsWt <- lapply(gtabs, reweightG, Gweights=cbind(c(40,1,1,1),c(3,1,1,1))) ## reweight them
    pmat$gVals <- unlist(lapply(gtabsWt, sumGs))
    ord <- order(pmat$gVals)
    pmat <- pmat[ord,]
    gtabs <- gtabs[ord]
}
Gsel <- pmat$gVals<cutf 
print(paste0(round(mean(Gsel)*100),'% of simulations pass G statistic threshold criteria'))
## Match number of couples
nprevsMin <- c(100,110,rep(120,3))
nprevsMax <- c(500,400,rep(300,3))
nincsMin <- c(20,rep(23,4))
nincsMax <- c(80,rep(60,4))
rightSizeNinc <- with(pmat, inc > nincsMin[batch] & inc < nincsMax[batch])
rightSizeNprev <- with(pmat, prev > nprevsMin[batch] & prev < nprevsMax[batch])
sum(rightSizeNinc & rightSizeNprev)
## Match prevalent couple hazards
WawPrevHaz <- with(wtab.rl$prevt, sum(i)/sum(n)) ## .0843
prevHazMin <- WawPrevHaz*c(.01, .7,.8,.85,.8)
prevHazMax <- WawPrevHaz/c(.01, .7,.8, .85,.8)
rightPrevHaz <- with(pmat, prevHazs > prevHazMin[batch] & prevHazs < prevHazMax[batch])
sum(rightPrevHaz)
## Match incident to prevalent couple hazard ratios
RHunivMin <- c(0,11*c(.7,.8,.83,.8)) ## Wawer univariate RH was 11, so try to hone in on this too
RHunivMax <- c(Inf,11/c(.7,.8,.83,.8)) ## Wawer univariate RH was 11, so try to hone in on this too
rightRHuniv <- with(pmat, univ > RHunivMin[batch] & univ < RHunivMax[batch])
sum(rightRHuniv)
## Combine criteria
choice <- with(pmat, Gsel & enoughHet & rightSizeNinc & rightSizeNprev & rightRHuniv & rightPrevHaz)
print(paste0(round(mean(choice)*100),'% of simulations pass all threshold criteria'))
print(sum(choice))

pmatChosen <- pmat[choice,]
gtabs <- gtabs[choice]
tail(pmatChosen,3)
head(pmatChosen,3)
gtabs[nrow(pmatChosen)] ## Look at last few Gtables to see how the worst simulations selected fare
priorParms <- addEHM(simParmSamp(4*10^4)) ## sample prior
pmatChosen <- addEHM(pmatChosen)

## Look at CIs nonlog
apply(pmatChosen[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
# apply(priorParms[,c(parnms,'EHMacute')], 2, function(x) quantile(x,c(.025,.5,.975)))
with(pmatChosen, cor.test(EHMacute, het.gen.sd))

## Compare prior distributions to intermediate distribution
if(batch==1) { ##Plot Priors for comparison (call it Post0 to make name switching easy)
    sbpairs(priorParms[, sel], file.path(fig.dir, 'Post', 'Post0'), do.jpeg=T, rgs=rgs) ## Post0 = Prior
    sbpairs(logtransParms(priorParms), file.path(fig.dir, 'logPost', 'logPost0'), do.jpeg=T, rgs=rgs)
    sbpairs(logtransParms(priorParms,T)[,sel], file.path(fig.dir,  'logPostEHM', 'logPostEHM0'), do.jpeg=T, rgs=rgs)
}

## Plot intermediate distribution
## Untransformed
sel <- c(ghazs,'acute.sc','dur.ac','het.gen.sd','EHMacute')
rgs <- apply(priorParms[,sel],2,range)
sbpairs(pmatChosen[,sel], file.path(fig.dir, 'Post', paste0('Post',batch)), do.jpeg=T, rgs=rgs)
## Log scale
rgs <- apply(logtransParms(priorParms),2,range) 
sbpairs(logtransParms(pmatChosen[,parnms]), file.path(fig.dir, 'logPost', paste0('logPost',batch)), do.jpeg=T, rgs=rgs)
## with EHMacute
sel <- c('acute.sc','dur.ac','EHMacute','het.gen.sd','bp')
rgs <- apply(logtransParms(priorParms,T)[,sel],2,range)
rgs[1,'EHMacute'] <- 0
sbpairs(logtransParms(pmatChosen[,parnms],T)[,sel], file.path(fig.dir, 'logPostEHM', paste0('logPostEHM',batch)), do.jpeg=T, rgs=rgs)

## Calculate std dev to use in particle perturbations
sdsNew <- sdPost(pmatChosen) ## log scale
gtabsNew <- gtabs

## Weight particles
if(batch==1) {
    pmatChosen$weight <- 1/nrow(pmatChosen) 
}else{
    pmatCNew <- pmatChosen ## new batch output
    load(file = file.path(out.dir,paste0('IntermedDistr',batch-1,'.Rdata')))
    pmatCOld <- pmatChosen ## last batch output
    sdsOld <- sds
    pmatCNew$weight <- weightParticles(pmatCNew,pmatCOld, sdsUse = sdsOld, browse=F)
    pmatCNew$weight <- pmatCNew$weight/sum(pmatCNew$weight)
}
pmatChosen <- pmatCNew
sds <- sdsNew ## kept name as New to avoid loading over it above
gtabs <- gtabsNew
save(pmatChosen, sds, gtabs, file=file.path(out.dir, paste0('IntermedDistr',batch,'.Rdata'))) ## Save particles & their sds to results directory
save(pmatChosen, sds, gtabs, file=file.path('abcSummary', paste0('IntermedDistr',batch,'.Rdata'))) ## Save particles & their sds to local directory in repository
