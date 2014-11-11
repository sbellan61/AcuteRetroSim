setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
rm(list=ls()); gc()
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
library(parallel)
ncores <- 12
fig.dir <- file.path('FiguresAndTables','abcFig')
out.dir <- file.path('results','abcSummary')
if(!file.exists(fig.dir))      dir.create(fig.dir) # create directory if necessary
if(!file.exists(out.dir))      dir.create(out.dir) # create directory if necessary
batch <- 1
in.dir <- paste0('results/abcBatch',batch)
#in.dir <- 'results/abcBatch1Old5good/'

fls <- list.files(in.dir, pattern='Rdata', full.names=T)
## tst <- collectSumStat(fls[16], browse=T, returnGtable=T, ncores=1)
## pmatLs <- lapply(fls, collectSumStat, returnGtable=T) ## not enough mem to do more than a few cores

## Collect results and put them in a matrix and a list of Wawer-type tables
pmatLs <- mclapply(fls, collectSumStat, returnGtable=T, mc.cores=12, ncores = 1)
pmat <- do.call(rbind.data.frame, lapply(pmatLs, '[[', 'pmat'))
gtabs <- do.call(c, lapply(pmatLs, '[[', 'Gtable'))
ord <- order(pmat$gVals)
pmat <- pmat[ord,]
gtabs <- gtabs[ord]
pmat <- within(pmat, {for(ll in logParms) assign(paste0('log',ll), log(get(ll))); rm(ll)})
save(gtabs,pmat, pmatLs, file = file.path(out.dir,paste0('pmatLs',batch,'.Rdata')))
nrow(pmat)
sum(duplicated(pmat[,parnms])) ## make sure we have unique parm draws

## Look at Gtables for the best ones to make sure they seem to be good matchest
head(pmat[,c('acute.sc','dur.ac','het.gen.sd','enoughHet','inc','prev','gVals')],10)
gtabs[1]
## load(file = file.path(out.dir,paste0('pmatLs',batch,'.Rdata')))
pmat <- pmat[!is.na(pmat$gVals),]


## Histogram of g value statistics
pdf(file.path(fig.dir, paste0('hist gval',batch,'.pdf')))
hist(pmat$gVals, breaks = 0:max(pmat$gVals+1), col = 'black', xlab = 'G stat', main='')
graphics.off()

## Look at hazard & number of couple parameters collinearity
sbpairs(pmat[,c('bb','inc')], file.path(fig.dir, paste0('bbInc',batch)), do.jpeg=T)
sbpairs(pmat[,c('be','inc')], file.path(fig.dir, paste0('beInc',batch)), do.jpeg=T)
sbpairs(pmat[,c('bp','inc')], file.path(fig.dir, paste0('bpInc',batch)), do.jpeg=T)
sbpairs(pmat[,c('bb','prev')], file.path(fig.dir, paste0('bbPrev',batch)), do.jpeg=T)
sbpairs(pmat[,c('be','prev')], file.path(fig.dir, paste0('bePrev',batch)), do.jpeg=T)
sbpairs(pmat[,c('bp','prev')], file.path(fig.dir, paste0('bpPrev',batch)), do.jpeg=T)


## Select intermediate distribution based on threshold criteria
qchisq(.95, 7) ## based on having 6-8ish gVal comparisons each of which is approximately chisq distr
cutfs <- c(15,8,5,2,1) ## cutoffs for each batch
cutf <- cutfs[batch]
print(paste0(round(mean(pmat$gVals<cutf)*100),'% of simulations pass threshold criteria'))
rightSizeNinc <- with(pmat, inc > 20 & inc < 80)
rightSizeNprev <- with(pmat, prev > 100 & prev < 500)
choice <- with(pmat, gVals<cutf & enoughHet & rightSizeNinc & rightSizeNprev)
print(paste0(round(mean(choice)*100),'% of simulations pass threshold criteria with # of couples restrictions'))
pmatChosen <- pmat[choice,]

gtabs[nrow(pmatChosen)+0:-3] ## Look at last few Gtables to see how the worst simulations selected fare

## Compare prior distributions to intermediate distribution
priorParms <- simParmSamp(10^5) ## sample prior
sel <- c(ghazs,'acute.sc','dur.ac','het.gen.sd')
rgs <- apply(priorParms[,sel],2,range)
sbpairs(pmatChosen[,sel], file.path(fig.dir, paste0('Post',batch)), do.jpeg=T, rgs=rgs)
sbpairs(priorParms[, sel], file.path(fig.dir, 'Prior'), do.jpeg=T, rgs=rgs)
graphics.off()

## Log scale
rgs <- apply(logtransParms(priorParms),2,range) 
sbpairs(logtransParms(pmatChosen[,parnms]), file.path(fig.dir, paste0('logPost',batch)), do.jpeg=T, rgs=rgs)
sbpairs(logtransParms(priorParms), file.path(fig.dir, 'logPrior'), do.jpeg=T, rgs=rgs)
graphics.off()

## Look at CIs
apply(logtransParms(pmatChosen[,parnms]), 2, function(x) quantile(x,c(.025,.5,.975)))
apply(logtransParms(priorParms[,parnms]), 2, function(x) quantile(x,c(.025,.5,.975)))

## Calculate std dev to use in particle perturbations
sds <- sdPost(pmatChosen) ## log scale

## Weight particles
if(batch==1) pmatChosen$weight <- 1/nrow(pmatChosen) else{
    pmatChosen$weight <- weightParticles(pmatChosen,pmatChosen, T)
}
save(pmatChosen, sds, file=file.path(out.dir, paste0('IntermedDistr',batch,'.Rdata'))) ## Save particles & their sds



## ## Check it's working (should always have same uniform density of getting to perturbed particle from 'from' particle
## frm <- pmatChosen[5,]; print(frm[,parnms])
## for(jj in 1:100) {
##     prt <- perturbParticle(frm, sds = sds, browse=F)
##     print(perturbParticle(prt, from = frm, sds = sds))
## }

## for(ii in 1:100) {prt <- perturbParticle(pmatChosen[ii,], sds = sds)
## sum(perturbParticle(prt, from = pmatChosen, sds = sds)>0)}
