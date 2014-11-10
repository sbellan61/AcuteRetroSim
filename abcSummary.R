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
in.dir <- 'results/abcBatch1Old5good/'

fls <- list.files(in.dir, pattern='Rdata', full.names=T)
tst <- collectSumStat(fls[2], browse=F, returnGtable=T, ncores=1)
##pmatLs <- lapply(fls, collectSumStat, returnGtable=T) ## not enough mem to do more than a few cores

## ncores goes to collectSumStat, use each core to load rather
pmatLs <- mclapply(fls, collectSumStat, returnGtable=T, mc.cores=12, ncores = 1)
unique(do.call(colnames, lapply(pmatLs, '[[', 'pmat')))
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
sel <- c('acute.sc','dur.ac','het.gen.sd','gVals')
head(pmat[,c('acute.sc','dur.ac','het.gen.sd','enoughHet','inc','prev','gVals')],10)
gtabs[1]

load(file = file.path(out.dir,paste0('pmatLs',batch,'.Rdata')))
pmat <- pmat[!is.na(pmat$gVals),]

pdf(file.path(fig.dir, 'hist gval.pdf'))
hist(pmat$gVals, breaks = 0:max(pmat$gVals+1), col = 'black', xlab = 'G stat', main='')
graphics.off()

qchisq(.95, 7) ## based on having 6-8ish gVal comparisons each of which is approximately chisq distr
nrow(pmat)
cutf <- 4 ## above cannot deal with small sample sizes, so let's be conservative
sum(pmat$gVals<cutf)
rightSizeNinc <- with(pmat, Ninc > 20 & Ninc < 80)
rightSizeNprev <- with(pmat, Nprev > 100 & Nprev < 500)
pmatChosen <- pmat[with(pmat, gVals<cutf & enoughHet),]
                                        #pmatChosen <- pmat[with(pmat, gVals<cutf & enoughHet & rightSizeNinc & rightSizeNprev),]
head(pmatChosen)
## Look at last few Gtables
gtabs[nrow(pmatChosen)+0:-3]

priorParms <- simParmSamp(10^5)
priorParms <- within(priorParms, {for(ll in logParms) assign(paste0('log',ll), log(get(ll))); rm(ll)})
head(priorParms)

sel <- c(ghazs,'acute.sc','dur.ac','het.gen.sd')
rgs <- apply(priorParms[,sel],2,range)
sbpairs(pmatChosen[,sel], file.path(fig.dir, paste0('Post',batch)), do.jpeg=T, rgs=rgs)
sbpairs(priorParms[, sel], file.path(fig.dir, 'Prior'), do.jpeg=T, rgs=rgs)
graphics.off()

logsel <- c(paste0('log',ghazs), 'logacute.sc','logdur.ac','het.gen.sd')
rgs <- apply(priorParms[,logsel],2,range)
sbpairs(pmat[with(pmat, gVals<cutf & enoughHet), logsel], file.path(fig.dir, paste0('logPost',batch)), do.jpeg=T, rgs=rgs)
sbpairs(priorParms[, logsel], file.path(fig.dir, 'logPrior'), do.jpeg=T, rgs=rgs)
graphics.off()


logsel <- c(paste0('log',ghazs), 'logacute.sc','logdur.ac','het.gen.sd')
sdPost <- function(pm) { ## get sd of posterior (on appropriate scale)
    pmUnTransf <- pm[,parnms]
    pmUnTransf[,logParms] <- log(pmUnTransf[,logParms])
    return(apply(pmUnTransf, 2, sd))
}

max(priorParms$bp)
## Look at CIs
apply(pmatChosen[,logsel], 2, function(x) quantile(x,c(.025,.5,.975)))
sds <- sdPost(pmatChosen)
apply(priorParms[,logsel], 2, function(x) quantile(x,c(.025,.5,.975)))

## Weight particles
if(batch==1) pmatChosen$weight <- 1/nrow(pmatChosen) else{
    


}

weightParticles(pmatChosen,pmatChosen, T)

               
weightParticles <- function(currentBatch, lastBatch, browse = F) {
    if(browse) browser()
    head(currentBatch[,parnms])
    dpriors <- simParmSamp(parms=currentBatch[,parnms])
    head(dpriors)
    dpi <- apply(dpriors, 1, prod)
    ## for each particle, calculate the probability it could have been gotten to from all previous
    ## particles, weighted by their weights
    for(jj in 1:length(dpriors)) { 
        lastBatch$weight * sapply(parmsChosen[pars], perturbParticle)
    }}

logtransParms <- function(parms) {
    Lparms <- parms
    Lparms[,logParms] <- log(parms[,logParms])
    return(Lparms)
}

unlogtransParms <- function(Lparms) {
    parms <- Lparms
    parms[,logParms] <- exp(Lparms[,logParms])
    return(parms)
}

## perturb a particle (must be done on appropriate scale). If from given, then calculute probability
## density of perturbing to parms from 'from' for kernel in weight denominator
perturbParticle <- function(parms, from=NULL, sds, browse=F) { ## from can be a matrix, to calculate probability of getting to parms from all particles in from
    newLparms <- Lparms <- logtransParms(parms[,parnms])
    if(is.null(from)) {
        zeroPriorDensity <- T
        ii <- 1
        if(browse) browser()
        while(zeroPriorDensity) { ## sample until a particle with nonzero prior density is found
            newLparms <- Lparms + runif(length(parms), - sds, sds) ## others
            ii <- ii+1
            #if(ii>50) browser()
            newparms <- unlogtransParms(newLparms)
            if(prod(simParmSamp(parms=newparms))>0) zeroPriorDensity <- FALSE
        }
        return(newparms)   
    }else{
        Lparms <- logtransParms(as.matrix(parms[,parnms]))
        dprior <- Lfrom <- logtransParms(as.matrix(from[,parnms]))
        if(browse) browser()
        for(pp in 1:length(parms)) dprior[,pp] <- dunif(Lparms[,pp], Lfrom[,pp] - sds[pp], Lfrom[,pp] + sds[pp])
        dpriorProd <- as.numeric(apply(dprior, 1, prod))
        return(dpriorProd)
    }}

## Check it's working (should always have same uniform density of getting to perturbed particle from 'from' particle
frm <- pmatChosen[5,]; print(frm[,parnms])
for(jj in 1:100) {
    prt <- perturbParticle(frm, sds = sds, browse=F)
    print(perturbParticle(prt, from = frm, sds = sds))
}

for(ii in 1:100) {prt <- perturbParticle(pmatChosen[ii,], sds = sds)
sum(perturbParticle(prt, from = pmatChosen, sds = sds)>0)}
