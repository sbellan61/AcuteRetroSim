####################################################################################################
## Create simulations for first ABC sample
####################################################################################################
library(parallel)
load("data files/ds.nm.all.Rdata") ## DHS country names
load('data files/pars.arr.ac.Rdata')   ## load hazards for each transmission route as fit to DHS data (in.arr[,,2])
load('FiguresAndTables/VL Profile/ehms.vl.Rdata') ## VL EHM estimates for prior
load("data files/pars.arr.ac.Rdata") ## Load pre-couple & extra-couple transmission coefficients estimated from DHS data.
out.dir <- file.path('results','abcSummary')
fig.dir <- file.path('FiguresAndTables','abcFig')
if(!file.exists(fig.dir))      dir.create(fig.dir) # create directory if necessary
if(!file.exists(out.dir))      dir.create(out.dir) # create directory if necessary
ncores <- 12                             ## cores per simulation
vfreq <- 200 ## how often to report on results
s.epic.ind <- s.epic.nm <- NA # not substituting epidemic curves, this will cause rcop() to use default country epidemic curves
## Get fitted transmission coefficients from DHS for prior
pars.arr <- out.arr[,,which(in.arr[,1,2]==7),] 
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") 
ghazs <- c('bb','be','bp') ## geometric mean hazards
mfs <- c('bm2f','em2f','pm2f')
logParms <- c(ghazs,mfs,'acute.sc','dur.ac')
notlogParms <- c('het.gen.sd')
gpars <- pars.arr[ghazs,,13]              # get transmission coefficients from base country
mfratios <- pars.arr[mfs,,13]

## Parameterize by geometric mean hazard b/w genders M2F ratio to get at correlation structure from DHS

## Prior on geometric mean within-couple transmission parameters.  Go more up than down
## because we know we're going to need more hazard if heterogeneity is added, but no reason to need
## less hazard since the estimates are from a homogoneous model which can't get any more homogenous
logGHazLoBound <- logGHazHiBound <- logGHazMeans <- log(gpars[,'50%'])
logGHazLoBound['bp'] <- log(gpars['bp','50%']/5) ## 5X smaller than DHS estimates
logGHazHiBound['bp'] <- log(gpars['bp','50%']*1000) ## 1000X bigger than DHS
## User tighter prior for pre-couple/extra-couple (can't be too big b/c couldn't have DHS data
## results if these were substantially different
logGHazLoBound[c('bb','be')] <- log(gpars[c('bb','be'),'50%']/5) ## 5X smaller than DHS estimates 
logGHazHiBound[c('bb','be')] <- log(gpars[c('bb','be'),'50%']*10) ## 15X bigger than DHS estimates
gpars[,'50%']
exp(rbind(logGHazLoBound,logGHazHiBound))
## parameter estimates on log scale monthly hazard = infections/(coital acts) * (coital acts)/month
## from Wawer. bump is because we know this needs to be bigger than Wawer's since prevalent
## couples have much lower hazard than the baseline in a hetrogeneous model

logGHazPrior <- function(n=1, parms = NULL, low=logGHazLoBound, high=logGHazHiBound) { ## less informative prior, uniform on wide range (on log scale)
    if(is.null(parms[1])){ ## RNG
        logGHazSamp <- runif(3*n, low, high)
        samp <- matrix(exp(logGHazSamp),nr=n, ncol = 3, byrow = T)
        colnames(samp) <- ghazs
        return(samp)
    }else{ ### give prior probability density
        dprior <- log(parms)
        for(ii in 1:ncol(parms)) dprior[,ii] <- dunif(dprior[,ii], low[ii], high[ii])
        return(dprior)
    }}
## Check to make sure it's working
logGHazPrior(1)
logGHazPrior(p=logGHazPrior(1))
apply(logGHazPrior(10^4), 2, function(x) quantile(x, c(.025,.5,.975)))
exp(rbind(logGHazLoBound, logGHazHiBound))
exp(logGHazMeans)
gpars[,'50%']


logM2Fsd <- rowMeans(log(mfratios[,c(2,3)])-log(mfratios[,c(1,2)]))/1.96 ## SD on log scale (using average of 97.5-50 & 50-2.5%)
logM2Fmean <- log(mfratios[,2])
## Priors on log(M2F) ratios for each route
logM2FPrior <- function(n=1, parms=NULL, mu=logM2Fmean, sig=logM2Fsd) {
    if(is.null(parms[1])) {
        logM2FSamp <- rnorm(3*n, mean = mu, sd = sig)
        samp <- matrix(exp(logM2FSamp), nr=n, ncol = 3, byrow = T)
        colnames(samp) <- mfs
        return(samp)
    }else{
        dprior <- log(parms)
        for(ii in 1:ncol(parms)) dprior[,ii] <- dnorm(dprior[,ii], mu[ii], sig[ii])
        return(dprior)
    }}
## Check to make sure it's working
logM2FPrior(3)
logM2FPrior(parms=logM2FPrior(3))
apply(logM2FPrior(10^4), 2, function(x) quantile(x, c(.025,.5,.975)))
mfratios ## looks like it's simulating well!

## Flat uniform prior from 1 to 3
hetGenSDPrior <- function(n, parms=NULL, lo = 0, hi = 3) {
    if(is.null(parms[1])) return(runif(n, lo, hi)) else return(dunif(parms, lo, hi))## RNG
}
x <- hetGenSDPrior(10)
hetGenSDPrior(parms=c(-.2,x,4)) # check its working

## uniform from 1/2X to 50X
rhAcutePrior <- function(n, parms = NULL, lo = .5, hi = 200) {
    if(is.null(parms[1])) {
        samp <- runif(n, log(lo), log(hi))
        return(exp(samp)) 
    }else{
        return(dunif(log(parms), log(lo), log(hi)))
    }}
rhAcutePrior(10)
quantile(rhAcutePrior(10^4), c(.025,.5,.975))

## uniform from 2 weeks to 8 months
durAcutePrior <- function(n, parms = NULL, lo = .5, hi = 8) {
    if(is.null(parms[1])) {
        samp <- runif(n, log(lo), log(hi))
        return(exp(samp)) 
    }else{
        return(dunif(log(parms), log(lo), log(hi)))
    }}
x <- durAcutePrior(10)
durAcutePrior(parms=c(.1,x,4)) # check its working

simParmSamp <- function(n, parms=NULL) {
    if(is.null(parms)) {
        samp <- data.frame(logGHazPrior(n), logM2FPrior(n),
                           acute.sc = rhAcutePrior(n), dur.ac = durAcutePrior(n), het.gen.sd = hetGenSDPrior(n))
        return(samp)
    }else{
        dprior <- parms
        dprior[,ghazs] <- logGHazPrior(parms = parms[,ghazs])
        dprior[,mfs] <- logM2FPrior(parms = parms[,mfs])
        dprior[,'acute.sc'] <- rhAcutePrior(parms = parms[,'acute.sc'])
        dprior[,'dur.ac'] <- durAcutePrior(parms = parms[,'dur.ac'])
        dprior[,'het.gen.sd'] <- hetGenSDPrior(parms = parms[,'het.gen.sd'])
        return(dprior)
    }}

simParmSamp(10)
x <- simParmSamp(10)
simParmSamp(parms=x)
parnms <- names(simParmSamp(1)) ## useful later

simParmConverter <- function(parms) {
    ghazSamp <- log(parms[,ghazs])
    M2FSamp <- log(parms[,mfs])
    mSamp <- ghazSamp + M2FSamp/2
    fSamp <- ghazSamp - M2FSamp/2
    colnames(mSamp) <- c('bmb','bme','bmp')
    colnames(fSamp) <- c('bfb','bfe','bfp')
    hazSamp <- exp(cbind(mSamp,fSamp)[,hazs]) ## order columns
    parms <- data.frame(parms,hazSamp)
    return(parms)
}
simParmConverter(simParmSamp(10)) ## looking good
tst <- apply(simParmConverter(simParmSamp(10^4)), 2, function(x) quantile(x, c(.025,.5,.975))) ## looking good
tst['2.5%',hazs]
t(pars.arr[hazs,2,13]*.2) ## 
tst['97.5%',hazs]
t(pars.arr[hazs,2,13]*20) ## all good going 5x lower & 20x higher

sdPost <- function(pm) { ## get sd of posterior (on appropriate scale)
    pmUnTransf <- pm[,parnms]
    pmUnTransf[,logParms] <- log(pmUnTransf[,logParms])
    return(apply(pmUnTransf, 2, sd))
}

addEHM <- function(x) within(x, EHMacute <- (x[,'acute.sc']-1)*x[,'dur.ac'])

logtransParms <- function(parms,includeEHM=F) {
    Lparms <- parms
    if(includeEHM) {
        Lparms <- cbind(parms, addEHM(parms))
        colnames(Lparms)[ncol(Lparms)] <- 'EHMacute'
        Lparms[,'EHMacute'] <- sapply(Lparms[,'EHMacute'],function(x) max(x,0)) ## careful with this later
        Lparms[,'EHMacute'] <- log(Lparms[,'EHMacute'])
    }
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
perturbParticle <- function(parms, ## matrix of last batch of particles from which we are perturbing
                            from=NULL, ## last batch of particles from which we are perturbing, used to calculate K() kernel PDFs for weights
                            sds, ## std dev of last batch on appropriate transformed scale
                            browse=F) { 
    newLparms <- Lparms <- logtransParms(parms[,parnms,drop=F])
    if(is.null(from)) {
        ii <- 1
        if(browse) browser()
        for(pp in 1:ncol(Lparms)) { ## perturb each parameter column
            newLparms[,pp] <- Lparms[,pp,drop=F] + runif(nrow(Lparms), - sds[pp], sds[pp]) ## others
        }
        newparms <- unlogtransParms(newLparms)
        PriorDensities <- simParmSamp(parms=newparms)
        PriorDensity <- apply(PriorDensities, 1, prod)
        newparms <- newparms[PriorDensity>0,] 
        return(newparms)   
    }else{
        Lparms <- logtransParms(as.matrix(parms[,parnms]))
        dprior <- Lfrom <- logtransParms(as.matrix(from[,parnms]))
        if(browse) browser()
        for(pp in 1:length(Lparms)) {
            dprior[,pp] <- dunif(Lparms[,pp], Lfrom[,pp] - sds[pp], Lfrom[,pp] + sds[pp])
        }
        dpriorProd <- as.numeric(apply(dprior, 1, prod))
        return(dpriorProd)
    }}

## ## Check that when perturbing and then asking for the probability of having been perturbed there
## ## from a parameter set, that we always get the same value (since doing uniform sampling on
## ## (-sds,+sds)
## frm <- pmatChosen[5,]; print(frm[,parnms])
## for(jj in 1:20) {
##     prt <- perturbParticle(frm, sds = sds, browse=F)
##     print(perturbParticle(prt, from = frm, sds = sds))
## }

## ## Check the number of particles from t-1 from which a perturbed particle at time t could have come
## ## from (equivalent to K() kernel since we're doing uniform K())
## for(ii in 1:100) {
##     prt <- perturbParticle(pmatChosen[ii,], sds = sds)
## #    print(sum(perturbParticle(prt, from = pmatChosen, sds = sds)>0))
##     print(sum(pmatChosen$weight * perturbParticle(prt, from = pmatChosen, sds = sds))) ## denominator
## } ## always at least 1 so that's good

weightParticles <- function(currentBatch, lastBatch, sdsUse, browse = F) {
    if(browse) browser()
    dpriors <- simParmSamp(parms=currentBatch[,parnms])
    numerator <- apply(dpriors, 1, prod)
    denominator <- rep(NA,length(numerator))
    ## for each particle, calculate the probability it could have been gotten to from all previous
    ## particles, weighted by their weights
    for(jj in 1:length(numerator)) { 
        Ks <- perturbParticle(currentBatch[jj,], from = lastBatch, sds = sdsUse, browse=F)
        denominator[jj] <- sum(lastBatch$weight * Ks)
    }
    return(numerator/denominator)
}

####################################################################################################
## Function that does simulation & gets wtab all at once
retroCohSim <- function(parms=simParmSamp(1), maxN=10000, seed=1, nc=12, browse=F) {
    startTime <- Sys.time()
    parms <- simParmConverter(parms)
    output <- with(parms, psrun(maxN = maxN, jobnum = seed, pars = parms[hazs], save.new = T, return.ts = T, returnFileNm=F, saveFile=F, seed = seed,
                                    acute.sc = acute.sc, dur.ac = dur.ac, het.gen=T, het.gen.sd = het.gen.sd, browse = F, nc = 12))
    cohsim <- rak.coh.fxn(output, #ts.ap = output$ts, dat = output$evout, dpars = output$rakpars,
                          ltf.prob=0.0287682072451781, 
                          rr.ltf.ff=1.5, rr.ltf.mm=1.5, rr.ltf.hh=1, rr.ltf.d=1, rr.inc.sdc=1.5,
                          verbose = F, browse = F)
                                        #    rm(output); gc() ## free up memory
    rcohsim <- rak.wawer(rak.coh = cohsim, excl.extram=T, het.gen.sd = with(parms,het.gen.sd),
                         cov.mods = F, verbose = F, fit.Pois=F, prop.controlled=c(NA,1), browse=F)
    rcohsim$pars <- parms
    minutesTaken <- as.numeric(difftime(Sys.time(), startTime, units='mins'))
    print(paste('minutes taken', round(minutesTaken,2)))
    return(rcohsim)
}

gStat <- function(x) { ## G test of independence
    nn <- sum(x)
    sr <- rowSums(x)
    sc <- colSums(x)
    E <- outer(sr, sc, "*")/nn
    gs <- x*log(x/E)
    gs[x==0] <- 0
    2*sum(gs)
}

gSumStat <- function(wtab) { ## compare wtab to wtab.rl
    if(!is.na(wtab[1])) {
        prevG <- incG <- numeric(4)
        ## inc
        if(nrow(wtab$inct)<4) { ## deal with less rows
            nr <- nrow(wtab$inct)
            for(rr in (nr+1):4) wtab$inct[rr,] <- c(NA,NA)
        }
        for(ii in 1:4) {
            incG[ii] <- gStat(rbind(contTabsRl$inct[ii,], wtab$inct[ii,]))
        }
        incG[is.nan(incG)] <- 0
        ## prev
        if(nrow(wtab$prevt)<4) { ## deal with less rows
            nr <- nrow(wtab$prevt)
            for(rr in (nr+1):4) wtab$prevt[rr,] <- c(NA,NA)
        }
        for(ii in 1:4) {
            prevG[ii] <- gStat(rbind(contTabsRl$prevt[ii,], wtab$prevt[ii,]))
        }
        prevG[is.nan(prevG)] <- 0
        Gtable <- abind(cbind(contTabsRl$inct, wtab$inct, incG), cbind(contTabsRl$prevt, wtab$prevt, prevG),  along = 3)
        gtab <- cbind(incG,prevG)
        return(list(Gval = sum(gtab,na.rm=T), Gtable))
    }else{    return(list(Gval = NA, Gtable=NA)) }
}

reweightG <- function(x, Gweights=cbind(rep(1,4),rep(1,4))) {
    x[[1]][,'prevG',] <- try(x[[1]][,'prevG',]*Gweights)
    return(x)
}

sumGs <- function(x) sum(x[[1]][,'prevG',], na.rm=T)

contTabFxn <- function(x) within(x, { ## turn wtab into only having # infected & # uninfected
    if(!is.na(as.vector(inct)[1])) {
        inct$ni <- with(inct, n-i)
        inct <- inct[,c('i','ni')]
        prevt$ni <- with(prevt, n-i)
        prevt <- prevt[,c('i','ni')]
    }else{ return(NA) }
})

propExcluded <- function(x) with(x, { ## proportion excluded from simulated couples (not conditioning on Rakai sample size)
    incSDCtotal <- sbmod.to.wdat(rakll, excl.by.err=F, condRakai=F, giveLate=F)$inct[1,'n']
    incSDCincl <- sbmod.to.wdat(rakll, excl.by.err=T, condRakai=F, giveLate=F)$inct[1,'n']
    propExcl <- (incSDCtotal-incSDCincl)/incSDCtotal
    return(propExcl)
})

numExcluded <- function(propEx) 23 / (1-propEx) - 23 ## based on that proportion from the total, how many were excluded to get down to 23

collectSumStat <- function(filenm, returnGtable = F, browse=F, ncores = ncores, rmNull=T) {
    if(browse) browser()
    print(filenm)
    load(filenm)
    ## if returning no output because of non-plausible (not enough couples in inc/prev or 0 infected
    ## inc)
    plausible <- !unlist(lapply(rcohsList, is.null)) 
    rcohsList <- rcohsList[plausible]
    ## Proportion seroincident excluded in ABC-SMC fits
    propExcl <- unlist(mclapply(rcohsList, propExcluded, mc.cores = ncores))    
    wtabSims <- mclapply(rcohsList, function(x) sbmod.to.wdat(x$rakll, excl.by.err = T, browse=F, giveLate=F, 
                                                              condRakai=T, giveProp=T, simpPois=T), mc.cores=ncores)
    ## Convert to wtab
    wtabSims <- mclapply(rcohsList, function(x) sbmod.to.wdat(x$rakll, excl.by.err = T, browse=F, giveLate=F, 
                                                              condRakai=T, giveProp=T, simpPois=T), mc.cores=ncores)
    prevHazs <- unlist(lapply(wtabSims, function(x) if(is.null(x$error)) with(x$prevt, sum(i)/sum(n)) else NA))
    nCplMat <- do.call(rbind.data.frame, mclapply(rcohsList, function(x) xtabs(~phase, x$rakll)))
    colnames(nCplMat) <- c('prev','inc','late')
    parmsMat <- matrix(unlist(lapply(rcohsList, '[[', 'pars')), nr = length(rcohsList), nc = 15, byrow=T)
    PoisRHsMat <- data.frame(matrix(unlist(lapply(wtabSims, '[[', 'PoisRHs')), nr = length(rcohsList), nc = 2, byrow=T))
    colnames(parmsMat) <- names(rcohsList[[1]]$pars)
    colnames(PoisRHsMat) <- c('univ','omn')
    PoisRHsMat <- within(PoisRHsMat, {
        RHreduction <- univ/omn
        enoughHet <- RHreduction > 11/7.25      ##  univariate unadjusted / multivariate from Wawer paper
    })
    rm(rcohsList); gc()
    ## just get # infected & # uninfected
    contTabsSim <- suppressWarnings(mclapply(wtabSims, contTabFxn, mc.cores=ncores)) ## suppressing length 1 for is.na(), not problematic
    ## Get G statistics
    gS <- mclapply(contTabsSim, gSumStat, mc.cores = ncores)
    gVals <- unlist(lapply(gS, '[',1))
    ## Create parameter matrix
    parmsMat <- data.frame(parmsMat, nCplMat, prevHazs, PoisRHsMat, gVals, propExcl = propExcl,
                           filenm = filenm, job = 1:nrow(parmsMat))
    if(returnGtable) Gtable <- lapply(gS, '[',2) else Gtable <- NA
    rm(contTabsSim, gS,gVals); gc()
    return(list(pmat=parmsMat, Gtable=Gtable))
}


## to delete a range of jobs
## qdel echo `seq -f "%.0f" 2282389 2282404`
