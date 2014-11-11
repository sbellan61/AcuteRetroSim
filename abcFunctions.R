####################################################################################################
## Create simulations for first ABC sample
####################################################################################################
load("data files/ds.nm.all.Rdata") ## DHS country names
load('data files/pars.arr.ac.Rdata')   ## load hazards for each transmission route as fit to DHS data (in.arr[,,2])
load('FiguresAndTables/VL Profile/ehms.vl.Rdata') ## VL EHM estimates for prior
load("data files/pars.arr.ac.Rdata") ## Load pre-couple & extra-couple transmission coefficients estimated from DHS data. 
nc <- 12                             ## cores per simulation
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

pdf('FiguresAndTables/abcFig/RHprior.pdf',5.5,4)
par(mar=c(5,4,1,.5), 'ps'=10)
xs <- exp(seq(log(.5),log(500),l=500))
ys <- rhAcutePrior(parms=xs)
plot(xs,ys, type = 'l', lwd = 2, log='x', ylab = 'probability density', xlab=expression(RH[acute]), main = '', las = 1, bty = 'n', xlim=c(.5,500))
graphics.off()

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

####################################################################################################
## Function that does simulation & gets wtab all at once
retroCohSim <- function(parms=simParmSamp(1), maxN=10000, seed=1, nc=12, browse=F) {
    startTime <- Sys.time()
    parms <- simParmConverter(parms)
    output <- with(parms, psrun(maxN = maxN, jobnum = seed, pars = parms[hazs], save.new = T, return.ts = T, returnFileNm=F, saveFile=F, seed = seed,
                                    acute.sc = acute.sc, dur.ac = dur.ac, het.gen=T, het.gen.sd = het.gen.sd, browse = F, nc = nc))
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
        return(list(Gval = sum(incG,prevG, na.rm=T), Gtable))
    }else{    return(list(Gval = NA, Gtable=NA)) }
}

contTabFxn <- function(x) within(x, { ## turn wtab into only having # infected & # uninfected
    if(!is.na(as.vector(inct)[1])) {
        inct$ni <- with(inct, n-i)
        inct <- inct[,c('i','ni')]
        prevt$ni <- with(prevt, n-i)
        prevt <- prevt[,c('i','ni')]
    }else{ return(NA) }
})

collectSumStat <- function(filenm, returnGtable = F, browse=F, ncores = 12, rmNull=T) {
    if(browse) browser()
    print(filenm)
    load(filenm)
    ## if returning no output because of non-plausible (not enough couples in inc/prev or 0 infected
    ## inc)
    plausible <- !unlist(lapply(rcohsList, is.null)) 
    rcohsList <- rcohsList[plausible]
    ## Convert to wtab
    wtabSims <- mclapply(rcohsList, function(x) sbmod.to.wdat(x$rakll, excl.by.err = T, browse=F, giveLate=F, 
                                                              condRakai=T, giveProp=T, simpPois=T), mc.cores=ncores)
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
    parmsMat <- data.frame(parmsMat, nCplMat, PoisRHsMat, gVals, filenm = filenm, job = 1:nrow(parmsMat))
    if(returnGtable) Gtable <- lapply(gS, '[',2) else Gtable <- NA
    rm(contTabsSim, gS,gVals); gc()
    return(list(pmat=parmsMat, Gtable=Gtable))
}


## to delete a range of jobs
## qdel echo `seq -f "%.0f" 2282389 2282404`
