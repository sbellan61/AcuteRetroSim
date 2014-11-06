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
spars <- pars.arr[hazs,,13]              # get transmission coefficients from base country

## Prior on pre-, extra-, and within-couple transmission parameters
bump <- 5
logHazMeans <- log(spars[,'50%']*bump)
logHazSDs <- rowMeans(abs(log(spars[,c('2.5%','97.5%')])-log(spars[,'50%']))) / 1.96 ## ~1.96*SD of
## parameter estimates on log scale monthly hazard = infections/(coital acts) * (coital acts)/month
## from Wawer. bump is because we know this needs to be bigger than Wawer's since prevalent
## couples have much lower hazard than the baseline in a hetrogeneous model
logHazPrior <- function(n=1, scl = 1, mean = logHazSDs, sd = logHazSDs, extraUncertaintyHazs = 1.5) {
    logHazSamp <- rnorm(6*n, logHazMeans, sd = logHazSDs*extraUncertaintyHazs)
    return(matrix(exp(logHazSamp),nr=n, ncol = 6, byrow = T))
}
tst <- logHazPrior(1000)
apply(tst, 2, function(x) quantile(x, c(.025,.5,.975)))
exp(logHazMeans)
spars[,'50%']

## Flat uniform prior from 1 to 3
hetGenSDPrior <- function(n) runif(n, 0, 3)

## Get RHacute prior centered on estimates from viral loads but with way more uncertainty
extraUncertaintyRH <- 3.5
log(ehms.vl['med'])-log(ehms.vl[c('lci','uci')]) ## How far away are bounds on a log scale? about .5
sd_RHacutePrior <- .5/1.96 * extraUncertaintyRH ## add a bunch of extra uncertainty since this is just VL
exp(qnorm(c(.025,.975), log(ehms.vl['med']), sd = sd_RHacutePrior)) ## 95% CI bounds: .97 - 32X as infectious
rhAcutePrior <- function(n) exp(rnorm(n, log(ehms.vl['med']), sd = sd_RHacutePrior))

## uniform from 2 weeks to 8 months
durAcutePrior <- function(n) runif(n, .5, 8)

simParmSamp <- function(n,...) {
    samp <- data.frame(logHazPrior(n), acute.sc = rhAcutePrior(n), dur.ac = durAcutePrior(n), het.gen.sd = hetGenSDPrior(n))
    names(samp)[1:6] <- hazs
    return(samp)
}

####################################################################################################
## Function that does simulation & gets wtab all at once
retroCohSim <- function(parms=simParmSamp(1), maxN=10000, seed=1, nc=12, browse=F) {
    startTime <- Sys.time()
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

abcSimSumStat <- function(rcohsim) {
    wtabSim <- sbmod.to.wdat(rcohsim$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=T)
    return(list(wtabSim=wtabSim, parms=parms, minutesTaken=))
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
}

collectSumStat <- function(filenm, returnGtable = F, browse=F, ncores = 12) {
    if(browse) browser()
    print(filenm)
    load(filenm)
    ## Convert to wtab
    wtabSims <- mclapply(rcohsList, function(x) sbmod.to.wdat(x$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=T, giveProp=T, simpPois=T),
                         mc.cores=ncores)
    parmsMat <- matrix(unlist(lapply(rcohsList, '[[', 'pars')), nr = length(rcohsList), nc = 9, byrow=T)
    PoisRHsMat <- data.frame(matrix(unlist(lapply(wtabSims, '[[', 'PoisRHs')), nr = length(rcohsList), nc = 2, byrow=T))
    colnames(parmsMat) <- names(rcohsList[[1]]$pars)
    colnames(PoisRHsMat) <- names(wtabSims[[1]]$PoisRHs)
    PoisRHsMat <- within(PoisRHsMat, {
        RHreduction <- univ.phaseinc/omn.phaseinc
        enoughHet <- RHreduction > 11/7.25}) ## univariate unadjusted / multivariate from Wawer paper
    rm(rcohsList); gc()
    ## just get # infected & # uninfected
    contTabsSim <- mclapply(wtabSims, function(x) {
        x <- within(x, {
            inct$ni <- with(inct, n-i)
            inct <- inct[,c('i','ni')]
            prevt$ni <- with(prevt, n-i)
            prevt <- prevt[,c('i','ni')]})}, mc.cores=ncores)
    ## Get G statistics
    gS <- mclapply(contTabsSim, gSumStat, mc.cores = ncores)
    gVals <- unlist(lapply(gS, '[',1))
    ## Create parameter matrix
    parmsMat <- data.frame(parmsMat, PoisRHsMat, gVals, filenm = filenm, job = 1:nrow(parmsMat))
    pmat <- parmsMat[order(gVals),]
    if(returnGtable) Gtable <- unlist(lapply(gS, '[',2)) else Gtable <- NA
    rm(parmsMat, contTabsSim, gS,gVals); gc()
    return(list(pmat=pmat, Gtable=Gtable))
}

