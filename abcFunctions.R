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
logHazMeans <- log(spars[,'50%']* 1.5)
logHazSDs <- rowMeans(abs(log(spars[,c('2.5%','97.5%')])-log(spars[,'50%']))) / 1.96 ## ~1.96*SD of
## parameter estimates on log scale monthly hazard = infections/(coital acts) * (coital acts)/month
## from Wawer. 1.3 bump is because we know this needs to be bigger than Wawer's since prevalent
## couples have much lower hazard than the baseline in a hetrogeneous model
logHazMeans[c('bmp','bfp')] <- log(36/48525 * 10 * 1.3) 
logHazPrior <- function(n=1, scl = 1, mean = logHazSDs, sd = logHazSDs, extraUncertaintyHazs = 1.5) {
    logHazSamp <- rnorm(6*n, logHazMeans, sd = logHazSDs*extraUncertaintyHazs)
    return(matrix(exp(logHazSamp),nr=n, ncol = 6, byrow = T))
}
tst <- logHazPrior(1000)
apply(tst, 2, function(x) quantile(x, c(.025,.5,.975)))
exp(logHazMeans)
spars[,'50%']

## Flat uniform prior from 1 to 3
hetGenSDPrior <- function(n) runif(n, 1, 3)

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
abcSimSumStat <- function(parms=simParmSamp(1), maxN=10000, seed=1, nc=12, browse=F,
                          condRakai = FALSE, RakSamp = c(inc = 23, prev = 161, late = 51)) {
    tempFileNm <- with(parms, psrun(maxN = maxN, jobnum = seed, pars = parms[hazs], save.new = T, return.ts = T,
                                    acute.sc = acute.sc, dur.ac = dur.ac, het.gen.sd = het.gen.sd, browse = F, nc = nc))
    load(tempFileNm)
    cohsim <- rak.coh.fxn(output, #ts.ap = output$ts, dat = output$evout, dpars = output$rakpars,
                          ltf.prob=0.0287682072451781, 
                          rr.ltf.ff=1.5, rr.ltf.mm=1.5, rr.ltf.hh=1, rr.ltf.d=1, rr.inc.sdc=1.5,
                          verbose = F, browse = F)
                                        #    rm(output); gc() ## free up memory
    rcohsim <- rak.wawer(rak.coh = cohsim, excl.extram=T, het.gen.sd = with(parms,het.gen.sd),
                         cov.mods = F, verbose = F, fit.Pois=F, browse=F)
    ## Condition on Rakai retrospecive cohort sample sizes
    if(condRakai) rcohsim <- within(rcohsim, {
        inc.wh <- which(rakll$phase=='inc')
        prev.wh <- which(rakll$phase=='prev')
        late.wh <- which(rakll$phase=='late')
        resamp <- c(sample(inc.wh, RakSamp['inc'], replace=T), sample(prev.wh, RakSamp['prev'], replace=T), 
                    sample(late.wh, RakSamp['late'], replace=T))
        rakll <- rakll[resamp,]})
    wtabSim <- sbmod.to.wdat(rcohsim$rakll, excl.by.err = T, browse=F)
    return(list(wtabSim=wtabSim, parms=parms))
}


