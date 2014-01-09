rm(list=ls())
library(mnormt);library(coda); library(stats4) #library(lme4)
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/') # setwd
source('PlotFunctions.R')                    # load functions to collect & plot results
## #load('results/Uganda-46900-321.Rdata')
## save.image('results/Uganda-46900-321.Rdata')
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions
##source('/home1/02413/sbellan/SDPSimulations/SimulationFunctions.R')                   # load old simulation functions

nc <- 12
interv <- 10
max.vis <- floor(40/interv) + 1

dpars <- c(acute.sc = 15, late.sc = 4, bp = .007,
             dur.ac = 2, dur.lt = 6, dur.aids = 4)
ldpars <- log(dpars)

nn <- 2000 # crashes for > 4000, not sure why
sim <- rbind(holl.mod(nn,nn,nn,dpars=dpars, verbose=F), holl.mod(nn,nn,nn,dpars=dpars, verbose=F))
xtabs(~inf + kk + phase, sim)

holl.lik(log(dpars), sim, excl.by.err = F, verbose = T, browse = F)

nc <- 12
niter <- 2*1000
nburn <- 500
sd.props <- c(.05, .1, .01, .05, .1, .05)*1.2
names(sd.props) <- names(ldpars)
seed.bump <- 0
d.out <- mclapply((seed.bump + 1:nc), holl.wrp, jit = .5,
                  sd.props = sd.props, force.inits = ldpars, rakdat = sim, excl.by.err = F,
                  multiv = F, covar = NULL, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = niter, nthin = 1, nburn = nburn, browse=F)
ndirnm <- 'results'
posts <- procpost(d.out)$posts
dmcmc <- procpost(d.out)$mcmc
hollsbpairs(posts, truepars = ldpars, show.lines = T, ## plot posterior correlations after adaptive phase
        file.nm = file.path(ndirnm,"Holl posterior pairs after adaptive phase"), width = 12, height = 12,
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
exposts <- as.data.frame(exp(posts))
exposts$atr.month.ac <- (exposts$acute.sc-1)*exposts$dur.ac
edpars <- dpars
edpars <- c(edpars, atr.month.ac = (edpars['acute.sc']-1)*edpars['dur.ac'])
hollsbpairs(exposts, truepars = edpars, show.lines = T, ## plot posterior correlations after adaptive phase
        file.nm = file.path(ndirnm,"exp Holl posterior pairs after adaptive phase"), width = 12, height = 12,
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
apply(exposts, 2, range)
apply(exposts, 2, function(x) signif(quantile(x, c(.025, .5, .975)),3))
edpars
sigma.ad <- cov.wt(posts)$cov ## posterior covariance matrix, then plot what proposal distr is gonna look like
sigma <- sigma.ad


####################################################################################################
## fitting phase with multivariate block sampling
niter <- 4*1000
nburn <- 1500
f.out <- mclapply((seed.bump + 1:nc), holl.wrp, force.inits = ldpars, jit = 1,
                  sd.props = sd.props,  rakdat = sim, excl.by.err = F,
                  multiv = T, covar = sigma, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = niter, nthin = 1, nburn = nburn, browse=F)
fposts <- procpost(f.out)$posts
fmcmc <- procpost(f.out)$mcmc.out
pdf('mcmc summ.pdf'); plot(fmcmc); dev.off()
ord <- c(3,1,4,2,5,6)
hollsbpairs(fposts[,ord], truepars = ldpars[ord], show.lines = T, ## plot posterior correlations after fit phase
        file.nm = file.path(ndirnm,"Holl posterior pairs after fit phase"), width = 12, height = 12,
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
exfposts <- as.data.frame(exp(fposts))
exfposts$atr.month.ac <- (exfposts$acute.sc-1)*exfposts$dur.ac
edpars <- dpars
edpars <- c(edpars, atr.month.ac = (edpars['acute.sc']-1)*edpars['dur.ac'])
ord <- c(3,7,1,4,2,5,6)
hollsbpairs(exfposts[,ord], truepars = edpars[ord], show.lines = T, ## plot posterior correlations after fit phase
        file.nm = file.path(ndirnm,"exp Holl posterior pairs after fit phase"), width = 12, height = 12,
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
apply(exfposts, 2, range)
apply(exfposts, 2, function(x) signif(quantile(x, c(.025, .5, .975)),3))
edpars

fsigma.ad <- cov.wt(fposts)$cov ## posterior covariance matrix, then plot what proposal distr is gonna look like
fsigma <- fsigma.ad
save.image('results/140108.Rdata')

## #ltf.prob <- NA
## rak.coh <- rak.coh.fxn(ts=output$ts, dat = output$evout, interv=interv, max.vis=max.vis, 
##                        ltf.prob = 0, rr.ltf.ff = 4.3, rr.ltf.mm = 2, rr.ltf.hh = 1, rr.ltf.d = 0,
##                        start.rak=start.rak,end.rak=end.rak,browse = F)
## nrow(rak.coh$dat)

## rkout <- rak.wawer(rak.coh = rak.coh, verbose = F, browse=F)
## rkout[[1]]; rkout[[2]]
## rkt <- rkout$rakll
## xtabs(~ phase + excl.by.err, rkt)
## xtabs(inf.trunc ~ phase, rkt)
## xtabs(pm.trunc ~ phase + inf.trunc + excl.by.err, rkt)

## holl.lik(ldpars, rakll = rkt, excl.by.err = F, verbose = T, browse = F)

## system.time(
##             for(ac in seq(1,8,by=.5)) {
##               print(ac)
##               dpars.t <- dpars
##               dpars.t[1] <- ac
##               ldpars.t <- log(dpars.t)
##               print(holl.lik(ldpars.t, rkt, excl.by.err = F, verbose = F, browse = F))
##             }
##             )


## ldpars2 <- jitter(ldpars, a = 2)
## opt <- optim(par = ldpars2, fn = holl.lik, method='SANN', rakll = sim, verbose=F, control=list(maxit=100, trace=5))
## opt <- optim(par = opt$par, fn = holl.lik, method='Nelder-Mead', rakll = sim, verbose=F, control=list(maxit=400, trace=3))
## opt$conv
## exp(opt$par)
## dpars
## arr(ldpars)
## arr(opt$par)

## opt <- optim(par = ldpars, fn = holl.lik, method='SANN', rakll = sim, verbose=F, control=list(maxit=200, trace=5))

## pfit <- c('acute.sc','dur.ac','bp')
## to.fit <- names(dpars) %in% pfit


## mlo <- mle(minuslogl = holl.lik.mle, start = as.list(ldpars[to.fit]),
##            fixed = append(as.list(ldpars[!to.fit]), list(rakll = sim, excl.by.err = F, verbose = F, browse = F)),
##            method = 'Nelder-Mead', control=list(maxit=300, trace=3))
## summary(mlo)
## exp(coef(mlo)[pfit])

## cis <- confint(mlo)
  


####################################################################################################

## ######################################################################
## ## look at age-haz & partnership-haz relationship without any
## ## controlling for covariates
## ##
## ## need to decide how to include age, let's do age in years.
## ## males
## ts.mage <- matrix(rep(dat$mage, each = nrow(ts)), nr = nrow(ts), nc = nrow(dat))
## ts.mage <- apply(ts.mage, 2, function(x) {x - (length(x)-1):0})
## ts.mage[which(is.na(ts), arr.ind=T)] <- NA
## ## females
## ts.fage <- matrix(rep(dat$fage, each = nrow(ts)), nr = nrow(ts), nc = nrow(dat))
## ts.fage <- apply(ts.fage, 2, function(x) {x - (length(x)-1):0})
## ts.fage[which(is.na(ts), arr.ind=T)] <- NA
## ## partnership duration in months
## ts.pdur <- matrix(rep(dat$mardur.mon, each = nrow(ts)), nr = nrow(ts), nc = nrow(dat))
## ts.pdur <- apply(ts.pdur, 2, function(x) {x - (length(x)-1):0})
## ts.pdur[which(is.na(ts), arr.ind=T)] <- NA
## ## round to years
## ts.mage <- (floor(ts.mage/12))#*12
## ts.fage <- (floor(ts.fage/12))#*12
## ts.pdur <- (floor(ts.pdur/12))#*12

##   ## build line list data
## for(ph in phs[2:5])
##   {
##     ## temporary stuff
##     selt <- get(paste('sel.',ph,sep=''))
##     pmt <-  get(paste('pm.',ph,sep=''))
##     inft <- get(paste(ph,'.inf',sep=''))
##     temp <- data.frame(uid = which(selt), inf = 0, mon = pmt[selt],
##                        pdur = NA, age = NA, sex = gend[selt], ph = ph) # gender of uninfected
##     temp$inf[temp$uid %in% inft] <- 1   # who was infected
##     fmonth <- apply(ts[,selt], 2, function(x) {min(which(grepl(ph, x)))}) # first month in that phase
##     temp$pdur <- ts.pdur[cbind(fmonth, which(selt))] # partnership duration in 1st month of phase
##     temp$age[temp$sex=='f'] <- ts.fage[cbind(fmonth[temp$sex=='f'], which(selt)[temp$sex=='f'])] # age of uninfected in 1st month of phase: female
##     temp$age[temp$sex=='m'] <- ts.fage[cbind(fmonth[temp$sex=='m'], which(selt)[temp$sex=='m'])] # male
##     temp$yr <- temp$mon/12
##     assign(paste('ldat.',ph,sep=''), temp)
##   }

## ldat <- rbind(ldat.ac, ldat.ch, ldat.lt, ldat.ad)
## ldat <- ldat[order(ldat$uid),]
## ldat$ph <- factor(ldat$ph, levels = c('ch','ac','lt','ad'))
## xtabs(inf ~ ph,ldat) / xtabs(yr ~ ph, ldat) * 100


## nn <- 10000
## enn <- c(.9*nn,.1*nn)
## x <- rep(c(1,15)*.005, enn)
## x.ph <- factor(rep(c('ch','ac'), enn),levels=c('ch','ac'))
## oo <- runif(nn, 1/12, 2)
## y <- rpois(nn, x*oo)
## exp(coef(glm(y ~ x.ph + offset(oo), family = poisson)))[2]
## z <- xtabs(y~x.ph) / xtabs(oo~x.ph)
## z/z[1]

## 1-(1-.01)^10
## 1-exp(-.01*10)
## xs <- rnorm(1000, sd = 3.25)
## clxs <- exp(xs)
## sd(clxs)

## cldat <- ldat                       # change ldat to see what's wronog
## ## something about being unbalanced?
## cldat$ph[1:nrow(cldat) %% 2 == 1] <- 'ac'
## cldat$ph[!1:nrow(cldat) %% 2 == 1] <- 'ch'
## cldat$yr[cldat$ph=='ac'] <- ceiling(runif(sum(cldat$ph=='ac'), 0,2))/12
## cldat$yr[!cldat$ph=='ac'] <- ceiling(runif(sum(cldat$ph!='ac'), 0,23))/12
## cldat$inf[cldat$ph=='ac'] <- rpois(sum(cldat$ph=='ac'), 90/1200 * cldat$yr[cldat$ph=='ac'])
## cldat$inf[!cldat$ph=='ac'] <- rpois(sum(!cldat$ph=='ac'), 10/1200 * cldat$yr[!ldat$ph=='ac'])
## mod1 <- glm(inf ~ offset(yr) + ph, family = poisson, data = cldat)
## summary(mod1)
## exp(coef(mod1))['phac']
## z <- xtabs(inf~ph,cldat) / xtabs(yr~ph,cldat)
## z/z[1]

## mod1 <- glm(inf ~ offset(yr) + ph, family = poisson, data = ldat)
## summary(mod1)
## exp(coef(mod1))['phac']
## z <- xtabs(inf~ph,ldat) / xtabs(yr~ph,ldat)
## z/z[1]

## mod1 <- glm(inf ~ offset(yr) + pdur + age + sex + ph, family = poisson, data = ldat)
## mod1 <- glm(inf ~ offset(yr) + ph, family = poisson, data = ldat)
## summary(mod1)
## exp(coef(mod1))
## mod2 <- lmer(inf ~ offset(yr) + pdur + age + sex + ph + (1 | uid), family = poisson, data = ldat)
## mod2 <- lmer(inf ~ offset(yr) +  ph + (1 | uid), family = poisson(link=log), data = ldat)
## summary(mod2)
## exp(fixef(mod2))

## ## left off here, trying to create line list data frame for poisson reg by age/pdur
## ac.ind <- which(grepl(ts[,sel.sdc]
## for(ii in 1:sum(sel.sdc))
##     {
##         ind <- which(sel.sdc)[ii]

##       }
        

## load('ac5het0.Rdata')
## rak.fxn(ac5het0$evout, ac5het0$ts, browse = F)

## ## add poisson regression looking at partner age & partnership duration
    
