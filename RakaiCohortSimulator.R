rm(list=ls())
library(stats4) #library(lme4)
#setwd('/home1/02413/sbellan/Rakai/SDPSimulations/') # setwd
source('PlotFunctions.R')                    # load functions to collect & plot results
load('results/Uganda-46900-321.Rdata')
#save.image('results/Uganda-46900-321.Rdata')
source('SimulationFunctions.R')                   # load simulation functions
##source('/home1/02413/sbellan/SDPSimulations/SimulationFunctions.R')                   # load old simulation functions

interv <- 10
max.vis <- floor(40/interv) + 1

#ltf.prob <- NA
rak.coh <- rak.coh.fxn(ts=output$ts, dat = output$evout, interv=interv, max.vis=max.vis, 
                       ltf.prob = 0, rr.ltf.ff = 4.3, rr.ltf.mm = 2, rr.ltf.hh = 1, rr.ltf.d = 0,
                       start.rak=start.rak,end.rak=end.rak,browse = F)
nrow(rak.coh$dat)

rkout <- rak.wawer(rak.coh = rak.coh, verbose = F, browse=F)
rkout[[1]]; rkout[[2]]
rkt <- rkout$rakll
xtabs(~ phase + excl.by.err, rkt)
xtabs(inf.trunc ~ phase, rkt)
xtabs(pm.trunc ~ phase + inf.trunc + excl.by.err, rkt)

dpars <- c(acute.sc = 7, late.sc = 1, aids.sc = 1, bp = sqrt(prod(spars[5:6])),
             dur.ac = 2, dur.lt = 9, dur.aids = 10)
ldpars <- log(dpars)

holl.lik(ldpars, rakll = rkt, excl.by.err = F, verbose = T, browse = F)

system.time(
            for(ac in seq(1,8,by=.5)) {
              print(ac)
              dpars.t <- dpars
              dpars.t[1] <- ac
              ldpars.t <- log(dpars.t)
              print(holl.lik(ldpars.t, rkt, excl.by.err = F, verbose = F, browse = F))
            }
            )

nn <- 1000
sim <- holl.mod(nn,nn,nn,dpars=dpars, verbose=F)
## dpars['late.sc'] <- 10
## dpars['aids.sc'] <- 0
xtabs(~inf + kk + phase, sim)
xtabs(~inf + kk + phase, rkt)

####################################################################################################
## left off here, having problems calculating likelihood for simulated data
holl.lik(log(dpars), sim, excl.by.err = F, verbose = F, browse = F)

ldpars2 <- jitter(ldpars, a = 2)
opt <- optim(par = ldpars2, fn = holl.lik, method='SANN', rakll = sim, verbose=F, control=list(maxit=100, trace=5))
opt <- optim(par = opt$par, fn = holl.lik, method='Nelder-Mead', rakll = sim, verbose=F, control=list(maxit=400, trace=3))
opt$conv
exp(opt$par)
dpars
arr(ldpars)
arr(opt$par)

opt <- optim(par = ldpars, fn = holl.lik, method='SANN', rakll = sim, verbose=F, control=list(maxit=200, trace=5))

pfit <- c('acute.sc','dur.ac','bp')
to.fit <- names(dpars) %in% pfit


mlo <- mle(minuslogl = holl.lik.mle, start = as.list(ldpars[to.fit]),
           fixed = append(as.list(ldpars[!to.fit]), list(rakll = sim, excl.by.err = F, verbose = F, browse = F)),
           method = 'Nelder-Mead', control=list(maxit=300, trace=3))
summary(mlo)
exp(coef(mlo)[pfit])

cis <- confint(mlo)
  


####################################################################################################

######################################################################
## look at age-haz & partnership-haz relationship without any
## controlling for covariates
##
## need to decide how to include age, let's do age in years.
## males
ts.mage <- matrix(rep(dat$mage, each = nrow(ts)), nr = nrow(ts), nc = nrow(dat))
ts.mage <- apply(ts.mage, 2, function(x) {x - (length(x)-1):0})
ts.mage[which(is.na(ts), arr.ind=T)] <- NA
## females
ts.fage <- matrix(rep(dat$fage, each = nrow(ts)), nr = nrow(ts), nc = nrow(dat))
ts.fage <- apply(ts.fage, 2, function(x) {x - (length(x)-1):0})
ts.fage[which(is.na(ts), arr.ind=T)] <- NA
## partnership duration in months
ts.pdur <- matrix(rep(dat$mardur.mon, each = nrow(ts)), nr = nrow(ts), nc = nrow(dat))
ts.pdur <- apply(ts.pdur, 2, function(x) {x - (length(x)-1):0})
ts.pdur[which(is.na(ts), arr.ind=T)] <- NA
## round to years
ts.mage <- (floor(ts.mage/12))#*12
ts.fage <- (floor(ts.fage/12))#*12
ts.pdur <- (floor(ts.pdur/12))#*12

  ## build line list data
for(ph in phs[2:5])
  {
    ## temporary stuff
    selt <- get(paste('sel.',ph,sep=''))
    pmt <-  get(paste('pm.',ph,sep=''))
    inft <- get(paste(ph,'.inf',sep=''))
    temp <- data.frame(uid = which(selt), inf = 0, mon = pmt[selt],
                       pdur = NA, age = NA, sex = gend[selt], ph = ph) # gender of uninfected
    temp$inf[temp$uid %in% inft] <- 1   # who was infected
    fmonth <- apply(ts[,selt], 2, function(x) {min(which(grepl(ph, x)))}) # first month in that phase
    temp$pdur <- ts.pdur[cbind(fmonth, which(selt))] # partnership duration in 1st month of phase
    temp$age[temp$sex=='f'] <- ts.fage[cbind(fmonth[temp$sex=='f'], which(selt)[temp$sex=='f'])] # age of uninfected in 1st month of phase: female
    temp$age[temp$sex=='m'] <- ts.fage[cbind(fmonth[temp$sex=='m'], which(selt)[temp$sex=='m'])] # male
    temp$yr <- temp$mon/12
    assign(paste('ldat.',ph,sep=''), temp)
  }

ldat <- rbind(ldat.ac, ldat.ch, ldat.lt, ldat.ad)
ldat <- ldat[order(ldat$uid),]
ldat$ph <- factor(ldat$ph, levels = c('ch','ac','lt','ad'))
xtabs(inf ~ ph,ldat) / xtabs(yr ~ ph, ldat) * 100


nn <- 10000
enn <- c(.9*nn,.1*nn)
x <- rep(c(1,15)*.005, enn)
x.ph <- factor(rep(c('ch','ac'), enn),levels=c('ch','ac'))
oo <- runif(nn, 1/12, 2)
y <- rpois(nn, x*oo)
exp(coef(glm(y ~ x.ph + offset(oo), family = poisson)))[2]
z <- xtabs(y~x.ph) / xtabs(oo~x.ph)
z/z[1]

1-(1-.01)^10
1-exp(-.01*10)
xs <- rnorm(1000, sd = 3.25)
clxs <- exp(xs)
sd(clxs)

cldat <- ldat                       # change ldat to see what's wronog
## something about being unbalanced?
cldat$ph[1:nrow(cldat) %% 2 == 1] <- 'ac'
cldat$ph[!1:nrow(cldat) %% 2 == 1] <- 'ch'
cldat$yr[cldat$ph=='ac'] <- ceiling(runif(sum(cldat$ph=='ac'), 0,2))/12
cldat$yr[!cldat$ph=='ac'] <- ceiling(runif(sum(cldat$ph!='ac'), 0,23))/12
cldat$inf[cldat$ph=='ac'] <- rpois(sum(cldat$ph=='ac'), 90/1200 * cldat$yr[cldat$ph=='ac'])
cldat$inf[!cldat$ph=='ac'] <- rpois(sum(!cldat$ph=='ac'), 10/1200 * cldat$yr[!ldat$ph=='ac'])
mod1 <- glm(inf ~ offset(yr) + ph, family = poisson, data = cldat)
summary(mod1)
exp(coef(mod1))['phac']
z <- xtabs(inf~ph,cldat) / xtabs(yr~ph,cldat)
z/z[1]

mod1 <- glm(inf ~ offset(yr) + ph, family = poisson, data = ldat)
summary(mod1)
exp(coef(mod1))['phac']
z <- xtabs(inf~ph,ldat) / xtabs(yr~ph,ldat)
z/z[1]

mod1 <- glm(inf ~ offset(yr) + pdur + age + sex + ph, family = poisson, data = ldat)
mod1 <- glm(inf ~ offset(yr) + ph, family = poisson, data = ldat)
summary(mod1)
exp(coef(mod1))
mod2 <- lmer(inf ~ offset(yr) + pdur + age + sex + ph + (1 | uid), family = poisson, data = ldat)
mod2 <- lmer(inf ~ offset(yr) +  ph + (1 | uid), family = poisson(link=log), data = ldat)
summary(mod2)
exp(fixef(mod2))

## left off here, trying to create line list data frame for poisson reg by age/pdur
ac.ind <- which(grepl(ts[,sel.sdc]
for(ii in 1:sum(sel.sdc))
    {
        ind <- which(sel.sdc)[ii]

      }
        

load('ac5het0.Rdata')
rak.fxn(ac5het0$evout, ac5het0$ts, browse = F)

## make sure that all partner infections within this period add up to # of inf's counted
sum(grepl('hh', ts[endmon,]) & (ev$mcoi=='p' | ev$fcoi=='p'))
sum(acute.n + chron.n + late.n + aids.n)

ww <-c(acute.w,chron.w,late.w,aids.w)
ww <- ww[order(ww)]

wwa <- which(grepl('hh', ts[endmon,]) & (ev$mcoi=='p' | ev$fcoi=='p'))

bz <- wwa[which(!wwa %in% ww)]
bzl <- 1:nrow(ev) %in% bz

ts[sttmon,bz]
ts[endmon,bz]

ev[bz[1],]
ev$fcoi.phase[bz]
ev$mcoi.phase[bz]

ev$mcoi[bzl & inc.inf]
ev$mcoi[bzl & chron.inf]
ev$mcoi.phase[bzl & late.inf]
ev$mcoi.phase[bzl & aids.inf]

ev[which(bzl&late.inf)[1],]


ts[sttmon,acute.ww]
ts[endmon,acute.ww]

## add poisson regression looking at partner age & partnership duration
    
