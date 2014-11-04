####################################################################################################
## Create simulations for first ABC sample
####################################################################################################
vlDir <- file.path('FiguresAndTables','VL Profile')
load("data files/ds.nm.all.Rdata")                      ## DHS country names
load('data files/pars.arr.ac.Rdata')   ## load hazards for each transmission route as fit to DHS data (in.arr[,,2])
load(file=file.path(vlDir, 'ehms.vl.Rdata'))
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') ##  transmission coefficient names, for convenience
nc <- 12                                       ## cores per simulation

## Prior on average hazard.
## Transmission rate from Wawer et al. 2005
lambdaMean <- 36/48525 * 10 ## monthly hazard = infections/(coital acts) * (coital acts)/month
## Want a prior around this that allows substantial uncertainty because we expect lambdaMean to be
## much bigger for more heterogeneity (since prevalent couples will represent the lowest risk
## couples). Use scl below to increase the prior so we sample reasonably.
lambdaPrior <- function(n, scl = 1, sd = 1) exp(rnorm(n, log(lambdaMean*scl), sd = sd))
exp(qnorm(c(.025,.975), mean = log(lambdaMean), sd = 1)) ## seems broad enough

hetGenSDPrior <- function(n) runif(n, 1, 3)

## Get RHacute prior centered on estimates from viral loads but with way more uncertainty
log(ehms.vl['med'])-log(ehms.vl[c('lci','uci')]) ## How far away are bounds on a log scale? about .5
sd_RHacutePrior <- .5/1.96*4 ## add a bunch of extra uncertainty since this is just VL
exp(qnorm(c(.025,.975), log(ehms.vl['med']), sd = sd_RHacutePrior)) ## 95% CI bounds: .75 - 41X as infectious
rhAcutePrior <- function(n) exp(rnorm(n, log(ehms.vl['med']), sd = sd_RHacutePrior))

durAcutePrior <- function(n) runif(n, .5, 8) ## uniform from 2 weeks to 8 months
