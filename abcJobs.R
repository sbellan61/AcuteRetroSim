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

n <- 10e4
simParms <- data.frame(lambda = lambdaPrior(n), acute.sc = rhAcutePrior(n), dur.ac = durAcutePrior(n), het.gen.sd = hetGenSDPrior(n),
                       dur.lt = 10, dur.aids = 10, late.sc = 5, aids.sc=0,## Fix these since they don't affect acute estimates much
                       het.gen.cor = 0, ## going with no intra-partner correlation in risk (conservative to this assumption)
                       group = which(ds.nm=='Uganda'), maxN = 10e5, each = 200) 
simParms$jobnum <- 1:nrow(simParms)
head(simParms)
nn <- nrow(simParms)
to.do <- 1:nrow(simParms) # which rows in simParms to do (default is all)

outdir <- file.path('results','RakABC1')
if(!file.exists(outdir))      dir.create(outdir) # create directory if necessary
if(!file.exists(file.path(outdir,'routs')))      dir.create(file.path(outdir, 'routs')) # setup directory to store Rout files

## #################################################################### 
sink("RakABC1.txt") ## create a control file to send to the cluster
for(ii in to.do) {
    cmd <- with(simParms, paste("R CMD BATCH '--args jobnum=", jobnum[ii], " simj=", ii, " outdir=\"", outdir, "\"", " nc=", nc,
                                " group.ind=", group[ii], " lambda=", lambda[ii],
                                " acute.sc=", acute.sc[ii], " late.sc=", late.sc[ii]," aids.sc=", aids.sc[ii], # acute phase varying throughout loop
                                " dur.ac=", dur.ac[ii], " dur.lt=", dur.lt[ii], " dur.aids=", dur.aids[ii],
                                " het.gen=T het.gen.sd=", het.gen.sd[ii], " het.gen.cor=", het.gen.cor[ii],
                                " maxN=", maxN[ii], " seed=1 each=", each[ii],
                                " return.ts=TRUE", ## cohort dates, return a time series as output
                                " one.couple=F", ## debugging with one couple replicated a bunch of times
                                ## interview date, call simulation script and specify output
                                " tint=100*12' SimulationStarter.R ", file.path(outdir, "routs", paste0(ds.nm[group[ii]], ii, ".Rout")), sep='') )
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
save(simParms, file = file.path(outdir,'simParms.Rdata')) # these are country-acute phase specific simParms


