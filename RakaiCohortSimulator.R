rm(list=ls())
library(mnormt);library(coda); library(stats4) #library(lme4)
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/') # setwd
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions


simul <- T ## Simulate or analyze Wawer et al. Table 1 data?
excl.by.err <- F ## exclude couples by error?
seed.bump <- 0 ## for changing simulation seed
nc <- 12 ## number of cores
interv <- 10 ## interval between survey visits
max.vis <- floor(40/interv) + 1 ## max visits
aniter <- 3*1000 ## adaptive MCMC iterations
anburn <- 1000 ## adaptive MCMC burn-in
niter <- 5*1000 ## MCMC iterations
nburn <- 1*1000 ## MCMC burn-in
## simulation paramters
dpars <- c(acute.sc = 26, late.sc = 15, bp = 10.6/1200, dur.ac = 2.9, dur.lt = 9, dur.aids = 10)
ldpars <- log(dpars)

if(simul) { ## simulate data
    nn <- 500 # crashes for > 4000, not sure why
    sim <- rbind(holl.mod(nn,nn,nn,dpars=dpars, verbose=F), holl.mod(nn,nn,nn,dpars=dpars, verbose=F))
    wtab.sim <- sbmod.to.wdat(sim, excl.by.err = excl.by.err)
    wtab.do <- wtab.sim
    print(wtab.do)
    nm <- ifelse(excl.by.err,'sim exclErr', 'sim no exclErr')
    sd.props <- c(.05, .1, .01, .05, .05, .05)*2.5
    if(excl.by.err) {
        ndirnm <- file.path('results','HollingsworthAn','SimExclbyErr')
    }else{
        ndirnm <- file.path('results','HollingsworthAn','SimNOExclbyErr')
    }
}else{ ## real data
    if(excl.by.err) {
        wtab.do <- wtab.rl
        ndirnm <- file.path('results','HollingsworthAn','RealExclbyErr')
    }else{
        wtab.do <- wtab.rl.no.err
        ndirnm <- file.path('results','HollingsworthAn','RealNOExclbyErr')
    }
    nm <- ifelse(excl.by.err,'real exclErr', 'real no exclErr')
    sd.props <- c(.1, .3, .01, .2, .5, .3)*1.2    
}
names(sd.props) <- names(ldpars)
if(!file.exists(ndirnm)) dir.create(ndirnm)
print(ndirnm)
print(wtab.do)
holl.lik(log(dpars), wtab.do, verbose = T, browse = F)

## Initialize a covariance matrix with -.9 correlation between acute/late scalar & duration.
sigma <- matrix(0, nr = length(sd.props), nc = length(sd.props))
rownames(sigma) <- names(ldpars)
colnames(sigma) <- names(ldpars)
diag(sigma) <- sd.props^2
sigma[1,4] <- -sd.props[1]*sd.props[4]*.9
sigma[4,1] <- -sd.props[1]*sd.props[4]*.9
sigma[2,5] <- -sd.props[2]*sd.props[5]*.9
sigma[5,2] <- -sd.props[2]*sd.props[5]*.9

d.out <- mclapply((seed.bump + 1:nc), holl.wrp, jit = .2,
                  sd.props = sd.props, force.inits = ldpars, wtab = wtab.do, 
                  multiv = T, covar = sigma, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = aniter, nthin = 1, nburn = 1, browse=F)
pout <- procpost(d.out, ldpars = ldpars, dirnm=ndirnm, nm = 'adapt1'); print(pout$gel)
sigma1 <- pout$sigma
## Adaptive phase #2, this time with sigma from previous run
d.out2 <- mclapply((seed.bump + 1:nc), holl.wrp, force.inits = ldpars, jit = .2,
                  sd.props = sd.props,  wtab = wtab.do,
                  multiv = T, covar = sigma1, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = aniter, nthin = 1, nburn = anburn, browse=F)
pout2 <- procpost(d.out2, ldpars = ldpars, dirnm=ndirnm, nm = 'adapt2'); print(pout2$gel)
sigma2 <- pout2$sigma

####################################################################################################
## fitting phase with multivariate block sampling
f.out <- mclapply((seed.bump + 1:nc), holl.wrp, force.inits = ldpars, jit = .2,
                  sd.props = sd.props,  wtab = wtab.do,
                  multiv = T, covar = sigma2, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = niter, nthin = 1, nburn = nburn, browse=F)
fout <- procpost(f.out, ldpars = ldpars, dirnm=ndirnm, nm = 'final'); print(fout$gel)
pdf(file.path(ndirnm,'mcmc summ.pdf')); plot(fmcmc); dev.off()
outtab <- fout$outtab
rownames(outtab)[4] <- 'Hollingsworth Estimates'
write.csv(outtab, file.path(ndirnm,'output table.csv'))
outtab

save.image(file.path(ndirnm,'workspace.Rdata'))
#load(file.path(ndirnm,'workspace.Rdata'))
