rm(list=ls())
library(mnormt);library(coda); library(stats4) #library(lme4)
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/') # setwd
source('PlotFunctions.R')                    # load functions to collect & plot results
## #load('results/Uganda-46900-321.Rdata')
## save.image('results/Uganda-46900-321.Rdata')
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions
##source('/home1/02413/sbellan/SDPSimulations/SimulationFunctions.R')                   # load old simulation functions

simul <- F
excl.by.err <- F
seed.bump <- 0
nc <- 12
interv <- 10
max.vis <- floor(40/interv) + 1
aniter <- 3*1000
anburn <- 1000
niter <- 8*1000
nburn <- 2*1000
dpars <- c(acute.sc = 26, late.sc = 7, bp = 10.6/1200, dur.ac = 2.9, dur.lt = 9, dur.aids = 10)
ldpars <- log(dpars)

if(simul) { ## simulate data
    nn <- 2000 # crashes for > 4000, not sure why
    sim <- rbind(holl.mod(nn,nn,nn,dpars=dpars, verbose=F), holl.mod(nn,nn,nn,dpars=dpars, verbose=F))
    wtab.sim <- hmod.to.wdat(sim, excl.by.err = excl.by.err)
    wtab.do <- wtab.sim
    nm <- ifelse(excl.by.err,'sim exclErr', 'sim no exclErr')
    sd.props <- c(.05, .1, .01, .05, .1, .05)*1.2
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

d.out <- mclapply((seed.bump + 1:nc), holl.wrp, jit = .5,
                  sd.props = sd.props, force.inits = ldpars, wtab = wtab.do, 
                  multiv = F, covar = NULL, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = aniter, nthin = 1, nburn = anburn, browse=F)
posts <- procpost(d.out)$posts
dmcmc <- procpost(d.out)$mcmc
parnames <- c('bp','acute.sc','dur.ac','late.sc','dur.lt','dur.aids')
hollsbpairs(posts[,parnames], truepars = ldpars[parnames], show.lines = T, ## plot posterior correlations after adaptive phase
        file.nm = file.path(ndirnm, paste0("adapt posterior (log)",nm)), width = 12, height = 12,
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
exposts <- posts
exposts[names(exposts)!='nll'] <- exp(posts[names(posts)!='nll'])
edpars <- dpars
edpars <- c(edpars, atr.month.ac = as.numeric((edpars['acute.sc']-1)*edpars['dur.ac']),
                    atr.month.lt = as.numeric((edpars['late.sc']-1)*edpars['dur.lt']))
tracenames <- c('bp','atr.month.ac','acute.sc','dur.ac','atr.month.lt','late.sc','dur.lt','dur.aids')
hollsbpairs(exposts[,tracenames], truepars = edpars[tracenames], show.lines = T, ## plot posterior correlations after adaptive phase
        file.nm = file.path(ndirnm, paste0("adapt posterior",nm)), width = 12, height = 12,
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
## apply(exposts, 2, range)
## apply(exposts, 2, function(x) signif(quantile(x, c(.025, .5, .975),na.rm=T),3))
## edpars
sigma <- cov.wt(posts[,names(ldpars)])$cov ## posterior covariance matrix, then plot what proposal distr is gonna look like
## Adaptive phase #2, this time with multivariate block sampling
niter <- 3*1000
nburn <- 1000
d.out2 <- mclapply((seed.bump + 1:nc), holl.wrp, force.inits = ldpars, jit = .5,
                  sd.props = sd.props,  wtab = wtab.do,
                  multiv = T, covar = sigma, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = niter, nthin = 1, nburn = nburn, browse=F)
posts2 <- procpost(d.out2, nc = 12)$posts
sigma <- cov.wt(posts2[,names(ldpars)])$cov ## update covariance matrix

####################################################################################################
## fitting phase with multivariate block sampling
niter <- 6*1000
nburn <- 1500
f.out <- mclapply((seed.bump + 1:nc), holl.wrp, force.inits = ldpars, jit = .5,
                  sd.props = sd.props,  wtab = wtab.do,
                  multiv = T, covar = sigma, 
                  verbose = T, verbose2 = F, tell = 100, 
                  niter = niter, nthin = 1, nburn = nburn, browse=F)
fposts <- procpost(f.out, nc = 12)$posts
fmcmc <- procpost(f.out, nc = 12)$mcmc.out
pdf(file.path(ndirnm,'mcmc summ.pdf')); plot(fmcmc); dev.off()
hollsbpairs(fposts[,parnames], truepars = ldpars[parnames], show.lines = T, ## plot posterior correlations after fit phase
        file.nm = file.path(ndirnm, paste0("final posterior (log)",nm)), width = 12, height = 12,
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
expfposts <- posts
expfposts[names(expfposts)!='nll'] <- exp(posts[names(posts)!='nll'])
edpars <- dpars
edpars <- c(edpars, atr.month.ac = as.numeric((edpars['acute.sc']-1)*edpars['dur.ac']),
                    atr.month.lt = as.numeric((edpars['late.sc']-1)*edpars['dur.lt']))
hollsbpairs(expfposts[,tracenames], truepars = edpars[tracenames], show.lines = T, ## plot posterior correlations after fit phase
        file.nm = file.path(ndirnm, paste0("final posterior",nm)), width = 12, height = 12, range0=c(2:length(edpars)),
        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
outtab <- rbind(apply(expfposts[tracenames], 2, function(x) signif(quantile(x, c(.025, .5, .975),na.rm=T),3)),
                edpars[tracenames])
rownames(outtab)[4] <- 'Hollingsworth Estimates'
write.csv(outtab, file.path(ndirnm,'output table.csv'))
outtab
gelman.diag(fmcmc)

#sigma <- cov.wt(fposts[,names(ldpars)])$cov ## update covariance matrix


f1 <- apply(expfposts, 2, median)
f1['acute.sc']*f1['bp']*12*100
f1['bp']*12*100
f1['late.sc']*f1['bp']*12*100

save.image(file.path(ndirnm,'workspace.Rdata'))
