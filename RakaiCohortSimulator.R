library(mnormt);library(coda); library(stats4);library(plyr)#library(lme4)
rm(list=ls(all=TRUE))
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions
load("data files/pars.arr.ac.Rdata")
load("data files/ds.nm.all.Rdata")
load("data files/epic.Rdata")     # infectious HIV prevalence
to.plot <- T
real.dat <- F ## Set this to be TRUE to run the analyses on the REAL data (just do it once manually isntead of calling from a control file)

## Real data values
if(real.dat) {
  jobnum=1;simj=1;batchdirnm="";nc=12;sim.nm=NULL;simul=F;sb.sim=F;nc=12;seed.bump=0;interv=10;max.vis=5;ltf.prob=0;rr.ltf.ff=1;rr.ltf.mm=1;rr.ltf.hh=1;rr.ltf.d=1;aniter=10000;anburn=2000;niter=50000;nburn=2500;init.jit=0.6;excl.extram=TRUE;decont=FALSE
}else{
  args=(commandArgs(TRUE)) ## load arguments from R CMD BATCH 
  if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    for(i in 1:length(args)) {
      eval(parse(text=args[[i]]))
    }  }
}

if(!simul) {
  ## Parameters as used in SB simulation
  country <- which(ds.nm=='Uganda')
  epic.ind <- which(colnames(epicf)=='Uganda')
  pars.arr <- out.arr[,,which(in.arr[,1,2]==7),] # median and credible intervals for transmission coefficient estimates for this acute relative hazard
  hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") # six gender-route specific transmission coefficients (*b*efore-, *e*xtra-, from-*p*artner-) for *m*ale & *f*emale
  spars <- pars.arr[hazs,2,country]              # get transmission coefficients from base country
  spars[5:6] <- .007
}

if(!is.null(sim.nm)) { ## if running from a qsub job
  load(sim.nm) ## load simulation
  dpars <- output$rakpars ## Hollingsworth Model Parametes
  output$infpm[,,3]/output$infpm[,2,3]
  infpm <- output$infpm
  init <- dpars
  ldpars <- log(dpars)
  het.gen.sd <- output$pars['het.gen.sd']
  cohsim <- rak.coh.fxn(ts.ap = output$ts, dat = output$evout, start.rak = 1994, end.rak = 1999.5,
                        interv = interv, max.vis = max.vis, ltf.prob = ltf.prob,
                        rr.ltf.ff = rr.ltf.ff, rr.ltf.mm = rr.ltf.mm, rr.ltf.hh = rr.ltf.hh, rr.ltf.d = rr.ltf.d, rr.inc.sdc = rr.inc.sdc,
                        verbose = T, browse = F)
  rm(output); gc() ## free up memory
  if(resamp) { ## Resample simulations of 100,000 couples, to yield sample sizes equivalent to the
    ## original Rakai retrospective cohort data reformat into wawer style line list & do
    ## Wawer style Poisson regressions, controlling for various amounts of heterogeneity
  }
  rcohsim <- rak.wawer(rak.coh = cohsim, excl.extram=excl.extram, decont=decont, start.rak = 1994,
                       het.gen.sd = het.gen.sd,
                       verbose = T, browse=F)
  rakll <- rcohsim$rakll
}else{
  ## sim.nm <- file.path('results','RakAcute','Uganda','Uganda-48100-50.Rdata')
  ## simulation paramters
  hpars <- c(acute.sc = 276/10.6, late.sc = 76/10.6, bp = 10.6/1200, dur.ac = 2.9, dur.lt = 9, dur.aids = 10) ## Hollingsworth 2008 ests
  dpars <- hpars
  ##dpars <- c(acute.sc = 276/10.6, late.sc = 3*76/10.6, bp = 10.6/1200, dur.ac = 2.9, dur.lt = 3, dur.aids = 10) ## other values to sim
  if(simul & !sb.sim)  {
    load(file.path(batchdirnm, 'blocksHollFitTest.Rdata'))
    simpars <- unlist(blocks[jobnum, names(dpars)])
    sim <- holl.mod(i.n,p.n,l.n,dpars=simpars, verbose=F)
    dpars <- simpars
  }
  ldpars <- log(dpars)
}
## simulate using Hollingsworth et al. 2008 model (only do once)


stepper <- 1
excl.by.errs <- c(F,T) ##ifelse(sb.sim, c(F,T), F)
for(excl.by.err in excl.by.errs) { ## Fit both with & without excluding -- to +- to LTF couples 
  if(simul) { ## simulated data
    if(!sb.sim) { ## simulate using Hollingsworth et al. 2008 model
      if(excl.by.err) stop(paste0("excl.by.err=T & !sb.sim: Hollingsworth simulator doesn't have any loss to follow-up"))
      wtab.do <- sbmod.to.wdat(sim, excl.by.err = excl.by.err, browse=F)
      print(wtab.do)
      nm <- ifelse(excl.by.err,'sim exclErr', 'sim no exclErr')
      ## ndirnm.bs <- file.path('results','HollingsworthAn','SimNOExclbyErr')
      ## while(file.exists(paste0(ndirnm.bs,'-',ss))) ss <- ss+1
      ## ndirnm <- paste0(ndirnm.bs,'-',ss)
      ndirnm <- file.path(batchdirnm,jobnum)
      ss <- 1
    }else{ ## SB cohort simulated data loaded from file & converted to rakll format abovve
      wtab.do <- sbmod.to.wdat(rakll, excl.by.err = excl.by.err, browse=F)
      nm <- paste0('sim-',jobnum, ' Xerr'[excl.by.err], 'Xext'[excl.extram])
      ndirnm <- file.path(batchdirnm,jobnum)
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
  }
  if(!file.exists(ndirnm)) dir.create(ndirnm)
  print(ndirnm)
  print(wtab.do)
  holl.lik(log(dpars), wtab.do, verbose = T, browse = F)
  ## Initialize a covariance matrix with -.9 correlation between acute/late scalar & duration.
  sd.props <- c(.1, .3, .01, .2, .5, .3)*1.2  ## Guess, works well for most fits & adapts anyways.
  names(sd.props) <- names(ldpars)
  sigma <- matrix(0, nr = length(sd.props), nc = length(sd.props))
  rownames(sigma) <- names(ldpars)
  colnames(sigma) <- names(ldpars)
  diag(sigma) <- sd.props^2
  sigma[1,4] <- -sd.props[1]*sd.props[4]*.9
  sigma[4,1] <- -sd.props[1]*sd.props[4]*.9
  sigma[2,5] <- -sd.props[2]*sd.props[5]*.9
  sigma[5,2] <- -sd.props[2]*sd.props[5]*.9
  gc()
  d.out <- mclapply((seed.bump + 1:nc), holl.wrp, jit = init.jit,
                    sd.props = sd.props, force.inits = ldpars, wtab = wtab.do, 
                    multiv = T, covar = sigma, 
                    verbose = T, verbose2 = F, tell = 100, 
                    niter = aniter, nthin = 1, nburn = anburn, browse=F)
  pout <- procpost(d.out, ldpars = ldpars, dirnm=ndirnm, nm = paste0('adapt1',nm), to.plot = to.plot)
  sigma1 <- pout$sigma
  ## Adaptive phase #2, this time with sigma from previous run
  d.out2 <- mclapply((seed.bump + 1:nc), holl.wrp, force.inits = ldpars, jit = init.jit,
                     sd.props = sd.props,  wtab = wtab.do,
                     multiv = T, covar = sigma1, 
                     verbose = T, verbose2 = F, tell = 100, 
                     niter = aniter, nthin = 1, nburn = anburn, browse=F)
  pout2 <- procpost(d.out2, ldpars = ldpars, dirnm=ndirnm, nm = paste0('adapt2',nm), to.plot = to.plot)
  sigma2 <- pout2$sigma
####################################################################################################
  ## fitting phase with multivariate block sampling
  f.out <- mclapply((seed.bump + 1:nc), holl.wrp, force.inits = ldpars, jit = init.jit,
                    sd.props = sd.props,  wtab = wtab.do,
                    multiv = T, covar = sigma2, 
                    verbose = T, verbose2 = F, tell = 100, 
                    niter = niter, nthin = 1, nburn = nburn, browse=F)
  fout <- procpost(f.out, ldpars = ldpars, dirnm=ndirnm, nm = paste0('final',nm), to.plot = to.plot)
  print(fout$gel)
  pdf(file.path(ndirnm,'mcmc summ.pdf')); plot(fout$mcmc); dev.off()
  if(stepper==1) {
    gel <- fout$gel
    outtab <- fout$outtab
    rownames(outtab)[4] <- ifelse(simul, 'true values', 'Hollingsworth Estimates')
    stepper <- stepper + 1
  }else{
    gel <- c(gel, fout$gel)
    outtab <- abind(outtab, fout$outtab, along = 3)
    dimnames(outtab)[[3]] <- c('base', 'XbErr')
  }
  rm(pout, pout2, fout)
  gc()
}

if(!simul) {
  write.csv(signif(outtab[,,'base'],3), file.path(ndirnm, 'outtab base.csv'))
  write.csv(signif(outtab[,,'XbErr'],3), file.path(ndirnm, 'outtab XbErr.csv'))
  save.image(file.path(ndirnm,'workspace.Rdata'))
}else{ ## Save results as list
  if(sb.sim) fitout <- list(jobnum = jobnum, outtab.h = outtab, outtab.w = rcohsim$armod, erhs.w = rcohsim$erhs,
                            wtab = wtab.do, rakll = rcohsim$rakll, gel = gel, infpm = infpm,
                            dpars = dpars, het.gen.sd = het.gen.sd, 
                            excl.extram = excl.extram, decont = decont,
                            interv = interv, max.vis = max.vis, ltf.prob = ltf.prob,
                            rr.ltf.ff = rr.ltf.ff, rr.ltf.mm = rr.ltf.mm, rr.ltf.hh = rr.ltf.hh, rr.ltf.d = rr.ltf.d)
  fitout <- list(jobnum = jobnum, outtab.h = outtab, ## Hollingsworth-model generated data
                 wtab = wtab.do, gel = gel, dpars = dpars, 
                 excl.extram = excl.extram, decont = decont,
                 interv = interv, max.vis = max.vis, ltf.prob = ltf.prob,
                 rr.ltf.ff = rr.ltf.ff, rr.ltf.mm = rr.ltf.mm, rr.ltf.hh = rr.ltf.hh, rr.ltf.d = rr.ltf.d)
  if(!file.exists(file.path(batchdirnm, 'fitouts'))) dir.create(file.path(batchdirnm, 'fitouts'))
  save(fitout, file = file.path(batchdirnm, 'fitouts', paste0('fitout-',jobnum,'-ltf',signif(ltf.prob,3),'.Rdata')))
}
#load(file.path(ndirnm,'workspace.Rdata'))


  ## excl.by.err <- F
  ## wtab.do <- sbmod.to.wdat(rcohsim$rakll, excl.by.err = excl.by.err, browse=F)
  ## print(wtab.do)
  ## rc <- rcohsim
  ## rc$erhs
  ## with(rc,    {print(armod[,'ehm.ac',,]); print(armod['true value','ehm.ac',1,1])})
  ## with(rc,    {print(armod['50 %','ehm.ac',,]); print(armod['true value','ehm.ac',1,1])})
  ## with(rc,    {print(armod['50 %','late.sc',,]); print(armod['true value','late.sc',1,1])})
  ## with(rc,    {print(armod['50 %','dur.lt',,]); print(armod['true value','dur.lt',1,1])})          
  ## with(rc,    {print(armod['50 %','ehm.ltaids',,]); print(armod['true value','ehm.ltaids',1,1])})
  ## with(rc,    {print(armod['50 %','ehm.lt',,]); print(armod['true value','ehm.lt',1,1])})      
  ## with(rc,    {print(armod['50 %','bp',,]); print(armod['true value','bp',1,1])})

    ## pdf(file.path(batchdirnm, 'test.pdf'))

    
