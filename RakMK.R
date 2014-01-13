####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################
#rm(list=ls())                                  # clear workspace
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
load("data files/ds.nm.all.Rdata") # country names
load('data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
#load('data files/CFJobsToDo.Rdata') ## for finishing up jobs from last run that didn't get finished due to cluster problems.
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
## source('RakMK.R')

####################################################################################################
####################################################################################################
##   COUNTERFACTUAL ANALYSIS
####################################################################################################
## Using an HIV transmission model fit to Demographic and Health Surveys in 18 sub-Saharan African
## countries, we used counter-factual simulations to examine how serodiscordant proportions are
## affected by AIDS mortality rates and pre-couple, extra-couple, and within-couple HIV transmission
## rates.
####################################################################################################
countries <- 1:length(ds.nm)
countries <- which(ds.nm=='Uganda')
each.val <- 100                          #  number of couples per couple formation (marital) cohort
counterf.betas <- F                       # change betas in counterfactuals? if not change beta_within & c's (so beta_within affects all routes)
sub.betas <- F                           # substitute betas? if not beta_within & c's
rtsc <- c(0, 1/10, 1/5, 1/2, 1, 2, 5, 10)  # transmission route scalars
nrtsc <- length(rtsc)                      # how many scalars?
acutes <- as.numeric(in.arr[,1,2])    # acute phase relative hazards we used to fit in fitting phase
nac <- length(acutes)                 # how many are there?
hsds <- c(.5,1,2,3)                     # standard deviation of log(hazards)
hsds.rts <- c(0,1,2)                      # same but use smaller subset when also scaling transmission routes
cors <- c(0,.4,.8)                      # inter-partner correlations
nn <- 400 # number of simulations per country-acute combination (must be bigger than max(sel) later but is fine to leave big
substitute <- F                         # not a substitution analysis
totn <- 0                               # total number of simulations (steps up to final value)
num.doing <- 0
sink("RakAcute.txt")         # create a control file to send to the cluster
outdir <- file.path('results','RakAcute')
if(!file.exists(outdir))      dir.create(outdir) # create directory if necessary
for(aa in acutes)  {                    # loop through acute phase relative hazard
  acdirnm <- file.path(outdir,paste0('Acute', aa)) # set up directory for each acute phase relative hazard
  if(!file.exists(acdirnm))      dir.create(acdirnm) # create directory if necessary
  for(cc in countries) {
    batchdirnm <- file.path(acdirnm, ds.nm[cc]) # setup directory for country
    if(!file.exists(batchdirnm))      dir.create(batchdirnm) # created if necessary
    if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) # setup directory to store Rout files
######################################################################
######################################################################           
    ## Set defaults for all parameters for each simulation, simulatin specific-values set later
######################################################################
    group <- rep(cc,nn)          # country group.
    ## set substitution country (donor country) indices to country:
    ## s.epic, s.demog, s.bmb...
    for(svar in c('epic','demog', hazs)) assign(paste0('s.',svar), rep(cc,nn))
    ## set all phases to have same infectivity (change below)
    for(ph in c('acute','late', 'aids')) assign(paste0(ph,'.sc'), rep(1,nn))
    death <- rep(T,nn)       # include death
    ## all haz scalars (bmb.sc,etc...) set to 1
    for(hh in hazs) assign(paste0(hh,'.sc'), rep(1,nn))
    ## set heterogeneity defaults (i.e. het.b, het.b.sd, het.b.cor,etc...)
    for(ht in c('b','e','p','gen','beh')) {
      assign(paste0('het.',ht),         rep(F,nn))
      assign(paste0('het.',ht,'.sd'),   rep(0,nn))
      assign(paste0('het.',ht,'.cor'),  rep(0,nn))
    }
    scale.by.sd <- rep(T,nn)     # scale by standard deviation?
    scale.adj <- rep(1,nn)       # arbitrary scalar if not doing that.
    infl.fac <- rep(200,nn)  # inflation factor for non-parametric couple pseudo-population builder
    maxN <- rep(10^5,nn)     # max pseudopopulation size
    sample.tmar <- rep(F,nn) # sample marital (couple formation) date from copulas?
    psNonPar <- rep(F,nn) #  use non-parametric couple pseudo-population builder?
    each <- rep(each.val, nn) # how many couples per marital (couple formation) cohort
######################################################################
######################################################################
    ## simulation-specific settings
    ## 
    ## as fitted
    sel <- 1
    blocks <- data.frame(start = min(sel), end = max(sel), lab = 'as fitted')
    for(hsd in hsds.rts)  {            # with varying amounts of within-couple heterogeneity
      ## no AIDS mortality
      sel <- sel[length(sel)] + 1
      death[sel] <- F
      het.gen[sel] <- hsd > 0
      het.gen.sd[sel] <- hsd
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('no AIDS mortality',paste0(' gen het=',hsd))))
      ## scale pre-couple
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + nrtsc)
      bmb.sc[sel] <- rtsc
      bfb.sc[sel] <- rtsc
      het.gen[sel] <- hsd > 0
      het.gen.sd[sel] <- hsd
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('scale pre-couple',paste0(' gen het=',hsd))))
      ## scale extra-couple
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + nrtsc)
      bme.sc[sel] <- rtsc
      bfe.sc[sel] <- rtsc
      het.gen[sel] <- hsd > 0
      het.gen.sd[sel] <- hsd
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('scale extra-couple',paste0(' gen het=',hsd))))  
      ## scale within-couple
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + nrtsc)
      bmp.sc[sel] <- rtsc
      bfp.sc[sel] <- rtsc
      het.gen[sel] <- hsd > 0
      het.gen.sd[sel] <- hsd
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('scale within-couple',paste0(' gen het=',hsd))))
    }
    for(cr in cors) {      # for each inter-partner correlation
      ## pre-couple heterogeneity (no heterogeneity for other routes)
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + length(hsds))
      het.b[sel] <- T
      het.b.sd[sel] <- hsds
      het.b.cor[sel] <- cr        
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('pre- heterogeneity cor=',cr)))
      ## extra-couple heterogeneity (no heterogeneity for other routes)
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + length(hsds))
      het.e[sel] <- T
      het.e.sd[sel] <- hsds
      het.e.cor[sel] <- cr        
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('extra- heterogeneity cor=',cr)))
      ## within-couple heterogeneity (no heterogeneity for other routes)
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + length(hsds))
      het.p[sel] <- T
      het.p.sd[sel] <- hsds
      het.p.cor[sel] <- cr        
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('within- heterogeneity cor=',cr)))
      ## genetic heterogeneity (same individual risk factor for all routes)
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + length(hsds))
      het.gen[sel] <- T
      het.gen.sd[sel] <- hsds
      het.gen.cor[sel] <- cr
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('genetic heterogeneity cor=',cr)))
      ## behavioral heterogeneity (same individual risk factor for pre-/extra- routes, no heterogeneity for within-)
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + length(hsds))
      het.beh[sel] <- T
      het.beh.sd[sel] <- hsds
      het.beh.cor[sel] <- cr        
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('behavioral heterogeneity cor=',cr)))
      ## all route heterogeneity but different risk factors
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + length(hsds))
      het.b[sel] <- T
      het.b.sd[sel] <- hsds
      het.b.cor[sel] <- cr        
      het.e[sel] <- T
      het.e.sd[sel] <- hsds
      het.e.cor[sel] <- cr        
      het.p[sel] <- T
      het.p.sd[sel] <- hsds
      het.p.cor[sel] <- cr        
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('all route heterogeneity cor=',cr)))
      ## pre-/extra- route heterogeneity but different risk factors
      sel <- (sel[length(sel)] + 1):(sel[length(sel)] + length(hsds))
      het.b[sel] <- T
      het.b.sd[sel] <- hsds
      het.b.cor[sel] <- cr        
      het.e[sel] <- T
      het.e.sd[sel] <- hsds
      het.e.cor[sel] <- cr        
      blocks <- rbind(blocks, data.frame(start = min(sel), end = max(sel), lab = paste0('pre extra route heterogeneity cor=',cr)))
    }  
    ## LEFT OFF HERE
    for(ii in 1:max(sel)) {
      jb <- ii                   # job num
      totn <- totn+1             # total jobs
      cmd <- paste("R CMD BATCH '--args jobnum=", totn, " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
                   " group.ind=", group[ii], " substitute=", substitute, " sub.betas=", sub.betas, " counterf.betas=", counterf.betas,
                   " s.epic=", s.epic[ii],  " s.demog=", s.demog[ii],
                   " s.bmb=", s.bmb[ii], " s.bfb=", s.bfb[ii],
                   " s.bme=", s.bme[ii], " s.bfe=", s.bfe[ii],
                   " s.bmp=", s.bmp[ii], " s.bfp=", s.bfp[ii], 
                   " death=", death[ii],
                   " acute.sc=", aa, " late.sc=", late.sc[ii]," aids.sc=", aids.sc[ii], # acute phase varying throughout loop
                   " bmb.sc=", bmb.sc[ii], " bfb.sc=", bfb.sc[ii],
                   " bme.sc=", bme.sc[ii], " bfe.sc=", bfe.sc[ii],
                   " bmp.sc=", bmp.sc[ii], " bfp.sc=", bfp.sc[ii],
                   " het.b=", het.b[ii], " het.b.sd=", het.b.sd[ii], " het.b.cor=", het.b.cor[ii],
                   " het.e=", het.e[ii], " het.e.sd=", het.e.sd[ii], " het.e.cor=", het.e.cor[ii],
                   " het.p=", het.p[ii], " het.p.sd=", het.p.sd[ii], " het.p.cor=", het.p.cor[ii],                     
                   " het.gen=", het.gen[ii], " het.gen.sd=", het.gen.sd[ii], " het.gen.cor=", het.gen.cor[ii],
                   " het.beh=", het.beh[ii], " het.beh.sd=", het.beh.sd[ii], " het.beh.cor=", het.beh.cor[ii],
                   " scale.by.sd=", scale.by.sd[ii], " scale.adj=", scale.adj[ii],
                   " infl.fac=", infl.fac[ii], " maxN=", maxN[ii], " sample.tmar=", sample.tmar[ii],
                   " psNonPar=", psNonPar[ii], " seed=1 tmar=(60*12):(99*12) each=", each[ii],
                   " start.rak=1994 end.rak=1999 return.ts=TRUE",
                   " one.couple=F",
                   " tint=99*12' SimulationStarter.R ", file.path(batchdirnm, "routs", paste0(ds.nm[group[ii]], ii, ".Rout")), sep='')
                                        #      if(ii > 158) {
  #     if(totn %in% jtd & ii %in% 89:92 & aa==7) { ## for finishing up jobs that didn't get properly submitted (cluster issues sometimes)
#      if(ii < 27 | ii <)) {
#          if((ii > 26 & ii < 77) & aa %in% c(1,7,25,50) & totn %in% jtd) { # acute/het sensitivity (called Ac7 still)
#          if((ii %in% c(77:84,89:92)) & aa %in% 7 & totn %in% jtd) { # het sensitivity at acute of 7
#      if(ii < 77 & aa%in% c(1,7,25,50) & totn %in% jtd) {
      if(ii==1) {
            num.doing <- num.doing+1
            cat(cmd)               # add command
            cat('\n')              # add new line
          }
        }
    } 
}
sink()
blocks
totn
save(blocks, file = file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
####################################################################################################
print(totn)
print(num.doing)
