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
cc <- which(ds.nm=='Uganda')
each.val <- 100                          #  number of couples per couple formation (marital) cohort
counterf.betas <- F                       # change betas in counterfactuals? if not change beta_within & c's (so beta_within affects all routes)
sub.betas <- F                           # substitute betas? if not beta_within & c's
blocks <- expand.grid(acute.sc = c(1:9,seq(10,50,by=5)),
                      dur.ac = seq(.5,8, by = .5),
                      het.gen.sd = 0:3)
blocks$het.gen <- blocks$het.gen.sd > 0
#acutes <- as.numeric(in.arr[,1,2])    # acute phase relative hazards we used to fit in fitting phase
substitute <- F                         # not a substitution analysis
num.doing <- 0
sink("RakAcute.txt")         # create a control file to send to the cluster
outdir <- file.path('results','RakAcute')
if(!file.exists(outdir))      dir.create(outdir) # create directory if necessary
batchdirnm <- file.path(outdir, ds.nm[cc]) # setup directory for country
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
## LEFT OFF HERE
for(ii in 1:max(sel)) {
  jb <- ii                   # job num
  totn <- totn+1             # total jobs
  cmd <- paste("R CMD BATCH '--args jobnum=", totn, " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
               " group.ind=", group[ii], " substitute=", substitute, " sub.betas=", sub.betas, " counterf.betas=", counterf.betas,
               " s.epic=", s.epic[ii],  " s.demog=", s.demog[ii],
               " s.bmb=", s.bmb[ii], " s.bfb=", s.bfb[ii], # country to substitute in for each input (if doing that)
               " s.bme=", s.bme[ii], " s.bfe=", s.bfe[ii],
               " s.bmp=", s.bmp[ii], " s.bfp=", s.bfp[ii], 
               " death=", death[ii],
               " acute.sc=", blocks$acute.sc[ii], " late.sc=", late.sc[ii]," aids.sc=", aids.sc[ii], # acute phase varying throughout loop
               " dur.ac=", blocks$dur.acs[ii], " dur.lt=", dur.lt[ii], "dur.aids=", dur.aids[ii],
               " bmb.sc=", bmb.sc[ii], " bfb.sc=", bfb.sc[ii],
               " bme.sc=", bme.sc[ii], " bfe.sc=", bfe.sc[ii],
               " bmp.sc=", bmp.sc[ii], " bfp.sc=", bfp.sc[ii],
               " het.b=", het.b[ii], " het.b.sd=", het.b.sd[ii], " het.b.cor=", het.b.cor[ii],
               " het.e=", het.e[ii], " het.e.sd=", het.e.sd[ii], " het.e.cor=", het.e.cor[ii],
               " het.p=", het.p[ii], " het.p.sd=", het.p.sd[ii], " het.p.cor=", het.p.cor[ii],                     
               " het.gen=", blocks$het.gen[ii], " het.gen.sd=", blocks$het.gen.sd[ii], " het.gen.cor=", het.gen.cor[ii],
               " het.beh=", het.beh[ii], " het.beh.sd=", het.beh.sd[ii], " het.beh.cor=", het.beh.cor[ii],
               " scale.by.sd=", scale.by.sd[ii], " scale.adj=", scale.adj[ii],
               " infl.fac=", infl.fac[ii], " maxN=", maxN[ii], " sample.tmar=", sample.tmar[ii],
               " psNonPar=", psNonPar[ii], " seed=1 tmar=(60*12):(100*12) each=", each[ii],
               " start.rak=1994 end.rak=2000 return.ts=TRUE",
               " one.couple=F",
               " tint=100*12' SimulationStarter.R ", file.path(batchdirnm, "routs", paste0(ds.nm[group[ii]], ii, ".Rout")), sep='')
                                        #      if(ii > 158) {
                                        #     if(totn %in% jtd & ii %in% 89:92 & aa==7) { ## for finishing up jobs that didn't get properly submitted (cluster issues sometimes)
                                        #      if(ii < 27 | ii <)) {
                                        #          if((ii > 26 & ii < 77) & aa %in% c(1,7,25,50) & totn %in% jtd) { # acute/het sensitivity (called Ac7 still)
                                        #          if((ii %in% c(77:84,89:92)) & aa %in% 7 & totn %in% jtd) { # het sensitivity at acute of 7
                                        #      if(ii < 77 & aa%in% c(1,7,25,50) & totn %in% jtd) {
                                        #      if(ii==1) {
  num.doing <- num.doing+1
  cat(cmd)               # add command
  cat('\n')              # add new line
                                        #         }
}
sink()
blocks
totn
save(blocks, file = file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
####################################################################################################
print(totn)
print(num.doing)
