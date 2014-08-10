####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster. For fitting SB cohort simulations across sensitivity analysis.
####################################################################################################
## rm(list=ls())                                  # clear workspace
load("data files/ds.nm.all.Rdata") # country names
load('data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
outdir <- file.path('results','RakAcute')
load(file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks describing simulations
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
## source('FitRakMK.R')


sb.sim <- T
## sim.nm <- 'Uganda-48100-'
sim.nm <- 'Uganda-96200-'
## Find those that have already been simulated
fls <- list.files(file.path('results','RakAcute','Uganda'), pattern = '.Rdata')
fls <- fls[grepl(sim.nm,fls)]
fls <- sub(sim.nm,'', fls)
fls <- sub('.Rdata','', fls)
fls <- as.numeric(fls)
fls <- fls[order(fls)]
sim.nm <- file.path('results','RakAcute','Uganda',sim.nm)
length(fls)

## Find those that have been already fit (i.e. have a fit folder)
fit.fls <- list.files(file.path('results','RakAcute','UgandaFits','fitouts'))
fit.fls <- sub('fitout-','', fit.fls)
fit.fls <- sub('-ltf0.0288.Rdata','', fit.fls)
fit.fls <- as.numeric(fit.fls)
fit.fls <- fit.fls[order(fit.fls)]

##to.do <- 1:nrow(blocks)
to.do <- fls
to.do <- with(blocks, which(dur.lt==10 & late.sc==5))
##to.do <- to.do[!to.do %in% fit.fls] ## which haven't been fit yet
##to.do <- jobnums.to.do

blocks[head(fls,10),]
blocks[with(blocks, which(acute.sc==7 & het.gen.sd==1.5 & late.sc==5 & dur.ac == 3 & dur.lt == 10)),]
blocks[with(blocks, which(acute.sc==7 & het.gen.sd==0 & late.sc==5 & dur.ac == 3 & dur.lt == 10)),]

####################################################################################################
####################################################################################################
##   MCMC Fitting of SB Cohort Simulations
####################################################################################################
## Fitting simulated couples data with Wawer & Hollingsworth models to explore potential biases.
####################################################################################################
batchdirnm <- file.path('results','RakAcute','UgandaFitsExtram')
if(!file.exists(batchdirnm))      dir.create(batchdirnm) # create directory if necessary
if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) # setup directory to store Rout files
max.vis <- 5
interv <- 10
ltf.prob <- -log(1-.25)/10 ## loss to follow-up probability of each couple, this yields 35% LOSS over 10 months as observed
rr.ltf.ff <- 1.5     ## yields loss to follow-up of 50% for SDCs
rr.ltf.mm <- 1.5     ## yields loss to follow-up of 50% for SDCs
rr.ltf.hh <- 1 ## if ++?
rr.ltf.d <- 1 ## if dead?
rr.inc.sdc <- 1.5 ## how much faster is the rate at which incident SDC (ie who were previously -- in a survey visit) are LTF than regular SDC?

excl.extram <- F ## exclude couples with 2nd partner infected extra-couply?
decont <- F      ## decontaminate? (remove prevalent couples with person-time exposed to acute/late partners etc)
seed.bump <- 0
aniter <- 5000
anburn <- 1000
niter <- 10000
nburn <- 1500
init.jit <- .45
nc <- 12

num.doing <- 0
totn <- 0
sink("FitRak.txt")         # create a control file to send to the cluster
## ####################################################################
for(ii in to.do) {
  jb <- ii                   # job num
  totn <- totn+1             # total jobs
  cmd <- paste("R CMD BATCH '--args jobnum=", blocks$jobnum[ii], " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
               " sim.nm=\"", sim.nm, blocks$jobnum[ii], ".Rdata\" simul=T sb.sim=", sb.sim, " nc=12 seed.bump=", seed.bump,
               " interv=", interv, " max.vis=", max.vis, " ltf.prob=", ltf.prob,
               " rr.ltf.ff=", rr.ltf.ff, " rr.ltf.mm=", rr.ltf.mm, " rr.ltf.hh=", rr.ltf.hh, " rr.ltf.d=", rr.ltf.d, " rr.inc.sdc=", rr.inc.sdc,
               " aniter=", aniter, " anburn=", anburn, " niter=", niter, " nburn=", nburn,
               " init.jit=", init.jit, " excl.extram=", excl.extram, " decont=", decont,
               "' RakaiCohortSimulator.R ", file.path(batchdirnm, "routs", paste0('Fit', ii, ".Rout")), sep='')
  num.doing <- num.doing+1
  cat(cmd)               ## add command
  cat('\n')              ## add new line
}
sink()
totn

####################################################################################################
print(totn)
print(num.doing)
