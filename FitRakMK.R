####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster. For fitting SB cohort simulations across sensitivity analysis.
####################################################################################################
## rm(list=ls())                                  # clear workspace
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
load("data files/ds.nm.all.Rdata") # country names
load('data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
outdir <- file.path('results','RakAcute')
load(file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
## source('FitRakMK.R')

sim.nm <- 'Uganda-96200-'
fls <- list.files(file.path('results','RakAcute','Uganda'), pattern = '.Rdata')
fls <- fls[grepl(sim.nm,fls)]
fls <- sub(sim.nm,'', fls)
fls <- sub('.Rdata','', fls)
fls <- as.numeric(fls)
fls <- fls[order(fls)]
sim.nm <- file.path('results','RakAcute','Uganda',sim.nm)

blocks[head(fls,10),]
blocks[with(blocks, which(acute.sc==7 & dur.ac==2 & het.gen.sd==0)),]


####################################################################################################
####################################################################################################
##   MCMC Fitting of SB Cohort Simulations
####################################################################################################
## Fitting simulated couples data with Wawer & Hollingsworth models to explore potential biases.
####################################################################################################
batchdirnm <- file.path('results','RakAcute','UgandaFits')
if(!file.exists(batchdirnm))      dir.create(batchdirnm) # create directory if necessary
if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) # setup directory to store Rout files
max.vis <- 5
interv <- 10
ltf.prob <- -log(1-.5)/10 ## loss to follow-up probability of each couple, this yields 50% LOSS over 10 months as observed
rr.ltf.ff <- 1  ## assume all these are equal for now
rr.ltf.mm <- 1
rr.ltf.hh <- 1
rr.ltf.d <- 1
excl.extram <- T ## exclude couples with 2nd partner infected extra-couply?
decont <- F      ## decontaminate? (remove prevalent couples with person-time exposed to acute/late partners etc)
seed.bump <- 0
aniter <- 5000
anburn <- 1000
niter <- 10000
nburn <- 1500
init.jit <- .6
nc <- 12


to.do <- 1:nrow(blocks)
to.do <- c(55:72,1783:1800)

num.doing <- 0
totn <- 0
sink("FitRak.txt")         # create a control file to send to the cluster
## ####################################################################
for(ii in to.do) {
  jb <- ii                   # job num
  totn <- totn+1             # total jobs
  cmd <- paste("R CMD BATCH '--args jobnum=", blocks$jobnum[ii], " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
               " sim.nm=\"", sim.nm, blocks$jobnum[ii], ".Rdata\" simul=T sb.sim=T nc=12 seed.bump=", seed.bump,
               " interv=", interv, " max.vis=", max.vis, " ltf.prob=", ltf.prob,
               " rr.ltf.ff=", rr.ltf.ff, " rr.ltf.mm=", rr.ltf.mm, " rr.ltf.hh=", rr.ltf.hh, " rr.ltf.d=", rr.ltf.d, 
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
