######################################################################
## Couples cohort transmission simulator
###################################################################### 
## Steve Bellan, 2014, steve.bellan@gmail.com
setwd('Rakai/AcuteRetroSim/')
rm(list=ls(all=T))                           # clear workspace
args=(commandArgs(TRUE))                # load arguments from R CMD BATCH 
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
      }  }
vfreq <- 200 ## how often to report on results
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
#set.seed(seed)

parms <- c(bmb=0.009259789, bfb=0.021657899, bme=0.003787767, bfe=0.004677159, bmp=36/48525*10, bfp=36/48525*10)
#parms[1:4] <- parms[1:4]*100
system.time(blah <- psrun(maxN = 4785, browse = F))

out.dir <- batchdirnm #  output directory
if(!file.exists(out.dir)) dir.create(out.dir) # create this directory if it doesn't exist
odat <- dat          # backup original data in this workspace before it gets manipulated
print(ds.nm[group.ind])                   # Print country name.

######################################################################
## Transmission coefficients
## Load pre-couple & extra-couple transmission coefficients estimated from DHS data. 
load("data files/pars.arr.ac.Rdata")
pars.arr <- out.arr[,,which(in.arr[,1,2]==7),] 
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") 
spars <- pars.arr[hazs,2,13]              # get transmission coefficients from base country
spars[5:6] <- lambda
s.epic.ind <- s.epic.nm <- NA # not substituting epidemic curves, this will cause rcop() to use default country epidemic curves

## nc <- 12
## each <- 100
## source("SimulationFunctions.R")                   # load simulation functions from script
## Simulate couples transmission model (calling psrun() from sim fxns3.R). Output (temp) is the name of the file that is produced.
temp <- psrun(country = group.ind,   # what country are we 'simulating'
              s.demog = group.ind,        # what country to use for demograhy
              maxN = maxN, vfreq = vfreq, # maximum # of couples in pseudo-population; Frequency with which to report on simulation progress.
              pars = spars,             # transmission coefficients
              last.int = F, sample.tmar = sample.tmar, # if non-parametric approach decide how to distribute marriage dates & interview times
              tmar = tmar, each = each, tint = tint, # parametric approach, couple pseudopop built parametrically from multivariate copulas
              start.rak = start.rak, end.rak = end.rak, # interval for Rakai style retrospective cohort analyses
              return.ts = return.ts,                    # return ts?
              acute.sc = acute.sc, late.sc = late.sc, aids.sc = aids.sc, # relative of infectivity of different phases vs chronic phase
              dur.ac=dur.ac, dur.lt=dur.lt, dur.aids=dur.aids,
              het.gen = het.gen, het.gen.sd = het.gen.sd, het.gen.cor = het.gen.cor,               ## Heterogeneity
              one.couple = one.couple, # for debugging, repeat one couple married in 1990, 100,000 times
              scale.by.sd = scale.by.sd, # adjust beta means to keep geometric mean constant with increasing heterogeneity (defaults to TRUE)
              scale.adj = scale.adj,   # adjust them arbitrarily
              out.dir = out.dir,       # output directory
              nc = nc,                 # number of cores
              make.jpgs = F,           # create jpgs of results
              browse = F)              # debug


load(temp)                              # load output of simulations
names(output)                           # list objects summarized below
## jobnum=which job; evout=line list of each couple; tss=time series of pseudopopulation, rakacRR=
## retrospective cohort estimate of acute infectivity; pars=input parameters; tmar=marriage cohort
## dates; each=# couples per marriage cohort
## 
head(output$evout)                      # columns explained below
## uid=unique couple identifier; ds=dat set; ser=couple serostatus (1:4::{++,+-,-+,--});
## tms/tfs=male/female sexual debuts (months since 1900); tmar=couple formation (i.e. marriage)
## date; tint=DHS interview/testing date; mardur.mon=couple duration in months; circ=circumciscion
## status of male; mage/fage=male/female age at interview; epic.ind=epidemic curve index in
## allepicm/f; epic.nm=epidemic curve name; group=analysis group; m/fser=m/f serostatus at
## interview; f/mdoi=m/f date of infection; m/fdod=m/f date of death; m/fcoi=m/f cause (i.e. route)
## of infection; m/f.het.gen/beh/b/e/p=lognormal risk deviate for m/f for genetic, behavioral or
## route-specific heterogeneity; m/fcoi.phase=partner's phase (acute/chronic/late) for infection;
## tend=end of couple-time (interview or death); dm/fage=date at which m/f age out of DHS cohort
## (50yr for f, 60yr for m); taend=end of couple-time (interview, aging out of DHS cohort, or
## death); m/falive=m/f alive at interview; tmsdc=time spent as M+F- couple; tfsdc=time spent as
## M-F+ couple; tccc=time spent ++; m/fsdc=ever M/F+ serodicsordant?; ccc=ever ++?
tail(output$tss)                        # columns explained below
## tt=months since 1900; yr=year; b.ss=#couples -- *b*efore couple formation; b.mm/b.ff=# couples
## +-/-+ before cf; b.hh=#couples ++ before cf; ss/ff/mm/hh=# couples in each serostatus formed &
## still alive; d.xx=# couples in each serostatus but >=1 partner dead; a.xx=# couples in each
## serostatus but >=1 partner aged out of DHS cohort (gender that's aged is appended); alive=#
## couples alive; inf.alive=# couples infected & alive; inf=# couples infected; cu.xx=#cumulative
## infections by gender-route combination
print(output$pars)                      # input parameters
print(spars)                            # inputted hazards

rm(list=ls(all=T))                           # clear workspace
