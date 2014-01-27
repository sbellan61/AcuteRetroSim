######################################################################
## Drivers of serodiscordance patterns and sub-Saharan Africa
###################################################################### 
## Steve Bellan, 2013
## steve.bellan@gmail.com
#################################################################################################### 
##  The goal of this analysis is to assess the extent to which pre-, extra-, & within-couple
##  transmission, AIDS mortality, and acute phase infectivity drive observed levels of
##  serodiscordance in sub-Saharan Africa as well as to determine which country-specific
##  characteristics can explain their variable serodiscordance proportions and epidemic intensities.
#################################################################################################### 
## This code is called by R CMD BATCH with input arguments (each line of a control file sent to a
## cluster), runs a simulation using functions from sim fxns3.R, and then saves the output in a
## specified directory structure
#################################################################################################### 
#rm(list=ls())                           # clear workspace

setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
args=(commandArgs(TRUE))                # load arguments from R CMD BATCH 
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
      }  }
set.seed(seed)
vfreq <- 200 ## how often to report on results
source("SimulationFunctions.R")                   # load simulation functions from script
source("RakFunctions.R") # load Rakai analysis simulation functions from script
load("data files/copula sigmas.Rdata")  # multivariate copula covariance matrix for simulating couple pseudopopulations
load("data files/epic.Rdata")     # infectious HIV prevalence
##  transmission parameters fit to DHS across the range of a acute phase relative hazards
load("data files/pars.arr.ac.Rdata")
## loads out.arr [parname, ci.l med ci.u, acute, country] and in.arr (describes inputs & gelman diagnostics)
load("data files/csurv.Rdata")    #  probability of survival (row) months by age (in months) at seroconversion (column)
load('data files/ds.nm.all.Rdata')        # load country names
load("data files/allDHSAIS.Rdata")         # DHS data
load("data files/dframe.Rdata")   # summary characteristics of DHS data
country <- group.ind              # set country to country-group index
out.dir <- batchdirnm #  output directory
if(!file.exists(out.dir)) dir.create(out.dir) # create this directory if it doesn't exist
odat <- dat          # backup original data in this workspace before it gets manipulated
print(ds.nm[country])                   # Print country name.

######################################################################
## Transmission coefficients
## Load parameters for acute phase assuming (acute.sc)
pars.arr <- out.arr[,,which(in.arr[,1,2]==acute.sc),] # median and credible intervals for transmission coefficient estimates for this acute relative hazard
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") # six gender-route specific transmission coefficients (*b*efore-, *e*xtra-, from-*p*artner-) for *m*ale & *f*emale
spars <- pars.arr[hazs,2,country]              # get transmission coefficients from base country
## If substituting, substitute parameters/epidemic curves out for those from donor country
if(substitute) {
  ##  For each country substitute things from others and see how close
  ##  it gets to the serodiscordance levels of the other country.
######################################################################
  ##  Epidemic curve
  if(s.epic!=country) {     # if substituting, replce it in data (otherwise it defaults to original)
    s.epic.nm <- ds.nm[s.epic]              # epidemic curve to use (country name)
    s.epic.ind <- which(colnames(epicf)==s.epic.nm) # epidemic curve to use (column index of epicf/m matrices)
    ## substitute epidemic curves into data
    dat$epic.ind <- s.epic.ind
    dat$epic.nm <- s.epic.nm
  }
  if(sub.betas) { ## substituting betas between countries
      ##  substitute parameters from substitution (donor) country. At most two of these substitutions will change spars.
      spars['bmb'] <- pars.arr['bmb',2,s.bmb] # pre-couple to male
      spars['bfb'] <- pars.arr['bfb',2,s.bfb] # pre-couple to female
      spars['bme'] <- pars.arr['bme',2,s.bme] # extra-couple to male  
      spars['bfe'] <- pars.arr['bfe',2,s.bfe] # extra-couple to female  
      spars['bmp'] <- pars.arr['bmp',2,s.bmp] # within-couple to male  
      spars['bfp'] <- pars.arr['bfp',2,s.bfp] # within-couple to female  
      ## substitute demography is just done by using parameters & epi curves from the opposite country in psrun() below.
    }else{ ## substituting HIV transmission rate & contact coefficients
      spars['bmb'] <- pars.arr['bmp',2,s.bmp] * pars.arr['cmb',2,s.bmb] # pre-couple to male
      spars['bfb'] <- pars.arr['bfp',2,s.bfp] * pars.arr['cfb',2,s.bfb] # pre-couple to female
      spars['bme'] <- pars.arr['bmp',2,s.bmp] * pars.arr['cme',2,s.bme] # extra-couple to male  
      spars['bfe'] <- pars.arr['bfp',2,s.bfp] * pars.arr['cfe',2,s.bfe] # extra-couple to female  
      spars['bmp'] <- pars.arr['bmp',2,s.bmp] # within-couple to male  
      spars['bfp'] <- pars.arr['bfp',2,s.bfp] # within-couple to female  
    }
}else{
  s.epic.nm <- NA # not substituting epidemic curves, this will cause rcop() to use default country epidemic curves
  s.epic.ind <- NA  
}

## Simulate couples transmission model (calling psrun() from sim fxns3.R). Output (temp) is the name of the file that is produced.
temp <- psrun(country = country,   # what country are we 'simulating'
              s.demog = s.demog,        # what country to use for demograhy
              maxN = maxN, vfreq = vfreq, # maximum # of couples in pseudo-population; Frequency with which to report on simulation progress.
              pars = spars,             # transmission coefficients
              psNonPar = psNonPar,      # nonparametric approach?
              infl.fac = infl.fac, # if non-parametric approach, how many pseudocouples to inflate each couple to
              last.int = F, sample.tmar = sample.tmar, # if non-parametric approach decide how to distribute marriage dates & interview times.
              tmar = tmar, each = each, tint = tint, # parametric approach, couple pseudopop built parametrically from multivariate copulas
              start.rak = start.rak, end.rak = end.rak, # interval for Rakai style retrospective cohort analyses
              return.ts = return.ts,                    # return ts?
              death = death,          # include HIV mortality in model?
              acute.sc = acute.sc, late.sc = late.sc, aids.sc = aids.sc, # relative of infectivity of different phases vs chronic phase
              ##  how to scale transmission coefficients (for counterfactual analysis of how
              ##  changing these routes affects SDP, defaults to 1 for substitution analyses)
              bmb.sc = bmb.sc, bfb.sc = bfb.sc, # pre-couple 
              bme.sc = bme.sc, bfe.sc = bfe.sc, # extra-couple
              bmp.sc = bmp.sc, bfp.sc = bfp.sc, # within-couple
              counterf.betas = counterf.betas, # change betas in counterfactuals? if not change beta_within & c's (so beta_within affects all routes)
              ## Heterogeneity
              het.gen = het.gen,        # do it?
              het.gen.sd = het.gen.sd, # standard deviation of genetic heterogeneity
              het.gen.cor = het.gen.cor, # inter-partner correlation of individual risk deviate for genetic heterogeneity
              ## next 3 are same as previous 3, but for behavioral heterogeneity (risk deviate
              ## amplifies only pre- and extra-couple transmission susceptibility)
              het.beh = het.beh, het.beh.sd = het.beh.sd, het.beh.cor = het.beh.cor,
              ## route-specific heterogeneity (i.e. individuals' hazards are multiplied by a different lognormal risk deviate for each route)
              het.b = het.b, het.b.sd = het.b.sd, het.b.cor = het.b.cor, # pre-couple   
              het.e = het.e, het.e.sd = het.e.sd, het.e.cor = het.e.cor, # extra-couple 
              het.p = het.p, het.p.sd = het.p.sd, het.p.cor = het.p.cor, # within-couple
              ##########
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
