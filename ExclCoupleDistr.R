library(mnormt);library(coda); library(stats4);library(plyr)#library(lme4)
rm(list=ls(all=TRUE))
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions
load("data files/pars.arr.ac.Rdata")
load("data files/ds.nm.all.Rdata")
load("data files/epic.Rdata")     # infectious HIV prevalence
to.plot <- T
real.dat <- F ## Set this to be TRUE to run the analyses on the REAL data (just do it once manually isntead of calling from a control file)

sim.nm="results/RakAcute/Uganda/Uganda-96200-911.Rdata"

args=(commandArgs(TRUE)) ## load arguments from R CMD BATCH 
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }  }

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
cov.mods <- T
rcohsim <- rak.wawer(rak.coh = cohsim, excl.extram=excl.extram, decont=decont, start.rak = 1994,
                     het.gen.sd = het.gen.sd,
                     verbose = T, browse=F)
rakll <- rcohsim$rakll

sim <- holl.mod(i.n,p.n,l.n,dpars=simpars, verbose=F)
if(excl.by.err) stop(paste0("excl.by.err=T & !sb.sim: Hollingsworth simulator doesn't have any loss to follow-up"))
wtab.base <- sbmod.to.wdat(sim, excl.by.err = F, browse=F)
wtab.XbErr <- sbmod.to.wdat(sim, excl.by.err = T, browse=F)

print(wtab.do)

####################################################################################################
## How many couples were excluded?
wtabs.ls <- lapply(alldat, function(x) { list(wtab = x$wtab,x$wtab$prevt[1,'n'],x$wtab$latet[nrow(x$wtab$latet),'n'])})
wtabs <- as.data.frame(wtabs)
colnames(wtabs) <- c('incn','prevn','laten')
wtabs$job <- jobnums
apply(wtabs[,-4],2, function(x) quantile(x, c(.025,.5,.975))) ## summary of sample sizes by couple classing
