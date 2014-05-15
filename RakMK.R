####################################################################################################
## Makes control files for each Ugandan couples cohort simulation within which each line giving one
## R CMD BATCH command line to run on a cluster.
####################################################################################################
## rm(list=ls())                                  # clear workspace
load("data files/ds.nm.all.Rdata")                      ## DHS country names
load('data files/pars.arr.ac.Rdata')   ## load hazards for each transmission route as fit to DHS data (in.arr[,,2])
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') ##  transmission coefficient names, for convenience
nc <- 12                                       ## cores per simulation

####################################################################################################
####################################################################################################
##   Simulate couples cohorts with different acute & late phase characteristics
####################################################################################################
cc <- which(ds.nm=='Uganda') ## to get Ugandan prameters from DHS
each.val <- 200 ##  equates to ~100,000 couples
blocks <- expand.grid(acute.sc = c(1,2,5,7, seq(10,50,by=5)),
                      dur.ac = seq(.5,5, by = .5),
                      het.gen.sd = seq(0,3, by = .5),
                      dur.lt = c(5,10), dur.aids = 10, late.sc = c(2,5,10), aids.sc=0)
## added a second batch with more late.sc values, later, kept in this form to keep track in analysis
blocks.add <- expand.grid(acute.sc = c(1,2,5,7, seq(10,50,by=5)), 
                      dur.ac = seq(.5,5, by = .5),
                      het.gen.sd = seq(0,3, by = .5),
                      dur.lt = c(5,10), dur.aids = 10, late.sc = 1, aids.sc=0)
blocks <- rbind(blocks,blocks.add)
blocks$het.gen <- blocks$het.gen.sd > 0
blocks$jobnum <- 1:nrow(blocks)
nn <- nrow(blocks)

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
## set all phases to have same infectivity (changed below)
for(ph in c('acute','late', 'aids')) assign(paste0(ph,'.sc'), rep(1,nn))
## all haz scalars (bmb.sc,etc...) set to 1  (scalars were used for another project & not needed here)
for(hh in hazs) assign(paste0(hh,'.sc'), rep(1,nn))
## set heterogeneity defaults (i.e. het.b, het.b.sd, het.b.cor,etc...; only het.gen.sd is used in this project)
for(ht in c('b','e','p','gen','beh')) {
  assign(paste0('het.',ht),         rep(F,nn))
  assign(paste0('het.',ht,'.sd'),   rep(0,nn))
  assign(paste0('het.',ht,'.cor'),  rep(0,nn))
}
scale.by.sd <- rep(T,nn)     # scale by standard deviation? (account for mean of lognormal distribution not being mu)
scale.adj <- rep(1,nn)       # arbitrary scalar if not doing that. (not doing)
infl.fac <- rep(200,nn)  # inflation factor for non-parametric couple pseudo-population builder (for building simulated couples)
maxN <- rep(10^5,nn)     # max pseudopopulation size (for building simulated couples)
sample.tmar <- rep(F,nn) # sample marital (couple formation) date from copulas? No, use marital cohorts
psNonPar <- rep(F,nn) #  use non-parametric couple pseudo-population builder? (no use copula models)
each <- rep(each.val, nn) # how many couples per marital (couple formation) cohort
to.do <- 1:nrow(blocks) # which rows in blocks to do (default is all)
 
num.doing <- 0 ## initiate coutners
totn <- 0
sink("RakAcute.txt") ## create a control file to send to the cluster
## ####################################################################
for(ii in to.do) {
  jb <- ii                   # job num
  totn <- totn+1             # total jobs
  cmd <- paste("R CMD BATCH '--args jobnum=", blocks$jobnum[ii], " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
               " group.ind=", group[ii], " substitute=FALSE sub.betas=FALSE counterf.betas=FALSE",
               " s.epic=", cc,  " s.demog=", cc, ## get parameters for appropriate country (Uganda in this case) from UNAIDS data or previous DHS fits
               " s.bmb=", cc, " s.bfb=", cc, ## get parameters for appropriate country (Uganda in this case) from UNAIDS data or previous DHS fits
               " s.bme=", cc, " s.bfe=", cc, ## get parameters for appropriate country (Uganda in this case) from UNAIDS data or previous DHS fits
               " s.bmp=", cc, " s.bfp=", cc, ## get parameters for appropriate country (Uganda in this case) from UNAIDS data or previous DHS fits
               " death=TRUE", ## model includes Weibull HIV mortality
               " acute.sc=", blocks$acute.sc[ii], " late.sc=", blocks$late.sc[ii]," aids.sc=", blocks$aids.sc[ii], # acute phase varying throughout loop
               " dur.ac=", blocks$dur.ac[ii], " dur.lt=", blocks$dur.lt[ii], " dur.aids=", blocks$dur.aids[ii],
               " bmb.sc=", bmb.sc[ii], " bfb.sc=", bfb.sc[ii], ## Not scaling hazards (this is from another project, all these = 1)
               " bme.sc=", bme.sc[ii], " bfe.sc=", bfe.sc[ii], ## Not scaling hazards (this is from another project, all these = 1)
               " bmp.sc=", bmp.sc[ii], " bfp.sc=", bfp.sc[ii], ## Not scaling hazards (this is from another project, all these = 1)
               " het.b=", het.b[ii], " het.b.sd=", het.b.sd[ii], " het.b.cor=", het.b.cor[ii], ## these heterogeneities not included in this project
               " het.e=", het.e[ii], " het.e.sd=", het.e.sd[ii], " het.e.cor=", het.e.cor[ii], ## these heterogeneities not included in this project
               " het.p=", het.p[ii], " het.p.sd=", het.p.sd[ii], " het.p.cor=", het.p.cor[ii], ## these heterogeneities not included in this project
               " het.gen=", blocks$het.gen[ii], " het.gen.sd=", blocks$het.gen.sd[ii], " het.gen.cor=", het.gen.cor[ii], ## INCLUDED
               " het.beh=", het.beh[ii], " het.beh.sd=", het.beh.sd[ii], " het.beh.cor=", het.beh.cor[ii], ## these heterogeneities not included in this project
               " scale.by.sd=", scale.by.sd[ii], " scale.adj=", scale.adj[ii],
               " infl.fac=", infl.fac[ii], " maxN=", maxN[ii], " sample.tmar=", sample.tmar[ii],
               " psNonPar=", psNonPar[ii], " seed=1 tmar=(60*12):(100*12) each=", each[ii],
               " start.rak=1994 end.rak=2000 return.ts=TRUE", ## cohort dates, return a time series as output
               " one.couple=F", ## debugging with one couple replicated a bunch of times
               ## interview date, call simulation script and specify output
               " tint=100*12' SimulationStarter.R ", file.path(batchdirnm, "routs", paste0(ds.nm[group[ii]], ii, ".Rout")), sep='') 
  num.doing <- num.doing+1
  cat(cmd)               # add command
  cat('\n')              # add new line
}
sink()
head(blocks)
totn
save(blocks, file = file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
####################################################################################################
print(totn) ## total jobs
print(num.doing) ## total doing in this control file (allows us to just do a few if we changed to.do above)
head(blocks[to.do,])
