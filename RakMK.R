####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################
## rm(list=ls())                                  # clear workspace
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
load("data files/ds.nm.all.Rdata") # country names
load('data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
## load('data files/CFJobsToDo.Rdata') ## for finishing up jobs from last run that didn't get finished due to cluster problems.
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
## source('RakMK.R')

fls <- list.files(file.path('results','RakAcute','Uganda'), pattern = '.Rdata')
fls <- fls[grepl('96200',fls)]
fls <- sub('Uganda-96200-','', fls)
fls <- sub('.Rdata','', fls)
fls <- as.numeric(fls)
fls <- fls[order(fls)]
length(fls)

####################################################################################################
####################################################################################################
##   Simulate couples cohorts with different acute & late phase characteristics
####################################################################################################
## 
####################################################################################################
cc <- which(ds.nm=='Uganda')
each.val <- 200 ##  number of couples per couple formation (marital) cohort
blocks <- expand.grid(acute.sc = c(1,2,5,7, seq(10,50,by=5)),
                      dur.ac = seq(.5,5, by = .5),
                      het.gen.sd = seq(0,3, by = .5),
                      dur.lt = c(5,10), dur.aids = 10, late.sc = c(2,5,10), aids.sc=0)
blocks.add <- expand.grid(acute.sc = c(1,2,5,7, seq(10,50,by=5)), ## added a second batch with more late.sc later
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
## set all phases to have same infectivity (change below)
for(ph in c('acute','late', 'aids')) assign(paste0(ph,'.sc'), rep(1,nn))
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

to.do <- 1:nrow(blocks)
to.do <- to.do[blocks$dur.lt==10 & blocks$late.sc==1]
to.do <- to.do[!to.do %in% fls]
## to.do <- to.do[!1:nrow(blocks) %in% fls]
 
num.doing <- 0
totn <- 0
sink("RakAcute.txt")         # create a control file to send to the cluster
## ####################################################################
for(ii in to.do) {
  jb <- ii                   # job num
  totn <- totn+1             # total jobs
  cmd <- paste("R CMD BATCH '--args jobnum=", blocks$jobnum[ii], " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
               " group.ind=", group[ii], " substitute=FALSE sub.betas=FALSE counterf.betas=FALSE",
               " s.epic=", cc,  " s.demog=", cc,
               " s.bmb=", cc, " s.bfb=", cc, # country to substitute in for each input (if doing that)
               " s.bme=", cc, " s.bfe=", cc,
               " s.bmp=", cc, " s.bfp=", cc, 
               " death=TRUE",
               " acute.sc=", blocks$acute.sc[ii], " late.sc=", blocks$late.sc[ii]," aids.sc=", blocks$aids.sc[ii], # acute phase varying throughout loop
               " dur.ac=", blocks$dur.ac[ii], " dur.lt=", blocks$dur.lt[ii], " dur.aids=", blocks$dur.aids[ii],
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
  num.doing <- num.doing+1
  cat(cmd)               # add command
  cat('\n')              # add new line
                                        #         }
}
sink()
head(blocks)
totn
save(blocks, file = file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
####################################################################################################
print(totn)
print(num.doing)
head(blocks[to.do,])
