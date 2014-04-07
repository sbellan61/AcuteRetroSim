library(plyr); library(data.table); library(abind); library(multicore)
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')
if(!file.exists(outdir)) dir.create(outdir)
load(file.path('results','RakAcute','blocks.Rdata')) # these are country-acute phase specific blocks
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation


## ##################################################################################################
## Load fit files from Hollingsworth-model-generated data
## ##################################################################################################
## Get file names
fls.holl <- list.files(file.path('results','RakAcute','HollModTestFits','fitouts'), pattern = 'fitout-', full.names=T)
length(fls.holl)
jobnums.holl <-  as.numeric(sapply(fls.holl, function(x) as.numeric(strsplit(x,'-')[[1]][2])))
fls.holl <- fls.holl[order(jobnums.holl)]
jobnums.holl <- jobnums.holl[order(jobnums.holl)]

## Load & label files by jobnum (fitout1, fitout2, etc...)
system.time(holldat <- mclapply(fls.holl, function(x) { jobnum <- as.numeric(strsplit(x,'-')[[1]][2]); load(x)
                                    assign(paste0('fitout',jobnum),fitout, env=.GlobalEnv)}))

## Each list element is one simulation output from a job
names(holldat[[1]])
holldat[[1]]$wtab ## raw simulation data

## Sample Sizes
wtabs.holl <- t(abind(lapply(holldat, function(x) { c(x$wtab$inct[1,'n'],x$wtab$prevt[1,'n'],x$wtab$latet[nrow(x$wtab$latet),'n'])}), along = 2))
wtabs.holl <- as.data.frame(wtabs.holl)
colnames(wtabs.holl) <- c('incn','prevn','laten')
wtabs.holl$job <- jobnums.holl
head(wtabs.holl)
apply(wtabs.holl[,-4],2, summary)

holldat[[1]]$outtab.h
## Hollingsworth Fits to Holl-gen data
hfits <- abind(lapply(holldat, function(x) { x$outtab.h}), along = 3)
dimnames(hfits)[[3]] <- jobnums.holl
dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(holldat, function(x) {as.data.frame(t(x$dpars))})))
cohvars$job <- jobnums.holl
## Pars
dpars <- blocks[jobnums.holl,]

## turn these into easy to work with data frames
hfh <- adply(hfits, 3:2, .parallel=T)
colnames(hfh)[1:2] <- c('job','var')
hfh$het.sd <- dpars$het.gen.sd[match(hfh$job, dpars$job)]
hfh <- merge(hfh, cohvars, by = 'job')
hfh <- hfh[order(hfh$var,hfh$job),]
head(hfh)
hfh <- within(hfh, {ehm.late <- (late.sc - 1)*dur.lt})
hfh <- within(hfh, {ehm.acute <- (acute.sc - 1)*dur.ac})

## ##################################################################################################
## Load fit files from SB-generated data
## ##################################################################################################
## Get file names
fls <- list.files(file.path('results','RakAcute','UgandaFits','fitouts'), pattern = 'fitout-', full.names=T)
length(fls)
jobnums <-  as.numeric(sapply(fls, function(x) as.numeric(strsplit(x,'-')[[1]][2])))
fls <- fls[order(jobnums)] ## sort on job number
jobnums <- jobnums[order(jobnums)]

## Load & label files by jobnum (fitout1, fitout2, etc...)
system.time(alldat <- mclapply(fls, function(x) { jobnum <- as.numeric(strsplit(x,'-')[[1]][2]); load(x)
                                    assign(paste0('fitout',jobnum),fitout, env=.GlobalEnv)}))

## Each list element is one simulation output from a job
names(alldat[[1]])
table(unlist(lapply(alldat, function(x) x$decont))) ## decontaminating between phases?
table(unlist(lapply(alldat, function(x) x$excl.extram)))## excluding extramarital transmission?
alldat[[1]]$wtab ## raw simulation data
       
## Sample Sizes
wtabs <- t(abind(lapply(alldat, function(x) { c(x$wtab$inct[1,'n'],x$wtab$prevt[1,'n'],x$wtab$latet[nrow(x$wtab$latet),'n'])}), along = 2))
wtabs <- as.data.frame(wtabs)
colnames(wtabs) <- c('incn','prevn','laten')
wtabs$job <- jobnums
head(wtabs)
apply(wtabs[,-4],2, summary)

## ## Hollingsworth Fits
## hfits <- abind(lapply(alldat, function(x) { x$outtab.h}), along = 4)
## dimnames(hfits)[[4]] <- jobnums
## dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## dim(alldat[[26]]$outtab.h)
## dim(alldat[[25]]$outtab.h)

## Hollingsworth Fits (without xbyerr fits) (REMOVE WHEN FIXED)
hfits <- abind(lapply(alldat, function(x) {
                                if(!is.na(dim(x$outtab.h)[3])) return(x$outtab.h[,,'base'])
                                if(is.na(dim(x$outtab.h)[3])) return(x$outtab.h)
                                 }), along = 3)
dimnames(hfits)[[3]] <- jobnums
dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## Wawer Fits (Poisson Regressions)
wfits <- abind(lapply(alldat, function(x) { x$outtab.w}), along = 6)
dimnames(wfits)[[6]] <- jobnums
#dimnames(wfits)
dimnames(wfits)[[1]] <- c('lci','med','uci','true')


## Empirical Wawer Hazards (no model fit)
ewfits <- t(abind(lapply(alldat, function(x) {x$erhs.w}), along = 2))
rownames(ewfits) <- jobnums

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(alldat, function(x) {as.data.frame(t(x$dpars))})))
cohvars$job <- jobnums
## Pars
dpars <- blocks[jobnums,]

## turn these into easy to work with data frames
hf <- adply(hfits, 4:2, .parallel=T)
colnames(hf)[1:3] <- c('job','err','var')
hf$het.sd <- dpars$het.gen.sd[match(hf$job, dpars$job)]
hf <- merge(hf, cohvars, by = 'job')
hf <- hf[order(hf$err,hf$var,hf$job),]
head(hf)


## turn these into easy to work with data frames (with xbyerr missing)
hf <- adply(hfits, 3:2, .parallel=T)
colnames(hf)[1:2] <- c('job','var')
hf$het.sd <- dpars$het.gen.sd[match(hf$job, dpars$job)]
hf <- merge(hf, cohvars, by = 'job')
hf <- hf[order(hf$err,hf$var,hf$job),]
head(hf)

## ## this can take a long time because each simulation had 24model fits
## wf <- adply(wfits, 6:2, .parallel=T)
## colnames(wf)[1:5] <- c('job','cov','hobs','err','var')
## wf$het.sd <- dpars$het.gen.sd[match(wf$job, dpars$job)]
## wf <- merge(wf, cohvars, by = 'job')
## wf <- wf[order(wf$err,wf$var,wf$job),]
## save(wf, hf, file=file.path(outdir, 'wf.Rdata'))
