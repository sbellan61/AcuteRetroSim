setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
rm(list=ls()); gc()
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
library(parallel)
ncores <- 12
fig.dir <- file.path('FiguresAndTables','abcFig')
out.dir <- file.path('results','abcSummary')
if(!file.exists(fig.dir))      dir.create(fig.dir) # create directory if necessary
if(!file.exists(out.dir))      dir.create(out.dir) # create directory if necessary
fls <- list.files('results/abcBatch1', pattern='Rdata', full.names=T)
#fls <- list.files('results/testDir', pattern='Rdata', full.names=T)

##tst <- collectSumStat(fls[2], browse=T, returnGtable=T, ncores=1)
#pmatLs <- lapply(fls, collectSumStat, returnGtable=T) ## not enough mem to do more than a few cores

sel <- c('acute.sc','dur.ac','het.gen.sd','gVals')
pmatLs <- mclapply(fls, collectSumStat, returnGtable=T, mc.cores=12, ncores = 1) ## ncores goes to collectSumStat, use each core to load rather
pmat <- do.call(rbind.data.frame, lapply(pmatLs, '[[', 'pmat'))
gtabs <- do.call(c, lapply(pmatLs, '[[', 'Gtable'))
ord <- order(pmat$gVals)
pmat <- pmat[ord,]
gtabs <- gtabs[ord]
save(gtabs,pmat, pmatLs, file = file.path(out.dir,'pmatLs.Rdata'))

## Look at Gtables for the best ones to make sure they seem to be good matchest
head(pmat[,c('acute.sc','dur.ac','het.gen.sd','enoughHet','gVals')],10)
gtabs[600]
pmat[600,]

load(file = file.path(out.dir,'pmatLs.Rdata'))

pdf(file.path(fig.dir, 'hist gval.pdf'))
hist(pmat$gVals, breaks = 0:max(pmat$gVals+1), col = 'black', xlab = 'G stat', main='')
graphics.off()

nrow(pmat)
cutf <- 4
sum(pmat$gVals<cutf)

pmat <- within(pmat, {logacute.sc <- log(acute.sc)})
sel <- c('logacute.sc','dur.ac','het.gen.sd','gVals','RHreduction')
rgs <- apply(pmat[,sel],2,range)
rgs[,'RHreduction'] <- c(.01,10)

sbpairs(pmat[with(pmat, gVals<cutf & enoughHet), sel], file.path(fig.dir, 'Post1'), do.jpeg=T, rgs=rgs)
sbpairs(pmat[, sel], file.path(fig.dir, 'Prior'), do.jpeg=T, rgs=rgs)
graphics.off()


sel <- c('bmp','bfp','acute.sc','dur.ac','het.gen.sd','gVals')
rgs <- apply(pmat[,sel],2,range)
sbpairs(pmat[with(pmat, gVals<cutf & enoughHet), sel], file.path(fig.dir, 'Post1haz'), do.jpeg=T, rgs=rgs)
sbpairs(pmat[, sel], file.path(fig.dir, 'Priorhaz'), do.jpeg=T, rgs=rgs)
graphics.off()

apply(pmat[,sel], 2, function(x) quantile(x,c(.025,.5,.975)))
apply(pmat[with(pmat, gVals<cutf & enoughHet),sel], 2, function(x) quantile(x,c(.025,.5,.975)))

cutf <- 2
pmatChosen <- pmat[with(pmat, gVals<cutf & enoughHet), parnms]
dim(pmatChosen)
head(pmatChosen)
ehms <- with(pmatChosen, (acute.sc-1)*dur.ac)
quantile(ehms,c(.025, .5, .975))
ehmsprior <- with(pmat, (acute.sc-1)*dur.ac)
quantile(ehmsprior,c(.025, .5, .975))
head(simParmSamp(parms=pmatChosen))

perturbParticle(pmatChosen[1,], sds)
