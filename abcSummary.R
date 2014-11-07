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

sel <- c('acute.sc','dur.ac','het.gen.sd','gVals')
pmatLs <- mclapply(fls, collectSumStat, returnGtable=T, mc.cores=8) ## not enough mem to do more than a few cores
pmat <- do.call(rbind.data.frame, lapply(pmatLs, '[[', 'pmat'))
gtabs <- do.call(c, lapply(pmatLs, '[[', 'Gtable'))
pmat <- pmat[order(pmat$gVals),]
gtabs <- gtabs[order(pmat$gVals)]
head(pmat[,c('acute.sc','dur.ac','het.gen.sd','enoughHet','gVals')],20)
gtabs[1:2]
save(gtabs,pmat, pmatLs, file = file.path(out.dir,'pmatLs.Rdata'))

sum(gtabs[1][[1]][[1]][,5,])
head(pmat)

load(file = file.path(out.dir,'pmatLs.Rdata'))

head(pmatLs[[1]]$Gtable,20)

length(gtabs)
dim(pmat)

pdf(file.path(fig.dir, 'hist gval.pdf'))
hist(pmat$gVals, breaks = 0:max(pmat$gVals+1), col = 'black', xlab = 'G stat', main='')
graphics.off()

nrow(pmat)
cutf <- 10
sum(pmat$gVals<cutf)

sel <- c('acute.sc','dur.ac','het.gen.sd','gVals','RHreduction')
rgs <- apply(pmat[,sel],2,range)
rgs[,'RHreduction'] <- c(.2,10)
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
