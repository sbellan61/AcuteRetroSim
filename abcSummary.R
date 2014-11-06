setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
rm(list=ls()); gc()
sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
library(parallel)
ncores <- 12
fig.dir <- file.path('FiguresAndTables','abcFig')
out.dir <- file.path('results','abcSummary')
if(!file.exists(fig.dir))      dir.create(fig.dir) # create directory if necessary
if(!file.exists(out.dir))      dir.create(out.dir) # create directory if necessary
fls <- list.files('results/abcBatch1Old', pattern='Rdata', full.names=T)

pmatLs <- mclapply(fls, collectSumStat, mc.cores=3) ## not enough mem to do more than a few cores
pmat <- do.call(rbind.data.frame, lapply(pmatLs, '[[', 'pmat'))
pmat <- pmat[order(pmat$gVals),]
head(pmat[,c('acute.sc','dur.ac','het.gen.sd','enoughHet','gVals')],80)
save(pmat, pmatLs, file = file.path(out.dir,'pmatLs.Rdata'))

pdf(file.path(fig.dir, 'hist gval.pdf'))
hist(pmat$gVals, breaks = 0:max(pmat$gVals+1), col = 'black', xlab = 'G stat', main='')
graphics.off()

sbpairs(pmat[pmat$gVals<2, c('acute.sc','dur.ac','het.gen.sd','gVals')], file.path(fig.dir, 'Post1'), do.jpeg=T)
sbpairs(pmat[, c('acute.sc','dur.ac','het.gen.sd','gVals')], file.path(fig.dir, 'Prior'), do.jpeg=T)

mostFrq <- names(table(pmat$acute.sc)[which.max(table(pmat$acute.sc))])

xtabs(~filenm, pmat[pmat$acute.sc==mostFrq,])
head(pmat)

load(fls[1])
wtab <- sbmod.to.wdat(rcohsList[[1]]$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=T, giveProp=T)
contTabsSim <- mclapply(list(wtab), function(x) {
    x <- within(x, {
        inct$ni <- with(inct, n-i)
        inct <- inct[,c('i','ni')]
        prevt$ni <- with(prevt, n-i)
        prevt <- prevt[,c('i','ni')]})}, mc.cores=12)
ct <- contTabsSim[[1]]
ct$inct[3,1] <- 9

gSumStat(ct)
