library(plyr); library(data.table); library(abind)
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')
if(!file.exists(outdir)) dir.create(outdir)
load(file.path('results','RakAcute','blocks.Rdata')) # these are country-acute phase specific blocks
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation

## Load fit files
fls <- list.files(file.path('results','RakAcute','UgandaFits'), pattern = 'fitout.Rdata', full.names=T, recursive=T)
length(fls)
jobnums <-  as.numeric(gsub("\\D", "", fls))
fls <- fls[order(jobnums)]
jobnums <- jobnums[order(jobnums)]

## Load & label files by jobnum (fitout1, fitout2, etc...)
alldat <- lapply(fls, function(x) { jobnum <- as.numeric(gsub("\\D", "", x)); load(x); assign(paste0('fitout',jobnum),fitout, env=.GlobalEnv)})

names(alldat[[1]])

## Hollingsworth Fits
hfits <- abind(lapply(alldat, function(x) { x$outtab.h}), along = 4)
dimnames(hfits)[[4]] <- jobnums
dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## Wawer Fits (Poisson Regressions)
wfits <- abind(lapply(alldat, function(x) { x$outtab.w}), along = 5)
dimnames(wfits)[[5]] <- jobnums
dimnames(wfits)
dimnames(wfits)[[1]] <- c('lci','med','uci','true')

## Empirical Wawer Hazards (no model fit)
ewfits <- t(abind(lapply(alldat, function(x) {x$erhs.w}), along = 2))
rownames(ewfits) <- jobnums

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(alldat, function(x) {as.data.frame(x[4:12])})))
cohvars$jobnum <- jobnums

## Pars
dpars <- blocks[jobnums,]

sel <- which(jobnums < 100)

hf <- adply(hfits, 4:2)
colnames(hf)[1:3] <- c('job','err','var')
hf$het.sd <- dpars$het.gen.sd[match(hf$job, dpars$job)]
head(hf)

cols <- rainbow(11)
pdf(file.path(outdir,'test.pdf'))
par(mfrow=c(1,3))
maxs <- c(100,60,10)
for(vv in 1:3) {
    var <- c('ehm.ac','acute.sc', 'dur.ac')[vv]
    plot(0,0, type='n', xlim = c(0, maxs[vv]), ylim = c(0,maxs[vv]), bty = 'n', xlab = 'EHM',
         ylab = expression(paste('estimated ',EHM[acute])), main = var)
    sel <- hf$var==var & hf$err=='base' & hf$het.sd==0
    ddply(hf[sel,], .(het.sd), with, {polygon(c(true,rev(true)), c(lci, rev(uci)), col = cols[het.sd+1])})
    ddply(hf[sel,], .(het.sd), with, {lines(true, med)}) #, col = cols[het.sd+1])})
    abline(a=0, b=1, col = 'black', lty = 3)
}
dev.off()

ddply(hf[hf$var=='ehm.ltaids',], .(var, het.sd, err), summarise, mean(med))
   
wf <- adply(wfits, 5:2)
colnames(wf)[1:4] <- c('job','cov','err','var')
wf$het.sd <- dpars$het.gen.sd[match(wf$job, dpars$job)]
head(wf)
ddply(wf[wf$var=='ehm.ltaids' & wf$cov=='0obs',], .(var, het.sd, err), summarise, mean(med))

ddply(wf[wf$var=='ehm.ac',], .(het.sd, err, cov), summarise, mean(med))
