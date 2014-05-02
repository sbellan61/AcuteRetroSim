rm(list=ls(all=T)); gc()
library(plyr); library(data.table); library(abind); library(multicore)
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','HollModTestFitsSummary')
if(!file.exists(outdir)) dir.create(outdir)
load(file = file.path('results','RakAcute','HollModTestFits', 'blocksHollFitTest.Rdata'))
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation

## Load fit files
fls <- list.files(file.path('results','RakAcute','HollModTestFits','fitouts'), pattern = 'fitout-', full.names=T)
length(fls)
jobnums <-  as.numeric(sapply(fls, function(x) as.numeric(strsplit(x,'-')[[1]][2])))
fls <- fls[order(jobnums)]
jobnums <- jobnums[order(jobnums)]

## Load & label files by jobnum (fitout1, fitout2, etc...)
system.time(alldat <- mclapply(fls, function(x) { jobnum <- as.numeric(strsplit(x,'-')[[1]][2]); load(x)
                                    assign(paste0('fitout',jobnum),fitout, env=.GlobalEnv)}))

## Each list element is one simulation output from a job
names(alldat[[1]])
alldat[[1]]$wtab ## raw simulation data

## Sample Sizes
wtabs <- t(abind(lapply(alldat, function(x) { c(x$wtab$inct[1,'n'],x$wtab$prevt[1,'n'],x$wtab$latet[nrow(x$wtab$latet),'n'])}), along = 2))
wtabs <- as.data.frame(wtabs)
colnames(wtabs) <- c('incn','prevn','laten')
wtabs$job <- jobnums
head(wtabs)
apply(wtabs[,-4],2, summary)

alldat[[1]]$outtab.h
## Hollingsworth Fits to Holl-gen data
hfits <- abind(lapply(alldat, function(x) { x$outtab.h}), along = 3)
dimnames(hfits)[[3]] <- jobnums
dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(alldat, function(x) {as.data.frame(t(x$dpars))})))
cohvars$job <- jobnums
## Pars
dpars <- blocks[jobnums,]

## turn these into easy to work with data frames
hfh <- adply(hfits, 3:2, .parallel=T)
colnames(hfh)[1:2] <- c('job','var')
hfh$het.sd <- dpars$het.gen.sd[match(hfh$job, dpars$job)]
hfh <- merge(hfh, cohvars, by = 'job')
hfh <- hfh[order(hfh$var,hfh$job),]
head(hfh)
hfh <- within(hfh, {ehm.late <- (late.sc - 1)*dur.lt})
hfh <- within(hfh, {ehm.acute <- (acute.sc - 1)*dur.ac})

## Compare these to fits of SB cohort data
load(file=file.path('results','RakAcute','UgandaFitSummaries', 'wf.Rdata'))
 
## Hollingsworth sensitivity analysis: showing estimates the elevated hazard months of the acute phase
ehm.lates <- unique(c(hf$ehm.late,hfh$ehm.late))
ehm.lates <- ehm.lates[order(ehm.lates)] ## which ehmls exist
ehmls <- c(10,40,90) ## which EHMlates to show (for clearer plots)
ehmlas <- ehmls-10 ## -10 missing chronic months from AIDS phase
xmax <- 100
ymax <- 100
var <- 'ehm.ac'
cols <- c('purple','red','orange')
pdf(file.path(outdir,'HollAn HollMod ehm acute.pdf'), w = 3.27, h = 3.5)
par(mfrow=c(1,1), mar = c(3.5,3.5,1,.5), 'ps'=10, mgp = c(2.2,1,0))
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = expression(paste('true ',EHM[acute])),
     ylab = expression(paste('estimated ',EHM[acute])), main = '', asp=1)
abline(a=0, b=1, col = gray(.4), lty = 1, lwd=2)
for(jj in 1:2) { ## for Holl,SB generated data
  dd <- c('hf','hfh')[jj]
  temp <- get(dd)
  temp <- temp[temp$ehm.late %in% ehmls,]
  temp <- temp[rev(order(temp$ehm.late)),]
  ## with(temp[sel,],  arrows(true, lci, true,uci, length=.05, angle = 90, code = 3,col = cols[ehm.late+1]), lwd = .7)
  ## with(temp[sel,],  points(true, med, col = cols[ehm.late+1], cex = .5, pch = 19))
  for(ee in length(ehmls):1) {
      if(dd=='hfh') sel <- temp$var==var & temp$ehm.late==ehmls[ee]
      if(dd=='hf') sel <- temp$var==var & temp$het.sd==0 & temp$err=='base' & temp$ehm.late==ehmls[ee]
      lines(with(temp[sel,], lowess(med ~ true)), col = cols[ee], lwd = 1, lty = jj)
  }
  ## segments(1, 75.4, xmax, 75.4, col = 'black', lty = 3, lwd = 2)
  ## text(30, 72.4, 'Hollingsworth Estimate', pos = 3, cex = ct)
}
legend(0,ymax, title= expression(paste(EHM[late+AIDS])), ncol=1,
       leg = ehmlas, lty = 1, col=cols, bty ='n', cex = .7)
legend('bottomright', leg = c('Hollingsworth', 'Bellan'), title='Data-Generation Model', lty = 1:2, cex = .7, bty = 'n')
dev.off()


## Hollingsworth sensitivity analysis: showing estimates the elevated hazard months of the late phase
xmax <- 100
ymax <- 150
cols <- colorRampPalette(c('purple','orange'))(max(hf$acute.sc))
var <- 'ehm.lt'
pdf(file.path(outdir,'HollAn HollMod ehm late.pdf'), w = 6.83, h = 3.5)
par(mfrow=c(1,2), mar = c(3.2,3,2,.5), cex.main = ct, cex.lab = ct*1.2, cex.axis = ct, oma = c(.3,0,0,0))
mains <- paste('Data from', c('Hollingsworth', 'SB'), '\nData-Generation Model')
for(jj in 1:2) { ## for Holl,SB generated data
  dd <- c('hfh','hf')[jj]
  temp <- get(dd)
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = expression(paste('true ',EHM[late])),
       ylab = expression(paste('estimated ',EHM[late])), main = mains[jj])
  if(dd=='hfh') sel <- temp$var==var
  if(dd=='hf') sel <- temp$var==var & temp$het.sd==0 & temp$err=='base'
  with(temp[sel,],  arrows(true+acute.sc/10, lci, true+acute.sc/10,uci, length=.05, angle = 90, code = 3,col = cols[acute.sc]))
  with(temp[sel,],  points(true+acute.sc/10, med, col = cols[acute.sc], cex = .7, pch = 19))
  ddply(temp[sel,], .(acute.sc), function(x) {lines(lowess(x$med ~ x$true), col = cols[x$acute.sc][1], lwd = 1)})
  abline(a=0, b=1, col = 'black', lty = 1)
  segments(0, (76/10.6-1)*9, xmax*.8, (76/10.6-1)*9, col = 'red', lty = 2, lwd = 2)
  text(xmax*.8, (76/10.6-1)*9, 'Hollingsworth \nEstimate', pos = 4, cex = .5)
  acute.scs <- unique(temp$acute.sc)[order(unique(temp$acute.sc))]
  if(dd=='hfh') legend(5,ymax, title= 'acute to chronic relative hazard \n(for months 20-11 prior to death)', ncol = 3,
       leg =acute.scs, pch = 19, col=cols[acute.scs], bty ='n', cex = .7)
}
dev.off()

## Hollingsworth sensitivity analysis: showing estimates the elevated hazard months of the *LATE* phase
ehm.acutes <- unique(c(hf$ehm.acute,hfh$ehm.acute))
ehm.acutes <- ehm.acutes[order(ehm.acutes)] ## which ehmls exist
ehmas <- c(0,10,30,48) ## which EHMlates to show (for clearer plots)
xmax <- 100
ymax <- 100
var <- 'ehm.lt'
cols <- c('purple','red','blue','orange')
pdf(file.path(outdir,'HollAn HollMod ehm LATE.pdf'), w = 3.27, h = 3.5)
par(mfrow=c(1,1), mar = c(3.5,3.5,1,.5), 'ps'=10, mgp = c(2.2,1,0))
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = expression(paste('true ',EHM[late])),
     ylab = expression(paste('estimated ',EHM[late])), main = '', asp=1)
abline(a=0, b=1, col = gray(.4), lty = 1, lwd=2)
for(jj in 1:2) { ## for Holl,SB generated data
  dd <- c('hf','hfh')[jj]
  temp <- get(dd)
  temp <- temp[temp$ehm.acute %in% ehmas,]
  temp <- temp[rev(order(temp$ehm.acute)),]
  ## with(temp[sel,],  arrows(true, lci, true,uci, length=.05, angle = 90, code = 3,col = cols[ehm.late+1]), lwd = .7)
  ## with(temp[sel,],  points(true, med, col = cols[ehm.late+1], cex = .5, pch = 19))
  for(ee in length(ehmas):1) {
      if(dd=='hfh') sel <- temp$var==var & temp$ehm.acute==ehmas[ee]
      if(dd=='hf') sel <- temp$var==var & temp$het.sd==0 & temp$err=='base' & temp$ehm.acute==ehmas[ee]
      lines(with(temp[sel,], lowess(med ~ true)), col = cols[ee], lwd = 1, lty = jj)
  }
  ## segments(1, 75.4, xmax, 75.4, col = 'black', lty = 3, lwd = 2)
  ## text(30, 72.4, 'Hollingsworth Estimate', pos = 3, cex = ct)
}
legend(0,ymax, title= expression(paste(EHM[acute])), ncol=1,
       leg = ehmas, lty = 1, col=cols, bty ='n', cex = .7)
legend('bottomright', leg = c('Hollingsworth', 'Bellan'), title='Data-Generation Model', lty = 1:2, cex = .7, bty = 'n')
dev.off()



