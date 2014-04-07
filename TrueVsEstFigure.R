library(plyr); library(data.table); library(abind); library(multicore)
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')
nc <- 12                                       # core per simulation
load(file=file.path(outdir, 'wf.Rdata'))
load(file = file.path('results','HollingsworthAn','RealExclbyErr','workspace.Rdata')) ## fit to real data

west <- 36.25
hest <- 65.4

####################################################################################################
## True vs est 4x4 plot
pdf(file.path(outdir, 'TrueVsEst.pdf'), w = 6.83, h = 6)
ct <- 12
madj <- .5
mln <- 2
mlnv <- 3
ptcol <- gray(.3)
arcol <- gray(.5)
xmax <- 100
ymax <- 110
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var & hf$ehm.late==40,] 
twf <- wf[wf$var==var & wf$ehm.late==40 & wf$cov=='',]
## Layout 4x4 plot
layout(matrix(c(1,11,7,2,12,8,3,5,9,4,6,10), 3, 4))
par(mar = c(3.2,3,1,.5), ps = ct, oma = c(1,6,3,0))
## 1) Waw, correct, het.gen = 0
sel <-  twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0'
temp <- twf[sel,]
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = expression(sigma[hazard]==0), las = 1, axes=F)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 2)
mtext('full inclusion', side = 3, adj = madj, line = mln)
mtext('Wawer:\nno covariates', side = 2, adj = madj, line = mlnv)
with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
abline(a=0, b=1, col = 'black', lty = 1)
segments(0, west, xmax, west, col = 'black', lty = 3, lwd = 2)
## 2-4) Waw, err, het.gen = 0,1,2
for(hh in 0:2) {
  sel <-  twf$err=='XbErr' & twf$het.sd==hh & twf$hobs=='obs0'
  temp <- twf[sel,]
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = bquote(sigma[hazard]==.(hh)), las = 1, axes = F)
  axis(2, seq(0,100, by = 20), las = 2)
  axis(1, seq(0,100, by = 20), las = 2)
  with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
  with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
  abline(a=0, b=1, col = 'black', lty = 1)
  segments(0, west, xmax, west, col = 'black', lty = 3, lwd = 2)
  if(hh==1) mtext('excluded incident SDCs lost to follow-up', side = 3, adj = madj, line = mln)
}
## 5-6) Waw, err, het.gen = 1,2, with 50% variance controlled
for(hh in 1:2) {
  sel <-  twf$err=='XbErr' & twf$het.sd==hh & twf$hobs=='obs0.7'
  temp <- twf[sel,]
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes = F)
  if(hh==1) mtext('Wawer:\n50% variance \ncontrolled \nthrough covariates', side = 2, adj = madj, line = mlnv+3)
  axis(2, seq(0,100, by = 20), las = 2)
  axis(1, seq(0,100, by = 20), las = 2)
  with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
  with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
  abline(a=0, b=1, col = 'black', lty = 1)
  segments(0, west, xmax, west, col = 'black', lty = 3, lwd = 2)
}
## 7) Holl, correct, het.gen = 0
sel <-  thf$err=='base' & thf$het.sd==0
temp <- thf[sel,]
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes = F)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 2)
with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
abline(a=0, b=1, col = 'black', lty = 1)
segments(0, hest, xmax, hest, col = 'black', lty = 3, lwd = 2)
mtext('Hollingsworth:\nno covariates', side = 2, adj = madj, line = mlnv)
## 8-10) Holl, err, het.gen = 0,1,2
for(hh in 0:2) {
  sel <-  thf$err=='XbErr' & thf$het.sd==hh
  temp <- thf[sel,]
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes = F)
  axis(2, seq(0,100, by = 20), las = 2)
  axis(1, seq(0,100, by = 20), las = 2)
  with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
  with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
  abline(a=0, b=1, col = 'black', lty = 1)
  segments(0, hest, xmax, hest, col = 'black', lty = 3, lwd = 2)
}
mtext(expression(paste('estimated ',EHM[acute])), side = 2, line = 4, adj = .5, outer = T, ps = 12)
mtext(expression(paste('true (simulated) ',EHM[acute])), side = 1, line = 0, adj = .5, outer = T, ps = 12)
dev.off()

####################################################################################################
## function to plot LOESS for a subset
ploes <- function(temp, err='XbErr', hsd=0, ehm.lt=40, hobs=NULL, x = 1:100, col='black', lty = 1, lwd = 1, ...) {
  if(is.null(hobs))     sel <- with(temp, err==err & het.sd==hsd & ehm.late==ehm.lt)
  if(!is.null(hobs))     sel <- with(temp, err==err & het.sd==hsd & hobs==hobs  & ehm.late==ehm.lt)
  with(temp[sel,], lines(x, predict(loess(med~true, ...), x), col = col, lty = lty, lwd = lwd))
}

####################################################################################################
## True vs est 4x4 plot: 2 panels with loess lines
##############################
pdf(file.path(outdir, 'TrueVsEst LOESS.pdf'), w = 6.83, h = 3)
hsds <- unique(c(twf$het.sd, thf$het.sd))
cols <- colorRampPalette(c('purple','orange','red'))(length(hsds))
## cols <- rainbow(length(hsds))
ct <- 8
par(mfrow=c(1,2), mar = c(3.2,3,1,.5), pointsize = ct, oma = c(0,0,0,0))#, cex.main = ct, cex.lab = ct, cex.axis = ct)
cex.leg <- .6
madj <- .5
mln <- 2
mlnv <- 3
xmax <- 100
ymax <- 100
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var & hf$ehm.late==40,] 
twf <- wf[wf$var==var & wf$ehm.late==40 & wf$cov=='',]
## Plot layout
## Wawer plot
par('ps'=ct)
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '(A) Wawer Model', axes = F)
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
ploes(twf, err='base', hsd=0, hobs='obs0', lty = 2, col=cols[1]) ## Base, het=0, Wawer
## Err, het=0:3, obs = 0,5 Wawer
for(hh in 1:length(hsds)) {
  for(bb in 1) { ## options to show multiple obs lines, do this on a separate figure though
    hobs <- c('obs0','obs0.7')[bb]
    ploes(twf, err='XbErr', hobs=hobs, hsd=hsds[hh], col =cols[hh], lty = bb)
  }}
##################################################
## Holl plot
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '(B) Hollingsworth Model', axes = F)
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
ploes(thf, err='base', hsd=0, lty = 2, col=cols[1]) ## Base, het=0, Wawer
## Err, het=0:3, obs = 0,5 Holl
hsds <- unique(thf$het.sd) ## seq(0,3,by=.5)
for(hh in 1:length(hsds)) {
  for(bb in 1) { ## options to show multiple obs lines, do this on a separate figure though
    hobs <- c('obs0','obs0.7')[bb]
    ploes(thf, err='XbErr', hsd=hsds[hh], col =cols[hh], lty = bb)
  }}
legend('bottomright', leg = hsds, col = cols, lty = 1, bty = 'n', cex = cex.leg, title = 'degree of heterogeneity\nsigma =')
legend('bottom', leg = c('included','excluded'), lty = 2:1, cex = cex.leg, bty = 'n', title='incident SDC \nlost to follow-up')
##################################################
mtext(expression(paste('estimated ',EHM[acute])), side = 2, line = -1, adj = .5, outer = T)#, ps = 12)
mtext(expression(paste('true (simulated) ',EHM[acute])), side = 1, line = -1, adj = .5, outer = T)#, ps = 12)
graphics.off()
 
####################################################################################################
## Show how controlling for error affects things
##############################
pdf(file.path(outdir, 'TrueVsEst LOESS het obs.pdf'), w = 3.27, h = 3)
## cols <- rainbow(length(hsds))
hsd <- 2 ## for this plot pick a het.gen
cov <- c(0,.5,.7,.9)
hobs <- paste0('obs', cov) #seq(0,1, by = .1))
cols <- colorRampPalette(c('purple','orange'))(length(hobs))
ct <- 8
par(mfrow=c(1,1), mar = c(3.2,3,1,.5), pointsize = ct, oma = c(0,0,0,0))#, cex.main = ct, cex.lab = ct, cex.axis = ct)
cex.leg <- .8
madj <- .5
mln <- 2
mlnv <- 3
xmax <- 100
ymax <- 100
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var & hf$ehm.late==40,] 
twf <- wf[wf$var==var & wf$ehm.late==40 & wf$cov=='',]
## Plot layout
## Wawer plot
par('ps'=ct)
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', axes = F)
     ## main = 'Removing Heterogeneity by \nControlling for Measured Confounders',
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
ploes(twf, err='XbErr', hobs='obs0', hsd=0, col =cols[1], lty = 2)
## Err, het=0:3, obs = 0,5 Wawer
  for(bb in 1:length(hobs)) { ## options to show multiple obs lines, do this on a separate figure though
    ploes(twf, err='XbErr', hobs=hobs[bb], hsd=hsd, col =cols[bb], lty = 1)
  }
legend('bottomright', leg = cov^2, col = cols, lty = 1, bty = 'n', cex = cex.leg, ncol=2, title = 'proportion variance \ncontrolled for \nvia covariates')
legend('topleft', leg = c('0','2'), lty = 2:1, cex = cex.leg, bty = 'n', title='sigma')
##################################################
mtext(expression(paste('estimated ',EHM[acute])), side = 2, line = -1, adj = .5, outer = T)#, ps = 12)
mtext(expression(paste('true (simulated) ',EHM[acute])), side = 1, line = -1, adj = .5, outer = T)#, ps = 12)
graphics.off()

####################################################################################################
## Show how late phase affects things
##############################
pdf(file.path(outdir, 'TrueVsEst LOESS by ehmlate.pdf'), w = 6.83, h = 3)
## cols <- rainbow(length(hsds))
ehm.lates <- unique(c(hf$ehm.late,wf$ehm.late))
ehm.lates <- ehm.lates[order(ehm.lates)]
cols <- colorRampPalette(c('purple','orange'))(length(ehm.lates))
ct <- 8
par(mfrow=c(1,1), mar = c(3.2,3,1,.5), pointsize = ct, oma = c(0,0,0,0))#, cex.main = ct, cex.lab = ct, cex.axis = ct)
cex.leg <- .8
madj <- .5
mln <- 2
mlnv <- 3
xmax <- 100
ymax <- 100
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var & hf$err=='base' ,] 
twf <- wf[with(wf, var==var & err=='base' & cov=='' & hobs=='obs0'),]
## Plot layout
par('ps'=ct, mfrow=c(1,2))
## Wawer plot
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', axes = F)
     ## main = 'Removing Heterogeneity by \nControlling for Measured Confounders',
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
## Wawer, Err, het = 0, late=0,10,40,90
  for(bb in 1:length(ehm.lates)) { ## options to show multiple obs lines, do this on a separate figure though
    ploes(twf, err='base', hobs='obs0', hsd=0, ehm.lt=ehm.lates[bb], col =cols[bb], lty = 1)
  }
##########
## Holl plot
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', axes = F)
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
## Holl, Err, het = 0, late=0,10,40,90
  for(bb in 1:length(ehm.lates)) { ## options to show multiple obs lines, do this on a separate figure though
    ploes(thf, err='base', hsd=0, ehm.lt=ehm.lates[bb], col =cols[bb], lty = 1)
  }
legend('bottomright', leg = ehm.lates, col = cols, lty = 1, bty = 'n', cex = cex.leg, ncol=2, title = expression(EHM['late']))
##################################################
mtext(expression(paste('estimated ',EHM[acute])), side = 2, line = -1, adj = .5, outer = T)#, ps = 12)
mtext(expression(paste('true (simulated) ',EHM[acute])), side = 1, line = -1, adj = .5, outer = T)#, ps = 12)
graphics.off()


####################################################################################################
## Get fake CI's for the real estimates using simulation

ehmacp <- fout$exposts$atr.month.ac
## qehms <- quantile(ehmacp, c(.025,seq(.05,.95,by=.05), .975)) ## posterior quantiles
ehm.q <- ecdf(ehmacp) ## empirical quantile function
## Assign a quantile to every simulation results
twf <- wf[with(wf, var=='ehm.ac' & err=='base' & cov=='' & hobs=='obs0' & ehm.late==40),]
head(twf)
twf <- ddply(twf, .(), transform, rl.qnt = ehm.q(med))
pal <- colorRamp(c('purple','orange','red'))
twf$col <- apply(pal(abs(twf$rl.qnt-.5)*2), 1, function(x) rgb(x[1],x[2],x[3], alpha = 100, maxColorValue=255))

pdf(file.path(outdir, 'fake cis.pdf'), w = 3.27, h = 3)
ct <- 8
par('ps'=ct, mar=c(3,3,1,1))
with(twf, plot(jitter(het.sd, a = .2), true, col = col, xlim = c(-.2,3), ylim = c(0,100), pch = 19, cex = .2,
               main='true values consistent with real estimates'))
mtext(expression(paste('true ',EHM[acute])), side = 2, line = -1, adj = .5, outer = T)#, ps = 12)
mtext(expression(paste(sigma['hazard'])), side = 1, line = -1, adj = .5, outer = T)#, ps = 12)
dev.off()


