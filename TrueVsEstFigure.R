library(plyr); library(data.table); library(abind); library(multicore)
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
## Creates Figures 3 & S4 from the manuscript
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')
nc <- 12                                       # core per simulation
load(file=file.path(outdir, 'wf.Rdata'))
load(file = file.path('results','HollingsworthAn','RealExclbyErr','workspace.Rdata')) ## fit to real data

west <- 36.25 ## Wawer median (if plotting horizontal line to show their estimate on the plots)
hest <- 65.4 ## Holl median (refit) (if plotting horizontal line to show their estimate on the plots)

####################################################################################################
## Figure S4
pdf(file.path(outdir, 'TrueVsEst.pdf'), w = 6.83, h = 4.5)
cols <- c('purple','red','blue','orange')
ct <- 12
mcex <- .8
madj <- .5
mln <- 2
mlnv <- 3
ptcol <- gray(.4)
arcol <- gray(.6)
xmax <- 100
ymax <- 110
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var ,] 
twf <- wf[wf$var==var & wf$cov=='',]
## Layout 4x4 plot
layout(t(matrix(1:18, 6, 3)))
par(mar = c(3.2,2.5,1,.1), ps = ct, oma = c(1,6,5.5,0))
## 1) Waw, base, het.gen = 0, no late phase
sel <-  twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==10
temp <- twf[sel,]
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes=F)
abline(a=0, b=1, col = 'black', lty = 1)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 2)
mtext('full inclusion', side = 3, adj = madj, line = mln, cex = mcex)
mtext('Wawer Model', side = 2, adj = madj, line = mlnv, cex = mcex)
with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
with(temp, lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 3, lwd = 1.6))
## 1) Waw, base, het.gen = 0, w/ late phase
sel <-  twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==40
temp <- twf[sel,]
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes=F)
mtext(bquote(sigma[hazard]==0), side = 3, line = .5, adj = .5)
abline(a=0, b=1, col = 'black', lty = 1)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), lab = NA)
with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
with(temp, lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 2, lwd = 1.6))
## segments(0, west, xmax, west, col = 'black', lty = 3, lwd = 2)
## 2-4) Waw, err, het.gen = 0,1,2,3
for(hh in 0:3) {
  sel <-  twf$err=='XbErr' & twf$het.sd==hh & twf$hobs=='obs0' & twf$ehm.late==40
  temp <- twf[sel,]
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '',
       main = '', las = 1, axes = F)
  abline(a=0, b=1, col = 'black', lty = 1)
  if(hh>0) mtext(bquote(sigma[hazard]==.(hh)), side = 3, line = .5, adj = .5)
  axis(2, seq(0,100, by = 20), lab = NA)
  axis(1, seq(0,100, by = 20), las = 2)
  with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
  with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
  with(temp, lines(1:100, predict(loess(med~true), 1:100), col = cols[hh+1], lty = 1, lwd = 1.3))
  ## segments(0, west, xmax, west, col = 'black', lty = 3, lwd = 2)
  if(hh==1) mtext('excluded incident SDCs lost to follow-up', side = 3, adj = madj, line = mln, cex = mcex)
}
mtext('elevated \nlate infectivity', side = 3, adj = .8, line = mln+1, cex = mcex, outer = T)
mtext('no elevated \nlate infectivity', side = 3, adj = .1, line = mln+1, cex = mcex, outer = T)
frame()
mtext('Wawer Model:\n50% variance \ncontrolled \nthrough covariates', side = 2, adj = madj, line = mlnv-3, cex = mcex)
frame()
frame()
## 5-6) Waw, err, het.gen = 1,2, with 50% variance controlled
for(hh in 1:3) {
  sel <-  twf$err=='XbErr' & twf$het.sd==hh & twf$hobs=='obs0.7' & twf$ehm.late==40
  temp <- twf[sel,]
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes = F)
  abline(a=0, b=1, col = 'black', lty = 1)
  if(hh==1) axis(2, seq(0,100, by = 20), las = 2) else axis(2, seq(0,100, by = 20), lab = NA)
  axis(1, seq(0,100, by = 20), las = 2)
  with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
  with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
  with(temp, lines(1:100, predict(loess(med~true), 1:100), col = cols[hh+1], lty = 3, lwd = 1.6))
  ## segments(0, west, xmax, west, col = 'black', lty = 3, lwd = 2)
}
## 7) Holl, base, het.gen = 0, no late
sel <-  thf$err=='base' & thf$het.sd==0 & thf$ehm.late==10
temp <- thf[sel,]
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes = F)
abline(a=0, b=1, col = 'black', lty = 1)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 2)
with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
with(temp, lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 3, lwd = 1.6))
mtext('Hollingsworth Model', side = 2, adj = madj, line = mlnv, cex = mcex)
## 7) Holl, base, het.gen = 0, with late
sel <-  thf$err=='base' & thf$het.sd==0 & thf$ehm.late==40
temp <- thf[sel,]
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes = F)
abline(a=0, b=1, col = 'black', lty = 1)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 2)
with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
with(temp, lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 2, lwd = 1.6))
## segments(0, hest, xmax, hest, col = 'black', lty = 3, lwd = 2)
## 8-10) Holl, err, het.gen = 0,1,2
for(hh in 0:3) {
  sel <-  thf$err=='XbErr' & thf$het.sd==hh
  temp <- thf[sel,]
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '', main = '', las = 1, axes = F)
  abline(a=0, b=1, col = 'black', lty = 1)
  axis(2, seq(0,100, by = 20), lab = NA)
  axis(1, seq(0,100, by = 20), las = 2)
  with(temp, arrows(true , lci, true ,uci, length=.02, angle = 90, code = 3, col = arcol))
  with(temp, points(true , med, cex = .7, pch = 19, col = ptcol))
  with(temp, lines(1:100, predict(loess(med~true), 1:100), col = cols[hh+1], lty = 1, lwd = 1.3))
  ## segments(0, hest, xmax, hest, col = 'black', lty = 3, lwd = 2)
}
mtext(expression(paste('estimated ',EHM[acute])), side = 2, line = 4, adj = .5, outer = T, ps = 12)
mtext(expression(paste('true (simulated) ',EHM[acute])), side = 1, line = 0, adj = .5, outer = T, ps = 12)
graphics.off()


####################################################################################################
## Figure 3
##############################
pdf(file.path(outdir, 'TrueVsEst LOESS.pdf'), w = 6.83, h = 4.5)
hsds <- 0:3 #unique(c(twf$het.sd, thf$het.sd))
## cols <- colorRampPalette(c('purple','red','orange'))(length(hsds))
## cols <- rainbow(length(hsds))
cols <- c('purple','red','blue','orange')
ct <- 10
par(mfrow=c(2,2), mar = c(4,4.5,2.5,1.5), mgp=c(2.5,1,0), 'ps'=ct, oma = c(0,0,0,8))#, cex.main = ct, cex.lab = ct, cex.axis = ct)
cex.leg <- 1
madj <- .5
mln <- 2
mlnv <- 3
####################################################################################################
## Log-hazard distribution figure
## pdf(file.path(outdir, 'log-hazard distributions.pdf'), w = 3, h = 3)
## par(mar = c(4,.5,1,.5), 'ps'=8)
plot(0,0, type = 'n', xlab = bquote(paste('multiple of ', lambda[hazard])), ylab='probability density', main='(A) Log-Hazard Distribution', bty = 'nn',
     xlim = c(10^-3,10^3), log='x',ylim = c(0,.6), axes=F) ##xaxt='n')
#axis(1, at = 10^(-3:3), lab = c(0.001, 0.01, 0.1, 1, 10, 100, 1000), las = 1)
labs <- c(expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0),expression(10^1), expression(10^2), expression(10^3))
axis(1, at = 10^(-3:3), lab = labs, las = 1)
axis(2, at = seq(0,.6, by = .2), las = 2)
for(hsd in 1:3) {
  curve(dnorm(log(x), 0,  hsd), from = 10^-3, to = 10^3, add = T, col = cols[hsd+1])
}
segments(1,0,1,1, col = cols[1])
legend('topleft', leg = 0:3, lty = 1, col = cols, title = expression(sigma[hazard]), bty = 'n', cex = cex.leg)
##################################################
## ehm plots
xmax <- 100
ymax <- 100
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var,] 
twf <- wf[wf$var==var &  wf$cov=='',]
## Plot layout
ylab <- expression(paste('estimated ',EHM[acute]))
xlab <- expression(paste('true (simulated) ',EHM[acute]))
## Wawer plot
par('ps'=ct)
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = xlab, ylab = ylab, main = '(B) Wawer Model', axes = F)
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
## seroincident ltf couples included EHM_late=0
sel <- twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==0
with(twf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 3, lwd = 1))
## seroincident ltf couples included
sel <- twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==40
with(twf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 2, lwd = 1))
## Err, het=0:3, obs = 0,5 Wawer
for(hh in 1:length(hsds)) {
  sel <- twf$err=='XbErr' & twf$het.sd==hsds[hh] & twf$hobs=='obs0' & twf$ehm.late==40
  with(twf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[hh], lty = 1, lwd = 1))
}
legend('bottomright', leg = c('included','excluded'), lty = 2:1, cex = cex.leg, title='incident SDC \nlost to follow-up', bty = 'n')
##################################################
## Holl plot
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = xlab, ylab = ylab, main = '(C) Hollingsworth Model', axes = F)
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
## seroincident ltf couples included EHM_late=0 
sel <- thf$err=='base' & thf$het.sd==0 & thf$ehm.late==0
with(thf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 3, lwd = 1))
## seroincident ltf couples included
sel <- thf$err=='base' & thf$het.sd==0 & thf$ehm.late==40
with(thf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 2, lwd = 1))
## Err, het=0:3, obs = 0,5 Holl
for(hh in 1:length(hsds)) {
  sel <- thf$err=='XbErr' & thf$het.sd==hsds[hh] & thf$ehm.late==40
  with(thf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[hh], lty = 1, lwd = 1))
}
####################################################################################################
## Show how controlling for error affects things Wawer model
cov <- c(0,.5,.7,.9)
hobs <- paste0('obs', cov) #seq(0,1, by = .1))
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = xlab, ylab = ylab, axes = F, main = '(D) Wawer Model (multivariate)')
     ## main = 'Removing Heterogeneity by \nControlling for Measured Confounders',
abline(a=0,b=1, lwd=2)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,100, by = 20), las = 1)
sel <- twf$err=='XbErr' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==40
with(twf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[1], lty = 1, lwd = 1))
## Err, het=0:3, obs = 0,5 Wawer
  for(bb in 1:length(hobs)) { ## options to show multiple obs lines, do this on a separate figure though
    sel <- twf$err=='XbErr' & twf$het.sd==3 & twf$hobs==hobs[bb] & twf$ehm.late==40 
    with(twf[sel,], lines(1:100, predict(loess(med~true), 1:100), col = cols[4], lty = bb, lwd = 1))
  }
legend('bottomright', leg = paste0(c(0,25,50,80),'%'), col = 'orange', lty = 1:4, bty = 'n', cex = cex.leg, ncol=1,
       title = 'variance controlled')
graphics.off()

