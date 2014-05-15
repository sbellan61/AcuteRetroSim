library(plyr)
####################################################################################################
## Viral Load Infectivity Plots, & expected EHM[acute] based on viral load trajectory
## Figures 3 and S5
####################################################################################################
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','VL Profile')
if(!file.exists(outdir)) dir.create(outdir)
## Donnell et al. (2010) , Attia et al. (2009), & Lingappa et al. (2010) data viral load-infectivity data
raw <- read.csv('data files/viral load infectivity results.csv')
## Fiebig et al. (2003) acute phase data (digitized from plot)
fbx <- read.csv('data files/Fiebig boxplot.csv')
names(fbx)[2] <- 'lvl'
## Robb (2012) acute phase data (digitized from plot)
rbaf <- read.csv('data files/Robb Africa summary.csv')

## Clean VL-infectivity data frame
names(raw)
dat <- raw[, 1:7]
## log-VL
dat$lvlow <- log(dat$Vlow,10)
dat$lvhi <- log(dat$Vhigh,10)
dat$lvlow[dat$lvlow==-Inf] <- 0
dat$lv <- (dat$lvlow+dat$lvhi)/2
## Studies by color
dat$cols <- c('purple','red','blue','orange')[as.numeric(dat$Study)]
dat$cols[1:3] <- c('purple','red','blue','orange')[4]
## Studies by Strata
dat$Strata <- as.character(dat$Strata)
dat$Strata[is.na(dat$Strata)] <- ''
dat$cat <- with(dat, paste0(Study, Strata))
dat$ltran <- with(dat, log(Transmission))
 
####################################################################################################
## Figure 3 with Robb 2012 data
pdf(file.path(outdir, 'Figure 3 - RH transmission profile Robb Fiebig.pdf'), w = 6.83, h = 3.5)
layout(matrix(c(1,1,2,3),2,2), w = c(.8,1))
par(mar=c(4,4,1,.5), 'ps'=8)#, mfrow = c(1,2))
with(dat[grepl('Ling',dat$Study),], {
  plot(1,1, type='n', xlim = c(2,7), ylim = c(.1,23), xlab = 'HIV RNA (copies/ml)', ylab = 'hazard (per = 100 person-years)',
       axes=F, bty='n', log='y', main = '(A)')#, main = cat[1])
  arrows(lv, LCI, lv, UCI, code=3, len = .15, angle = 90)
  segments(lvlow, Transmission, lvhi, Transmission, lwd = 2)
  mod <- lm(ltran ~ lv)
  xs <- seq(2,7, l = 100)
  ys <- predict(mod, newdata=data.frame(lv = xs, ltran = NA))
  lines(xs, exp(ys), col='blue', lty = 1, lwd = 1)
  vl.ticks <- c(expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7))
  axis(1, at = 2:7, lab = vl.ticks)
  axis(2, at = 10^c(-1:1), las = 2)
  axis(2, at = c(seq(.1, .9, by = .1), seq(1, 10, by = 1), 20), lab = NA)
  legfr <- unique(cbind(cols,cat))
})
par(mar=c(3,4,1,.5), mgp = c(3,1,0))
plot(rbaf$days, rbaf$logvl, ylab='', xlab = '', type = 'l', bty = 'n', ylim = c(2,7), main = '(B)', las = 1, yaxt = 'n')
vl.ticks <- c(expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7))
axis(2, at=2:7, lab = vl.ticks, las = 2)
title(xlab='days since first RNA positive', line = 2)
title(ylab='HIV RNA (copies/ml)', line = 2)
mod <- with(dat[grepl('Ling',dat$Study),],lm(ltran ~ lv)) # Lingappa model
ys.t <- exp(predict(mod, newdata=data.frame(lv = rbaf$logvl, Transmission = NA)))
ys.t <- ys.t/ys.t[length(ys.t)] ## relative to chronic
rbaf$rh <- ys.t
plot(rbaf$days, ys.t, ylab='', xlab = '', type = 'n', bty = 'n', ylim = c(0,10), main = '(C)', las = 1)
title(xlab='days since first RNA positive', line = 2)
title(ylab='relative hazard \n(versus chronic phase)', line = 2)
sel <- nrow(rbaf)
polygon(c(rbaf$days[-sel],rev(rbaf$days[-sel])), c(ys.t[-sel], rep(0, length(ys.t)-1)), col = 'purple', border = NA)
sel <- (nrow(rbaf)-1):nrow(rbaf)
polygon(c(rbaf$days[sel],rev(rbaf$days[sel])), c(rep(ys.t[nrow(rbaf)],2),0,0), col = 'orange', border = NA)
polygon(c(rbaf$days[1], rbaf$days[sel][1], rbaf$days[sel][1], rbaf$days[1]),c(rep(ys.t[nrow(rbaf)],2),0,0), col = 'orange', density = 70, angle = 45, border=NA)
## Trapezoidal rule for getting area under acute phase curve
durs <- rbaf$days[-1]-rbaf$days[-nrow(rbaf)]
hms <- .5*(rbaf$rh[-sel]+rbaf$rh[-c(1,nrow(rbaf))])*durs[-length(durs)]/30 ## hazard months
hms <- sum(hms)
hms.expected <- 1*sum(durs[-length(durs)])/30
ehms <- hms - hms.expected ## excess hazard-months
ehms
text(45, 3.5, paste0('excess hazard months = ',signif(ehms,3)), pos = 4)
graphics.off()

####################################################################################################
## Figure 3 with Fiebig 2003 data
pdf(file.path(outdir, 'Figure 3 (version 2) - RH transmission profile Lingappa Fiebig.pdf'), w = 6.83, h = 4)
layout(matrix(c(1,1,2,3),2,2))
par(mar=c(4,4,1,1), 'ps'=8)#, mfrow = c(1,2))
with(dat[grepl('Ling',dat$Study),], {
  plot(1,1, type='n', xlim = c(2,8.2), ylim = c(.005,23), xlab = 'viral load', ylab = 'hazard (per = 100 person-years)',
       axes=F, bty='n', log='y', main = '(A)')#, main = cat[1])
  arrows(lv, LCI, lv, UCI, code=3, len = .15, angle = 90)
  segments(lvlow, Transmission, lvhi, Transmission, lwd = 2)
  mod <- lm(ltran ~ lv)
  xs <- seq(2,7, l = 100)
  ys <- predict(mod, newdata=data.frame(lv = xs, ltran = NA))
  lines(xs, exp(ys), col='blue', lty = 1, lwd = 1)
  xs <- sapply(paste0('10^',2:7), function(x) as.expression((x)))
  axis(1, at = 2:7, lab = xs)
  axis(2, at = 10^c(-2:1), las = 2)
  axis(2, at = c(seq(.01,.09, by = .01), seq(.1, .9, by = .1), seq(1, 10, by = 1), 20), lab = NA)
  legfr <- unique(cbind(cols,cat))
  ## Fiebig
  polygon(c(1.8, 8.3, 8.3, 1.8), c(.005, .005, .075, .075), col = gray(.7))
  ats <- exp(rev(seq(log(.009), log(.06), l=6)))
  cols <- c('orange', rep('red',3),'purple', gray(.25))
  boxplot(lvl ~ stage, fbx, add=T, at = ats, horizontal = T, axes = F, boxwex=.23, lty = 1, col = cols)
  stages <- c(paste0(c('I: 5 ', 'II: 5.3 ', 'III: 3.2 ', 'IV: 5.6 ', 'V: 69.5 '),'days'), 'VI: chronic')
  for(ii in 1:length(ats))  text(6.5, ats[ii] , stages[ii], pos = 4, ps = 8)
})
lvl <- c(0, ddply(fbx, .(stage), with, median(lvl))[,2])
lvl <- c(lvl, lvl[length(lvl)])
stepfbx <- data.frame(day = c(0, 11, 11+c(5.0, 10.3,  13.5,  19.1,  88.6, 150)), lvl) ## eclipse is 11 days
## Function to make polygons for each piece
cols <- c('orange', rep('red',3),'purple', gray(.25))
polyf <- function(ii, col, stepfbx) with(stepfbx, polygon(c(day[ii],day[ii+1],day[ii+1],day[ii]), c(0,0, lvl[ii], lvl[ii]), col = col, border = NA))
plot(0,0, type ='n', xlim = c(0, 170), ylim = c(0, 6), xlab = 'days from infection', ylab = 'HIV RNA (cp/ml)', bty = 'n', main = '(B)')
for(ii in 2:7) polyf(ii, col = c('black',cols)[ii], stepfbx)
## Transform to relative hazard scale
mod <- with(dat[grepl('Ling',dat$Study),],lm(ltran ~ lv)) # Lingappa model
ys.t <- exp(predict(mod, newdata=data.frame(lv = stepfbx[,2], Transmission = NA)))
ys.t <- ys.t/ys.t[7] # relative to chronic
ys.t[1] <- 0 # set to 0
plot(0,0, type ='n', xlim = c(0, 170), ylim = c(0, 3), xlab = 'days from infection', ylab = 'relative hazard', bty = 'n', yaxt='n', main = '(C)')
axis(2, 0:3)
ehmfr <- data.frame(iday = c(11,5.0, 5.3,  3.2,  5.6, 69.5), iday.low = c(11, 3.1, 3.7, 2.1, 3.8, 39.7), iday.hi = c(11, 8.1, 13.5, 17, 22.9, 129.8),rhs = ys.t[1:6])
ehmfr$hm <- ehmfr$iday/30 * ehmfr$rhs
ehmfr$ehm <- ehmfr$iday/30 * (ehmfr$rhs-1)
ehmfr$ehm.low <- ehmfr$iday.low/30 * (ehmfr$rhs-1)
ehmfr$ehm.hi <- ehmfr$iday.hi/30 * (ehmfr$rhs-1)
signif(colSums(ehmfr[,c('ehm','ehm.low','ehm.hi')]),3)
for(ii in 2:7) polyf(ii, col = c('black',cols)[ii], data.frame(day = stepfbx[,1], lvl =ys.t))
polygon(c(0,99.6,99.6,0),c(0,0,1,1), lty = 2)#, col = cols[6])
text(45, 2, paste0('excess hazard days = ',signif(sum(ehmfr$ehm)*30,3)), pos = 4)
text(45, 1.5, paste0('excess hazard months = ',signif(sum(ehmfr$ehm),3)), pos = 4)
graphics.off()

####################################################################################################
## Figure S5
pdf(file.path(outdir, 'Figure S5 - published VL-infectivity relationships.pdf'), w = 3.23, h = 3)
par(mar=c(4,4,0,1), 'ps'=8, mgp = c(2.5,1,0))
plot(1,1, type='n', xlim = c(1,8), ylim = c(.01,23), xlab = 'viral load', ylab = 'hazard (per 100 person-years)',
     axes=F, bty='n', log='y')
with(dat[grepl('Ling',dat$Study),],  {
                                        #browser()
  mod <- lm(ltran ~ lv)
  xs <- seq(1,8, l = 100)
  ys <- predict(mod, newdata=data.frame(lv = xs, ltran = NA))
  lines(xs, exp(ys), col = dat$cols[grepl("Ling", dat$cat)][1])
})
vl.ticks <- c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7),expression(10^8))
axis(1, at = 1:8, lab = vl.ticks)
axis(2, at = 10^c(-2:1), c(.01,.01,'1','10'),las = 2)
axis(2, at = c(seq(.01,.09, by = .01), seq(.1, .9, by = .1), seq(1, 10, by = 1), 20), lab = NA)
ddply(dat, .(Study, Strata), with, points(lv, Transmission, type = 'p', pch = 19, col = cols))
##ddply(dat, .(Study, Strata), with, arrows(lv, LCI, lv, UCI, code=3, len = .35, angle = 90, col = cols))
##ddply(dat, .(Study, Strata), with, segments(lvlow, Transmission, lvhi, Transmission, , col = cols, lwd = 2))
legfr <- unique(dat[,c('cols','cat')])
leg <- legfr$cat
leg[1:2] <- c("Donnell 2010: CD4 200-349/ml", "Donnell 2010: CD4>350/ml")
legend('bottomright', leg = leg, col = legfr$cols, pch = 19, cex = .9, bty ='o')
graphics.off()
