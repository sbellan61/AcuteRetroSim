library(plyr)
####################################################################################################
## Viral Load Infectivity Plots
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','VL Profile')
if(!file.exists(outdir)) dir.create(outdir)
raw <- read.csv('data files/viral load infectivity results.csv')
fbx <- read.csv('data files/Fiebig boxplot.csv')
names(fbx)[2] <- 'lvl'

####################################################################################################
## Fiebig stages
fbx$stage <- round(fbx$stage)
pdf(file.path(outdir,'fiebig.pdf'))
boxplot(lvl~ stage, data = fbx, ylim = c(0,8), ylab = 'log HIV RNA (cp/ml)', xlab = 'Fiebig stage')
graphics.off()

names(raw)

dat <- raw[, 1:7]
dat$lvlow <- log(dat$Vlow,10)
dat$lvhi <- log(dat$Vhigh,10)
dat$lvlow[dat$lvlow==-Inf] <- 0
dat$lv <- (dat$lvlow+dat$lvhi)/2
dat$cols <- rainbow(4)[as.numeric(dat$Study)]
dat$cols[1:3] <- rainbow(4)[4]
dat$Strata <- as.character(dat$Strata)
dat$Strata[is.na(dat$Strata)] <- ''
dat$cat <- with(dat, paste0(Study, Strata))
dat$ltran <- with(dat, log(Transmission))



stepfun.polygon<-function(x,y,miny=min(y),col=NULL,...) {
  plot(x,y,type="n",...)
  polyx<-c(x[1],rep(x[2:length(x)],each=2),x[1])
  polyy<-c(rep(y,each=2))
  polygon(c(polyx, rev(polyx)),c(polyy, rep(0, length(polyy))),col=col)
}

mtrans<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

pdf(file.path(outdir, 'transmission profile Lingappa Fiebig.pdf'), w = 6.83, h = 4)
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
#cols <- mtrans(cols, alpha = 180)
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

     
pdf(file.path(outdir, 'transmission profile.pdf'), w = 3.23, h = 3)
par(mar=c(4,4,1,1), 'ps'=8)
plot(1,1, type='n', xlim = c(1,8), ylim = c(.01,23), xlab = 'viral load', ylab = 'hazard (per = 100 person-years)',
     axes=F, bty='n', log='y')
with(dat[grepl('Ling',dat$Study),],  {
                                        #browser()
  mod <- lm(ltran ~ lv)
  xs <- seq(1,8, l = 100)
  ys <- predict(mod, newdata=data.frame(lv = xs, ltran = NA))
  lines(xs, exp(ys), col = dat$cols[grepl("Ling", dat$cat)][1])
})
xs <- sapply(paste0('10^',1:8), function(x) as.expression((x)))
axis(1, at = 1:8, lab = xs)
axis(2, at = 10^c(-2:1), las = 2)
axis(2, at = c(seq(.01,.09, by = .01), seq(.1, .9, by = .1), seq(1, 10, by = 1), 20), lab = NA)
ddply(dat, .(Study, Strata), with, points(lv, Transmission, type = 'p', pch = 19, col = cols))
##ddply(dat, .(Study, Strata), with, arrows(lv, LCI, lv, UCI, code=3, len = .35, angle = 90, col = cols))
##ddply(dat, .(Study, Strata), with, segments(lvlow, Transmission, lvhi, Transmission, , col = cols, lwd = 2))
legfr <- unique(dat[,c('cols','cat')])
legend('bottomright', leg = legfr$cat, col = legfr$cols, pch = 19, cex = .7, bty ='n')
graphics.off()
