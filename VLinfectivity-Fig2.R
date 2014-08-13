library(plyr)
####################################################################################################
## Viral Load Infectivity Plots, & expected EHM[acute] based on viral load trajectory
## Figures 3 and S5
####################################################################################################
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
outdir <- file.path('FiguresAndTables','VL Profile')
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
pdf(file.path(outdir, 'Figure 3 - RH transmission profile Robb Lingappa.pdf'), w = 6.83, h = 3.5)
layout(matrix(c(1,1,2,3),2,2), w = c(.8,1))
par(mar=c(4,4,1,.5), 'ps'=10)#, mfrow = c(1,2))
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
## Redo Lingappa's model for plotting purposes
mod <- with(dat[grepl('Ling',dat$Study),],lm(ltran ~ lv)) # Lingappa model
ys.t <- exp(predict(mod, newdata=data.frame(lv = rbaf$logvl, Transmission = NA)))
ys.t <- ys.t/ys.t[length(ys.t)] ## relative to chronic
rbaf$rh <- ys.t
## Use Lingappa's actual results to get the EHM calculation (since it accounts for different #'s of individual in each VL strata)
## Get CI's on EHM
vl.inf50 <- c(.97,.74,.6) ## Lingappa's necessary reduction of logVL to achieve a 50% reduction in infectivity
vl.inf <- exp(log(2)/vl.inf50) ## to get increase in infectivity per log10 increase in VL
## Now use these to get CI on EHM.acute
rbaf$logvlabove.ch <- rbaf$logvl-rbaf$logvl[nrow(rbaf)]
rbaf$rh.lci <- vl.inf[1]^rbaf$logvlabove.ch
rbaf$rh.med <- vl.inf[2]^rbaf$logvlabove.ch
rbaf$rh.uci <- vl.inf[3]^rbaf$logvlabove.ch
plot(rbaf$days, rbaf$rh.med, ylab='', xlab = '', type = 'n', bty = 'n', ylim = c(0,10), main = '(C)', las = 1)
title(xlab='days since first RNA positive', line = 2)
title(ylab='relative hazard \n(versus chronic phase)', line = 2)
sel <- nrow(rbaf)
polygon(c(rbaf$days[-sel],rev(rbaf$days[-sel])), c(ys.t[-sel], rep(0, length(ys.t)-1)), col = 'purple', border = NA)
sel <- (nrow(rbaf)-1):nrow(rbaf)
polygon(c(rbaf$days[sel],rev(rbaf$days[sel])), c(rep(ys.t[nrow(rbaf)],2),0,0), col = 'orange', border = NA)
polygon(c(rbaf$days[1], rbaf$days[sel][1], rbaf$days[sel][1], rbaf$days[1]),c(rep(ys.t[nrow(rbaf)],2),0,0), col = 'orange', angle = 45, border=NA)
## Trapezoidal rule for getting area under acute phase curve
ehms <- numeric(4)
names(ehms) <- c('digitized','lci','med','uci')
for(ii in 1:4) { ## median and CI on ehm acute
    rh.temp <- rbaf[,c('rh','rh.lci','rh.med','rh.uci')[ii]]
    durs <- rbaf$days[-1]-rbaf$days[-nrow(rbaf)]
    hms <- .5*(rh.temp[-sel]+rh.temp[-c(1,nrow(rbaf))])*durs[-length(durs)]/30 ## hazard months
    hms <- sum(hms)
    hms.expected <- 1*sum(durs[-length(durs)])/30
    ehms[ii] <- hms - hms.expected ## excess hazard-months
}
ehms.vl <- ehms
ehms.vl
text(45, 3.5, paste0('excess hazard months = ',signif(ehms['med'],3),' (',signif(ehms['lci'],3),'-', signif(ehms['uci'],3),')'), pos = 4)
graphics.off()

save(ehms.vl, file=file.path(outdir, 'ehms.vl.Rdata'))

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
