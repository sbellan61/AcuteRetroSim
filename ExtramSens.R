library(mgcv);library(plyr); library(data.table); library(abind); library(multicore)
rm(list=ls(all=T)); gc()
####################################################################################################
## Summarize Hollingsworth & Wawer style fits to simulated data.
## Create wf.Rdata with simulation results stored for plotting.
####################################################################################################
outdir <- file.path('FiguresAndTables','UgandaFitSummaries')
if(!file.exists(outdir)) dir.create(outdir)
load(file.path('results','RakAcute','blocks.Rdata')) ## these are country-acute phase specific blocks
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') ## transmission coefficient names, for convenience
nc <- 12                                       ## core per simulation

## Load fit files (output from RakFitMK.R files)
fls <- list.files(file.path('results','RakAcute','UgandaFitsExtram','fitouts'), pattern = 'fitout-', full.names=T)
length(fls)
jobnums <-  as.numeric(sapply(fls, function(x) as.numeric(strsplit(x,'-')[[1]][2])))
fls <- fls[order(jobnums)]
jobnums <- jobnums[order(jobnums)]

to.do <- with(blocks, which(dur.lt==10 & late.sc==c(1,5)))
jtd <- to.do[!to.do %in% jobnums]

## Load & label files by jobnum (fitout1, fitout2, etc...)
system.time(alldat <- mclapply(fls, function(x) { jobnum <- as.numeric(strsplit(x,'-')[[1]][2]); load(x)
                                    assign(paste0('fitout',jobnum),fitout, env=.GlobalEnv)}))

## Each list element is one simulation output from a job
names(alldat[[1]])
table(unlist(lapply(alldat, function(x) x$decont))) ## decontaminating between phases?
table(unlist(lapply(alldat, function(x) x$excl.extram)))## excluding extramarital transmission?
alldat[[1]]$wtab ## raw simulation data
       
## Sample Sizes
wtabs <- t(abind(lapply(alldat, function(x) { c(x$wtab$inct[1,'n'],x$wtab$prevt[1,'n'],x$wtab$latet[nrow(x$wtab$latet),'n'])}), along = 2))
wtabs <- as.data.frame(wtabs)
colnames(wtabs) <- c('incn','prevn','laten')
wtabs$job <- jobnums
apply(wtabs[,-4],2, function(x) quantile(x, c(.025,.5,.975))) ## summary of sample sizes by couple classing

dimnames(alldat[[1]]$outtab.h)

## Hollingsworth Fits
hfitsExtram <- abind(lapply(alldat, function(x) { x$outtab.h}), along = 4)
dimnames(hfitsExtram)[[4]] <- jobnums
dimnames(hfitsExtram)[[1]] <- c('lci','med','uci','true')

## Wawer Fits (Poisson Regressions)
wfitsExtram <- abind(lapply(alldat, function(x) { x$outtab.w}), along = 6)
dimnames(wfitsExtram)[[6]] <- jobnums
dimnames(wfitsExtram)[[1]] <- c('lci','med','uci','true')

## Empirical Wawer Hazards (no model fit)
ewfitsExtram <- t(abind(lapply(alldat, function(x) {x$erhs.w}), along = 2))
rownames(ewfitsExtram) <- jobnums

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(alldat, function(x) {as.data.frame(t(x$dpars))})))
cohvars$job <- jobnums
dpars <- blocks[jobnums,]

## turn these into easy to work with data frames
hfExtram <- adply(hfitsExtram, 4:2, .parallel=T)
colnames(hfExtram)[1:3] <- c('job','err','var')
hfExtram$het.sd <- dpars$het.gen.sd[match(hfExtram$job, dpars$job)]
hfExtram <- merge(hfExtram, cohvars, by = 'job')
hfExtram <- hfExtram[order(hfExtram$err,hfExtram$var,hfExtram$job),]
hfExtram <- within(hfExtram, {ehm.late <- (late.sc - 1)*dur.lt})
hfExtram <- within(hfExtram, {ehm.acute <- (acute.sc - 1)*dur.ac})
head(hfExtram)

## this can take a long time because each simulation had 24model fits
obs.sh <- paste0('obs',c(0,0.3,.5,.7,1))
wfitsExtram.small <- wfitsExtram[,c('ehm.ac','acute.sc','dur.ac','ehm.lt','late.sc','dur.lt'),,,,]
wfExtram <- adply(wfitsExtram.small, 6:2, .parallel=T)
colnames(wfExtram)[1:5] <- c('job','cov','hobs','err','var')
wfExtram$het.sd <- dpars$het.gen.sd[match(wfExtram$job, dpars$job)]
wfExtram <- merge(wfExtram, cohvars, by = 'job')
wfExtram <- wfExtram[order(wfExtram$err,wfExtram$var,wfExtram$job),]
wfExtram <- within(wfExtram, {ehm.late <- (late.sc - 1)*dur.lt})
wfExtram <- within(wfExtram, {ehm.acute <- (acute.sc - 1)*dur.ac})
save(wfExtram, hfExtram, file=file.path(outdir, 'wfExtram.Rdata'))

load(file=file.path(outdir, 'wfExtram.Rdata'))
load(file=file.path(outdir, 'wf.Rdata'))
head(wfExtram)

names(hfExtram)
xtabs(~acute.sc + dur.ac + het.sd, hfExtram[hfExtram$var=='ehm.ac',])
td <- which(with(blocks, acute.sc==5 & dur.ac == 3 & late.sc == 1 & het.gen.sd==1 & dur.lt==10))
blocks[td,]

td %in% to.do


####################################################################################################
## Figure SX Extra-couple effect
##############################
pdf(file.path(outdir, 'Figure SX - extracouple effect.pdf'), w = 6.83, h = 4)
asp <- 1
hsds <- 1
cols <- c('black','red','blue','orange')
bias2col <- gray(.4)
bias1col <- gray(.6)
yx.col <- gray(.8)
ct <- 12
par(mfrow = c(1,2), mar = c(4,4.5,2.5,1.5), mgp=c(2.5,1,0), 'ps'=ct, oma = c(0,0,0,0)) #, cex.main = ct, cex.lab = ct, cex.axis = ct)
cex.leg <- 1
lwd.abline <- 5
lwd <- 2
madj <- .5
mln <- 2
mlnv <- 3
##################################################
## ehm plots
xmax <- 80
xmax.loess <- 80
ymax <- 100
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var,] 
twf <- wf[wf$var==var &  wf$cov=='',]
thfExtram <- hfExtram[hfExtram$var==var,] 
twfExtram <- wfExtram[wfExtram$var==var &  wfExtram$cov=='',]
## Plot layout
ylab <- expression(paste('estimated ',EHM[acute]))
xlab <- expression(paste('true ',EHM[acute]))
## Wawer plot
par('ps'=ct)

plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = ylab, main = '(A) Wawer et al. Model', axes = F, asp=asp)
clip(0, xmax, 0, ymax)
abline(a=0,b=1, lwd=lwd.abline, col=yx.col)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,xmax, by = 20), las = 1)
## seroincident ltf couples included EHM_late=00
sel <- twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==0
with(twf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 1, lwd = lwd))
sel <- twfExtram$err=='base' & twfExtram$het.sd==0 & twfExtram$hobs=='obs0' & twfExtram$ehm.late==0
##with(twfExtram[sel,], points(true, med, col = bias1col, pch = 16))
with(twfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 2, lwd = lwd))
## seroincident ltf couples included EHM_late=40
sel <- twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==40
with(twf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias2col, lty = 1, lwd = lwd))
sel <- twfExtram$err=='base' & twfExtram$het.sd==0 & twfExtram$hobs=='obs0' & twfExtram$ehm.late==40
##with(twfExtram[sel,], points(true, med, col = bias2col, pch =16))
with(twfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias2col, lty = 2, lwd = lwd))
## Err, het=0:3, obs = 0,5 Wawer
hsds <- 0:3
for(hh in 1:length(hsds)) {
  sel <- twf$err=='XbErr' & twf$het.sd==hsds[hh] & twf$hobs=='obs0' & twf$ehm.late==40
  with(twf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = cols[hh], lty = 1, lwd = lwd))
}
## Err, het=0:3, obs = 0,5 Wawer
for(hh in 1:length(hsds)) {
  sel <- twfExtram$err=='XbErr' & twfExtram$het.sd==hsds[hh] & twfExtram$hobs=='obs0' & twfExtram$ehm.late==40
##  with(twfExtram[sel,], points(true, med, col = cols[hh], pch = 16))
  with(twfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = cols[hh], lty = 2, lwd = lwd))
}

##################################################
## Holl plot
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = ylab, main = '(B) Hollingsworth et al. Model', axes = F, asp=asp)
clip(0, xmax, 0, ymax)
abline(a=0,b=1, lwd=lwd.abline, col=yx.col)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,xmax, by = 20), las = 1)
## seroincident ltf couples included EHM_late=00
sel <- thf$err=='base' & thf$het.sd==0 &  thf$ehm.late==0
with(thf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 1, lwd = lwd))
sel <- thfExtram$err=='base' & thfExtram$het.sd==0 & thfExtram$ehm.late==0
with(thfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 2, lwd = lwd))
## seroincident ltf couples included EHM_late=40
sel <- thf$err=='base' & thf$het.sd==0 & thf$ehm.late==40
with(thf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias2col, lty = 1, lwd = lwd))
sel <- thfExtram$err=='base' & thfExtram$het.sd==0 & thfExtram$ehm.late==40
with(thfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias2col, lty = 2, lwd = lwd))
## Err, het=0:3, obs = 0,5 Wawer
hsds <- 0:3
for(hh in 1:length(hsds)) {
  sel <- thf$err=='XbErr' & thf$het.sd==hsds[hh] & thf$ehm.late==40
  with(thf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = cols[hh], lty = 1, lwd = lwd))
}
## Err, het=0:3, obs = 0,5 Wawer
for(hh in 1:length(hsds)) {
  sel <- thfExtram$err=='XbErr' & thfExtram$het.sd==hsds[hh] & thfExtram$ehm.late==40
  with(thfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = cols[hh], lty = 2, lwd = lwd))
}
mtext(xlab, outer =T, side = 1, line = -1)
graphics.off()

####################################################################################################

####################################################################################################
## Figure SX Extra-couple effect
##############################
pdf(file.path(outdir, 'Figure SX - extracouple effect (last bias?).pdf'), w = 6.83, h = 4)
asp <- 1
hsds <- 1
cols <- c('black','red','blue','orange')
bias2col <- gray(.4)
bias1col <- gray(.6)
yx.col <- gray(.8)
ct <- 12
par(mfrow = c(1,2), mar = c(4,4.5,2.5,1.5), mgp=c(2.5,1,0), 'ps'=ct, oma = c(0,0,0,0)) #, cex.main = ct, cex.lab = ct, cex.axis = ct)
cex.leg <- 1
lwd.abline <- 5
lwd <- 2
madj <- .5
mln <- 2
mlnv <- 3
##################################################
## ehm plots
xmax <- 80
xmax.loess <- 80
ymax <- 100
var <- 'ehm.ac'
## just get subsets we're interested in
thf <- hf[hf$var==var,] 
twf <- wf[wf$var==var &  wf$cov=='',]
thfExtram <- hfExtram[hfExtram$var==var,] 
twfExtram <- wfExtram[wfExtram$var==var &  wfExtram$cov=='',]
## Plot layout
ylab <- expression(paste('estimated ',EHM[acute]))
xlab <- expression(paste('true ',EHM[acute]))
## Wawer plot
par('ps'=ct)
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = ylab, main = '(A) Wawer et al. Model', axes = F, asp=asp)
clip(0, xmax, 0, ymax)
abline(a=0,b=1, lwd=lwd.abline, col=yx.col)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,xmax, by = 20), las = 1)
## seroincident ltf couples included EHM_late=00
sel <- twf$err=='base' & twf$het.sd==0 & twf$hobs=='obs0' & twf$ehm.late==0
with(twf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 1, lwd = lwd))
sel <- twfExtram$err=='base' & twfExtram$het.sd==0 & twfExtram$hobs=='obs0' & twfExtram$ehm.late==0
##with(twfExtram[sel,], points(true, med, col = bias1col, pch = 16))
with(twfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 2, lwd = lwd))
##################################################
## Holl plot
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = ylab, main = '(B) Hollingsworth et al. Model', axes = F, asp=asp)
clip(0, xmax, 0, ymax)
abline(a=0,b=1, lwd=lwd.abline, col=yx.col)
axis(2, seq(0,100, by = 20), las = 2)
axis(1, seq(0,xmax, by = 20), las = 1)
## seroincident ltf couples included EHM_late=00
sel <- thf$err=='base' & thf$het.sd==0 &  thf$ehm.late==0
with(thf[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 1, lwd = lwd))
sel <- thfExtram$err=='base' & thfExtram$het.sd==0 & thfExtram$ehm.late==0
with(thfExtram[sel,], lines(1:xmax.loess, predict(loess(med~true), 1:xmax.loess), col = bias1col, lty = 2, lwd = lwd))
mtext(xlab, outer =T, side = 1, line = -1)
graphics.off()
