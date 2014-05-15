library(mgcv);library(plyr); library(data.table); library(abind); library(multicore)
rm(list=ls(all=T)); gc()
####################################################################################################
## Summarize Hollingsworth & Wawer style fits to simulated data.
## Create wf.Rdata with simulation results stored for plotting.
####################################################################################################
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     ## setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')
if(!file.exists(outdir)) dir.create(outdir)
load(file.path('results','RakAcute','blocks.Rdata')) ## these are country-acute phase specific blocks
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') ## transmission coefficient names, for convenience
nc <- 12                                       ## core per simulation

## Load fit files (output from RakFitMK.R files)
fls <- list.files(file.path('results','RakAcute','UgandaFits','fitouts'), pattern = 'fitout-', full.names=T)
length(fls)
jobnums <-  as.numeric(sapply(fls, function(x) as.numeric(strsplit(x,'-')[[1]][2])))
fls <- fls[order(jobnums)]
jobnums <- jobnums[order(jobnums)]
load(file=file.path(outdir, 'jtd.Rdata'))

fls <- fls[!jobnums %in% jobnums.to.do]
jobnums <- jobnums[!jobnums %in% jobnums.to.do]

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

## Hollingsworth Fits
hfits <- abind(lapply(alldat, function(x) { x$outtab.h}), along = 4)
dimnames(hfits)[[4]] <- jobnums
dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## Wawer Fits (Poisson Regressions)
wfits <- abind(lapply(alldat, function(x) { x$outtab.w}), along = 6)
dimnames(wfits)[[6]] <- jobnums
dimnames(wfits)[[1]] <- c('lci','med','uci','true')

## Empirical Wawer Hazards (no model fit)
ewfits <- t(abind(lapply(alldat, function(x) {x$erhs.w}), along = 2))
rownames(ewfits) <- jobnums

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(alldat, function(x) {as.data.frame(t(x$dpars))})))
cohvars$job <- jobnums
dpars <- blocks[jobnums,]

## turn these into easy to work with data frames
hf <- adply(hfits, 4:2, .parallel=T)
colnames(hf)[1:3] <- c('job','err','var')
hf$het.sd <- dpars$het.gen.sd[match(hf$job, dpars$job)]
hf <- merge(hf, cohvars, by = 'job')
hf <- hf[order(hf$err,hf$var,hf$job),]
hf <- within(hf, {ehm.late <- (late.sc - 1)*dur.lt})
hf <- within(hf, {ehm.acute <- (acute.sc - 1)*dur.ac})
head(hf)

## this can take a long time because each simulation had 24model fits
obs.sh <- paste0('obs',c(0,0.3,.5,.7,1))
wfits.small <- wfits[,c('ehm.ac','acute.sc','dur.ac','ehm.lt','late.sc','dur.lt'),,,,]
wf <- adply(wfits.small, 6:2, .parallel=T)
colnames(wf)[1:5] <- c('job','cov','hobs','err','var')
wf$het.sd <- dpars$het.gen.sd[match(wf$job, dpars$job)]
wf <- merge(wf, cohvars, by = 'job')
wf <- wf[order(wf$err,wf$var,wf$job),]
wf <- within(wf, {ehm.late <- (late.sc - 1)*dur.lt})
wf <- within(wf, {ehm.acute <- (acute.sc - 1)*dur.ac})
save(wf, hf, file=file.path(outdir, 'wf.Rdata'))

load(file=file.path(outdir, 'wf.Rdata'))
head(wf)

####################################################################################################
## Hollingsworth Model: True vs Estimated by EHM.late
####################################################################################################
ct <- .7
xmax <- 100
ymax <- 130
cols <- rep(NA,100)
## cols[unique(hf$late.sc)] <- c('blue','orange','red')
ehm.lates <- unique(c(hf$ehm.late))
ehm.lates <- ehm.lates[order(ehm.lates)][-1]
cols[ehm.lates+1] <- colorRampPalette(c('purple','orange','red'))(length(ehm.lates))
cols <- paste0(cols,'8f')
pdf(file.path(outdir,'Hollingsworth Analysis of Couples Cohort Simulated Data.pdf'), w = 6.83, h = 3)
for(vv in (unique(hf$het.sd))) {
  var <- 'ehm.ac'
  par(mfrow=c(1,2), mar = c(3.2,3,2,.5), cex.main = ct, cex.lab = ct*1.2, cex.axis = ct, oma = c(0,1,1,0))
  for(ee in 1:2) {
    err <- c('base','XbErr')[ee]
    err.lab <- paste0('(',LETTERS[1:2], ') ', c('full inclusion','excluded incident SDCs lost to follow-up'))
    plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '',#expression(paste('true ',EHM[acute])),
         ylab = '', main = err.lab[ee], las = 1)
    sel <- hf$var==var & hf$err==err & hf$het.sd==vv & hf$late.sc!=1# & hf$dur.ac==2 #
    ddply(hf[sel,], .(het.sd), with, {arrows(true+late.sc/10, lci, true+late.sc/10,uci, length=.05, angle = 90, code = 3,col = cols[ehm.late+1])})    
    ddply(hf[sel,], .(het.sd), with, {points(true+late.sc/10, med, col = cols[ehm.late+1], cex = .7, pch = 19)}) #, col = cols[het.sd+1])})
    abline(a=0, b=1, col = 'black', lty = 1)
    segments(0, 75.4, xmax, 75.4, col = 'black', lty = 3, lwd = 2)
    ## text(xmax*.25, 72.4, 'Hollingsworth \nEstimate', pos = 3, cex = .5)
    ## axis(2, 75.4, 75.4, las = 2, line = -2, tick=F)
    ## segments(0, 36.25, 100, 36.25, col = 'red', lty = 3, lwd = 2)
    ## text(100, 36.25, 'Wawer \nEstimate', pos = 4, cex = .5)    
    if(ee==1) legend('topleft', title= expression(EHM['late']), ##'late to chronic relative hazard \n(for months 20-11 prior to death)',
         leg = ehm.lates, pch = 19, col=cols[ehm.lates+1], cex = .7,  bty ='n')
  }
  mtext(paste('indivual heterogeneity sd=',vv), side = 3, outer = T, line = -.5, adj = .5, cex = ct)
  mtext(expression(paste('estimated ',EHM[acute])), side = 2, line = 0, adj = .5, outer = T, cex = ct*1.2)
  mtext(expression(paste('true ',EHM[acute])), side = 1, line = -1, adj = .5, outer = T, cex = ct*1.2)
}
dev.off()


####################################################################################################
## Hollingsworth Model: True vs Estimated: by acute phase duration
####################################################################################################
cols <- colorRampPalette(c('purple','orange'))(10)
lt.sc.temp <- 5 ## only plot for one late scalar 
pdf(file.path(outdir,'Hollingsworth by acute duration.pdf'), w = 8, h = 5)
for(vv in (unique(hf$het.sd))) {
  var <- 'ehm.ac'
  par(mfrow=c(1,2), oma = c(0,0,1,0))
  for(ee in 1:2) {
    err <- c('base','XbErr')[ee]
    err.lab <- c('correct inclusion','excluded couples \nlost to follow-up by error')
    plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = expression(paste('true ',EHM[acute])),
         ylab = expression(paste('estimated ',EHM[acute])), main = err.lab[ee])
    sel <- hf$var==var & hf$err==err & hf$het.sd==vv & hf$late.sc==lt.sc.temp ## & hf$dur.ac==2 #
    ddply(hf[sel,], .(het.sd), with, {arrows(true+late.sc/10, lci, true+late.sc/10,uci, length=.05, angle = 90, code = 3,col = cols[dur.ac*2])})    
    ddply(hf[sel,], .(het.sd), with, {points(true+late.sc/10, med, col = cols[dur.ac*2], pch = 19, cex = .7)}) #, col = cols[het.sd+1])})
    abline(a=0, b=1, col = 'black', lty = 3)
    segments(0, 75.4, 100, 75.4, col = 'red', lty = 2, lwd = 2)
    text(100, 75.4, 'Hollingsworth \nEstimate', pos = 4, cex = .5)
    ## segments(0, 36.25, 100, 36.25, col = 'red', lty = 3, lwd = 2)
    ## text(100, 36.25, 'Wawer \nEstimate', pos = 4, cex = .5)    
    if(ee==1) legend(5,ymax, title= 'acute phase duration (months)', ncol = 2,
         leg = unique(hf$dur.ac), pch = 19, col=cols[unique(hf$dur.ac)*2], bty ='n', cex = .7)
  }
  mtext(paste('indivual heterogeneity sd=',vv, 'late RH=',lt.sc.temp), side = 3, outer = T, line = -.5, adj = .5, cex = 1.5)
}
dev.off()
 
####################################################################################################
## Hollingsworth Model: True vs Estimated ********EHM_LATE* (by acute.sc)
####################################################################################################
xmax <- 130
ymax <- 130
axis.cex <- 1.5
cols <- rep(NA,100)
acute.scalars <- unique(hf$acute.sc)[order(unique(hf$acute.sc))]
cols[acute.scalars] <- colorRampPalette(c('purple','orange'))(13)
pdf(file.path(outdir,'Hollingsworth EHM.LATE estimates of Couples Cohort Simulated Data.pdf'), w = 8, h = 8)
for(vv in (unique(hf$het.sd))) {
  var <- 'ehm.lt'
  par(mfrow=c(4,2), mar = c(3,3,2,.5), oma = c(4,4,4,0))
  for(ad in c(.5,1,2,4)) {
    for(ee in 1:2) {
      err <- c('base','XbErr')[ee]
      err.lab <- c('correct inclusion','excluded couples \nlost to follow-up by error')
      plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '',
           ylab = '', main = ifelse(ad==.5,err.lab[ee],''))
      sel <- hf$var==var & hf$err==err & hf$het.sd==vv & hf$dur.ac==ad #& hf$late.sc==10
      temp <- hf[sel,]
      temp$xs <- with(temp, jitter(true+acute.sc/4, a = .3))
      ddply(temp, .(het.sd), with, {arrows(xs, lci, xs,uci, length=.05, angle = 90, code = 3,col = cols[acute.sc])})    
      ddply(temp, .(het.sd), with, {points(xs, med, col = cols[acute.sc], cex = .7, pch = 19)}) #, col = cols[het.sd+1])})
      abline(a=0, b=1, col = 'black', lty = 3)
      ## segments(0, 75.4, 100, 75.4, col = 'red', lty = 2, lwd = 2)
      ## text(100, 75.4, 'Hollingsworth \nEstimate', pos = 4, cex = .5)
      ## segments(0, 36.25, 100, 36.25, col = 'red', lty = 3, lwd = 2)
      ## text(100, 36.25, 'Wawer \nEstimate', pos = 4, cex = .5)    
      if(ee==1 & ad==.5) legend('topleft', title= 'acute to chronic relative hazard',
                   leg = acute.scalars, ncol = 3, pch = 19, col=cols[acute.scalars], bty ='n', cex = .7)
      if(ee==1)  mtext(paste(ad,'month acute phase'), side = 2, line = 5)
    }
  }
  mtext(paste('indivual heterogeneity sd=',vv), side = 3, outer = T, line = 2, adj = .5, cex = 1.5)
  mtext(expression(paste('true ',EHM[late], ' (jittered')), side = 1, outer = T, line = 1.5, cex = axis.cex)
  mtext(expression(paste('estimated ',EHM[late])), side = 2, outer = T, line = -.5, cex = axis.cex)
}
dev.off()


####################################################################################################
## Wawer Model: True vs Estimated
####################################################################################################
## Plot of Wawer-type analyses of simulated data across different amounts of heterogeneity,
## different acute/late parameters, & w/ & w/o exclusion error.
xmax <- 130
ymax <- 150
axis.cex <- 1.5
leg <- unique(wf$late.sc)
cols <- rep(NA,10)
cols[unique(wf$late.sc)] <- c('blue','orange','red')
legcols <- cols[unique(wf$late.sc)]
err.lab <- c('correct inclusion','excluded couples \nlost to follow-up by error')
pdf(file.path(outdir,'Wawer Analysis of Couples Cohort Simulated Data.pdf'), w = 10, h = 6)
for(vv in unique(wf$het.sd)) {
  var.an <- 'ehm.ac'
  par(mfrow=c(2,5), mar=c(2.5,3,1,.5), oma = c(3,7.5,7,0))
  #layout(matrix(1:8,4,2))
  for(ee in 1:2) {
    for(bb in c(0,0.3,.5,.7,1)) { ## amount of correlation between measured confounders and underlying confounders
      obs <- paste0('obs',bb)
      err <- c('base','XbErr')[ee]
      plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '')
      sel <- wf$var==var.an & wf$err==err & wf$het.sd==vv & wf$hobs==obs & wf$cov==''# & wf$dur.ac==2 #& wf$late.sc==10
      temp <- wf[sel,]
      ddply(temp, .(het.sd), with, {arrows(true+late.sc/10, lci, true+late.sc/10,uci, length=.05, angle = 90, code = 3,col = cols[late.sc])})    
      ddply(temp, .(het.sd), with, {points(true+late.sc/10, med, cex = .7, col = cols[late.sc], pch = 19)}) #, col = cols[het.sd+1])})
      abline(a=0, b=1, col = 'black', lty = 1, lwd = 1)
      ## segments(0, 75.4, 100, 75.4, col = 'red', lty = 2, lwd = 2)
      ## text(100, 75.4, 'Hollingsworth \nEstimate', pos = 4, cex = .5)
      segments(0, 36.25, 100, 36.25, col = 'red', lty = 3, lwd = 2)
      ## text(100, 36.25, 'Wawer \nEstimate', pos = 4, cex = .5)
      if(ee==1) mtext(bb, side = 3, line = 1, cex = 1.2, las = 1)
      if(ee==1 & bb==0) legend(5,ymax, title= 'late to chronic relative hazard \n(for months 20-11 prior to death)',
                   leg = leg, pch = 19, col=legcols, bty ='n', cex = 1)
      if(bb==0)  mtext(paste(err.lab[ee]), side = 2, line = 7)
    }
  }
  mtext(paste('indivual heterogeneity sd=',vv), side = 3, outer = T, line = 5, adj = .5, cex = axis.cex)
  mtext(expression(paste('true ',EHM[acute])), side = 1, outer = T, line = 1.5, cex = axis.cex)
  mtext(expression(paste('estimated ',EHM[acute])), side = 2, outer = T, line = .5, cex = axis.cex)
  mtext(expression(underline('proportion of heterogeneity observed and controlled for')), side = 3, outer = T, line = 2, cex = 1.2)
}
dev.off()

####################################################################################################
## Log-hazard distribution figure
pdf(file.path(outdir, 'log-hazard distributions.pdf'), w = 3.27, h = 3)
par(mar = c(4,.5,1,.5), 'ps'=8)
plot(0,0, type = 'n', xlab = bquote(paste('multiple of ', lambda[hazard])), ylab='', main='', bty = 'nn',
     xlim = c(10^-3,10^3), log='x',ylim = c(0,1), axes = F)
axis(1, at = 10^(-3:3), lab = c(0.001, 0.01, 0.1, 1, 10, 100, 1000))
cols <- rainbow(4)
cols[1] <- 'orange'
for(hsd in 1:3) {
  curve(dnorm(log(x), 0,  hsd), from = 10^-3, to = 10^3, add = T, col = cols[hsd+1])
}
segments(1,0,1,1, col = 'orange')
legend('topleft', leg = 0:3, lty = 1, col = cols, title = expression(sigma[hazard]))
graphics.off()


