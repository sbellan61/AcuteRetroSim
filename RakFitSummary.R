library(plyr); library(data.table); library(abind); library(multicore)
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')
if(!file.exists(outdir)) dir.create(outdir)
load(file.path('results','RakAcute','blocks.Rdata')) # these are country-acute phase specific blocks
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation

## Load fit files
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

## Get rid of those with the wrong dimensions (one batch screwed up)
dims <- unlist(lapply(alldat, function(x) { dim(x$outtab.h)[3]}))
## dims
fls <- fls[!is.na(dims)]
jobnums.to.do <- jobnums[is.na(dims)]
length(jobnums.to.do)
save(jobnums.to.do, file=file.path(outdir, 'jtd.Rdata'))
jobnums <- jobnums[!is.na(dims)]
rm(alldat); gc()

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
head(wtabs)
apply(wtabs[,-4],2, summary)

## Hollingsworth Fits
hfits <- abind(lapply(alldat, function(x) { x$outtab.h}), along = 4)
dimnames(hfits)[[4]] <- jobnums
dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## Wawer Fits (Poisson Regressions)
wfits <- abind(lapply(alldat, function(x) { x$outtab.w}), along = 6)
dimnames(wfits)[[6]] <- jobnums
#dimnames(wfits)
dimnames(wfits)[[1]] <- c('lci','med','uci','true')

## Empirical Wawer Hazards (no model fit)
ewfits <- t(abind(lapply(alldat, function(x) {x$erhs.w}), along = 2))
rownames(ewfits) <- jobnums

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(alldat, function(x) {as.data.frame(t(x$dpars))})))
cohvars$job <- jobnums
## Pars
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

## ## this can take a long time because each simulation had 24model fits
## obs.sh <- paste0('obs',c(0,0.3,.5,.7,1))
## wfits.small <- wfits[,c('ehm.ac','acute.sc','dur.ac','ehm.lt','late.sc','dur.lt'),,,,]
## wf <- adply(wfits.small, 6:2, .parallel=T)
## colnames(wf)[1:5] <- c('job','cov','hobs','err','var')
## wf$het.sd <- dpars$het.gen.sd[match(wf$job, dpars$job)]
## wf <- merge(wf, cohvars, by = 'job')
## wf <- wf[order(wf$err,wf$var,wf$job),]
## wf <- within(wf, {ehm.late <- (late.sc - 1)*dur.lt})
## wf <- within(wf, {ehm.acute <- (acute.sc - 1)*dur.ac})
## save(wf, hf, file=file.path(outdir, 'wf.Rdata'))

load(file=file.path(outdir, 'wf.Rdata'))
head(wf)


####################################################################################################
## Are any of the age or marital duration covariates significant predictors of hazard were included
## in a Poisson model? THESE ARE ALREADY IN YEARS
## all duration covariates are labeled secp.pdsa, but dimensions are labeled appropriately
alldat[[1]]$outtab.w[,'secp.pdsa',,,] ## one simulation run
covests <- abind(lapply(alldat, function(x) { x$outtab.w[,'secp.pdsa',,,]}), along = 5)
dim(covests)
dimnames(covests)[1:4]
dimnames(covests)[[5]] <- jobnums
cef <- adply(covests, 5:2, .parallel=TRUE)
colnames(cef)[1:4] <- c('job','cov','var','err')
cef$het.sd <- dpars$het.gen.sd[match(cef$job, dpars$job)]
cef <- merge(cef, cohvars, by = 'job')
tail(cef)
## turn it into 5 year RH
cef[,c('lci','med','uci')] <- (cef[,c('lci','med','uci')])^5
range(cef$med, na.rm=T)


####################################################################################################
## Are the empiricaly hazards tallied during simulation equivalent to the parameters that are going
## in?
alldat[[1]]$infpm ## tallied infections & person-months at risk during event-driven model
## Look at raw hazards, tabulated during simulation
rawhaz <- as.data.frame(t(abind(lapply(alldat, function(x) { colMeans(x$infpm[,,1])/colMeans(x$infpm[,,2])}), along = 2)))
rawrhs <- rawhaz/rawhaz$ch
rawrhs$job <- jobnums
rawrhs <- merge(rawrhs, cohvars, by = 'job')
## throw out simulations with fractional acute periods since that complicates calculations (could
## figure it out algebraically later)
rawrhs <- rawrhs[rawrhs$dur.ac==round(rawrhs$dur.ac),]
rawrhs$het.sd <- dpars$het.gen.sd[match(rawrhs$job, dpars$job)]
## Plot this 
pdf(file.path(outdir,'check simulated hazard.pdf'))
cols <- rep(NA,100)
dur.acs <- unique(rawrhs$dur.ac)[order(unique(rawrhs$dur.ac))]
cols[dur.acs] <- colorRampPalette(c('purple','orange'))(length(dur.acs))
with(rawrhs[rawrhs$het.sd==0,], plot(acute.sc, ac, col = cols[dur.ac], pch = 19, cex = .5, bty = 'n', ylim = c(0,50), xlim = c(0,50),
                  xlab = expression(paste('true ',RH[acute])), ylab = expression(paste('simulated ',RH[acute]))))
abline(a = 0, b = 1)
##  for large acute phase relative hazards, simulated relative hazards (calculated as number of
##  infections divided by number of person-month exposed) are consistently smaller than the actual
##  parameters that are used to simulate them. I think this is because of the discretization at the
##  monthly level. For very high acute phase relative hazards, individuals are likely to have been
##  infected very early in a month of exposure. Thus, when we tally it up by calculating person
##  months in monthly units we overestimate the person time at risk in the denominator.
## 
##  There an analytical relationship that we can use to correct for this.
hzp <- function(x, bp = .007) { (1-exp(-bp*x)) / (1-exp(-bp)) }
curve(hzp, from = 1, to = 50, add = T)
dev.off()
## However, the reason Hollingsworth is underestimating isn't because of this since their likelihood
## model isn't discretizing at all. It has to do with the uniformity assumption.

 
## Hollingsworth sensitivity analysis: showing estimates the elevated hazard months of the acute phase
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


## Colored by acute duration to see if that affects estimates of
## EHM. Shown for one late scalar only.
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
    sel <- hf$var==var & hf$err==err & hf$het.sd==vv & hf$late.sc==lt.sc.temp # & hf$dur.ac==2 #
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
 
## LATE PHASE Hollingsworth sensitivity analysis: showing estimates the elevated hazard months of the late phase
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

## Does age become a spurious predictor of susceptibility via a cohort survival effect?
axis.cex <- 1.5
leg <- unique(wf$late.sc)
cols <- rep(NA,10)
cols[unique(wf$late.sc)] <- c('blue','orange','red')
legcols <- cols[unique(wf$late.sc)]
err.lab <- c('correct inclusion','excluded couples \nlost to follow-up by error')
pdf(file.path(outdir,'Wawer Analysis with Age etc.pdf'), w = 6, h = 6)
for(vv in unique(wf$het.sd)) {
  var.an <- 'ehm.ac'
  par(mfrow=c(2,2), mar=c(2.5,3,1,.5), oma = c(3,7.5,7,0))
  #layout(matrix(1:8,4,2))
  for(ee in 1:2) {
    for(cc in 1:2) { ## amount of correlation between measured confounders and underlying confounders
      cc.var <- c('','secp.age')[cc]
      cc.lab <- c('no covariates', 'secondary partner age \nat first visit')
      err <- c('base','XbErr')[ee]
      plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = '', ylab = '')
      sel <- wf$var==var.an & wf$err==err & wf$het.sd==vv & wf$hobs=='obsNA' & wf$cov==cc.var# & wf$dur.ac==2 #& wf$late.sc==10
      temp <- wf[sel,]
      ddply(temp, .(het.sd), with, {arrows(true+late.sc/10, lci, true+late.sc/10,uci, length=.05, angle = 90, code = 3,col = cols[late.sc])})    
      ddply(temp, .(het.sd), with, {points(true+late.sc/10, med, cex = .5, col = cols[late.sc], pch = 19)}) #, col = cols[het.sd+1])})
      abline(a=0, b=1, col = 'black', lty = 1, lwd = 1)
      ## segments(0, 75.4, 100, 75.4, col = 'red', lty = 2, lwd = 2)
      ## text(100, 75.4, 'Hollingsworth \nEstimate', pos = 4, cex = .5)
      segments(0, 36.25, 100, 36.25, col = 'red', lty = 3, lwd = 2)
      ## text(100, 36.25, 'Wawer \nEstimate', pos = 4, cex = .5)
      if(ee==1) mtext(cc.lab[cc], side = 3, line = 1, cex = 1.2, las = 1)
      if(ee==1 & cc==0) legend(5,ymax, title= 'late to chronic relative hazard \n(for months 20-11 prior to death)',
                   leg = leg, pch = 19, col=legcols, bty ='n', cex = 1)
      if(cc==1)  mtext(paste(err.lab[ee]), side = 2, line = 8)
    }
  }
  mtext(paste('indivual heterogeneity sd=',vv), side = 3, outer = T, line = 5, adj = .5, cex = axis.cex)
  mtext(expression(paste('true ',EHM[acute])), side = 1, outer = T, line = 1.5, cex = axis.cex)
  mtext(expression(paste('estimated ',EHM[acute])), side = 2, outer = T, line = .5, cex = axis.cex)
  ##mtext(expression(underline('proportion of heterogeneity observed and controlled for')), side = 3, outer = T, line = 2, cex = 1.2)
}
dev.off()

####################################################################################################
## Are any of the age or marital duration covariates significant predictors of hazard were included
## in a Poisson model?
leg <- c('index partner age', 'secondary partner age', 'couple duration',
         'secondary partner total duration of sexual activity',
         'secondary partner pre-couple duration of sexual activity')
pdf(file.path(outdir,'Wawer Analysis covariate estimates.pdf'), w = 8, h = 6)
plot(0,0, type='n', xlim = c(0, 3.5), ylim = c(3/5,5/3), bty = 'n', xlab = expression(paste('hazard heterogeneity ',sigma, ' (jittered)')), ylab = '', log='y')
abline(h=1, lwd = 3)
sel <- cef$cov!='' & cef$var=='obsNA' & cef$err=='base' & cef$late.sc==2
## sel <- cef$var==var.an & cef$err==err & cef$het.sd==vv & cef$hobs=='obsNA' & cef$cov==cc.var# & cef$dur.ac==2 #& cef$late.sc==10
temp <- cef[sel,]
## ddply(temp, .(het.sd), with, {arrows(true+late.sc/10, lci, true+late.sc/10,uci, length=.05, angle = 90, code = 3,col = cols[late.sc])})
xs <- jitter(temp$het.sd + as.numeric(temp$cov)/35, a = .05)
arrows(xs, temp$lci, xs, temp$uci, col = cols[as.numeric(temp$cov)], length=.05, angle = 90, code = 3)
points(xs, temp$med, col = 'black', pch = 19, cex = .3)
abline(, col = 'red')
##ddply(temp, .(het.sd), with, {points(true+late.sc/10, med, cex = .5, col = cols[late.sc], pch = 19)}) #, col = cols[het.sd+1])})
legend('topleft', leg = leg, pch = 19, col = cols[2:nlevels(temp$cov)], bty = 'n')
mtext('relative hazard per 5 years', side = 2, outer = T, line = 3, cex = axis.cex)
dev.off()
 
cols <- rep(NA,100)
cols[unique(cef$acute.sc)] <- colorRampPalette(c('purple','orange'))(length(unique(cef$acute.sc)))
## Now do on different plots
axis.cex <- 1.5
cols <- rainbow(7)
err.lab <- c('correct inclusion','excluded couples \nlost to follow-up by error')
labs <- c('index partner age', '2nd partner age', 'couple duration',
          '2nd partner \ntotal dur. of sex. activity',
          '2nd partner \npre-couple dur. of sex. activity')
pdf(file.path(outdir,'Wawer Analysis covariate estimates (panels).pdf'), w = 25, h = 12)
var.an <- 'ehm.ac'
par(mar=c(2.5,3,3,.5), oma = c(3,6,5,0))
layout(matrix(1:25,5,5))
for(bb in c(0,0.3,.5,.7,1)) { ## amount of correlation between measured confounders and underlying confounders
  obs <- paste0('obs',bb)
  for(cc in 1:5) { ## amount of correlation between measured confounders and underlying confounders
    cc.var <- levels(cef$cov)[-1][cc]
    err <- c('base','XbErr')[ee]
    plot(0,0, type='n', xlim = c(0, 3.3), ylim = c(3/5,5/3), bty = 'n', xlab = '', ylab = '', log='y',
         main = '')
    sel <- cef$cov==cc.var & cef$var==obs & cef$err=='base' & cef$late.sc==5 & cef$dur.ac==2 # &  cef$acute.sc==7
    temp <- cef[sel,]
    xs <- jitter(temp$het.sd, a = .1)
    arrows(xs, temp$lci, xs, temp$uci, col = cols[temp$acute.sc], length=.05, angle = 90, code = 3)
    points(xs, temp$med, col = cols[temp$acute.sc], pch = 19, cex = .8)
    abline(h=1, lwd = 2)
    ## legend('topleft', leg = leg, pch = 19, col = cols[2:nlevels(temp$cov)], bty = 'n')
    if(ee==1) mtext(labs[cc], side = 3, line = 1, cex = 1.2, las = 1)
    if(bb==0) mtext(labs[cc], side = 2, line = 3)
    if(ee==1 & cc==0) legend(5,150, title= 'late to chronic relative hazard \n(for months 20-11 prior to death)',
                 leg = c(2,5), pch = 19, col=cols[c(2,5)], bty ='n', cex = 1)
    if(cc==1) mtext(bb, side = 3, line = 3, cex = 1.2, las = 1)
  }
}
mtext(expression(paste(sigma, ' of inter-individual heterogeneity of hazard (jittered)')), side = 1, outer = T, line = 1.5, cex = axis.cex)
mtext('relative hazard per 5 years of:', side = 2, outer = T, line = 3, cex = axis.cex)
mtext(expression(underline('proportion of heterogeneity observed and controlled for')), side = 3, outer = T, line = 2, cex = 1.2)
dev.off()

## WTF is going on with age patterns not showing up??
ddply(cef[cef$err=='base' &cef$cov=='secp.pdsa' & cef$late.sc==5 & cef$dur.ac==2 & cef$acute.sc==7,], .(het.sd, var), summarise,lci = mean(lci),  med = mean(med), uci = mean(uci))

blocks[with(blocks, which(acute.sc==7 & dur.ac ==2 & het.gen.sd==2)),]



load(file.path('results','RakAcute','Uganda','Uganda-96200-1473.Rdata'))

rcohsim <- rak.wawer(rak.coh = cohsim, excl.extram=excl.extram, decont=decont, start.rak = 1994,
                     het.gen.sd = het.gen.sd,
                     verbose = T, browse=F)

dtt <- cohsim$dat

head(dtt)
mean(log(dtt$m.het.gen))
mean(log(dtt$f.het.gen))

## look at late phase CI's ranges
head(hf)
head(hf[hf$var=='late.sc',])
apply(hf[hf$var=='late.sc',c('lci','med','uci')],2, range)
pdf(file.path(outdir,'late test.pdf'))
temp <- hf[hf$var=='late.sc',]
hist(log(temp$lci))
hist(log(temp$med))
hist(log(temp$uci))
dev.off()

####################################################################################################
## What are plausible values for the true EHMacute given real estimates?
realan <- read.csv(file.path('results','HollingsworthAn','RealExclbyErr','output table.csv'))
rownames(realan) <- realan[,1]
realan <- realan[,-1]
realan$ehm.ac


## by CI's from analyzing real data
hf$plaus.rl <- NA
hf$plaus.rl[hf$var=='ehm.ac'] <- hf$med[hf$var=='ehm.ac'] < realan['97.5%','ehm.ac'] & hf$med[hf$var=='ehm.ac'] > realan['2.5%','ehm.ac']
## by CI's from analyzing simulated data
hf$plaus.sim <- NA
hf$plaus.sim[hf$var=='ehm.ac'] <- hf$lci[hf$var=='ehm.ac'] < realan['50%','ehm.ac'] & hf$uci[hf$var=='ehm.ac'] > realan['50%','ehm.ac']

## Plot plausible range of EHM.ac by heterogeneity assumption.
pdf(file.path(outdir,'plausible range (simCI).pdf'))
temp <- hf[hf$var=='ehm.ac' & hf$err=='XbErr',]
## plot(0,0, type = 'n', xlim = c(0,3), ylim = c(0,100), xlab = expression(paste('hazard heterogeneity ',sigma)), bty = 'n',
##      ylab = expression(paste('true ',EHM[acute])), main = expression(paste('true ',EHM[acute], ' consistent with Rakai Data')))
## points(temp$het.sd[temp$plaus.sim], temp$true[temp$plaus.sim])
with(temp[temp$plaus.sim,],
     boxplot(true ~ het.sd, ylim = c(0,100),xlab = expression(paste('hazard heterogeneity ',sigma)), bty = 'n',
      ylab = expression(paste('true ',EHM[acute])), main = expression(paste('true ',EHM[acute], ' consistent with Rakai Data')))
)
abline(h=realan['50%','ehm.ac'], col = 'red', lty = 2)
dev.off()

## Plot plausible range of EHM.ac by heterogeneity assumption.
pdf(file.path(outdir,'plausible range (realCI).pdf'))
temp <- hf[hf$var=='ehm.ac' & hf$err=='XbErr',]
## plot(0,0, type = 'n', xlim = c(0,3), ylim = c(0,100), xlab = expression(paste('hazard heterogeneity ',sigma)), bty = 'n',
##      ylab = expression(paste('true ',EHM[acute])), main = expression(paste('true ',EHM[acute], ' consistent with Rakai Data')))
## points(temp$het.sd[temp$plaus.sim], temp$true[temp$plaus.sim])
with(temp[temp$plaus.rl,],
     boxplot(true ~ het.sd, ylim = c(0,100),xlab = expression(paste('hazard heterogeneity ',sigma)), bty = 'n',
      ylab = expression(paste('true ',EHM[acute])), main = expression(paste('true ',EHM[acute], ' consistent with Rakai Data')))
   )
abline(h=realan['50%','ehm.ac'], col = 'red', lty = 2)
dev.off()

exp(qnorm(.975, 0, sd = 1))
exp(qnorm(.975, 0, sd = 2))
exp(qnorm(.975, 0, sd = 3))
