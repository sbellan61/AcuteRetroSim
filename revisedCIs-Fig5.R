library(plyr); library(data.table); library(abind); library(multicore); library(coda); library(grid); library(utils)
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')
nc <- 12                                       # core per simulation
load(file=file.path(outdir, 'wf.Rdata'))
load(file = file.path('results','HollingsworthAn','RealExclbyErr','workspace.Rdata')) ## fit to real data
names(fout)

## Get summary estimates of ehm.ac from our refit of Hollingsworth's model to Rakai data
head(fout$exposts)
fexsum <- apply(fout$expost, 2, function(x) quantile(x, c(.025,.5,.975)))

## Get summary estimates of ehm.ac from Wawer raw estimates with no heterogeneity
infs <- c(10, 36) ## infections
cas <- c(1221, 48525)  ## coital acts
prev <- infs/cas
rh.waw <- prev[1]/prev[2]

## Wawer CIs
## Multivariate from their paper
wexsum <- (c(3.05, 7.25, 17.25)-1)*5
## Univariate from binomial regression, ORs are approximately equal to prevalence ratios because
## probabilities are so small so we can use logistic regression.
wraw <- data.frame(ph = c('acute','chronic'), infs = infs, cas = cas)
wraw$ph <- relevel(wraw$ph, levels = c('acute','chronic'), ref = 'chronic')
mod <- glm(infs ~ offset(log(cas)) + ph, wraw, family= poisson)
wexsum.univ <- c(exp(confint(mod))['phacute',1],
               exp(coef(summary(mod))['phacute','Estimate']),
               exp(confint(mod))['phacute',2])
wexsum.univ <- (wexsum.univ-1)*5

####################################################################################################
## Do LOESS for each sigma_hazard for Waw & Holl with XbErr, then find where true CI's hit the Loess
## just get subsets we're interested in
thf <- hf[hf$var=='ehm.ac' & hf$ehm.late==40,] 
twf <- wf[wf$var=='ehm.ac' & wf$ehm.late==40 & wf$cov=='',]
hsds <- seq(0,3, by = .5)
  
## Wawer
xs <- seq(0,250, by = .1)
wpred <- data.frame(xs = xs)
for(hsd in hsds) {
  sel <- twf$err=='XbErr' & twf$het.sd==hsd & twf$hobs=='obs0'
  wpred <- cbind(wpred, with(twf[sel,], predict(loess(med~true), xs)))
  names(wpred)[ncol(wpred)] <- paste0('sig',hsd)
}
head(wpred)

## Holl
hpred <- data.frame(xs = xs)
for(hsd in hsds) {
  sel <- thf$err=='XbErr' & thf$het.sd==hsd
  hpred <- cbind(hpred, with(thf[sel,], predict(loess(med~true), xs)))
  names(hpred)[ncol(hpred)] <- paste0('sig',hsd)
}
head(hpred)

findCI.h <- function(x) {
  lci <- which.min(abs(x-fexsum['2.5%','ehm.ac']))
  med <- which.min(abs(x-fexsum['50%','ehm.ac']))
  uci <- which.min(abs(x-fexsum['97.5%','ehm.ac']))
  return(xs[c(lci, med, uci)])
}

findCI.w <- function(x) {
  lci <- which.min(abs(x-wexsum[1]))
  med <- which.min(abs(x-wexsum[2]))
  uci <- which.min(abs(x-wexsum[3]))
  return(xs[c(lci, med, uci)])
}

findCI.w.univ <- function(x) {
  lci <- which.min(abs(x-wexsum.univ[1]))
  med <- which.min(abs(x-wexsum.univ[2]))
  uci <- which.min(abs(x-wexsum.univ[3]))
  return(xs[c(lci, med, uci)])
}

hcis <- signif(apply(hpred[,-1], 2, findCI.h),3)
wcis <- signif(apply(wpred[,-1], 2, findCI.w),3)
wcis.univ <- signif(apply(wpred[,-1], 2, findCI.w.univ),3)

ci.frame <- data.frame(hsds, t(hcis),t(wcis), t(wcis.univ))
names(ci.frame) <- c('sigma','hlci','hmed','huci','wlci.m','wmed.m','wuci.m','wlci.u','wmed.u','wuci.u')
write.csv(ci.frame, file = file.path(outdir, 'CI data frame.csv'))

## Literature estimate table
tab1 <- rbind(c(NA, 4.3, NA), c(20.9, 50.2,102), c(11.9, 36.3, 96), wexsum, fexsum[,'ehm.ac'], c(NA, (30.3-1)*4.8, NA), c(NA, 238, NA))
tab1 <- data.frame(study = NA, tab1)
colnames(tab1)[-1] <- c('lci','med','uci')
tab1$study <- c('based on VL', 'Wawer (raw)', 'Wawer (coital acts)', 'Wawer (full model)', 'Hollingsworth (refit)', 'Powers et al. 2011',
                'Rasmussen et al. 2014')

hsds <- 0:3
add.leg <- F
x.offset <- 4
modlab <- c('Multivar')
cis <- wcis
wsum <- wexsum
pdf(file.path(outdir, paste0('Figure 5 - CIs Wawer ', modlab[mm], '.pdf')), w = ifelse(add.leg, 6.83, 7.5), h = ifelse(add.leg, 2.5, 4.5))
cols <- c('purple','red','blue','orange')
## cols <- rainbow(length(hsds))
ylab <- expression(paste('estimated ',EHM[acute]))
xlab <- expression(paste('true ',EHM[acute]))
if(!add.leg) { ## presentations
  to.do <- 1:6
  mains <- c('Wawer Model', expression(paste(EHM[acute], ' estimates')))
  par(mar=c(4,4,2,2), mgp = c(2.7,1,0),'ps'=12, oma = c(3.5,0,0,.5))
  layout(matrix(1:2, nr=1, nc=2), w = c(1,1.3))
  ymax <- 150
  spc <- 1.2
}else{ ## publication
  to.do <- c(1,4,5,6,7)
  mains <- paste0('(',LETTERS[1:2],')')
  par(mar=c(4,4,2,2), mgp = c(2.7,1,0),'ps'=13, oma = c(0,0,0,0))
  layout(matrix(1:3, nr=1, nc=3), w = c(1,1,.6))
  ymax <- 250
}
##################################################
xmax <- 100
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = xlab, ylab = ylab, main = mains[1], las = 2)
## polygon(c(0,xmax,xmax,0), rep(fexsum[c('2.5%','97.5%'),'ehm.ac'], each = 2), col = gray(.5), border = NA) ## Holl
## segments(0,fexsum['50%','ehm.ac'], xmax, fexsum['50%','ehm.ac'], col = gray(.3))
polygon(c(0,xmax,xmax,0), rep(wsum[c(1,3)], each = 2), col = gray(.5), border = NA) ## Waw
segments(0,wsum[2], xmax, wsum[2], col = gray(.3))
abline(a=0,b=1, lwd=2)
for(hh in 1:length(hsds)) {
  if(hsds[hh] %in% 0:3) {
    ## sel <- thf$err=='XbErr' & thf$het.sd==hsds[hh] ##2
    ## with(thf[sel,], lines(xs, predict(loess(med~true), xs), col = cols[hh], lty = 1, lwd = 1)) # holl
    ## points(hcis[,paste0('sig',hsds[hh])], fexsum[c('2.5%','50%','97.5%'),'ehm.ac'], pch = 19, col = cols[hh], cex = .4)
    sel <- twf$err=='XbErr' & twf$het.sd==hsds[hh] & twf$hobs=='obsNA'##2
    with(twf[sel,], lines(xs, predict(loess(med~true), xs), col = cols[hh], lty = 1, lwd = 1)) # waw
    points(cis[,paste0('sig',hsds[hh])], wsum, pch = 19, col = cols[hh], cex = .4)
  }
}
par('xpd'=T)
arrows(xmax*1.1, wsum[3], xmax*1.1, wsum[1], code = 3, len = .05, angle = 90, col = gray(.3))
points(xmax*1.1, wsum[2], pch = 19, col = gray(.3))
if(add.leg)   text(xmax*1.1, wsum[3], '(II)', pos = 3) ## 
## text(xmax+2, wsum[3]+1, 'point estimate', pos = 3)
## text(xmax, wsum[2], 'point estimate', pos = 4)
## text(xmax, wsum[1], '95% CI lower bound', pos = 4)
## text(xmax, wsum[3], '95% CI upper bound', pos = 4)
## Show all CIs
par(mar=c(4,5,2,0.5))
plot(0,0, type = 'n', xlab = '', ylab = xlab, bty = 'n',
     main = mains[2], xlim = c(-.2,x.offset + (length(to.do)*spc)), ylim = c(0,ymax), axes=F)
mtext(bquote(sigma['hazard']), side =  1, line = 3, adj = .15)
axis(1, hsds)
par('xpd'=F)
abline(h=seq(0,ymax,by=10), col = gray(.9))
par('xpd'=T)
axis(2, seq(0,ymax, by = 50), las = 2)
axis(2, seq(0,ymax, by = 10), lab=NA)
## Holl
## for(ii in 1:length(hsds)) arrows(hsds[ii], hcis[1,ii], hsds[ii], hcis[3,ii], code = 3, len = .05, angle = 90, col = cols[ii])
## for(ii in 1:length(hsds)) points(hsds[ii], hcis[2,ii], pch = 19, col = cols[ii])
## Wawer
for(ii in 1:length(hsds)) arrows(hsds[ii], cis[1,paste0('sig',hsds[ii])], hsds[ii], cis[3,paste0('sig',hsds[ii])], code = 3, len = .05, angle = 90, col = cols[ii])
for(ii in 1:length(hsds)) points(hsds[ii], cis[2,paste0('sig',hsds[ii])], pch = 19, col = cols[ii])
for(ii in 1:length(to.do))
  {
    dd <- to.do[ii]
    iid <- ii*spc
    points(iid+x.offset, tab1[dd,'med'], pch = 19, col = gray(.3))    
    arrows(iid+x.offset, tab1[dd,'lci'], iid+x.offset, tab1[dd,'uci'], code = 3, len = .05, angle = 90, col = gray(.3))
    if(add.leg) text(iid+x.offset, max(tab1[dd,-1],na.rm=T), paste0('(',as.roman(ii), ')'), pos = 3)
  }
if(add.leg) {
  axis(1, c(x.offset+1, x.offset + length(to.do)), lab = NA)
  text(x.offset + length(to.do)/2 + .5, -30, 'other estimates')
}else{
  axis(1, x.offset+ 1:length(to.do)*spc, lab = tab1[to.do,'study'], las = 2)
}
if(add.leg) {
  par(mar=rep(.5,4))
  plot(0,0, type = 'n', xlab ='', ylab='', axes =F, bty = 'n')
  box()
  labs <- paste0('(',as.roman(1:length(to.do)), ') ', tab1$study[to.do])
  mtext('Legend', side = 3, line = -1, adj = .5, cex = .6)
  for(ii in 1:6) mtext(labs[ii], side = 3, line = -ii*2.8, adj = .5, cex = .6)
}
graphics.off()

  
####################################################################################################
## EHM literature estimates for presentation
to.do <- 1:6
pdf(file.path(outdir, paste0('published EHMs.pdf')), w = 6, h = 5)
par(mar=c(10,4,.5,.5), mgp = c(2.4,1,0),'ps'=12)
##################################################
xmax <- 100
ymax <- ceiling(max(tab1[to.do,-1],na.rm=T)/50)*50
plot(0,0, type = 'n', xlab = '', ylab = xlab, bty = 'n',
     main = '', xlim = c(1,length(to.do)), ylim = c(0,ymax), axes=F)
axis(1, 1:length(to.do), lab = tab1[to.do,'study'], las = 2)
axis(2, seq(0,ymax, by = 50), las = 2)
axis(2, seq(0,ymax, by = 10), lab=NA)
par('xpd'=F)
abline(h=seq(0,ymax,by=10), col = gray(.9))
par('xpd'=T)
for(ii in 1:length(to.do))
  {
    dd <- to.do[ii]
    points(ii, tab1[ii,'med'], pch = 19, col = gray(.3))    
    arrows(ii, tab1[ii,'lci'], ii, tab1[ii,'uci'], code = 3, len = .05, angle = 90, col = gray(.3))
    ##text(ii, tab1[ii,'uci'], tab1[ii,'study'], pos = 3)
  }
graphics.off()

