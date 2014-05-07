library(plyr); library(data.table); library(abind); library(multicore); library(coda); library(grid)
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

modlab <- c('Univar','Multivar')
for(mm in 1:2) {
  if(mm==1) {## univar
    cis <- wcis.univ
    wsum <- wexsum.univ
  }else{ ## multivar
    cis <- wcis
    wsum <- wexsum
  }
  pdf(file.path(outdir, paste0('CIs Wawer ', modlab[mm], '.pdf')), w = 6.83, h = 2.5)
  cols <- c('purple','red','blue','orange')
  ## cols <- rainbow(length(hsds))
  hsds <- 0:3
  ylab <- expression(paste('estimated ',EHM[acute]))
  xlab <- expression(paste('true ',EHM[acute]))
  par(mar=c(4,4,2,1.5), mgp = c(2.4,1,0),'ps'=12)
  layout(matrix(1:3, nr=1, nc=3), w = c(1,1,.6))
##################################################
  xmax <- 100
  ymax <- 250
  plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = xlab, ylab = ylab, main = '(A)', las = 2)
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
  text(xmax*1.1, wsum[3], '(i)', pos = 3) ## 
  ## text(xmax+2, wsum[3]+1, 'point estimate', pos = 3)
  ## text(xmax, wsum[2], 'point estimate', pos = 4)
  ## text(xmax, wsum[1], '95% CI lower bound', pos = 4)
  ## text(xmax, wsum[3], '95% CI upper bound', pos = 4)
  ## Show all CIs
  par(mar=c(4,4,2,0))
  plot(0,0, type = 'n', xlab = bquote(sigma['hazard']), ylab = xlab, bty = 'n',
       main = '(B)', xlim = c(-.2,6.2), ylim = c(0,ymax), axes=F)
  axis(1, hsds)
  axis(2, seq(0,ymax, by = 50), las = 2)
  axis(2, seq(0,ymax, by = 10), lab=NA)
  par('xpd'=F)
  abline(h=seq(0,ymax,by=10), col = gray(.9))
  par('xpd'=T)
  ## Holl
  ## for(ii in 1:length(hsds)) arrows(hsds[ii], hcis[1,ii], hsds[ii], hcis[3,ii], code = 3, len = .05, angle = 90, col = cols[ii])
  ## for(ii in 1:length(hsds)) points(hsds[ii], hcis[2,ii], pch = 19, col = cols[ii])
  ## Wawer
  for(ii in 1:length(hsds)) arrows(hsds[ii], cis[1,ii], hsds[ii], cis[3,ii], code = 3, len = .05, angle = 90, col = cols[ii])
  for(ii in 1:length(hsds)) points(hsds[ii], cis[2,ii], pch = 19, col = cols[ii])
  ## Wawer 2005 multivariate
  arrows(3.5, wsum[3], 3.5, wsum[1], code = 3, len = .05, angle = 90, col = gray(.3))
  points(3.5, wsum[2], pch = 19, col = gray(.3))
  text(3.5, wsum[3], '(i)', pos = 3)    ## "Wawer \n(multivariate)",
  ## Holl our refit
  arrows(4, fexsum['97.5%','ehm.ac'], 4, fexsum['2.5%','ehm.ac'], code = 3, len = .05, angle = 90, col = gray(.3))
  points(4, fexsum['50%','ehm.ac'], pch = 19, col = gray(.3))
  text(4, fexsum['97.5%','ehm.ac'], "(ii)", pos = 3) ##"Hollingsworth \n(our refit)
  ## Holl 2008 original
  points(4.5, (26-1)*2.9, pch = 15)
  text(4.5, (26-1)*2.9, '(iii)', pos = 3) ## "Hollingsworth \n(original estimate)"
  ## Powers 2011
  points(5, (30.3-1)*4.8, pch = 15)
  text(5, (30.3-1)*4.8, '(iv)', pos = 3) ## "Powers"
  ## Rasmussen 2014
  points(5.5, (20-1)*12, pch = 15)
  text(5.5, (20-1)*12, '(v)', pos = 3) ## "Rasmussen"
  axis(1, c(3.25,5.65), lab = c('','previous \nestimates'))
  ## Expected based on viral load
  ##arrows(3.5, wsum.univ[3], 3.5, wsum.univ[1], code = 3, len = .05, angle = 90, col = gray(.3))
  points(6, 4.3, pch = 15, col = gray(.3))
  text(6, 4.3, '(vi)', pos = 3)    ## "Wawer \n(multivariate)",
#### legend
  par(mar=rep(0,4))
  plot(0,0, type = 'n', xlab ='', ylab='', axes =F)
  labs <- c('(i) Wawer et al. 2005 \nmultivariate', '(ii) Hollingsworth et al. 2008 \n(our refit)',
            '(ii) Hollingsworth et al. 2008 \n(original)', '(iv) Powers et al. 2011',
            '(v) Rasmussen et al. 2014', '(vi) expected based on viral loads')
  mtext('Legend', side = 3, line = -1, adj = .5, cex = .7)
  for(ii in 1:6) mtext(labs[ii], side = 3, line = -ii*3, adj = .5, cex = .7)
  graphics.off()
}
