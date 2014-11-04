library(plyr); library(data.table); library(abind); library(multicore); library(coda); library(grid); library(utils)
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
outdir <- file.path('FiguresAndTables','UgandaFitSummaries')
nc <- 12                                       # core per simulation
load(file='results/RakAcute/wf.Rdata') ## big file, stored outside of git repo
load(file = file.path('results','HollingsworthAn','RealExclbyErr','workspace.Rdata')) ## fit to real data
names(fout)
outdir <- file.path('FiguresAndTables','UgandaFitSummaries') ## incase workspace load changed it

## Get summary estimates of ehm.ac from our refit of Hollingsworth's model to Rakai data
head(fout$exposts)
fexsum <- apply(fout$expost, 2, function(x) quantile(x, c(.025,.5,.975)))

## Get summary estimates of ehm.ac from Wawer raw estimates with no heterogeneity
infs <- c(10, 36) ## infections
cas <- c(1221, 48525)  ## coital acts
prev <- infs/cas
rh.waw <- prev[1]/prev[2]

## Wawer CIs
## Multivariate from Wawer paper, direct from Table 2
wexsum <- (c(3.05, 7.25, 17.25)-1)*5
## Univariate estimates from Wawer data using binomial regression, ORs are approximately equal to
## prevalence ratios because probabilities are so small so we can use logistic regression.
wraw <- data.frame(ph = c('acute','chronic'), infs = infs, cas = cas)
wraw$ph <- relevel(wraw$ph, levels = c('acute','chronic'), ref = 'chronic')
mod <- glm(infs ~ offset(log(cas)) + ph, wraw, family= poisson)
wexsum.univ <- c(exp(confint(mod))['phacute',1],
               exp(coef(summary(mod))['phacute','Estimate']),
               exp(confint(mod))['phacute',2])
wexsum.univ <- (wexsum.univ-1)*5

####################################################################################################
## Do LOESS for each sigma_hazard for Waw & Holl with XbErr, then find where previously published
## CI's hit the Loess.
thf <- hf[hf$var=='ehm.ac' & hf$ehm.late==40,] ## First, just get subsets we're interested in.
twf <- wf[wf$var=='ehm.ac' & wf$ehm.late==40 & wf$cov=='',]
hsds <- seq(0,3, by = .5)
  
## Wawer
xs <- seq(0,250, by = .1)
wpred <- data.frame(xs = xs)
for(hsd in hsds) {
  sel <- twf$err=='XbErr' & twf$het.sd==hsd & twf$hobs=='obs0'
  ## Add a column that has the predicted estimated EHM for a true EHM (do it at all xs, insted of the subset of EHMs simulated)
  wpred <- cbind(wpred, with(twf[sel,], predict(loess(med~true), xs)))
  names(wpred)[ncol(wpred)] <- paste0('sig',hsd)
}
head(wpred)

## Holl
hpred <- data.frame(xs = xs)
for(hsd in hsds) {
  sel <- thf$err=='XbErr' & thf$het.sd==hsd
  ## Add a column that has the predicted estimated EHM for a true EHM (do it at all xs, insted of the subset of EHMs simulated)
  hpred <- cbind(hpred, with(thf[sel,], predict(loess(med~true), xs)))
  names(hpred)[ncol(hpred)] <- paste0('sig',hsd)
}
head(hpred)

findCI.h <- function(x) { ## find value closest to median & 95CI bounds
  lci <- which.min(abs(x-fexsum['2.5%','ehm.ac']))
  med <- which.min(abs(x-fexsum['50%','ehm.ac']))
  uci <- which.min(abs(x-fexsum['97.5%','ehm.ac']))
  return(xs[c(lci, med, uci)])
}

findCI.w <- function(x) {               ## wawer multivariate
  lci <- which.min(abs(x-wexsum[1]))
  med <- which.min(abs(x-wexsum[2]))
  uci <- which.min(abs(x-wexsum[3]))
  return(xs[c(lci, med, uci)])
}

findCI.w.univ <- function(x) { ## wawer univariate
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

load('FiguresAndTables/VL Profile/ehms.vl.Rdata')
## Literature estimate table
tab1 <- rbind(ehms.vl[-1], c(20.9, 50.2,102), c(11.9, 36.3, 96), wexsum, fexsum[,'ehm.ac'], c(NA, (30.3-1)*4.8, NA), c(NA, 238, NA))
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
pdf(file.path(outdir, paste0('Figure 5 - CIs Wawer ', modlab, '.pdf')), w = ifelse(add.leg, 6.83, 7.5), h = ifelse(add.leg, 2.5, 4.5))
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


## Literature estimate table
tab1 <- rbind(ehms.vl[-1], c(20.9, 50.2,102), c(11.9, 36.3, 96), wexsum, fexsum[,'ehm.ac'], c(NA, (30.3-1)*4.8, NA), c(NA, 238, NA))
tab1 <- data.frame(study = NA, tab1)
colnames(tab1)[-1] <- c('lci','med','uci')
tab1$study <- c('viral load estimate', 'Wawer et al. (unadjusted)',
                'Wawer et al. (coital acts)', 'Wawer et al. (full model)',
                'Hollingsworth et al.', 'Powers et al.',
                'Rasmussen et al.')
tab1 <- rbind(tab1, c(NA, cis[,'sig1']))
tab1$study[8] <- 'Rakai re-estimate'
tab1 <- tab1[c(1,8,2:7),]

####################################################################################################
## Simple Fig 5
for(do.log in 0:1) {
    to.do <- c(1:2,5,3,6:7)
    pdf(file.path(outdir, paste0('Fig 5 simple', 'log'[do.log], '.pdf')), w = 6, h = 5)
    par(mar=c(10,5,.5,.5), mgp = c(2.4,1,0),'ps'=12)
##################################################
    xmax <- 100
    ymax <- 150
    ylab <- expression(EHM[acute])
    ylab <- 'excess hazard months \nattributable to the acute phase'
    if(do.log) {
        plog <- 'y'
        ylim <- c(1,ymax)
    }else{
        plog <- ''
        ylim <- c(0,ymax)
    }
    plot(0,0, type = 'n', xlab = '', ylab = ylab, bty = 'n', log = plog,
         main = '', xlim = c(1,length(to.do)), ylim = ylim, axes=F)
    axis(1, 1:length(to.do), lab = paste0('(',LETTERS[1:length(to.do)], ') ', tab1[to.do,'study']), las = 2)
    if(!do.log) {
        axis(2, seq(0,ymax, by = 50), las = 2)
        axis(2, seq(0,ymax, by = 10), lab=NA)
    }else{
        axis(2, c(1,5, 10, 50, 100, 150), las = 2)
        axis(2, c(1:9, seq(10,100,10)), lab=NA)
    }
    par('xpd'=F)
    abline(h=seq(0,ymax,by=10), col = gray(.9))
    par('xpd'=T)
    for(ii in 1:length(to.do))  {
        dd <- to.do[ii]
        if(ii<3) col <- 'black' else col <- gray(.6)
        points(ii, tab1[dd,'med'], pch = 15, col = col, cex = ifelse(do.log, 1.5,1))
        arrows(ii, tab1[dd,'lci'], ii, tab1[dd,'uci'], code = 3, len = .05, angle = 90, col = col, lwd = 2)
        ##text(ii, tab1[ii,'uci'], tab1[ii,'study'], pos = 3)
    }
    graphics.off()
}

####################################################################################################
## 95% contour for RH[acute] & d[acute]
## Which values produced things within Wawer's 95% CI?

sel <- with(twf, err=='XbErr' & het.sd==1 & hobs=='obs0', ehm.late==40)
sum(sel)
temp <- twf[sel,]
temp <- temp[temp$med > wexsum[1] & temp$med < wexsum[3],c('acute.sc','dur.ac')]

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor, poly=F){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  if(!poly) lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  if(poly) polygon(c(xcoord[hpts], rev(xcoord[hpts])), c(ycoord[hpts], rev(ycoord[hpts])), col = lcolor)
}  

library(plotrix)

pdf(file=file.path(outdir,"95CI region.pdf"), w = 3.27, h = 3)
ct <- .8
layout(matrix(c(1:2),1,2), w = c(1,.3))
par(mar=c(4,4,1,1))
xs <- seq(-.3,7.2, by = .05) ##log(seq(0.1,8, by = 1))
ys <- seq(-2.4,2.7, by = .05) ##log(seq(0.1, 11, by = .5))
levels <- c(0,2,5,10,25,50,70,100,500,1000)
cols <- colorRampPalette(c('purple','orange','red'))(length(levels)-1)
image(xs, ys, outer(exp(xs)-1,exp(ys)), breaks = levels, xlim = c(-.3,log(60)), ylim = c(-2.4,2.35), mgp = c(3,0,0),
               col = cols, axes = F, xlab = expression(paste(RH['acute'])), ylab = expression(d[acute]))
points(log(temp[,1]), log(temp[,2]), pch = 16)
xts <- c(1:9, seq(10, 100, by = 10))
xsh <- c(1,10,50)
xls <- xts
xls[!xls %in% xsh] <- NA
axis(1, at = log(xts), lab = xls)
yts <- c(seq(.1,.9, by = .1),1:10)
ysh <- c(.1,1,10)
yls <- yts
yls[!yls %in% ysh] <- NA
axis(2, at = log(yts), lab = yls, las = 2)
## Palette legend
par(mar=rep(0,4), 'ps'=12)
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
color.legend(.09,.1,.15,.8, levels, rect.col = cols, gradient = "y", cex = ct*.8)
text(.025, .8, expression(paste(EHM['acute'])), pos = 3, cex = ct*1.2)
graphics.off()


write.csv(temp, file = file.path(outdir,'RHacute & dacute combos with 95CI.csv'))
