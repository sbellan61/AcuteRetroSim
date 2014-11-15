library(plyr); library(data.table); library(abind); library(multicore); library(coda); library(grid); library(utils)
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
outdir <- file.path('FiguresAndTables','UgandaFitSummaries')
nc <- 12                                       # core per simulation
load(file='results/RakAcute/wf.Rdata') ## big file, stored outside of git repo
load(file = file.path('results','HollingsworthAn','RealExclbyErr','workspace.Rdata')) ## fit to real data
names(fout)
outdir <- file.path('FiguresAndTables','UgandaFitSummaries') ## incase workspace load changed it
finalbatch <- 4 ## final abc batch to use as posterior
abcdir <- 'results/abcSummary'
load(file=file.path(abcdir, paste0('IntermedDistr',finalbatch,'.Rdata'))) ## Load last distribution (already filtered)

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

wexsum
wexsum.univ

load('FiguresAndTables/VL Profile/ehms.vl.Rdata')
## Literature estimate table
tab1 <- rbind(ehms.vl[-1], c(20.9, 50.2,102), c(11.9, 36.3, 96), wexsum, fexsum[,'ehm.ac'], c(NA, (30.3-1)*4.8, NA), c(NA, 238, NA))
tab1 <- data.frame(study = NA, tab1)
colnames(tab1)[-1] <- c('lci','med','uci')
tab1$study <- c('based on VL', 'Wawer (raw)', 'Wawer (coital acts)', 'Wawer (full model)', 'Hollingsworth (refit)', 'Powers et al. 2011',
                'Rasmussen et al. 2014')

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

tab1 <- rbind(tab1, c('Bellan et al. Rakai re-estimate', quantile(pmatChosen$EHMacute, c(.025,.5,.975))))
tab1 <- tab1[c(1,8,2:7),]
tab1[,-1] <- apply(tab1[,-1], 2, as.numeric)
tab1[,-1] <- signif(tab1[,-1],2)

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
    plot(1,1, type = 'n', xlab = '', ylab = ylab, bty = 'n', log = plog,
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
        ys <- as.numeric(tab1[dd,c('lci','uci')])
        if(do.log==1) ys[ys<0] <- .01
        arrows(ii, ys[1], ii, ys[2], code = 3, len = .05, angle = 90, col = col, lwd = 2)
        ##text(ii, tab1[ii,'uci'], tab1[ii,'study'], pos = 3)
    }
    graphics.off()
}

####################################################################################################
## 95% contour for RH[acute] & d[acute]
## Which values produced things within Wawer's 95% CI?

library(plyr); library(data.table); library(abind); library(multicore); library(emdbook);library(coda); library(plotrix)

## Get 95% posterior contour using volume of points in that region
test <- HPDregionplot(pmatChosen, vars = c('logacute.sc','logdur.ac'), prob = c(.95), lims = c(-1.5,10,-2,11), n = 27, h = c(.3,.1),
                      bty = 'n', xlab = expression(RH['acute']), ylab = 'acute phase duration (months)')
graphics.off()
####################################################################################################
## Figure 1. Show EHM diagram and RH[acute]/d[acute] collinearity plot with MCMC posteriors from refit of Holl Mod.
for(ff in 1:2) {
    if(ff==1) pdf(file.path(outdir,'Figure SX - fitted RHacute dacute.pdf'), w = 3.27, h = 3)
    if(ff==2) png(file.path(outdir,'Figure SX - fitted RHacute dacute.png'), w = 3.27, h = 3, units='in',res=200)
    par('ps'=10)
    layout(matrix(c(1:2),1,2), w = c(1,.3))
##################################################
    ## contour plot of EHM
    par(mar=c(4,4,.3,0))
    xs <- seq(-1,7.2, by = .05) ##log(seq(0.1,8, by = 1))
    ys <- seq(-2.4,2.7, by = .05) ##log(seq(0.1, 11, by = .5))
    levels <- c(0,2,5,10,25,50,70,100,200,500,1000,10000,10^5)
    cols <- colorRampPalette(c('purple','orange','red'))(length(levels)-1)
    blevels <- c(-25,-10,-2,0)
    bcols <- colorRampPalette(c('dark green','purple'))(length(blevels))
    cols <- c(bcols[-length(bcols)], cols)
    levels <- c(blevels[-length(blevels)],levels)
    image(xs, ys, outer(exp(xs)-1,exp(ys)), breaks = levels, xlim = c(-1,6.2), ylim = c(-2.4,2.6), mgp = c(3,0,0),
          col = cols, axes = F, xlab = expression(paste(RH['acute'])), ylab = expression(d[acute]))
    xts <- c(1:9, seq(10, 90, by = 10), seq(100, 500, by = 100))
    xsh <- c(1,10,100)
    xls <- xts
    xls[!xls %in% xsh] <- NA
    axis(1, at = log(xts), lab = xls)
    yts <- c(seq(.1,.9, by = .1),1:10)
    ysh <- c(.1,1,10)
    yls <- yts
    yls[!yls %in% ysh] <- NA
    axis(2, at = log(yts), lab = yls, las = 2)
                                        #lines(test[[1]]$x, test[[1]]$y)
    with(pmatChosen, points(log(acute.sc), log(dur.ac), cex = .3, col = 'black', pch = 16))
    ## text(log(meds['acute.sc']), log(meds['dur.ac']), 'median', pos=1, cex=ct)
    ## with(data.frame(outtab), axis(1, at = log(acute.sc[c(1,3)]), lab = NA, line = -1.5))
    ## with(data.frame(outtab), axis(2, at = log(dur.ac[c(1,3)]), lab = NA, line = -2)) 
    ## Palette legend
    par(mar=rep(0,4))
    plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
    color.legend(.09,.1,.15,.8, levels, rect.col = cols, gradient = "y", cex = .8)
    text(.025, .8, expression(paste(EHM['acute'])), pos = 3, cex = 1.2)
}
graphics.off()
####################################################################################################



pdf(file.path(outdir,'Figure SX - fitted RHacute dacute.pdf'), w = 4, h = 4)
with(pmatChosen, plot(acute.sc, dur.ac, cex = .3, log='xy', xlim = c(.5,200), ylim = c(.5,10), las = 3))
graphics.off()
