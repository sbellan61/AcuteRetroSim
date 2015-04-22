library(plyr); library(data.table); library(abind); library(coda); library(plotrix); library(sp); library(emdbook); library(multicore); 
rm(list=ls(all=T)); gc()
source('RakFunctions.R')
####################################################################################################
## Show that duration & infectivity of acute phase are not identifiable
## Create Figure 1.
####################################################################################################
## Load our MCMC refit of Holl Model to the real data (with exclusion criteria error)
load(file = file.path('HollModRefit','HollModRefit.Rdata'))
outtab <- outtab[,,'XbErr']
print(outtab) ## Holl Mod refit with MCMC
fexp <- (fout$exposts) ## posteriors
fexp$bp.ac <- fexp$bp * fexp$acute.sc ## acute hazard = chronic hazard * RH[ac]
quantile(fexp$bp.ac, c(0.025, .5, .975))*12 ## acute hazad per peson-year
fmcmc <- as.mcmc(fout$posts)
meds <- apply(fout$exposts[,c('acute.sc','dur.ac')], 2, median) ## medians of these parameters
apply(fout$exposts[,c('acute.sc','dur.ac')], 2, function(x) quantile(x, c(.025,.5,.975)))

outdir <- file.path('FiguresAndTables','HollingsworthAn')
if(!file.exists(outdir)) dir.create(outdir)

## Get 95% posterior contour using volume of points in that region. USING 90% because HPD tends to
## overshoot a bit, true values for EHMacute CI's given in text
HollCont <- HPDregionplot(fmcmc, vars = c('acute.sc','dur.ac'), prob = c(.95), lims = c(-1.5,10,-2,11), #n = 28, h = c(1,.08),
                      bty = 'n', xlab = expression(RH['acute']), ylab = 'acute phase duration (months)')

graphics.off()

HollXY <- with(HollCont[[1]], cbind(x=x,y=y))
HollXY <- rbind(HollXY, tail(HollXY,1))
mean(point.in.polygon(fmcmc[,'acute.sc'],fmcmc[,'dur.ac'],HollXY[,'x'], HollXY[,'y']))

####################################################################################################
## ABC-SMC results
####################################################################################################
abcdir <- 'results/abcSummary'
finalbatch <- 5
ABCoutdir <- 'FiguresAndTables/abcFig'
load(file=file.path(abcdir, paste0('IntermedDistr',finalbatch,'.Rdata'))) ## Load last distribution (already filtered)
## Get 95% posterior contour using volume of points in that region, .89 chosen based on point.in.polygon below and to visually match the true EHMacute quantile (which has 95% upper bound of 63, so contour can't go over the 70 color boundary
abcCont <- HPDregionplot(pmatChosen, vars = c('logacute.sc','logdur.ac'), prob = c(.87), lims = c(-1.5,10,-2,11), n = 28, h =c(.8,1),# h = c(.3,.1),
                      bty = 'n', xlab = expression(RH['acute']), ylab = 'acute phase duration (months)')
abcCIs <- log(quants(pmatChosen[,c('acute.sc','dur.ac')]))
abcEHM <- quants(pmatChosen[,'EHMacute',drop=F])
graphics.off()

abcXY <- with(abcCont[[1]], cbind(x=x,y=y))
abcXY <- rbind(abcXY, tail(abcXY,1))
mean(point.in.polygon(pmatChosen$logacute.sc,pmatChosen$logdur.ac,abcXY[,'x'], abcXY[,'y']))

drawHazProf <- function(yt, acute.sc = 101, dur.ac = .75, examnum = 1, lab=T) {
    polygon(c(0,dur.ac,dur.ac,0), yt + c(0,0,acute.sc,acute.sc), col = col.ac, border=NA) ## acute
    polygon(c(100,110,110,100), yt + c(0,0,7,7),col = col.lt, border=NA) ## late
    polygon(c(0,110,110,0), yt + c(0,0,1,1), col = col.ch, border=NA) ## chronic
    segments(110,yt,120,yt, lty = 3, col = col.aids)
    points(-8,yt+5, pch = 21, cex = 2.5)
    points(-8,yt+5, pch = as.character(examnum))
    if(lab) {
        text(1.5, yt + 12, bquote(RH[acute]==paste(.(acute.sc), ' for ', .(dur.ac), ' months')),pos = 4, col = col.ac)
    }}

hazLabs <- function(yt=30)  {   text(5,yt, 'acute', col = col.ac, pos=1)
                                text(60,yt, 'chronic', col = col.ch, pos=1)
                                text(103,yt, 'late', col = col.lt, pos=1)
                                text(118,yt, 'AIDS', col = col.aids, pos=1)}

## Contour/image plot settings
xs <- seq(-1,7.2, by = .05) ##log(seq(0.1,8, by = 1))
ys <- seq(-2.4,2.7, by = .05) ##log(seq(0.1, 11, by = .5))
zs <- outer(exp(xs)-1,exp(ys))
range(zs)
which(zs==min(zs),T)
zs[1,103]
ys[103]
levels <- c(-10,-2,0,2,10,20,50,100,200,500,1000,10000,10^5)
trfxn <- function(x) log(x+10)
utrfxn <- function(x) exp(x)-10
levels <- utrfxn(seq(trfxn(min(zs)),trfxn(max(zs)), l = 20))
levels[which.min(abs(levels))] <- 0
levels <- c(-10, levels[levels>=0])
levels <- round(levels, -1)
levels <- c(-10,0,10,20,40,70,130,250,500,1000,2000,4000,7000,13000,20000)
cols <- colorRampPalette(c('skyblue2','purple','orange','red'))(length(levels)-1)
## blevels <- c(-25,-10,-2,0)
## bcols <- colorRampPalette(c('dark green','purple'))(length(blevels)-1)
## cols <- c(bcols[-length(bcols)], cols)
## levels <- c(blevels[-length(blevels)],levels)
mar1 <- c(4,2,2.5,1)
mar2 <- c(4,5,2.5,0)
options("scipen"=999)
col.ac <- 'brown'
col.ch <- 'black'
col.lt <- gray(.3)
col.aids <- gray(.6)
xlab1 <- expression(RH[acute])
ylab1 <- expression(d[acute])
xlab2 <- 'relative infectivity of acute phase vs. chronic phase'
ylab2 <- 'acute phase duration in months'

####################################################################################################
## Figure 1. Show EHM diagram and RH[acute]/d[acute] collinearity plot with MCMC posteriors from refit of Holl Mod.

showpts <- FALSE
for(ff in 1:3) {
    if(ff==1) pdf(file.path(ABCoutdir,'Figure 1 - EHM diagram and RH_ac d_ac collinearity.pdf'), w = 6.83, h = 5)
    if(ff==2) png(file.path(ABCoutdir,'Figure 1 - EHM diagram and RH_ac d_ac collinearity.png'), w = 6.83, h = 5, units='in',res=200)
    if(ff==3) tiff(file.path(ABCoutdir,'Figure 1 - EHM diagram and RH_ac d_ac collinearity.tiff'), w = 6.83, h = 5, units='in',res=300)
    par('ps'=12)
    layout(t(matrix(c(1:2,4,1,3,4),3,2)), w = c(.8,1,.27))
    ## ################################################
    ## conceptual diagram
    par(mar=mar1)
    plot(0,0, type = 'n', xlim = c(-10,125), ylim = c(0, 600), bty = 'n', axes=F, xlab='', ylab='')
    title(main='A', adj = 0)
    title(xlab='years since infection', mgp = c(2.5,0,0))
    title(ylab='hazard profile', mgp=c(.5,1,0))
    axis(1, at = seq(0,120, by = 12), 0:10)
    ## hazard profiles
    drawHazProf(600, 16, 5,1)
    drawHazProf(510, 26, 3,2)
    drawHazProf(370, 101, .75,3)
    hazLabs(360)
  #  drawHazProf(400, 42, 1.5,3)
    rhmn <- signif(exp(abcCIs['50%','acute.sc']),2)
    dmn <- signif(exp(abcCIs['50%','dur.ac']),2)
    drawHazProf(100, rhmn, dmn,4)
    hazLabs(90)

    ## ################################################
    ## contour plot of EHM from Hollingsworth et al. variable hazard survival model refit with MCMC
    par(mar=mar2)
    image(xs, ys, zs, breaks = levels, xlim = c(-1,6.2), ylim = c(-2.4,2.6), mgp = c(3,0,0),
          col = cols, axes = F, xlab='',ylab=ylab2)
    if(showpts) with(fout$posts, points(acute.sc, dur.ac, pch=16, cex = .4, col = gray(.6)))
    title(main='B', adj = 0)
    title(xlab=xlab2, mgp=c(2,0,0))
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
    lines(HollCont[[1]]$x, HollCont[[1]]$y)
    text(log(38), log(5.8),  '95% credible contour', pos = 4, srt=-45)
    points(log(c(101,26,16)),log(c(.75,3,5)), pch=as.character(3:1))
    points(log(c(101,26,16)),log(c(.75,3,5)), pch = 21, cex = 2.5)
    text(log(outtab['50%','acute.sc']), log(.14), '95% CrI', pos = 3)
    text(log(4), log(outtab['50%','dur.ac']),  '95% CrI', pos = 2, srt=90)
    with(data.frame(outtab), arrows(log(acute.sc[c(1)]), log(.15), log(acute.sc[c(3)]), log(.15), len = .05, angle = 90, code = 3))
    with(data.frame(outtab), arrows(log(4), log(dur.ac[c(1)]), log(4), log(dur.ac[c(3)]), len = .05, angle = 90, code = 3))

    ## contour plot of EHM from Bellan et al. couple transmission model fit with ABC SMC
    par(mar=mar2)
    image(xs, ys, zs, breaks = levels, xlim = c(-1,6.2), ylim = c(-2.4,2.6), mgp = c(3,0,0),
          col = cols, axes = F, xlab='', ylab=ylab2)
    if(showpts) with(pmatChosen, points(logacute.sc, logdur.ac, pch=16, cex = .4, col = gray(.6)))
    title(main='C', adj = 0)
    title(xlab=xlab2, mgp=c(2,0,0))
    points(abcCIs['50%','acute.sc'], abcCIs['50%','dur.ac'], pch = as.character(4))
    points(abcCIs['50%','acute.sc'], abcCIs['50%','dur.ac'], pch = 21, cex = 2.5)
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
    lines(abcCont[[1]]$x, abcCont[[1]]$y) #
    arrows(abcCIs['2.5%','acute.sc'], log(.2), abcCIs['97.5%','acute.sc'], log(.2), len = .05, angle = 90, code = 3)
    arrows(log(.45), abcCIs['2.5%','dur.ac'], log(.45), abcCIs['97.5%','dur.ac'], len = .05, angle = 90, code = 3)
    ## Palette legend
    par(mar=rep(0,4))
    plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
    legseq <- levels[levels<=1000] # c(0,10,20,40,70,130,240,1000)
    show <- levels<=1000 #& levels > -5
    sbcolor.legend(.09,.1,.15,.8, legend=levels[show], legseq=legseq, rect.col = cols[show], gradient = "y", cex = .8, browse=F)
    text(.025, .85, expression(paste(EHM['acute'])), pos = 3, cex = 1.2)
    dev.off()
}

####################################################################################################

####################################################################################################
## Same for late phase
####################################################################################################
fout$exposts <- ddply(fout$exposts, .(), transform, atr.month.ltaids = (late.sc-1)*(dur.lt) -dur.aids)
head(fout$exposts)
ltp <- c('late.sc','dur.lt')
cont.lt <- HPDregionplot(fmcmc, vars = ltp, prob = c(.95), lims = c(-1.5,10,-2,11), n = 80, h = c(.5,.1),
                      bty = 'n', xlab = expression(RH['late']), ylab = 'late phase duration (months)')
meds.lt <- apply(fout$exposts[,ltp], 2, median)
apply(fout$exposts[,ltp], 2, function(x) quantile(x, c(.025,.975)))
adp <- c('late.sc','dur.aids')
cont.ad <- HPDregionplot(fmcmc, vars = adp, prob = c(.95), lims = c(-1.5,10,-2,11), n = 80, h = c(.5,.1),
                      bty = 'n', xlab = expression(RH['late']), ylab = 'aids phase duration (months)')
meds.ad <- apply(fout$exposts[,adp], 2, median)
apply(fout$exposts[,adp], 2, function(x) quantile(x, c(.025,.975)))
## ehm.late plot
pdf(file.path(outdir,'late contour (LOG).pdf'), w = 6.83, h = 3)
ct <- .7
layout(matrix(c(1:3),1,3), w = c(1,.3,1))
## RHlate vs dur.lt
par(cex.lab = ct, cex.axis = ct, cex.main = ct)
par(mar=c(3,3,.3,0))
xs <- seq(-1,7.2, by = .05) ##log(seq(0.1,8, by = 1))
ys <- seq(-2.4,6, by = .05) ##log(seq(0.1, 11, by = .5))
levels <- c(0,2,5,10,25,50,70,100,200,500,1000,10000,10^5)
cols <- colorRampPalette(c('purple','orange','red'))(length(levels)-1)
image(xs, ys, outer(exp(xs)-1,exp(ys)), breaks = levels, xlim = c(-1,5.7), ylim = c(-2.4,5), mgp = c(2,0,0),
               col = cols, axes = F, xlab = expression(paste(RH['late'])), ylab = 'late phase duration (months)')
xts <- c(1:9, seq(10, 90, by = 10), seq(100, 300, by = 100))
xsh <- c(1,10,100,1000)
xls <- xts
xls[!xls %in% xsh] <- NA
axis(1, at = log(xts), lab = xls)
yts <- c(seq(.1,.9, by = .1),1:9, seq(10, 100, by = 10))
ysh <- c(.1,1,10,100)
yls <- yts
yls[!yls %in% ysh] <- NA
axis(2, at = log(yts), lab = yls, las = 2)
lines(cont.lt[[1]]$x, cont.lt[[1]]$y)
## points(fout$posts$late.sc, fout$posts$dur.lt, cex = .3)
points(log(outtab[4,'late.sc']), log(outtab[4,'dur.lt']), pch = 15)#, col = 'points')
## text(log(outtab[4,'late.sc']), log(outtab[4,'dur.lt']), 'Hollingsworth', pos=3, cex=.5)
points(log(meds.lt['late.sc']), log(meds.lt['dur.lt']), pch = 19)#, col = 'points')
## text(log(meds.lt['late.sc']), log(meds.lt['dur.lt']), 'median', pos=1, cex=.5)
with(data.frame(outtab), axis(1, at = log(late.sc[c(1,3)]), lab = NA, line = -1.5))
with(data.frame(outtab), axis(2, at = log(dur.lt[c(1,3)]), lab = NA, line = -1))
## Palette legend
par(mar=rep(0,4), cex.lab = ct, cex.axis = ct, cex.main = ct)
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
color.legend(.09,.1,.15,.8, levels, rect.col = cols, gradient = "y", cex = ct)
text(.025, .8, expression(paste(EHM['late'])), pos = 3, cex = ct)
####################################################################################################
## RHlate vs dur.aids
par(cex.lab = ct, cex.axis = ct, cex.main = ct)
par(mar=c(3,4,.3,1))
xs <- seq(-1,7.2, by = .05) ##log(seq(0.1,8, by = 1))
ys <- seq(-2.4,6, by = .05) ##log(seq(0.1, 11, by = .5))
levels <- c(0,2,5,10,25,50,70,100,200,500,1000,10000,10^5)
cols <- rep(NA, length(levels)-1) #colorRampPalette(c('purple','orange','red'))(length(levels)-1)
## filled.contour(xs, ys, outer(exp(xs),exp(ys)), level = levels, xlim = c(.5,7.2), ylim = c(-2,3),
##                col = cols, plot.axes=F)
image(xs, ys, outer(exp(xs),exp(ys)), breaks = levels, xlim = c(-1,5.7), ylim = c(-2.4,5), mgp = c(2,0,0),
               col = cols, axes = F, xlab = expression(paste(RH['late'])), ylab = 'AIDS phase duration (months)')
xts <- c(1:9, seq(10, 90, by = 10), seq(100, 200, by = 100))
xsh <- c(1,10,100,1000)
xls <- xts
xls[!xls %in% xsh] <- NA
axis(1, at = log(xts), lab = xls)
yts <- c(seq(.1,.9, by = .1),1:9, seq(10, 100, by = 10))
ysh <- c(.1,1,10,100)
yls <- yts
yls[!yls %in% ysh] <- NA
axis(2, at = log(yts), lab = yls, las = 2)
lines(cont.ad[[1]]$x, cont.ad[[1]]$y)
## points(fout$posts$late.sc, fout$posts$dur.ad, cex = .3)
points(log(outtab[4,'late.sc']), log(outtab[4,'dur.aids']), pch = 15)#, col = 'points')
## text(log(outtab[4,'late.sc']), log(outtab[4,'dur.aids']), 'Hollingsworth', pos=3, cex=.5)
points(log(meds.ad['late.sc']), log(meds.ad['dur.aids']), pch = 19)#, col = 'points')
## text(log(meds.ad['late.sc']), log(meds.ad['dur.aids']), 'median', pos=1, cex=.5)
with(data.frame(outtab), axis(1, at = log(late.sc[c(1,3)]), lab = NA, line = -1.5))
with(data.frame(outtab), axis(2, at = log(dur.aids[c(1,3)]), lab = NA, line = -1))
dev.off()
####################################################################################################

####################################################################################################
## Same for late+AIDS phase
####################################################################################################
fout$exposts <- ddply(fout$exposts, .(), transform, atr.month.ltaids = (late.sc-1)*(dur.lt) -dur.aids, dur.ltaids = dur.lt + dur.aids,
                      ltaids.sc = (late.sc*dur.lt + 0*dur.aids)/(dur.lt+dur.aids))
outtab <- ddply(outtab, .(), transform, atr.month.ltaids = (late.sc-1)*(dur.lt) -dur.aids, dur.ltaids = dur.lt + dur.aids,
                      ltaids.sc = (late.sc*dur.lt + 0*dur.aids)/(dur.lt+dur.aids))
head(fout$exposts)
posts <- fout$exposts
posts[,c('ltaids.sc','dur.ltaids','atr.month.ltaids')] <- log(posts[,c('ltaids.sc','dur.ltaids','atr.month.ltaids')])
ltp <- c('ltaids.sc','dur.ltaids')
cont.ltaids <- HPDregionplot(posts, vars = ltp, prob = c(.95), lims = c(-1.5,10,-2,11), n = 80, h = c(.5,.1))
meds.ltaids <- apply(posts[,ltp], 2, median)
cis <- apply(exp(posts[,ltp]), 2, function(x) quantile(x, c(.025,.975)))
apply(exp(posts[,ltp]), 2, range)
quantile(fout$exposts$atr.month.ltaids, c(.025, .5, .975))
## EHM.ltaids plot
pdf(file.path(outdir,'EHM lateaids (LOG).pdf'), w = 3.27, h = 3)
ct <- .6
layout(matrix(c(1:2),1,2), w = c(1,.25))
## (RHlate*dur.lt + RHaids*dur.aids)/(dur.lt + dur.aids)  vs (dur.lt+dur.aids)
par(cex.lab = ct, cex.axis = ct, cex.main = ct)
par(mar=c(3,3,.3,0))
xs <- seq(-1,3, by = .05) ##log(seq(0.1,8, by = 1))
ys <- seq(0,4, by = .05) ##log(seq(0.1, 11, by = .5))
pal <- colorRamp(c('purple','orange','red'))
levels <- c(-25,-10,-5,-2,0,2,5,10,25,50,70,100,200,500,1000)
rg <- range(levels)
cols <- colorRampPalette(c('purple','orange','red'))(length(levels)-1)
## filled.contour(xs, ys, outer(exp(xs),exp(ys)), level = levels, xlim = c(.5,7.2), ylim = c(-2,3),
##                col = cols, plot.axes=F)
zs <- outer(exp(xs)-1,exp(ys))
image(xs, ys, zs, breaks = levels, xlim = c(-.4,2.5), ylim = c(0,4), mgp = c(2,0,0),
               col = cols, axes = F, xlab = expression(paste((RH['late']*d['late']+RH['AIDS']*d['AIDS'])/(d['late']+d['AIDS']))),
      ylab = expression(paste((d['late']+d['AIDS']))))
xts <- c(1:10)
xsh <- c(1,10,100,1000)
xls <- xts
xls[!xls %in% xsh] <- NA
axis(1, at = log(xts), lab = xls)
yts <- c(1:9, seq(10, 50, by = 10))
ysh <- c(.1,1,10,100)
yls <- yts
yls[!yls %in% ysh] <- NA
axis(2, at = log(yts), lab = yls, las = 2)
lines(cont.ltaids[[1]]$x, cont.ltaids[[1]]$y)
#with(posts, points(ltaids.sc, dur.ltaids, cex = .08, pch = 19))
points(log(outtab[4,'ltaids.sc']), log(outtab[4,'dur.ltaids']), pch = 15, cex = .8, col = gray(.2))#, col = 'points')
## text(log(outtab[4,'late.sc']), log(outtab[4,'dur.lt']), 'Hollingsworth', pos=3, cex=.5)
points(meds.ltaids['ltaids.sc'], meds.ltaids['dur.ltaids'], pch = 19, cex = .8, col=gray(.3))#, col = 'points')
## text(log(meds.lt['late.sc']), log(meds.lt['dur.lt']), 'median', pos=1, cex=.5)
with(data.frame(cis), axis(1, at = log(ltaids.sc[c(1,2)]), lab = NA, line = -1.5))
with(data.frame(cis), axis(2, at = log(dur.ltaids[c(1,2)]), lab = NA, line = -1))
## Palette legend
par(mar=rep(0,4), cex.lab = ct, cex.axis = ct, cex.main = ct)
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
ticks <- levels#[-length(levels)]#ticks[-length(ticks)]
#cols <- apply(pal(trf(ticks)), 1, function(x) rgb(x[1],x[2],x[3], maxColorValue=255))
color.legend(.09,.1,.15,.8, ticks, rect.col = cols, gradient = "y", cex = ct)
text(.025, .8, expression(paste(EHM['late'+'AIDS'])), pos = 3, cex = ct)
graphics.off()



## pdf(file.path(outdir,'Figure SX - fitted RHacute dacute.pdf'), w = 4, h = 4)
## with(pmatChosen, plot(acute.sc, dur.ac, cex = .3, log='xy', xlim = c(.5,200), ylim = c(.5,10), las = 3))
## graphics.off()
