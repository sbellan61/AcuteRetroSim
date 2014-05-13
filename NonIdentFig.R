library(plyr); library(data.table); library(abind); library(multicore); library(emdbook);library(coda); library(plotrix)
rm(list=ls(all=T)); gc()
## Show that duration & infectivity of acute phase are not identifiable
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','HollingsworthAn','Figures')
if(!file.exists(outdir)) dir.create(outdir)
load(file = file.path('results','HollingsworthAn','RealExclbyErr','workspace.Rdata'))


fexp <- (fout$exposts)
fexp$bp.ac <- fexp$bp * fexp$acute.sc
quantile(fexp$bp.ac, c(0.025, .5, .975))*12
outtab <- outtab[,,'XbErr']

fmcmc <- as.mcmc(fout$posts)
meds <- apply(fout$exposts[,c('acute.sc','dur.ac')], 2, median)
apply(fout$exposts[,c('acute.sc','dur.ac')], 2, function(x) quantile(x, c(.025,.975)))

test <- HPDregionplot(fmcmc, vars = c('acute.sc','dur.ac'), prob = c(.95), lims = c(-1.5,10,-2,11), n = 28, h = c(1,.08),
                      bty = 'n', xlab = expression(RH['acute']), ylab = 'acute phase duration (months)')

pdf(file.path(outdir,'acute contour.pdf'))
xs <- seq(0.1,150, by = 1)
ys <- seq(0.1, 11, by = .5)
filled.contour(log(xs), log(ys), outer(xs-1,ys), nlevels=10, xlim = c(0,200), ylim = c(-.3, 11))
lines(exp(test[[1]]$x), exp(test[[1]]$y))
points(meds['acute.sc'], meds['dur.ac'], pch = 19)
points(fout$exposts$acute.sc, fout$exposts$dur.ac, cex = .3)
dev.off()


pdf(file.path(outdir,'acute contour (LOG).pdf'), w = 3.27, h = 3)
ct <- .7
par(cex.lab = ct, cex.axis = ct, cex.main = ct)
layout(matrix(c(1,2),1,2), w = c(1,.3))
par(mar=c(3,3,.3,0))
xs <- seq(-.3,7.2, by = .05) ##log(seq(0.1,8, by = 1))
ys <- seq(-2.4,2.7, by = .05) ##log(seq(0.1, 11, by = .5))
levels <- c(0,2,5,10,25,50,70,100,200,500,1000,10000,10^5)
cols <- colorRampPalette(c('purple','orange','red'))(length(levels)-1)
image(xs, ys, outer(exp(xs)-1,exp(ys)), breaks = levels, xlim = c(-.3,7.2), ylim = c(-2.4,2.7), mgp = c(2,0,0),
               col = cols, axes = F, xlab = expression(paste(RH['acute'])), ylab = 'acute phase duration (months)')
xts <- c(1:9, seq(10, 90, by = 10), seq(100, 1000, by = 100))
xsh <- c(1,5,10,50,100,500)
xls <- xts
xls[!xls %in% xsh] <- NA
axis(1, at = log(xts), lab = xls)
yts <- c(seq(.1,.9, by = .1),1:10)
ysh <- c(.1,.5,1,5,10)
yls <- yts
yls[!yls %in% ysh] <- NA
axis(2, at = log(yts), lab = yls, las = 2)
lines(test[[1]]$x, test[[1]]$y)
#points(fout$posts$acute.sc, fout$posts$dur.ac, cex = .3)
points(log(outtab[4,'acute.sc']), log(outtab[4,'dur.ac']), pch = 15)#, col = 'points')
text(log(outtab[4,'acute.sc']), log(outtab[4,'dur.ac']), 'Hollingsworth', pos=3, cex=.5)
points(log(meds['acute.sc']), log(meds['dur.ac']), pch = 19)#, col = 'points')
text(log(meds['acute.sc']), log(meds['dur.ac']), 'median', pos=1, cex=.5)
with(data.frame(outtab), axis(1, at = log(acute.sc[c(1,3)]), lab = NA, line = -1.5))
with(data.frame(outtab), axis(2, at = log(dur.ac[c(1,3)]), lab = NA, line = -2)) 
## Palette legend
par(mar=rep(0,4), cex.lab = ct, cex.axis = ct, cex.main = ct)
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
color.legend(.09,.1,.15,.8, levels, rect.col = cols, gradient = "y", cex = ct)
text(.025, .8, expression(paste(EHM['acute'])), pos = 3, cex = ct)
dev.off()

####################################################################################################
## Late phase
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

####################################################################################################
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
## Late+AIDS phase
fout$exposts <- ddply(fout$exposts, .(), transform, atr.month.ltaids = (late.sc-1)*(dur.lt) -dur.aids, dur.ltaids = dur.lt + dur.aids,
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

####################################################################################################
## EHM.ltaids plot
pdf(file.path(outdir,'EHM lateaids JD (LOG).pdf'), w = 3.27, h = 3)
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
               col = cols, axes = F, xlab = expression(paste((RH['L']*d['L']+RH['AIDS']*d['AIDS'])/(d['L']+d['AIDS']))),
      ylab = expression(paste((d['L']+d['AIDS']))))
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
text(.025, .8, expression(paste(EHM['L'+'AIDS'])), pos = 3, cex = ct)
graphics.off()

####################################################################################################
## With conceptual plot too
pdf(file.path(outdir,'acute contour & EHM diagram (LOG).pdf'), w = 6.83, h = 2.7)
ct <- 1
par(cex.lab = ct, cex.axis = ct, cex.main = ct)
layout(matrix(c(1:3),1,3), w = c(.8,1,.3))
##################################################
## conceptual diagram
par(mar= c(4,.5,1,1))
plot(0,0, type = 'n', xlim = c(-5,125), ylim = c(0, 250), bty = 'n', axes=F, xlab='years since infection',ylab='')
axis(1, at = seq(0,120, by = 12), 0:10)
mtext('relative hazard (vs chronic phase)', side = 2, line = -.5, adj = .5, cex = .7)
## example 1
##axis(1, at = c(0,.5,110,120), pos = 0)
yt <- 0
polygon(c(0,.75,.75,0), yt + c(0,0,101,101), col = gray(.4), border=NA) ## acute
polygon(c(100,110,110,100), yt + c(0,0,7,7),col = gray(.6), border=NA) ## late
polygon(c(0,110,110,0), yt + c(0,0,1,1), col = gray(.3), border=NA) ## chronic
segments(110,yt,120,yt, lty = 3)
## example 2
##axis(1, at = c(0,2,110,120), pos = yt)
yt <- 170
polygon(c(0,3,3,0), yt +c(0,0,26,26), col = gray(.4), border=NA) ## acute
polygon(c(100,110,110,100), yt +c(0,0,7,7),col = gray(.6), border=NA) ## late
polygon(c(0,110,110,0), yt +c(0,0,1,1), col = gray(.3), border=NA) ## chronic
segments(110,yt,120,yt, lty = 3)
## example 3
##axis(1, at = c(0,2,110,120), pos = yt)
yt <- 215
polygon(c(0,5,5,0), yt + c(0,0,16,16), col = gray(.4), border=NA) ## acute
polygon(c(100,110,110,100), yt + c(0,0,7,7),col = gray(.6), border=NA) ## late
polygon(c(0,110,110,0), yt + c(0,0,1,1), col = gray(.3), border=NA) ## chronic
segments(110,yt,120,yt, lty = 3)
##################################################
## contour plot of EHM
par(mar=c(4,4,.3,0))
xs <- seq(-.3,7.2, by = .05) ##log(seq(0.1,8, by = 1))
ys <- seq(-2.4,2.7, by = .05) ##log(seq(0.1, 11, by = .5))
levels <- c(0,2,5,10,25,50,70,100,200,500,1000,10000,10^5)
cols <- colorRampPalette(c('purple','orange','red'))(length(levels)-1)
image(xs, ys, outer(exp(xs)-1,exp(ys)), breaks = levels, xlim = c(-.3,6.2), ylim = c(-2.4,2.35), mgp = c(3,0,0),
               col = cols, axes = F, xlab = expression(paste(RH['acute'])), ylab = 'acute phase duration (months)')
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
lines(test[[1]]$x, test[[1]]$y)
#points(fout$posts$acute.sc, fout$posts$dur.ac, cex = .3)
points(log(outtab[4,'acute.sc']), log(outtab[4,'dur.ac']), pch = 15)#, col = 'points')
## text(log(outtab[4,'acute.sc']), log(outtab[4,'dur.ac']), 'Hollingsworth', pos=3, cex=ct)
points(log(meds['acute.sc']), log(meds['dur.ac']), pch = 19)#, col = 'points')
## text(log(meds['acute.sc']), log(meds['dur.ac']), 'median', pos=1, cex=ct)
## with(data.frame(outtab), axis(1, at = log(acute.sc[c(1,3)]), lab = NA, line = -1.5))
## with(data.frame(outtab), axis(2, at = log(dur.ac[c(1,3)]), lab = NA, line = -2)) 
with(data.frame(outtab), arrows(log(acute.sc[c(1)]), log(.2), log(acute.sc[c(3)]), log(.2), len = .05, angle = 90, code = 3))
with(data.frame(outtab), arrows(log(5), log(dur.ac[c(1)]), log(5), log(dur.ac[c(3)]), len = .05, angle = 90, code = 3))
## Palette legend
par(mar=rep(0,4), cex.lab = ct, cex.axis = ct, cex.main = ct)
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
color.legend(.09,.1,.15,.8, levels, rect.col = cols, gradient = "y", cex = ct*.8)
text(.025, .8, expression(paste(EHM['acute'])), pos = 3, cex = ct*1.2)
dev.off()
