setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','UgandaFitSummaries')

pdf(file.path(outdir,'EHM figure.pdf'), w = 3, h = 3)
par(mar= rep(0,4))
plot(0,0, type = 'n', xlim = c(-5,125), ylim = c(0, 250), bty = 'n', axes=F)
## example 1
##axis(1, at = c(0,.5,110,120), pos = 0)
yt <- 0
polygon(c(0,.5,.5,0), yt + c(0,0,151,151), col = gray(.7), border=NA) ## acute
polygon(c(100,110,110,100), yt + c(0,0,7,7),col = gray(.9), border=NA) ## late
polygon(c(0,110,110,0), yt + c(0,0,1,1), col = gray(.3), border=NA) ## chronic
segments(110,yt,120,yt, lty = 3)
## example 2
##axis(1, at = c(0,2,110,120), pos = yt)
yt <- 170
polygon(c(0,3,3,0), yt +c(0,0,26,26), col = gray(.7), border=NA) ## acute
polygon(c(100,110,110,100), yt +c(0,0,7,7),col = gray(.9), border=NA) ## late
polygon(c(0,110,110,0), yt +c(0,0,1,1), col = gray(.3), border=NA) ## chronic
segments(110,yt,120,yt, lty = 3)
## example 3
##axis(1, at = c(0,2,110,120), pos = yt)
yt <- 215
polygon(c(0,5,5,0), yt + c(0,0,16,16), col = gray(.7), border=NA) ## acute
polygon(c(100,110,110,100), yt + c(0,0,7,7),col = gray(.9), border=NA) ## late
polygon(c(0,110,110,0), yt + c(0,0,1,1), col = gray(.3), border=NA) ## chronic
segments(110,yt,120,yt, lty = 3)
dev.off()
 

