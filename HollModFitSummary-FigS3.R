rm(list=ls(all=T)); gc()
library(plyr); library(data.table); library(abind); library(multicore)
####################################################################################################
## Summarize Hollingsworth Model fits to Hollingsworth simulated & Bellan simulated data.
## Create Figure S3.
####################################################################################################
outdir <- file.path('FiguresAndTables','HollModTestFitsSummary')
if(!file.exists(outdir)) dir.create(outdir)
load(file = file.path('results','RakAcute','HollModTestFits', 'blocksHollFitTest.Rdata'))
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') ##  transmission coefficient names, for convenience
nc <- 12                                       ## core per simulation
## Compare these to fits of SB cohort data, load those (from SummarizeFits.R output) for comparison
load(file='results/RakAcute/wf.Rdata')

## Load fit files from Hollingsworth-Model-generated data
fls <- list.files(file.path('results','RakAcute','HollModTestFits','fitouts'), pattern = 'fitout-', full.names=T)
length(fls)
jobnums <-  as.numeric(sapply(fls, function(x) as.numeric(strsplit(x,'-')[[1]][2])))
fls <- fls[order(jobnums)]
jobnums <- jobnums[order(jobnums)]

## Load & label files by jobnum (fitout1, fitout2, etc...)
system.time(alldat <- mclapply(fls, function(x) { jobnum <- as.numeric(strsplit(x,'-')[[1]][2]); load(x)
                                    assign(paste0('fitout',jobnum),fitout, env=.GlobalEnv)}))

## Each list element is one simulation output from a job
names(alldat[[1]])
alldat[[1]]$wtab ## raw simulation data

## Sample Sizes
wtabs <- t(abind(lapply(alldat, function(x) { c(x$wtab$inct[1,'n'],x$wtab$prevt[1,'n'],x$wtab$latet[nrow(x$wtab$latet),'n'])}), along = 2))
wtabs <- as.data.frame(wtabs)
colnames(wtabs) <- c('incn','prevn','laten')
wtabs$job <- jobnums
head(wtabs)
apply(wtabs[,-4],2, summary)

## Hollingsworth Fits to Holl-gen data
alldat[[1]]$outtab.h
hfits <- abind(lapply(alldat, function(x) { x$outtab.h}), along = 3)
dimnames(hfits)[[3]] <- jobnums
dimnames(hfits)[[1]] <- c('lci','med','uci','true')

## Variables use to create and fit the retrospective cohort
cohvars <- data.frame(rbindlist(lapply(alldat, function(x) {as.data.frame(t(x$dpars))})))
cohvars$job <- jobnums
dpars <- blocks[jobnums,]

## turn these into easy to work with data frames
hfh <- adply(hfits, 3:2, .parallel=T)
colnames(hfh)[1:2] <- c('job','var')
hfh$het.sd <- dpars$het.gen.sd[match(hfh$job, dpars$job)]
hfh <- merge(hfh, cohvars, by = 'job')
hfh <- hfh[order(hfh$var,hfh$job),]
head(hfh)
hfh <- within(hfh, {ehm.late <- (late.sc - 1)*dur.lt})
hfh <- within(hfh, {ehm.acute <- (acute.sc - 1)*dur.ac})

####################################################################################################
## Figure S3: Hollingsworth sensitivity analysis: showing estimates the elevated hazard months of the acute phase
unique(hf$ehm.late)
unique(hfh$ehm.late)
ehm.lates <- unique(c(hf$ehm.late,hfh$ehm.late))
ehm.lates <- ehm.lates[order(ehm.lates)] ## which ehmls exist
ehmls <- c(0,40,90) ## which EHMlates to show (for clearer plots)
xmax <- 60
ymax <- 60
var <- 'ehm.ac'
cols <- c('purple','red','orange')
pdf(file.path(outdir,'Figure S3- HollMod of SB & HollMod generated data true vs est ehmacute.pdf'), w = 3.27, h = 3)
par(mfrow=c(1,1), mar = c(3.5,3.5,1,.5), 'ps'=10, mgp = c(2.2,1,0))
plot(0,0, type='n', xlim = c(0, xmax), ylim = c(0,ymax), bty = 'n', xlab = expression(paste('true ',EHM[acute])),
     ylab = expression(paste('estimated ',EHM[acute])), main = '', asp=1)
clip(0,xmax,0,ymax)
abline(a=0, b=1, col = gray(.4), lty = 1, lwd=2)
for(jj in 1:2) { ## for Holl,SB generated data
  dd <- c('hf','hfh')[jj]
  temp <- get(dd)
  ##temp <- temp[temp$ehm.late %in% ehmls,]
  temp <- temp[rev(order(temp$ehm.late)),]
  for(ee in length(ehmls):1) {
      if(dd=='hfh') sel <- temp$var==var & temp$ehm.late==ehmls[ee]
      if(dd=='hf') sel <- temp$var==var & temp$het.sd==0 & temp$err=='base' & temp$ehm.late==ehmls[ee]
      ## with(temp[sel,],  arrows(true, lci, true, uci, length=.05, angle = 90, code = 3,col = cols[ee]), lwd = .7)
      ## with(temp[sel,],  points(true, med, col = cols[ee], cex = .5, pch = 19))
      clip(0,xmax,0,ymax)
      lines(1:100, with(temp[sel,], predict(loess(med ~ true), 1:100)), col = cols[ee], lwd = 1, lty = jj)
  }
  ## segments(1, 75.4, xmax, 75.4, col = 'black', lty = 3, lwd = 2)
  ## text(30, 72.4, 'Hollingsworth Estimate', pos = 3, cex = ct)
}
clip(0,xmax*1.5,0,ymax)
legend(0,ymax, title= expression(paste(EHM[late])), ncol=1,
       leg = ehmls, lty = 1, col=cols, bty ='n', cex = .7)
legend('bottomright', leg = c('Hollingsworth', 'Bellan'), title='Data-Generation Model', lty = 2:1, cex = .7, bty = 'n')
dev.off()
