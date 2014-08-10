library(mnormt);library(coda); library(stats4);library(plyr)#library(lme4)
rm(list=ls(all=TRUE))
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions
load(file.path('results','RakAcute','blocks.Rdata')) ## these are country-acute phase specific blocks

fls <- list.files("results/RakAcute/UgandaFits/fitouts", full.names=T)
fls.sh <- list.files("results/RakAcute/UgandaFits/fitouts")
jobnums <- as.numeric(sapply(fls.sh, function(x) strsplit(x, '-')[[1]][2]))

wtab.base.ls <- list(NA)
wtab.XbErr.ls <- list(NA)

for(ii in 1:length(jobnums)) {
    if(ii %% 10 == 0)   print(paste('on', ii, 'of', length(jobnums)))
    load(fls[ii])
    wtab.base.ls[[ii]] <- sbmod.to.wdat(fitout$rakll, excl.by.err = F, browse=F)
    wtab.XbErr.ls[[ii]] <- sbmod.to.wdat(fitout$rakll, excl.by.err = T, browse=F)
}

names(wtab.base.ls) <- jobnums
names(wtab.XbErr.ls) <- jobnums

inc.sdp.base <- as.numeric(lapply(wtab.base.ls, function(x) x$inct[1,"n"]))
inc.sdp.XbErr <- as.numeric(lapply(wtab.XbErr.ls, function(x) x$inct[1,"n"]))
head(inc.sdp.base)
head(inc.sdp.XbErr)

cpls.excl <- inc.sdp.base - inc.sdp.XbErr
cpls.excl.prop <- (inc.sdp.base - inc.sdp.XbErr)/ inc.sdp.base

cdat <- data.frame(job = jobnums, cpls = cpls.excl, cpls.prop = cpls.excl.prop)
cdat$acute.sc <- blocks$acute.sc[match(cdat$job, blocks$jobnum)]
cdat$dur.ac <- blocks$dur.ac[match(cdat$job, blocks$jobnum)]
cdat$het.gen.sd <- blocks$het.gen.sd[match(cdat$job, blocks$jobnum)]
cdat$ehm.ac <- with(cdat, (acute.sc-1)*dur.ac)

## Only plot for plausible parameters (since this is what we should be comparing to the real data censorship rates)
plaus.par <- read.csv("results/RakAcute/UgandaFitSummaries/RHacute & dacute combos with 95CI.csv")
cdat$plaus <- with(cdat, paste0(acute.sc,'-',dur.ac)) %in% with(plaus.par, paste0(acute.sc, '-', dur.ac))

## Save it!
save(cdat, wtab.base.ls, wtab.XbErr.ls, file = file.path('results/RakAcute/UgandaFitSummaries/','wtabs.Rdata'))

head(cdat)

pdf(file.path("results/RakAcute/UgandaFitSummaries/", "Incident SD couples excluded distribution (ALL).pdf"), w = 3.27, h = 3)
par('ps'=8, mar = c(4,4,.5,.5))
hsds <- 0:3
cols <- c("black", 'red','purple','orange')
plot(0, 0, bty='n', xlab = expression(EHM[acute]), ylab = 'proportion seroincident couples excluded',
         xlim = c(0, 100), ylim = c(.2, .7), type = 'n', las = 1)
for(ii in 4:1) {
    hh <- hsds[ii]
    with(cdat[cdat$het==hh,], points(ehm.ac, cpls.prop, col = cols[ii], pch = 21, cex = .7))
}
legend('topright', leg = 0:3, title = expression(sigma[hazard]), pch = 21, col = cols, ncol=4)
dev.off()

pdf(file.path("results/RakAcute/UgandaFitSummaries/", "Incident SD couples excluded distribution (plausible parameters).pdf"),
    w = 3.27, h = 3)
par('ps'=8, mar = c(4,4,.5,.5))
hsds <- 0:3
plot(0, 0, bty='n', xlab = expression(EHM[acute]), ylab = 'proportion seroincident couples excluded',
         xlim = c(0, 50), ylim = c(.4, .6), type = 'n', las = 1)
for(ii in 2) {
    hh <- hsds[ii]
    with(cdat[cdat$plaus & cdat$het==hh,], points(ehm.ac, cpls.prop, col = cols[ii], pch = 21, cex = .7))
}
dev.off()



## hist(cpls.excl, xlab = "# couples excluded by \n following Rakai exclusion criteria", ylab = "frequency", las = 1, col = 'black',
##      xlim = c(0,1000), main = '')

