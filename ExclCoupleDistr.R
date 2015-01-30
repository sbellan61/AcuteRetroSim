library(mnormt);library(coda); library(stats4);library(plyr)#library(lme4)
rm(list=ls(all=TRUE))
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions
source('abcFunctions.R')                   # 
load(file.path('results','RakAcute','blocks.Rdata')) ## these are country-acute phase specific blocks
outdir <- file.path('FiguresAndTables/UgandaFitSummaries')

## Take this from ABC-SMC fits
batch <- 5
load(file=file.path(out.dir, paste0('IntermedDistr',batch,'.Rdata'))) ## Load last distribution (already filtered)
pmatChosen$numExcluded <- numExcluded(pmatChosen$propExcl)

exclSummary <- signif(rbind(
    100*c(mean = mean(pmatChosen$propExcl), quantile(pmatChosen$propExcl, c(.025,.5,.975), na.rm=T)),
    c(mean = mean(pmatChosen$numExcluded), quantile(pmatChosen$numExcluded, c(.025,.5,.975), na.rm=T))),3)
rownames(exclSummary) <- c('proportion excluded', '# excluded')
write.csv(exclSummary, file.path(outdir, 'excluded Summary.csv'))

pdf(file.path(outdir, "Incident SD couples excluded distribution (ABC-SMC fit).pdf"),
    w = 3.27, h = 5)
par('ps'=12, mar = c(4,4,1.5,.5), mfrow = c(2,1))
plot(1, 1, bty='n', xlab = expression(EHM[acute]), ylab = 'proportion excluded',
     xlim = c(-10, 100), ylim = c(0, .7), type = 'n', las = 1)#, log='x')
    with(as.data.frame(pmatChosen), points(EHMacute, propExcl, pch = 16, cex = .7))
title('A', adj = 0)
abline(h=.47, col = 'red')
    with(as.data.frame(pmatChosen), hist(numExcluded, breaks = 0:max(numExcluded+1), col='black',main='', xlab = '# excluded', ylab = 'frequency'))
title('B', adj = 0)
dev.off()

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
plaus.par <- read.csv(file.path(outdir, 'RHacute & dacute combos with 95CI.csv'))
cdat$plaus <- with(cdat, paste0(acute.sc,'-',dur.ac)) %in% with(plaus.par, paste0(acute.sc, '-', dur.ac))

## Save it!
save(cdat, blocks, wtab.base.ls, wtab.XbErr.ls, file = file.path(outdir,'wtabs.Rdata'))

head(cdat)

pdf(file.path(outdir, "Incident SD couples excluded distribution (ALL).pdf"), w = 3.27, h = 3)
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

pdf(file.path(outdir, "Incident SD couples excluded distribution (plausible parameters).pdf"),
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

## Look at distribution of # of couples
load(file.path(outdir,'wtabs.Rdata'))
plausParms <- read.csv(file = file.path(outdir,'RHacute & dacute combos with 95CI.csv'))

infprob <- function(x) data.frame(inct = x$inct$i[1:4] / x$inct$n[1:4],  prevt = x$prevt$i[1:4] / x$prevt$n[1:4], latet = x$latet$i[1:4] / x$latet$n[1:4])

lsSelect <- names(wtab.XbErr.ls) %in% plausParms$job
names(wtab.XbErr.ls)[lsSelect] ## jobs that match

lsSelect <- 1:length(wtab.XbErr.ls)

wtabXbErrPlaus <- wtab.XbErr.ls[lsSelect]
wtabBasePlaus <- wtab.base.ls[lsSelect]

wtabXbErrPlaus[1][[1]]$inct
wtabBasePlaus[1][[1]]$inct

tail(fls)

infProbs95 <- mclapply(wtabXbErrPlaus, infprob)
IncProbs95 <- t(as.data.frame(lapply(infProbs95, '[', 'inct')))
PrevProbs95 <- t(as.data.frame(lapply(infProbs95, '[', 'prevt')))
rownames(IncProbs95) <- rownames(PrevProbs95)<- names(wtabXbErrPlaus)

IncProbs95
PrevProbs95

apply(IncProbs95, 2, range)
apply(PrevProbs95, 2, range)
apply(PrevProbs95, 2, mean)
infprob(wtab.rl)

-log(1-apply(PrevProbs95, 2, mean)) ## rate


load(fls[jobnums==3058])
fitout$dpars
fitout$het.gen.sd
names(fitout)
plausParms[plausParms$job==3058,]

head(twf)
with(twf, twf[job==3058 & hobs=='obsNA' & err == 'XbErr',])

with(twf, twf[acute.sc == 5 & dur.ac == 3 & hobs=='obsNA' & err == 'XbErr',])
## 2798 is het = 0, 3058 is hete = 1
