####################################################################################################
## Explore plasuibility of hypothesis that SDP variation is driven by differences in couple
## dissolution rates between countries. If this was the case we'd expect SDP (particularly in
## couples with partners who are both in their first stable couple) to be greatest in countries with
## the shortest relationship durations, because ++ couples would be condissipating quickly and +-
## couples when dissolving are mmore likely to form another +- couple than a ++ couple since
## prevalence < 50%.
####################################################################################################
rm(list=ls())                           # clear workspace
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
source('SimulationFunctions.R')
load('data files/ds.nm.all.Rdata')        # country names
#load("data files/allDHSAIS.Rdata")         # DHS data
load("data files/allDHSAISdis.Rdata")      # raw dhsais with extra mfun & ffun variables
load('data files/dframe.Rdata') # country summary data (peak prevalence, country-prevalence, etc...)
load('data files/dframe.s.Rdata')# country summary data (peak prevalence, country-prevalence, etc...) by DHS survey
load('data files/draw.s.Rdata')# country summary data (peak prevalence, country-prevalence, etc...) by DHS survey (unfiltered data)
load('data files/pars.arr.ac.Rdata') # fitted transmission coefficients & other parameters for each assumed acute phase
load("data files/epic.Rdata")       # infectious HIV prevalence in women
load("data files/csurv.Rdata")    #  probability of survival (row) months by age (in months) at seroconversion (column)
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")   ## transmission coefficient names
outdir <- file.path('results','DissolutionFigs') ## create output directory
if(!file.exists(outdir))        dir.create(outdir)
dat <- dat.dis
odat <- dat.dis
dat$bd <- apply(cbind(dat$tmar-dat$tms,dat$tmar-dat$tfs), 1, max) # before couple duration of sexual activity (max across M/F)
dat$cd <- dat$tint - dat$tmar                                     # couple duration

####################################################################################################
## histogram of relationship durations for each country (1st unions only)
pdf(file.path(outdir,'Marital Duration Histograms (1st unions).pdf'))
brks <- seq(0,500/12, by = 1)
par(mfrow=c(4,4), mar = c(4,4,3,.5))
for(ii in 1:length(ds.nm)) {
    hist(dat$mardur.mon[dat$group==ds.nm[ii] & dat$mfun==1 & dat$ffun==1]/12, breaks = brks,
         xlim = c(0,500/12), main = ds.nm[ii], xlab='', col = 'black')
    abline(v=mean(dat$mardur.mon[dat$group==ds.nm[ii] & dat$mfun==1 & dat$ffun==1]/12, na.rm=T), col = 'red')
}
mtext('couple duration (yrs)', side = 1, outer = T, line = -2, adj = .5)
dev.off()

####################################################################################################
## Add mean, and 2.5, 97.5% quantile relationship duration by country-group to dframe both for (1)
## unweighted couples & (2) couples weighted by their inverse probability of survival to account for
## the fact that longer-lasting couples are more likely to have died of AIDS just from being older
## and around longer.
dframe$mreldur <- NA
for(ii in 1:nrow(dframe)) {
    print(paste0('appending relationship durations to dframe:', ds.nm[ii]))
    dframe$mreldur[ii] <- mean(dat$mardur.mon[dat$group==dframe$country[ii]])
    dframe$lreldur[ii] <- quantile(dat$mardur.mon[dat$group==dframe$country[ii]], .025)
    dframe$ureldur[ii] <- quantile(dat$mardur.mon[dat$group==dframe$country[ii]], .975)
    temp <- dat[dat$group==ds.nm[ii],]
    pars <- out.arr[hazs,2,which(in.arr[,1,2]==7), ii] # use survival inflator from fits with acute RH=7
    sim <- pcalc(pars, temp, sim = T, survive = T, browse = F, trace = F)
    ## infl fator = 1 / survival probabilities
    surp <- rowSums(sim$pser.a)
    infl <- 1/surp
    ## do weighted average of relationship durations to get survival-corrected values
    frel <- temp$mfun==1 & temp$ffun==1
    frel[is.na(frel)] <- F
    dframe$psdc.f[ii] <- sum(temp$ser[frel] %in% 2:3) / sum(temp$ser[frel] %in% 1:3)
    dframe$mreldur.s[ii] <- sum(temp$mardur.mon * infl) / (sum(infl)) # weighted
    dframe$mreldur.f[ii] <- mean(temp$mardur.mon[frel]) # first union
    dframe$mreldur.fs[ii] <- sum((temp$mardur.mon * infl)[frel]) / sum(infl[frel]) # weighted first union
  }

####################################################################################################
## Plot SDP vs mean relationship duration for survival inflated estimates & for all couples fitted
## in model & for all couples fitted in which both partners were in their first union)
dframe$col <- rainbow(nrow(dframe))
rdtab <- data.frame(matrix(NA,nr=4,nc=4))
colnames(rdtab) <- c('survival \ninflated', 'couples', 'effect (CI)', 'P')
rdtab[,1] <- rep(c('yes','no'),2)
rdtab[,2] <- rep(c('all','first \nunions'),each=2)
pdf(file.path(outdir, 'SDP vs Relationship Duration.pdf'), w = 7, h = 3.5) # initialize PDF
par(mar=c(4,4,1.5,0), mfrow = c(1,2))                            # set up margins & panels
ys <- dframe$psdc
mains <- c('all unions', 'first unions only')
## mains <- LETTERS ## if lettering-naming plots
ylabs <- c("serodiscordant proportion",'')
for(fu in 1:2) { ## for all unions & those in first union only
    if(fu==1) {
        xs <- dframe$mreldur/12         # all
        xs.s <- dframe$mreldur.s/12     # all weighted
    }else{
        xs <- dframe$mreldur.f/12       # first
        xs.s <- dframe$mreldur.fs/12    # first weighted
    }
    plot(0,0, type = 'n', xlim = c(8,13), ylim = c(0,1), bty = 'n', # initialize plot
         xlab = '', ylab=ylabs[fu], yaxt = 'n', main = mains[fu])
    axis(2, at = seq(0,1,l=5), las = 2)     # add y-axis
    points(xs, ys, pch = 19, col = dframe$col) # SDP vs un-weighted mean relationship duration
    points(xs.s, ys, pch = 21, col = dframe$col) # SDP vs weighted mean relationship duration
    ## arrows(dframe$lreldur/12, ys, dframe$ureldur/12, ys) ## add CI's   
    lmod <- lm(ys ~ xs) ## build linear model for both weighted & unweighted
    lmod.s <- lm(ys ~ xs.s)
    rdtab[1+2*(fu-1),3] <- paste0(signif(summary(lmod)[[4]][2,1],3), ' (', # adding effect (CI) to table
                          signif(confint(lmod)[2,][1],3), ', ', signif(confint(lmod)[2,][2],3), ')')
    rdtab[2+2*(fu-1),3] <- paste0(signif(summary(lmod.s)[[4]][2,1],3), ' (',
                            signif(confint(lmod.s)[2,][1],3), ', ', signif(confint(lmod.s)[2,][2],3), ')')
    rdtab[1+2*(fu-1),4] <- signif(summary(lmod)[[4]][2,4],3) # adding P value to table
    rdtab[2+2*(fu-1),4] <- signif(summary(lmod.s)[[4]][2,4],3)
    lfun <- function(x) {coef(lmod)[1] + coef(lmod)[2]*x} # make model into y(x) functions for flexible plotting
    lfun.s <- function(x) {coef(lmod.s)[1] + coef(lmod.s)[2]*x}
    curve(lfun(x), from = 8, to = 13, add = T) # add curves to plot
    curve(lfun.s(x), from = 8, to = 13, add = T, lty = 2)
    leg.ord <- rev(order(ys))      # add legend (order as SDP)
}
mtext("mean relationship duration (years)", outer = T, side = 1, line = -1.5, adj = .6)
dev.off()
rdtab ## linear model table, save to file
write.csv(rdtab, file = file.path(outdir, 'SDP vs relationship duration LMs.csv'))
## Make separate file with a legend to put in the plot where wanted with Illustrator
pdf(file.path(outdir,'SDP vs reldur legend.pdf'), w = 3, h = 3)
par(mar=c(0,0,0,0))
plot(0,0, type = 'n', axes = F, xlab = '', ylab='', bty='n')
legend('topleft', ds.nm[leg.ord], pch = 19, col = dframe$col[leg.ord], ncol = 1, cex = .6)
legend('bottom', c('DHS data', 'AIDS-survival inflated DHS data'), pch = c(19,21), lty = 1:2, cex = .6)
dev.off()
####################################################################################################


####################################################################################################
## Plot SDP vs mean relationship duration for survival inflated estimates for just couples in their
## first union
pdf(file.path(outdir, 'SDP vs Relationship Duration (first unions only).pdf'), w = 4, h = 3.5) # initialize PDF
par(mar=c(4,4,1.5,0)) # set up margins & panels
ys <- dframe$psdc
mains <- c('all unions', 'first unions only')
## mains <- LETTERS ## if lettering-naming plots
ylabs <- c("serodiscordant proportion",'')
xs <- dframe$mreldur.f/12       # first
xs.s <- dframe$mreldur.fs/12    # first weighted
    plot(0,0, type = 'n', xlim = c(8,13), ylim = c(0,1), bty = 'n', # initialize plot
         xlab = '', ylab=ylabs[fu], yaxt = 'n', main = mains[fu])
    axis(2, at = seq(0,1,l=5), las = 2)     # add y-axis
    points(xs, ys, pch = 19, col = dframe$col) # SDP vs un-weighted mean relationship duration
    points(xs.s, ys, pch = 21, col = dframe$col) # SDP vs weighted mean relationship duration
    ## arrows(dframe$lreldur/12, ys, dframe$ureldur/12, ys) ## add CI's   
    lmod <- lm(ys ~ xs) ## build linear model for both weighted & unweighted
    lmod.s <- lm(ys ~ xs.s)
    lfun <- function(x) {coef(lmod)[1] + coef(lmod)[2]*x} # make model into y(x) functions for flexible plotting
    lfun.s <- function(x) {coef(lmod.s)[1] + coef(lmod.s)[2]*x}
    curve(lfun(x), from = 8, to = 13, add = T) # add curves to plot
    curve(lfun.s(x), from = 8, to = 13, add = T, lty = 2)
    leg.ord <- rev(order(ys))      # add legend (order as SDP)
mtext("mean relationship duration (years)", outer = T, side = 1, line = -1.5, adj = .6)
dev.off()
####################################################################################################

####################################################################################################
## Now explore how these figures would change if we only look at coules with a relationship duration
## shorter/longer than a specified cutoff (ct) and do a sensitivity analysis to the cutoff.
dnames <- c('dfa.sh','dfa.ln') ## shorter/longer than array names
cts <- 1:30
for(ss in 1:2) {                        ## loop over shorter/longer
    for(ct in 1:length(cts)) {
        cutoff <- cts[ct]
        if(ss==1) dat <- odat[odat$mardur.mon/12 <= cutoff,] else dat <- odat[odat$mardur.mon/12 >= cutoff,]
        temp.dfrs <- data.frame(ds = unique(dat$ds), prev = NA, group = NA, 
                                psdc = NA, pmsdc = NA, pfsdc = NA,
                                psdc.f = NA, pmsdc.f = NA, pfsdc.f = NA,
                                yr = NA)
        for(cc in 1:nrow(temp.dfrs)) {              # order so WA first so plot looks best
            temp <- dat[dat$ds==temp.dfrs$ds[cc],]
            temp.dfrs$group[cc] <- temp$group[1]
            bfmar <- temp$mfun == 1 & temp$ffun == 1
            temp.dfrs$psdc[cc] <- sum(temp$ser %in% 2:3) / sum(temp$ser %in% 1:3)
            temp.dfrs$pmsdc[cc] <- sum(temp$ser %in% 2) / sum(temp$ser %in% 1:3)
            temp.dfrs$pfsdc[cc] <- sum(temp$ser %in% 3) / sum(temp$ser %in% 1:3)
            ## first marriage only
            temp.dfrs$psdc.f[cc] <- sum(temp$ser[bfmar] %in% 2:3) / sum(temp$ser[bfmar] %in% 1:3)
            temp.dfrs$pmsdc.f[cc] <- sum(temp$ser[bfmar] %in% 2) / sum(temp$ser[bfmar] %in% 1:3)
            temp.dfrs$pfsdc.f[cc] <- sum(temp$ser[bfmar] %in% 3) / sum(temp$ser[bfmar] %in% 1:3)
            temp.dfrs$yr[cc] <- mean(temp$tint)/12+1900 
            ## temp.dfrs$country[cc] <- temp$epic.nm[1]
            ## add hiv prevalence at survey date
            temp.dfrs$prev[cc] <- epicf[mean(temp$tint),temp$epic.ind[1]]
        }
        if(ct==1) tempout <- temp.dfrs else tempout <- abind(tempout, temp.dfrs, along = 3)
    }
    assign(dnames[ss], tempout)
}
###################################################################################################
## plot results
pdf(file.path(outdir,'SDP by reldur (after duration cutting filter).pdf'), w = 6, h = 4)
layout(matrix(1:3,1,3), w = c(1,1,.3), h = rep(1,3))
par(mar = c(5,4,.5,.5))
xlabs <- c('relationships shorter than (yrs)','relationships longer than (yrs)')
for(ss in 1:2) {
    tempout <- get(dnames[ss])
    plot(0,0, type = 'n', xlim = c(0,max(cts)), ylim = c(0,1), xlab = xlabs[ss], ylab = 'SDP', xaxt = 'n', bty='n')
    axis(1, at = pretty(cts,4))
    for(cc in 1:nrow(tempout)) lines(cts,tempout[cc,'psdc',], type = 'l', col = dframe.s$col[cc])
    }
par(mar=rep(0,4))
plot(0,0,type = 'n', axes = F, bty ='n', xlab ='',ylab='')
sel <- !duplicated(dframe.s$country)
legend('topleft', dframe.s$country[sel], col = dframe.s$col[sel], pch = 19, bty = 'n', ncol = 1, cex = .7)
dev.off()
w
