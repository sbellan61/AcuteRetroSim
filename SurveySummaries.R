######################################################################
#  Create a data frame with one row per country (dframe) or per
##  survey (dframe.s) with accompanying summary characteristics.
######################################################################
rm(list=ls())                           # clear workspace
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
#load("data files/ds.name.Rdata")        # country names
load("data files/ds.nm.all.Rdata") # country names
load("data files/allDHSAIS.Rdata")         # DHS data
#load("data files/alldhs.Rdata")         # DHS data
load("data files/epic.Rdata")       # Infectious HIV prevalence in men
head(dat,2)
outdir <- file.path('results','PrevFigs')
if(!file.exists(outdir)) dir.create(outdir)
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")

####################################################################################################
## initialize dframe: one row per country: name; serodiscordant proportion (SDP; psdc), male SDP, female SDP
dframe <- data.frame(matrix(NA,length(ds.nm), 14))
colnames(dframe) <- c('country','psdc','pmsdc','pfsdc','curprev','peak','peak.nart','tpeak','tpeak.nart',
                      'pms','pfs','rms','rfs', 'col')
dframe$col <- rainbow(nrow(dframe))
dframe$country <- levels(dat$group)
####################################################################################################
## plot epidemic curves for each country, get epidemic peaks (of prevalence & infectious
## prevalence--the latter assumes the proportion infected but on ART are not infectious). Also fill
## in dframe as we go.
pdf(file.path(outdir,'epics.pdf'))
plot(0,0,type='n',xlim = c(1975,2015), ylim = c(0,.5), bty='n', xlab='', ylab = 'female prevalence')
cols <- rainbow(10)
ord <- c(which(ds.nm=='WA'), which(ds.nm!='WA')) # order so WA first so plot looks best
for(cc in ord) {## for each country
    temp <- dat[dat$group==dframe$country[cc],] ## select that country's data
    dframe$psdc[cc] <- sum(temp$ser %in% 2:3) / sum(temp$ser %in% 1:3) ## SDP
    dframe$pmsdc[cc] <- sum(temp$ser %in% 2) / sum(temp$ser %in% 1:3) ## M SDP
    dframe$pfsdc[cc] <- sum(temp$ser %in% 3) / sum(temp$ser %in% 1:3) ## F SDP
    dframe$curprev[cc] <- mean(prev.inf[temp$tint,temp$epic.ind])        ## infectious prevalence at mean interview time
    ## add epidemic info (female curves)
    if(dframe$country[cc] !='WA') { # non-WA countries
        lines((1:nrow(prev.inf))/12+1900, prev.inf[,temp$epic.ind[1]], col = cols[cc], lty = 2) # infectious prevalence
        lines((1:nrow(prev.inf))/12+1900, prev.all[,temp$epic.ind[1]], col = cols[cc], lty = 1) # all prevalence
        ## peak height        
        dframe$peak[cc] <- max(prev.inf[,temp$epic.ind[1]])
        dframe$peak.nart[cc] <- max(prev.all[,temp$epic.ind[1]])        
        ## months between epidemic peak & mean(survey time)
        dframe$tpeak[cc] <- mean(temp$tint) - which.max(prev.inf[,temp$epic.ind[1]]) #infectious
        dframe$tpeak.nart[cc] <- mean(temp$tint) - which.max(prev.all[,temp$epic.ind[1]]) #all
        points(which.max(prev.inf[,temp$epic.ind[1]])/12 + 1900, dframe$peak[cc], col = cols[cc], pch = 21) # add points
        points(which.max(prev.all[,temp$epic.ind[1]])/12 + 1900, dframe$peak.nart[cc], col = cols[cc], pch = 19)        
      }else{ ## get weighted average of epidemic peak for W Africa
        maxes <- apply(prev.inf[,unique(temp$epic.ind)], 2, max) # max infectious prevalence for countries
        tpeaks <- apply(prev.inf[,unique(temp$epic.ind)], 2, which.max) # peak infectious prevalence timing for all
        maxes.nart <- apply(prev.all[,unique(temp$epic.ind)], 2, max) # ditto for all prevalence
        tpeaks.nart <- apply(prev.all[,unique(temp$epic.ind)], 2, which.max)
        weights <- table(temp$epic.nm)/nrow(temp) # weight these by # of couples from each country
        dframe$peak[cc] <- sum(maxes*weights)
        dframe$tpeak[cc] <- mean(temp$tint) - sum(tpeaks*weights)
        dframe$peak.nart[cc] <- sum(maxes.nart*weights)
        dframe$tpeak.nart[cc] <- mean(temp$tint) - sum(tpeaks.nart*weights)
        for(cc.wa in unique(temp$epic.ind)) { ## lines for each country
            lines((1:nrow(prev.inf))/12+1900, prev.inf[,cc.wa], col = cols[cc])
            lines((1:nrow(prev.inf))/12+1900, prev.all[,cc.wa], col = cols[cc], lty = 2)            
          }
        points(sum(tpeaks*weights)/12+1900, dframe$peak[cc], col = cols[cc], pch = 19) # points for all countries
        points(sum(tpeaks.nart*weights)/12+1900, dframe$peak.nat[cc], col = cols[cc], pch = 21)        
      }
    ## now add demography summaries to the data
    dframe$pms[cc] <- mean(temp$tmar - temp$tms)/12    ##  premarital duration of sexual activity
    dframe$pfs[cc] <- mean(temp$tmar - temp$tfs)/12 
    dframe$rms[cc] <- mean(temp$tint - temp$tmar)/12    ##  marital duration of sexual activity
    dframe$rfs[cc] <- mean(temp$tint - temp$tmar)/12 
  }
dframe$tpeak <- dframe$tpeak/12 ## transform to years
dframe$tpeak.nart <- dframe$tpeak.nart/12
dframe$psdc.m <- dframe$pmsdc / (dframe$pmsdc + dframe$pfsdc) # proportion of SDC couples that are M+
legend('topleft', leg = dframe$country, col = cols, bty = 'n', pch = 19) 
dev.off()
## log-ratio of time spent in couple to time spent sexually active before
dframe$logmrelpre <- log(dframe$rms / dframe$pms) 
dframe$logfrelpre <- log(dframe$rfs / dframe$pfs)
save(dframe, file=file.path("data files/dframe.Rdata"))

####################################################################################################
## initialize dframe.s: one row per survey
dframe.s <- data.frame(matrix(NA,length(unique(dat$ds)), 19))
colnames(dframe.s) <- c('ds', 'country','group','yr',
                        'psdc','pmsdc','pfsdc',
                        'psdc.f','pmsdc.f','pfsdc.f',
                        'curprev','peak','peak.nart','tpeak','tpeak.nart',
                      'pms','pfs','rms','rfs')
dframe.s$ds <- unique(dat$ds)
dframe.s <- dframe.s[order(dframe.s$ds),]
for(cc in 1:nrow(dframe.s)) {
    temp <- dat[dat$ds==dframe.s$ds[cc],]
    dframe.s$group[cc] <- temp$group[1]
    bfmar <- temp$m.fun == 1 & temp$f.fun == 1 ## both partners in first union?
    dframe.s$psdc[cc] <- sum(temp$ser %in% 2:3) / sum(temp$ser %in% 1:3) # SDP
    dframe.s$pmsdc[cc] <- sum(temp$ser %in% 2) / sum(temp$ser %in% 1:3)  # M SDP
    dframe.s$pfsdc[cc] <- sum(temp$ser %in% 3) / sum(temp$ser %in% 1:3)  # F SDP
    ## ditto above but first marriage only
    dframe.s$psdc.f[cc] <- sum(temp$ser[bfmar] %in% 2:3) / sum(temp$ser[bfmar] %in% 1:3) # SDP  
    dframe.s$pmsdc.f[cc] <- sum(temp$ser[bfmar] %in% 2) / sum(temp$ser[bfmar] %in% 1:3)  # M SDP
    dframe.s$pfsdc.f[cc] <- sum(temp$ser[bfmar] %in% 3) / sum(temp$ser[bfmar] %in% 1:3)  # F SDP
    dframe.s$yr[cc] <- mean(temp$tint)/12+1900 
    dframe.s$country[cc] <- temp$epic.nm[1]
    ## add hiv prevalence at mean survey date
    dframe.s$prev[cc] <- prev.inf[mean(temp$tint),temp$epic.ind[1]]
  }
dframe.s
 
####################################################################################################
## Plot SDP by survey year for all couples & for those only with both in their first union (to see
## if there are SDP changes over time).
cols <- rainbow(length(unique(dframe.s$country)))
dframe.s$col <- NA
pdf(file.path(outdir, 'sdp by survey (after filter).pdf'), w = 8, h = 8)
plot(0,0, type = 'n', xlim = c(2000,2013), ylim = c(0,1), xlab = 'year', ylab = 'SDP')
for(cc in 1:length(unique(dframe.s$country))) {
    tcount <- unique(dframe.s$country)[cc] # country name
    dframe.s$col[dframe.s$country==tcount] <- cols[cc] # color
    points(dframe.s$yr[dframe.s$country == tcount], dframe.s$psdc[dframe.s$country == tcount], # all
           type = 'b', pch = 21, col = dframe.s$col[dframe.s$country == tcount])
    points(dframe.s$yr[dframe.s$country == tcount], dframe.s$psdc.f[dframe.s$country == tcount], # both in first union
           type = 'b', pch = 19, col = dframe.s$col[dframe.s$country == tcount])
  }
legend('topleft', unique(dframe.s$country), col = cols, pch = 19, bty = 'n')
legend('bottomright', c('all','first marriage'), pch = c(21,19), bty = 'n')
dev.off()
save(dframe.s, file=file.path("data files","dframe.s.Rdata"))
###################################################################### 

####################################################################################################
## Now make same types of summary plots but with all the data before pre-processing to remove
## couples with missing/inconsistent data for model fitting.
load("data files/allRawDHSAIS.Rdata")       # raw data
##################################################
## First need to fix Ethiopian survey dates due to their calendar shift.
## Ethiopia 2005: "The survey was fielded from April 27 to August 30, 2005." (p. 11, Ethiopia DHS
## 2005 report) difference in dates is then
diff2005 <- 105*12 + 4 - min(allraw[allraw$ds=="Ethiopia 2005","tint"])
allraw[allraw$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] <-
    allraw[allraw$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] + diff2005
## Ethiopia 2011:"All data collection took place over a five-month period from 27 December 2010 to 3
## June 2011." (p. 10, Ethiopia DHS 2011 report)
diff2011 <- 110*12 + 12 - min(allraw[allraw$ds=="Ethiopia 2011","tint"])
allraw[allraw$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] <-
    allraw[allraw$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] + diff2011
##################################################
##  get discordance rates by countries for raw data (draw)
draw <- data.frame(matrix(NA,length(ds.nm), 13))
colnames(draw) <- c('country','psdc','pmsdc','pfsdc','curprev','peak','peak.nart','tpeak','tpeak.nart',
                      'pms','pfs','rms','rfs')
draw$country <- unique(dat$group)
draw <- draw[order(draw$country),]
head(draw)
for(cc in 1:nrow(draw)) {
    temp <- allraw[allraw$group==draw$country[cc],]
    draw$psdc[cc] <- sum(temp$ser %in% 2:3, na.rm=T) / sum(temp$ser %in% 1:3, na.rm=T) # SDP  
    draw$pmsdc[cc] <- sum(temp$ser %in% 2, na.rm=T) / sum(temp$ser %in% 1:3, na.rm=T)  # M SDP
    draw$pfsdc[cc] <- sum(temp$ser %in% 3, na.rm=T) / sum(temp$ser %in% 1:3, na.rm=T)  # F SDP
    ##  premarital duration of sexual activity
    draw$pms[cc] <- mean(temp$tmar - temp$tms, na.rm=T)/12
    draw$pfs[cc] <- mean(temp$tmar - temp$tfs, na.rm=T)/12 
    ##  marital duration of sexual activity
    draw$rms[cc] <- mean(temp$tint - temp$tmar, na.rm=T)/12
    draw$rfs[cc] <- mean(temp$tint - temp$tmar, na.rm=T)/12 
  }
draw$psdc.m <- draw$pmsdc / (draw$pmsdc + draw$pfsdc)
draw <- draw[order(draw$country),]
head(draw,3)
save(draw, file="data files/draw.Rdata")

####################################################################################################
##  get discordance rates by surveys for raw data (draw.s), do for first mar only too.
draw.s <- data.frame(matrix(NA,length(unique(allraw$ds)), 17))
colnames(draw.s) <- c('country','ds', 'tint.yr',
                      'psdc','pmsdc','pfsdc',
                      'psdc.f','pmsdc.f','pfsdc.f',
                      'pms','pfs','rms','rfs',
                      'pms.f','pfs.f','rms.f','rfs.f')
draw.s$ds <- unique(allraw$ds)
draw.s <- draw.s[order(draw.s$ds),]
show <- c('Mnumber.of.unions', 'Fnumber.of.unions')
for(cc in 1:nrow(draw.s)) {
    temp <- allraw[allraw$ds==draw.s$ds[cc],]
    draw.s$country[cc] <- temp$epic.nm[1]
    draw.s$psdc[cc] <- sum(temp$ser %in% 2:3, na.rm=T) / sum(temp$ser %in% 1:3, na.rm=T) # SDP  
    draw.s$pmsdc[cc] <- sum(temp$ser %in% 2, na.rm=T) / sum(temp$ser %in% 1:3, na.rm=T)  # M SDP
    draw.s$pfsdc[cc] <- sum(temp$ser %in% 3, na.rm=T) / sum(temp$ser %in% 1:3, na.rm=T)  # F SDP
    ## for both partners in first marriage only
    bfmar <- temp$Mnumber.of.unions == 'Once' & temp$Fnumber.of.unions == 'Once' 
    draw.s$psdc.f[cc] <- sum(temp$ser[bfmar] %in% 2:3, na.rm=T) / sum(temp$ser[bfmar] %in% 1:3, na.rm=T) # SDP  
    draw.s$pmsdc.f[cc] <- sum(temp$ser[bfmar] %in% 2, na.rm=T) / sum(temp$ser[bfmar] %in% 1:3, na.rm=T)  # M SDP
    draw.s$pfsdc.f[cc] <- sum(temp$ser[bfmar] %in% 3, na.rm=T) / sum(temp$ser[bfmar] %in% 1:3, na.rm=T)  # F SDP
    ##  premarital duration of sexual activity
    draw.s$pms[cc] <- mean(temp$tmar - temp$tms, na.rm=T)/12
    draw.s$pfs[cc] <- mean(temp$tmar - temp$tfs, na.rm=T)/12 
    ##  marital duration of sexual activity
    draw.s$rms[cc] <- mean(temp$tint - temp$tmar, na.rm=T)/12
    draw.s$rfs[cc] <- mean(temp$tint - temp$tmar, na.rm=T)/12
    ##  premarital duration of sexual activity (first union)
    draw.s$pms.f[cc] <- mean(temp$tmar[bfmar] - temp$tms[bfmar], na.rm=T)/12
    draw.s$pfs.f[cc] <- mean(temp$tmar[bfmar] - temp$tfs[bfmar], na.rm=T)/12 
    ##  marital duration of sexual activity (first union)
    draw.s$rms.f[cc] <- mean(temp$tint[bfmar] - temp$tmar[bfmar], na.rm=T)/12
    draw.s$rfs.f[cc] <- mean(temp$tint[bfmar] - temp$tmar[bfmar], na.rm=T)/12
    ## mean interview date
    draw.s$tint.yr[cc] <- mean(temp$tint, na.rm = T) / 12 + 1900
    ## add hiv prevalence at mean survey date
    draw.s$prev[cc] <- prev.all[mean(temp$tint),temp$epic.ind[1]]
  }
draw.s$psdc.m <- draw.s$pmsdc / (draw.s$pmsdc + draw.s$pfsdc)
## add country colors to data frame
cols <- rainbow(length(unique(draw.s$country)))
draw.s$col <- NA
for(cc in 1:length(unique(draw.s$country))) {
    tcount <- unique(draw.s$country)[cc] # country name
    draw.s$col[draw.s$country==tcount] <- cols[cc] # color
}
head(draw.s,5)
save(draw.s, file="data files/draw.s.Rdata")


####################################################################################################
## Plot SDP by survey year for all couples & for those only with both in their first union (to see
## if there are SDP changes over time). *RAW DATA*
cols <- rainbow(length(unique(draw.s$country)))
draw.s$col <- NA
pdf(file.path(outdir,'sdp by survey (unfiltered data).pdf'), w = 8, h = 8)
plot(0,0, type = 'n', xlim = c(2000,2013), ylim = c(0,1), xlab = 'year', ylab = 'SDP')
for(cc in 1:length(unique(draw.s$country)))  {
    tcount <- unique(draw.s$country)[cc] # temp country
    draw.s$col[draw.s$country==tcount] <- cols[cc] # assign color
    points(draw.s$tint.yr[draw.s$country == tcount], draw.s$psdc[draw.s$country == tcount],
           type = 'b', pch = 21, col = draw.s$col[draw.s$country == tcount])
    points(draw.s$tint.yr[draw.s$country == tcount], draw.s$psdc.f[draw.s$country == tcount],
           type = 'b', pch = 19, col = draw.s$col[draw.s$country == tcount])
  }
legend('topleft', unique(draw.s$country), col = cols, pch = 19, bty = 'n')
legend('bottomright', c('all','first marriage'), pch = c(21,19), bty = 'n')
dev.off()

####################################################################################################
## Miscellaneous serostatus plots for presentations with raw data
####################################################################################################
for(col in c('black','white')) {
    presdir <- file.path('results', 'PresentationFigures')
    if(!file.exists(presdir))       dir.create(presdir)
###################################################################### 
    pdf(file.path(presdir,paste0("serostatus breakdown ",col,".pdf")), width = 5.5, height = 4)
    par(mar=c(5,6,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
    tab1 <- xtabs(~group + ser, allraw)
    tab1 <- tab1[nrow(tab1):1,]
    total <- rowSums(tab1)
    total <- as.matrix(total)[,rep(1,4)]
    tab1 <- tab1/total
    tab1
    rownames(tab1)[rownames(tab1)=="WA"] <- "West Africa"
    cols <- c("yellow","green","purple","dark gray")
    barplot(t(tab1), beside = F, names.arg = rownames(tab1), horiz = T, las = 2,
            col = cols, xlab = "proportion of couples in serogroup", border = NA)
    dev.off()
###################################################################### 
    pdf(file.path(presdir,paste0("serostatus breakdown leg ",col,".pdf")), width = 6, height = 2.5)
    par(mar=rep(0,4), fg = col, col.axis = col, col.main = col, col.lab = col)
    plot(0,0, type = "n", bty = "n", xlab = "", ylab = "", axes = F)
    legend("top", c("M+ F+", "M+ F-", "M- F+", "M- F-"),
           pch = 15, bty = "n", col = cols, cex = 1.5)
    dev.off()
###################################################################### 
    pdf(file.path(presdir,paste0("serostatus breakdown neg only ",col,".pdf")), width = 8, height = 5)
    par(mar=c(5,6,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
    tab1 <- xtabs(~group + ser, allraw)
    tab1 <- tab1[nrow(tab1):1,]
    total <- rowSums(tab1)
    total <- as.matrix(total)[,rep(1,4)]
    tab1 <- tab1/total
    tab1
    rownames(tab1)[rownames(tab1)=="WA"] <- "West Africa"
    cols <- c(NA,"green","purple","dark gray")
    barplot(t(tab1), beside = F, names.arg = rownames(tab1), horiz = T, las = 2,
            col = cols, xlab = "proportion of couples in serogroup", border = NA)
    dev.off()
###################################################################### 
###################################################################### 
    pdf(file.path(presdir,paste0("serostatus breakdown pos only ",col,".pdf")), width = 5.5, height = 4)
    par(mar=c(5,6,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
    tab1 <- xtabs(~group + ser, allraw)
    tab1 <- tab1[nrow(tab1):1,]
    total <- rowSums(tab1)
    total <- as.matrix(total)[,rep(1,4)]
    tab1 <- tab1/total
    tab1
    rownames(tab1)[rownames(tab1)=="WA"] <- "West Africa"
    cols <- c("yellow","green","purple",NA)
    barplot(t(tab1), beside = F, names.arg = rownames(tab1), horiz = T, las = 2,
            col = cols, xlab = "proportion of couples in serogroup", border = NA)
    dev.off()
###################################################################### 
###################################################################### 
    pdf(file.path(presdir,paste0("serostatus breakdown pos only perc ",col,".pdf")), width = 5.5, height = 4)
    par(mar=c(5,6,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
    tab1 <- xtabs(~group + ser, allraw)
    tab1 <- tab1[nrow(tab1):1,]
    total <- rowSums(tab1)
    total <- as.matrix(total)[,rep(1,4)]
    tab1 <- tab1/total
    tab1
    rownames(tab1)[rownames(tab1)=="WA"] <- "West Africa"
    cols <- c("yellow","green","purple",NA)
    bb <- barplot(t(tab1), beside = F, names.arg = rownames(tab1), horiz = T, las = 2,
                  col = cols, xlab = "proportion of couples in serogroup", border = NA)
    text(rev(1-tab1[,4]),rev(bb), paste(signif(dframe$psdc*100,3),'%',sep=''), pos = 4)
    dev.off()
###################################################################### 
###################################################################### 
    pdf(file.path(presdir,paste0("serostatus proportion ",col,".pdf")), width = 5.5, height = 4)
    par(mar=c(5,6,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
    tab2 <- xtabs(~group + ser, allraw)
    tab2 <- tab2[nrow(tab2):1,]
    total <- rowSums(tab2[,-4])
    total <- as.matrix(total)[,rep(1,3)]
    tab2 <- tab2[,-4]/total
    tab2 <- tab2[,3:1]
    rownames(tab2)[rownames(tab2)=="WA"] <- "West Africa"
    cols <- rev(c("yellow","green","purple"))
    bb <- barplot(t(tab2), beside = F, names.arg = rownames(tab2), horiz = T, las = 2,
                  col = cols, xlab = "serodiscordance proportion (SDP)", border = NA)
    dev.off()
###################################################################### 
###################################################################### 
    pdf(file.path(presdir,paste0("serostatus proportion only sdc ",col,".pdf")), width = 5.5, height = 4)
    par(mar=c(5,6,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
    tab2 <- xtabs(~group + ser, allraw)
    tab2 <- tab2[nrow(tab2):1,]
    total <- rowSums(tab2[,-4])
    total <- as.matrix(total)[,rep(1,3)]
    tab2 <- tab2[,-4]/total
    tab2 <- tab2[,3:1]
    rownames(tab2)[rownames(tab2)=="WA"] <- "West Africa"
    cols <- rev(c(NA,"green","purple"))
    bb <- barplot(t(tab2), beside = F, names.arg = rownames(tab2), horiz = T, las = 2,
                  col = cols, xlab = "serodiscordance proportion (SDP)", border = NA)
    dev.off()
###################################################################### 
    ## lognormal risk distributions
    for(sd in c(.5, 1, 2, 3))  {
        pdf(file.path(presdir,paste('log normal sd', sd, col, '.pdf')), w = 3, h = 3)
        par(mar=c(5,.5,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
        curve(dnorm(x, 0, sd), from = -8, 8, axes = F, ylab = '', xlab = expression(Z[i]), ylim = c(0,.8))
        axis(1, at = log(10^c(-3:3)), label = 10^c(-3:3), las = 2)
        dev.off()
    }
}


####################################################################################################
## prevalence vs SDP
col <- 'white'
cols <- rainbow(length(ds.nm))
pdf(file.path(presdir,paste('prevalence vs SDP', col, '.pdf')), w = 7, h = 4.4)
layout(matrix(c(1,2),nr=1,nc=2), w = c(1,.38))
par(mar=c(5,4,.5,.5), fg = col, col.axis = col, col.main = col, col.lab = col)
plot(dframe$peak, dframe$psdc, pch = 19, col = cols, xlim = c(0,.3), ylim = c(0,1), xlab = 'peak HIV prevalence', ylab = 'serodiscordant proportion', main ='',
     bty='n', axes = F)
mod <- lm(psdc~peak, dframe)
abline(mod)
fx <- function(x) { 2*x*(1-x) / ( 2*x*(1-x) + x^2 ) }
#curve(fx, from = 0, to = .3, col = 'yellow', add=T)
axis(1, at = seq(0,.3, b = .05))
axis(2, at = seq(0,1,b=.25), las = 2)
ord <- rev(order(dframe$psdc))
par(mar=c(.2,0,0,0))
plot.new()
legend('top', ds.nm[ord], col=cols[ord], title = 'top to bottom', cex = .8, inset = .1, pch = 19) #, bty = 'n')
graphics.off()
