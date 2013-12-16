####################################################################################################
## Plot HIV prevalence & SDP vs fit transmission coefficients, and fit contact mixing coefficients.
rm(list=ls())                           # clear workspace
library(metatest);library(coda);library(faraway); library(hier.part); library(AICcmodavg); library(abind)
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
load("data files/allDHSAIS.Rdata")         # DHS data
load("data files/ds.nm.all.Rdata") # country names
load('data files/dframe.Rdata') # country summary data (peak prevalence, country-prevalence, etc...)
load("data files/draw.Rdata") # raw country summary data (peak prevalence, country-prevalence, etc...)
load('data files/dframe.s.Rdata')# country summary data (peak prevalence, country-prevalence, etc...) by DHS survey
load('data files/draw.s.Rdata')# country summary data (peak prevalence, country-prevalence, etc...) by DHS survey (unfiltered data)
load('data files/pars.arr.ac.Rdata') # fitted transmission coefficients & other parameters for each assumed acute phase
source('SimulationFunctions.R')

## add couple prevalence to dframe
for(gg in 1:nrow(dframe)) {             # for each country
    sel <- dat$group==dframe$country[gg] # select DHS data % of
    # coupled individuals that are HIV+ (1/discordant, 2/ ++ couple):
    dframe$cprev[gg] <- (sum(dat$ser[sel] %in% 2:3) + 2*sum(dat$ser[sel] == 1)) / (sum(sel)*2) 
  }


acutes <- as.numeric(in.arr[,1,2])
## make a list of data frame with betas, and contact coefficients (log-ed) for each acute fit
acs.to.do <- which(rowSums(!is.na(in.arr[,,2]))>0)
rdatl <- list()
for(aa in acs.to.do)  { ## for each acute fit
    pars.arr <- out.arr[,,aa,]          # select fitted parameters
    ## create data frame with couple prevalence, SDP, betas, & sexual mixing rates
    rdat <- data.frame(dframe$cprev, dframe$psdc, # DHS HIV prevalence in coupled indivs; proportion serodiscordant;
                       ## b=before (pre-); e=extra; p=partner (within-)
                       log(t(12*100*pars.arr['bp',c(2,1,3),])), log(t(12*100*pars.arr['be',c(2,1,3),])), # log(beta_within); log(beta_extra)
                       log(t(12*100*pars.arr['bmb',c(2,1,3),])), log(t(12*100*pars.arr['bfb',c(2,1,3),])), # log(beta male pre-); log(beta female pre-)
                       log(t(pars.arr['rr.ep',c(2,1,3),])), # log(extra- / within-) = log(extra-couple mixing coefficient) geometric mean across genders
                       ## below: log(pre- / within-) = log(pre-couple mixing coefficient) (by gender)
                       log(t(pars.arr['rr.bp.m',c(2,1,3),])), log(t(pars.arr['rr.bp.f',c(2,1,3),]))) 
    names(rdat) <- c('cprev','psdc',    # DHS HIV prevalence in coupled indivs; proportion serodiscordant;
                     'lbp','lbp.l','lbp.u', # log(beta_within) median & credible intervals (lower, upper)
                     'lbe','lbe.l','lbe.u', # ditto for log(beta_extra)
                     'lbmb','lbmb.l','lbmb.u', # log(male beta_pre)
                     'lbfb','lbfb.l','lbfb.u', # log(female beta_pre)
                     'lrr.ep','lrr.ep.l','lrr.ep.u', # log(extra-couple mixing coefficient)
                     'lrr.bp.m','lrr.bp.m.l','lrr.bp.m.u',# log(male pre-couple mixing coefficient)
                     'lrr.bp.f','lrr.bp.f.l','lrr.bp.f.u')# log(female pre-couple mixing coefficient)
    rdat$lgt.cprev <- logit(rdat$cprev) # logistic(DHS HIV prevalence in coupled individuals)
    rdat$lgt.pprev <- logit(dframe$peak.nart) # logistic(peak UNAIDS HIV prevalence in country)
    rdat$lgt.psdc <- logit(rdat$psdc) # logistic(DHS serodiscordant proportion)
    rdat$lgt.fD <- logit(rdat$psdc / (2*rdat$cprev * (1-rdat$cprev) / (2*rdat$cprev * (1-rdat$cprev) + rdat$cprev^2) )) # fancy D (normalized SDP)
    rdat$country <- dframe$country
    dat$fysa <- dat$tmar-dat$tfs
    dat$mysa <- dat$tmar-dat$tms
    rdat$fysa <- aggregate(dat$fysa, list(dat$group), mean)[,2]/12
    rdat$mysa <- aggregate(dat$mysa, list(dat$group), mean)[,2]/12
    rdat$mardur <- aggregate(dat$mardur.mon, list(dat$group), mean)[,2]/12
    rdat$lmpret <- rdat$lrr.bp.m + log(rdat$mysa) # pre * time
    rdat$lfpret <- rdat$lrr.bp.f + log(rdat$fysa) # pre * time
    rdat$lextrat <- rdat$lrr.ep + log(rdat$mardur) # extra * time
    rdat$linft <- rdat$lbp + log(rdat$mardur)      # infectivity * time
    rdat$lmysa <- log(rdat$mysa)
    rdat$lfysa <- log(rdat$fysa)
    rdat$lmardur <- log(rdat$mardur)
    rdat$lbp.w <- 1/((rdat$lbp.u - rdat$lbp.l)/(1.96*4))^2 # weights for each variable
    rdat$lbe.w <- 1/((rdat$lbe.u - rdat$lbe.l)/(1.96*4))^2
    rdat$lbmb.w <- 1/((rdat$lbmb.u - rdat$lbmb.l)/(1.96*4))^2
    rdat$lbfb.w <- 1/((rdat$lbfb.u - rdat$lbfb.l)/(1.96*4))^2
    rdat$geom.w <- (rdat$lbp.w * rdat$lbe.w * rdat$lbmb.w * rdat$lbfb.w)^.25 # geometric average of weights
    rdatl[[aa]] <- rdat               # add to array over acute fits
  }

wts <- rdatl[[3]][,c('country','lbp.w','lbe.w','lbmb.w','lbfb.w','geom.w')] # weights
wts
outdir <- file.path('results','PrevFigs')
sbpairs(wts[,-1], file.nm = file.path(outdir, 'wt correlations'))

## set plot parameters
las <- 1                                # axis text direction
leg.ord <- order(dframe$cprev)          # have legend order be by coupled prevalence so it matches top to bottom order of points
cols <- rainbow(length(ds.nm))[leg.ord]            # order colors similarly
ylab <- ''                              # no ylabels
r2 <- T                                 # show p values on regression plots
yaxt <- 'n'                             # no default axes
xaxt <- 'n'
leg.ord <- order(dframe$cprev)          # countries ordered by HIV coupled prevalence in legend
cols <- rainbow(length(ds.nm))[leg.ord]            # same for colors

## plot axes, logistic y, log x, for different limits & labels (pre-determined)
axis.fxn <- function(ylogit.axis, xlog.axis, ytype, xtype) {
  if(ylogit.axis) {
    if(ytype=='prev') {
      y.tcks <- c(seq(.01, .1, by = .01), seq(.1, .5, by = .1))
      y.labels <- y.tcks
      y.labels[!y.labels %in% c(.01,.1,.5)] <- NA
    }else{ ## serodiscordant proportion (ranges .4-.8)
      y.tcks <- c(seq(.3, .9, by = .1), seq(.91,.95, by = .01))
      y.labels <- y.tcks
      y.labels[!y.labels %in% c(.5,.9,.99)] <- NA
    }        
    axis(2, at = logit(y.tcks), label = y.labels, las = 2)
  }else{
    axis(2) ## not logistic
  }
  if(xlog.axis) {
    if(xtype=='within') { ## within-couple beta = HIV infectivity
      x.tcks <- c(1:9, seq(10,50, by = 10)) #.4,.9, by = .1),
      x.labels <- x.tcks
      x.labels[!x.labels %in% c(1,10,50)] <- NA
      axis(1, at = log(x.tcks), label = x.labels)
    }else{
      if(xtype=='prev') { ## prevalence
        x.tcks <- c(seq(.01, .1, by = .01), seq(.1, .5, by = .1))
        x.labels <- x.tcks
        x.labels[!x.labels %in% c(.01,.1,.5)] <- NA
      }else{ ## pre- & extra- contact mixing coefficients or betas (transmission coefficients)
        x.tcks <- c(seq(.1,.9, by = .1), 1:9, seq(10,100, by = 10))
        x.labels <- x.tcks
        x.labels[!x.labels %in% c(.01,.1,1,10,100)] <- NA
      }
    }
    axis(1, at = log(x.tcks), label = x.labels)
  }
}

outdir <- file.path('results','PrevFigs')
if(!file.exists(outdir)) dir.create(outdir)
yvars <- c('lgt.cprev','lgt.pprev','lgt.psdc','lgt.fD')
ytypes <- c('prev','prev','SDP','fancyD')
ynames <- c('DHS prev','peak prev','SDP','fancyD')
ytexts <- c('DHS HIV prevalence in couples', 'peak HIV prevalence (UNAIDS)','DHS serodiscordant proportion', 'fancy D')
xclasses <- c('beta','contact')
ylims <- list(c(.01, .55), c(.01, .55), c(.3, .99), c(.3, .99))
for(yv in 1:4) { ##  for each Y variable
  for(xc in 1:2) { ##  and each X class
    xclass <- xclasses[xc]
    if(xclass=='contact') { # name HIV infectivity & contact coefficients
      xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f')
      xlabs <- c('HIV infectivity (per 100 person-years)', expression(c['extra-couple']),
                 expression(c['pre-couple,male']), expression(c['pre-couple,female']))    
    }else{                    # name HIV infectivity & other betas
      xvars <- c('lbp','lbe','lbmb','lbfb')
      xlabs <- c('HIV infectivity (per 100 person-years)', expression(beta['extra-couple']),
                 expression(beta['pre-couple,male']), expression(beta['pre-couple,female']))
    }
    xtypes <- c('within', rep('contact',3)) # for axis function
    xlim <- list(c(.3,50),c(.08,100),c(.08,100),c(.08,100)) # xlimits for each x variable, same regardless of class
    ytype <- ytypes[yv]                                     # yvar type
    yvar <- rdat[,yvars[yv]]                                # yvar
    pdf(file.path(outdir,paste0(ynames[yv],' vs ', xclass,'.pdf')), w = 6.5, h = 4.5) # initialize PDF
    for(aa in acs.to.do) {     # for each assumed acute phase relative hazard
      rdat <- rdatl[[aa]]           # pick fitted parameters from array for that acute phase relative hazard
      par(mar = c(4,3,2,.5), mfrow = c(2,2), oma = c(0,1,1,0), cex.lab = .8)
      for(xv in 1:4) { ## for each plot panel (each x variable)
        xvar   <- rdat[,xvars[xv]] # get xvariable
        xvar.l <- rdat[,paste0(xvars[xv],'.l')] # get xvariable lower credible limit
        xvar.u <- rdat[,paste0(xvars[xv],'.u')] # get xvariable upper credible limit      
        plot(0,0, type = 'n', xlab = xlabs[xv], ylab = ylab, las = las, yaxt = yaxt, xaxt = xaxt, # initialize plot
             xlim = log(xlim[[xv]]), ylim = logit(ylims[[yv]]), bty = 'n')
        axis.fxn(T,T, ytype=ytype, xtype=xtypes[xv]) # add axes using function above
        points(xvar, yvar, pch = 19, col = cols) # medians
        arrows(xvar.l, yvar, xvar.u, yvar, code = 3, col = cols, len = .05, angle = 90) # credible intervals
        ## weighted regression model (variance in xvars, from fitting transmission coefficients)
        mod <- lm(yvar ~ xvar, rdat, weights = 1/((xvar.u - xvar.l)/(1.96*2))^2) 
        newx<- seq(min(xvar.l,na.rm=T), max(xvar.u,na.rm=T), l = 120) # sequence of x variables over which to predict
        prd<-predict(mod,newdata=data.frame(xvar = newx),interval = c("confidence"), level = 0.95,type="response") # predict
        ## plot thick solid regression line with thin dashed CI lines
        for(pl in 1:3)    lines(newx,prd[,pl], lty=ifelse(pl==1,1,2), lwd = ifelse(pl==1,2,1)) 
        if(r2)  mtext(paste(expression(P), '=',signif(summary(mod)[[5]][2,4],2)),
                      side = 3, line = 0, adj = .95, cex = .7) # add p value to plot
      }
      legend(.95*log(xlim[[xv]][1]), 1*logit(ylims[[yv]][2]), ncol = 2, # country legend
             rev(ds.nm[leg.ord]), col = rev(cols[leg.ord]), pch = 19, cex = .5, bg = 'white', title = 'from top to bottom')
      mtext(ytexts[yv], side = 2, outer = T, adj = .6, line = -.3) # add one y label
      mtext(paste('acute RR =', in.arr[aa,1,2],'during fit'), side = 3, line = -.5, outer=T) # show assumed acute relative hazard
    }
    dev.off()
  }
}

acs.to.do <- which(rowSums(!is.na(in.arr[,,2]))>0)
fg.col <- 'black'
wid <- 6.5
hei <- 4.5
cex <- 2
####################################################################################################
## New version of figures with color scale for SDP & no arrows for CI's
rmp <- colorRamp(c("red","yellow"))     #create color ramp
outdir <- file.path('results','PrevFigsNew')
if(!file.exists(outdir)) dir.create(outdir)
yvars <- c('lgt.cprev','lgt.pprev','lgt.psdc','lgt.fD')
ytypes <- c('prev','prev')#,'SDP','fancyD')
ynames <- c('DHS prev','peak prev')#,'SDP','fancyD')
ytexts <- c('DHS HIV prevalence in couples', 'peak HIV prevalence (UNAIDS)')#,'DHS serodiscordant proportion', 'fancy D')
xclasses <- c('beta','contact')
ylims <- list(c(.01, .55), c(.01, .55), c(.3, .99), c(.3, .99))
for(yv in 1:2) { ##  for each Y variable
    for(xc in 1:2) { ##  and each X class
        xclass <- xclasses[xc]
        if(xclass=='contact') { # name HIV infectivity & contact coefficients
            xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f')
            xlabs <- c('HIV infectivity (per 100 person-years)', expression(c['extra-couple']),
                       expression(c['pre-couple,male']), expression(c['pre-couple,female']))
            xlabs <- c('HIV infectivity (per 100 person-years)', 'extra-couple mixing coefficient',
                       'male pre-couple mixing coefficient', 'female pre-couple mixing coefficient')
            xlim <- list(c(1,50),c(.08,50),c(.08,50),c(.08,50)) # xlimits for each x variable, same regardless of class            
        }else{                    # name HIV infectivity & other betas
            xvars <- c('lbp','lbe','lbmb','lbfb')
            xlabs <- c('HIV infectivity (per 100 person-years)', expression(beta['extra-couple']),
                       expression(beta['pre-couple,male']), expression(beta['pre-couple,female']))
            xlabs <- c('HIV infectivity (per 100 person-years)', 'extra-couple transmission coefficient',
                       'male pre-couple transmission coefficient', 'female pre-couple transmission coefficient')
            xlim <- list(c(.3,50),c(.08,100),c(.08,100),c(.08,100)) # xlimits for each x variable, same regardless of class
        }
        xtypes <- c('within', rep('contact',3)) # for axis function
        ytype <- ytypes[yv]                                     # yvar type
        yvar <- rdat[,yvars[yv]]                                # yvar
        zvar <- rdat[,'psdc']                                     # zvar for point color
        zord <- order(order(zvar))
        rzvar <- max(zvar)-min(zvar)
        zvar <- (zvar - min(zvar))/rzvar                     # scale so it's between 0 & 1
        cols <- rgb(rmp(zvar), alpha = 200,max = 255)
        pdf(file.path(outdir,paste0(ynames[yv],' vs ', xclass,'-',fg.col,'.pdf')), w = wid, h = hei) # initialize PDF
        for(aa in acs.to.do) {     # for each assumed acute phase relative hazard
            rdat <- rdatl[[aa]]           # pick fitted parameters from array for that acute phase relative hazard
            layout(matrix(c(1,3,2,4,5,5),2,3), widths = c(1,1,.6), h = rep(1, 3))
            par(mar = c(4,4,4,.5), cex.lab = .8, fg = fg.col, col.axis=fg.col, col.lab=fg.col)
            for(xv in 1:4) { ## for each plot panel (each x variable)
                xvar   <- rdat[,xvars[xv]] # get xvariable
                xvar.l <- rdat[,paste0(xvars[xv],'.l')] # get xvariable lower credible limit
                xvar.u <- rdat[,paste0(xvars[xv],'.u')] # get xvariable upper credible limit      
                plot(0,0, type = 'n', ylab = ylab, xlab='', las = las, yaxt = yaxt, xaxt = xaxt, # initialize plot
                     xlim = log(xlim[[xv]]), ylim = logit(ylims[[yv]]), bty = 'n')
                mtext(xlabs[xv], side = 1, line = 2.3, adj = .5, cex = .75)
                axis.fxn(T,T, ytype=ytype, xtype=xtypes[xv]) # add axes using function above
                ## arrows(xvar.l, yvar, xvar.u, yvar, code = 3, col = cols, len = .05, angle = 90) # credible intervals
                ## weighted regression model (variance in xvars, from fitting transmission coefficients)
                mod <- lm(yvar ~ xvar, rdat, weights = 1/((xvar.u - xvar.l)/(1.96*2))^2) 
                newx<- seq(min(xvar.l,na.rm=T), max(xvar.u,na.rm=T), l = 120) # sequence of x variables over which to predict
                prd<-predict(mod,newdata=data.frame(xvar = newx),interval = c("confidence"), level = 0.95,type="response") # predict
                ## plot thick solid regression line with thin dashed CI lines
                for(pl in 1:3)    lines(newx,prd[,pl], lty=ifelse(pl==1,1,2), lwd = ifelse(pl==1,2,1))
                points(xvar, yvar, cex = cex, pch = 19, col = cols) # medians
                text(xvar, yvar, c(16:1)[zord], cex = .5)
                if(r2)  mtext(paste(expression(P), '=',signif(summary(mod)[[5]][2,4],2)),
                              side = 3, line = 0, adj = .95, cex = .7) # add p value to plot
            }
            par(mar=rep(0,4))
            plot(0,0, xlim = c(0, 10), ylim = c(0, 100), type = 'n', bty = 'n', axes = F, xlab = '', ylab ='')
            points(rep(.5,16), seq(40, 90, l = 16), pch = 19, cex = cex, col = cols[order(zvar)])
            text(.5, seq(40, 90, l = 16), 16:1, cex = cex*.3)
            rdat$country[rdat$country=='WA'] <- 'West Africa'
            text(1, seq(40, 90, l = 16), paste0(rdat$country[order(zvar)],
                                    ' (',signif(rdat$psdc[order(zvar)],2)*100,'%)'), pos = 4)
            text(5, 100, 'country \n(serodiscordant \nproportion as %)')
            ## legend('topleft', ncol = 1, paste0(rdat$country, ' (',signif(zvar,2)*100,'%)'),
            ##        pch = 19, col = cols, pt.cex = 1.5, cex = .8, bg = 'white')
            ## legend('topleft', ncol = 1, paste0(rdat$country, ' (',signif(zvar,2)*100,'%)'),
            ##        pch = as.character(1:16), col = 'black', pt.cex = .5, cex = .8, bty = 'n')
            mtext(ytexts[yv], side = 2, outer = T, adj = .6, line = -1.5) # add one y label
            mtext(paste('acute RR =', in.arr[aa,1,2],'during fit'), side = 3, line = -2, outer=T) # show assumed acute relative hazard
        }
        dev.off()
    }
}
####################################################################################################

####################################################################################################
####################################################################################################
####################################################################################################
## for PPT
r2 <- T                                 # show p values on regression plots
acs.to.do <- which(in.arr[,2,2]==7)
fg.col <- 'black'
wid <- 5.5
hei <- 3.5
cex <- 2
####################################################################################################
## New version of figures with color scale for SDP & no arrows for CI's
rmp <- colorRamp(c("red","yellow"))     #create color ramp
outdir <- file.path('results','PrevFigsNew')
if(!file.exists(outdir)) dir.create(outdir)
yvars <- c('lgt.cprev','lgt.pprev','lgt.psdc','lgt.fD')
ytypes <- c('prev','prev')#,'SDP','fancyD')
ynames <- c('DHS prev','peak prev')#,'SDP','fancyD')
ytexts <- c('HIV prevalence', 'peak HIV prevalence')#,'DHS serodiscordant proportion', 'fancy D')
xclasses <- c('contact')
ylims <- list(c(.01, .5), c(.01, .5), c(.3, .99), c(.3, .99))
for(yv in 2) { ##  for each Y variable
  for(xc in 1) { ##  and each X class
    xclass <- xclasses[xc]
    if(xclass=='contact') { # name HIV infectivity & contact coefficients
      xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f')
      xlabs <- c('intrinsic HIV transmission rate \n(per 100 person-years)', expression(c['extra-couple']),
                 expression(c['pre-couple,male']), expression(c['pre-couple,female']))
      xlabs <- c('intrinsic HIV transmission rate \n(per 100 person-years)', 'extra-couple mixing coefficient',
                 'male pre-couple mixing coefficient', 'female pre-couple mixing coefficient')
      xlim <- list(c(1,50),c(.08,50),c(.08,50),c(.08,50)) # xlimits for each x variable, same regardless of class            
    }else{                    # name HIV infectivity & other betas
      xvars <- c('lbp','lbe','lbmb','lbfb')
      xlabs <- c('HIV transmission rate (per 100 person-years)', expression(beta['extra-couple']),
                 expression(beta['pre-couple,male']), expression(beta['pre-couple,female']))
      xlabs <- c('HIV transmission rate (per 100 person-years)', 'extra-couple transmission coefficient',
                 'male pre-couple transmission coefficient', 'female pre-couple transmission coefficient')
      xlim <- list(c(.3,50),c(.08,100),c(.08,100),c(.08,100)) # xlimits for each x variable, same regardless of class
    }
    xtypes <- c('within', rep('contact',3)) # for axis function
    ytype <- ytypes[yv]                                     # yvar type
    yvar <- rdat[,yvars[yv]]                                # yvar
    zvar <- rdat[,'psdc']                                     # zvar for point color
    zord <- order(order(zvar))
    rzvar <- max(zvar)-min(zvar)
    zvar <- (zvar - min(zvar))/rzvar                     # scale so it's between 0 & 1
    cols <- rgb(rmp(zvar), alpha = 200,max = 255)
    for(aa in acs.to.do) {     # for each assumed acute phase relative hazard
      rdat <- rdatl[[aa]]           # pick fitted parameters from array for that acute phase relative hazard
      for(xv in 1:4) { ## for each plot panel (each x variable)
        pdf(file.path(outdir,paste0(ynames[yv],' vs ', xclass,xv,'-',fg.col,'.pdf')), w = wid, h = hei) # initialize PDF
        layout(matrix(c(1,2),1,2), widths = c(1,.4)) 
        par(mar = c(4,4,1,.5), oma = c(0,0,0,0), cex.lab = .8, fg = fg.col, col.axis=fg.col, col.lab=fg.col)
        xvar   <- rdat[,xvars[xv]] # get xvariable
        xvar.l <- rdat[,paste0(xvars[xv],'.l')] # get xvariable lower credible limit
        xvar.u <- rdat[,paste0(xvars[xv],'.u')] # get xvariable upper credible limit      
        plot(0,0, type = 'n', ylab = ytexts[yv], xlab='', las = las, yaxt = yaxt, xaxt = xaxt, # initialize plot
             xlim = log(xlim[[xv]]), ylim = logit(ylims[[yv]]), bty = 'n')
        mtext(xlabs[xv], side = 1, line = 2.3, adj = .5, cex = .75)
        axis.fxn(T,T, ytype=ytype, xtype=xtypes[xv]) # add axes using function above
        ## arrows(xvar.l, yvar, xvar.u, yvar, code = 3, col = cols, len = .05, angle = 90) # credible intervals
        ## weighted regression model (variance in xvars, from fitting transmission coefficients)
        mod <- lm(yvar ~ xvar, rdat, weights = geom.w)
        newx<- seq(min(xvar.l,na.rm=T), max(xvar.u,na.rm=T), l = 120) # sequence of x variables over which to predict
        prd<-predict(mod,newdata=data.frame(xvar = newx),interval = c("confidence"), level = 0.95,type="response") # predict
        ## plot thick solid regression line with thin dashed CI lines
        for(pl in 1:3)    lines(newx,prd[,pl], lty=ifelse(pl==1,1,2), lwd = ifelse(pl==1,2,1))
        points(xvar, yvar, cex = cex, pch = 19, col = cols) # medians
        text(xvar, yvar, c(16:1)[zord], cex = .5)
        if(r2)  mtext(paste(expression(P), '=',signif(summary(mod)[[5]][2,4],2)),
                      side = 3, line = -2, adj = .25, cex = .7) # add p value to plot
        par(mar=rep(0,4), cex = .8)
        plot(0,0, xlim = c(0, 10), ylim = c(0, 100), type = 'n', bty = 'n', axes = F, xlab = '', ylab ='')
        points(rep(.5,16), seq(0, 85, l = 16), pch = 19, cex = cex, col = cols[order(zvar)])
        text(.5, seq(0, 85, l = 16), 16:1, cex = cex*.3)
        rdat$country[rdat$country=='WA'] <- 'West Africa'
        text(1, seq(0, 85, l = 16), paste0(rdat$country[order(zvar)],  ' (',signif(rdat$psdc[order(zvar)],2)*100,'%)'),
                            pos = 4, cex = cex*.5)
        text(5, 95, 'country \n(serodiscordant \nproportion as %)', cex = .9)
        par(xpd=NA)
        arrows(8,0,8,85, code = 2, len = .05, lwd = 3)
        mtext('serodiscordant proportion', side = 4, line = -1.3, cex = .8, adj = .4)
        dev.off()
      }
    }
  }
}

####################################################################################################
## Full model & univariate models
acs.to.do <- 1:8
prev.type <- 'lgt.pprev'
for(aa in acs.to.do)  {
    modlist <- list(NA)
    rdat <- rdatl[[aa]]
    xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f')
    form4c <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    modlist[[5]] <- lm(form4c, dat = rdat,  weights = geom.w)
    modtab <- cbind(summary(modlist[[5]])[[5]], confint(modlist[[5]]))
    univtab <- modtab[-1,]                  # get rid of intercept row
    univtab[,] <- NA
    for(xx in 1:4) {
        xvar <- xvars[xx]
        temp.form <- as.formula(paste(prev.type, ' ~', xvar))
        modlist[[xx]] <- lm(temp.form, dat = rdat, weights = geom.w)
        univtab[xx,] <- cbind(summary(modlist[[xx]])[[5]], confint(modlist[[xx]]))[-1,]
    } 
    tab <- aictab(modlist, second.ord=T, modnames = c(xvars,'mult'))
    tab[,-1] <- signif(tab[,-1],3)
    write.csv(tab, file.path(outdir, paste0('aicc table Ac',acutes[aa],'.csv')))
    un.est <- c(NA,paste0(signif(univtab[,'Estimate'],2), ' (', signif(univtab[,'2.5 %'],2), ', ', signif(univtab[,'97.5 %'],2), ')'))
    un.p <-  c(NA,signif(univtab[,'Pr(>|t|)'],2))
    mult.est <- paste0(signif(modtab[,'Estimate'],2), ' (', signif(modtab[,'2.5 %'],2), ', ', signif(modtab[,'97.5 %'],2), ')')
    mult.p <-  signif(modtab[,'Pr(>|t|)'],2)
    outtab <- data.frame('univariate'=un.est, 'P'=un.p, 'multivariate'=mult.est, 'P'=mult.p)[-1,]
    rownames(outtab) <- c('transmission rate','extra-couple contact coefficient',
                          'male pre-couple contact coefficient', 'female pre-couple contact coefficient')
    outtab
    write.csv(outtab, file.path(outdir, paste0('final model table Ac',acutes[aa],'.csv')))
    print(paste('R2 = ', summary(modlist[[1]])$r.squared, 'for acute = ', acutes[aa]))
  }

####################################################################################################
## Show evolution of {SDP,Prev} ove time
ds <- function(x) (2*x*(1-x))/(2*x*(1-x) + x^2) # Null SDP model based on mixing SDP = (2p(1-p)) / (2p(1-p) + p^2)
col <- 'black'
cc <- F
xlim <- c(0,.3)
ylim <- c(0,1)
pdf(file.path(outdir, paste('SDP vs prevalence by survey over time','curve'[cc],'lm.pdf')), w = 7, h=4)
par(mar=c(4,5,1,1), col = col, fg = col, col.axis = col, col.main = col, col.lab=col, bty = 'l')
plot(0,0, type = 'n', xlim=xlim, ylim = ylim, xlab = 'HIV prevalence', ylab = 'serodiscordance proportion', axes=F)
if(cc) curve(ds(x), from = 0, to = .3, col = 'red', lty = 2, add=T)
axis(2, seq(0,1,l=5), las = 2)
axis(1, seq(0,.3,by=.05), las = 1)
w.africa <- c("Burkina Faso", "Cameroon", 'Cote dIvoire', "Ghana", "Guinea", "Liberia", "Mali", "Niger", "Senegal",
              "Sierra Leone")
temp <- draw.s[!draw.s$country %in% w.africa & !is.na(draw.s$country),]
temp$country <- factor(temp$country)
temp$col <- rainbow(nlevels(temp$country), alpha = .7)[as.numeric(temp$country)]
points(temp$prev, temp$psdc, pch = 19, col = temp$col, cex = 2.5)
leg <- as.character(unique(temp$country))
for(cc in 1:nlevels(temp$country)) {
  sel <- temp$country==levels(temp$country)[cc]
  leg[cc] <- paste(leg[cc], paste(round(temp$tint.yr[sel]), collapse = ', '))
  lines(temp$prev[sel], temp$psdc[sel], col = temp$col[sel][1])
  for(jj in which(sel))  text(temp$prev[sel], temp$psdc[sel], substr(round(temp$tint.yr[sel]),3,4), cex = .6) #, col = temp$col[sel])
}
abline(lm(temp$psdc ~ temp$prev), col = 'black')
legend('topright', leg = leg, col = unique(temp$col), pch = 19, bty = 'n', cex = .7, ncol=2)
dev.off()


