
## Full 9 variable AIC combinations
xvars.all <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f', 'lmysa', 'lfysa', 'lmardur','lmpret','lfpret')
nvar <- length(xvars.all)
cbn <- combn(nvar,nvar)
for(ii in (nvar-1):1) {
  to.add <- combn(nvar,ii)
  for(jj in 1:(nvar-ii)) to.add <- rbind(to.add, NA)
  cbn <- cbind(cbn, to.add)
}
acs.to.do <- which(rowSums(!is.na(in.arr[,,2]))>0)
typ <- c('pprev','cprev')
modtab <- list(NA)
modtab2 <- list(NA)
for(pt in 1:2) {
  prev.type <- paste0('lgt.',typ[pt])
####################################################################################################
  ## Multivariate model: contact * time at risk for all 4 variables
  for(aa in acs.to.do) {
    rdat <- rdatl[[which(in.arr[,1,2]==acutes[aa])]] # 
    modlist <- list(NA)
    for(cb in 1:ncol(cbn)) {
      xvars <- xvars.all[cbn[,cb]]
      xvars <- xvars[!is.na(xvars)]
      form.temp <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
      modlist[[cb]] <- lm(form.temp, dat = rdat,  weights = geom.w)
    }
    tab <- head(aictab(modlist, second.ord=T, modnames = paste0('mod',1:ncol(cbn))))
    ## look at best model
    if(pt==1) {
      if(aa==1) tab.arr <- tab else tab.arr <- abind(tab.arr, tab, along = 3)
      modtab[[aa]] <- summary(modlist[[as.numeric(sub('mod','',tab[1,1]))]])[[5]]
    }else{
      if(aa==1) tab.arr2 <- tab else tab.arr2 <- abind(tab.arr2, tab, along = 3)
      modtab2[[aa]] <- summary(modlist[[as.numeric(sub('mod','',tab[1,1]))]])[[5]]      
    }
  }
}
tab.arr <- abind(tab.arr,tab.arr2, along = 4)
tab.arr[,,3,]
tab.arr[,,,]
modtab
modtab2

acs.to.do <- which(rowSums(!is.na(in.arr[,,2]))>0)
typ <- c('pprev','cprev')
for(pt in 1:2) {
  prev.type <- paste0('lgt.',typ[pt])
####################################################################################################
  ## Multivariate model: contact * time at risk for all 4 variables
  for(aa in acs.to.do) {
    rdat <- rdatl[[which(in.arr[,1,2]==acutes[aa])]] # 
    ## contact*time for all four variables
    xvars <- c('linft','lextrat','lmpret','lfpret')
    form4t <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod4t <- lm(form4t, dat = rdat,  weights = geom.w)
    ##summary(mod4t) ## all * time
    ##hier.part(rdat[,prev.type], rdat[,xvars])
    ## contact + time for all 3 time variables
    xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f', 'lmysa', 'lfysa', 'lmardur')
    form7t <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod7t <- lm(form7t, dat = rdat,  weights = geom.w)
    ##summary(mod7t) ## all * time
    ##hier.part(rdat[,prev.type], rdat[,xvars])
    ## contact + time for all 2 time variables
    xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f', 'lmysa', 'lfysa')
    form6t <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod6t <- lm(form6t, dat = rdat,  weights = geom.w)
    ##summary(mod6t) ## all * time
    ##hier.part(rdat[,prev.type], rdat[,xvars])
    ## contacts only
    xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f')
    form4c <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod4c <- lm(form4c, dat = rdat,  weights = geom.w)
    ##summary(mod4c) ## all * time
    ##hier.part(rdat[,prev.type], rdat[,xvars])
    ## contact*time for pre-couple only
    xvars <- c('lbp','lrr.ep','lmpret','lfpret')
    form2t <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod2t <- lm(form2t, dat = rdat,  weights = geom.w)
    ## transmission only
    xvars <- c('lbp')
  #  xvars <- c('lbp','lrr.ep','lmpret','lfpret')    
    form1b <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod1b <- lm(form1b, dat = rdat,  weights = geom.w)
    ##summary(mod2t) ## pre * time, others not
    ##hier.part(rdat[,prev.type], rdat[,xvars], gof='Rsqu') 
    tab <- aictab(list(mod4c, mod4t, mod2t, mod7t, mod6t, mod1b), second.ord=T, modnames = c('mod4c', 'mod4t', 'mod2t', 'mod7t', 'mod6t','mod1b')) #, modnames = nms)
    tab[,3:8] <- signif(tab[,3:8],3)
    ##    tab <- tab[,-c(5,7,8)]
    ## Get names of X variables from best model
    xvars <- names(summary(get(as.character(tab[1,1])))[[5]][,1])[-1]
    modtab <- summary(get(as.character(tab[1,1])))[[5]]
    if(length(xvars)>1) {
      htab <- hier.part(rdat[,prev.type], rdat[,xvars], gof='Rsqu')$IJ
      phtab <- hier.part(rdat[,prev.type], rdat[,xvars], gof='Rsqu')$I.perc
    }else{
      htab <- summary(get(as.character(tab[1,1])))$r.squared
      phtab <- 1
    }
    if(pt==1) {
      if(aa==1) modtab.arr <- modtab else modtab.arr <- abind(modtab.arr, modtab, along = 3)          
      if(aa==1) tab.arr <- tab else tab.arr <- abind(tab.arr, tab, along = 3)    
      if(aa==1) htab.arr <- htab else htab.arr <- abind(htab.arr, htab, along = 3)
      if(aa==1) phtab.arr <- phtab else phtab.arr <- abind(phtab.arr, phtab, along = 3)
    }else{
      if(aa==1) modtab.arr2 <- modtab else modtab.arr2 <- abind(modtab.arr2, modtab, along = 3)                
      if(aa==1) tab.arr2 <- tab else tab.arr2 <- abind(tab.arr2, tab, along = 3)    
      if(aa==1) htab.arr2 <- htab else htab.arr2 <- abind(htab.arr2, htab, along = 3)
      if(aa==1) phtab.arr2 <- phtab else phtab.arr2 <- abind(phtab.arr2, phtab, along = 3)
    }
  }
}
modtab.arr <- abind(modtab.arr, modtab.arr2, along = 4)
tab.arr <- abind(tab.arr, tab.arr2, along = 4)
htab.arr <- abind(htab.arr, htab.arr2, along = 4)
phtab.arr <- abind(phtab.arr, phtab.arr2, along = 4)


modtab.arr
tab.arr
htab.arr
phtab.arr

modtab.arr[,,3,]
tab.arr[,,3,]

typ <- c('pprev','cprev')
  prev.type <- paste0('lgt.',typ)[1]
for(aa in acs.to.do)
  {
    rdat <- rdatl[[aa]]
    xvars <- c('lbp','lrr.ep','lrr.bp.m','lrr.bp.f')
    form4c <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod4c <- lm(form4c, dat = rdat,  weights = geom.w)
    summary(mod4c)
    xvars <- c('lrr.ep','lrr.bp.m','lrr.bp.f')
    form3c <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod3c <- lm(form3c, dat = rdat,  weights = geom.w)
    summary(mod3c)
    xvars <- c('lbp')
    form1b <- as.formula(paste(prev.type, ' ~', paste(xvars, collapse = '+')))
    mod1b <- lm(form1b, dat = rdat,  weights = geom.w)
    summary(mod1b)
    tab <- aictab(list(mod4c, mod3c, mod1b), second.ord=T, modnames = c('mod4c', 'mod3c', 'mod1b'))
    print(tab)
  }

    ## sel <- !rdat$country %in% c('Burundi','Rwanda','Ethiopia')
    ## mod2t.ol <- lm(form2t, dat = rdat[sel,],  weights = geom.w)
    ## summary(mod2t.ol) ## pre * time, others not
    ## hier.part(rdat$lgt.cprev[sel], rdat[sel,xvars], gof='Rsqu')

    names(tab) <- c('model', 'df', 'AICc', 'deltaAICc', 'AICc weight')
    write.csv(tab, file.path(outdir, 'aicc table.csv'))

    ## This is the best model so make it into a table output
    names(summary(mod2t))
    modtab <- cbind(summary(mod2t)[[5]], confint(mod2t))
    ## Add univariate estimates
    univtab <- modtab[-1,]                  # get rid of intercept row
    for(ii in 1:4) {
        xv <- xvars[ii]
        xvar   <- rdat[,xv] # get xvariable
        ## xvar.l <- rdat[,paste0(xv,'.l')] # get xvariable lower credible limit
        ## xvar.u <- rdat[,paste0(xv,'.u')] # get xvariable upper credible limit      
        mod.temp <- lm(as.formula(paste('lgt.cprev ~', xv)), dat = rdat,  weights = 1/((lbp.u - lbp.l)/(1.96*2))^2)
        univtab[ii,] <- cbind(summary(mod.temp)[[5]], confint(mod.temp))[xv,]
    }
    un.est <- c(NA,paste0(signif(univtab[,'Estimate'],2), ' (', signif(univtab[,'2.5 %'],2), ', ', signif(univtab[,'97.5 %'],2), ')'))
    un.p <-  c(NA,signif(univtab[,'Pr(>|t|)'],2))
    mult.est <- paste0(signif(modtab[,'Estimate'],2), ' (', signif(modtab[,'2.5 %'],2), ', ', signif(modtab[,'97.5 %'],2), ')')
    mult.p <-  signif(modtab[,'Pr(>|t|)'],2)
    outtab <- data.frame('univariate'=un.est, 'P'=un.p, 'multivariate'=mult.est, 'P'=mult.p)
    htab <- rbind(NA, signif(hier.part(rdat$lgt.cprev, rdat[,xvars], gof='Rsqu')$I.perc,2))
    outtab <- data.frame(outtab, htab)
    names(outtab)[5] <- 'independent variance \ncontributed'
    outtab
    write.csv(outtab, file.path(outdir, 'final model table.csv'))


####################################################################################################
## Multivariate mod3el: contact
rdat <- rdatl[[which(in.arr[,1,2]==7)]] # acute==7

mod3 <- lm(lgt.cprev ~ lbp + lrr.ep + lrr.bp.m + lrr.bp.f, rdat,  weights = geom.w)
#mod2 <- lm(lgt.cprev ~ lbp + lrr.ep, rdat,  weights = geom.w)
summary(mod3)
hier.part(rdat$lgt.cprev, rdat[,xvars])
hier.part(rdat$lgt.cprev, rdat[,xvars[1:2]])

####################################################################################################
## Multivariate model: betas
rdat <- rdatl[[which(in.arr[,1,2]==7)]] # acute==7
mod2 <- lm(lgt.cprev ~ lbp + lbe + lbmb + lbfb, rdat,  weights = geom.w)
xvars <- c('lbp','lbe','lbmb','lbfb')
summary(mod2)
sel <- !is.na(rdat[,xvars[1]])
hier.part(rdat$lgt.cprev, rdat[,xvars])
hier.part(rdat$lgt.cprev, rdat[,xvars[1:2]])

####################################################################################################
## Do we expect SDP to be affected by prevalence? Remake Chemaitelly et al. (2012)
## no points null mixing model
## sdp w/o w-c infection
ds <- function(x) (2*x*(1-x))/(2*x*(1-x) + x^2) # Null SDP model based on mixing SDP = (2p(1-p)) / (2p(1-p) + p^2)
col <- 'black'
for(cc in c(T,F)) {
    xlim <- c(0,.3)
    ylim <- c(0,1)
    pdf(file.path(outdir, paste('SDP vs prevalence','curve'[cc],'lm.pdf')), w = 7, h=4)
    par(mar=c(4,5,1,1), col = col, fg = col, col.axis = col, col.main = col, col.lab=col, bty = 'l')
    plot(0,0, type = 'n', xlim=xlim, ylim = ylim, xlab = 'HIV prevalence', ylab = 'serodiscordance proportion', axes=F)
    if(cc) curve(ds(x), from = 0, to = .3, col = 'red', lty = 2, add=T)
    axis(2, seq(0,1,l=5), las = 2)
    axis(1, seq(0,.3,by=.05), las = 1)
    temp <- draw.s
    points(temp$prev, temp$psdc, pch = 19, col = temp$col)
    abline(lm(temp$psdc ~ temp$prev), col = 'black')
    legend('topright', leg = ds.nm, col = cols, pch = 19, bty = 'n', cex = .7)
    dev.off()
  }
