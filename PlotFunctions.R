####################################################################################################
## Functions for collecting, analyzing, & visualizing simulation results
####################################################################################################
## Steve Bellan, 2013
####################################################################################################

####################################################################################################
## Collect simulation results named by 'files' into arrays/dataframes
coll <- function(files, give.ev = T,    # files to collect; return 'evout' (line list) array?
                 nc = 12, browse = F)   # number cores; debug
  {
    if(browse) browser()              # debug
    ## make batches to divide between cores to increase speed
    batchs <- rep(1:nc, l = length(files))
    wrp <- function(bb, files = files, browse) { # wrapper function to feed to mclapply() for parallel processing
        if(browse) browser()
      files <- files[batchs == bb] # select batch
      for(ff in 1:length(files)) { # loop through files
        if(ff%%5 == 0) print(paste('collecting result file', ff, 'of', length(files))) # update on progress
        load(files[ff]) 
        if(ff==1) {     # initialize cframe (description of jobs), t.arr (time
                                        # serieses in array), e.arr (line lists in array)
          t.arr.temp <- output$tss[output$tss$tt > 62*12, ]
          if(max(t.arr.temp$tt) < 1356) { ## add rows up to 2013, do avoid doing all simulations that just did to 2011
              temp <- t.arr.temp[1:(1356-max(t.arr.temp$tt)),] ## copy df structure
              temp[,] <- NA
              temp$tt <- (max(t.arr.temp$tt)+1):1356
              t.arr.temp <- rbind(t.arr.temp, temp)
          }
         t.arr.temp$sdp <- (t.arr.temp[,'mm'] + t.arr.temp[,'ff']) / t.arr.temp[,'inf.alive']          
          t.arr <- t.arr.temp
          cframe <- data.frame(t(c(output$jobnum, output$simj, rakacRR = output$rakacRR, output$pars, output$hours,
                                   sdp08 = t.arr.temp[which(t.arr.temp$yr==2008),'sdp'])))
          ## causes of infection (m/f): before = 1, extra = 2, partner = 3, not infected = 4.
          output$ev$mcoi <- as.numeric(output$ev$mcoi) 
          output$ev$fcoi <- as.numeric(output$ev$fcoi)
          if(give.ev) { # if returning line list
            ## Next if statement: if parametrically sampling couple formation (marital) dates,
            ## select 98% of first line list, and sample that from the
            ## populations to account fo the fact that each simulation has
            ## slightly different # of couples (otherwise can't store in
            ## array due to different sizes.
            if(output$pars['sample.tmar']==1 & output$pars['psNonPar']==0) {
              maxrow <- round(.98*nrow(output$ev))
              samp <- sample(1:nrow(output$ev), maxrow )
            }else{ samp <- 1:nrow(output$ev)} # return all data
            col.nms <- c( "uid", "ser", "tms", "tfs", "tmar", # variables to return
                         "tint", "mardur.mon","mage", "fage", "epic.ind", 
                         "mser", "fser", "mdoi", "fdoi", "mdod", 
                         "fdod", "mcoi", "fcoi",'dage','dmage','dfage',
                         'm.het.b','f.het.b','m.het.e','f.het.e','m.het.p','f.het.p',
                         "m.het.gen", "f.het.gen",
                         "m.het.beh", "f.het.beh", "malive", "falive", "alive", 
                         "tend", "tmsdc", "tfsdc", "tccc", "msdc", 
                         "fsdc", "ccc")       
            e.arr <- output$ev[samp,col.nms] # select rows & columns to return
          }
        }else{      # append cframe, t.arr, e.arr
          t.arr.temp <- output$tss[output$tss$tt > 62*12, ]
         t.arr.temp$sdp <- (t.arr.temp[,'mm'] + t.arr.temp[,'ff']) / t.arr.temp[,'inf.alive']
          cframe <- rbind(cframe, data.frame(t(c(output$jobnum, output$simj, rakacRR = output$rakacRR,
                                                 output$pars, output$hours, sdp08 = t.arr.temp[which(t.arr.temp$yr==2008),'sdp'])))) # append
          if(max(t.arr.temp$tt) < 1356) { ## add rows up to 2013, do avoid doing all simulations that just did to 2011
              temp <- t.arr.temp[1:(1356-max(t.arr.temp$tt)),] ## copy df structure
              temp[,] <- NA
              temp$tt <- (max(t.arr.temp$tt)+1):1356
              t.arr.temp <- rbind(t.arr.temp, temp)
          }
          t.arr <- abind(t.arr, t.arr.temp, along = 3) # append
          if(give.ev) { # must be numeric since it's an array
            output$ev$mcoi <- as.numeric(output$ev$mcoi)
            output$ev$fcoi <- as.numeric(output$ev$fcoi)
            if(output$pars['sample.tmar']==1 & output$pars['psNonPar']==0) { 
              samp <- sample(1:nrow(output$ev), maxrow )# downsample as above
            }else{ samp <- 1:nrow(output$ev)}             # complete sample
            e.arr <- abind(e.arr, output$ev[samp,col.nms], along = 3) # append
          }
        }
      }
      if(give.ev) {           
        return(list(cframe = cframe, t.arr = t.arr, e.arr = e.arr))
      }else{
        return(list(cframe = cframe, t.arr = t.arr))
      }
    } ## end multicore wrapper
    temp <- mclapply(1:nc, wrp, files = files, browse = browse) # send wrp tasks to cores
#    browser()
    for(ii in 1:nc) {                          ## collect results from single core outputs
      if(ii==1) {                        ## initialize cframe, t.arr, e.arr
        cframe <- temp[[ii]]$cframe
        t.arr <- temp[[ii]]$t.arr
        if(give.ev)     e.arr <- temp[[ii]]$e.arr
      }else{              # append cframe, t.arr, e.arr
        cframe <- rbind(cframe, temp[[ii]]$cframe)
        t.arr <- abind(t.arr, temp[[ii]]$t.arr, along = 3)
        if(give.ev)     e.arr <- abind(e.arr, temp[[ii]]$e.arr, along = 3)
      }
    }                                                      
    names(cframe)[1] <- 'job'
    jord <- order(cframe$job)
    names(cframe)[2] <- 'simj'
    if(give.ev) { ## return results.
      list(cframe = cframe[jord,], t.arr = t.arr[,,jord], e.arr = e.arr[,,jord])
    }else{ ## return results without line lists.
      list(cframe = cframe[jord,], t.arr = t.arr[,,jord])
    }
  }

####################################################################################################
## Plot time series array file from psrun()
####################################################################################################
plot.tss <- function(js, # jobs to plot (indexing cframe)
                     leg, # legend for plot
                     title = NA,        # legend title
                     pdf.nm,   # pdf name
                     dens = T, # histogrames or densities?
                     cols = NA, ltys = NA, lwds = NA, # colors, line types, line widths
                     make.pdf = T,                    # make pdf? if not plots to default R graphics (quartz)
                     ymax = .35,                       # ymax for time series plots
                     ymax2 = .35, # ymax for density plots
                     rmp = colorRampPalette(c("blue","orange")), # color ramp for different jobs
                     fitted.col = 'green',                       # color for 'as fitted' simulation (non-counterfactual)
                     early.yr = 1985,                            # earliest year to plot in time series
                     browse = F)                                 # debut
  {
    if(browse) browser()                # debug
    tss <- cfs$t.arr[,,js]              # select time series for selected jobs
    tss <- tss[tss[,"yr",1]>early.yr,,] # select time series after earliest year
    if(make.pdf) {                      # make pdf
        pdf(paste(pdf.nm, ".pdf", sep = ""), w = 8, h = 8)
        par(mfrow=c(4,3), oma = c(1,3,2,0), mar = c(3,3,1,1)) # set plot panels, inner/outer margins
      }
    ## Plot prevalence of infected couples
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## set colors, line types/widths from ramp if not specified, 1st job (as fitted to DHS defaults
    ## to green)
    if(is.na(cols[1]))      cols <- rmp(length(js)); cols[js==1] <- fitted.col
    if(is.na(ltys))     ltys <- rep(1, length(js))
    if(is.na(lwds))     lwds <- rep(1, length(js))
    ## for each infected couple type, plot the prevalence of that couple type over time
    mains <- c('M+F-','M-F+','M+F+')          # panel titles
    cts <- c('mm','ff','hh')
    for(ct.i in 1:length(cts)) {
        ct <- cts[ct.i]                 # which couple type to work on
        plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), ## initialize plot
             bty = "n", main = mains[ct.i], xlab="", ylab = "")
        for(ii in 1:dim(tss)[3]) {      # loop through jobs to plot
            lines(tss[,"yr",ii], tss[,ct,ii]/tss[,'alive',ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
        }
        legend("topleft", leg = leg, col = cols, lwd = lwds, lty = ltys, bty = "n")
    }
    ## for each infected couple type, plot the prevalence of that couple type amongst *infected*
    ## couples over time
    for(ct.i in 1:length(cts)) {
        ct <- cts[ct.i]                 # which couple type to work on
        plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), ## initialize plot
             bty = "n", xlab="", ylab = "")
        for(ii in 1:dim(tss)[3]) {      # loop through jobs to plot
            lines(tss[,"yr",ii], tss[,ct,ii]/tss[,'inf.alive',ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
        }
    }
    ## plot proportion of serodiscordant couples that are M+ & F+ over time
        for(ct.i in length(cts[-3])) {
        ct <- cts[ct.i]                 # which couple type to work on
        plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, 1), bty = "n", main = "", xlab="", ylab = "") 
        for(ii in 1:dim(tss)[3]) {      # loop through jobs to plot
            lines(tss[,"yr",ii], tss[,ct,ii]/(tss[,"mm",ii]+tss[,"ff",ii]), col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
        }
    }
    ## Plot country-wide infectious prevalence (from UNAIDS model)
    country <- cfs$e.arr[1,'epic.ind',1] # which country?
    eyr <- early.yr-1900                 # earliest year in yrs since 1900, for CMC indexing
    if(length(unique(cframe$s.epic.ind[js]))==1) {                          # only one epidemic curve used amongst jobs
        plot(1900+1/12*(eyr*12):(110*12), epicf[(eyr*12):(110*12),country], # initialize plot, plot female curve
             col = "black", main = "country-wide infectious prevalence",
             xlim = range(tss[,"yr",1]), ylim = c(0,.5), bty = "n", xlab = "", ylab = "", type = "l")
        lines(1900+1/12*(eyr*12):(110*12), epicm[(eyr*12):(110*12),country], col = "dark gray") # males
        legend("topleft", c("males", "females"), col = c("black","dark gray"), lwd = 1, bty = "n")
      }else{                            # multiple epidemic curves used amongst jobs
        plot(0, 0, type = 'n', col = "black", main = "country-wide infectious prevalence", # initialize plot
             xlim = range(tss[,"yr",1]), ylim = c(0,.5), bty = "n", xlab = "", ylab = "")
        for(ii in 1:dim(tss)[3]) { # for each country, plot m/f epidemic curves in differentcolors with different lty
            e.ind <- which(colnames(epicf) == ds.nm[cframe$s.epic[ii]])
            lines(1900+1/12*(eyr*12):(110*12), epicf[(eyr*12):(110*12),e.ind], col = cols[ii], lty = 1)
            lines(1900+1/12*(eyr*12):(110*12), epicm[(eyr*12):(110*12),e.ind], col = cols[ii], lty = 2)
          }
        legend("topleft", c("males", "females"), lwd = 1, lty = 1:2, bty = "n")        
      }
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## histograms of duration spent in each state
    ev <- cfs$e.arr[,,js] ## select line list subset
    xmax <- 60
    tcs <- c('tmsdc','tfsdc','tccc')
    for(tc.i in 1:length(tcs)) {
        tc <- tcs[tc.i]                 # which couple type to work on
        plot(0,0, type = 'n', bty = 'n', xlim = c(0,30), ylim = c(0,ymax2), xlab="years spent in state", main = "")
        for(ff in 1:length(js)) {
            if(dens) {
                lines(density(ev[,tc,ff]/12, na.rm = T, from = 0), col = cols[ff])
            }else{
                hm <- hist(ev[,tc,ff]/12, breaks = seq(0,ceiling(xmax), by = .5), plot=F)
                lines(hm$mids, hm$density, type = 'l', col = cols[ff])
            }
        }
    }
    ## Add text to outer margins
    mtext("prop. of couples", 2, 0, outer = T, at = .87)
    mtext("prop. of infected couples", 2, 0, outer = T, at = .63)
    mtext("prop. of SDC couples", 2, 0, outer = T, at = .38)
    mtext("density", 2, 0, outer = T, at = .12)
    mtext("time spent in state", 1, 0, outer = T, adj = .5)    
    mtext(ds.nm[cframe$group.ind[1]], side = 3, outer = T, at = .5, line = 0) # country name
    if(make.pdf)        dev.off()                                             # finish pdf if making one
  }

####################################################################################################
## Explore how time spent in a serostatus changes with things
durhet <- function(js,                  # jobs
                   het.type = 'gen',    # type of heterogeneity examined
                   rmp = colorRamp(c("yellow","purple")), # color ramp
                   browse = F)          # debug
  {
    if(browse)        browser()         # debug
    ev <- cfs$e.arr[,,js]               # select subset of line lists
    mhet <- paste('m.het','.',het.type,sep='') # individual risk factors to be examined: m.het.gen, m.het.beh, m.het.b...
    fhet <- paste('f.het','.',het.type,sep='') # ditto for female
    het.sd <- get(paste('cframe$het.',het.type,'.sd',sep=''))[js] # job het.sds for this type of heterogeneity
    het.cor <- get(paste('cframe$het.',het.type,'.cor',sep=''))[js] # inter-partner correlations
    pdf.nm <- paste(ds.nm[group],' het',het.type,'-',sd.temp,'-cor',cframe[js,het.cor],'.jpg',sep='') # make pdf
    jpeg(pdf.nm, w = 1000, h = 1000, quality = 100, pointsize = 18)
    par(mfrow=c(2,2))
    ## mSDC vs fHet vs m.age
    sel <- ev[,'msdc']==1
    pchs <- rep(21, sum(sel))
    pchs[ev[sel,'ccc']==1] <- 19
    mage.inf <- ev[sel,'mage'] - (ev[sel,'tint'] - ev[sel,'mdoi'])
    mage.norm <- mage.inf / (12*100)        # fraction of 100 yrs old
    plot(log(ev[sel,fhet]), ev[sel,'tmsdc']/12, col = rgb(rmp(mage.norm), max = 255), pch = pchs,
         ylab = 'years spent mSDC', xlab = 'log(female risk scalar)', cex = .5, xlim = c(-8,8), ylim = c(0,25), bty = 'n')
    ## fSDC vs mHet vs f.age
    sel <- ev[,'fsdc']==1
    pchs <- rep(21, sum(sel))
    pchs[ev[sel,'ccc']==1] <- 19
    fage.inf <- ev[sel,'fage'] - (ev[sel,'tint'] - ev[sel,'fdoi'])
    fage.norm <- fage.inf / (12*100)        # fraction of 100 yrs old
    plot(log(ev[sel,mhet]), ev[sel,'tfsdc']/12, col = rgb(rmp(fage.norm), max = 255), pch = pchs,
         ylab = 'years spent fSDC', xlab = 'log(male risk scalar)', cex = .5, xlim = c(-8,8), ylim = c(0,25), bty = 'n')
    ## CCC vs fHet vs f.age
    sel <- ev[,'ccc']==1
    pchs <- rep(21, sum(sel))
    pchs[ev[sel,'fsdc']==1] <- 19
    plot(log(ev[sel,fhet]), ev[sel,'tccc']/12, col = rgb(rmp(fage.norm), max = 255), pch = pchs,
         ylab = 'years spent CCC', xlab = 'log(female risk scalar)', cex = .5, xlim = c(-8,8), ylim = c(0,25), bty = 'n')
    ## CCC vs mHet vs m.age
    sel <- ev[,'ccc']==1
    pchs <- rep(21, sum(sel))
    pchs[ev[sel,'msdc']==1] <- 19
    plot(log(ev[sel,mhet]), ev[sel,'tccc']/12, col = rgb(rmp(mage.norm), max = 255), pch = pchs,
         ylab = 'years spent CCC', xlab = 'log(male risk scalar)', cex = .5, xlim = c(-8,8), ylim = c(0,25), bty = 'n')
    dev.off()
  }

coitr <- function(js, # jobs to plot
                  leg, # legend for them
                  pdf.nm = 'test',
                  cols = NA,
                  ltys = NA,
                  lwds = NA,
                  rmp = colorRampPalette(c("blue","orange")),
                  early.yr = 1985, 
                  browse = F)
  {
    if(browse) browser()
    tss <- cfs$t.arr[,,js]
    tss <- tss[tss[,"yr",1]>early.yr,,]
    pdf(paste(pdf.nm, ".pdf", sep = ""), w = 7, h = 7)
    par(mfrow=c(3,2), oma = c(0,2,2,0), mar = c(5,5,2,1))
######################################################################
    ymax <- max(tss[,c('cu.mb','cu.me','cu.mp','cu.fb','cu.fe','cu.fp'),])
    if(is.na(cols))      cols <- rmp(length(js)); cols[js==1] <- 'green'    
    if(is.na(ltys))     ltys <- rep(1, length(js))
    if(is.na(lwds))     lwds <- rep(1, length(js))
    ## cumulative premarital infections    
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "Male", xlab="", mgp = c(3,2,0), ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.mb",ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
      }
    legend("topleft", leg = leg, col = cols, lwd = lwds, lty = ltys, bty = "n")
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "Female", xlab="", ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.fb",ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
      }
    ## cumulative extramarital infections    
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "", xlab="", mgp = c(3,2,0), ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.me",ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
      }
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "", xlab="", ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.fe",ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
      }
    ## cumulative marital infections    
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "", xlab="", mgp = c(3,2,0), ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.mp",ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
      }
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "", xlab="", ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.fp",ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
      }
    mtext('pre-couple # inf', 2, outer = T, line = .5, adj = .93)
    mtext('extra-couple # inf', 2, outer = T, line = .5, adj = .53)
    mtext('within-couple # inf', 2, outer = T, line = .5, adj = .12)    
    dev.off()
  }


coitr2 <- function(js, # jobs to plot
                   leg, # legend for them
                   pdf.nm = 'test',
                   cols = NA,
                   ltys = NA,
                   lwds = NA,
                   rmp = colorRampPalette(c("blue","orange")),
                   early.yr = 1985, 
                   browse = F)
  {
    if(browse) browser()
    tss <- cfs$t.arr[,,js]
    tss <- tss[tss[,"yr",1]>early.yr,,]
    pdf(paste(pdf.nm, ".pdf", sep = ""), w = 7, h = 7)
    par(mfrow=c(3,1), oma = c(0,0,2,0), mar = c(5,5,2,1))
######################################################################
    ymax <- max(tss[,c('cu.mb','cu.me','cu.mp','cu.fb','cu.fe','cu.fp'),])
    if(is.na(cols))      cols <- rmp(length(js)); cols[js==1] <- 'green'    
    if(is.na(ltys))     ltys <- rep(1, length(js))
    if(is.na(lwds))     lwds <- rep(1, length(js))
    ## cumulative premarital infections    
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", ylab = "", xlab="", main = "pre-couple infections")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.mb",ii], col = cols[ii], lty = 1, lwd = lwds[ii])
        lines(tss[,"yr",ii], tss[,"cu.fb",ii], col = cols[ii], lty = 2, lwd = lwds[ii])
      }
    legend("topleft", leg = c(leg,"male",'female'), col = c(cols,rep('black',2)), lwd = c(lwds,1,1), lty = c(ltys,1:2), bty = "n")
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", ylab = "", xlab="", main = "extra-couple infections")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.me",ii], col = cols[ii], lty = 1, lwd = lwds[ii])
        lines(tss[,"yr",ii], tss[,"cu.fe",ii], col = cols[ii], lty = 2, lwd = lwds[ii])
      }
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", ylab = "", xlab="", main = "within-couple infections")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.mp",ii], col = cols[ii], lty = 1, lwd = lwds[ii])
        lines(tss[,"yr",ii], tss[,"cu.fp",ii], col = cols[ii], lty = 2, lwd = lwds[ii])
      }
    dev.off()
  }



## Plot time series array file from psrun()
plot.sdp <- function(js, # jobs to plot
                     leg, # legend for them
                     show.leg = F, cex.leg = 1, leg.col = 1,
                     repl = NA, ## coding for show.sdp
                     show.sdp = T, show.all = T, show.all.pts = T, to.show = NA,
                     remove.wa = F,
                     pdf.nm = 'test',
                     make.pdf = T,
                     dens = T, # histogrames or densities?
                     cols = NA,
                     ltys = NA,
                     lwds = NA,
                     main = NA,
                     ymax2 = .35, # for density plotss
                     rmp = colorRampPalette(c("blue","orange")),
                     ylab = "serodiscordance proportion",
                     early.yr = 1985,
                     col.pl = 'black',
                     sep.leg = F,
                     yaxt = T, xaxt = T,
                     browse = F)
  {
    if(browse) browser()
    if(remove.wa)       js <- js[cframe[js,repl]!=which(ds.nm=='WA')]
    tss <- cfs$t.arr[,,js]
    if(length(js)>1)
      {
        tss <- tss[tss[,"yr",1]>early.yr,,]
        ns <- dim(tss)[3]
        xlim <- range(tss[,"yr",1], na.rm=T)
      }else{
        tss <- tss[tss[,"yr"]>early.yr,]
        ns <- 1
        xlim <- range(tss[,"yr"], na.rm=T)
      }
    xlim[2] <- 2014.5
    if(make.pdf)
      {
        pdf(paste(pdf.nm, ".pdf", sep = ""), w = 3.5, h = 3.5)
        par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(3,4,1,1), col = col.pl, col.lab = col.pl, col.main = col.pl, fg = col.pl)
      }
######################################################################
    ## prevalence of infected couples
    ymax <- .35 # max(tss[,c("mm","ff","hh"),])       #
    if(is.na(cols[1]))
      {
        cols <- rainbow(length(ds.nm))
        cols <- data.frame(country = 1:length(ds.nm), cols = cols)
        cols$cols[cols$country==cframe$group.ind[js[1]]]
      }
    if(is.na(ltys))     ltys <- rep(1, length(js))
    if(is.na(lwds))     lwds <- rep(1, length(js)); lwds[1] <- 2
###################################################################### 
    ## SDC proportion
    plot(0,0, type = "n", xlim = xlim, ylim = c(0, 1), bty = "n", main = main, xlab="", ylab = ylab,
         las = 1, axes = F)
    if(!xaxt) axis(1, at = seq(early.yr, 2010, by = 5), col.axis = col.pl, las = 2, labels = NA, col.ticks = col.pl)
    if(xaxt) axis(1, at = seq(early.yr, 2010, by = 5), col.axis = col.pl, las = 2, col.ticks = col.pl)
    if(!yaxt)    axis(2, col.axis = col.pl, las = 2, labels = NA, col.ticks = col.pl)
    if(yaxt)    axis(2, col.axis = col.pl, las = 2, col.ticks = col.pl)
    if(show.all) {
      to.do <- ns:1
    }
    if(!is.na(to.show)) to.do <- to.show
    for(ii in to.do) {
      ## donor or recipient plot types to assign colors (i.e. if all substituted countries are the same then it's donor)
      cols.ind <- ifelse(length(unique(cframe[js,repl]))==1, cframe[js,'group.ind'][ii], cframe[js,repl][ii])     
      lines(tss[,"yr",ii], (tss[,"mm",ii] + tss[,"ff",ii])/tss[,"inf.alive",ii],
            col = rainbow(length(ds.nm))[ii], lty = ltys[ii], lwd = lwds[ii])
    }
    lines(tss[,"yr",1], (tss[,"mm",1] + tss[,"ff",1])/tss[,"inf.alive",1], col = 'black', lty = ltys[1], lwd = lwds[1])
    if(show.leg)               legend("bottomleft", leg = leg, col = cols$col, lwd = lwds, lty = ltys, bty = "n", cex = cex.leg, ncol = leg.col)
    if(show.sdp) {
      for(cc in 1:length(js)) {
        ## donor country simulations have 'group.ind' as c.temp
        c.temp <- ifelse(length(unique(cframe[js,repl]))==1, cframe[js[cc],'group.ind'], cframe[js[cc],repl])
        sel <- dframe.s$group == c.temp
        if(show.all.pts & c.temp !=which(ds.nm=='WA')) points(dframe.s$yr[sel], dframe.s$psdc[sel],
                                           col = rainbow(length(ds.nm))[c.temp], pch = 19, lty = 1)
      }
      if(!is.na(to.show))    {
        sel <- dframe.s$group == to.show
        points(dframe.s$yr[sel], dframe.s$psdc[sel], col = rainbow(length(ds.nm))[to.show], pch = 19, lty = 1)
      }
      points(dframe.s$yr[dframe.s$group ==cframe[js[1],'group.ind']], dframe.s$psdc[dframe.s$group==cframe[js[1],'group.ind']],
             col = 'black', pch = 19, lty = 1)
    }
    if(make.pdf)
      {
        dev.off()
        if(sep.leg)
          {
            pdf(paste(pdf.nm, "leg .pdf", sep = ""), w = 2.5, h = 3.5)
            par(mar = rep(0,4))
            plot(0,0,type = 'n', bty = 'n', axes=F, xlab='', ylab = '')
            legend("topleft", leg = leg, col = cols, lwd = lwds, lty = ltys, bty = "n")
            dev.off()
          }
      }
  }


## Plot time series array file from psrun()
plot.sdp.nsub <- function(js, # jobs to plot
                          leg, # legend for them
                          js1=1,        # which simulation is as fitted? (make black & big)
                          show.pts = T, pts.group = NA, # show DHS pts from simulation group county (dframe.s)
                          show.leg = F, cex.leg = .5, title = NA,
                          pdf.nm = 'test',
                          make.pdf = T,
                          dens = T, # histogrames or densities?
                          cols = NA,
                          show = NA,     # which lines to not plot
                          ltys = NA,
                          lwds = NA,
                          main = NA,
                          ymax2 = .35, # for density plotss
                          rmp = colorRampPalette(c("yellow","red")),
                          ylab = "serodiscordance proportion",
                          early.yr = 1985,
                          col.pl = 'black',
                          sep.leg = T,
                          yaxt = T,
                          xaxt = T,
                          browse = F)
  {
    if(browse) browser()
    tss <- cfs$t.arr[,,js]
    if(length(js)>1)
      {
        tss <- tss[tss[,"yr",1]>early.yr,,]
        ns <- dim(tss)[3]
        xlim <- range(tss[,"yr",1], na.rm=T)
      }else{
        tss <- tss[tss[,"yr"]>early.yr,]
        ns <- 1
        xlim <- range(tss[,"yr"], na.rm=T)
      }
    xlim[2] <- 2012
    if(make.pdf)
      {
        pdf(paste(pdf.nm, ".pdf", sep = ""), w = 3.5, h = 3.5)
        par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(3,4,1,1), col = col.pl, col.lab = col.pl, col.main = col.pl, fg = col.pl)
      }
######################################################################
    ## prevalence of infected couples
    ymax <- .35 # max(tss[,c("mm","ff","hh"),])       #
    if(is.na(cols[1]))      cols <- rmp(length(js))
    cols[-show] <- NA
    cols[cframe$job[js]==js1] <- 'black'
    ## print(cols)
    if(is.na(ltys[1]))     ltys <- rep(1, length(js))
    if(is.na(lwds[1]))     lwds <- rep(1, length(js)); lwds[cframe$job[js]==js1] <- 2
###################################################################### 
    ## SDC proportion
    plot(0,0, type = "n", xlim = xlim, ylim = c(0, 1), bty = "n", main = main, xlab="", ylab = ylab,
         las = 1, axes = F)
    if(xaxt) axis(1, col.axis = col.pl, las = 2, col.ticks = col.pl)
    if(!xaxt) axis(1, labels = NA, col.axis = col.pl, las = 2, col.ticks = col.pl)
    if(!yaxt)    axis(2, at = seq(0,1, l = 5), col.axis = col.pl, las = 2, labels = NA, col.ticks = col.pl)
    if(yaxt)    axis(2, at = seq(0,1, l = 5),  col.axis = col.pl, las = 2, col.ticks = col.pl)
    if(show.pts) {
            points(dframe.s$yr[dframe.s$group==pts.group], dframe.s$psdc[dframe.s$group==pts.group], pch = 19, col = 'black', cex = 1.1)
               }
    for(ii in 1:ns)
      {
        if(ns>1) lines(tss[,"yr",ii], (tss[,"mm",ii] + tss[,"ff",ii])/tss[,"inf.alive",ii], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])
        if(ns==1) lines(tss[,"yr"], (tss[,"mm"] + tss[,"ff"])/tss[,"inf.alive"], col = cols[ii], lty = ltys[ii], lwd = lwds[ii])        
      }
    if(show.leg)               legend("bottomleft", leg = leg, col = cols, lwd = 2, lty = ltys, bty = "n", cex = cex.leg, ncol = 2, title = title)
    if(make.pdf)
      {dev.off()
       if(sep.leg)
         {
           pdf(paste(pdf.nm, "leg .pdf", sep = ""), w = 2.5, h = 3.5)
           par(mar = rep(0,4), col = col.pl, col.lab = col.pl, col.main = col.pl, fg = col.pl)
           plot(0,0,type = 'n', bty = 'n', axes=F, xlab='', ylab = '')
           legend("topleft", leg = leg, col = cols, lwd = lwds, lty = ltys, bty = "n")
           dev.off()
         }
     }
  }

## dplot(35, pdf.nm = 'het dot.pdf', browse=F)

## dot plot at 2010
dplot <- function(js,
                  pdf.nm = 'test.pdf',
                  cex = .8,
                  cex.ax = .7,
                  maxnum = 400,
                  xlim = c(1975, 2012),
                  ylim = c(0, 35),
                  ylab = "years sexually active before couple formation",
                  xlab = "couple formation",
                  mains = c("M+F+","M+F-","M-F+","M-F-"),
                  cols = c('grey','red','blue',NA),
                  browse = F)
  {
    if(browse) browser()
    dat <- as.data.frame(e.arr[,,js])
    pdf(pdf.nm)
    par(mar = c(.5,1,1,.5), oma = c(5,5,0,9))
    cc <- dat$epic.ind[1]   #find epidemic curve for that data set
    epic.col <- "blue"
    layout(t(matrix(1:4,2,2)))
    ## mains <- c("B","A","D","C")
    sel <-  dat$alive==1 & dat$mage < 60*12 & dat$fage < 50*12 & dat$ser != 4
    dat <- dat[sel,]
    dat <- dat[sample(1:nrow(dat), maxnum),]
    for(ii in 2:1)
      {
        plot(1900 + 1/12*(dat$tint-dat$mardur.mon)[dat$ser==ii],
             1/12*(dat$tmar-dat$tms)[dat$ser==ii ],
             col = cols[dat$mcoi[dat$ser==ii ]], xlab = "",
             pch = 19, cex = cex, axes = F, ylab = ylab,
             xlim = xlim, ylim = ylim, bty = "n",
             main = mains[ii])
        xs <- epicm[,1]/12 + 1900
        if(ii==2) axis(2, at = seq(0, 30, by = 10), las = 2)
        if(ii==1) axis(2, at = seq(0, 30, by = 10), labels = NA, las = 2)
        axis(1, at = seq(1980, 2010, by = 10), labels = NA)
        lines(xs[xs>1975], epicf[xs>1975, cc]*max(ylim), col = epic.col, lwd = 1)
        ## if(ii==1) axis(4, at = seq(0, max(ylim), l=5), seq(0, 100, l = 5), col = epic.col, las = 2)
      }
    for(ii in c(3,1))
      {
        plot(1900 + 1/12*(dat$tint-dat$mardur.mon)[dat$ser==ii ],
             1/12*(dat$tmar-dat$tfs)[dat$ser==ii ],
             col = cols[dat$fcoi[dat$ser==ii ]], pch = 19, las = 2, axes = F,
             xlim = xlim, ylim = ylim, bty = "n", xlab = "", cex = cex, ylab = ylab,
             main = mains[ii])
        if(ii==3) axis(2, at = seq(0, 30, by = 10), las = 2)
        if(ii==1) axis(2, at = seq(0, 30, by = 10), labels = NA, las = 2)
        axis(1, at = seq(1980, 2010, by = 10), las = 2)
        lines(xs[xs>1975], epicm[xs>1975, cc]*max(ylim), col = epic.col, lwd = 1)
        ## if(ii==1) axis(4, at = seq(0, max(ylim), l=5), seq(0, 100, l = 5), col = epic.col, las = 2)
      }
        mtext("male YSA before couple formation", side = 2, outer = T, line = 2, cex = cex.ax, adj = .82)
    mtext("female YSA before couple formation", side = 2, outer = T, line = 2, cex = cex.ax, adj = .08)
    dev.off()
  }

coitrsex <- function(js, # jobs to plot
                  leg, # legend for them
                  pdf.nm = 'test',
                  cols = NA,
                  ltys = NA,
                  lwds = NA,
                  rmp = colorRampPalette(c("blue","orange")),
                  early.yr = 1985, 
                  browse = F)
  {
    if(browse) browser()
    tss <- cfs$t.arr[,,js]
    tss <- tss[tss[,"yr",1]>early.yr,,]
    pdf(paste(pdf.nm, ".pdf", sep = ""), w = 7, h = 7)
    par(mfrow=c(2,1), oma = c(0,2,2,0), mar = c(5,5,2,1))
######################################################################
    ymax <- max(tss[,c('cu.mb','cu.me','cu.mp','cu.fb','cu.fe','cu.fp'),])
    if(is.na(cols))      cols <- rmp(length(js)); cols[js==1] <- 'green'    
    if(is.na(ltys))     ltys <- rep(1, length(js))
    if(is.na(lwds))     lwds <- rep(1, length(js))
    ## cumulative premarital infections    
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "Male", xlab="", ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.mb",ii], col = cols[ii], lty = 1, lwd = lwds[ii])
        lines(tss[,"yr",ii], tss[,"cu.me",ii], col = cols[ii], lty = 2, lwd = lwds[ii])
        lines(tss[,"yr",ii], tss[,"cu.mp",ii], col = cols[ii], lty = 3, lwd = lwds[ii])        
      }
    legend("topleft", leg = c('pre-','extra-','within'), col = c(rep('black',3)), lwd = c(rep(1,3)), lty = c(1:3), bty = "n")
    plot(0,0, type = "n", xlim = range(tss[,"yr",1]), ylim = c(0, ymax*1.2), las = 2, bty = "n", main = "Female", xlab="", ylab = "")
    for(ii in 1:dim(tss)[3])     
      {
        lines(tss[,"yr",ii], tss[,"cu.mb",ii], col = cols[ii], lty = 1, lwd = lwds[ii])
        lines(tss[,"yr",ii], tss[,"cu.me",ii], col = cols[ii], lty = 2, lwd = lwds[ii])
        lines(tss[,"yr",ii], tss[,"cu.mp",ii], col = cols[ii], lty = 3, lwd = lwds[ii])        
      }
    legend("topleft", leg = leg, col = cols, lwd = lwds, lty = ltys, bty = "n")    
    mtext('cumulative infections', side = 2, line = .5, outer=T, adj=.5)
    dev.off()
  }
