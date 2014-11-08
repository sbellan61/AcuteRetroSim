## These functions take the output of couples simulations and create a retrospective cohort of
## couples that were observed as serodiscordant to simulate the analyses of Wawer et al. 2005 and
## Hollingsworth et al. 2008 performed to estimate the relative infectivity of the acute phase to
## the chronic phase.

## Calculate hazard ratios for acute, chronic, late & aids phases from data
empir.arh <- function(dat, ts, start.rak = 1994, end.rak = 1999, browse = F)
  {
    if(browse) browser()
################################################## 
    ## find all couples for which they were serodiscordant for at least one month in the [start, end) interval
    sttmon <- (start.rak-1900)*12       # start month in CMC
    endmon <- (end.rak-1900)*12         # end month in CMC
    ## Select couples that were serodiscordant during this time period.
    sel <- apply(ts, 2, function(x) {sum(x[sttmon:(endmon -1)] %in% c('mm','ff')) > 0 })
    ts <- ts[,sel]   ## subset cohort frame
    dat <- dat[sel,] ## subset couple line list    
    si <- dat$ser %in% 1:3 ## couples in which someone is infected (redundant with sel above but don't feel like changing below code yet)
    idoi <- rep(NA, nrow(dat))       ## date of index partner infection
    idoi[si] <- apply(dat[si,c('mdoi','fdoi')], 1, function(x) {min(x,na.rm=T)})
##################################################
    ## Now, add stuff to the cohort frame that identifies what phase each couple is in if they are serodiscordant.
    ## 1st month of acute phase: append 'ac' to 'ff'/'mm' whenever they're in acute phase if they're SDC still ff/mm
    ## (i.e. haven't become hh), even if before marriage:
    ## is.sdc: are they SDC at date of index infection?
    is.sdc <- grepl('mm', ts[cbind(idoi[si], which(si))]) | grepl('ff', ts[cbind(idoi[si], which(si))])
    si.temp <- which(si)[is.sdc] ## indices of those that are SDC @ index infection date
    ts[cbind(idoi[si.temp], si.temp)] <- paste(ts[cbind(idoi[si.temp], si.temp)], '.ac', sep = '') # first month, will always still be SDC
    ## 2nd month of acute phase: append 'ac' to 'ff'/'mm' whenever they're in acute phase if they're SDC still (ff/mm, even if before marriage
    si.temp <- si & !(idoi==nrow(ts))       # get rid of infections at the end of the simulation with no 2nd month
    is.sdc <- grepl('mm', ts[cbind(idoi[si.temp]+1, which(si.temp))]) | grepl('ff', ts[cbind(idoi[si.temp]+1, which(si.temp))])
    si.temp <- which(si.temp)[is.sdc]
    ts[cbind(idoi[si.temp]+1, si.temp)] <- paste(ts[cbind(idoi[si.temp]+1, si.temp)], '.ac', sep = '')
    ## date of index partner's death
    idod <- rep(NA, nrow(dat))
    idod[si] <- apply(dat[si, c('mdod','fdod')], 1, function(x) {min(x, na.rm = T)})
    for(mm in -19:-10) {## find late stage SDC's
      si.temp <- si & !(idod + mm > nrow(ts)) # who was late infected before sim ended
      is.sdc <- grepl('mm', ts[cbind(idod[si.temp]+mm, which(si.temp))]) | grepl('ff', ts[cbind(idod[si.temp]+mm, which(si.temp))])    
      si.temp <- which(si.temp)[is.sdc]
      ts[cbind(idod[si.temp]+mm, si.temp)] <- paste(ts[cbind(idod[si.temp]+mm, si.temp)], '.lt', sep = '')
    }
    for(mm in -9:-1) {## find aids stage SDC's
      si.temp <- si & !(idod + mm > nrow(ts)) # who was late infected before sim ended
      is.sdc <- grepl('mm', ts[cbind(idod[si.temp]+mm, which(si.temp))]) | grepl('ff', ts[cbind(idod[si.temp]+mm, which(si.temp))])    
      si.temp <- which(si.temp)[is.sdc]
      ts[cbind(idod[si.temp]+mm, si.temp)] <- paste(ts[cbind(idod[si.temp]+mm, si.temp)], '.ad', sep = '')
    }
    ## rest are chronic
    ts[which(ts=='mm',arr.ind=T)] <- 'mm.ch'
    ts[which(ts=='ff',arr.ind=T)] <- 'ff.ch'
##################################################
    ## tally up person-months of exposure by phase for each couple
    pm.ac <- apply(ts, 2, function(x) {sum(x[sttmon:(endmon -1)] %in% c('mm.ac','ff.ac'))})
    pm.ch <- apply(ts, 2, function(x) {sum(x[sttmon:(endmon -1)] %in% c('mm.ch','ff.ch'))})
    pm.lt <- apply(ts, 2, function(x) {sum(x[sttmon:(endmon -1)] %in% c('mm.lt','ff.lt'))})
    pm.ad <- apply(ts, 2, function(x) {sum(x[sttmon:(endmon -1)] %in% c('mm.ad','ff.ad'))})
    pm.sdc <- pm.ac + pm.ch + pm.lt + pm.ad
    pm.nac <- pm.sdc - pm.ac                # all phases but acute *n*ot *ac*ute
    ## for each couple identify gender of partner at risk
    gend <- rep(NA,nrow(dat))
    gend[apply(ts, 2, function(x) {sum(grepl('mm', x[sttmon:(endmon -1)])) > 0})] <- 'f'
    gend[apply(ts, 2, function(x) {sum(grepl('ff', x[sttmon:(endmon -1)])) > 0})] <- 'm'
    ## Person Months at Risk
    gens <- c('m','f')              ## genders
    phs <- c('sdc','ac','ch','lt','ad','nac') ## phases (all sdcs, acute, chronic, late, aids, not acute)
    for(gen in gens) { ## calculate person months of exposure by gender exposed & phase of their partner
      for(pp in phs)  {
        assign(paste0('sel.',pp,'.',gen), get(paste0('pm.',pp))>0 & gend==gen) # create logical seletors
        temp <- get(paste0('pm.',pp))
        temp[gend!=gen] <- 0
        assign(paste0('pm.',pp,'.',gen), temp)
      }
    }
##################################################
    ## initialize vector of names: how many couples by exposure category
    num.couples <- sapply(paste0('sel.',phs, '.', rep(gens, each = length(phs))),function(x) {sum(get(x))})
    ## initialize vector of names: how many person-months by exposure category
    num.pms <- sapply(paste0('pm.',phs, '.', rep(gens, each = length(phs))),function(x) {sum(get(x))}) 
    for(gen in gens) { ## for each gender and phase, calculate number of infections
      ## acute
      tsel <- get(paste0('sel.ac.',gen))
      temp.inf <- which(grepl('hh', ts[endmon,tsel]) & (dat$mcoi.phase[tsel] == 'a' | dat$fcoi.phase[tsel] == 'a'))
      assign(paste0('ac.inf.',gen), which(tsel)[temp.inf])
      ## chronic
      tsel <- get(paste0('sel.ch.',gen))
      temp.inf <- which(grepl('hh', ts[endmon,tsel]) & (dat$mcoi.phase[tsel] == 'c' | dat$fcoi.phase[tsel] == 'c'))
      assign(paste0('ch.inf.',gen), which(tsel)[temp.inf])
      ## late
      tsel <- get(paste0('sel.lt.',gen))
      temp.inf <- which(grepl('hh', ts[endmon,tsel]) & (dat$mcoi.phase[tsel] == 'l' | dat$fcoi.phase[tsel] == 'l'))
      assign(paste0('lt.inf.',gen), which(tsel)[temp.inf])
      ## aids
      tsel <- get(paste0('sel.ad.',gen))
      temp.inf <- which(grepl('hh', ts[endmon,tsel]) & (dat$mcoi.phase[tsel] == 'ad' | dat$fcoi.phase[tsel] == 'ad'))
      assign(paste0('ad.inf.',gen), which(tsel)[temp.inf])
      ## anything after acute
      tsel <- get(paste0('sel.nac.',gen))
      temp.inf <- which(grepl('hh', ts[endmon,tsel]) & (dat$mcoi.phase[tsel] %in% c('c','l','ad') | dat$fcoi.phase[tsel] %in% c('c','l','ad')))
      assign(paste0('nac.inf.',gen), which(tsel)[temp.inf])
      ## anything
      tsel <- get(paste0('sel.sdc.',gen))
      temp.inf <- which(grepl('hh', ts[endmon,tsel]) & (dat$mcoi.phase[tsel] %in% c('a','c','l','ad') | dat$fcoi.phase[tsel] %in% c('a','c','l','ad')))
      assign(paste0('sdc.inf.',gen), which(tsel)[temp.inf])
    }
    ## calculate raw hazards per 100 person years stratified by gender & phase
    for(gen in gens) {
      temp.hzs <- rep(NA,length(phs))
      names(temp.hzs) <- phs
      for(ph in phs) {
        temp.inf <- length(get(paste0(ph,'.inf.',gen)))
        temp.pm <- sum(get(paste0('pm.',ph,'.',gen)))
        temp.haz <- temp.inf / temp.pm ## infections/person-months
        temp.hzs[ph] <- -log(1-temp.haz) ## adjustment from risk to hazard: prob = 1 - exp(-rate*time); rate = -log(1-prob)/time
      }
      assign(paste0('hzs.',gen), temp.hzs)
    }
    ## genders combined
    hzs.both <- rep(NA,length(phs))
    names(hzs.both) <- phs
    for(ph in phs) {
      temp.inf <- length(get(paste0(ph,'.inf.',gens[1]))) + length(get(paste0(ph,'.inf.',gens[2])))
      temp.pm <- sum(get(paste0('pm.',ph,'.',gens[1]))) + sum(get(paste0('pm.',ph,'.',gens[2])))
      temp.haz <- temp.inf / temp.pm ## infections/person-months
      hzs.both[ph] <- -log(1-temp.haz) ## adjustment from risk to hazard: prob = 1 - exp(rate*time); rate = log(1-prob)/time
    }
    hzs <- cbind(hzs.both, hzs.m, hzs.f)
    colnames(hzs) <- c('all','m','f')
    rownames(hzs) <- phs
    hrs <- hzs / hzs[rep(which(phs=='ch'),length(phs)),]
    return(list(hzs, hrs))
  }


## Calculate hazards from a ts type object for observed person months 
truehrs <- function(tst, evt, pars, sl=NULL) {
  ##  browser()
  if(!is.null(sl)) {
    tst <- tst[,sl] #cohsim$ts.rak.all
    evt <- evt[sl,] ##cohsim$dat.rak
  }
  pms.ac <- sum(grepl('ac', as.vector(tst))) ## person-months acute
  infs.ac <- sum(apply(tst,2,function(x) sum(grepl('hh',x))>0) & (evt$mcoi.phase=='a' | evt$fcoi.phase=='a'), na.rm=T)
  pms.ch <- sum(as.vector(tst) %in% c('mm','ff'))
  infs.ch <- sum(apply(tst,2,function(x) sum(grepl('hh',x))>0) & (evt$mcoi.phase=='c' | evt$fcoi.phase=='c'), na.rm=T)
  pms.lt <- sum(grepl('lt', as.vector(tst)))
  infs.lt <- sum(apply(tst,2,function(x) sum(grepl('hh',x))>0) & (evt$mcoi.phase=='l' | evt$fcoi.phase=='l'), na.rm=T)
  pms.aids <- sum(grepl('aids', as.vector(tst)))
  infs.aids <- sum(apply(tst,2,function(x) sum(grepl('hh',x))>0) & (evt$mcoi.phase=='ad' | evt$fcoi.phase=='ad'), na.rm=T)
  for(ph in c('ac','ch','lt','aids')) assign(paste0('hz.',ph),get(paste0('infs.',ph))/get(paste0('pms.',ph)))
  hzs <- c(ac=hz.ac,ch=hz.ch, lt=hz.lt)
  hrs <- hzs/hz.ch
  ehms <- c(ac = as.numeric((hz.ac/hz.ch-1)*pars['dur.ac']),
            lt = as.numeric((hz.lt/hz.ch-1)*pars['dur.lt']),
            ltaids = as.numeric((hz.lt/hz.ch-1)*pars['dur.lt']  +   (hz.aids/hz.ch-1)*pars['dur.aids']))
  rm(list=setdiff(ls(), c("hzs","hrs","ehms"))) ## remove everything but output
  out <- lapply(list(hzs=hzs, hrs = hrs, ehms = ehms), function(x) signif(x, 3))
  return(out)
  gc()
}

sdcs <- c('mm','mm.ac','mm.lt','mm.aids','ff','ff.ac','ff.lt','ff.aids')
sers.ap <- list(ss='ss', mm =c('mm','mm.ac','mm.lt','mm.aids'), ff = c('ff','ff.ac','ff.lt','ff.aids'), hh = 'hh')
sers <- list('ss','mm','ff','hh')

## create a retrospective cohort at interv monthly intervals from start.rak to end.rak, only keep
## couples that were observed more than once, and which were serodiscordant at some point during the
## observed time (they could go -- to ++ in one visit interval though).
rak.coh.fxn <- function(output, interv = 10, max.vis = 5, start.rak=1994, end.rak=1999.5, ## interv is interval in months between visits
                        ## ltf.prob = monthly probability of loss to follow-up, rr = +- or -+ vs -- or ++
                        ltf.prob = NA, rr.ltf.ff = 1, rr.ltf.mm = 1, rr.ltf.hh = 1, rr.ltf.d = 0, rr.inc.sdc = 1, # .d is ltf when dead
                        verbose = F, browse = F) 
  {
      ts.ap <- output$ts
      dat <- output$evout
      dpars <- output$rakpars
    if(browse) browser()
    if(verbose)     {
      print('from full simulation:')
      sel <- apply(ts.ap, 2, function(x) { sum(x %in% sdcs)>0 } )
      print(truehrs(ts.ap[,sel], dat[sel,], pars = dpars))
    }
    sttmon <- (start.rak-1900)*12       # start month in CMC
    endmon <- (end.rak-1900)*12         # end month in CMC
    ## Select all couples that were serodiscordant at some point during the cohort (don't have to worry about ss
    ## -> hh transition in one month, because in our model an individual cannot be infected & infect their
    ## partner in the same month)
    sel <- apply(ts.ap[sttmon:endmon,], 2, function(x) { sum(x %in% sdcs)>0 } )
    ts.sdc <- ts.ap[,sel]
    dat.sdc <- dat[sel,]
    ## reduce data frame to visit months
    vis.mon <- seq(sttmon, endmon, by = interv)
    vis.mon.all <- min(vis.mon):max(vis.mon)
    ts.vm <- ts.sdc[vis.mon,]
    ts.vm.all <- ts.sdc[vis.mon.all,]
    row.names(ts.vm) <- vis.mon
    row.names(ts.vm.all) <- vis.mon.all
    ## replace time points of couples with partners aged out of cohort with NA (>49 for f, >59 for m; note actual Rakai
    ## cohort is >59 for both, but we stick with DHS criteria, shouldn't make a difference)
    ts.vm <- apply(ts.vm, 2, function(x) { x[grepl('a\\.',x)] <- NA; return(x) } )
    ts.vm.all <- apply(ts.vm.all, 2, function(x) { x[grepl('a\\.',x)] <- NA; return(x) } )        
    ## replace pre-couple time points with NAs
    ts.vm <- apply(ts.vm, 2, function(x) { x[grepl('b\\.',x)] <- NA; return(x) } )
    ts.vm.all <- apply(ts.vm.all, 2, function(x) { x[grepl('b\\.',x)] <- NA; return(x) } )
    ## Remove any couples with only NA's now
    no.obs <- apply(ts.vm,2, function(x) sum(!is.na(x)))==0
    dat.vm <- dat.sdc[!no.obs,]
    ts.vm <- ts.vm[,!no.obs]
    ts.vm.all <- ts.vm.all[,!no.obs]
    if(verbose)     {
      print(paste('after subsetting to SDCs from', start.rak, 'to', end.rak))
      print(truehrs(ts.vm.all, dat.vm, pars = dpars))
    }
    ## Censor data for couples visited more than 5 times, (only 5 rounds were done in Rakai study, 40 months total)
    num.vis <- apply(ts.vm, 2, function(x) sum(!is.na(x)))
    ## number of visits each couple makes before being lost to follow-up
    if(!is.na(ltf.prob)) {
      ltf.rt <- -log(1-ltf.prob) ## monthly rate         
      ltfp.ss <- 1 - exp(-ltf.rt * interv) ## probability of loss to follow up in each interval (ss)
      ltfp.mm <- 1 - exp(-rr.ltf.mm * ltf.rt * interv) ## probability of loss to follow up in each interval (mm)
      ltfp.ff <- 1 - exp(-rr.ltf.ff * ltf.rt * interv) ## probability of loss to follow up in each interval (ff)
      ltfp.hh <- 1 - exp(-rr.ltf.hh * ltf.rt * interv) ## probability of loss to follow up in each interval (hh)
      ltfp.d <- 1 - exp(-rr.ltf.d * ltf.rt * interv) ## probability of loss to follow up in each interval (hh)
      ## Identify incident serodiscordant visits as high risk for LTF
      inc.sdc.wh <- which(apply(ts.vm, 2, function(x) sum(x[-nrow(ts.vm)]=='ss' & x[-1] %in% sdcs))==1)
      inc.sdc.wh.tt <- apply(ts.vm[,inc.sdc.wh], 2, function(x) which(x[-nrow(ts.vm)]=='ss' & x[-1] %in% sdcs)) + 1
      ## ts.vm[,head(inc.sdc.wh)] ## couples that were seroincident
      ## head(inc.sdc.wh.tt)      ## which visit were they first serodiscordant?
      ## create ltfp matrix
      ltfps <- ts.vm
      ## Set serostatus dependent LTFps, we are marking the probability that a couple is gone at the NEXT visit (so first visits are always seen)
      for(ser in sers)        ltfps[which(ts.vm%in%sers.ap[[ser]], arr.ind=T)] <- get(paste0('ltfp.',ser))
      ## Set death dependent LTFps
      ltfps[apply(ltfps, 2, function(x) grepl('d\\.', x))] <- ltfp.d
      ## Increase relative rate of first seroconversion LTFps
      inc.sdc.temp.ltfp <- as.numeric(ltfps[cbind(inc.sdc.wh.tt, inc.sdc.wh)]) ## extract these probabilities
      inc.sdc.temp.ltfr <- -log(1-inc.sdc.temp.ltfp)/interv * rr.inc.sdc ## convert to rates & multiply by relative rate
      inc.sdc.temp.ltfp <- 1-exp(-inc.sdc.temp.ltfr * interv) ## convert back to probailities
      ltfps[cbind(inc.sdc.wh.tt, inc.sdc.wh)] <- inc.sdc.temp.ltfp ## replace back in matrix logical
      ## matrix of ltf censoring where serostatus at each visit affects the probability of being
      ## observed at the next. This line determines the probability of loss to follow-up after any
      ## visit, the next lines then find the earliest one that was lost.
      ltf.log <- apply(ltfps,2, function(x) { x <- as.numeric(x) ## currently a character matrix
                                              x[!is.na(x)] <- rbinom(sum(!is.na(x)), 1, x[!is.na(x)]) ## only do it for !NA
                                              return(x) })
      rand <- sample(1:ncol(ts.vm),8)
      ltfps[,rand]
      ts.vm[,rand]
      ltf.log[,rand]
      ## first censorship
      ltf.vis <- rep(NA,ncol(ltf.log))
      ltf.vis1 <- apply(ltf.log, 2, function(x) sum(x, na.rm=T)) > 0
      ltf.vis[ltf.vis1] <- apply(ltf.log[,ltf.vis1], 2, function(x) min(which(x==1),na.rm=T))
    }else{ ## no loss to follow-up
      ltf.vis <- NA
    }
    for(ii in 1:length(num.vis)) {  # for each couple, determine censorship
      if(!is.na(ltf.vis[ii])) { ## if lost to follow-up ever
        ts.vm[1:nrow(ts.vm) > ltf.vis[ii], ii] <-  NA ## censor all visits *AFTER* loss to follow up
      }
      temp.num.vis <- sum(!is.na(ts.vm[,ii])) ## how many observations are left?
      if(temp.num.vis>max.vis) { ## if more observations than max observations
        temp.visits <- which(!is.na(ts.vm[,ii])) # which visits were observed         
        show.visits <- temp.visits[1:max.vis] # visits to show
        ts.vm[!1:nrow(ts.vm) %in% show.visits, ii] <- NA # censor others
      }
      ## Remove person-time months from ts.vm.all for censored time
      vis.obs <- which(!is.na(ts.vm[,ii]))
      obs.range <- range(rownames(ts.vm)[vis.obs])
      unseen <- rownames(ts.vm.all) < obs.range[1] | rownames(ts.vm.all) > obs.range[2]
      ts.vm.all[unseen,ii]<- NA
    }
    ## now, again reduce cohort data frame to all couples that were observed serodiscordant at some
    ## point during the cohort, or that were observed -- and then ++ between visit intervals.
    sel.coh <- apply(ts.vm, 2, function(x) { sum(x %in% sdcs, na.rm=T)>0 | (sum(x=='ss', na.rm=T)>0 & sum(x=='hh', na.rm=T)>0) } )
    ts.vm <- ts.vm[, sel.coh]
    ts.vm.all <- ts.vm.all[, sel.coh]        
    dat.vm <- dat.vm[sel.coh,]
    if(verbose) {
      print('EHMs after censoring due to LTF or max.vis:')
      print(truehrs(ts.vm.all, dat.vm, pars = dpars))
    }
    ## remove any couples that were only observed once because of aging out
    nage <- colSums(!is.na(ts.vm))>1
    ts.vm <- ts.vm[,nage]
    ts.vm.all <- ts.vm.all[,nage]        
    dat.vm <- dat.vm[nage,]
    ## remove any couples that were only observed once because of death (is.na clause makes sure to exclude pre-couple times)
    death1 <- apply(ts.vm, 2, function(x) { sum(!grepl('d\\.',x) & !is.na(x)) ==1 })
    ts.vm <- ts.vm[,!death1]
    ts.vm.all <- ts.vm.all[,!death1]        
    dat.vm <- dat.vm[!death1,]
    ## Find CMC month of first survey visit to add to dat.vm
    dat.vm$fvis.tt <- as.numeric(rownames(ts.vm)[apply(ts.vm, 2, function(x) min(which(x %in% c(sdcs,'hh'))))])
    ## Return results
    rak.coh <- list(dat.rak = dat.vm, ts.rak = ts.vm, ts.rak.all = ts.vm.all, interv = interv, dpars = dpars)
    rm(list=setdiff(ls(), c('rak.coh'))) ## remove everything but output
    gc() # clear memory
    return(rak.coh)
  }

## Calculate person-months at risk for second partner in each couple group.
## Used in rak.wawer below
make.rakll <- function(dat.vm, ts.vm, cov.mods=F, interv=10, verbose2=F) { 
    rakll <- data.frame(uid = dat.vm$uid, phase = NA, pm = NA, inf = 0, pm.trunc = NA, inf.trunc = NA, excl.by.err = F,
                        mcoi = dat.vm$mcoi, fcoi = dat.vm$fcoi, mcoi.phase = dat.vm$mcoi.phase, fcoi.phase = dat.vm$fcoi.phase,
                        secp = NA, secp.lhet = NA, secp.age = NA, indp.age = NA, mardur = NA,
                        secp.tdsa = NA, secp.pdsa = NA, indp.duri = NA, secp.hazm = NA, secp.thazm = NA) ## second partner infected
    secp.m <- which(dat.vm$mdoi>dat.vm$fdoi | is.na(dat.vm$mdoi)) ## second partner is male is infected after her or never gets infected
    secp.f <- which(dat.vm$mdoi<dat.vm$fdoi | is.na(dat.vm$fdoi)) ##
    rakll$secp[secp.m] <- 'm'
    rakll$secp[secp.f] <- 'f'
    rakll$secp.lhet[secp.m] <- log(dat.vm$m.het.gen[secp.m])
    rakll$secp.lhet[secp.f] <- log(dat.vm$f.het.gen[secp.f])
    if(cov.mods)  { ## if we are going to run multivariate regression models
        ## age of secondary partner at first interval followed in cohort study
        rakll$secp.age[secp.m] <- with(dat.vm, mage[secp.m] - (tint[secp.m] - fvis.tt[secp.m]))
        rakll$secp.age[secp.f] <- with(dat.vm, fage[secp.f] - (tint[secp.f] - fvis.tt[secp.f]))
        ## age of index partner at first interval followed in cohort study
        rakll$indp.age[secp.m] <- with(dat.vm, fage[secp.m] - (tint[secp.m] - fvis.tt[secp.m]))
        rakll$indp.age[secp.f] <- with(dat.vm, mage[secp.f] - (tint[secp.f] - fvis.tt[secp.f]))
        ## partnership duration at first interval followed up
        rakll$mardur <- with(dat.vm, mardur.mon - (tint-fvis.tt))
        ## total sexual activity time of secondary partner at first interval followed in cohort study
        rakll$secp.tdsa[secp.m] <- with(dat.vm, fvis.tt[secp.m] - (tms[secp.m]))
        rakll$secp.tdsa[secp.f] <- with(dat.vm, fvis.tt[secp.f] - (tfs[secp.f]))    
        ## pre-couple sexual activity time of secondary partner
        rakll$secp.pdsa[secp.m] <- with(dat.vm, tmar[secp.m] - (tms[secp.m]))
        rakll$secp.pdsa[secp.f] <- with(dat.vm, tmar[secp.f] - (tfs[secp.f]))
        ## duration of primary partner's infection
        rakll$indp.duri[secp.m] <- with(dat.vm, fvis.tt[secp.m] - fdoi[secp.m])
        rakll$indp.duri[secp.f] <- with(dat.vm, fvis.tt[secp.f] - mdoi[secp.f])
        ## ############################
        ## hazard months 2nd partner has already been exposed to by first visit
        rakll$secp.hazm[secp.m] <- with(dat.vm[secp.m,], fvis.tt - apply(cbind(tmar, fdoi),1, max))
        rakll$secp.hazm[secp.f] <- with(dat.vm[secp.f,], fvis.tt - apply(cbind(tmar, mdoi),1, max))
        ## add excess hazard months due to acute phase: those exposed to both months
        secp.m.ac2 <- with(dat.vm[secp.m,], fdoi-tmar) >= -1
        secp.f.ac2 <- with(dat.vm[secp.f,], mdoi-tmar) >= -1
        rakll$secp.hazm[secp.m[secp.m.ac2]] <- rakll$secp.hazm[secp.m[secp.m.ac2]] + 2*(dpars['acute.sc']-1)
        rakll$secp.hazm[secp.f[secp.f.ac2]] <- rakll$secp.hazm[secp.f[secp.f.ac2]] + 2*(dpars['acute.sc']-1)
        ## those exposed to one month
        secp.m.ac1 <- with(dat.vm[secp.m,], fdoi-tmar) == -2
        secp.f.ac1 <- with(dat.vm[secp.f,], mdoi-tmar) == -2
        rakll$secp.hazm[secp.m[secp.m.ac1]] <- rakll$secp.hazm[secp.m[secp.m.ac1]] + 1*(dpars['acute.sc']-1)
        rakll$secp.hazm[secp.f[secp.f.ac1]] <- rakll$secp.hazm[secp.f[secp.f.ac1]] + 1*(dpars['acute.sc']-1)
        ## Add total hazard exposed to by first visit, this equals equivalent hazard-months exposed to
        ## in chronic phase plus those pre- & extra-couple
        rakll$secp.thazm <- rakll$secp.hazm
        rakll$secp.thazm[secp.m] <- spars['bmp'] * rakll$secp.thazm[secp.m] + 
            spars['bmb'] * apply(dat.vm[secp.m,c('tms','tmar')], 1, function(x) sum(epicf[x['tms']:(x['tmar']-1),epic.ind])) +
                spars['bme'] * apply(dat.vm[secp.m,c('fvis.tt','tmar')], 1, function(x) sum(epicf[x['tmar']:(x['fvis.tt']-1),epic.ind]))
        rakll$secp.thazm[secp.f] <- spars['bfp'] * rakll$secp.thazm[secp.f] + 
            spars['bfb'] * apply(dat.vm[secp.f,c('tfs','tmar')], 1, function(x)     sum(epicm[x['tfs']: (x['tmar']-1),epic.ind])) +
                spars['bfe'] * apply(dat.vm[secp.f,c('fvis.tt','tmar')], 1, function(x) sum(epicm[x['tmar']:(x['fvis.tt']-1),epic.ind]))
    }
    ## Excluded by error are couples that were observed during at least 2 visits, but with the
    ## last observed visit being serodiscordant. Based on Wawer et al.'s methods description and
    ## the fact that the exact same # of incident couples were followed for 1 interval as for 2
    ## (indicating that they excluded any just followed for 1).      
##################################################
    ## Incident infections
    inc.wh <- which(apply(ts.vm, 2, function(x) sum(grepl('ss',x))>0))
    rakll$phase[inc.wh] <- 'inc'
    last.sus <- rep(NA, ncol(ts.vm))
    last.sus[inc.wh] <- apply(ts.vm[,inc.wh,drop=F], 2, function(x) max(which(x=='ss')))
    ## those that went to ++ at some point
    inc.wh.hh <- inc.wh[which(apply(ts.vm[,inc.wh,drop=F], 2, function(x) sum(grepl('hh',x))>0))]
    rakll$inf[inc.wh.hh] <- 1
    ## those that went -- to ++ in one interval
    inc.wh.hh1 <- inc.wh[ts.vm[cbind(last.sus[inc.wh]+1,inc.wh)]=='hh'] 
    rakll$pm[inc.wh.hh1] <- interv/4
    ## those that went -- to ++ but not in one interval
    inc.wh.hh2 <- inc.wh.hh[!inc.wh.hh %in% inc.wh.hh1]
    rakll$pm[inc.wh.hh2] <- apply(ts.vm[,inc.wh.hh2, drop=F], 2, function(x) sum(x %in% sdcs)) * interv
    ## those that never went to ++
    inc.wh.nhh <- inc.wh[!inc.wh %in% inc.wh.hh]
    ## those that never went to ++ that were only observed SDC once and were consequently probably excluded from the Wawer study
    inc.wh.nhh.exl.err <- inc.wh.nhh[apply(ts.vm[,inc.wh.nhh,drop=F], 2, function(x) sum(x %in% sdcs)==1)]
    rakll$excl.by.err[inc.wh.nhh.exl.err] <- T
    ## person months = (# times observed SDC -1)*interv + interv/2
    rakll$pm[inc.wh.nhh] <- apply(ts.vm[,inc.wh.nhh,drop=F], 2, function(x) sum(x %in% sdcs)-1)*interv + interv/2
    ## ################################################
    ## Chronic infections
    ch.wh <- which(apply(ts.vm, 2, function(x) { sum(grepl('ss',x) + grepl('d\\.',x))==0  | (sum(grepl('ss',x))==0 & sum(grepl('d\\.hh\\.m',x))>0 & sum(x %in% sers.ap$ff)>0) | (sum(grepl('ss',x))==0 & sum(grepl('d\\.hh\\.f',x))>0 & sum(x %in% sers.ap$mm)>0) }))
    rakll$phase[ch.wh] <- 'prev'
    ## those that became ++
    ch.wh.hh <- ch.wh[which(apply(ts.vm[,ch.wh,drop=F], 2, function(x) sum(grepl('hh',x))>0))]
    ##person months = (# times observed SDC -1)*interv + interv/2 for ++
    rakll$pm[ch.wh.hh] <- apply(ts.vm[,ch.wh.hh,drop=F], 2, function(x) sum(x %in% sdcs)-1)*interv + interv/2
    rakll$inf[ch.wh.hh] <- 1
    ## those that stayed +-
    ch.wh.nhh <- ch.wh[!ch.wh %in% ch.wh.hh]
    rakll$pm[ch.wh.nhh] <- apply(ts.vm[,ch.wh.nhh,drop=F], 2, function(x) sum(x %in% sdcs)-1)*interv
    ## ################################################
    ## Late infections ## male death after male SDC, or vice versa
    lt.wh <- which(apply(ts.vm, 2, function(x) {sum(grepl('d\\.mm',x) | grepl('d\\.ff',x)) > 0 | (sum(grepl('d\\.hh\\.m',x))>0 & sum(x %in% sers.ap$mm)>0) | (sum(grepl('d\\.hh\\.f',x))>0 & sum(x %in% sers.ap$ff)>0) } ))
    ## **************************************************???
    ## What to do with couples that are both EARLY & LATE?? for now leave them as early only
    lt.inc.wh <- lt.wh[lt.wh %in% inc.wh]
    if(verbose2) print(paste(length(lt.inc.wh), 'couples were classified as both early & late. We keep them as early for the analysis, though if decont=T, they are completey excluded later'))
    lt.wh <- lt.wh[!lt.wh %in% lt.inc.wh]
    rakll$phase[lt.wh] <- 'late'
    ## those that became ++ and were seen ++ at a visit *including* if only first seen ++ at the first visit after a partner's death
    lt.wh.hh <- lt.wh[which(apply(ts.vm[,lt.wh,drop=F], 2, function(x) sum(grepl('hh',x) )>0))]
    rakll$inf[lt.wh.hh] <- 1
    rakll$pm[lt.wh.hh] <- apply(ts.vm[,lt.wh.hh,drop=F], 2, function(x) sum(x %in% sdcs)-1)*interv + interv/2
    ## those that stayed +- up until last visit
    lt.wh.nhh <- lt.wh[!lt.wh %in% lt.wh.hh]
    rakll$pm[lt.wh.nhh] <- apply(ts.vm[,lt.wh.nhh,drop=F], 2, function(x) sum(x %in% sdcs))*interv
    rakll$phase <- factor(rakll$phase)
    rakll$phase <- relevel(rakll$phase, ref = 'prev')
    ## ################################################
    ## one final adjustment, Wawer compare hazards from first '5 months' post incident couples'
    ## index partner's infection to the hazards in prevalent SDCs. So we need to truncate the
    ## person-months observe in incident couples to 5 months for those observed for longer in
    ## this analysis, and also need to exclude any infections that occurred after.
    ##
    rakll$pm.trunc <- rakll$pm
    rakll$pm.trunc[rakll$phase=='inc' & rakll$pm.trunc>interv/2] <- interv/2
    rakll$inf.trunc <- rakll$inf
    rakll$inf.trunc[inc.wh.hh2] <- 0
    ## For late couples, they excluded any interval right before death in calculations, only
    ## analyzing the 2nd & 3rd intervals before death (also ignoring the 4th) so subtract interv
    ## person-months from these
    rakll$pm.trunc[rakll$phase=='late' & rakll$inf.trunc==0] <- rakll$pm.trunc[rakll$phase=='late' & rakll$inf.trunc==0] - interv
    ## ##################################################################################################
    ## Add variables needed for Hollingsworth et al. style analysis
    ## ##################################################################################################
    rakll$kk <- NA ## intervals followed (or interval of infection for late couples)
    rakll$kkt <- NA ## for late couples total intervals between first observation & index partner death
    ## ########
    ## Incident Couples {ss->hh->hh} k=1; {ss->mm} k = 1; {ss->mm->mm} k =2; {ss->mm->hh->hh} k = 2:
    ## If infected: All non-susceptible visits minus all ++ visits except one.
    rakll$kk[inc.wh.hh] <- apply(ts.vm[,inc.wh.hh,drop=F], 2, function(x) sum(x!='ss',na.rm=T) - (sum(x=='hh',na.rm=T)-1))
    ## If uninfected:All non-susceptible visits.        
    rakll$kk[inc.wh.nhh] <- apply(ts.vm[,inc.wh.nhh,drop=F], 2, function(x) sum(x!='ss',na.rm=T))
    ## Check that calcultions are working
    if(verbose2) {rnd <- sample(inc.wh,10); print(ts.vm[,rnd]); print(rakll[rnd,])}
    ## ########
    ## Prevalent couples {mm->hh} k=1; {mm->hh->hh} k=1; {mm->mm->hh->hh} k = 2; {mm->mm->mm} k=2
    ## If infected: All SDC visits
    rakll$kk[ch.wh.hh] <- apply(ts.vm[,ch.wh.hh,drop=F], 2, function(x) sum(x%in%sdcs, na.rm=T))
    ## If uninfected: All SDC visits - 1
    rakll$kk[ch.wh.nhh] <- apply(ts.vm[,ch.wh.nhh,drop=F], 2, function(x) sum(x%in%sdcs, na.rm=T)) - 1
    ## Check that calcultions are working
    if(verbose2) {rnd <- sample(ch.wh,10); print(ts.vm[,rnd]); print(rakll[rnd,])}
    ## ########
    ## Late couples: more complicated because we have to account for both total follow-up
    ## intervals & interval of infection (if second partner gets infected) which may not be the
    ## same since we are also keeping track of time until index partner's death.
    ## {mm->mm->hh->d.hh} k=2 (interval of infection), kkt=3 (intervals followed before death)
    ## {mm->mm->mm->d.m} k=NA (no infection), kkt = 3
    ## {mm->d.m} k=NA, kkt=1
    ## {mm->mm->d.hh} k = 1 (last inteval before death), kkt=2
    rakll$kkt[lt.wh] <- apply(ts.vm[,lt.wh,drop=F], 2, function(x) sum(x %in%c(sdcs,'hh'), na.rm=T))
    ## If infected: k = last time SDC, kkt = all non-dead observations
    if(length(lt.wh.hh)>0) {
        rakll$kk[lt.wh.hh] <- apply(ts.vm[,lt.wh.hh,drop=F], 2, function(x) min(which(grepl('d\\.',x))) - max(which(x %in% sdcs)))
        ## If uninfected: k = 0, because NA's screw up code later on, but this is meaningless
        rakll$kk[lt.wh.nhh] <- 0}
    ## If uninfected or infected: kkt = all non-dead observations
    ## remove first interval of observation for all individuals watched all 4 intervals
    ## before death (a la Wawer's 6-25 month assumption in Table 2)
    lt.wh.log <- 1:nrow(rakll) %in% lt.wh
    rakll$pm.trunc[lt.wh.log & rakll$kkt==4] <- rakll$pm.trunc[lt.wh.log & rakll$kkt==4] - interv
    ## remove individuals infected in that 4th interval before death
    lt.wh.hh.log <- 1:nrow(rakll) %in% lt.wh.hh
    rakll$inf.trunc[lt.wh.hh.log & rakll$kk==4] <- 0
    return(list(rakll=rakll, inc.wh=inc.wh, inc.wh.hh=inc.wh.hh, inc.wh.nhh=inc.wh.nhh,
                ch.wh=ch.wh, ch.wh.hh=ch.wh.hh, ch.wh.nhh=ch.wh.nhh,
                lt.wh=lt.wh, lt.wh.hh=lt.wh.hh, lt.wh.nhh=lt.wh.nhh))
}

####################################################################################################
## Poisson Model with specified correlation with true heterogeneous variables (for feeding into mclapply)
do.hetmod <- function(het) {
    if(is.na(het)) { ## if not controlling for covariates
        formul <- formula(paste('inf.trunc ~ offset(log(pm.trunc)) + phase', '+'[hps>1], hetproxies[hps]))
    }else{ ## controlling for covariates, create a random covariate with het amount of correlation with true underlying individual risk factors
        temp <- rnorm(nrow(rtrunc), mean = het*rtrunc$secp.lhet, sd = sqrt(het.gen.sd^2 - het^2*het.gen.sd^2))
        formul <- formula(paste('inf.trunc ~ offset(log(pm.trunc)) + phase + temp', '+'[hps>1], hetproxies[hps]))
    }
    temp.arr <- abind(poismod.to.tab(glm(formul, family = "poisson", data = rtrunc)),
                      poismod.to.tab(glm(formul, family = "poisson", data = rtrunc, subset = !excl.by.err)),
                      along = 3)
    dimnames(temp.arr)[[3]] <- c('base', 'XbErr')
    rm(list=setdiff(ls(), "temp.arr")) ## remove everything but output
    gc() ## clean up memory
    return(temp.arr)
}


####################################################################################################
## Poisson Regression (ignoring any source of heterogeneity, i.e. no coital acts, GUD, age, etc)
## assume coital acts are a function of person-months and use that as the offset
## Used in rak.wawer below
poismod.to.tab <- function(mod) {
    poistab <- data.frame(t(cbind(coef(mod), suppressMessages(confint(mod)))))
    if(sum(colnames(poistab) %in% hetproxies)>0) { ## convert to years
        poistab[,colnames(poistab) %in% hetproxies] <- poistab[,colnames(poistab) %in% hetproxies]*12 ## convert het proxy month variables to yars
    }
    poistab <- exp(poistab)
    rownames(poistab) <- c('med','lci','uci')
    poistab <- poistab[c('lci','med','uci'),]
    poistab <- poistab[,!colnames(poistab)=='temp']
    mainpars <- c('bp', 'acute.sc', 'late.sc')
    colnames(poistab)[1:3] <- mainpars
    poistab$dur.ac <- interv/2 ## Wawer et al. assumption that they're seeing 2nd partner at risk for 1/2 of interval
    poistab$dur.lt <- interv 
    poistab$dur.aids <- interv
    if(sum(colnames(poistab) %in% hetproxies)>0) {
        hp <- which(colnames(poistab) %in% hetproxies)
        nord <- c(c(1:ncol(poistab))[-hp], hp)
        poistab <- poistab[, nord]
    }else{
        poistab <- data.frame(poistab, empty = NA)
    }
    poistab <- cbind(poistab, ehm.ac = as.numeric((poistab[,'acute.sc']-1)*poistab[,'dur.ac']),
                     ehm.lt = as.numeric((poistab[,'late.sc']-1)*2*interv), ## poistab[,'dur.lt']), ## assumed to be 2-3rd intervals before death
                     ehm.ltaids = as.numeric((poistab[,'late.sc']-1)*2*interv +(0-1)*0.5*interv)) ## assumed to be last half interval before death
    tdpars <- cbind(t(dpars[parnames]), empty = 1, ehm.ac = as.numeric((dpars['acute.sc']-1)*dpars['dur.ac']),
                    ehm.lt = as.numeric((dpars['late.sc']-1)*dpars['dur.lt']),
                    ehm.ltaids = as.numeric((dpars['late.sc']-1)*dpars['dur.lt'] +(0-1)*dpars['dur.aids']))
    hetp.nm <- colnames(poistab)[colnames(tdpars)=='empty']
    colnames(tdpars)[colnames(tdpars)=='empty'] <- hetp.nm
    poistab <- rbind(poistab, tdpars)
    rownames(poistab)[4] <- 'true'
    tracenames <- c(tracenames, hetp.nm)
    return(poistab[,tracenames])
}


####################################################################################################
## Wawer et al. style analysis of Rakai retrospective cohort
rak.wawer <- function(rak.coh, verbose = F, verbose2=F, browse = F, excl.extram = T, decont=F, start.rak=1994, het.gen.sd, late.ph,
                      resamp=F, cov.mods=T, fit.Pois=T, 
                      prop.controlled = c(NA,seq(0, 1, by = .1)), hetproxies = '') { ## amount of heteroeneity controlled for, other covariates to add
    if(browse) browser()
    ts.vm <- rak.coh$ts.rak
    ts.vm.all <- rak.coh$ts.rak.all
    dat.vm <- rak.coh$dat.rak
    dpars <- rak.coh$dpars
    if(verbose) {
        print('EHMs as inputted:')
        print(truehrs(ts.vm.all, dat.vm, pars  = dpars))
    }
    interv <- rak.coh$interv
    rm(rak.coh) ## to release memory
    ## Deal with extra-couply infected 2n partners
    sel <- which(apply(ts.vm, 2, function(x) sum(grepl('hh',x))>0))
    sel.m2e <- sel[dat.vm$mcoi[sel]=='e' & dat.vm$mdoi[sel] > dat.vm$fdoi[sel]] # male 2nd
    sel.f2e <- sel[dat.vm$fcoi[sel]=='e' & dat.vm$fdoi[sel] > dat.vm$mdoi[sel]] # female 2nd
    sel.b2e <- sel[dat.vm$mcoi[sel]=='e' & dat.vm$fcoi[sel]=='e' & dat.vm$mdoi[sel] == dat.vm$fdoi[sel]] # both same time
    extram <- c(sel.m2e, sel.f2e, sel.b2e)
    if(excl.extram) { ## Remove couples where 2nd partner was infected extra-couly
        rem <- 1:nrow(dat.vm) %in% extram
        dat.vm <- dat.vm[!rem,]
        ts.vm <- ts.vm[,!rem]
        ts.vm.all <- ts.vm.all[,!rem]
        if(verbose) {
            print('EHMs after excluding couples with 2nd partner infected extra-couply:')
            print(truehrs(ts.vm.all, dat.vm, pars  = dpars))
        }
    }else{ ## Censor couples where 2nd partner was infected extra-couply from infection forward
        which(ts.vm[,extram]=='hh',arr.in=T)
        ## change all ++ to NA's (this doesn't censor perfectly since we should assume infection midpoint in the interval.
        ts.vm[,extram][which(grepl('hh', ts.vm[,extram]),arr.in=T)] <- NA
    }
    ncpl <- ncol(ts.vm)
    ## Create line list
    rakllout <- make.rakll(dat.vm = dat.vm, ts.vm = ts.vm, cov.mods=cov.mods, interv=interv, verbose2=verbose2)
    ## ##################################################################################################
    ## Look at contamination between phases
    if(decont) { ## remove couples that are in wrong grouping (e.g. prevalent  couples with person-time spent in acute/late/aids phase)
        contam <- with(rakllout, {
            inc.contam <- inc.wh[which(apply(ts.vm.all[,inc.wh], 2, function(x) sum(grepl('lt',x) | grepl('aids',x))>0))]
            ch.contam.inc <- ch.wh[which(apply(ts.vm.all[,ch.wh], 2, function(x) sum(grepl('ac',x))>0))]
            ch.contam.lt <- ch.wh[which(apply(ts.vm.all[,ch.wh], 2, function(x) sum(grepl('lt',x) | grepl('aids',x))>0))]
            ch.contam <- unique(c(ch.contam.inc, ch.contam.lt))
            lt.contam <- lt.wh[which(apply(ts.vm.all[,lt.wh], 2, function(x) sum(grepl('ac',x))>0))]
            contam <- 1:ncol(ts.vm) %in% c(inc.contam, ch.contam, lt.contam)
            return(contam)
        })
        ts.vm <- ts.vm[,!contam]
        ts.vm.all <- ts.vm.all[,!contam]
        dat.vm <- dat.vm[!contam,]
        if(verbose) {
            print('EHMs after decontaminating person-time between phases (e.g. late phase person-time in prevalent couples):')
            print(truehrs(ts.vm.all, dat.vm, pars  = dpars))
        }
        rakllout <- make.rakll(dat.vm = dat.vm, ts.vm = ts.vm, cov.mods=cov.mods, interv=interv, verose2=verbose2) 
    }
    ## ##################################################################################################
    ## empirical hazards
    erhs <- with(rakllout, {   
        ehz.inc <- sum(rakll$inf.trunc[inc.wh]) / sum(rakll$pm.trunc[inc.wh])
        ehz.ch <- sum(rakll$inf.trunc[ch.wh]) / sum(rakll$pm.trunc[ch.wh])
        ehz.lt <- sum(rakll$inf.trunc[lt.wh]) / sum(rakll$pm.trunc[lt.wh])
        ## empirical hazard ratio
        e.arh <- ehz.inc/ehz.ch
        e.lrh <- ehz.lt/ehz.ch
        ## empirical hazards excluding couples only observed +- at one point (& excluded by error)
        sel <- inc.wh[!inc.wh %in% which(rakll$excl.by.err)]
        ehz.inc.err <- sum(rakll$inf.trunc[sel]) / sum(rakll$pm.trunc[sel])
        sel <- ch.wh[!ch.wh %in% which(rakll$excl.by.err)]        
        ehz.ch.err <- sum(rakll$inf.trunc[sel]) / sum(rakll$pm.trunc[sel])
        sel <- lt.wh[!lt.wh %in% which(rakll$excl.by.err)]        
        ehz.lt.err <- sum(rakll$inf.trunc[sel]) / sum(rakll$pm.trunc[sel])
        ## empirical hazard ratio
        e.arh.err <- ehz.inc.err/ehz.ch.err
        e.lrh.err <- ehz.lt.err/ehz.ch.err
        ## store empirical hazards & their ratios as a vector for output
        erhs <- c(e.arh = e.arh, e.lrh = e.lrh, e.arh.err = e.arh.err, e.lrh.err = e.lrh.err,
                  ehz.inc = ehz.inc, ehz.ch = ehz.ch, ehz.lt = ehz.lt,
                  ehz.inc.err = ehz.inc.err, ehz.ch.err = ehz.ch.err, ehz.lt.err = ehz.lt.err)
        return(erhs)
    })
    ## ##################################################################################################
    parnames <- c('bp','acute.sc','dur.ac','late.sc','dur.lt','dur.aids')
    tracenames <- c('bp','ehm.ac','acute.sc','dur.ac','ehm.ltaids', 'ehm.lt','late.sc','dur.lt','dur.aids')
    ## ##################################################################################################
    ## remove late couples for which we zero-ed out their person months (i.e. only observed in 4th or
    ## 1st interval before death)
    rtrunc <- rakllout$rakll 
    rtrunc <- rtrunc[rtrunc$pm.trunc>0,]
    dat.vm <- dat.vm[rtrunc$pm.trunc>0,]
    ts.vm.all <- ts.vm.all[,rtrunc$pm.trunc>0]
    if(verbose) {
        print('EHMs after removing late couples not seen in the 2nd-3rd intervals prior to death')
        print(truehrs(ts.vm.all, dat.vm, pars  = dpars))
        print('EHMs when tabulating within *assigned* phase (i.e. omnitient tallying of person-months exposed to a given phase, and infections really due to that phase but *within* assigned phases')
        inc.ac.hz <- truehrs(ts.vm.all, dat.vm, pars  = dpars, rakllout$inc.wh)$hzs['ac']
        prev.ch.hz <- truehrs(ts.vm.all, dat.vm, pars  = dpars, rakllout$ch.wh)$hzs['ch']
        late.lt.hz <- truehrs(ts.vm.all, dat.vm, pars  = dpars, rakllout$lt.wh)$hzs['lt']
        ehm.ac.phase <- (inc.ac.hz/prev.ch.hz - 1)* dpars['dur.ac']
        ehm.lt.phase <- (late.lt.hz/prev.ch.hz - 1)* dpars['dur.lt']
        ehm.ltaids.phase <- (late.lt.hz/prev.ch.hz - 1)* dpars['dur.lt'] + (0-1)*dpars['dur.aids']    
        print(signif(c(ehm.ac.phase, ehm.lt.phase, ehm.ltaids.phase),3))
        print('infections by phase')
        print(xtabs(inf.trunc ~ phase, rtrunc))
        print('person-months by phase')    
        print(xtabs(pm.trunc ~ phase, rtrunc))
        print('empirical hazards by phase (not omnitient)')
        print(xtabs(inf.trunc ~ phase, rtrunc) / xtabs(pm.trunc ~ phase, rtrunc))
    }
    ## Do several models
    obs.hets <- paste0('obs',prop.controlled) ## name variables
    ## should we include any other proxy of heterogeneity in the model?
    gc()
####################################################################################################
    if(fit.Pois) {
        for(hps in 1:length(hetproxies)) { ## for each covariate that could be included which might be indicative of heterogeneity
            print(paste('fitting Poisson models with heterogeneity &', hetproxies[hps]))
            tout <- mclapply(prop.controlled, do.hetmod)
            if(hps==1) {
                armod <- abind(tout,along = 4)
            }else{
                armod <- abind(armod, abind(tout,along = 4), along = 5)
            }
        }
        dimnames(armod)[[4]] <- obs.hets
        dimnames(armod)[[5]] <- hetproxies    
        armod['med','ehm.ac','base',,]
        armod['true','ehm.ac',1,1,1]
        armod['med','ehm.ltaids','base',,]
        armod['true','ehm.ltaids',1,1,1]
        if(verbose) { ## examine het proxies a bit
            print(ddply(rtrunc, .(phase), summarise, mean.indiv.RH = exp(mean(secp.lhet))))
            rtrunc <- ddply(rtrunc, .(), transform, secp.age.cat = cut(secp.age/12, seq(0,65, by = 5)),
                            indp.age.cat = cut(indp.age/12, seq(0,65, by = 5)),
                            mardur.cat = cut(mardur/12+1, seq(0,65, by = 5)),
                            secp.tdsa.cat = cut(secp.tdsa/12+1, seq(0,65, by = 5)),
                            secp.pdsa.cat = cut(secp.pdsa/12+1, seq(0,65, by = 5)))
            print(ddply(rtrunc, .(secp.age.cat), summarise, mean.indiv.RH = exp(mean(secp.lhet))))
            print(ddply(rtrunc, .(indp.age.cat), summarise, mean.indiv.RH = exp(mean(secp.lhet))))
            print(ddply(rtrunc, .(mardur.cat), summarise, mean.indiv.RH = exp(mean(secp.lhet))))
            print(ddply(rtrunc, .(secp.tdsa.cat), summarise, mean.indiv.RH = exp(mean(secp.lhet))))
            print(ddply(rtrunc, .(secp.pdsa.cat), summarise, mean.indiv.RH = exp(mean(secp.lhet))))        
        }
    }else{armod <- NA}
    ## Check that calcultions are working
    if(verbose2) {with(rakllout, {rnd <- sample(lt.wh,10); print(ts.vm[,rnd]); print(rakll[rnd,])})}
    rakll <- rakllout$rakll
    rm(list=setdiff(ls(), c("erhs","rakll","armod","PoisRHs"))) ## remove everything but output
  return(list(erhs = erhs, rakll = rakll, armod = armod))
  gc()
}

####################################################################################################
####################################################################################################
## Functions for Hollingsworth style Likelihood
####################################################################################################
## beta(t) for incident couples (function of time since index partner infection tt)
b.inc <- function(tt, dpars)  {
   for(nm in names(dpars))      assign(nm, dpars[nm])
    return(ifelse(tt <= dur.ac, bp* acute.sc, bp))
}
## beta(t) for late couples (function of time until index partner death td)
b.lt.uv <- function(td,dpars) {
  for(nm in names(dpars))      assign(nm, dpars[nm])
  if(td < dur.aids) { ## if within dur.aids of death, no transmission
    bet <- 0
  }else{ ## if td is in (dur.aids+dur.lt) to dur.aids, it's late phase
    if(td <= (dur.aids + dur.lt)) {
      bet <- bp * late.sc
    }else{ ## anything earlier than (dur.aids+dur.lt) is still chronic
      bet <- bp
    }
  }
  return(bet)
}
b.lt <- Vectorize(b.lt.uv,'td')
####################################################################################################
## 1) Conditional probability of transmitting in each inteval given it didn't happen previously
##############################
## Incident Couples
##############################
## probability of transmitting to second partner before end of 1st interval (function of
## when in (-- to ++) interval index partner is infected = tt)
cp.inc.1 <- function(tt, dpars) {
    for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
    ## print(1 - exp(-integrate(b.inc, lower=0, upper=interv - tt, dpars = dpars)$val)) ## function should return something equivalent to this line    
    if(dur.ac <= interv - tt) return(1 - exp( - bp*acute.sc*dur.ac - bp*(interv-tt-dur.ac) ))
    if(dur.ac >  interv - tt) return(1 - exp( - bp*acute.sc*(interv-tt))) ## acute phase is full first interval
}
 
## probability of transmitting to second partner before end of kth interval (function of
## when in (-- to +-) interval index partner is infected = tt) STARTS at 2nd interval
cp.inc.k <- function(tt,kk,dpars) {
  for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
  ## print(1 - exp(-integrate(b.inc, lower=interv*(kk-1) - tt, upper=interv*kk - tt, dpars = dpars)$val)) ## function should return something equivalent to this line
  if(interv*(kk-1) >= tt + dur.ac) return(1-exp(-bp*interv)) ## chronic phase only
  if(interv*(kk-1) < tt + dur.ac)  { ## some acute phase still in interal
    dur.ac.temp <- min(dur.ac - (interv*(kk-1) - tt), interv)
    dur.ch.temp <- interv-dur.ac.temp
    return(1-exp(-bp*acute.sc*dur.ac.temp - bp*dur.ch.temp))
  }
}

##############################
## Prevalent Couples (assumed to be the same in all intervals)
##############################
cp.prev.k <- function(kk=NA,dpars)  {
       for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
        1 - exp(-bp*interv)
    }
##############################
## Late Couples
##############################
cp.lt.1 <- function(td, dpars) {
  ## print(1 - exp(-integrate(b.lt, lower=0, upper=interv - td, dpars = dpars)$val)) ##  function should return something equivalent to this line
  for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
  if(interv - td > dur.lt + dur.aids) { ## interval before death also includes some chronic infectivity
    dur.ch <- interv - td - dur.lt - dur.aids
    cp <- 1 - exp(-bp*dur.ch - bp*late.sc*dur.lt - bp*0*dur.aids)
  }else{
    if(interv - td > dur.aids) { ## interval before death includes some but not all late phase infectivity
      dur.lt.temp <- interv - td - dur.aids
      cp <- 1 - exp(-bp*late.sc*dur.lt.temp - bp*0*dur.aids)
    }else{ ## interval before death only includes AIDS phase infectivity
      dur.aids.temp <- interv - td
      cp <- 1 - exp(-bp*0*dur.aids.temp)
    } }
  return(cp)
}
  
## probability of transmitting to second partner in k-th interval BEFORE death (starts at
## interval 2, since previous line is for 1st interval before index partner death)
cp.lt.k <- function(td,kk,dpars) {
#  cpi <- 1 - exp(-integrate(b.lt, lower=interv*(kk-1) - td, upper=interv*kk - td, dpars = dpars)$val)
  for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
  if(interv*(kk-1) - td >= dur.aids + dur.lt) { ## interval only includes chronic phase
    cp <-  1 - exp(-bp*interv)
  }else{ ##  interval includes some but not all late phase
    if(interv*(kk-1) - td >= dur.aids) {
      dur.lt.temp <- min(interv, dur.lt - (interv*(kk-1) - td - dur.aids))
      dur.ch <- max(0, interv - dur.lt.temp) ## the rest is chronic
      cp <- 1 - exp(-bp*late.sc*dur.lt.temp - bp*dur.ch)
    }else{ ## interval includes all aids phase & some but not all late phase
      dur.aids.temp <- min(dur.aids - (interv*(kk-1) - td),interv)
      dur.lt.temp <- min(interv - dur.aids.temp, dur.lt,interv)
      dur.ch <- max(0, interv - dur.aids.temp - dur.lt.temp)
      cp <- 1 - exp(-bp*dur.ch - bp*late.sc*dur.lt.temp - bp*0*dur.aids.temp)
    } }
#  if(abs(cp-cpi)>10^-10) browser() #stop('analytic integral error')
return(cp)
} 
## debug(ucp.lt)
## ucp.lt(kk = kk.temp, kkt=kkt.temp, inf = 1, dpars=dpars)
## debug(cp.lt.k)
## cp.lt.k(td=4, kk=jj, dpars)

####################################################################################################
## 2) Unconditional probability of transmitting in a particular inteval
####################################################################################################
## Incident Couples: k = interval since -- (1st, 2nd, 3rd, 4th, etc..)
ucp.inc.2 <- function(inct,dpars) {
  for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
  ## numerical integrals (to check calculus)
  pps <- integrate(Vectorize(cp.inc.1,'tt'), lower=0, upper=interv, dpars = dpars)$val / interv ## 
    for(jj in 2:nrow(inct)) {
      ## numerical integrals
      pps <- c(pps, integrate(Vectorize(cp.inc.k,'tt'), lower=0, upper=interv, kk=jj, dpars = dpars)$val / interv)
    }
  lls <- inct$i * log(pps) + (inct$n-inct$i) * log(1-pps)
  nll.inc <- -sum(lls)
  return(nll.inc)
}

####################################################################################################
##############################
## Prevalent Couples (assumed to be the same in all intervals)
##############################
ucp.prev.2 <- function(prevt,dpars) { 
    for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
    pp <- 1-exp(-bp*interv)
    lls <- prevt$i * log(pp) + (prevt$n-prevt$i) * log(1-pp)
    nll.prev <- -sum(lls)
    return(nll.prev)
}

####################################################################################################
##############################
## Late Couples (assumed to be the same in all intervals)
##############################
## probability of seroconversion in interval kk BEFORE death, given not in previous
## intervals, **when observed for a total of kkt intervals before death**
ucp.lt.2 <- function(latet, dpars, browse=F) {
  for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
  ## i1p
#  print(dpars)
#  browser()
  if(browse) browser()
  max.int <- max(0,(interv-dpars['dur.aids']))
  i1p <-   try(integrate(Vectorize(cp.lt.1,'td'), lower=0, upper=max.int, dpars=dpars)$val / interv, silent = T)
  w.step <- 0
  while(inherits(i1p, 'try-error')){
    if(w.step>30) ucpN <- 0 # give up on this integral eventually & just reject this proposal
    max.int <- max.int*.98
    i1p <-   try(integrate(Vectorize(cp.lt.1,'td'), lower=0, upper=max.int, dpars=dpars)$val / interv, silent = T)
    w.step <- w.step+1
  }
  pps <- i1p
  #browser()
  #ll.lt <- latet$i[latet$int==1]*log(i1p) + (latet$n[latet$int==1]-latet$i[latet$int==1])*log(1-i1p)
  #print(i1p)
  #ucN <- 1-i1p
  #temp <- (latet$n[latet$int==1]-latet$i[latet$int==1])*log(1-i1p)
  for(vv in 2:(max.vis-1)) { # remaining intervals observed (2nd before death, 3rd, etc..)
    ## ivvp
    max.int <- min(10,max(0,(interv*vv-dpars['dur.aids'])))
    ivvp <- try(integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=vv, dpars = dpars)$val / interv, silent=T)
    w.step <- 0
    while(inherits(ivvp, 'try-error')){
      if(w.step>30) ucpN <- 0 # give up on this integral eventually & just reject this proposal
      max.int <- max.int*.98
      ivvp <- try(integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=vv, dpars = dpars)$val / interv, silent=T)
      w.step <- w.step+1
    }
    #ucN <- ucN*(1-ivvp)
    pps <- c(pps, ivvp)
    ## i1p^i1 * (1-i1p)^(n1-i1) etc...
    #ll.lt <- ll.lt + latet$i[latet$int==vv]*log(ivvp) + (latet$n[latet$int==vv]-latet$i[latet$int==vv])*log(1-ivvp)
    #temp <- temp + (latet$n[latet$int==1]-latet$i[latet$int==1])*log(1-ivvp)
  }
  pps <- rev(pps) ## since 1st is last row in latet
  pps <- pps[(max.vis-1):1 %in% latet$int] ## get rid of any where there were no individuals even observed
  ## !=0 indicing below is to avoid 0*log(0) returning NaN. Instead it
  ## !should not be added in (0 events obesrved)
  lls.inf <- latet$i[latet$i!=0] * log(pps[latet$i!=0]) # those infected
  lls.ninf <- (latet$n-latet$i)[(latet$n-latet$i)!=0] * log(1-pps[(latet$n-latet$i)!=0]) # those uninfected
  nll.lt <- -sum(lls.inf)-sum(lls.ninf)
  return(nll.lt)
}
 
 
## turn SB simulation into Wawer style table
sbmod.to.wdat <- function(sim, browse=F, excl.by.err = F, giveLate = T, giveProp = F, condRakai=F, RakSamp = c(inc = 23, prev = 161, late=51), simpPois=F) {    
    ## Excluding incident couples seen serodiscordant once & then never again as in Wawer 2005?
    if(browse) browser()
    if(excl.by.err) sim <- sim[!sim$excl.by.err,]
    if(!sum(sim$inf[sim$phase=='inc'])==0 & !sum(sim$inf[sim$phase=='prev'])==0 & !sum(sim$phase=='inc')==0) { ##otherwise just return NA
    ## Condition on Rakai sample sizes by sampling with replacement
    if(condRakai) {
        inc.wh <- which(sim$phase=='inc')
        prev.wh <- which(sim$phase=='prev')
        if(giveLate) {prev.wh <- which(sim$phase=='prev')
                      resamp <- c(sample(inc.wh, RakSamp['inc'], replace=T), sample(prev.wh, RakSamp['prev'], replace=T),sample(late.wh, RakSamp['late'], replace=T))
                  }else{ ## only inc & prev
                      resamp <- c(sample(inc.wh, RakSamp['inc'], replace=T), sample(prev.wh, RakSamp['prev'], replace=T))
                  }
        sim <- sim[resamp,]
    }
        ## early
        sim$inf <- factor(sim$inf, levels=c(0,1))
        inctab <- xtabs(~inf+kk, sim, subset=phase=='inc')
        inct <- try(data.frame(int = 1, n = sum(sim$phase=='inc'), i = inctab[2,1]))
        if(inherits(inct, 'try-error')) {print(inctab); print(sum(sim$phase=='inc'))}
        if(ncol(inctab)>1) {
            for(ii in 2:min(4,ncol(inctab))) {
                temp <- data.frame(int = ii,
                                   n =  inct$n[ii-1] - inct$i[ii-1] - inctab[1,ii-1],
                                   i = inctab[2,ii])
                inct <- rbind(inct, temp)
            }}
        ## prev
        prevtab <- xtabs(~inf+kk, sim, subset=phase=='prev')
        prevt <- data.frame(int = 1, n = sum(sim$phase=='prev'), i = prevtab[2,1])
        for(ii in 2:ncol(prevtab)) {
            temp <- data.frame(int = ii,
                               n =  prevt$n[ii-1] - prevt$i[ii-1] - prevtab[1,ii-1],
                               i = prevtab[2,ii])
            prevt <- rbind(prevt, temp)
        }
        ## just get an unadj Pois & an omnitient Pois reg (controlling for all heterogeneity), ignoring late phase
        ## only do if there aren't 0's in any of the inf/phase categories (neither all infected, or none infected in a phase), otherwise return NA (commented for now
        if(simpPois) { # & sum(xtabs(~inf + phase, temprcoh$rakll)[,c('prev','inc')]==0)==0) {
            formul.uni <- formula('inf.trunc ~ offset(log(pm.trunc)) + phase')
            acuteRH.uni <- exp(coef(glm(formul.uni, family = "poisson", data = sim[sim$phase!='late',])))['phaseinc']
            formul.mult <- formula('inf.trunc ~ offset(log(pm.trunc)) + phase + secp.lhet')
            acuteRH.mult <- exp(coef(glm(formul.mult, family = "poisson", data = sim[sim$phase!='late',])))['phaseinc']
            PoisRHs <- c(univ = acuteRH.uni, omn = acuteRH.mult)
        }else{ PoisRHs <- c(NA,NA)}
        if(giveLate) {
            ## late
            latetab <- xtabs(~inf+kk, sim, subset=phase=='late')
            maxints <- max(sim$kkt, na.rm=T)
            if(maxints==-Inf) maxints <- 1
            latet <- data.frame(int = 1:1, n = 0, i = 0)
            for(ii in 1:nrow(latet)) {
                int <- latet$int[ii]
                ##wh.tm <- sim$phase=='late' & (sim$kkt == int | (sim$kkt>int & sim$kk<=int))
                ##print(head(sim[wh.tm,],20))
                latet$n[latet$int==int] <- sum(sim$phase=='late' & (sim$kkt == int | (sim$kkt>int & sim$kk<=int)) )
                latet$i[latet$int==int] <- sum(sim$phase=='late' & sim$kk==int)
            }              
            head(sim[sim$phase=='late',],10)
            if(giveProp) {
                inct$p <- with(inct, i/n)
                prevt$p <- with(prevt, i/n)
                if(giveLate) latet$p <- with(latet, i/n)
            }
            return(list(inct=inct, prevt=prevt, latet= latet))
        }else{
            if(giveProp) {
                inct$p <- with(inct, i/n)
                prevt$p <- with(prevt, i/n)
            }
            return(list(inct=inct, prevt=prevt, PoisRHs=PoisRHs))}
    }else{
        return(list(inct=NA, prevt=NA, PoisRHs=c(NA,NA), error='No infections in incident couples'))
    }
}

####################################################################################################
## Calculate Hollingsworth et al. style likelihood
holl.lik <- function(ldpars, wtab, ## log-parms, Wawer Table 1 type tables
                     range.dur.ac = c(.25,10), ##  dur.ac make flat on [.25,10] months to keep things bounded
                     range.dur.lt = c(.5,36), ## flat bounded prior on late phase
                     range.dur.aids = c(.5, 36),  ## flat bounded prior on aids phase
                     late.ph = T, ## include late phase couples in likelihood or just compare incident/prevalent couples
                     verbose = T, browse = F)
  {
    if(browse) browser()
    ## exponentiate ldpars
    dpars <- exp(ldpars)
    ## Incident Couples
    inc.nll <- ucp.inc.2(wtab$inct, dpars)
    if(verbose) print(inc.nll)
    ## Prevalent Couples
    prev.nll <- ucp.prev.2(wtab$prevt, dpars)
    if(verbose) print(prev.nll)    
    if(late.ph) { ## if including late phase
      ## Late Couples
      lt.nll <- ucp.lt.2(wtab$latet, dpars = dpars, browse=F)
      if(verbose) print(lt.nll)
      nll <- inc.nll + prev.nll + lt.nll
    }else{
      nll <- inc.nll + prev.nll
    }
    nll <- as.numeric(nll)
    if(length(range.dur.ac)>0) { ## give 0 posterior probability (nll - log(prior)=Inf) if outside range of dur.ac
      if(dpars['dur.ac'] < range.dur.ac[1] | dpars['dur.ac'] > range.dur.ac[2]) nll <- Inf
    }
    if(length(range.dur.lt)>0) { ## give 0 posterior probability (nll - log(prior)=Inf) if outside range of dur.lt
      if(dpars['dur.lt'] < range.dur.lt[1] | dpars['dur.lt'] > range.dur.lt[2]) nll <- Inf
    }
    if(length(range.dur.aids)>0) { ## give 0 posterior probability (nll - log(prior)=Inf) if outside range of dur.aids
      if(dpars['dur.aids'] < range.dur.aids[1] | dpars['dur.aids'] > range.dur.aids[2]) nll <- Inf
    }
    if(verbose) print('total nll:')
    return(nll)
  }

## function to calculate attributable chronic-infection-months due to elevated acute phase infectivity
arr <- function(ldp) as.numeric((exp(ldp['acute.sc'])-1)*exp(ldp['dur.ac']))

## Simulate dat using exactly Hollingsworth's model to make sure we're fitting it correctly
holl.mod <- function(i.n = 50, p.n = 50, l.n = 50, interv = 10,  max.vis = 4, dpars, verbose=F, browse = F) {
  for(nm in names(dpars))      assign(nm, as.numeric(dpars[nm])) ## loading 'bp'
  if(browse) browser()
  dat <- data.frame(uid=1:sum(i.n,p.n,l.n), phase = c(rep('inc',i.n), rep('prev',p.n), rep('late',l.n)), inf = NA, kk = NA, kkt = NA)
  ## Inc
  for(ii in 1:i.n) {
    if(verbose) print(paste('inc.ii',ii))
    t.inf1 <- runif(1,0,interv) ## time of index partner infection
    inf <- 0 ## Set secondary partner infection status 0
    ac.exp <- min(interv-t.inf1, dur.ac) ## duration of acute phase experienced in first interval
    cp1 <- 1-exp(-(bp*acute.sc*ac.exp + bp*(interv-t.inf1-ac.exp))) ##  cumulative probability of infection in first interval
    inf <- rbinom(1,1,cp1)                                         ##  Bernoulli infection variable
    vis <- 1 ## current visit number for while loop
    if(verbose)       print(paste(vis,signif(cp1,3)))
    while(inf==0 & vis<max.vis) {       # while and infected and before last visit
      ac.exp <- max(t.inf1+dur.ac - vis*interv, 0) ## acute phase expeienced in vis-th intrval
      cp <- 1-exp(-bp*acute.sc*ac.exp - bp*(interv-ac.exp)) ## cumulative probability of infection in this interval
      inf <- rbinom(1,1,cp)## Bernoulli infection variable
      vis <- vis+1 ## increment visit
      if(verbose)       print(paste(vis,signif(cp,3)))
    }
    dat$inf[ii] <- inf
    dat$kk[ii] <- vis
  }
  ## Prev
  for(ii in 1:p.n) {
    if(verbose) print(paste('prev.ii',ii))
    inf <- 0 ## Bernoulli infection variable
    vis <- 0 ##  visit number
    while(inf==0 & vis<max.vis) {
      cp <- 1-exp(-bp*interv) ## cumulative probability of infection within each interval
      if(verbose)       print(paste(vis,signif(cp,3)))
      inf <- rbinom(1,1,cp)
      vis <- vis+1 ## increment visit number
    }
    dat$inf[i.n + ii] <- inf
    dat$kk[i.n + ii] <- vis
  }
  ## Late
  for(ii in 1:l.n) {
    if(verbose) print(paste('late.ii',ii))
    inf <- 0 ##  Bernoulli infection variable
    vis <- 0 ## visit number
    t.dth <- runif(1,0,interv) ## time between death & visit after death (visit 5)
    while(inf==0 & vis<max.vis) {
      cp <- ifelse((max.vis-vis)>1, cp.lt.k(t.dth,(max.vis-vis),dpars), cp.lt.1(t.dth,dpars)) # use cumulative probability functions written above
      if(verbose)       print(paste(vis,signif(cp,3)))
      inf <- rbinom(1,1,cp)
      vis <- vis+1
    }
    dat$inf[i.n+p.n+ii] <- inf
    dat$kkt[i.n+p.n+ii] <- max.vis
    if(inf==1)  {
      dat$kk[i.n+p.n+ii] <- max.vis - vis + 1
    }else{
      dat$kk[i.n+p.n+ii] <- 0
    }
  }
  dat$phase <- factor(dat$phase, levels = c('inc','prev','late'))
  return(dat)
}

## Fit Hollingsworth model with MCMC
## MCMC SAMPLER
hollsampler <- function(sd.props, inits,
                        wtab, fixed = NULL,
                        multiv = F, covar = NULL, # if multiv, sample from multivariate distribution (calculated during adaptive phase)
                        late.ph = F, # include late phase?
                        verbose = T, verbose2 = F, tell = 100, seed = 1, 
                        niter = 6*1000, nthin = 5, nburn = 1000, browse=F)
  {
    if(browse)  browser()
    set.seed(seed)
    pars <- inits
    parnames <- names(pars)
    if(length(fixed)>0) {
      for(iii in 1:length(fixed))       pars[names(fixed)[iii]] <- fixed[iii]
      to.fit <- !names(pars) %in% names(fixed)
    }else{
      to.fit <- names(pars)
    }
    vv <- 2
    accept <- 0                  #track each parameters acceptance individually
    lprob.cur <- holl.lik(pars, wtab, verbose = F, browse = F)
    out <- as.data.frame(c(pars, nll =lprob.cur))
    last.it <- 0
    start <- Sys.time()
    while(vv < niter + 1) {
        if(verbose & vv%%tell+1==1) print(paste("on iteration",vv,"of",last.it + niter + 1))
        pars.prop <- pars              #initialize proposal parameterr vector
        ## propose new parameter vector
        if(multiv) {
            if(!vv%%5==0) { ## block sample
                pars.prop[to.fit] <- pars.prop[to.fit] + rmnorm(1, mean = 0, varcov = covar)
                pars.prop[to.fit] <- as.vector(pars.prop[to.fit]) #otherwise is a matrix
               names(pars.prop) <- parnames
            }else{ ## every 5 iterations just sample dur.aids since it appears to get stuck in weird
                   ## places during block sampling
                pars.prop['dur.aids'] <- pars.prop['dur.aids'] + rnorm(1, mean = 0, sd = covar[6,6])
                names(pars.prop) <- parnames
            }
        }else{
          pars.prop[to.fit] <- pars.prop[to.fit] + rnorm(length(pars[to.fit]), mean = 0, sd = sd.props[to.fit])
        }
        ## trace = T if in non-thinned iteration, or the previous one (in case of rejection)
        ## calculate proposal par log probability
        if(verbose2)    print(vv)
        lprob.prop <- holl.lik(pars.prop, wtab, verbose = F, browse = F)
        ##if(is.na(lprob.prop)) browser()
        if(verbose2) print(lprob.prop)
        ## holl.lik gives -log(L), so negative of that gives the log probability (assuming flat improper priors)
        lmh <- -lprob.prop + lprob.cur       # log Metropolis-Hastings ratio
        if(is.na(lmh)) print(paste0('lmh=',lmh,', pars.prop=',paste(exp(pars.prop), collapse=', '), ', vv=',vv,', seed=',seed))
        if(verbose2) print(lmh)
        ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
        if(lmh >= 0 | runif(1,0,1) <= exp(lmh)) {
          pars <- pars.prop
          if(vv>nburn) accept <- accept + 1 #only track acceptance after burn-in
          lprob.cur <- lprob.prop
        }
        if(vv%%nthin + 1 ==1)   {
          out <- cbind(out, as.data.frame(c(pars, nll =lprob.cur))) ## if not thinning, then record it
          vv <- vv+1
        }
      }
    colnames(out) <- 1:ncol(out)
    rownames(out) <- c(names(pars),'nll')
    if(verbose) print(paste("took", difftime(Sys.time(),start, units = "mins"),"mins"))
    aratio <- accept/((vv-nburn))    
    return(list(out = out[,1:ncol(out)>(nburn+1)/nthin], aratio = aratio, inits = inits))
  }

holl.init.fxn <- function(seed = 1) {
  set.seed(seed)
    ldpars <- c(acute.sc = runif(1, log(1), log(50)),
                late.sc = runif(1, log(1), log(50)),
                bp = runif(1, -6, -3),
                dur.ac = runif(1, log(1), log(6)),
                dur.lt = runif(1, log(3), log(12)),
                dur.aids = runif(1, log(3), log(12)))
    ldpars
  }

holl.wrp <- function(seed=1, sd.props, wtab, force.inits=NULL,
                multiv=F, covar=NULL, jit = .8,
                verbose = T, verbose2 = F, tell = 50, browse = F,
                niter, nthin, nburn)
  { ## new version of wrp needs to save progress for long chains (WA)
    if(length(seed)>1) browse <- F ## can't browse within mclapply
    if(length(force.inits)==0) {
        inits.temp <- holl.init.fxn(seed = seed) # initial conditions different for each seed
    }else{
        inits.temp <- jitter(force.inits, a = jit)
    }
    hollsampler(sd.props = sd.props, inits = inits.temp, wtab = wtab,
                multiv = multiv, covar = covar, 
                verbose = verbose, verbose2 = verbose2, tell = tell, seed = seed, 
                niter = niter, nthin = nthin, nburn = nburn, browse=browse)
  }


## Plot posterior pairwise-correlations & histograms
hollsbpairs <- function(posts, truepars = NULL, file.nm, width = 10, height = 10, show.lines = T, range0 = NULL,
                        inits.tp = NULL, # data frame of inits, columns matching posts
                        cex = 1, col = "black", nrpoints = 200, do.pdf = F, do.jpeg = T, ranges = NULL,
                        cex.axis = 1.5, cex.nm = 2.5, greek = T, show.cis = T, browse = F)
  {
    if(browse) browser()
    if(do.pdf) pdf(paste(file.nm, ".pdf", sep=""), width = width, height = height)
    if(do.jpeg) jpeg(paste(file.nm, ".jpeg", sep=""), width = width*100, height = height*100)
    if(length(ranges)==0)   ranges <- apply(rbind(posts,truepars), 2, function(x) range(x,na.rm=T)) ## not putting inits in here, they'll only show up if in the region.
    if(length(range0)>1) ranges[1,range0] <- 0
    cis <- apply(rbind(posts,truepars), 2, function(x) quantile(x, c(.025, .975), na.rm=T))
    par(mar=rep(3,4),oma=rep(2,4),mfrow=rep(ncol(posts),2))
    parnames <- colnames(posts)
    for(ii in 1:ncol(posts)) {
        for(jj in 1:ncol(posts)) {
            if(ii==jj) {
                hist(posts[,ii], breaks = 40, col = "black", main = "", xlab = "", ylab = "",
                     xlim = ranges[,ii], las = 2, cex.axis = cex.axis)
                if(show.lines) abline(v=truepars[ii], col = "red", lwd = 3)
                if(show.cis) abline(v=cis[,ii], col = "yellow", lwd = 3)
                mtext(parnames[ii], side = 3, line = -2, adj = .98, col = "red", cex = cex.nm)
              }else{
                smoothScatter(posts[,jj],posts[,ii], cex = cex, col = col, las = 2, cex.axis = cex.axis,
                              main = "", xlab = "", ylab = "", nrpoints = nrpoints,
                              xlim = ranges[,jj], ylim = ranges[,ii])
                if(!is.null(inits.tp)) points(inits.tp[,jj], inits.tp[,ii], cex = 1, pch = 19, col = 'purple')
                if(show.lines) abline(v=truepars[jj], col = "red", lwd = 2)
                if(show.lines) abline(h=truepars[ii], col = "red", lwd = 2)              
              }
          }
      }
    if(do.pdf | do.jpeg) dev.off()
  }
 
## get posteriors in friendly data frame
procpost <- function(d.out, nll = F, nc = 12, to.plot=T, to.plot.exp=T, ldpars, dirnm, nm) # show nll?
    {
        mcmc.d.out <- list(NA)
        d.aratio <- 0
        for(ii in 1:nc) {
                                        #  rownames(d.out[[ii]]$out) <- c(names(ldpars),'nll')
            mcmc.d.out[[ii]] <- as.mcmc(t(d.out[[ii]]$out))
            d.aratio <- d.aratio + d.out[[ii]]$aratio
            if(ii==1) { ## initialize
                inits <- d.out[[ii]]$inits
            }else{ ## append
                inits <- rbind(inits, d.out[[ii]]$inits)
            }
        }
        mcmc.d.out <- mcmc.list(mcmc.d.out) # as mcmc.list
        d.aratio <- d.aratio/nc             # acceptance ratio
        print(paste("acceptance ratio is", round(d.aratio,2)))
        parnames <- names(ldpars)
        if(nll) parnames <- c(parnames, 'nll')
        for(pp in parnames) assign(paste0(pp,'.vec'), unlist(mcmc.d.out[,pp])) ## pull out all chains into vectors
        for(pp in parnames) {
            if(pp==parnames[[1]]) {
                posts <- get(paste0(pp,'.vec'))
            }else{
                posts <- cbind(posts, get(paste0(pp,'.vec')))
            }
        }
        posts <- as.data.frame(posts)
        if(nll) colnames(posts) <- c(parnames,'nll') else colnames(posts) <- parnames
        if(nll) {
            posts <- posts[,c('bp','acute.sc','dur.ac','late.sc','dur.lt','dur.aids','nll')]
        }else{
            posts <- posts[,c('bp','acute.sc','dur.ac','late.sc','dur.lt','dur.aids')] 
        }
        exposts <- posts
        exposts[names(exposts)!='nll'] <- exp(posts[names(posts)!='nll'])
        exposts$ehm.ac <- (exposts$acute.sc-1)*exposts$dur.ac  ## excess hazard months during acute pha
        exposts$ehm.lt <- (exposts$late.sc-1)*exposts$dur.lt  ## excess hazard months during late phase
        ## excess hazard months during late & AIDS phase with RH=0
        exposts$ehm.ltaids <- (exposts$late.sc-1)*exposts$dur.lt + (0-1)*exposts$dur.aids
        edpars <- exp(ldpars)
        edpars <- c(edpars, ehm.ac = as.numeric((edpars['acute.sc']-1)*edpars['dur.ac']),
                    ehm.lt = as.numeric((edpars['late.sc']-1)*edpars['dur.lt']),
                    ehm.ltaids = as.numeric((edpars['late.sc']-1)*edpars['dur.lt'] + (0-1)*edpars['dur.aids']))
        einits <- exp(inits)
        einits <- cbind(einits, ehm.ac = as.numeric((einits[,'acute.sc']-1)*einits[,'dur.ac']),
                        ehm.lt = as.numeric((einits[,'late.sc']-1)*einits[,'dur.lt']),
                        ehm.ltaids = as.numeric((einits[,'late.sc']-1)*einits[,'dur.lt'] +(0-1)*einits[,'dur.aids']))
        parnames <- c('bp','acute.sc','dur.ac','late.sc','dur.lt','dur.aids')
        tracenames <- c('bp','ehm.ac','acute.sc','dur.ac','ehm.ltaids', 'ehm.lt','late.sc','dur.lt','dur.aids')
        ## Plotting
        if(to.plot) {
            hollsbpairs(posts[,parnames], truepars = ldpars[parnames], show.lines = T, ## plot posterior correlations after adaptive phase
                        inits.tp = inits[,parnames], file.nm = file.path(dirnm, paste0(nm,' LOG')), width = 12, height = 12,
                        cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)             
            if(to.plot.exp) {
                hollsbpairs(exposts[,tracenames], truepars = edpars[tracenames], show.lines = T, 
                            inits.tp = einits[,tracenames], file.nm = file.path(dirnm, nm),
                            width = 12, height = 12,
                            cex = 1, col = "black", nrpoints = 100, do.jpeg = T, browse=F)
            }}
        sigma <- cov.wt(posts[,names(ldpars)])$cov ## posterior covariance matrix, then plot what proposal distr is gonna look like
        print(gel <- gelman.diag(mcmc.d.out))
        gc()
        outtab <- rbind(apply(exposts[tracenames], 2, function(x) quantile(x, c(.025, .5, .975),na.rm=T)),
                        edpars[tracenames])
        rm(list=setdiff(ls(), c("posts","exposts","tracenames","mcmc.d.out","inits",'sigma','gel','outtab'))) ## remove everything but output
        gc()
        return(list(posts = posts, exposts = exposts[,tracenames], mcmc.out=mcmc.d.out, inits = inits, sigma = sigma, gel = gel,
                    outtab = outtab))
    }


####################################################################################################
## Real Wawer data
## Incident couples
## 10/23 in 1st interval (probably missing 40% loss to follow-up due to error, could really be  10/38)
## 2/13 in 2nd
## 1/7 in 3rd-4th (NOT SURE WHICH INTERVAL IT WAS, but let's make it 4th for now, since that keeps the denominator the same in both)
inct <- data.frame(int = 1:4, n = c(23,13,7,7), i = c(10,2,0,1))
inct.no.err <- data.frame(int = 1:4, n = c(round(23/.6),13,7,7), i = c(10,2,0,1))
## Prevalent couples
## 14/161 in 1st
## 9/129 in 2nd
## 10/92 in 3rd
## 3/45 in 4th
prevt <- data.frame(int = 1:4, n = c(161,129,92,45), i = c(14,9,10,3))
## Late couples
## 2/22 in 4th to last
## 9/35 in 3rd to last, 35 - (22-2) = 15 new couples followed 
## 8/31 in 2nd to last  31 - (35-9) = 5 new couples followed
## 0/13 in last before death, 23 - 31 - 8 = 0 new couple followed (23 couples were available last interval but only 13 had the live partner followed up)
## 19/51 total, though that's missing the 13, which would yield 64 late couples
latet <- data.frame(int = 4:1, n = c(22,35,31,13), i = c(2,9,8,0))
wtab.rl <- list(inct=inct, prevt=prevt, latet=latet)
wtab.rl.no.err <- list(inct=inct.no.err, prevt=prevt, latet=latet)

wtab.rlp <- within(wtab.rl, {
    inct$p <- with(inct, i/n)
    prevt$p <- with(prevt, i/n)
    latet$p <- with(latet, i/n)
    inct$ni <- with(inct, n-i)
    prevt$ni <- with(prevt, n-i)
    latet$ni <- with(latet, n-i)
})
contTabsRl <- lapply(wtab.rlp, '[', c('i','ni'))
