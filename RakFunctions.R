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


## create a retrospective cohort at interv monthly intervals from start.rak to end.rak, only keep
## couples that were observed more than once, and which were serodiscordant at some point during the
## observed time (they could go -- to ++ in one visit interval though).
rak.coh.fxn <- function(ts, dat, interv = 10, max.vis = 5, start.rak, end.rak, ## interv is interval in months between visits
                        ## ltf.prob = monthly probability of loss to follow-up, rr = +- or -+ vs -- or ++
                        ltf.prob = NA, rr.ltf.ff = 1, rr.ltf.mm = 1, rr.ltf.hh = 1, rr.ltf.d = 0, # .d is ltf when dead
                        browse = F) 
    {
        if(browse) browser()
        sttmon <- (start.rak-1900)*12       # start month in CMC
        endmon <- (end.rak-1900)*12         # end month in CMC
        ## Select all couples that were serodiscordant at some point during the cohort (don't have to worry about ss
        ## -> hh transition in one month, because in our model an individual cannot be infected & infect their
        ## partner in the same month)
        sel <- apply(ts[sttmon:endmon,], 2, function(x) { sum(x %in% c('mm','ff','hh'))>0 } )
        ts.sdc <- ts[,sel]
        dat.sdc <- dat[sel,]
        ## reduce data frame to visit months
        vis.mon <- seq(sttmon, endmon, by = interv) 
        ts.vm <- ts.sdc[vis.mon,]
        ## replace time points of couples with partners aged out of cohort with NA (>49 for f, >59 for m; note actual Rakai
        ## cohort is >59 for both, but we stick with DHS criteria, shouldn't make a difference)
        ts.vm <- apply(ts.vm, 2, function(x) { x[grepl('a.',x)] <- NA; return(x) } )
        ## replace pre-couple time points with NAs
        ts.vm <- apply(ts.vm, 2, function(x) { x[grepl('b.',x)] <- NA; return(x) } )
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
            ## create ltfp matrix
            ltfps <- ts.vm
            sers <- c('ss','mm','ff','hh')
            for(ser in sers)        ltfps[which(ts.vm==ser, arr.ind=T)] <- get(paste0('ltfp.',ser))
            ltfps[apply(ltfps, 2, function(x) grepl('d.', x))] <- ltfp.d
            ## logical matrix of ltf censoring
            ltf.log <- apply(ltfps,2, function(x) { x <- as.numeric(x); x[!is.na(x)] <- rbinom(sum(!is.na(x)), 1, x[!is.na(x)]); return(x)})
            ## first censorship
            ltf.vis <- rep(NA,ncol(ltf.log))
            ltf.vis1 <- apply(ltf.log, 2, function(x) sum(x, na.rm=T)) > 0
            ltf.vis[ltf.vis1] <- apply(ltf.log[,ltf.vis1], 2, function(x) min(which(x==1),na.rm=T))
        }else{
            ltf.vis <- NA
        }
        for(ii in 1:length(num.vis)) {  # for each couple, determine censorship
            if(!is.na(ltf.vis[ii])) { ## if lost to follow-up ever
                ts.vm[1:nrow(ts.vm) > ltf.vis[ii], ii] <-  NA ## censor all visits after loss to follow up
            }
            temp.num.vis <- sum(!is.na(ts.vm[,ii])) ## how many observations are left?
            if(temp.num.vis>max.vis) { ## if more observations than max observations
                temp.visits <- which(!is.na(ts.vm[,ii])) # which visits were observed         
                show.visits <- temp.visits[1:max.vis] # visits to show
                ts.vm[! 1:nrow(ts.vm) %in% show.visits, ii] <- NA # censor others
            }
        }           
        ## now reduce cohort data frame to all couples that were observed serodiscordant at some
        ## point during the cohort, or that were observed -- and then ++ between visit intervals.
        sel.coh <- apply(ts.vm, 2, function(x) { sum(x %in% c('mm','ff'), na.rm=T)>0 | (sum(x=='ss', na.rm=T)>0 & sum(x=='hh', na.rm=T)>0) } )
        ts.vm <- ts.vm[, sel.coh]
        dat.vm <- dat.sdc[sel.coh,]
        ## remove any couples that were only observed once because of aging out
        nage <- colSums(!is.na(ts.vm))>1
        ts.vm <- ts.vm[,nage]
        dat.vm <- dat.vm[nage,]
        ## remove any couples that were only observed once because of death (is.na clause makes sure to exclude pre-couple times)
        death1 <- apply(ts.vm, 2, function(x) { sum(!grepl('d.',x) & !is.na(x)) ==1 })
        ts.vm <- ts.vm[,!death1]
        dat.vm <- dat.vm[!death1,]      
        ## ts.vm[,sample(1:ncol(ts.vm),10)] # look at some couples
        rak.coh <- list(dat.rak = dat.vm, ts.rak = ts.vm, interv = interv)
        return(rak.coh)
    }
 
####################################################################################################
## Wawer et al. style analysis of Rakai retrospective cohort
rak.wawer <- function(rak.coh, verbose = F, browse = F)
    {
        if(browse) browser()
        ts.vm <- rak.coh$ts.rak
        dat.vm <- rak.coh$dat.rak
        interv <- rak.coh$interv
        ncpl <- ncol(ts.vm)
        ## Calculate person-months at risk for second partner in each couple group.
        ## Create line list
        rakll <- data.frame(uid = dat.vm$uid, phase = NA, pm = NA, inf = 0, pm.trunc = NA, inf.trunc = NA, excl.by.err = F)
        ## Excluded by error are couples that were observed during at least 2 visits, but with the
        ## last observed visit being serodiscordant. Based on Wawer et al.'s methods description and
        ## the fact that the exact same # of incident couples were followed for 1 interval as for 2
        ## (indicating that they excluded any just followed for 1).      
        first.obs <- apply(ts.vm, 2, function(x) min(which(!is.na(x))))
        ##################################################
        ## Incident infections
        inc.wh <- which(apply(ts.vm, 2, function(x) sum(grepl('ss',x))>0))
        rakll$phase[inc.wh] <- 'inc'
        last.sus <- rep(NA, ncpl)
        last.sus[inc.wh] <- apply(ts.vm[,inc.wh], 2, function(x) max(which(x=='ss')))
        ## those that went to ++ at some point
        inc.wh.hh <- inc.wh[which(apply(ts.vm[,inc.wh], 2, function(x) sum(grepl('hh',x))>0))]
        rakll$inf[inc.wh.hh] <- 1
        ## those that went -- to ++ in one interval
        inc.wh.hh1 <- inc.wh[ts.vm[cbind(last.sus[inc.wh]+1,inc.wh)]=='hh'] 
        rakll$pm[inc.wh.hh1] <- interv/4
        ## those that went -- to ++ but not in one interval
        inc.wh.hh2 <- inc.wh.hh[!inc.wh.hh %in% inc.wh.hh1]
        rakll$pm[inc.wh.hh2] <- apply(ts.vm[,inc.wh.hh2], 2, function(x) sum(grepl('mm',x)) + sum(grepl('ff',x))) * interv
        ## those that never went to ++
        inc.wh.nhh <- inc.wh[!inc.wh %in% inc.wh.hh]
        ## those that never went to ++ that were only observed SDC once and were consequently probably excluded from the Wawer study
        inc.wh.nhh.exl.err <- inc.wh.nhh[apply(ts.vm[,inc.wh.nhh], 2, function(x) sum(grepl('mm',x) + grepl('ff',x))==1)]
        rakll$excl.by.err[inc.wh.nhh.exl.err] <- T
        ## person months = (# times observed SDC -1)*interv + interv/2
        rakll$pm[inc.wh.nhh] <- apply(ts.vm[,inc.wh.nhh], 2, function(x) sum(grepl('mm',x) + grepl('ff',x))-1)*interv + interv/2
        ##################################################
        ## Chronic infections
       ch.wh <- which(apply(ts.vm, 2, function(x) sum(grepl('ss',x) + grepl('d.',x))==0))
        rakll$phase[ch.wh] <- 'prev'
        ## those that became ++
       ch.wh.hh <- ch.wh[which(apply(ts.vm[,ch.wh], 2, function(x) sum(grepl('hh',x))>0))]
        ##person months = (# times observed SDC -1)*interv + interv/2 for ++
        rakll$pm[ch.wh.hh] <- apply(ts.vm[,ch.wh.hh], 2, function(x) sum(grepl('mm',x) + grepl('ff',x))-1)*interv + interv/2
        rakll$inf[ch.wh.hh] <- 1
        ## those that stayed +-
        ch.wh.nhh <- ch.wh[!ch.wh %in% ch.wh.hh]
        rakll$pm[ch.wh.nhh] <- apply(ts.vm[,ch.wh.nhh], 2, function(x) sum(grepl('mm',x) + grepl('ff',x))-1)*interv
        ##################################################
        ## Late infections
        lt.wh <- which(apply(ts.vm, 2, function(x) sum(grepl('d.',x))>0))
        # **************************************************???
        ## What to do with couples that are both EARLY & LATE?? for now leave them as early only
        lt.wh <- lt.wh[!lt.wh %in% inc.wh]
        rakll$phase[lt.wh] <- 'late'
        ## those that became ++ and were seen ++ at a visit *including* if only first seen ++ at the first visit after a partner's death
        lt.wh.hh <- lt.wh[which(apply(ts.vm[,lt.wh], 2, function(x) sum(grepl('hh',x) )>0))]
        rakll$inf[lt.wh.hh] <- 1
        rakll$pm[lt.wh.hh] <- apply(ts.vm[,lt.wh.hh], 2, function(x) sum(grepl('mm',x) + grepl('ff',x))-1)*interv + interv/2
        ## those that stayed +- up until last visit
        lt.wh.nhh <- lt.wh[!lt.wh %in% lt.wh.hh]
        rakll$pm[lt.wh.nhh] <- apply(ts.vm[,lt.wh.nhh], 2, function(x) sum(grepl('mm',x) + grepl('ff',x))-1)*interv
        rakll$phase <- factor(rakll$phase)
        rakll$phase <- relevel(rakll$phase, ref = 'prev')
        ##################################################
        ## one final adjustment, Wawer compare hazards from first '5 months' post incident couples'
        ## index partner's infection to the hazards in prevalent SDCs. So we need to truncate the
        ## person-months observe in incident couples to 5 months for those observed for longer in
        ## this analysis, and also need to exclude any infections that occurred after.
        rakll$pm.trunc <- rakll$pm
        rakll$pm.trunc[rakll$phase=='inc' & rakll$pm.trunc>interv/2] <- interv/2
        rakll$inf.trunc <- rakll$inf
        rakll$inf.trunc[inc.wh.hh2] <- 0
        ####################################################################################################
        ## Poisson Regression (ignoring any source of heterogeneity, i.e. no coital acts, GUD, age, etc)
        ## assume coital acts are a function of person-months and use that as the offset
        ## formul <- formula(inf.trunc ~ offset(pm.trunc) + phase)
        ## mod <- glm(formul, family = poisson, data = rakll)
        ## summary(mod)
        ## exp(coef(mod))
        ################################################## ***************
        ## left off here, it's underestimating ARH & not running the Poisson regession... Try making
        ## a better function of interv & max.vis, it should converge to empir.arh for small interv &
        ## big max.vis
        ####################################################################################################
        ## empirical hazards
        ehz.inc <- sum(rakll$inf.trunc[inc.wh]) / sum(rakll$pm.trunc[inc.wh])
        ehz.ch <- sum(rakll$inf.trunc[ch.wh]) / sum(rakll$pm.trunc[ch.wh])
        ehz.lt <- sum(rakll$inf.trunc[lt.wh]) / sum(rakll$pm.trunc[lt.wh])
        ## empirical hazard ratio
        e.arh <- ehz.inc/ehz.ch
        e.arh
        ## empirical hazards excluding couples only observed +- at one point (& excluded by error)
        sel <- inc.wh[!inc.wh %in% which(rakll$excl.by.err)]
        ehz.inc <- sum(rakll$inf.trunc[sel]) / sum(rakll$pm.trunc[sel])
        sel <- ch.wh[!ch.wh %in% which(rakll$excl.by.err)]        
        ehz.ch <- sum(rakll$inf.trunc[sel]) / sum(rakll$pm.trunc[sel])
        sel <- lt.wh[!lt.wh %in% which(rakll$excl.by.err)]        
        ehz.lt <- sum(rakll$inf.trunc[sel]) / sum(rakll$pm.trunc[sel])
        ## empirical hazard ratio
        e.arh.err <- ehz.inc/ehz.ch
        e.arh.err
        ##################################################
        ## Add variables needed for Hollingsworth et al. style analysis     
        rakll$kk <- NA ## intervals followed (or interval of infection for late couples)
        rakll$kkt <- NA ## for late couples total intervals between first observation & index partner death
        ##########
        ## Incident Couples {ss->hh->hh} k=1; {ss->mm} k = 1; {ss->mm->mm} k =2; {ss->mm->hh->hh} k = 2:
        ## If infected: All non-susceptible visits minus all ++ visits except one.
        rakll$kk[inc.wh.hh] <- apply(ts.vm[,inc.wh.hh], 2, function(x) sum(x!='ss',na.rm=T) - (sum(x=='hh',na.rm=T)-1))
        ## If uninfected:All non-susceptible visits.        
        rakll$kk[inc.wh.nhh] <- apply(ts.vm[,inc.wh.nhh], 2, function(x) sum(x!='ss',na.rm=T))
        ## Check that calcultions are working
        if(verbose) {rnd <- sample(inc.wh,10); print(ts.vm[,rnd]); print(rakll[rnd,])}
        ##########
        ## Prevalent couples {mm->hh} k=1; {mm->hh->hh} k=1; {mm->mm->hh->hh} k = 2; {mm->mm->mm} k=2
        ## If infected: All SDC visits
        rakll$kk[ch.wh.hh] <- apply(ts.vm[,ch.wh.hh], 2, function(x) sum(x%in%c('mm','ff'), na.rm=T))
        ## If uninfected: All SDC visits - 1
        rakll$kk[ch.wh.nhh] <- apply(ts.vm[,ch.wh.nhh], 2, function(x) sum(x%in%c('mm','ff'), na.rm=T)) - 1
        ## Check that calcultions are working
        if(verbose) {rnd <- sample(ch.wh,10); print(ts.vm[,rnd]); print(rakll[rnd,])}
        ##########
        ## Late couples: more complicated because we have to account for both total follow-up
        ## intervals & interval of infection (if second partner gets infected) which may not be the
        ## same since we are also keeping track of time until index partner's death.
        ## {mm->mm->hh->d.hh} k=2 (interval of infection), kkt=3 (intervals followed before death)
        ## {mm->mm->mm->d.m} k=NA (no infection), kkt = 3
        ## {mm->d.m} k=NA, kkt=1
        ## {mm->mm->d.hh} k = 1 (last inteval before death), kkt=2
        ## If uninfected or infected: kkt = all non-dead observations
        rakll$kkt[lt.wh] <- apply(ts.vm[,lt.wh], 2, function(x) sum(x %in%c('mm','ff','hh'), na.rm=T))
        ## If infected: k = last time SDC, kkt = all non-dead observations
        rakll$kk[lt.wh.hh] <- apply(ts.vm[,lt.wh.hh], 2, function(x) min(which(grepl('d.',x))) - max(which(x %in% c('mm','ff'))))
        ## If uninfected: k = 0, because NA's screw up code later on, but this is meaningless
        rakll$kk[lt.wh.nhh] <- 0
        ## Check that calcultions are working
        if(verbose) {rnd <- sample(lt.wh,10); print(ts.vm[,rnd]); print(rakll[rnd,])}
        return(list(e.arh = e.arh, e.arh = e.arh.err, rakll = rakll))  
    }


####################################################################################################
## Hollingsworth et al. style analysis of Rakai retrospective cohort
rak.holl <- function(rak.coh, incl.exl.by.err = F, browse = F)
    {
        if(browse) browser()

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
## NOTE that while we replicate Hollingsworth et al.'s analysis, there is a Monte Hall like problem
## here: the probability of the index partner being infected at each time point T is not uniform on
## the interval and depends on the second partner's infection status. Couples going from -- to ++
## are more likely to have incident infections occuring earlier in the inteveral. Couples going from
## -- to +- are likely to have incident infections occuring later in the interval.
## **************************************************
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
##############################
## Need to integrate the conditional probability of 2ndary seroconversion (or not) in each interval
## over all possible timing of index partner seroconversion in the first interval.

####################################################################################################
## Incident Couples: k = interval since -- (1st, 2nd, 3rd, 4th, etc..)
ucp.inc <- function(kk,inf,dpars) {
  for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
  ## numerical integrals (to check calculus)
  ucpS <- integrate(Vectorize(cp.inc.1,'tt'), lower=0, upper=interv, dpars = dpars)$val / interv ## 
  ucpNS <- 1-ucpS
  if(kk >1) { ## otherwise, if k = 2,3,...
    for(jj in 2:kk) {
      ## numerical integrals
      ucpS <- ucpNS * integrate(Vectorize(cp.inc.k,'tt'), lower=0, upper=interv, kk=jj, dpars = dpars)$val / interv
      ucpNS <- ucpNS * (1- integrate(Vectorize(cp.inc.k,'tt'), lower=0, upper=interv, kk=jj, dpars = dpars)$val / interv)
    }
  }
  if(is.na(log(ucpS))) browser()
  if(is.na(log(ucpNS))) browser()  
  if(inf==1) return(as.numeric(ucpS)) #list(ucpA,ucpS))
  if(inf==0) return(as.numeric(ucpNS)) #list(ucpNA,ucpNS))
}
## Test code to make sure function is working right
## ucp.inc(1,0,dpars) + ucp.inc(1,1,dpars)                ## should be 1
## ucp.inc(2,1,dpars) + ucp.inc(2,0,dpars) - ucp.inc(1,0,dpars)  ## should be 0 (w/ rounding error,dpars)
## ucp.inc(3,1,dpars) + ucp.inc(3,0,dpars) - ucp.inc(2,0,dpars)  ## should be 0 (w/ rounding error,dpars)
## ucp.inc(4,1,dpars) + ucp.inc(4,0,dpars) - ucp.inc(3,0,dpars)  ## should be 0 (w/ rounding error)

####################################################################################################
##############################
## Prevalent Couples (assumed to be the same in all intervals)
##############################
ucp.prev <- function(kk,inf,dpars) {  ## probability of seroconversion in interval k, given not in previous intervals
  if(inf==1) ret <- cp.prev.k(dpars=dpars) * (1-cp.prev.k(dpars=dpars))^(kk-1)
  if(inf==0) ret <- (1-cp.prev.k(dpars=dpars))^kk
  return(as.numeric(ret))
}
## Test code to make sure function is working right
## ucp.prev(1,0) + ucp.prev(1,1)                ## should be 1
## ucp.prev(2,1) + ucp.prev(2,0) - ucp.prev(1,0)  ## should be 0 (w/ rounding error)
## ucp.prev(3,1) + ucp.prev(3,0) - ucp.prev(2,0)  ## should be 0 (w/ rounding error)
####################################################################################################
##############################
## Late Couples (assumed to be the same in all intervals)
##############################
## probability of seroconversion in interval kk BEFORE death, given not in previous
## intervals, **when observed for a total of kkt intervals before death**

ucp.lt <- function(kk, kkt=NA, inf,dpars) {
    if(inf==1) {
        if(kk>kkt) stop('error: kk>kkt')
        ## if kkt = 1, then it's the same as the conditional probability in the last month
        if(kkt==1) { ## only need to integrate up until (interv-dur.aids) because if dying more than
            ## dur.aids before last interval, there's no transmission in this interval
            max.int <- max(0,(interv-dpars['dur.aids']))
            ucp <- try(integrate(Vectorize(cp.lt.1,'td'), lower=0, upper=max.int, dpars=dpars)$val / interv,
                       silent = T)
            w.step <- 0
            while(inherits(ucp ,'try-error')){ ## if there's a divergent integral, then integrate over a
                ## little smaller of an interval to avoid the problematic area
                if(w.step>30) ucp <- 0 # give up on this integral eventually & just reject this proposal
                max.int <- max.int *.98
                ucp <- try(integrate(Vectorize(cp.lt.1,'td'), lower=0, upper=max.int, dpars=dpars)$val / interv,
                           silent = T)
                w.step <- w.step + 1
            }
        }
        if(kkt >1) { ## otherwise, if observed more than one interval before death
            if(kk==kkt) { ## if infected in first interval of observation
                ## min/max trickery is to avoid divergent integrals occurring when most of the function takes the value of 0 (due to aids phase)
                max.int <- min(10,max(0,(interv*kk-dpars['dur.aids'])))
                ucp <- try(integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=kk, dpars=dpars)$val / interv, silent = T)
                w.step <- 0
                while(inherits(ucp ,'try-error')){ ## if there's a divergent integral, then integrate over smaller region
                    if(w.step>30) ucp <- 0 # give up on this integral eventually & just reject this proposal
                    max.int <- max.int*.98
                    ucp <- try(integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=kk, dpars=dpars)$val / interv, silent = T)
                    w.step <- w.step + 1
                }              
            }else{ ## if not infected in at least one interval before infected (kk<kkt)
                ucp <- 1 ## initialize ucp for below loop
                for(jj in kkt:(kk+1)) { ## calculate probability of not being infected between
                    ## visits kkt to (kk+1) before death (where 1 = visit at dath)
                    max.int <- min(10,max(0,(interv*jj-dpars['dur.aids'])))
                    ucp <- try(ucp * (1 - integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=jj, dpars=dpars)$val /interv), silent = T)
                    w.step <- 0
                    while(inherits(ucp, 'try-error')){
                        if(w.step>30) ucp <- 0 # give up on this integral eventually & just reject this proposal
                        max.int <- max.int*.98
                        ucp <- try(ucp * (1 - integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=jj, dpars=dpars)$val /interv), silent = T)
                        w.step <- w.step + 1
                    }
                }
                ## times the probability of being infected in the kk-th interval
                max.int <- max.int*.98
                ucp <- try(ucp * integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=kk, dpars=dpars)$val / interv, silent = T)
                w.step <- 0
                while(inherits(ucp, 'try-error')){
                    if(w.step>30) ucp <- 0 # give up on this integral eventually & just reject this proposal
                    max.int <- max.int*.98
                    ucp <- try(ucp * integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=kk, dpars=dpars)$val / interv, silent = T)
                    w.step <- w.step + 1
                }                
            }
        }
        return(as.numeric(ucp))
    }
    if(inf==0) {
        max.int <- max(0,(interv-dpars['dur.aids']))
        ucpN <- try(1 - integrate(Vectorize(cp.lt.1,'td'), lower=0, upper=max.int, dpars=dpars)$val / interv, silent = T)
        w.step <- 0
        while(inherits(ucpN ,'try-error')){ ## if there's a divergent integral, then integrate over a
            if(w.step>30) ucpN <- 0 # give up on this integral eventually & just reject this proposal
            ## little smaller of an interval to avoid the problematic area
            max.int <- max.int*.98
            ucpN <- try(1 - integrate(Vectorize(cp.lt.1,'td'), lower=0, upper=max.int, dpars=dpars)$val / interv, silent = T)
            w.step <- w.step+1
        }
        if(kk >1) { ## otherwise, if observed more than one interval before death, multiply this by not being observed any of the other kk-2 intervals before death
            for(jj in kk:2) { ## calculate probability of not being infected between
                ## intervals kk to last before death (last is interval 1 before death)
                max.int <- min(10,max(0,(interv*jj-dpars['dur.aids'])))
                ucpN <- try(ucpN * (1 - integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=jj, dpars=dpars)$val / interv), silent=T)
                w.step <- 0
                while(inherits(ucpN, 'try-error')){
                    if(w.step>30) ucpN <- 0 # give up on this integral eventually & just reject this proposal
                    max.int <- max.int*.98
                    ucpN <- try(ucpN * (1 - integrate(Vectorize(cp.lt.k,'td'), lower=0, upper=max.int, kk=jj, dpars=dpars)$val / interv), silent=T)
                    w.step <- w.step+1
                }
            }
        }
        return(as.numeric(ucpN))        
    }
}

####################################################################################################
####################################################################################################
## Test code to make sure function is working right
## ucp.lt(1,1,1,dpars) + ucp.lt(1,NA,0,dpars)                ## should be 1 {mm->d.mm} or {mm->d.hh}
## ucp.lt(2,2,1,dpars) + ucp.lt(1,2,1,dpars) + ucp.lt(2,NA,0,dpars) ## should be 1 {mm->hh->d.hh} or {mm->mm->d.hh} or {mm->mm->d.mm}
## ucp.lt(3,3,1,dpars) + ucp.lt(2,3,1,dpars) + ucp.lt(1,3,1,dpars) + ucp.lt(3,NA,0,dpars) ## similarly, all possible options for 3 intervals...
 
## wrapper for mle()
holl.lik.mle <-  function(acute.sc, late.sc, bp, dur.ac, dur.lt, dur.aids,
                          rakll, excl.by.err, verbose = T, browse = F)
  {
    ldpars <- c(acute.sc = acute.sc, late.sc = late.sc, bp = bp,
                dur.ac = dur.ac, dur.lt = dur.lt, dur.aids = dur.aids)
    nll <- holl.lik(ldpars = ldpars, rakll = rakll, excl.by.err = excl.by.err, verbose = verbose, browse = browse)
    return(nll)
  }

####################################################################################################
## Calculate Hollingsworth et al. style likelihood
holl.lik <- function(ldpars, rakll, # dpars has disease progression/infectiousnes parameters
                     excl.by.err = F,         ## exclude {ss->mm/ff} couples as in Wawer?
                     verbose = T, browse = F)
  {
    if(browse) browser()
    ## exponentiate ldpars
    dpars <- exp(ldpars)
    ## Excluding incident couples seen serodiscordant once & then never again as in Wawer 2005?
    if(excl.by.err) rakll <- rakll[!rakll$excl.by.err,]
    ## Create data frame for Hollingsworth style analysis
    ## Incident Couples
    inc.kkinf <- xtabs(~kk + inf, data = rakll, subset = phase=='inc')
    inc.nll <- 0
    ## for incident couples observed jj intervals, calculate negative log likelihood of that
    ## many couples seroconverting or not seroconverting in exactly that interval.
    if(verbose) print('adding up incident nlls')
    for(jj in 1:nrow(inc.kkinf)) { 
      inc.nll <- inc.nll - inc.kkinf[jj,'0']*log(ucp.inc(jj,0, dpars=dpars)) -
        inc.kkinf[jj,'1']*log(ucp.inc(jj,1, dpars=dpars))
      if(verbose) print(inc.nll)
    }
    ## Prevalent Couples
    if(verbose) print('adding up prevalent nlls')
    prev.kkinf <- xtabs(~kk + inf, data = rakll, subset = phase=='prev')
    prev.nll <- 0
    ## for prevalent couples observed jj intervals, calculate negative log likelihood of that
    ## many couples seroconverting or not seroconverting in exactly that interval.
    for(jj in 1:nrow(prev.kkinf)) { 
      prev.nll <- prev.nll - prev.kkinf[jj,'0']*log(ucp.prev(jj,0, dpars=dpars)) -
        prev.kkinf[jj,'1']*log(ucp.prev(jj,1, dpars=dpars))
      if(verbose) print(prev.nll)
    }
    ## Late Couples
    if(verbose) print('adding up late nlls')
    lt.kkinf <- xtabs(~kk + kkt + inf, data = rakll, subset = phase=='late')
    lt.nll <- 0
    ## late couples without infected 2ndary partners observed for kk intervals before death (called kkt in rakll)
                                        #browser()
    for(jj in 1:ncol(as.matrix(lt.kkinf[,,'0']))) {
      lt.nll <- lt.nll - lt.kkinf['0',jj,'0'] * log(ucp.lt(kk = as.numeric(colnames(lt.kkinf)[jj]), kkt=NA, inf = 0, dpars=dpars))
      if(verbose) print(lt.nll)
    }
    ## late couples with infected 2ndary partners observed kkt intervals before death that with
    ## 2ndary seroconversion occuring in kk-th interval before death
    wh.row <- which(rownames(lt.kkinf)!='0')
    for(ii in wh.row) { ## values of kk
      for(jj in 1:ncol(as.matrix(lt.kkinf[,,'1']))) { ## values of kkt
        kk.temp <- as.numeric(rownames(lt.kkinf)[ii])
        kkt.temp <- as.numeric(colnames(lt.kkinf)[jj])
        if(! kk.temp > kkt.temp ) { ## if !kk>kkt
          lt.nll <- lt.nll - lt.kkinf[ii,jj,'1'] * log(ucp.lt(kk = kk.temp, kkt=kkt.temp, inf = 1, dpars=dpars))
          ## if(is.na(lt.nll))  print('lt.nll=NA'); browser()
          ## if(lt.nll==Inf) print('lt.nll=Inf'); browser()
          if(verbose) print(lt.nll)
        }
      }
    }
    nll <- inc.nll + prev.nll + lt.nll
    nll <- as.numeric(nll)
    if(verbose) print('total nll:')
    return(nll)
  }

## function to calculate attributable chronic-infection-months due to elevated acute phase infectivity
arr <- function(ldp) as.numeric((exp(ldp['acute.sc'])-1)*exp(ldp['dur.ac']))

## Simulate dat using exactly Hollingsworth's model to make sure we're fitting it correctly
holl.mod <- function(i.n = 50, p.n = 50, l.n = 50, interv = 10,  max.vis = 4, dpars, verbose=F, browse = F) {
  for(nm in names(dpars))      assign(nm, dpars[nm]) ## loading 'bp'
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
      inf <- rbinom(1,1,cp) ## Bernoulli infection variable
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
                        rakdat, excl.by.err,
                        fixed = NULL, 
                        multiv = F, covar = NULL, # if multiv, sample from multivariate distribution (calculated during adaptive phase)
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
    lprob.cur <- holl.lik(pars, rakdat, excl.by.err = excl.by.err, verbose = F, browse = F)
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
        if(vv==3311) browser()
        lprob.prop <- holl.lik(pars.prop, rakdat, excl.by.err = excl.by.err, verbose = F, browse = F)
        if(verbose2) print(lprob.prop)
        ## holl.lik gives -log(L), so negative of that gives the log probability (assuming flat improper priors)
        lmh <- -lprob.prop + lprob.cur       # log Metropolis-Hastings ratio
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

holl.wrp <- function(seed=1, sd.props, rakdat, excl.by.err, force.inits=NULL,
                multiv=F, covar=NULL, 
                verbose = T, verbose2 = F, tell = 50, browse = F,
                niter, nthin, nburn)
  { ## new version of wrp needs to save progress for long chains (WA)
    if(length(seed)>1) browse <- F ## can't browse within mclapply
    if(length(force.inits)==0) {
        inits.temp <- holl.init.fxn(seed = seed) # initial conditions different for each seed
    }else{
        inits.temp <- jitter(force.inits, a = .5)
    }
    hollsampler(sd.props = sd.props, inits = inits.temp, rakdat = sim, excl.by.err = F,
                multiv = multiv, covar = covar, 
                verbose = verbose, verbose2 = verbose2, tell = tell, seed = seed, 
                niter = niter, nthin = nthin, nburn = nburn, browse=browse)
  }



## Plot posterior pairwise-correlations & histograms
hollsbpairs <- function(posts, truepars = NULL, file.nm, width = 10, height = 10, show.lines = T,
                        cex = 1, col = "black", nrpoints = 200, do.pdf = F, do.jpeg = T, ranges = NULL,
                        cex.axis = 1.5, cex.nm = 2.5, greek = T, show.cis = T, browse = F)
  {
    if(browse) browser()
    if(do.pdf) pdf(paste(file.nm, ".pdf", sep=""), width = width, height = height)
    if(do.jpeg) jpeg(paste(file.nm, ".jpeg", sep=""), width = width*100, height = height*100)
    if(length(ranges)==0)   ranges <- apply(rbind(posts,truepars), 2, range)
    cis <- apply(rbind(posts,truepars), 2, function(x) quantile(x, c(.025, .975)))
    par(mar=rep(3,4),oma=rep(2,4),mfrow=rep(ncol(posts),2))
    parnames <- colnames(posts)
    for(ii in 1:ncol(posts))
      {
        for(jj in 1:ncol(posts))
          {
            if(ii==jj)
              {
                hist(posts[,ii], breaks = 40, col = "black", main = "", xlab = "", ylab = "",
                     xlim = ranges[,ii], las = 2, cex.axis = cex.axis)
                if(show.lines) abline(v=truepars[ii], col = "red", lwd = 3)
                if(show.cis) abline(v=cis[,ii], col = "yellow", lwd = 3)
                mtext(parnames[ii], side = 3, line = -2, adj = .98, col = "red", cex = cex.nm)
              }else{
                smoothScatter(posts[,jj],posts[,ii], cex = cex, col = col, las = 2, cex.axis = cex.axis,
                              main = "", xlab = "", ylab = "", nrpoints = nrpoints,
                              xlim = ranges[,jj], ylim = ranges[,ii])
                if(show.lines) abline(v=truepars[jj], col = "red", lwd = 2)
                if(show.lines) abline(h=truepars[ii], col = "red", lwd = 2)              
              }
          }
      }
    if(do.pdf | do.jpeg) dev.off()
  }

## get posteriors in friendly data frame
procpost <- function(d.out)
    {
        mcmc.d.out <- list(NA)
        d.aratio <- 0
        for(ii in 1:nc) {
                                        #  rownames(d.out[[ii]]$out) <- c(names(ldpars),'nll')
            mcmc.d.out[[ii]] <- as.mcmc(t(d.out[[ii]]$out))
            d.aratio <- d.aratio + d.out[[ii]]$aratio
            if(ii==1) { ## initialize
                init.adapt <- d.out[[ii]]$inits
            }else{ ## append
                init.adapt <- rbind(init.adapt, d.out[[ii]]$inits)
            }
        }
        mcmc.d.out <- mcmc.list(mcmc.d.out) # as mcmc.list
        d.aratio <- d.aratio/nc             # acceptance ratio
        print(paste("adaptive phase aratio is", round(d.aratio,2)))
        parnames <- names(ldpars)
        for(pp in parnames) assign(paste0(pp,'.vec'), unlist(mcmc.d.out[,pp])) ## pull out all chains into vectors
        for(pp in parnames) {
            if(pp==parnames[[1]]) {
                posts <- get(paste0(pp,'.vec'))
            }else{
                posts <- cbind(posts, get(paste0(pp,'.vec')))
            }
        }
        colnames(posts) <- parnames
        posts
    }
