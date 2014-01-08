library(graphics); library(abind); library(mvtnorm); library(multicore); library(lme4) # load necessary libraries


## Age-at-seroconversion dependent Weibull survival times fit to
## CASCADE 2000 data. See Bellan et al. (2013) Supp Info for details.
ageweib <- function(age, death = T)
  {
    if(death)
      {
        shp <- 2.3                      # Weibull shape parameter
        scl <- 2000/shp/(age/12)^.53    # Weibull scale parameter
        return(round(rweibull(length(age), shape = shp, scale = scl),0))
      }else{
        return(200*12)               # if there's no death arbitrarily return 200 year survival time
      }
  }

##  Event-driven simulation function
event.fn <- function(pars, dat, browse = F, # transmission coefficients to use for simulation, DHS couples data, debug
                     death = T, # have aids death in model?
                     acute.sc = 26, # Acute to chronic phase relative hazard (2 months after infection)
                     late.sc = 1, # Late to chronic phase relative hazard (20-11 months before death)
                     aids.sc = 0, # AIDS to chronic phase relative hazard (2 months after infection) (10-0 months before death)
                     ##  route-specific heterogeneity
                     het.b = F, het.b.sd = 0, het.b.cor = 0, # pre-couple
                     het.e = F, het.e.sd = 0, het.e.cor = 0, # extra-couple
                     het.p = F, het.p.sd = 0, het.p.cor = 0, # within-couple
                     ## genetic heterogeneity (create individual lognormal risk deviate determining
                     ## an individual susceptibility that is consistent throughout life)
                     het.gen = F,     # do it?                                                                           
                     het.gen.sd = 1,  # standard deviation of genetic heterogeneity                                     
                     het.gen.cor = 0, # inter-partner correlation of individual risk deviate for genetic heterogeneity
                     ## individual suscepbitility constant through life is either HIGH or LOW risk,
                     ## proportion high risk (phigh) are at rrhigh times higher risk
                     hilo = F, phigh.m = .2, phigh.f = .2, rrhigh.m = 10, rrhigh.f = 10,
                     ## next 3 are same as previous 3, but for behavioral heterogeneity (risk
                     ## deviate amplifies only pre- and extra-couple transmission susceptibility)
                     het.beh = F, het.beh.sd = 1, het.beh.cor = 0, 
                     scale.by.sd = T, # adjust beta means to keep geometric mean constant with increasing het
                     scale.adj = 1,   # adjust them arbitrarily
                     end.at.int = T, # end infections at interviews (otherwise go til first partner is 60)
                     nc = 12,        # number of cores
                     vfreq = 200)    # how often to show progress
  {
    if(hilo & (het.gen | het.b | het.e | het.p | het.beh)) stop("Can't have both continuous and discrete heterogeneity")
    K <- nrow(dat)
    mser <- rep(0, nrow(dat))    ## serostatus at end
    fser <- rep(0, nrow(dat))
    mdod <- rep(NA, nrow(dat))    ## date of death (CMC)
    fdod <- rep(NA, nrow(dat))
    mdoi <- rep(NA, nrow(dat))    ## date of infection (CMC)
    fdoi <- rep(NA, nrow(dat))
    mcoi <- factor(rep(NA, nrow(dat)), levels = c("b","e","p","-"))    ## route of infection
    fcoi <- factor(rep(NA, nrow(dat)), levels = c("b","e","p","-"))
    hpars <- pars                       # backup transmission coefficients
    ##################################################################
    ## For each type of heterogeneity, sample random lognormal deviates
    for(het.t in c('b','e','p','gen','beh')) {     # pre-, -extra, -within, genetic (bep), behavioral (be)
        het.log <- get(paste('het.',het.t,sep='')) ## temporary logical indicator of whether this type of heterogeneity is being simulated
        het.sd <- get(paste('het.',het.t,'.sd',sep='')) # corresponding standard deviation
        het.cor <- get(paste('het.',het.t,'.cor',sep='')) # corresponding inter-partner correlation
        if(het.log) { #  include this type of heterogeneity?
          if(het.cor==0) {              #   uncorrelated
            assign(paste('m.het.',het.t,sep=''), exp( rnorm(K, mean = 0, sd = het.sd) ))
            assign(paste('f.het.',het.t,sep=''), exp( rnorm(K, mean = 0, sd = het.sd) ))            
          }else{                        #  correlated between partners
            b.corSig <- matrix(c(het.sd,rep(het.cor*het.sd,2),het.sd),2,2)
            het.both <- rmvnorm(K, mean = rep(0,2), sigma = b.corSig)
            assign(paste('m.het.',het.t,sep=''), exp(het.both[,1]))
            assign(paste('f.het.',het.t,sep=''), exp(het.both[,2]))            
          }
          if(het.t %in% c('b','e','p')) { # which hazards to rescale to keep geometric mean of risk deviates = 1
              het.hazs <- paste0('b',c('m','f'),het.t)
          }else{
              if(het.t=='gen') het.hazs <- hazs
              if(het.t=='beh') het.hazs <- hazs[1:4]
          }
          if(scale.by.sd) { ##  keep geometric mean transmission coefficient the same
            hpars[het.hazs] <-  hpars[het.hazs] / exp(het.sd^2 / 2)  # correct for increased geometric mean w/ hetogeneity
          }else{
            hpars[het.hazs] <-  hpars[het.hazs] / scale.adj # scale arbitrarily
          }
        }else{ ## no  pre-couple-specific heterogeneity
          assign(paste('m.het.',het.t,sep=''), rep(1, K))
          assign(paste('f.het.',het.t,sep=''), rep(1, K))          
        }
      } # end heterogeneity risk deviate assignment loop
    ## Discrete (binary) heterogeneity: High risk & low risk groups
    m.het.hilo <- rep(1,K)
    f.het.hilo <- rep(1,K)
    if(hilo) {
      m.high <- rbinom(K, 1, phigh.m)
      f.high <- rbinom(K, 1, phigh.f)
      m.het.hilo[m.high] <- rrhigh.m
      f.het.hilo[f.high] <- rrhigh.f      
    }
    ## make data frame for storing output
    dat <- data.frame(dat, mser, fser, mdoi, fdoi, mdod, fdod, mcoi, fcoi, m.het.gen, f.het.gen, m.het.beh, f.het.beh,
                      m.het.b, f.het.b, # pre-couple   
                      m.het.e, f.het.e, # extra-couple 
                      m.het.p, f.het.p, # within-couple
                      m.het.hilo, f.het.hilo) # genetic, binary (hi vs lo)
    ## track if infections from partner were due to acute phase infectiousness
    dat$mcoi.phase <- NA
    dat$fcoi.phase <- NA
    if(browse) browser()
    ## Break couples into nc batches to split between cores. Divide between cores by couple
    ## formation date strata since early couples will always take longest.
    breaks <- rep(1:nc, length.out = nrow(dat))
    ## Call couple loop function below in multicore
    multi.out <- mclapply(1:nc, cloop, dat = dat, breaks = breaks, vfreq = vfreq, death = death, acute.sc = acute.sc, late.sc = late.sc, aids.sc = aids.sc,
                          pars = hpars, browse=F)
    temp <- multi.out[[1]]    ## Combine mclapply output from each core back into data frame
    if(nc>1) {
        for(ii in 2:nc)     temp <- rbind(temp, multi.out[[ii]])
      }
    dat <- temp
    dat <- dat[order(dat$uid),]     #  reorder by unique couple identifier
    ## End time is the minimum of dates of death and interview time.
    dat$tend <- ord(dat[, c("tint","mdod","fdod")], 1)
    ## create variable for month couple ages out of DHS couple cohort, if applicable
    dmage <- dat$tint - (dat$mage - 60*12) # males age out at age 60
    dmage[dmage<0] <- NA
    dfage <- dat$tint - (dat$fage - 50*12) # females age out at age 50
    dfage[dfage<0] <- NA
    dage <- apply(cbind(dmage,dfage),1, function(x) min(x, na.rm = T)) # find earliest date at which a partner ages out
    ## create taend which is endtime but also incorporating aging out
    taend <- apply(cbind(dat$tend,dage),1, function(x) min(x, na.rm = T)) 
    dat <- data.frame(dat, dmage = dmage, dfage = dfage, dage = dage, taend = taend) # append new variables to data frame
    ## alive at time of aging out/interview/AIDS death (ie did they die or did their partner)
    dat$malive <- rep(T, nrow(dat))
    dat$malive[dat$mdod <= dat$taend] <- F
    dat$falive <- rep(T, nrow(dat))
    dat$falive[dat$fdod <= dat$taend] <- F
    dat$alive <- dat$malive & dat$falive
    dat$mcoi[is.na(dat$mcoi)] <- "-"    # replace NA with -  for cause of infections
    dat$fcoi[is.na(dat$fcoi)] <- "-"
    dat$ser[dat$mser == 1 & dat$fser == 1] <- 1 # create ser (couple serostatus variable)
    dat$ser[dat$mser == 1 & dat$fser == 0] <- 2
    dat$ser[dat$mser == 0 & dat$fser == 1] <- 3
    dat$ser[dat$mser == 0 & dat$fser == 0] <- 4
######################################################################
    ## Truncated survival times
######################################################################
    ## Survival times of couple serostatus, as truncated by both death and
    ## interviews.
######################################################################
    ## Couples that were ever mSDC: Would have been either
    ## been CCC at end,  not have died before forming, mdoi < fdoi, and fdoi < tmar.
    msdc <- (dat$taend > dat$tmar & dat$ser==1 & dat$mdoi < dat$fdoi & dat$fdoi > dat$tmar)
    ## or mSDC at end,  not have died before forming.
    msdc <- msdc | (dat$taend > dat$tmar & dat$ser==2)
    msdc[is.na(msdc)] <- F
###################################################################### 
    ## same for females
    fsdc <- (dat$taend > dat$tmar & dat$ser==1 & dat$fdoi < dat$mdoi & dat$mdoi > dat$tmar)
    ## or mSDC at end,  not have died before forming.
    fsdc <- fsdc | (dat$taend > dat$tmar & dat$ser==3)
    fsdc[is.na(fsdc)] <- F
###################################################################### 
    ## Couples that were ever CCC: Would have been CC at
    ## the end (interview or time of death) and not have died before forming.
    ccc <- dat$ser == 1 & dat$taend > dat$tmar
    ## mSDC survival times: Time between couple formation/male infection and female infection/end time.
    tmsdc <- rep(NA, nrow(dat))
    tmsdc[msdc & ccc] <- dat$fdoi[msdc & ccc] -ord(dat[msdc & ccc, c("tmar","mdoi")], 2)
    tmsdc[msdc & !ccc] <- dat$taend[msdc & !ccc] -ord(dat[msdc & !ccc, c("tmar","mdoi")], 2)
    ## for female
    tfsdc <- rep(NA, nrow(dat))
    tfsdc[fsdc & ccc] <- dat$mdoi[fsdc & ccc] -ord(dat[fsdc & ccc, c("tmar","fdoi")], 2)
    tfsdc[fsdc & !ccc] <- dat$taend[fsdc & !ccc] -ord(dat[fsdc & !ccc, c("tmar","fdoi")], 2)
    ## CCC survival time: Time between couple formation/2nd infection and endtime
    tccc <- rep(NA, nrow(dat))
    tccc[ccc] <- dat$taend[ccc] - ord(dat[ccc, c("tmar","fdoi","mdoi")], 3)
    ## add to data frame.
    dat$tmsdc <- tmsdc
    dat$tfsdc <- tfsdc
    dat$tccc <- tccc
    dat$msdc <- msdc
    dat$fsdc <- fsdc
    dat$ccc <- ccc    
    return(dat)
  }

###################################################################### 
## couple loop for event.fn (to allow parallelization
## batch ii does couples in the range breaks[ii,1:2]
######################################################################
cloop <- function(batch, dat, pars, breaks, vfreq, browse = F, death, acute.sc, late.sc, aids.sc = 0)
  {                                     
    set.seed(batch)                     # in case doing nonparametric inflation, we want each inflated couple to have different trajectories
    if(browse) browser()
    dat <- dat[breaks==batch,]          # only simulate couples in the batch
    ## attach parameters
    bmb <- as.numeric(pars["bmb"])
    bfb <- as.numeric(pars["bfb"])
    bme <- as.numeric(pars["bme"])
    bfe <- as.numeric(pars["bfe"])
    bmp <- as.numeric(pars["bmp"])
    bfp <- as.numeric(pars["bfp"])
    for(ii in 1:nrow(dat))              # loop over couples
      {
        if(ii %% vfreq == 0) print(paste("On couple", ii, "of", nrow(dat))) # show progress
        if(substitute & exists(as.character(substitute(s.epic.ind)))) { # if substitution analysis and there is a substituted epidemic curve (otherwise not)
          epic.ind.temp <- s.epic.ind   # substitued epidemic curve
        }else{
          epic.ind.temp <- dat$epic.ind[ii]                                   # choose epidemic curve
        }
        ## keep simulating until date of death > date of couple
        ## formation, i.e. condition on couple forming.
        finished <- FALSE
        while(!finished) 
          {
            ## male before marriage: get bernoulli probability of infection in each month of sexual
            ## activity before marriage
            m.inf.bef <- rbinom(dat$tmar[ii] - dat$tms[ii],1, # hets: b,beh,gen
                                prob = 1 - exp(-bmb*dat$m.het.hilo[ii]*dat$m.het.b[ii]*dat$m.het.beh[ii]*dat$m.het.gen[ii]*epicf[dat$tms[ii]:(dat$tmar[ii]-1), epic.ind.temp]))
            if(sum(m.inf.bef)>0)            ## if he gets infected in 1 or more months
              {
                ## change serostatus to HIV+ & use earliest infection as date of infection. Then
                ## figure out date of death from age-at-seroconversion dependent Weibull survival
                ## times.
                temp.mdoi <- dat$tms[ii] + min(which(m.inf.bef==1)) - 1 # date of infection
                temp.mdod <- temp.mdoi + ageweib(dat$mage[ii] - (dat$tint[ii] - temp.mdoi), death = death) # date of death
                if(temp.mdod > dat$tmar[ii]) # if death is after couple formation, set values & finished
                  {
                    dat$mser[ii] <- 1         # sero-positive
                    dat$mdoi[ii] <- temp.mdoi # dat$tms[ii] + min(which(m.inf.bef==1)) - 1
                    dat$mdod[ii] <- temp.mdod # mdoi[ii] + ageweib(dat$mage[ii] - (dat$tint[ii] - mdoi[ii]), death = death)
                    dat$mcoi[ii] <- "b"       #  before couple formation (pre-couple)
                    finished <- TRUE # if death is after couple formation, finished
                  }
              }else{                # if no infections, end while loop
                finished <- TRUE
              }
          }
        ## keep simulting until date of death > date of couple
        ## formation, i.e. condition on couple forming.
        finished <- FALSE
        while(!finished) 
          {
            ## female before marriage: get bernoulli probability of
            ## infection in each month of sexual activity before
            ## marriage
            f.inf.bef <- rbinom(dat$tmar[ii] - dat$tfs[ii],1, # hets: b,beh,gen
                                prob = 1 - exp(-bfb*dat$f.het.hilo[ii]*dat$f.het.b[ii]*dat$f.het.gen[ii]*dat$f.het.beh[ii]*epicm[dat$tfs[ii]:(dat$tmar[ii]-1), epic.ind.temp]))
            if(sum(f.inf.bef)>0) { ## if she gets infected in 1 or more months
              ## change serostatus to HIV+ & use earliest infection as date of infection. Then
              ## figure out date of death from age-at-seroconversion dependent Weibull survival
              ## times.
                temp.fdoi <- dat$tfs[ii] + min(which(f.inf.bef==1)) - 1 # date of infection
                temp.fdod <- temp.fdoi + ageweib(dat$fage[ii] - (dat$tint[ii] - temp.fdoi), death = death) # date of death
                if(temp.fdod > dat$tmar[ii]) # if death is after couple formation, set values & finished
                  {
                    dat$fser[ii] <- 1   # seropositive
                    dat$fdoi[ii] <- temp.fdoi 
                    dat$fdod[ii] <- temp.fdod 
                    dat$fcoi[ii] <- "b"
                    finished <- TRUE # if death is after couple formation, finished
                  }
              }else{                # if no infections, end while loop
                finished <- TRUE
              }
          }
        if(dat$mser[ii] + dat$fser[ii] < 2) { # if both aren't infected already by date of couple
                                            # formation, simulate transmission during couple
                                            # duration
          ## only simulate up to interview time or date of death, whichever is earliest
            end <- min(dat$mdod[ii], dat$fdod[ii], dat$tint[ii], na.rm = TRUE) 
            tt <- dat$tmar[ii] ##  couple formation date
            ## While both are alive and both aren't infected, run infection loop up until they age out of cohort (f<50yrs, m<60yrs).
            while((tt <= end) & (dat$mser[ii] + dat$fser[ii] < 2) & (dat$mage[ii] - (dat$tint[ii]-tt)) < (60*12) & (dat$fage[ii] - (dat$tint[ii]-tt)) < (50*12))
              {
                ## infections from partners, length 2 vector,
                ## bernoulli, c(male inf by part, female inf by part);
                ## if(((tt - dat$mdoi[ii]) <= 2 & dat$mser[ii]) | ((tt - dat$fdoi[ii]) <= 2& dat$fser[ii])) browser()
                ###################################################################### 
                ## ACUTE phase scalars
                if(dat$mser[ii]) {
                    m.ac.sc <- ifelse((tt - dat$mdoi[ii]) <= 2, acute.sc, 1) # male acutely infected for 2 months after infection
                  }else{
                    m.ac.sc <- 1
                  }
                if(dat$fser[ii]) {
                    f.ac.sc <- ifelse((tt - dat$fdoi[ii]) <= 2, acute.sc, 1) # female acutely infected for 2 months after infection
                  }else{
                    f.ac.sc <- 1
                  }
                ###################################################################### 
                ## LATE phase scalars
                if(dat$mser[ii]) {
                    m.lt.sc <- ifelse(tt > (dat$mdod[ii] - 20) & tt <= (dat$mdod[ii] - 10), late.sc, 1) # male late-stage from 20-10 months before death
                  }else{
                    m.lt.sc <- 1
                  }
                if(dat$fser[ii]) {
                    f.lt.sc <- ifelse(tt > (dat$fdod[ii] - 20) & tt <= (dat$fdod[ii] - 10), late.sc, 1) # female late-stage from 20-10 months before death
                  }else{
                    f.lt.sc <- 1
                  }
                ###################################################################### 
                ## AIDS phase scalars
                if(dat$mser[ii]) {
                    m.aids.sc <- ifelse(tt > (dat$mdod[ii] - 10), aids.sc, 1) # male late-stage from 20-10 months before death
                  }else{
                    m.aids.sc <- 1
                  }
                if(dat$fser[ii]) {
                    f.aids.sc <- ifelse(tt > (dat$fdod[ii] - 10), aids.sc, 1) # female late-stage from 20-10 months before death
                  }else{
                    f.aids.sc <- 1
                  }
                ## Bernoulli within-couple transmission probabilities in current month; hets: p, gen
                temp.prob <- 1 - exp(-c(f.ac.sc*f.lt.sc*f.aids.sc*bmp*dat$m.het.hilo[ii]*dat$m.het.p[ii]*dat$m.het.gen[ii],
                                        m.ac.sc*m.lt.sc*m.aids.sc*bfp*dat$f.het.hilo[ii]*dat$f.het.p[ii]*dat$f.het.gen[ii]))
                from.part <- rbinom(2, 1, temp.prob)                   # Bernoulli random variables
                from.part <- from.part * c(dat$fser[ii], dat$mser[ii]) # only counts if partner is infected
                ## if new male infection from partner
                if(dat$mser[ii]==0 & from.part[1] == 1) {              
                    dat$mser[ii] <- 1
                    dat$mdoi[ii] <- tt
                    dat$mdod[ii] <- tt + ageweib(dat$mage[ii] - (dat$tint[ii] - tt), death = death)
                    dat$mcoi[ii] <- "p" # within-couple
                    ## track phase of infector
                    if((tt - dat$fdoi[ii])<=2) {
                        dat$mcoi.phase[ii] <- 'a' #acute
                      }else{
                        if(tt <= (dat$fdod[ii] - 20)) {
                            dat$mcoi.phase[ii] <- 'c' #chronic
                          }else{
                            if(tt <= (dat$fdod[ii] - 10)) {
                                dat$mcoi.phase[ii] <- 'l' #late
                              }else{
                            if(tt > (dat$fdod[ii] - 10) & tt <= (dat$fdod[ii])) dat$mcoi.phase[ii] <- 'ad' #aids
                              }
                          }
                      }
                  }
                ## if new female infection from partner
                if(dat$fser[ii]==0 & from.part[2] == 1) {
                    dat$fser[ii] <- 1
                    dat$fdoi[ii] <- tt
                    dat$fdod[ii] <- tt + ageweib(dat$fage[ii] - (dat$tint[ii] - tt), death = death)
                    dat$fcoi[ii] <- "p" # within-couple
                    if((tt - dat$mdoi[ii])<=2) {
                        dat$fcoi.phase[ii] <- 'a' #acute
                      }else{
                        if(tt <= (dat$mdod[ii] - 20)) {
                            dat$fcoi.phase[ii] <- 'c' #chronic
                          }else{
                            if(tt <= (dat$mdod[ii] - 10)) {
                                dat$fcoi.phase[ii] <- 'l' #late
                              }else{
                                if(tt > (dat$mdod[ii] - 10) & tt <= (dat$mdod[ii])) dat$fcoi.phase[ii] <- 'ad'
                              }
                          }
                      }
                  }
                if(sum(from.part)==0) { ## if there wasn't a within-couple infection (implying now concordant positive) (added Dec 6, 2013) to fix
                                        ## calculations of within-couple hazards for Rakai analysis. Previously we were slightly underestimating
                                        ## estimating hazards because some 'p' infections got replaced by 'e' infections.  extra-couple infections
                                        ## c(m,f); temporary transmission probability; hets: e, gen, beh
                  exc <- rbinom(2, 1, prob = c(1 - exp(-bme*dat$m.het.hilo[ii]*dat$m.het.e[ii]*dat$m.het.gen[ii]*dat$m.het.beh[ii]*epicf[tt, epic.ind.temp]),
                                               1 - exp(-bfe*dat$f.het.hilo[ii]*dat$f.het.e[ii]*dat$f.het.gen[ii]*dat$f.het.beh[ii]*epicm[tt, epic.ind.temp])))
                  ## if new m inf from extracouple
                  if(dat$mser[ii]==0 & exc[1] == 1) {
                    dat$mser[ii] <- 1
                    dat$mdoi[ii] <- tt
                    dat$mdod[ii] <- tt + ageweib(dat$mage[ii] - (dat$tint[ii] - tt), death = death)
                    dat$mcoi[ii] <- "e" # extra-couple
                  }
                  ## if new f inf from extracouple
                  if(dat$fser[ii]==0 & exc[2] == 1) {
                    dat$fser[ii] <- 1
                    dat$fdoi[ii] <- tt
                    dat$fdod[ii] <- tt + ageweib(dat$fage[ii] - (dat$tint[ii] - tt), death = death)
                    dat$fcoi[ii] <- "e" # extra-couple
                  }
                }
                end <- min(dat$mdod[ii], dat$fdod[ii], dat$tint[ii], na.rm = TRUE)
                tt <- tt+1
              } ## while still stuff to do
          }     ## end marriage if statement
      }         ## loop over couples
    return(dat)
  }

######################################################################
## Convert outut from event.fn into a time series (single core)
## ts: rows are months from 1900-2011, columns are couples, cell gives couple serostatus/death/aged out state
## ts: rows are months from 1900-2011, columns give # of couples in each serostatus/death/aged out state
## rakacRR: calcuates retrospective cohort style relative risk of acute phase infectiousness a la Rakai (see function below)
######################################################################
ts.fxn <- function(dat, nc = 12, verbose = T, vfreq = 500, return.ts = F, browse = F, do.rak = T, start.rak = 1994, end.rak = 1999)
    {
        if(browse)  browser()
        start <- min(dat[,c("tms","tfs")])
        end <- max(dat[,c("tint")])
        ## divide data into breaks for each core
        breaks <- rep(1:nc, length.out = nrow(dat))
        ## keep track of reordering for later
        uid.brk <- dat$uid[order(breaks)]
        ## use tsloop within each core
        multi.out <- mclapply(1:nc, tsloop, dat = dat, breaks = breaks, start = start, end = end, verbose = verbose, vfreq = vfreq)
        ## combine back into data frame
        tss <- multi.out[[1]][['tss']]              # each is a list(ts = ts, tss = tss) object
        ts <- multi.out[[1]][['ts']]                                               
        if(nc>1) {
            for(ii in 2:nc)     tss[,-(1:2)] <- tss[,-(1:2)] + multi.out[[ii]][['tss']][,-(1:2)]
            for(ii in 2:nc)     ts <- cbind(ts, multi.out[[ii]][['ts']])
            ts <- ts[,order(uid.brk)]
            if(do.rak) {
                hrout <- empir.arh(dat = dat, ts = ts, start.rak = start.rak, end.rak = end.rak)
                print(paste('Rakai style acute phase relative risk:', signif(hrout[[2]][2,1],5)))
            }else{ hrout <- NA }
        }
        gc() ## clear memory since this is memory intensive.
        if(!return.ts) { # ts is a months X couples matrix, and is extremeley big for large simulations, default to not return this.
            return(list(tss = tss, hrout = hrout))
        }else{
            return(list(ts = ts, tss = tss, hrout = hrout))
        }
    }

## Single core function to work on a batch of couples and turn it into a couple population time series
tsloop <- function(batch, dat, breaks, start, end, verbose, vfreq, browse = F)
  {
    if(browse)  browser()
    dat <- dat[breaks==batch,]          # batch to work on
    ## initialize time series array
    tseq <-start:end
    ts <- matrix(NA, nr = max(tseq), nc = nrow(dat))
    for(ii in 1:nrow(dat)) {            # for each couple
        if(verbose & ii%%vfreq==0) print(paste("on couple",ii,"of",nrow(dat))) # update on progress
        if(dat$ser[ii]==4)                #  if the couple is --
          {
            ts[min(dat[ii,c("tms","tfs")]):(dat$tmar[ii]-1),ii] <- "b.ss" # pre-couple concordant negative
            ts[dat$tmar[ii]:dat$tint[ii],ii] <- "ss"                      # couple concordant negative
          }
        ###################################################################### 
        if(dat$ser[ii]==3)              #  if the couple is F+
          {
            if(dat$fdoi[ii] < dat$tmar[ii]) # inf before mar
              {
                ts[min(dat[ii,c("tms","tfs")]):(dat$fdoi[ii]-1),ii] <- "b.ss" # pre-couple concordant negative
                ts[dat$fdoi[ii]:(dat$tmar[ii]-1),ii] <- "b.ff"                # pre-couple female positive
                if(dat$fdod[ii] < dat$tint[ii]) # if death before interview
                  {
                    ts[dat$tmar[ii]:(dat$fdod[ii]-1),ii] <- "ff" # female positive
                    ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.ff" # dead female positive
                  }else{                    # if death after interview
                    ts[dat$tmar[ii]:(dat$tint[ii]),ii] <- "ff" # female positive
                  }
              }else{                        # inf after mar
                ts[min(dat[ii,c("tms","tfs")]):(dat$tmar[ii]-1),ii] <- "b.ss" # pre-couple concordant negative
                if(dat$fdod[ii] < dat$tint[ii]) # if death before interview
                  {
                    ts[dat$tmar[ii]:(dat$fdoi[ii]-1),ii] <- "ss" # concordant negative
                    ts[dat$fdoi[ii]:(dat$fdod[ii]-1),ii] <- "ff" # female positive
                    ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.ff" # dead female positive
                  }else{                     # if death after interview
                    ts[dat$tmar[ii]:(dat$fdoi[ii]-1),ii] <- "ss" # concordant negative
                    ts[dat$fdoi[ii]:(dat$tint[ii]),ii] <- "ff"   # female positive
                  }
              }                             # end inf bef/after mar if/else
          }                                 # end F+
        ###################################################################### 
        if(dat$ser[ii]==2)              #  if the couple is M+ 
          {
            if(dat$mdoi[ii] < dat$tmar[ii]) # inf before mar
              {
                ts[min(dat[ii,c("tms","tfs")]):(dat$mdoi[ii]-1),ii] <- "b.ss"
                ts[dat$mdoi[ii]:(dat$tmar[ii]-1),ii] <- "b.mm"
                if(dat$mdod[ii] < dat$tint[ii]) # if death before interview
                  {
                    ts[dat$tmar[ii]:(dat$mdod[ii]-1),ii] <- "mm"
                    ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.mm" 
                  }else{                    # if death after interview
                    ts[dat$tmar[ii]:(dat$tint[ii]),ii] <- "mm" 
                  }
              }else{                        # inf after mar
                ts[min(dat[ii,c("tms","tfs")]):(dat$tmar[ii]-1),ii] <- "b.ss"
                if(dat$mdod[ii] < dat$tint[ii]) # if death before interview
                  {
                    ts[dat$tmar[ii]:(dat$mdoi[ii]-1),ii] <- "ss"
                    ts[dat$mdoi[ii]:(dat$mdod[ii]-1),ii] <- "mm"
                    ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.mm" 
                  }else{                     # if death after interview
                    ts[dat$tmar[ii]:(dat$mdoi[ii]-1),ii] <- "ss"
                    ts[dat$mdoi[ii]:(dat$tint[ii]),ii] <- "mm"
                  }
              }                             # end inf bef/after mar if/else
          }                                 # end M =    
        ###################################################################### 
        if(dat$ser[ii]==1)              #  if the couple is ++
          {                             # infected on the same date
            if(dat$mdoi[ii] == dat$fdoi[ii]) # mdoi=fdoi ##########
              {
                if(dat$mdoi[ii] < dat$tmar[ii]) # infections before marriage
                  {
                    ts[min(dat[ii,c("tms","tfs")]):(dat$mdoi[ii]-1),ii] <- "b.ss"
                    ts[dat$mdoi[ii]:(dat$tmar[ii]-1),ii] <- "b.hh"
                    if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                      {
                        if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                          {
                            ts[dat$tmar[ii]:(dat$mdod[ii]-1),ii] <- "hh" # concordant positive
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m" # male-dead-first dead concordant positive
                          }
                        if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                          {
                            ts[dat$tmar[ii]:(dat$fdod[ii]-1),ii] <- "hh" # concordant positive
                            ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f" # female-dead-first dead concordant positive
                          }
                        if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                          {
                            ts[dat$tmar[ii]:(dat$mdod[ii]-1),ii] <- "hh" # concordant positive
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1) # randomly choose who died first
                          }
                      }else{                # interview as a ++ couple
                        ts[dat$tmar[ii]:(dat$tint[ii]),ii] <- "hh"                 
                      }
                  ##### 
                  }else{                    # infections after marriage
                    ts[min(dat[ii,c("tms","tfs")]):(dat$tmar[ii]-1),ii] <- "b.ss"
                    ts[dat$tmar[ii]:(dat$mdoi[ii]-1),ii] <- "ss"
                    if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                      {
                        if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                          {
                            ts[dat$mdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh" 
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m"
                          }
                        if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                          {
                            ts[dat$mdoi[ii]:(dat$fdod[ii]-1),ii] <- "hh" 
                            ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f"
                          }
                        if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                          {
                            ts[dat$mdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh"                        
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1)
                          }
                      }else{                # interview as a ++ couple
                        ts[dat$mdoi[ii]:(dat$tint[ii]),ii] <- "hh"                 
                      }
                  }
              }
            ##############################
            ##  male infected first
            if(dat$mdoi[ii] < dat$fdoi[ii]) # M first ##########
              {
                if(dat$mdoi[ii] < dat$tmar[ii]) # if mdoi before mar
                  {
                    ts[min(dat[ii,c("tms","tfs")]):(dat$mdoi[ii]-1),ii] <- "b.ss"
                    if(dat$fdoi[ii] < dat$tmar[ii]) # and fdoi befor mar
                      {
                        ts[dat$mdoi[ii]:(dat$fdoi[ii]-1),ii] <- "b.mm"
                        ts[dat$fdoi[ii]:(dat$tmar[ii]-1),ii] <- "b.hh"
                        if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                          {
                            if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                              {
                                ts[dat$tmar[ii]:(dat$mdod[ii]-1),ii] <- "hh" 
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m"
                              }
                            if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                              {
                                ts[dat$tmar[ii]:(dat$fdod[ii]-1),ii] <- "hh" 
                                ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f"
                              }
                            if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                              {
                                ts[dat$tmar[ii]:(dat$mdod[ii]-1),ii] <- "hh"                        
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1)
                              }
                          }else{            # they interview as ++/alive
                            ts[dat$tmar[ii]:(dat$tint[ii]),ii] <- "hh"
                          }
                      }else{                # mdoi before mar but fdoi after
                        ts[dat$mdoi[ii]:(dat$tmar[ii]-1),ii] <- "b.mm"
                        ts[dat$tmar[ii]:(dat$fdoi[ii]-1),ii] <- "mm"
                        if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                          {
                            if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                              {
                                ts[dat$fdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh" 
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m"
                              }
                            if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                              {
                                ts[dat$fdoi[ii]:(dat$fdod[ii]-1),ii] <- "hh" 
                                ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f"
                              }
                            if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                              {
                                ts[dat$fdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh"                        
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1)
                              }
                          }else{            # they interview as ++/alive
                            ts[dat$fdoi[ii]:(dat$tint[ii]),ii] <- "hh"
                          }
                      }
                  }else{                        # both mdoi & fdoi after marriage
                    ts[min(dat[ii,c("tms","tfs")]):(dat$tmar[ii]-1),ii] <- "b.ss"
                    ts[dat$tmar[ii]:(dat$mdoi[ii]-1),ii] <- "ss"
                    ts[dat$mdoi[ii]:(dat$fdoi[ii]-1),ii] <- "mm"
                    if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                      {
                        if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                          {
                            ts[dat$fdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh" 
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m"
                          }
                        if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                          {
                            ts[dat$fdoi[ii]:(dat$fdod[ii]-1),ii] <- "hh" 
                            ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f"
                          }
                        if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                          {
                            ts[dat$fdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh"                        
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1)
                          }
                      }else{            # they interview as ++/alive
                        ts[dat$fdoi[ii]:(dat$tint[ii]),ii] <- "hh"
                      }
                  }            
              }
            ##################################################
            ##  female infected first
            if(dat$fdoi[ii] < dat$mdoi[ii]) # F first ##########
              {
                if(dat$fdoi[ii] < dat$tmar[ii]) # if mdoi before mar
                  {
                    ts[min(dat[ii,c("tms","tfs")]):(dat$fdoi[ii]-1),ii] <- "b.ss"
                    if(dat$mdoi[ii] < dat$tmar[ii]) # and mdoi befor mar
                      {
                        ts[dat$fdoi[ii]:(dat$mdoi[ii]-1),ii] <- "b.ff"
                        ts[dat$mdoi[ii]:(dat$tmar[ii]-1),ii] <- "b.hh"
                        if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                          {
                            if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                              {
                                ts[dat$tmar[ii]:(dat$mdod[ii]-1),ii] <- "hh" 
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m"
                              }
                            if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                              {
                                ts[dat$tmar[ii]:(dat$fdod[ii]-1),ii] <- "hh" 
                                ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f"
                              }
                            if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                              {
                                ts[dat$tmar[ii]:(dat$mdod[ii]-1),ii] <- "hh"                        
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1)
                              }
                          }else{            # they interview as ++/alive
                            ts[dat$tmar[ii]:(dat$tint[ii]),ii] <- "hh"
                          }
                      }else{                # fdoi before mar but mdoi after
                        ts[dat$fdoi[ii]:(dat$tmar[ii]-1),ii] <- "b.ff"
                        ts[dat$tmar[ii]:(dat$mdoi[ii]-1),ii] <- "ff"
                        if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                          {
                            if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                              {
                                ts[dat$mdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh" 
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m"
                              }
                            if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                              {
                                ts[dat$mdoi[ii]:(dat$fdod[ii]-1),ii] <- "hh" 
                                ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f"
                              }
                            if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                              {
                                ts[dat$mdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh"                        
                                ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1)
                              }
                          }else{            # they interview as ++/alive
                            ts[dat$mdoi[ii]:(dat$tint[ii]),ii] <- "hh"
                          }
                      }
                  }else{                        # both mdoi & fdoi after marriage
                    ts[min(dat[ii,c("tms","tfs")]):(dat$tmar[ii]-1),ii] <- "b.ss"
                    ts[dat$tmar[ii]:(dat$fdoi[ii]-1),ii] <- "ss"
                    ts[dat$fdoi[ii]:(dat$mdoi[ii]-1),ii] <- "ff"
                    if(min(dat[ii,c("mdod","fdod")]) < dat$tint[ii]) # if first death is before interview
                      {
                        if(dat$mdod[ii] < dat$fdod[ii]) # if male dies first
                          {
                            ts[dat$mdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh" 
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- "d.hh.m"
                          }
                        if(dat$fdod[ii] < dat$mdod[ii]) # if female dies first
                          {
                            ts[dat$mdoi[ii]:(dat$fdod[ii]-1),ii] <- "hh" 
                            ts[dat$fdod[ii]:(dat$tint[ii]),ii] <- "d.hh.f"
                          }
                        if(dat$mdod[ii] == dat$fdod[ii]) # if same date of death
                          {
                            ts[dat$mdoi[ii]:(dat$mdod[ii]-1),ii] <- "hh"                        
                            ts[dat$mdod[ii]:(dat$tint[ii]),ii] <- sample(c("d.hh.m","d.hh.f"),1)
                          }
                      }else{            # they interview as ++/alive
                        ts[dat$mdoi[ii]:(dat$tint[ii]),ii] <- "hh"
                      }
                  }            
              }
          }                                 # end ser==1
      }                                     # end couples for loop
    ## Now replace all live couple states with a.ss, a.mm, a.ff,
    ## a.hh, a.hh to reflect that they aged out of the cohort and
    ## should no longer be considered in SDP or other (i.e.,
    ## prevalence) calculations
    ## add postfix .m or .f to reflect who aged out
    for(cc in 1:nrow(dat))
      {
        if(!is.na(dat$dage[cc]))        # if they age out
          {
            if(dat$dage[cc] < dat$tend[cc]) # if they aged out before interview time or AIDS related death
              {
                ## then replace their couple state after they aged out with the state before aging out preceded by 'a.'
                ts[1:dat$tint > dat$dage[cc],cc] <- paste('a.',ts[1:dat$tint == dat$dage[cc],cc], sep = '')
                if(dat$dmage[cc]<=dat$dfage[cc]) # if male is first
                  {
                    ts[1:dat$tint > dat$dage[cc],cc] <- paste(ts[1:dat$tint > dat$dage[cc],cc], '.m', sep = '') # add m
                  }else{                # otherwise
                    ts[1:dat$tint > dat$dage[cc],cc] <- paste(ts[1:dat$tint > dat$dage[cc],cc], '.f', sep = '') # add f
                  }
              }
          }
      }
    ##  set names of columns
    nms <- c("b.ss","b.mm","b.ff","b.hh","ss","mm","ff","hh","d.mm","d.ff","d.hh.m","d.hh.f","a.ss.m","a.ss.f","a.mm.m","a.mm.f","a.ff.m","a.ff.f","a.hh.m","a.hh.f")
    tss <- as.data.frame(matrix(NA, nr = length(tseq), nc = length(nms) + 2))
    names(tss) <- c("tt","yr",nms)
    tss$tt <- tseq                      # months since 1900
    tss$yr <- tseq/12+1900              # year
    for(tt in 1:length(tseq)) {
        for(nn in nms) {
            tss[tt,nn] <- sum(ts[tseq[tt],]==nn, na.rm=T)
          }
      }
    tss$alive <- rowSums(tss[,c("ss","mm","ff","hh")]) # number formed couples alive at any given time
    tss$inf.alive <- rowSums(tss[,c("mm","ff","hh")])  # number of infected, formed live couples at any given time
    tss$sdp <- (tss[,'mm']+tss[,'ff'])/tss[,'inf.alive'] # serodiscordant proportion (SDP)
    tss$inf <- rowSums(tss[,c("mm","ff","hh","d.mm","d.ff","d.hh.m","d.hh.f")]) # number of infected couples
    ## cumulative infections by route
    tss$cu.mb <- NA
    tss$cu.me <- NA
    tss$cu.mp <- NA
    tss$cu.fb <- NA
    tss$cu.fe <- NA
    tss$cu.fp <- NA
    for(tt in tss[,'tt'])               # cumulative infections by route
      {
        tss[tss[,'tt']==tt, 'cu.mb'] <- sum(dat[dat$mser==1 & dat$mdoi<= tt, 'mcoi']=='b')
        tss[tss[,'tt']==tt, 'cu.me'] <- sum(dat[dat$mser==1 & dat$mdoi<= tt, 'mcoi']=='e')
        tss[tss[,'tt']==tt, 'cu.mp'] <- sum(dat[dat$mser==1 & dat$mdoi<= tt, 'mcoi']=='p')    
        tss[tss[,'tt']==tt, 'cu.fb'] <- sum(dat[dat$fser==1 & dat$fdoi<= tt, 'fcoi']=='b')
        tss[tss[,'tt']==tt, 'cu.fe'] <- sum(dat[dat$fser==1 & dat$fdoi<= tt, 'fcoi']=='e')
        tss[tss[,'tt']==tt, 'cu.fp'] <- sum(dat[dat$fser==1 & dat$fdoi<= tt, 'fcoi']=='p')
      }
    ## change ts to include acute/chronic/late/aids phases
    return(list(ts = ts, tss = tss))
  }

## Find order function. Finds the ind-th value of all the rows a matrix x.
ord <- function(x,ind) {
  unlist(apply(x, 1, function(x) x[which(x == x[order(x)][ind])][1]))
}



## Create a pseudo-population and do it event driven simulation with that then compiling it into a
## population timeseries.
psrun <- function(country, s.demog = NA, # country to simulate;  country whose relationship patterns are to be used in simulation
                  seed = 1,      # seed to set for random number generation
                  ## population size
                  maxN = 10^4, # maximum pseudo-population size (sample inflated pop) to avoid memory problems
                  ## non-parametric approach to generate couple pseudo-population
                  infl.fac = 4,      #  pseudo-population inflation factor
                  last.int = F,      #  set interview dates to be the latest one for all individuals
                  psNonPar = F,      # use non-parametric pseudo-population generator?
                  pars,              # transmission parameters to use
                  sample.tmar = T,   # sample marriage time too?
                  ## for parametric approach to generate couple pseudo-population
                  tmar, each, tint,    # marriage times, # at each marriage time & interview time
                  start.rak, end.rak, # years of retrospective cohort analysis for Rakai style acute phase study
                  return.ts = F,
                  ## sensitivity analyses
                  death = T,             #  include AIDS mortality?
                  acute.sc, late.sc, aids.sc, #  relative hazards of HIV phases compared to chronic phase
                  ##  next three lines scale transmission coefficients by the following values (used
                  ## for counterfactual simulations where different routes are amplified or
                  ## diminished)
                  bmb.sc = 1, bfb.sc = 1, #  pre-couple 
                  bme.sc = 1, bfe.sc = 1, #  extra-couple                 
                  bmp.sc = 1, bfp.sc = 1, #  within-couple
                  counterf.betas = F, # change betas in counterfactuals? if not change beta_within & c's (so beta_within affects all routes)
                  ##  Route-specific heterogeneity parameters {do this type of heterogeneity?;
                  ##  standard deviation of lognormal risk deviate; inter-partner correlation}
                  het.b = F, het.b.sd = 0, het.b.cor = 0, #  pre-couple   
                  het.e = F, het.e.sd = 0, het.e.cor = 0, #  extra-couple 
                  het.p = F, het.p.sd = 0, het.p.cor = 0, #  within-couple
                  het.gen = F, het.gen.sd = 0, het.gen.cor = 0, # same for genetic heterogeneity ( same lognormal risk deviate for all three routes)
                  het.beh = F, het.beh.sd = 0, het.beh.cor = 0, # same for behavioral heterogeneity ( same lognormal risk deviate for pre- & extra-couple routes)
                  ## individual suscepbitility constant through life is either HIGH or LOW risk,
                  ## proportion high risk (phigh) are at rrhigh times higher risk
                  hilo = F, phigh.m = .2, phigh.f = .2, rrhigh.m = 10, rrhigh.f = 10,
                  scale.by.sd = T, # adjust beta means to keep geometric mean constant with increasing heterogeneity
                  scale.adj = 1,   # adjust betas arbitrarily
                  one.couple = F, # for debugging, repeat one couple married in 1990, 100,000 times
                  ## processing & visualizations
                  nc = 12,               # number of cores
                  out.dir,               # output directory
                  make.jpgs = T,         # make some pictures of results
                  early.yr = 1985,       # earliest year to show in timeseries plots
                  plot.pdfs = T,         # make some PDFs of results
                  save.new = F,          # always save as a new file (don't overwrite)
                  vfreq = 100,           # progress update frequency
                  mdcol = "dark green",  # male discordant color
                  fdcol = "purple",      # female discordant color
                  sscol = "dark gray",   # concordant negative color
                  cccol = "red",         # concordant positive color
                  browse = F)            # debug
  {
    set.seed(seed)                      # set RNG seed
    start.time1 <- Sys.time()           # track computation time
    odat <- dat                         # store original data set workspace
    if(browse) browser()                # debug
    ## calculate pre-couple duration (earliest sexual debut to couple formation) 
    odat$bd <- apply(cbind(odat$tmar-odat$tms,odat$tmar-odat$tfs), 1, max) 
    ##  calculate couple duration
    odat$cd <- odat$tint - odat$tmar
    dat <- odat[odat$group==ds.nm[s.demog],]
    ##  male and female durations of sex before marriage, just to see how
    ##  they affect output later on.
    dat$mds <- dat$tmar - dat$tms
    dat$fds <- dat$tmar - dat$tfs
######################################## 
    ## Generate Pseudo-Population
########################################
    if(psNonPar) { ## NON-Parametric APPROACH:Inflate observed couples times their inverse
                   ## probability of having survived to be observed. First, simulate with state
                   ## probability model to get inflation factors for each couple based on
                   ## probability of survival
        sim <- pcalc(pars, dat, sim = T, survive = T, browse = F, trace = F) # Use Markov State Probability model to get survival probabilities
        surp <- rowSums(sim$pser.a) # survival probabilities =  probability couples are in one of the alive states at DHS interview
        infl <- 1/surp              # infl fator = 1 / survival probabilities
        ## Generate pseudocouple data set based on infl.fac X # of couples per infl couple (infl.fac
        ## just allows us to create bigger pseudopopulations).
        numb <- rpois(rep(1,nrow(dat)),infl.fac*infl) # poisson RNG for # of pseudocouples per real couple
        psdat <- dat[rep(1:nrow(dat), numb),]
        ## Reduce pseudopopulation to maximum psuedopopulation size desired
        if(nrow(psdat) > maxN) {
            smp <- sample(1:nrow(psdat), maxN)
            smp <- smp[order(smp)]
            psdat <- psdat[smp,]
          }
        ## make interview dates all the last interview date so we don't truncate some couple's simulations
        psdat$tint <- max(psdat$tint)
      }else{ # Parametric Approach: use normal copulas fitted to relationship history patterns for each country to simulate pseudo-population.
        psdat <- rcop(s.demog, NN = maxN, sample.tmar = sample.tmar, tmar = tmar, each = each, tint = tint, one.couple = one.couple, browse = F)
      }
    ## Plots pseudo-population distributions for exploration later.
    psdat.add <- add.vars(psdat) # Append extra variables (see psdat.add() below) to pseudo-population
    real <- add.vars(dat)       # Append extra variables (see psdat.add() below) to real population
    rgs <- apply(rbind(real[,vars],psdat.add[,vars]), 2 , range) # get full range for all variables
    sz <- 800                                                    # JPG dimensions
    if(make.jpgs) {                                              ##  show copula fits
        ## Secular trends of simulated data
        jpeg(file.path(out.dir,paste(ds.nm[country], "sec trend SIM.jpg", sep = "")), width = sz, height = sz)
        par(mfrow=c(2,2))
        for(vv in 1:4) {                # for each variable
            scatter.smooth(psdat.add$tmar, psdat.add[,vars[vv]], xlab = "couple formation date",
                           ylab = labs[vv], bty = "n", pch = 19, cex = .5, ylim = rgs[,vv])
            abline(lm(psdat.add[,vars[vv]]~psdat.add$tmar), col = "red", lwd = 3)
            mtext(paste(ds.nm[country],"nonpar"[psNonPar],"smpTmar"[sample.tmar], "synthCoh"[!psNonPar & !sample.tmar]),
                  side = 3, line = 0, outer = T, cex = 2)        
          }
        legend("topleft", c("linear fit", "loess fit"), col = c("red","black"), lty = 1, lwd = 1, bty = "n")    
        dev.off()
        ## Secular trends of real data for comparison
        jpeg(file.path(out.dir,paste(ds.nm[country], "sec trend Real.jpg", sep = "")), width = sz, height = sz)
        par(mfrow=c(2,2))
        for(vv in 1:4) {                # for each variable
            scatter.smooth(real$tmar, real[,vars[vv]], xlab = "couple formation date", ylab = labs[vv], bty = "n", pch = 19, cex = .5, ylim = rgs[,vv])
            abline(lm(real[,vars[vv]]~real$tmar), col = "red", lwd = 3)
            mtext(paste(ds.nm[country],"real data"), side = 3, line = 0, outer = T, cex = 2)        
          }
        legend("topleft", c("linear fit", "loess fit"), col = c("red","black"), lty = 1, lwd = 1, bty = "n")        
        dev.off()
        ## 5D multivariate distribution
        sz <- 1200
        ## Simulated data
        jpeg(file.path(out.dir,paste(ds.nm[country], "N-", maxN,"-distr-SIM.jpg", sep = "")), w = sz, h = sz)
        ## lims & breaks chosen from multi-countyr comparison code (pscor.R)
        ylim <- data.frame(l = rep(0,5), u = c(.02, .04, .25, .3, .01))
        breaks <- list(seq(0,60*12, by = 12), seq(0,60*12, by = 12), seq(0,60*12, by = 1), seq(0,60*12, by = 1), # for histograms below
                       seq(600, 1600, by = 12))
        ##  pairwise scatterplots of simulated data
        sbpairs(psdat.add[,vars], browse = F, do.pdf = F, do.jpeg = F, rgs = rgs, hist.freq = F, yrg = ylim, breaks = breaks)
        mtext(paste(ds.nm[country],"nonpar"[psNonPar],"smpTmar"[sample.tmar], "synthCoh"[!psNonPar & !sample.tmar]), side = 3, line = 0, outer = T, cex = 2)
        dev.off()
        ## pairwise scatterplots of real data
        jpeg(file.path(out.dir,paste(ds.nm[country], "N-", nrow(dat),"-distr-REAL.jpg", sep = "")), w = sz, h = sz)
        sbpairs(real[,vars], browse = F, do.pdf = F, do.jpeg = F, rgs = rgs, hist.freq = F, yrg = ylim, breaks = breaks)
        mtext(paste(ds.nm[country],"actual data"), side = 3, line = 0, outer = T, cex = 2)
        dev.off()
      }
    ######################################################################
    ## Simulate couple transmission model with pseudo-population
    ######################################################################
    ## Scale extra-couple parameters
    print('Running simulation')
    cpars <- pars
    if(counterf.betas) { ## change betas in counterfactuals
    cpars["bmb"] <- cpars["bmb"]*bmb.sc
    cpars["bfb"] <- cpars["bfb"]*bfb.sc   
    cpars["bme"] <- cpars["bme"]*bme.sc
    cpars["bfe"] <- cpars["bfe"]*bfe.sc   
    cpars["bmp"] <- cpars["bmp"]*bmp.sc
    cpars["bfp"] <- cpars["bfp"]*bfp.sc
  }else{ ## change HIV transmission rate (beta_within) & contact coefficients (c's) in counterfactuals
    cpars["bmb"] <- cpars["bmb"]*bmb.sc*bmp.sc
    cpars["bfb"] <- cpars["bfb"]*bfb.sc*bfp.sc
    cpars["bme"] <- cpars["bme"]*bme.sc*bmp.sc
    cpars["bfe"] <- cpars["bfe"]*bfe.sc*bfp.sc
    cpars["bmp"] <- cpars["bmp"]*bmp.sc
    cpars["bfp"] <- cpars["bfp"]*bfp.sc
  }
    ##  call event-driven simulation
    evout <- event.fn(cpars, psdat, vfreq = 100, death = death, acute.sc = acute.sc, late.sc = late.sc, aids.sc = aids.sc, nc = nc,
                      het.b = het.b, het.b.sd = het.b.sd, het.b.cor = het.b.cor,
                      het.e = het.e, het.e.sd = het.e.sd, het.e.cor = het.e.cor,
                      het.p = het.p, het.p.sd = het.p.sd, het.p.cor = het.p.cor, 
                      het.gen = het.gen, het.gen.sd = het.gen.sd, het.gen.cor = het.gen.cor,
                      het.beh = het.beh, het.beh.sd = het.beh.sd, het.beh.cor = het.beh.cor,
                      hilo = hilo, phigh.m = phigh.m, phigh.f = phigh.f, rrhigh.m = rrhigh.m, rrhigh.f = rrhigh.f,
                      scale.by.sd = scale.by.sd, scale.adj = scale.adj,
                      browse=F)
    ts <- NA # so there's something to return if not keeping it in next line
    temptsout <- ts.fxn(evout, return.ts = return.ts, start.rak = start.rak, end.rak = end.rak)     ##  convert line lists to population timeseries
    tss <- temptsout[['tss']]      # grab tss
    if(return.ts)       ts <- temptsout[['ts']]
    hrout <- temptsout[['hrout']]   # grab HIV phase object
    hours <- round(as.numeric(difftime(Sys.time(), start.time1, unit = "hour")),3) # runtime
    if(!exists(as.character(substitute(jobnum)))) jobnum <- NA # jobnum is set globally when on cluster
    if(!exists(as.character(substitute(simj)))) simj <- NA # simj is set globally when on cluster (job # within country/acute batch)
    ##  produce output list
    output <- list(jobnum = jobnum, simj = simj, evout = evout, tss = tss, ts = ts, hrout = hrout,
                   pars = c(bmb.sc = bmb.sc, bfb.sc = bfb.sc, bme.sc = bme.sc, bfe.sc = bfe.sc, bmp.sc = bmp.sc, bfp.sc = bfp.sc, 
                     death = death, acute.sc = acute.sc, late.sc = late.sc, aids.sc = aids.sc,
                     s.epic = s.epic, s.demog = s.demog,
                     s.bmb = s.bmb, s.bfb = s.bfb,
                     s.bme = s.bme, s.bfe = s.bfe,
                     s.bmp = s.bmp, s.bfp = s.bfp,                     
                     het.b = het.b, het.b.sd = het.b.sd, het.b.cor = het.b.cor,
                     het.e = het.e, het.e.sd = het.e.sd, het.e.cor = het.e.cor,
                     het.p = het.p, het.p.sd = het.p.sd, het.p.cor = het.p.cor, 
                     het.gen = het.gen, het.gen.sd = het.gen.sd, het.gen.cor = het.gen.cor,
                     het.beh = het.beh, het.beh.sd = het.beh.sd, het.beh.cor = het.beh.cor,
                     hilo = hilo, phigh.m = phigh.m, phigh.f = phigh.f, rrhigh.m = rrhigh.m, rrhigh.f = rrhigh.f,
                     group.ind = group.ind, infl.fac = infl.fac, maxN = maxN,
                     psNonPar = psNonPar, last.int = F, sample.tmar = sample.tmar, tint = tint,
                     hours = hours),
                   tmar = tmar, each = each) # these are vectors & so must be given separately
    sim.nm <- paste(ds.nm[country], "-", nrow(evout), '-', jobnum, sep = "") # name simulation results (country-#couples-jobum)
    ## create file name making sure not to save over old files
    output.nm.base <- file.path(out.dir, sim.nm) #paste(sim.nm, ".Rdata", sep = ""))
    stepper <- 1
    output.nm <- output.nm.base
    if(save.new) { ## don't save over old files
        while(file.exists(paste0(output.nm,'.Rdata'))) { 
            stepper <- stepper + 1
            output.nm <- paste0(output.nm.base, '-', stepper)
        }
    }
    output.nm <- paste0(output.nm, '.Rdata')
    print(paste('saving file',output.nm))
    save(output, file = output.nm)
    return(output.nm)
  }
## Returns results in a couples line list format, a time series format, and also gives all information on simulation parameters


## Create a function that takes the output from event.fn and makes it look like the sample
## proportions matching the output of pcalc$allstates to make sure that both event-driven and Markov
## State probability simulations give same input for same transmission coefficients
transf <- function(stuff)
  {
    with(stuff,
         {
           K <- nrow(stuff)
           s.. <- sum(mser+fser==0) / K
           mb. <- sum(mcoi == "b" & fser==0) / K
           me. <- sum(mcoi == "e" & fser==0) / K
           f.b <- sum(fcoi == "b" & mser==0) / K
           f.e <- sum(fcoi == "e" & mser==0) / K
           hb1b2 <- sum(mcoi == "b" & fcoi == "b" & (mdoi < fdoi)) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi)) / K
           hb2b1 <- sum(mcoi == "b" & fcoi == "b" & (mdoi > fdoi)) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi)) / K
           hbe <- sum(mcoi == "b" & fcoi == "e") / K
           heb <- sum(mcoi == "e" & fcoi == "b") / K
           hbp <- sum(mcoi == "b" & fcoi == "p") / K
           hpb <- sum(mcoi == "p" & fcoi == "b") / K
           hep <- sum(mcoi == "e" & fcoi == "p") / K
           hpe <- sum(mcoi == "p" & fcoi == "e") / K
           he1e2 <- sum(mcoi == "e" & fcoi == "e" & (mdoi < fdoi)) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi)) / K
           he2e1 <- sum(mcoi == "e" & fcoi == "e" & (mdoi > fdoi)) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi)) / K
           ## And alive
           mb.A <- sum(mcoi == "b" & fser==0 & alive) / K
           me.A <- sum(mcoi == "e" & fser==0 & alive) / K
           f.bA <- sum(fcoi == "b" & mser==0 & alive) / K
           f.eA <- sum(fcoi == "e" & mser==0 & alive) / K
           hb1b2A <- sum(mcoi == "b" & fcoi == "b" & (mdoi < fdoi) & alive) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi) & alive) / K
           hb2b1A <- sum(mcoi == "b" & fcoi == "b" & (mdoi > fdoi) & alive) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi) & alive) / K
           hbeA <- sum(mcoi == "b" & fcoi == "e" & alive) / K
           hebA <- sum(mcoi == "e" & fcoi == "b" & alive) / K
           hbpA <- sum(mcoi == "b" & fcoi == "p" & alive) / K
           hpbA <- sum(mcoi == "p" & fcoi == "b" & alive) / K
           hepA <- sum(mcoi == "e" & fcoi == "p" & alive) / K
           hpeA <- sum(mcoi == "p" & fcoi == "e" & alive) / K
           he1e2A <- sum(mcoi == "e" & fcoi == "e" & (mdoi < fdoi) & alive) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi) & alive) / K
           he2e1A <- sum(mcoi == "e" & fcoi == "e" & (mdoi > fdoi) & alive) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi) & alive) / K
           out <- c(s.., mb., me., f.b, f.e, hb1b2, hb2b1, hbe, heb, hbp, hpb, hep, hpe, he1e2, he2e1, mb.A,
                    me.A,  f.bA,  f.eA,  hb1b2A, hb2b1A, hbeA,  hebA,  hbpA, hpbA,  hepA,  hpeA,  he1e2A, he2e1A)
           return(out)
         })
  }


## Markov State probability transmission model formulation: each couple has a probability of being
## in each serostatus from the early sexual debut to the time of interview.  This is the model that
## is used to estimate transmission coefficients from DHS data.
pcalc <- function(pars, dat, browse = F, compars = NULL, # compare lprob for true pars with other pars (fitted)
                  give.pis = F,              # return individual pi values for couples (for outside mcmc)
                  sim = F,                   # just use this to simulate data given parameters? (outputs pser.a)
                  survive = T,               # account for survival in analysis
                  cond.sim = F,              # only simulate individuals that will live
                  trace = T) # only do certain calculations when tracing parameters (i.e. for non-thinned versions)
  {
    if(browse) browser()                # debug
    K <- nrow(dat)                      # number of couples
    if(!sim) {                          # using real data? Then assign logical indices so couples can be easily selected by serostatus
        hh.log <- dat$ser==1            #M+F+
        mm.log <- dat$ser==2            #M+F-
        ff.log <- dat$ser==3            #M-F+
        ss.log <- dat$ser==4            #M-F-
      }
    if(sum(pars[1:5]<0)>0)                   #if any parameters are <0 then the model must be rejected so we return logprob =-Inf
      {                                      # i.e., this is a boundary condition
        probs <- NA                     
        lprob <- -Inf
        pop.avs <- NA
        proj12 <- NA
        pser.a <- NA
        pser <- NA
        rrs <- NA
      }else{                            # attached transmission coefficients (prevalence-standardized hazards)
        bmb <- pars[["bmb"]]            #pre-couple to m
        bfb <- pars[["bfb"]]            #pre-couple to f
        bme <- pars[["bme"]]            #extra-couple to m
        bfe <- pars[["bfe"]]            #extra-couple to f
        bmp <- pars[["bmp"]]            #within-couple to m
        if("lrho" %in% names(pars))     # bfp/bmp = rho; i.e.  ratio of male to female to female to male transmission within a couple's
          {
            rho <- exp(pars[["lrho"]]) # feeding in log(rho)
            bfp <- bmp * rho           #within-couple to f
          }else{
            bfp <- pars[["bfp"]]
          }
        # L stands for *L*ast month iteration (always need to store state probabilities at time t-1 when calculating those that time t)
        s..L <- rep(1,K)                # s: concordant negative
        mb.L <- rep(0, K)               # m: male positive, pre-couply infected
        me.L <- rep(0, K)               # m: male positive, extra-couply infected
        f.bL <- rep(0, K)               # f: female positive, pre-couply infected  
        f.eL <- rep(0, K)               # f: female positive, extra-couply infected
        hb1b2L <- rep(0, K)      # h: concordant positive, both infected pre-couply but male first
        hb2b1L <- rep(0, K)      # h: concordant positive, both infected pre-couply but female first
        hbeL <- rep(0, K) # h: concordant positive,  male infected pre-couply,  female infected extra-couply
        hebL <- rep(0, K) # etc..., 2nd & 3rd character give route of transmission for M & F, respectively
        hbpL <- rep(0, K) 
        hpbL <- rep(0, K)
        hepL <- rep(0, K)
        hpeL <- rep(0, K)
        he2e1L <- rep(0, K)
        he1e2L <- rep(0, K)
        ## initiate vectors to update based on *L*ast state (i.e. time t)
        s.. <- rep(1,K)                 
        mb. <- rep(0, K)               
        me. <- rep(0, K)        
        f.b <- rep(0, K)               
        f.e <- rep(0, K)        
        hb1b2 <- rep(0, K)
        hb2b1 <- rep(0, K)                       
        hbe <- rep(0, K)
        heb <- rep(0, K)               
        hbp <- rep(0, K) 
        hpb <- rep(0, K)
        hep <- rep(0, K)
        hpe <- rep(0, K)
        he1e2 <- rep(0, K)
        he2e1 <- rep(0, K)        
        # i.e., hbp is a ++ couple in which the male was inf *b*efore
        # couple formation & the female by her *p*artner
        if(survive)
          {
            # A stands for *A*live, i.e. joint probability of serostatus
            # and both partners being alive at sampling
            mb.AL <- rep(0, K)
            me.AL <- rep(0, K)        
            f.bAL <- rep(0, K)
            f.eAL <- rep(0, K)        
            hb1b2AL <- rep(0, K)
            hb2b1AL <- rep(0, K)
            hbeAL <- rep(0, K)
            hebAL <- rep(0, K)
            hbpAL <- rep(0, K)
            hpbAL <- rep(0, K)
            hepAL <- rep(0, K)
            hpeAL <- rep(0, K)
            he1e2AL <- rep(0, K)
            he2e1AL <- rep(0, K)
            ## initiate vectors to update based on *L*ast state
            mb.A <- rep(0, K)
            me.A <- rep(0, K)        
            f.bA <- rep(0, K)
            f.eA <- rep(0, K)        
            hb1b2A <- rep(0, K)
            hb2b1A <- rep(0, K)
            hbeA <- rep(0, K)
            hebA <- rep(0, K)
            hbpA <- rep(0, K)
            hpbA <- rep(0, K)
            hepA <- rep(0, K)
            hpeA <- rep(0, K)
            he1e2A <- rep(0, K)
            he2e1A <- rep(0, K)                    
          }
        for(tt in 1:max(dat$bd))
          {
            ## probabilities are non-zero only for times after started having sex and before couple formation
            m.sex <- dat$tmar-dat$bd+tt-1 >= dat$tms & dat$tmar-dat$bd+tt-1 < dat$tmar
            f.sex <- dat$tmar-dat$bd+tt-1 >= dat$tfs & dat$tmar-dat$bd+tt-1 < dat$tmar
            e.sex <- m.sex|f.sex           # either are active
            ## probability infected in month tt
            p.m.bef <- rep(0,K)
            p.f.bef <- rep(0,K)    
            p.m.bef[m.sex] <- (1 - exp(-bmb * epicf[cbind(dat$tmar[m.sex]-dat$bd[m.sex]+tt-1, dat$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * epicm[cbind(dat$tmar[f.sex]-dat$bd[f.sex]+tt-1, dat$epic.ind[f.sex])]))
            ## probability infected in month tt and alive at sampling
            if(survive)
              {
                p.m.bef.a <- rep(0,K)
                p.f.bef.a <- rep(0,K)
                ## csurv[time til interview, age in months in this month]
                p.m.bef.a[m.sex] <- p.m.bef[m.sex] * csurv[cbind(dat$mage[m.sex]-dat$cd[m.sex]-dat$bd[m.sex]+tt-1, dat$cd[m.sex]+dat$bd[m.sex]-tt+1)]
                p.f.bef.a[f.sex] <- p.f.bef[f.sex] * csurv[cbind(dat$fage[f.sex]-dat$cd[f.sex]-dat$bd[f.sex]+tt-1, dat$cd[f.sex]+dat$bd[f.sex]-tt+1)]
              }
            ## iterate probabilities based on previous values for only cases where it needs updating
            s..[e.sex] <- s..L[e.sex]*(1-p.m.bef[e.sex])*(1-p.f.bef[e.sex])
            mb.[e.sex] <- mb.L[e.sex]*(1 - p.f.bef[e.sex]) + s..L[e.sex]*p.m.bef[e.sex]*(1-p.f.bef[e.sex])
            f.b[e.sex] <- f.bL[e.sex]*(1 - p.m.bef[e.sex]) + s..L[e.sex]*p.f.bef[e.sex]*(1-p.m.bef[e.sex])
            ## for individuals infected in the same month, assign
            ## the order of infection based on competing risks
            ## formula, but if the denominator is 0, replace both
            ## with 0 to avoid errors.
            p.mfirst <- p.m.bef[e.sex] / (p.m.bef[e.sex]+p.f.bef[e.sex])
            p.ffirst <- 1-p.mfirst
            p.mfirst[is.na(p.mfirst)] <- 0
            p.ffirst[is.na(p.ffirst)] <- 0                
            hb1b2[e.sex] <- hb1b2L[e.sex] + p.mfirst * s..L[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
                                            mb.L[e.sex]*p.f.bef[e.sex]
            hb2b1[e.sex] <- hb2b1L[e.sex] + p.ffirst * s..L[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
                                            f.bL[e.sex]*p.m.bef[e.sex]
            ## iterate joint probabilities with survival
            if(survive)
              {
                mb.A[e.sex] <- mb.AL[e.sex]*(1 - p.f.bef[e.sex]) + s..L[e.sex]*p.m.bef.a[e.sex]*(1-p.f.bef[e.sex])
                f.bA[e.sex] <- f.bAL[e.sex]*(1 - p.m.bef[e.sex]) + s..L[e.sex]*p.f.bef.a[e.sex]*(1-p.m.bef[e.sex])
                ## for individuals infected in the same month, assign
                ## the order of infection based on competing risks
                ## formula, but if the denominator is 0, replace both
                ## with 0 to avoid errors.
                p.mfirst.a <- p.m.bef.a[e.sex] / (p.m.bef.a[e.sex]+p.f.bef.a[e.sex])
                p.ffirst.a <- 1-p.mfirst.a
                p.mfirst.a[is.na(p.mfirst.a)] <- 0
                p.ffirst.a[is.na(p.ffirst.a)] <- 0                
                hb1b2A[e.sex] <- hb1b2AL[e.sex] + p.mfirst.a * s..L[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
                                                  mb.AL[e.sex] * p.f.bef.a[e.sex]
                hb2b1A[e.sex] <- hb2b1AL[e.sex] + p.ffirst.a * s..L[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
                                                  f.bAL[e.sex] * p.m.bef.a[e.sex]
                ## Update *A*live *L*ast states
                mb.AL[e.sex] <- mb.A[e.sex]
                f.bAL[e.sex] <- f.bA[e.sex]
                hb1b2AL[e.sex] <- hb1b2A[e.sex]
                hb2b1AL[e.sex] <- hb2b1A[e.sex]                
              }
            ## Update other *L*ast states
            s..L[e.sex] <- s..[e.sex]
            mb.L[e.sex] <- mb.[e.sex]
            f.bL[e.sex] <- f.b[e.sex]
            hb1b2L[e.sex] <- hb1b2[e.sex]
            hb2b1L[e.sex] <- hb2b1[e.sex]                        
          }
        ## probability of being infected by partner (constant, used inside loop)
        p.m.part <- 1 - exp(-bmp)
        p.f.part <- 1 - exp(-bfp)
        ## Now loop through marriage
        for(tt in 1:max(dat$cd-1))
          {
            ## are partners formed in a couple?
            fmd <- dat$cd >= tt              
            ######################################################################
            ## everything below is automatically sum(fmd) length except p.m/f.part which are length 1
            ## probability infected extracouply in the ttc-th month of couple
            p.m.exc <- (1 - exp(-bme*epicf[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe*epicm[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            if(survive)
              {
                ## Survival probabilities
                s.p.m <- csurv[cbind(dat$mage[fmd]-dat$cd[fmd]+tt-1, dat$tint[fmd] - (dat$tmar[fmd] + tt - 1))]
                s.p.f <- csurv[cbind(dat$fage[fmd]-dat$cd[fmd]+tt-1, dat$tint[fmd] - (dat$tmar[fmd] + tt - 1))]
                ## Transmission probabilities from partner (jointly with survival)
                p.m.part.a <- p.m.part * s.p.m
                p.f.part.a <- p.f.part * s.p.f
                p.m.exc.a <- p.m.exc * s.p.m
                p.f.exc.a <- p.f.exc * s.p.f
              }
            ######################################################################
            ## iterate probabilities
            s..[fmd] <- s..L[fmd]*(1-p.m.exc)*(1-p.f.exc)
            mb.[fmd] <- mb.L[fmd]*(1-p.f.exc)*(1-p.f.part)
            me.[fmd] <- me.L[fmd]*(1-p.f.exc)*(1-p.f.part) + s..L[fmd]*p.m.exc*(1-p.f.exc)
            f.b[fmd] <- f.bL[fmd]*(1-p.m.exc)*(1-p.m.part)
            f.e[fmd] <- f.eL[fmd]*(1-p.m.exc)*(1-p.m.part) + s..L[fmd]*p.f.exc*(1-p.m.exc)
##          hb1b2[fmd] <- hb1b2L[fmd] # Doesn't change during couple duration
##          hb2b1[fmd] <- hb2b1L[fmd] # Doesn't change during couple duration                
            hbe[fmd] <- hbeL[fmd] + mb.L[fmd]*(1-p.f.part)*p.f.exc
            heb[fmd] <- hebL[fmd] + f.bL[fmd]*(1-p.m.part)*p.m.exc
            hbp[fmd] <- hbpL[fmd] + mb.L[fmd]*p.f.part
            hpb[fmd] <- hpbL[fmd] + f.bL[fmd]*p.m.part
            hep[fmd] <- hepL[fmd] + me.L[fmd]*p.f.part
            hpe[fmd] <- hpeL[fmd] + f.eL[fmd]*p.m.part
            ## for individuals infected in the same month, assign
            ## the order of infection based on competing risks
            ## formula, but if the denominator is 0, replace both
            ## with 0 to avoid errors.
            p.mfirst <- p.m.exc / (p.m.exc+p.f.exc)
            p.ffirst <- 1-p.mfirst
            p.mfirst[is.na(p.mfirst)] <- 0
            p.ffirst[is.na(p.ffirst)] <- 0                
            he1e2[fmd] <- he1e2L[fmd] + p.mfirst * s..L[fmd]*p.m.exc*p.f.exc +
                                        me.L[fmd]*(1-p.f.part)*p.f.exc
            he2e1[fmd] <- he2e1L[fmd] + p.ffirst * s..L[fmd]*p.m.exc*p.f.exc +
                                        f.eL[fmd]*(1-p.m.part)*p.m.exc
            ######################################################################
            ## Iterate probabilities jointly with survival until survey.
            ## Note for probabilities of not being infected, we don't
            ## use the joint probability with being alive at sampling.
            if(survive)
              {
                mb.A[fmd] <- mb.AL[fmd]*(1-p.f.exc)*(1-p.f.part)
                me.A[fmd] <- me.AL[fmd]*(1-p.f.exc)*(1-p.f.part) + s..L[fmd]*p.m.exc.a*(1-p.f.exc)
                f.bA[fmd] <- f.bAL[fmd]*(1-p.m.exc)*(1-p.m.part)
                f.eA[fmd] <- f.eAL[fmd]*(1-p.m.exc)*(1-p.m.part) + s..L[fmd]*p.f.exc.a*(1-p.m.exc)
##              hb1b2A[fmd] <- hb1b2AL[fmd] # Doesn't change during couple duration
##              hb2b1A[fmd] <- hb2b1AL[fmd] # Doesn't change during couple duration                
                hbeA[fmd] <- hbeAL[fmd] + mb.AL[fmd]*(1-p.f.part)*p.f.exc.a
                hebA[fmd] <- hebAL[fmd] + f.bAL[fmd]*(1-p.m.part)*p.m.exc.a
                hbpA[fmd] <- hbpAL[fmd] + mb.AL[fmd]*p.f.part.a
                hpbA[fmd] <- hpbAL[fmd] + f.bAL[fmd]*p.m.part.a
                hepA[fmd] <- hepAL[fmd] + me.AL[fmd]*p.f.part.a
                hpeA[fmd] <- hpeAL[fmd] + f.eAL[fmd]*p.m.part.a
                ## for individuals infected in the same month, assign
                ## the order of infection based on competing risks
                ## formula, but if the denominator is 0, replace both
                ## with 0 to avoid errors.
                p.mfirst.a <- p.m.exc.a / (p.m.exc.a+p.f.exc.a)
                p.ffirst.a <- 1-p.mfirst.a
                p.mfirst.a[is.na(p.mfirst.a)] <- 0
                p.ffirst.a[is.na(p.ffirst.a)] <- 0                
                he1e2A[fmd] <- he1e2AL[fmd] + p.mfirst.a * s..L[fmd]*p.m.exc.a*p.f.exc.a +
                                          me.AL[fmd]*(1-p.f.part)*p.f.exc.a
                he2e1A[fmd] <- he2e1AL[fmd] + p.ffirst.a * s..L[fmd]*p.m.exc.a*p.f.exc.a +
                                          f.eAL[fmd]*(1-p.m.part)*p.m.exc.a
                ## update *L*ast month states for *A*live states
                mb.AL[fmd] <- mb.A[fmd]
                me.AL[fmd] <- me.A[fmd]
                f.bAL[fmd] <- f.bA[fmd]
                f.eAL[fmd] <- f.eA[fmd]
##              hb1b2AL[fmd] <- hb1b2A[fmd]
##              hb2b1AL[fmd] <- hb2b1A[fmd]                
                hbeAL[fmd] <- hbeA[fmd]
                hebAL[fmd] <- hebA[fmd]
                hbpAL[fmd] <- hbpA[fmd]
                hpbAL[fmd] <- hpbA[fmd]
                hepAL[fmd] <- hepA[fmd]
                hpeAL[fmd] <- hpeA[fmd]              
                he1e2AL[fmd] <- he1e2A[fmd]
                he2e1AL[fmd] <- he2e1A[fmd]                              
              }
            ## update other *L*ast month states
            s..L[fmd] <-  s..[fmd]
            mb.L[fmd] <-  mb.[fmd]
            me.L[fmd] <-  me.[fmd]
            f.bL[fmd] <-  f.b[fmd]
            f.eL[fmd] <-  f.e[fmd]
##          hb1b2L[fmd] <-  hb1b2[fmd]
##          hb2b1L[fmd] <-  hb2b1[fmd]            
            hbeL[fmd] <-  hbe[fmd]
            hebL[fmd] <-  heb[fmd]
            hbpL[fmd] <-  hbp[fmd]
            hpbL[fmd] <-  hpb[fmd]
            hepL[fmd] <-  hep[fmd]
            hpeL[fmd] <-  hpe[fmd]         
            he1e2L[fmd] <-  he1e2[fmd]
            he2e1L[fmd] <-  he2e1[fmd]                     
          }
        allstates <- data.frame(s..,mb.,me.,f.b,f.e,hb1b2,hb2b1,hbe,heb,hbp,hpb,hep,hpe,he1e2,he2e1,
                                mb.A,me.A,f.bA,f.eA,hb1b2A,hb2b1A,hbeA,hebA,hbpA,hpbA,hepA,hpeA,he1e2A,he2e1A)
        ss <- s..
        mm <- mb. + me.
        ff <- f.b + f.e
        hh <- hb1b2 + hb2b1 + hbe + heb + hbp + hpb + hep + hpe + he1e2 + he2e1
        ## Calculate probability of data given parameters * priors of paramters
        if(survive)
          {
            mmA <- mb.A + me.A
            ffA <- f.bA + f.eA
            hhA <- hb1b2A + hb2b1A + hbeA + hebA + hbpA + hpbA + hepA + hpeA + he1e2A + he2e1A
            pser.a <- cbind(hhA, mmA, ffA, ss)
          }
        pser <- cbind(hh, mm, ff, ss)
        if(trace & !sim)
          {
            ######################################################################
            ## Route of transmission breakdowns for observed couples
            ## (conditional on survival)
            ######################################################################
            ## male breakdown amongst observed M+F- couples, partner *N*egative
            pibNA <- sum(mb.A[mm.log] / mmA[mm.log]) / sum(mm.log)
            pieNA <- sum(me.A[mm.log] / mmA[mm.log]) / sum(mm.log)
            ## female breakdown amongst observed M-F+ couples, partner *N*egative
            piNbA <- sum(f.bA[ff.log] / ffA[ff.log]) / sum(ff.log)
            piNeA <- sum(f.eA[ff.log] / ffA[ff.log]) / sum(ff.log)
            ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
            pibPA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hbeA[hh.log] + hbpA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piePA <- sum((hebA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
            pipPA <- sum((hpbA[hh.log] + hpeA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
            piPbA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hebA[hh.log] + hpbA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piPeA <- sum((hbeA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piPpA <- sum((hbpA[hh.log] + hepA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
            pibUA <- (pibNA*sum(mm.log) + pibPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            pieUA <- (pieNA*sum(mm.log) + piePA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            pipUA <- (pipPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
            piUbA <- (piNbA*sum(mm.log) + piPbA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            piUeA <- (piNeA*sum(mm.log) + piPeA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            piUpA <- (piPpA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            ######################################################################
            ## Give pieNA, piNeA, piePA, piPeA, for each couple
            if(give.pis)
              {
                ## probability infection was extracouple given ser
                piCe.A <- rep(NA, K)
                piC.eA <- rep(NA, K)
                piCe.A[mm.log] <- me.A[mm.log] / mmA[mm.log]
                piCe.A[hh.log] <- (hebA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
                piC.eA[ff.log] <- f.eA[ff.log] / ffA[ff.log]
                piC.eA[hh.log] <- (hbeA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
                pis <- data.frame(piCe.A, piC.eA)
              }
            ######################################################################
            ## Route of transmission breakdowns for inferred
            ## pseudopopulation (unconditional on survival)
            ######################################################################
            ######################################################################
            ## Index infections, do with estimators summing over all
            ## infected couples and over all couples
            ## version 1 - all infected couples
            mb1. <- rowSums(allstates[!ss.log, c("mb.","hb1b2","hbe","hbp")])
            me1. <- rowSums(allstates[!ss.log, c("me.","he1e2","hep")])
            ## female
            f.b1 <- rowSums(allstates[!ss.log, c("f.b","hb2b1","heb","hpb")])
            f.e1 <- rowSums(allstates[!ss.log, c("f.e","he2e1","hpe")])
            all.infA <- rowSums(allstates[!ss.log,names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))]])
            ## number of inflated male before-couple index infections (in any couple)
            mb1.Infl <- sum( mb1. / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            me1.Infl <- sum( me1. / all.infA)
            ## nufber of inflated male before-couple index infections (in any couple)
            f.b1Infl <- sum( f.b1 / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            f.e1Infl <- sum( f.e1 / all.infA)
            ## number of inflated index infections (1 in each couple with an infection)
            IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
            ######################################################################
            ## Proportion of index infections pooling across gender
            piGb1.sumI <- mb1.Infl / IndInfl
            piGe1.sumI <- me1.Infl / IndInfl
            piG.b1sumI <- f.b1Infl / IndInfl
            piG.e1sumI <- f.e1Infl / IndInfl
            ## ## version 2 - sum over all couples
            mb1. <- rowSums(allstates[, c("mb.","hb1b2","hbe","hbp")])
            me1. <- rowSums(allstates[, c("me.","he1e2","hep")])
            ## female
            f.b1 <- rowSums(allstates[, c("f.b","hb2b1","heb","hpb")])
            f.e1 <- rowSums(allstates[, c("f.e","he2e1","hpe")])
            all.infA <- rowSums(allstates[,c("s..", names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))])])
            ## number of inflated male before-couple index infections (in any couple)
            mb1.Infl <- sum( mb1. / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            me1.Infl <- sum( me1. / all.infA)
            ## nufber of inflated male before-couple index infections (in any couple)
            f.b1Infl <- sum( f.b1 / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            f.e1Infl <- sum( f.e1 / all.infA)
            ## number of inflated index infections (1 in each couple with an infection)
            IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
######################################################################
            ## Proportion of index infections pooling across gender
            piGb1.sumIS <- mb1.Infl / IndInfl
            piGe1.sumIS <- me1.Infl / IndInfl
            piG.b1sumIS <- f.b1Infl / IndInfl
            piG.e1sumIS <- f.e1Infl / IndInfl
######################################################################
            ## put them all in a dataframe
            pop.avs <- data.frame(
                                  ## conditional on survival
                                  pibNA, pieNA, # b/e in +- given A
                                  piNbA, piNeA, # b/e in -+ given A
                                  pibPA, piePA, pipPA, # b/e/p in male in ++ given A
                                  piPbA, piPeA, piPpA, # b/e/p in female in ++ given A
                                  pibUA, pieUA, pipUA, # b/e/p in male in any given A
                                  piUbA, piUeA, piUpA, # b/e/p in female in any given A
                                  ## unconditional on survival version 1
                                  piGb1.sumI, piGe1.sumI, # b/e index in males amongst all infected 
                                  piG.b1sumI, piG.e1sumI, # b/e index in females amongst all infected
                                  piGb1.sumIS, piGe1.sumIS, # b/e index in males amongst all infected 
                                  piG.b1sumIS, piG.e1sumIS) # b/e index in females amongst all infected
            ######################################################################
            ## Project incidence forward 12 months for each of the 3
            ## couple types (ss, mh, fh) for each country in the data
            ## set (because they have different population prevalences
            num.country <- length(unique(dat$epic.ind))
            cc.inds <- unique(dat$epic.ind)
            ##  initialize proportion in each serostatus stratified by initial serostatus
            ## 
            ## proportion of concordant negative couples that will end up in each serostatus
            ss12.ssL <- rep(1, num.country)
            mm12.ssL <- rep(0, num.country)
            ff12.ssL <- rep(0, num.country)
            hh12.ssL <- rep(0, num.country)            
            ## proportion of male positive discordant couples that will end up in each serostatus
            mm12.mmL <- rep(1, num.country)
            hh12.mmL <- rep(0, num.country)            
            ## proportion of female positive discordant couples that will end up in each serostatus
            ff12.ffL <- rep(1, num.country)
            hh12.ffL <- rep(0, num.country)
            ## initialize pis
            pi.m.part12.ss <- 0 #   probability man will be infected by his female partner given couple started 12 months --
            pi.f.part12.ss <- 0 # ditto for woman
            pi.m.exc12.ss <- 0 # probability man will be infected extra-couply given couple started 12 months --
            pi.f.exc12.ss <- 0 # ditto for woman
            pi.f.part12.mm <- 0 # probability woman will be infected by her male partner given couple started 12 months M+F-
            pi.f.exc12.mm <- 0 # probability woman will be infected extra-couply given couple started 12 months M+F-           
            pi.m.part12.ff <- 0         # ditto last two for male
            pi.m.exc12.ff <- 0            # ditto last two for male
            for(tt in 1:12)             # for the next 12 months
              {
                ######################################################################
                ## Transmission probabilities
                ## probability infected extracouply in various months of 2011
                p.m.exc <- 1 - exp(-bme*epicf[1332+tt-1, cc.inds])
                p.f.exc <- 1 - exp(-bfe*epicm[1332+tt-1, cc.inds])
                ## concordant negative couples
                ss12.ss <- ss12.ssL*(1-p.m.exc)*(1-p.f.exc)
                mm12.ss <- mm12.ssL*(1-p.f.exc)*(1-p.f.part) + ss12.ssL*p.m.exc*(1-p.f.exc)
                ff12.ss <- ff12.ssL*(1-p.m.exc)*(1-p.m.part) + ss12.ssL*p.f.exc*(1-p.m.exc)
                hh12.ss <- hh12.ssL + ss12.ssL* p.m.exc*p.f.exc +
                    mm12.ssL*(p.f.part + (1-p.f.part)*p.f.exc) +
                    ff12.ssL*(p.m.part + (1-p.m.part)*p.m.exc)
                pi.m.part12.ss <- pi.m.part12.ss + ff12.ssL*p.m.part
                pi.f.part12.ss <- pi.f.part12.ss + mm12.ssL*p.f.part        
                pi.m.exc12.ss <- pi.m.exc12.ss + (ss12.ssL + ff12.ssL*(1-p.m.part))*p.m.exc
                pi.f.exc12.ss <- pi.f.exc12.ss + (ss12.ssL + mm12.ssL*(1-p.f.part))*p.f.exc
                ## male positive couples & female seroconversion
                mm12.mm <- mm12.mmL*(1-p.f.exc)*(1-p.f.part)
                hh12.mm <- hh12.mmL + mm12.mmL*(p.f.part + (1-p.f.part)*p.f.exc) 
                pi.f.part12.mm <- pi.f.part12.mm + mm12.mmL*p.f.part        
                pi.f.exc12.mm <- pi.f.exc12.mm + mm12.mmL*(1-p.f.part)*p.f.exc
                ## female positive couples & male seroconversion                  
                ff12.ff <- ff12.ffL*(1-p.m.exc)*(1-p.m.part)
                hh12.ff <- hh12.ffL + ff12.ffL*(p.m.part + (1-p.m.part)*p.m.exc)
                pi.m.part12.ff <- pi.m.part12.ff + ff12.ffL*p.m.part
                pi.m.exc12.ff <- pi.m.exc12.ff + ff12.ffL*(1-p.m.part)*p.m.exc
                ss12.ssL <- ss12.ss
                mm12.ssL <- mm12.ss
                ff12.ssL <- ff12.ss
                hh12.ssL <- hh12.ss
                ## male positive discordant
                mm12.mmL <- mm12.mm
                hh12.mmL <- hh12.mm
                ## female positive discordant
                ff12.ffL <- ff12.ff
                hh12.ffL <- hh12.ff
              }
            n.m.part.dc <- 0
            n.m.part.cc <- 0
            n.f.part.dc <- 0
            n.f.part.cc <- 0
            n.m.exc.dc <- 0
            n.m.exc.cc <- 0
            n.f.exc.dc <- 0
            n.f.exc.cc <- 0
            ## add all the different countries incidence by scaling by
            ## serotype (for West Africa pooled country analysis where
            ## multiple countries are analyzed together
            for(cc in 1:num.country)
              {
                n.m.part.dc <- n.m.part.dc + pi.m.part12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
                n.f.part.dc <- n.f.part.dc + pi.f.part12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
                n.m.part.cc <- n.m.part.cc + pi.m.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
                n.f.part.cc <- n.f.part.cc + pi.f.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
                n.m.exc.dc <- n.m.exc.dc + pi.m.exc12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
                n.f.exc.dc <- n.f.exc.dc + pi.f.exc12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
                n.m.exc.cc <- n.m.exc.cc + pi.m.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
                n.f.exc.cc <- n.f.exc.cc + pi.f.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
              }
            n.m.part.tot <- n.m.part.dc + n.m.part.cc #  total infections from partners to men
            n.f.part.tot <- n.f.part.dc + n.f.part.cc            
            n.m.exc.tot <- n.m.exc.dc + n.m.exc.cc # total infections from extra couple to men
            n.f.exc.tot <- n.f.exc.dc + n.f.exc.cc
            proj12 <- data.frame(n.m.part.dc, n.f.part.dc,
                                 n.m.part.cc, n.f.part.cc,
                                 n.m.exc.dc, n.f.exc.dc,
                                 n.m.exc.cc, n.f.exc.cc,
                                 n.m.part.tot, n.f.part.tot,
                                 n.m.exc.tot, n.f.exc.tot) / sum(!hh.log) * 1000  # incidence per 1000
            prop.exc.m <- n.m.exc.tot / (n.m.exc.tot + n.m.part.tot) # proportion of total infections that are extra couple
            prop.exc.f <- n.f.exc.tot / (n.f.exc.tot + n.f.part.tot)
            proj12 <- data.frame(proj12, prop.exc.m, prop.exc.f)            
            ## relative rate of transmission coefficient extracouply vs before relationship
            rr.m.out <- bme/bmb
            rr.f.out <- bfe/bfb
            rr.m.in <- bmp/bme
            rr.f.in <- bfp/bfe
            rr.m.pbef <- bmp/bmb        #partner to before
            rr.f.pbef <- bfp/bfb
            ## relative rate of transmission coefficient extracouply and before relationship between males and females
            rr.mf.bef <- bmb/bfb
            rr.mf.exc <- bme/bfe
            ## rho is the last one
            ## relative rate of contact/risk paramter (i.e. accounting
            ## for difference in per coital act probability as estimated
            ## from within partnership transmission.
            rr.mf.bef.cont <- rr.mf.bef * rho
            rr.mf.exc.cont <- rr.mf.exc * rho
            rrs <- data.frame(rr.m.out = rr.m.out, rr.f.out = rr.f.out,
                              rr.m.in = rr.m.in, rr.f.in = rr.f.in,
                              rr.m.pbef = rr.m.pbef, rr.f.pbef = rr.f.pbef,
                              rr.mf.bef = rr.mf.bef, rr.mf.exc = rr.mf.exc,
                              rr.mf.bef.cont = rr.mf.bef.cont, rr.mf.exc.cont = rr.mf.exc.cont)
          }
        if(sim) # if simulating data
          {
            probs <- NA
            lprob <- NA
            ## create couple state probability *A*live & *D*ead
            sim.probs <- data.frame(s..A = s..,
                                    mb.A,   mb.D   =  mb.- mb.A,
                                    me.A,   me.D   =  me. - me.A,
                                    f.bA,   f.bD   =  f.b - f.bA,
                                    f.eA,   f.eD   =  f.e - f.eA,
                                    hb1b2A, hb1b2D = hb1b2 - hb1b2A, # note some of the dead cases were infected by dead partners, so can only use the index cases in any h couple
                                    hb2b1A, hb2b1D = hb2b1 - hb2b1A,
                                    hbeA,   hbeD   =  hbe - hbeA,
                                    hebA,   hebD   =  heb - hebA,
                                    hepA,   hepD   =  hep - hepA,
                                    hpeA,   hpeD   =  hpe - hpeA,
                                    hbpA,   hbpD   =  hbp - hbpA,
                                    hpbA,   hpbD   =  hpb - hpbA,
                                    he1e2A, he1e2D = he1e2 - he1e2A,
                                    he2e1A, he2e1D = he2e1 - he2e1A)
            for(ii in 1:nrow(dat))
              {                         # sample from categorical distribution based on serostatus probabilities above
                        dat$cat[ii] <- which(rmultinom(1, 1, sim.probs[ii,])==1)
              }
            dat$cat.nm <- names(sim.probs)[dat$cat]
            dat$cat.nm <- factor(dat$cat.nm, levels = names(sim.probs))
            K <- nrow(dat)
            if(!survive) pser.a <- NA
          }else{ ## if not simulating data calculate likelihood p(data|pars)
            if(survive)
              {                         # must NORMALIZE probabilities to 1 for likelihood!
                probs <- pser.a[cbind(1:K,dat$ser)] / rowSums(pser.a) # accounting for survival
              }else{
                probs <- pser[cbind(1:K,dat$ser)] # if not accounting for survival
                pser.a <- NA
              }
            if(sum(probs==0)==0) # if non of the serotatuses occur with 0 probability in the current model
              {
                lprob <- sum(log(probs)) + dnorm(log(rho), log(trans.ratio), 1/2, log = T)
                if(length(compars)>0) clprob <- sum(log(cprobs)) + dnorm(as.numeric(compars["lrho"]),
                                                                         log(trans.ratio), 1/2, log = T)
              }else{ # if some of the serostatuses are 0, then the current parameters have 0 probability
                lprob <- -Inf
              }
          }
      }
    if(length(compars)==0)
      {
        clprob <- NA
        cprobs <- NA
      }
    if(sim) {
        if(trace) {
            if(give.pis) {
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs,  proj12=proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs))
              }else{ 
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs))
              }
          }else{
            return(list(lprob = lprob, pser.a = pser.a, pser = pser, dat = dat, sim.probs, allstates = allstates,
                        clprob = clprob, probs = probs, cprobs = cprobs))
          }
      }else{                            # if not simulating
        if(trace) {
            if(give.pis)
              {
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs,  proj12=proj12, pis = pis, allstates = allstates,
                            pser.a = pser.a, pser = pser, probs = probs))
              }else{ 
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12,
                            pser.a = pser.a, pser = pser, probs = probs))
              }
          }else{
            return(list(lprob = lprob, pser.a = pser.a, pser = pser))
          }
      }
  }

######################################################################
## Conditional multivariate normal sampling
## https://stat.ethz.ch/pipermail/r-devel/2012-March/063453.html
######################################################################
condNormal <- function(x.given, mu, sigma, given.ind, req.ind){
                                        # Returns conditional mean and variance of x[req.ind] 
                                        # Given x[given.ind] = x.given
                                        # where X is multivariate Normal with
                                        # mean = mu and covariance = sigma
                                        # 
  B <- sigma[req.ind, req.ind]
  C <- sigma[req.ind, given.ind, drop=FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
  cVar <- B - CDinv %*% t(C)
  return(list(condMean=cMu, condVar=cVar))
}

## RNG function from sigmas & marginal distribution functions with multivariate normal copulas
rcop <- function(country,          # country to simulate
                 s.epic = NA,           # epidemic curve index to use
                 sample.tmar = T,       # sample couple formation times?
                 NN,                    # total sample size if sampling tmar also
                 tmar,              # marriage month in CMC (months since 1900)
                 each,      # number of couples in each year (either one # or one # per year)
                 tint,      # interview year
                 one.couple = F, # for debugging, repeat one couple married in 1990, 100,000 times
                 browse = F)            # Debug?
  {
    if(browse)  browser()               # Debug
    sigma <- cop.sigmas[,,country]          # select country-specific normal copula covariance matrix (previously fitted in pscor.R)
    nvar <- 5                           # number of variables
    if(is.na(s.epic)) s.epic.nm <- NA
    if(sample.tmar)  { ## if sampling tmar, sample from it's distribution as well
        print(paste("Creating pseudo-population of size", NN, "for", ds.nm[country]))        
        rdat <- rmvnorm(NN, mean = rep(0,nvar), sigma = sigma) # correlated multivariate normal sample for all five variables
      }else{ ## if using inputted values in tmar & each (couple formation cohort dates, number of couples per cohort)
        ##  check that each is either length 1 (all cohorts have same number of couples) or the same length as tmar
        if(!length(each) %in% c(1,length(tmar)))    stop("Incorrect length(each)")
        if(length(each)==1) each <- rep(each, length(tmar))
        print("not sampling marriage dates so ignoring maxN and simulating cohorts")
        NN <- sum(each)                 # number of couples
        print(paste("Creating pseudo-population of size", NN, "for", ds.nm[country]))        
        tmar.cdf <- ctmar[match(tmar,ctmar$mm), country+1] ## transform tmar to [0,1] using CDF created in pscor.R
        rtmar <- qnorm(tmar.cdf)                           ## transform tmar to the normal deviate
        ## then sample from conditional normal, conditional on marital cohort
        temp <- condNormal(x.given = rtmar, mu = rep(0,5), sigma = sigma, given.ind = 5, req.ind = 1:4)
        tsig <- temp$condVar
        tmu <- t(matrix(temp$condMean, nr = 4)) # different tmu for each tmar, but tsig doesn't change
        for(tt in 1:length(tmar)) {             # for each marital cohort
            rdat.temp <- rmvnorm(each[tt], mean = tmu[tt,], sigma = tsig) # take conditional normal sample of remaining variables
            colnames(rdat.temp) <- c("ams","afs","mdur","fdur")
            if(tt==1) {                 
                rdat <- data.frame(rdat.temp, tmar = rep(rtmar[tt], each[tt]))
              }else{ rdat <- rbind(rdat, data.frame(rdat.temp, tmar = rep(rtmar[tt], each[tt])))
              }
          }                             # end marital cohort loop
      }                                 # end conditional sampling else statement
    qdat <- apply(rdat,2, pnorm)    ## get univariate normal quantiles for each relationship variable
    ## transform back using marginal copula distributions
    out <- data.frame(qdat)
    out[,] <- NA
    colnames(out) <- vars
    for(jj in 1:5) {
        print(paste("copulating", vars[jj]))
        temp <- get(cnames[jj])[,c(1,country + 1)]
        out[,jj] <- as.numeric(sapply(qdat[,jj], function(x) { temp$mm[which.min(abs(x - temp[,2]))] }))
      }
    ## Choose country epidemic curves
    if(ds.nm[country]!="WA") {
            epic.nm <- ds.nm[country]
            epic.ind <- which(colnames(epicm) == ds.nm[country])
            if(!is.na(s.epic.nm)) {
                epic.nm <- s.epic.nm
                epic.ind <- s.epic.ind
              }            
      }else{
        ## For West Africa just assume that all the marriage
        ## demographics are the same amongst countries and take each
        ## couple's country randomly proportional to the sampling
        ## breakdown in the real DHS data.
        epic.nm <- sample(odat$epic.nm[odat$group=="WA"], NN, replace = TRUE)
        epic.ind <- match(epic.nm, colnames(epicm))
      }
    ##  set up output data frame
    dat <- data.frame(uid = 1:nrow(out), ds = paste(country,"sim"), ser = NA, # unique couple identifier, data set, couple serostatus set interview
                      tms = out$tmar - out$mdur,                              # male sexual debut
                      tfs = out$tmar - out$fdur,                              # female sexual debut
                      tmar = out$tmar,                                        # couple formation month
                      tint = tint,                                            # interview month
                      mardur.mon = tint - out$tmar,                           # couple duration in months
                      circ = NA,                                              # male circumcision status
                      mage = tint - out$tmar + out$mdur + out$ams, # male age at interview
                      fage = tint - out$tmar + out$fdur + out$afs, # female age at interview
                      epic.ind = epic.ind,                         # epidemic curve index
                      epic.nm = epic.nm,                           # epidemic curve name
                      group = ds.nm[country])                      # country group name
    dat <- dat[dat$tmar<=tint,] # can't have marriage after interview time, otherwise not a couple
    if(one.couple) { ## for debugging, just take one couple that gets married in 1990 and replicate them 10^5 times
      dat <- dat[rep(which(dat$tmar==90*12)[1],nrow(dat)),]
    }
    return(dat)
  }

## Plot posterior pairwise-correlations & histograms to compare distributions
sbpairs <- function(dat, file.nm, width = 10, height = 10,
                    rgs = NA,           # ranges
                    yrg = NA,
                    breaks = NA,
                    hist.freq = T,
                    cex = 1, col = "black", nrpoints = 200, do.pdf = F, do.jpeg = T, big.ranges = F,
                    cex.axis = 1.5, cex.nm = 2.5, browse = F)
  {
    if(browse) browser()
    if(do.pdf) pdf(paste(file.nm, ".pdf", sep=""), width = width, height = height)
    if(do.jpeg) jpeg(file = paste(file.nm, ".jpg", sep=""), width = width*100, height = height*100)
    if(!is.na(rgs[1]))
      {
        ranges <- rgs
      }else{
        ranges <- apply(dat, 2, range)        
        if(big.ranges)
          {
            ranges[1,1:5] <- 0
            ranges[2,1:5] <- .06 #1.5* ranges[2,1:5]
            ranges[1,6] <- ranges[1,6] - .5
            ranges[2,6] <- ranges[2,6] + .5
          }
      }
    par(mar=rep(3,4),oma=rep(2,4),mfrow=rep(ncol(dat),2))
    parnames <- colnames(dat)
    for(ii in 1:ncol(dat))
      {
        for(jj in 1:ncol(dat))
          {
            if(ii==jj)
              {
                if(length(dim(yrg))==0)
                  {
                    hist(dat[,ii], breaks = 40, col = "black", main = "", xlab = "", ylab = "",
                     xlim = ranges[,ii], las = 2, cex.axis = cex.axis, freq = hist.freq)
                  }else{
                    hist(dat[,ii], breaks = breaks[[jj]], col = "black", main = "", xlab = "", ylab = "",
                     xlim = ranges[,ii], las = 2, cex.axis = cex.axis, freq = hist.freq, ylim = as.numeric(yrg[jj,]))
                  }
                mtext(parnames[ii], side = 3, line = -2, adj = .98, col = "red", cex = cex.nm)
              }else
            {
              smoothScatter(dat[,jj],dat[,ii], cex = cex, col = col, las = 2, cex.axis = cex.axis,
                            main = "", xlab = "", ylab = "", nrpoints = nrpoints,
                            xlim = ranges[,jj], ylim = ranges[,ii])
            }
          }
      }
    if(do.pdf | do.jpeg) dev.off()
  }

## Add ages @ sexual debut, duration of sex before marriage, & date of marriage to a data set
add.vars <- function(dat)
  {
    dat$ams <- dat$mage - (dat$tint - dat$tms) # age at sexual debut male
    dat$afs <- dat$fage - (dat$tint - dat$tfs) # age at sexual debut female
    dat$mdur <- dat$tmar - dat$tms       # sex dur before marriage male
    dat$fdur <- dat$tmar - dat$tfs       # sex dur before marriage female
    dat$ymar <- dat$tmar/12 + 1900       # year
    return(dat)
  }


## take a continuous density function & discretize it (deprecated, not used anymore)
discret <- function(dfun, range = c(0, 12*60),
                    leftz = -1,
                    rightz = 200*12,
                    redge.smooth = F, vals, # smooth right edge (needs vals to get max val)
                    browse = F) 
  {
    if(browse)  browser()
    xmin <- range[1]                    
    xmax <- range[2]
    rframe <- data.frame(l = xmin:xmax, u = (xmin+1):(xmax+1))
    nn <- nrow(rframe)
    out <- data.frame(mm = xmin:xmax, pp = NA)
    for(ii in 1:nn)
      {
        out$pp[ii] <- integrate(dfun, lower = rframe$l[ii], upper = rframe$u[ii])$value
      }
    out$pp[out$mm<leftz | out$mm>rightz] <- 0 # truncate
    if(redge.smooth)
      {
        redges <- max(vals) + 0:1
        out$pp[out$mm %in% redges] <- out$pp[out$mm == max(vals) - 1]
      }
    out$pp <- out$pp / sum(out$pp)      # normalize to 1
    return(out)
  }
