####################################################################################################
## Show transmission by time since infection using simulations that take into account increasing
## amounts of concurrency and heterogeneity.
####################################################################################################
library(plyr); library(parallel)
####################################################################################################
## Viral Load Infectivity Plots
rm(list=ls(all=T)); gc()
## Summarize Hollingsworth & Wawer style fits to simulated data
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
outdir <- file.path('results','RakAcute','TransmContr')
if(!file.exists(outdir)) dir.create(outdir)

## Age-at-seroconversion dependent Weibull survival times fit to
## CASCADE 2000 data. See Bellan et al. (2013) Supp Info for details.
ageweib <- function(age, death = T) {
  if(death) {
    shp <- 2.3                      # Weibull shape parameter
    scl <- 2000/shp/(age/12)^.53    # Weibull scale parameter
    return(round(rweibull(length(age), shape = shp, scale = scl),0))
  }else{
    return(200*12) # if there's no death arbitrarily return 200 year survival time
  }
}

trc <- function(sig.inf, vol = 1/12, mu.cont = 1, sigma.cont = .1, N = 10000, months = 12*20, sig.sus = 0, 
                rh.acute = 7, d.acute = 3, rh.late = 10, d.late = 10, d.aids = 10, bp = .007/12, pop.prev.inf = .15, browse = F) {
  if(browse)    browser()
  ## initialize het id data.frame
  if(sig.inf>0) {
    rr.inf <- exp(rnorm(N, mean = -sig.inf^2/2, sd = sig.inf))
  }else{ rr.inf <- rep(1, N) }
  if(sig.sus>0) {
    rr.sus <- exp(rnorm(N, mean = -sig.sus^2/2, sd = sig.sus))
  }else{ rr.sus <- rep(1, N) }
  idfr <- data.frame(id = 1:N, rr.inf = rr.inf, rr.sus = rr.sus)
  tts <- 1:months
  inf.arr <- array(NA, dim = c(N, length(tts), 6))
  ## ################################################…
  ## current partner contact rate, infected?, partner id, transmitted to partner? phase of infection
  dimnames(inf.arr)[[3]] <- c('cont', 'inf', 'pid', 'transm', 'ph','ph.haz') 
  inf.arr[,,c('inf','transm')] <- 0
  shp <- (mu.cont/sigma.cont)^2 ## gamma shape par
  scl <- sigma.cont^2 / mu.cont ## gamma scale parameter
  ## ################################################
  ## Set up contact rate using volatility
  ## ################################################
  ## Get risk period times using the geometric distribution
  pms <- N*months
  ## total person months divided by average risk period duration gives us approximate number of risk
  ## periods in this whole sample. Multiply X1.3 to make sure we have enough to cover random variation
  total.pers <- round(pms*vol*1.3)
  rand.pers <- 1+rgeom(total.pers, prob = 1 - exp(-vol)) ## probability of changing risk in each period is 1-exp(-volatility hazard * 1)
  rand.conts <- rgamma(length(rand.pers), shape = shp, scale = scl) ## resample contact
  rh.cont <- rep.int(rand.conts, rand.pers) ## repeat each contact for the interval  
  ## 
  ## String these together into a matrix (don't worry about dependence between last contact rate of
  ## previous individual and first of next since this is still consistent with distributional
  ## assumptions and because the infection process below resamples from such individuals
  inf.arr[,,'cont'] <- t(matrix(rh.cont[1:(N*months)], nr = months, nc = N))
  ## ################################################
  ## when do they get infected? a function of their current contact rate & their relative susceptibility
  inf.arr[,,'inf'] <- matrix(rbinom(months*N, 1, 1-exp(-bp * as.matrix(rr.sus)[,rep(1,months)] * inf.arr[,,'cont']*pop.prev.inf)), nr = N, nc = months)
  inf.log <- rowSums(inf.arr[,,'inf'])>0
  inf.t <- apply(inf.arr[inf.log,,'inf'], 1, function(x) min(which(x>0)))
  for(ii in 1:length(inf.t))     inf.arr[which(inf.log)[ii], tts> inf.t[ii],'inf'] <- 1
  ## ################################################
  ## Restrict to individuals that were infected to create transmitter array
  trarr <- inf.arr[inf.log,,]
  rr.inf <- rr.inf[inf.log]
  rr.sus <- rr.sus[inf.log]
  idfr <- idfr[inf.log,]  
  s.times <- ageweib(rep(25*12,length(inf.t)), death = TRUE) ## give them survival times based on 25 yrs old
  s.times[s.times<24] <- 24 ## make 2 year the minimum survival time
  for(ii in 1:length(inf.t))    { ## add phase data
    ## name phases
    trarr[ii,tts >= inf.t[ii] & tts < inf.t[ii] + d.acute,'ph'] <- 'acute'
    trarr[ii,tts >= inf.t[ii] + d.acute & tts < (inf.t[ii] + s.times[ii] - (d.late+d.aids)),'ph'] <- 'chronic'
    trarr[ii,tts >= (inf.t[ii] + s.times[ii] - (d.late+d.aids)) & tts < (inf.t[ii] + s.times[ii] - d.aids),'ph'] <- 'late'
    trarr[ii,tts >= (inf.t[ii] + s.times[ii] - d.aids) & tts < (inf.t[ii] + s.times[ii]),'ph'] <- 'AIDS'
    trarr[ii,tts >= (inf.t[ii] + s.times[ii]),'ph'] <- 'dead'
    ## assign hazards
    trarr[ii,tts >= inf.t[ii] & tts < inf.t[ii] + d.acute,'ph.haz'] <- rh.acute
    trarr[ii,tts >= inf.t[ii] + d.acute & tts < (inf.t[ii] + s.times[ii] - (d.late+d.aids)),'ph.haz'] <- 1
    trarr[ii,tts >= (inf.t[ii] + s.times[ii] - (d.late+d.aids)) & tts < (inf.t[ii] + s.times[ii] - d.aids),'ph.haz'] <- rh.late
    trarr[ii,tts >= (inf.t[ii] + s.times[ii] - d.aids) & tts < (inf.t[ii] + s.times[ii]),'ph.haz'] <- 0
    trarr[ii,tts >= (inf.t[ii] + s.times[ii]),'ph.haz'] <- 0
    trarr[ii,tts >= (inf.t[ii] + s.times[ii]),'inf'] <- NA
  }
  rhphs.mat <- trarr[,,'ph.haz']
  cont.mat <- trarr[,,'cont']  
  mode(rhphs.mat) <- 'numeric'
  mode(cont.mat) <- 'numeric'  
  ## phase relative hazards * monthly average hazard * individual infectivity * individual contact rate
  hazs.mat <- rhphs.mat * bp* as.matrix(rr.inf)[,rep(1,months)] * cont.mat
  ## Do they infect in a given interval?
  trarr[,,'transm'] <- matrix(rbinom(prod(dim(hazs.mat)), 1, 1-exp(-hazs.mat*1)), nr = nrow(hazs.mat), nc = ncol(hazs.mat))
  tr.wh <- which(rowSums(trarr[,,'transm']==1, na.rm=T)>0)
  ind.R <- rowSums(trarr[,,'transm']==1, na.rm=T)
  temp <- apply(trarr[tr.wh,,'transm'], 1, function(x) which(x==1)) ## day of transmission (on calendar time)
  m.trans <- mapply('-',temp, inf.t[tr.wh]) ## month of transmission (since infection date)
  m.trs <- unlist(m.trans)
  pars <- c(sig.inf=sig.inf, sig.sus=sig.sus, vol=vol, mu.cont=mu.cont, sigma.cont=sigma.cont, bp=bp)
  output <- list(m.trs = m.trs, m.trans = m.trans, ind.R = ind.R, pars = pars) #r, trarr = trarr)
  gc()
  return(output)
}

## out <- trc(N=10^4, months = 12*65, sig.inf = 0, vol = 1/24, bp=.007*20, mu.cont = 1, sigma.cont = 2, browse=F)
## range(out$m.trs)
## length(out$ind.R)
## mean(out$ind.R)
## ## Average total chronic phase equivalent hazard-months contributed
## hms <- out$trarr[,,'ph.haz']
## mode(hms) <- 'numeric'
## mean(rowSums(hms,na.rm=T))
  
louts0.1 <- mclapply(X=seq(0,3, by=.5), trc, N=10^4, months = 12*65, vol = 1/6, bp=.007*20, mu.cont = 1, sigma.cont = .1, browse=F, mc.cores = 12)
louts1 <- mclapply(X=seq(0,3, by=.5), trc, N=10^4, months = 12*65, vol = 1/6, bp=.007*20, mu.cont = 1, sigma.cont = 1, browse=F, mc.cores = 12)
louts2 <- mclapply(X=seq(0,3, by=.5), trc, N=10^4, months = 12*65, vol = 1/6, bp=.007*20, mu.cont = 1, sigma.cont = 2, browse=F, mc.cores = 12)

for(sc in c(.1,1,2)){
  nm <- paste0('test',sc,'.pdf')
  lout.temp <- get(paste0('louts',sc))
  sig.infs <- unlist(lapply(lout.temp, function(x) x$pars[1]))
  cols <- rainbow(length(lout.temp))
  pdf(file.path(outdir,nm), w = 8, h = 5)
  par(mfrow=c(1,2))
  breaks <- seq(0,500, by = 3)
  xmax <- 20
  plot(0,0, xlab = 'years since infection', ylab = 'proportion of transmission', type = 'n', xlim = c(0,xmax), ylim=c(0,.05), xaxt='n', bty='n')
  axis(1, at = seq(0,xmax, by = 5))#, 0:xmax)
  for(ii in 1:length(lout.temp)) {
    hs <- hist(lout.temp[[ii]]$m.trs, breaks = breaks, plot = F)
    lines(hs$mids/12, hs$dens, col = cols[ii])
    assign(paste0('h',ii), hs)
  }
  legend('topright', leg = signif(sig.infs, 2), col = cols, lty = 1, title = expression(sigma['infectiousness']), cex = .7, ncol=2)
  xmax <- 15
  plot(0, 0, type = 'l', ylim = c(0,1), xlim = c(0,xmax), bty = 'n', xaxt='n',
       xlab = 'years since infection', ylab = 'cumulative proportion of transmission')
  for(ii in 1:length(lout.temp)) {
    hs <- get(paste0('h',ii))
    lines(hs$mids/12, cumsum(hs$dens)*3, col = cols[ii])
  }
  axis(1, at = seq(0,xmax, by = 5))#, 0:xmax)
  dev.off()
}
