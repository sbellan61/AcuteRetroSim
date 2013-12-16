####################################################################################################
## Collect, summarize, & visualize results from counterfactual simulations.
####################################################################################################
rm(list=ls())                           # clear workspace
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
source('PlotFunctions.r')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
library(abind)                          # array binding
load('data files/ds.nm.all.Rdata')        # country names
load('data files/dframe.s.Rdata')        # SDP by survey
do.again <- F                           # collect results again (otherwise load cfs.Rdata)
show.pts <- T                           # show observed SDP in real DHS data
## source('CounterFactualSummaries.R')

dir.results <- file.path('results','CounterFactual') # results locations
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures
if(!file.exists(dir.figs)) dir.create(dir.figs)      # create it
## load 'blocks' which gives info on all simulations within each country-acute group
load(file.path(dir.results, 'blocks.Rdata'))
## Load all results files in the results directory (all Rdata files except blocks & cfs)
fs <- list.files(pattern = '.Rdata', path = file.path(dir.results), recursive = T, full.names = T)
fs <- fs[!grepl('blocks',fs) & !grepl('cfs',fs)]
#fs <- fs[grepl('Acute7',fs)]
print(length(fs))
    
## coll: collects all results into data frame of input parameters, and array of time series
if(!file.exists(file.path(dir.results, 'cfs.Rdata')) | do.again) {
  cfs <- coll(fs, nc = 12, give.ev = F)
  cfs$cframe <- cfs$cframe[order(cfs$cframe$job),] # order by jobs
  cfs$t.arr <- cfs$t.arr[,,order(cfs$cframe$job)]          # ditto
  attach(cfs) # for convenience, be careful later!
  save(cfs, file = file.path(dir.results, 'cfs.Rdata'))
}else{
  load(file.path(dir.results, 'cfs.Rdata'))
  attach(cfs)
}
dim(cframe) ## dimensions & check for duplicate runs
print(paste(sum(duplicated(cframe$job)), 'duplicate jobs'))
cframe <- cframe[!duplicated(cframe$job),]
t.arr <-  t.arr[,,!duplicated(cframe$job)]
acutes <- unique(cframe$acute.sc)       # which acute phase relative hazards (RHs) were simulated
nac <- 8 #length(acutes)                   # how many
col.pl <- 'black'                       # base plot color
mend <- max(blocks$end)                 # last job of each block
countries <- unique(cframe$group.ind)   # countries simulated
countries <- countries[order(countries)]
ngroup <- length(countries)             # how many
jtd <- which(!1:(mend*nac*ngroup) %in% cframe$job)    # jobs undone
std <- which(is.na(t.arr[612,4,])) # jobs that weren't run up to 2013 yet
print(paste("didn't do jobs:",paste(head(jtd,50), collapse=','))) # check to see if any jobs didn't complete
save(jtd, file='data files/CFJobsToDo.Rdata')
save(std, file='data files/CFShort2011ToDo.Rdata')

## Set labels & indices for plotting below
set.labs <- function(bb, js) {
  if(grepl('scale', blocks$lab[bb])) { # set up legend for scaled transmission routes
    route.ind <<- which(colSums(cframe[js,grepl('b...sc',names(cframe))]!=1)>0)[1]
    route <<- paste0(rep(c('pre','extra','within'),each=2),'-couple')[route.ind]
    col.ind <<- rep(c('bmb.sc','bme.sc','bmp.sc'),each=2)[route.ind]
    legtitle <<- paste(route, '\ntransmission X') # route-specific transmission scaled by
    leg <<- cframe[js,col.ind]                    # scalar values
  }
  if(grepl('heterogeneity', blocks$lab[bb])) { # set up legend for different heterogeneity
    route.ind <<- which(colSums(cframe[js,grepl('het.*sd',names(cframe))]!=0, na.rm=T)>0) # which het sd is not 0 
    route <<- paste0(c('pre','extra','within','pre, extra, within','pre & extra'),'-couple heterogeneity')[route.ind]
    col.ind <<- paste0('het.',c('b','e','p','gen','beh'),'.sd')[route.ind]
    if(length(route.ind)==3) { ## all 3 routes (individuals have different risk deviates for each route though)
      route <<- 'pre, extra, & within-couple heterogeneity \n(different risk deviates)'
      col.ind <<- col.ind[1]   # just pick first one to get legend since the sd's are the same for all routes in each sim
    }
    if(length(route.ind)==2) { ## pre-/extra- route (individuals have different risk deviates for each route though)
      route <<- 'pre- & extra-couple heterogeneity \n(different risk deviates)'
      col.ind <<- col.ind[1]   # just pick first one to get legend since the sd's are the same for all routes in each sim
    }                    
    legtitle <<- paste(route, '\nstd dev =') # route-specific standard deviation
    leg <<- cframe[js,col.ind]               # what were the sds?
  }
  if(grepl('mortality', blocks$lab[bb])) { # set up legend for AIDS mortality counterfactual
    legtitle <<- ''
    leg <<- c('as fitted', 'no AIDS mortality \ncounterfactual')
  }else{ ## below: call plotting function (from plot fxns.R)
    leg <<- ''
  }
}

ac.to.do <- c(1,7,25,50)
nac <- length(ac.to.do)
for(cc in countries) {   # make summary figures for each country
####################################################################################################
  ## Summary of all blocks: 1 page per block, 1 row per acute phase.
  print(paste0('visualizing results from ', ds.nm[cc], ' counterfactual simulations'))
  pdf(file.path(dir.figs,paste0(ds.nm[cc], 'SDP summary.pdf')), w = 8, h = 4)
  for(bb in 2:(nrow(blocks)-1)) { ## for each block
    par(mfrow = c(1,nac), oma = c(0,0,1.5,0))   # one panel per acute phase
    for(aa in 1:nac) { 
      jst <- c(1, blocks$start[bb]:blocks$end[bb])
      js <- which(cframe$simj %in% jst)
      js <- js[cframe$acute.sc[js]==ac.to.do[aa] & cframe$group.ind[js]==cc] # select sims for this acute phase RH & country
      set.labs(bb, js)
      js1 <- cframe$job[js[cframe$simj[js]==1]] # which line was as fitted? always simj=1
      plot.sdp.nsub(js = js, leg = leg, js1 = js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc,
                    main = paste('acute phase \nrelative hazard =',ac.to.do[aa]), cex.leg = .8,
                    title = legtitle, browse=F, col.pl = col.pl, show.leg = T, sep.leg = F)
      mtext(blocks$lab[bb], side = 3, outer = T, line = 0, adj = .5) # add label describing block to top of figure
    }
  }
  dev.off()
####################################################################################################
  ## Summary of transmission route-scaling counterfactuals, columns = acute RH, rows = route, page
  ## = genetic heterogeneity std dev
  pdf(file.path(dir.figs,paste0(ds.nm[cc], 'SDP scale summary big panels.pdf')), w = 8, h = 8)
  for(bb in 2:13) { ## for each block
    if(bb %in% c(2,6,10))    par(mfrow = c(4,nac), mar = c(4,4,1,.5), oma = c(0,2,4,0))   # one panel per acute phase
    for(aa in 1:nac) { 
      jst <- c(1, blocks$start[bb]:blocks$end[bb])
      js <- which(cframe$simj %in% jst)
      js <- js[cframe$acute.sc[js]==ac.to.do[aa] & cframe$group.ind[js]==cc] # select sims for this acute phase RH & country
      set.labs(bb, js)
      main <- ifelse(bb %in% c(2,6,10), paste('acute =',ac.to.do[aa]), '')
      js1 <- cframe$job[js[cframe$simj[js]==1]] # which line was as fitted? always simj=1
      plot.sdp.nsub(js = js, leg = leg, js1 = js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc,
                    main = main, cex.leg = .5,
                    title = legtitle, browse=F, col.pl = col.pl, show.leg = T, sep.leg = F)
      if(bb %in% c(2,6,10) & aa ==1) {
        mtext('AIDS death', side = 2, outer = T, line = 0, adj = .95) # add label describing block to top of figure          
        mtext('pre-couple', side = 2, outer = T, line = 0, adj = .65) # add label describing block to top of figure
        mtext('extra-couple', side = 2, outer = T, line = 0, adj = .4) # add label describing block to top of figure
        mtext('within-couple', side = 2, outer = T, line = 0, adj = .12) # add label describing block to top of figure
        if(bb==2)   mtext('scaling routes, no heterogeneity', side = 3, outer = T, line = 1, adj = .5, cex = 2) 
        if(bb==6)   mtext('scaling routes, genetic heterogeneity sd = 1', side = 3, outer = T, line = 1, adj = .5, cex = 2) 
        if(bb==10)   mtext('scaling routes, genetic heterogeneity sd = 2', side = 3, outer = T, line = 1, adj = .5, cex = 2) 
      }
    }
  }
  dev.off()
####################################################################################################
  ## Summary of heterogeneity counterfactuals, columns = heterogeneity type, rows = acute phase
  ## RH, page = genetic heterogeneity inter-partner correlation
  pdf(file.path(dir.figs,paste0(ds.nm[cc],'SDP heterogeneity summary big panels.pdf')), w = 14, h = 9)
  for(bb in 14:34) { ## for each block
    ## one panel per acute phase
    if(bb %in% c(14,21,28)) {
      layout(matrix(1:(7*nac),nac,7))
      par(mar = c(3,4,1,.5), oma = c(.5,6,5,0))
    }
    for(aa in 1:nac) { 
      jst <- c(1, blocks$start[bb]:blocks$end[bb])
      js <- which(cframe$simj %in% jst)
      js <- js[cframe$acute.sc[js]==ac.to.do[aa] & cframe$group.ind[js]==cc] # select sims for this acute phase RH & country
      set.labs(bb, js)
      js1 <- cframe$job[js[cframe$simj[js]==1]] # which line was as fitted? always simj=1
      plot.sdp.nsub(js = js, leg = leg, js1 =js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc,
                    main = '', cex.leg = .8,
                    title = legtitle, browse=F, col.pl = col.pl, show.leg = T, sep.leg = F)
      if(aa==1) {
        if(bb%%7==5) mtext('pre-couple', side = 3,  line = 0, adj = .5) # add label describing block to top of figure
        if(bb%%7==6) mtext('extra-couple', side = 3,  line = 0, adj = .5) # add label describing block to top of figure
        if(bb%%7==0) mtext('within-couple', side = 3,  line = 0, adj = .5) # add label describing block to top of figure
        if(bb%%7==1) mtext('all route \n(same)', side = 3,  line = 0, adj = .5) # add label describing block to top of figure
        if(bb%%7==2) mtext('pre/extra \n(same)', side = 3,  line = 0, adj = .5) # add label describing block to top of figure
        if(bb%%7==3) mtext('all route \n(different)', side = 3,  line = 0, adj = .5) # add label describing block to top of figure
        if(bb%%7==4) mtext('pre/extra \n(different)', side = 3,  line = 0, adj = .5) # add label describing block to top of figure
      }
      if(bb %in% c(12,19,26)) {
        if(aa==1) {
          if(bb==12)   mtext('different heterogeneity std devs, inter-partner correlation=0 ',
               side = 3, outer = T, line = 2, adj = .5, cex = 2)
          if(bb==19)   mtext('different heterogeneity std devs, inter-partner correlation=0.4',
               side = 3, outer = T, line = 2, adj = .5, cex = 2)
          if(bb==26)   mtext('different heterogeneity std devs, inter-partner correlation=0.8',
               side = 3, outer = T, line = 2, adj = .5, cex = 2)
        }
        ytext <- ifelse(bb %in% c(12,19,26), paste('acute phase \nrelative hazard =',ac.to.do[aa]), '')
        mtext(ytext, side = 2, line = 6, adj = .5)
      }
    }
  }
  dev.off()
}

####################################################################################################
## Figure 2 for the manuscript
####################################################################################################
col.pl <- 'white'
ac <- 7 ## acute phase RH to use in Figure
for(one.file in c(F,T)) {
  leg.cex <- .57
  rmp <- colorRampPalette(c("yellow","red"))
  hazm <- c('bmb.sc','bme.sc','bmp.sc')   ## for each route get simjob with as fitted, 0, 10 scalars; acuteRH=7, no heterogeneity
  inds <- matrix(NA, 4,3)
  for(bb in 1:2) inds[bb,] <- c(1,blocks$start[bb+2],blocks$end[bb+2])
  for(bb in 3) inds[bb,] <- c(1,blocks$start[bb+2]+1,blocks$end[bb+2])
  inds[4,] <- c(1,90:91)
  inds <- rbind(c(1,2,NA), inds)
  bltd <- c(2:5,17) ## blocks to show in summary figure (no AIDS mortality, 3 routes, within-couple, heterogeneity)
  tbl <- list(lab = blocks$lab[bltd], seq = list(inds[1,], inds[2,], inds[3,], inds[4,], inds[5,])) ## 5 blocks for figure
  mains <- c('A','B','C','D','E')
  mains <- c('mortality','pre-couple \ncontact coefficient','extra-couple \ncontact coefficient',
             'intrinsic \ntransmission rate', 'heterogeneity \nin transmission')
  ltys <- 1:3
  lwds <- c(3,1.2,1.2)
  cx <- .9
  if(one.file) pdf(file.path(dir.figs,paste0('Figure X - Counterfactual SummaryAc',ac,'.pdf')), w = 6.5, h = 3)
  for(cc1 in countries) {
    if(!one.file)  pdf(file.path(dir.figs,paste0('Figure X - Counterfactual Summary ', ds.nm[cc1],' Ac',ac,'.pdf')), w = 6.5, h = 3)
    layout(matrix(c(1,6,2,7,3,8,4,9,5,10,11,11),2,6),w = c(rep(1,5),.85))
    par(mar = c(3,1,2,0), oma = c(1,3,0,0), cex.lab = cx, cex.axis = cx, cex.main = cx, fg = col.pl, col.axis = col.pl,
        col.lab = col.pl, col = col.pl, col.main = col.pl)
    for(bb in 1:5) {
      jst <- c(1, tbl$seq[[bb]])
      js <- which(cframe$simj %in% jst)
      js <- js[cframe$acute.sc[js]==ac & cframe$group.ind[js]==cc1] # select sims for this acute phase RH & country
      main <- mains[bb]
      js1 <- cframe$job[js[cframe$simj[js]==1]] # which line was as fitted? always simj=1
      set.labs(bb,js)
      yaxt <- ifelse(bb==1, T, F)
      if(bb==1) cols <- c('black', 'gray')
      if(bb>1) cols <- c('black', 'gray', 'gray')
      ltys <- c(1,1,2)
      plot.sdp.nsub(js = js, leg = leg, js1 = js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc1,
                    main = main, cex.leg = .8, yaxt = yaxt, ylab = '', ltys = ltys, lwds = lwds, cols = cols,
                    title = legtitle, browse=F, col.pl = col.pl, show.leg = F, sep.leg = F)
      ## if(bb==1) legend('bottomleft', c('as fitted', 'no AIDS mortality'), lwd = 2:1, lty = 1:2, col = c('black', rmp(2)[2]), bty = 'n', cex = leg.cex)
      ## if(bb%in%2:3) legend('bottomleft', c('set to 0', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(2,1,3),
      ##                      col = c(rmp(3)[2], 'black', rmp(3)[3]), bty = 'n', cex = leg.cex)
      ## if(bb%in%4) legend('bottomleft', c('scaled X 1/10', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(2,1,3),
      ##                    col = c(rmp(3)[2], 'black', rmp(3)[3]), bty = 'n', cex = leg.cex)
      ## if(bb==5) legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = c(2,1,1), lty = 1:3,
      ##      col = c('black', rmp(3)[2:3]), bty = 'n', cex = leg.cex)
      if(bb==1) legend('bottomleft', c('as fitted', 'no AIDS mortality'), lwd = 2:1, lty = 1, col = c('black', 'gray'), bty = 'n', cex = leg.cex)
      if(bb%in%2:3) legend('bottomleft', c('set to 0', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(1,1,2),
                           col = c( 'gray', 'black', 'gray'), bty = 'n', cex = leg.cex)
      if(bb%in%4) legend('bottomleft', c('scaled X 1/10', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(1,1,2),
                         col = c( 'gray', 'black',  'gray'), bty = 'n', cex = leg.cex)
      if(bb==5) legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = c(2,1,1), lty = c(1,1,2),
           col = c('black',  'gray', 'gray'), bty = 'n', cex = leg.cex)
    }
####################################################################################################
    ## Extract SDP's from 2008 for all scenarios to plot in 2nd row of figure
    par(mar = c(3,1,.5,0))
    cols <- rainbow(length(ds.nm))
    ylab <- '' #ifelse(rr==1,'SDP','')
    xmax <- c(1,10,10,10,2)
    for(rr in 1:5) {
      yaxt <- ifelse(rr==1, T, F)
      if(rr %in% 2:4) {
        logd <- 'x'
        xlim <- c(.04, 10.5)
      }else{
        logd <- ''
        if(rr==1)  xlim <- c(-.2,1.2)     else    xlim <- c(-.2,3.2)
      }   
      plot(1,1, type = 'n', xlim = xlim, ylim = c(0,1), bty = 'n', axes = F,
           xlab = '', ylab = ylab, main = '', log = logd)
      if(rr==1)     axis(1, at = c(0,1), c('no AIDS \nmortality','as fitted'), las =1, padj = 1)
      if(rr %in% c(2:4))    {
        axis(1, at = c(.1,.2,.5,1,2,5,10), label = c('0.1','0.2','0.5','1','2','5','10'), las = 2)
        if(rr<4) axis(1, at = c(.05), '0', las = 2)
      }
      if(rr==5)             axis(1, at = 0:3)
      if(rr ==5) {
        grcol <- rgb(t(col2rgb(gray(.6), alpha = .5)),max=255)
        segments(0,0,0,1, col = grcol, lwd = 3)
      }else{
        grcol <- rgb(t(col2rgb(gray(.6), alpha = .5)),max=255)
        segments(1,0,1,1, col = grcol, lwd = 3)
      }
      if(yaxt) axis(2, at = seq(0,1,l=5), las = 2) else axis(2, at = seq(0,1,l=5), labels = NA)
      yrind <- which(t.arr[,2,1]==2008)
      for(cc in countries)   {
        if(rr==1) { 
          sel <- which(cframe$group.ind==cc & cframe$simj %in% 1:2 & cframe$acute.sc==ac)
          lines(c(1,0), (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel], col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr==2) {
          sel <- which(cframe$group.ind==cc & cframe$simj %in% blocks$start[bltd[rr]]:blocks$end[bltd[rr]] & cframe$acute.sc==ac)
          xs <- cframe$bmb.sc[sel]
          xs[cframe$bmb.sc[sel]==0] <- .05 # since on a log scale
          lines(xs, (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr==3) {
          sel <- which(cframe$group.ind==cc & cframe$simj %in% blocks$start[bltd[rr]]:blocks$end[bltd[rr]] & cframe$acute.sc==ac)
          xs <- cframe$bme.sc[sel]
          xs[cframe$bme.sc[sel]==0] <- .05 # since on a log scale
          lines(xs, (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr %in% 4) {
          sel <- which(cframe$group.ind==cc & cframe$simj %in% (blocks$start[bltd[rr]]+1):blocks$end[bltd[rr]] & cframe$acute.sc==ac) 
          lines(cframe$bmp.sc[sel], (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel], col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr == 5) {
          sel <- which(cframe$group.ind==cc & cframe$simj %in% c(1,90:92) & cframe$acute.sc==ac)
          lines(cframe$het.gen.sd[sel], (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                col = cols[cc], type = 'l', pch = 19, lty = 1)
        } 
      }
    }
    mtext('serodiscordant proportion (SDP)', side = 2, line = 2, adj = .5, cex = leg.cex, outer = T)
    mtext('scalar multiple of fitted parameter used', side = 1, line = -.3, adj = .35, cex = leg.cex, outer = T)
    mtext('standard deviation \nof risk distribution', side = 1, line = 3, adj = .5, cex = leg.cex, outer = F)      
    plot.new()
    par(mar=c(0,0,0,0))
    sel <- which(cframe$simj==1 & cframe$acute.sc==ac)
    ord <- rev(order((t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel]))
    save(ord, file='data files/sdp ord.Rdata')
    par(xpd=NA)
    legend('bottomleft', ds.nm[ord], col=rainbow(length(sel))[ord], lwd = 1, title = 'top to bottom', cex = .8, inset = -.1)
    if(one.file) mtext(ds.nm[cc1], adj = 0.5, line = -2, side = 3, outer = F, cex = .75)
    if(!one.file) dev.off()
  }
  if(one.file)    dev.off()
}
graphics.off()



####################################################################################################
## Figure for appendix showing how SDP sensitivity to pre- & extra-couple is affected by heterogeneity
####################################################################################################
## 3 row by 5 column figure, 1st column shows distribution of within-couple risk factors
## other 4 columns have mortality, pre, extra, transm rate with different lines for each acute phase
## 3 rows have 'genetic heterogeneity' with std dev of 0, 1, 2
###################################################################### 
## Just do for Zambia
##
## plotting parameters
leg.cex <- .57
rmp <- colorRampPalette(c("orange","purple"))
cols <- rmp(4)
hazm <- c('bmb.sc','bme.sc','bmp.sc') ## for each route get simjob with as fitted, 0, 10 scalars; acuteRH=7, no heterogeneity
mains <- LETTERS
mains <- c('AIDS mortality','pre-couple \ncontact coefficient','extra-couple \ncontact coefficient','intrinsic \ntransmission rate')
cx <- .9
acutes.to.do <- c(1,7,25,50)
sds <- c(0,1,2)
col.pl <- 'black'
for(one.file in c(F,T)) {
####################################################################################################
    if(one.file) pdf(file.path(dir.figs,paste0('Figure AX - SDP sensitivity vs acute & het.pdf')), w = 6.5, h = 4)  
    for(cc in 1:length(ds.nm)) { #which(ds.nm=='Swaziland')) { ## each country has its own page
        if(!one.file) pdf(file.path(dir.figs,paste0('Figure AX - ',ds.nm[cc],' SDP sensitivity vs acute & het.pdf')), w = 6.5, h = 4)  
        ## initialize array telling which simulations to plot
                                        # starts <- array(NA, dim= c(length(sds),length(acutes.to.do),length(mains)))    # sds by acutes by route
                                        #ends <- starts
        starts <- matrix(blocks$start[2:13], nr = length(mains), nc = length(sds)) # routes by het sd
        starts[1,] <- c(1, 32, 57) # just take scalar 1 from pre-couple for sd = 1,2 (simj=32,57) for as fitted but with het
        ends <- matrix(blocks$end[2:13], nr = length(mains), nc = length(sds)) # routes by het sd
####################################################################################################
        ## Extract SDP's from 2008 for all scenarios to plot in 2nd row of figure
        par(mfrow=c(length(sds),length(mains)), mar = c(1.5,2,.5,0), oma = c(3,3,2,4))
        ylab <- '' #ifelse(rr==1,'SDP','')
        xmax <- c(1,10,10,10,2)
        for(sd.i in 1:length(sds)) { ## sd
            for(rr in 1:length(mains)) { ## routes
                ## initialize plot
                if(rr %in% 2:4) {
                    logd <- 'x'
                    xlim <- c(.04, 10.5)
                }else{
                    logd <- ''
                    if(rr==1)  xlim <- c(-.2,1.2)     else    xlim <- c(-.2,3.2)
                }   
                plot(1,1, type = 'n', xlim = xlim, ylim = c(0,1), bty = 'n', axes = F,
                     xlab = '', ylab = ylab, main = '', log = logd)
                main <- ifelse(sd.i==1, mains[rr],'')
                mtext(main, side = 3, adj = .6, line = 0, cex = .75)
                if(rr==length(mains)) mtext(paste0(sds[sd.i]), side = 4, line = 1.5, outer = F, las = 2)
                if(sd.i==length(sds)) {            # if last row
                    if(rr==1)     axis(1, at = c(0,1), c('no AIDS \nmortality','as fitted'), las =1, padj = 1)
                    if(rr %in% c(2:4))    {
                        axis(1, at = c(.1,.2,.5,1,2,5,10), label = c('0.1','0.2','0.5','1','2','5','10'), las = 2)
                        if(rr<4) axis(1, at = c(.05), '0', las = 2)
                    }
                }else{                            # otherwise no labels
                    if(rr==1)     axis(1, at = c(0,1), label = NA)
                    if(rr %in% c(2:4))    {
                        axis(1, at = c(.1,.2,.5,1,2,5,10), label = NA)
                        if(rr<4) axis(1, at = c(.05), label = NA)
                    }
                }
                if(rr==1) {                       # if first column
                    axis(2, at = seq(0,1,l=5), label = c(0,NA,0.5,NA,1),las = 2)
                }else{
                    axis(2, at = seq(0,1,l=5), label = NA)
                }  
                for(aa.i in 1:length(acutes.to.do)) { ## each acute phase
                    if(rr==1) {
                                        # sel <- which(cframe$job %in% 
                        sel <- which(cframe$simj %in% c((starts[rr,sd.i]),ends[rr,sd.i]) & cframe$group.ind==cc & cframe$acute.sc==acutes.to.do[aa.i])
                        if(sd.i>1)      sel <- sel[rev(order(sel))]
                        lines(c(1,0), cframe$sdp08[sel], col = cols[aa.i], type = 'l', pch = 19, lty = 1)
                    }
                    if(rr==2) {
                        sel <- which(cframe$simj %in% c(starts[rr,sd.i]:ends[rr,sd.i]) & cframe$group.ind==cc & cframe$acute.sc==acutes.to.do[aa.i])
                        xs <- cframe$bmb.sc[sel]
                        xs[cframe$bmb.sc[sel]==0] <- .05 # since on a log scale
                        lines(xs, cframe$sdp08[sel], col = cols[aa.i], type = 'l', pch = 19, lty = 1)
                    }
                    if(rr==3) {
                        sel <- which(cframe$simj %in% c(starts[rr,sd.i]:ends[rr,sd.i]) & cframe$group.ind==cc & cframe$acute.sc==acutes.to.do[aa.i])
                        xs <- cframe$bme.sc[sel]
                        xs[cframe$bme.sc[sel]==0] <- .05 # since on a log scale
                        lines(xs, cframe$sdp08[sel], col = cols[aa.i], type = 'l', pch = 19, lty = 1)
                    }
                    if(rr %in% 4) {
                        sel <- which(cframe$simj %in% c(starts[rr,sd.i]:ends[rr,sd.i]) & cframe$group.ind==cc & cframe$acute.sc==acutes.to.do[aa.i])
                        lines(cframe$bmp.sc[sel], cframe$sdp08[sel], col = cols[aa.i], type = 'l', pch = 19, lty = 1)
                    }
                } 
            }
        }
        legend('bottomleft', leg= acutes.to.do, lty = 1, title = 'acute phase \nrelative hazard', bty = 'n', cex = .7, col = cols, ncol=2)
        mtext('serodiscordant proportion (SDP)', side = 2, line = 1.5, outer = T, adj = .5, cex = .8)
        mtext('scalar multiple of fitted parameter used', side = 1, line = 1.5, adj = .75, cex = leg.cex, outer = T)
        mtext('std dev', side = 3, outer = T, line = .5, adj = 1.09, cex = .8)
        if(one.file)    mtext(ds.nm[cc], side = 3, adj=-.05, outer = T, line = .5)
        if(!one.file) dev.off()
    }
    graphics.off()
}


####################################################################################################
## Figure for appendix showing how SDP sensitivity is affected by different kinds of heterogeneity
####################################################################################################
leg.cex <- .57
rmp <- colorRampPalette(c("yellow","red"))
hazm <- c('bmb.sc','bme.sc','bmp.sc') ## for each route get simjob with as fitted, 0, 10 scalars; acuteRH=7, no heterogeneity
mains <- LETTERS
mains <- c('pre-couple \ncontact coefficient','extra-couple \ncontact coefficient','intrinsic \ntransmission rate')
ltys <- 1:3
lwds <- c(2,1,1)
cx <- .9
ac <- 7 ## acute phase RH to use in Figure
col.pl <- 'black'
for(one.file in c(F)) {
  leg.cex <- .57
  rmp <- colorRampPalette(c("yellow","red"))
  inds <- matrix(NA, 3,2)
  for(bb in 1:3) inds[bb,] <- c(c(77,81,89)[bb], c(80,84,92)[bb])
  ## bltd <- c(2:5,17) ## blocks to show in summary figure (no AIDS mortality, 3 routes, within-couple, heterogeneity)
  ## tbl <- list(lab = blocks$lab[bltd], seq = list(inds[1,], inds[2,], inds[3,], inds[4,], inds[5,])) ## 5 blocks for figure
  ltys <- 1:3
  lwds <- c(3,1.2,1.2)
  cx <- .9
  col.pl <- 'black'
  if(one.file) pdf(file.path(dir.figs,paste0('Figure X - Counterfactual Heterogeneity SummaryAc',ac,'.pdf')), w = 6.5, h = 3)
  for(cc1 in which(ds.nm=='Zambia')) {
    if(!one.file)  pdf(file.path(dir.figs,paste0('Figure X - Counterfactual Heterogeneity Summary ', ds.nm[cc1],' Ac',ac,'.pdf')), w = 6.5, h = 3)
    layout(matrix(c(1,4,2,5,3,6,7,7), nr = 2, nc = 4), ,w = c(rep(1,3),.55))
    par(mar = c(3,1,2,0), oma = c(1,3,0,0), cex.lab = cx, cex.axis = cx, cex.main = cx)
    for(bb in 1:3) {
      jst <- c(1, inds[bb,1]+1,inds[bb,2]-1)
      js <- which(cframe$simj %in% jst)
      js <- js[cframe$acute.sc[js]==ac & cframe$group.ind[js]==cc1] # select sims for this acute phase RH & country
      main <- mains[bb]
      js1 <- cframe$job[js[cframe$simj[js]==1]] # which line was as fitted? always simj=1
      set.labs(bb,js)
      yaxt <- ifelse(bb==1, T, F)
      cols <- c('black', 'gray', 'gray')
      ltys <- c(1,1,2)
      plot.sdp.nsub(js = js, leg = leg, js1 = js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc1,
                    main = main, cex.leg = .8, yaxt = yaxt, ylab = '', ltys = ltys, lwds = lwds, cols = cols,
                    title = legtitle, browse=F, col.pl = col.pl, show.leg = F, sep.leg = F)
      ## if(bb==1) legend('bottomleft', c('as fitted', 'no AIDS mortality'), lwd = 2:1, lty = 1:2, col = c('black', rmp(2)[2]), bty = 'n', cex = leg.cex)
      ## if(bb%in%2:3) legend('bottomleft', c('set to 0', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(2,1,3),
      ##                      col = c(rmp(3)[2], 'black', rmp(3)[3]), bty = 'n', cex = leg.cex)
      ## if(bb%in%4) legend('bottomleft', c('scaled X 1/10', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(2,1,3),
      ##                    col = c(rmp(3)[2], 'black', rmp(3)[3]), bty = 'n', cex = leg.cex)
      ## if(bb==5) legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = c(2,1,1), lty = 1:3,
      ##      col = c('black', rmp(3)[2:3]), bty = 'n', cex = leg.cex)
      if(bb==3) legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = c(2,1,1), lty = c(1,1,2),
           col = c('black',  'gray', 'gray'), bty = 'n', cex = leg.cex)
    }
####################################################################################################
    ## Extract SDP's from 2008 for all scenarios to plot in 2nd row of figure
    par(mar = c(3,1,.5,0))
    cols <- rainbow(length(ds.nm))
    ylab <- '' #ifelse(rr==1,'SDP','')
    xmax <- c(1,10,10,10,2)
    for(rr in 1:3) {
      het.type <- paste0('het.',c('b','e','gen'),'.sd')[rr]
      yaxt <- ifelse(rr==1, T, F)
      xlim <- c(-.2,3.2)
      plot(1,1, type = 'n', xlim = xlim, ylim = c(0,1), bty = 'n', axes = F,
           xlab = '', ylab = ylab, main = '')
      axis(1, at = 0:3)
      grcol <- rgb(t(col2rgb(gray(.6), alpha = .5)),max=255)
      segments(0,0,0,1, col = grcol, lwd = 3)
      if(yaxt) axis(2, at = seq(0,1,l=5), las = 2) else axis(2, at = seq(0,1,l=5), labels = NA)
      for(cc in countries)   {
        sel <- which(cframe$group.ind==cc & cframe$simj %in% c(1,inds[rr,]) & cframe$acute.sc==ac)
        lines(cframe[sel,het.type], cframe$sdp08[sel], col = cols[cc], type = 'l', pch = 19, lty = 1)
      } 
    }
    mtext('standard deviation of risk distribution', side = 1, line = -.5, adj = .4, cex = leg.cex, outer = T)      
    plot.new()
    par(mar=c(0,0,0,0))
    sel <- which(cframe$simj==1 & cframe$acute.sc==ac)
    ord <- rev(order(cframe$sdp08[sel]))
    save(ord, file='data files/sdp ord.Rdata')
    par(xpd=NA)
    legend('bottomleft', ds.nm[ord], col=rainbow(length(sel))[ord], lwd = 1, title = 'top to bottom', cex = .8, inset = -.1)
    mtext('serodiscordant proportion (SDP)', side = 2, line = 1.8, outer = T, adj = .5, cex = .75)
    if(one.file) mtext(ds.nm[cc1], adj = 0.5, line = -2, side = 3, outer = F, cex = .75)
    if(!one.file) dev.off()
  }
  if(one.file)    dev.off()
}
graphics.off()
