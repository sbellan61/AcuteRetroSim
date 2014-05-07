
cptime <- function(output, ## simulation output file
                   ncpl = 12,    ## of couples sub-sampled for cohort-style figure
                   inf.only = F, ## only show infected couples
                   yrmin = 1965, yrmax = 2013, ## time range
                   lwd = 1, lwd.cs = 4.4, ## line widths for outside couple & inside
                   col.mpre = 'black', col.fpre = 'black', ## pre-couple M & F SD couples colors
                   col.m = 'dark green', col.f = 'purple', ## M & F SD couples colors
                   col.ccn = 'wheat2', col.msdc = 'dark green', col.fsdc = 'purple', col.ccp = 'orange',
                   mort.pch = 4, age.pch = '|',## symbol for death
                   inf.pch.pre = 17, inf.pch.extra = 17, inf.pch.within = 17, ## symbols for outside & from-partner infections
                   inf.col.pre = 'dark gray', inf.col.extra = 'red', inf.col.within = 'blue', 
                   ylab1 = 'example \nrelationship histories', ylab2 = 'serodiscordance proportion (SDP)',
                   stacked.col = F, # depcrecated (polygons for SDP between M/F)
                   n.inf = NA, # number infected couples to have in sub-sample
                   examp = 8,   # example couple
                   tmsd = 83*12, tfsd = 88*12, tmar = 90*12, tint = 113*12, # for example couple diagram
                   show.bigsim = T, # show SDP line for full simulation (not just subsampled couples)
                   cex = 1, # symbol size
                   seed = 1, # random seed for sub-sampling couples
                   browse = F) # debug
  {
    if(browse)  browser() # debug
    step <- 0
    set.seed(seed) # set random seed
    evout <- output$evout # load line list from event-driven simulation
    if(is.na(n.inf)) { ## if not selectin # of infected couples
      if(!inf.only) { ## & not only showing infected couples
        sel <- sample(1:nrow(evout), ncpl)
      }else{ ## select all infected couples
        sel <- which(evout$msdc | evout$fsdc | evout$ccc) 
        sel <- sample(sel,ncpl)
      }
    }else{ ## select both infected & uninfected couples 
      sel.i <- which(evout$msdc | evout$fsdc | evout$ccc) 
      sel.i <- sample(sel.i,n.inf)
      sel.ni <- which(evout$ser == 4)
      sel.ni <- sample(sel.ni, ncpl - n.inf)
      sel <- c(sel.ni,sel.i)
      sel.ord <- order(evout[sel,'tmar'])
      sel <- sel[sel.ord]
    }
    ev <- evout[sel,]                   ## sub-sample data
    layout(t(matrix(1:8,nr=2,nc=4)), hei= rep(1,4)) #c(.35,.5,.35,2.5)) ## panel layout (3 rows 1 column)
####################################################################################################
    ## Prevalence figures
    epic.ind <- output$evout$epic.ind[1]
    if(!is.na(examp)) {
#      browser()
      tmsd <- ev[examp,'tms']
      tfsd <- ev[examp,'tfs']
      tmar <- ev[examp,'tmar']
      tint <- ev[examp,'tint']
    }
    ev <- rbind(ev[-examp,], ev[examp,])
    ## female
    for(ss in 1:2) {
      par(mar = c(0,5,0,.5))
#      par(xpd=F)
      sex <- c('f','m')[ss]
      osex <- c('m','f')[ss]
      ymax <- c(.2, .2)
##      sexlab <- paste('HIV prevalence in \n', c('women','men'))
            sexlab <- c('women','men')
      var1 <- get(paste0('epic',sex,'.all')) # no ART
      var2 <- get(paste0('epic',sex))        # infectious prevalence
      plot(0,0, type = 'n', xlim = c((yrmin-1900)*12, (yrmax-1900)*12), ylim = c(0, ymax[ss]),
           axes = F, bty = 'n', xlab = '', ylab = sexlab[ss])
      tts <- ((yrmin-1900)*12):((yrmax-1900)*12)
      sd <- get(paste0('t',osex,'sd'))
      predur <- sd:tmar
      mardur <- tmar:(tint-1)
      polygon(c(predur,rev(predur)), c(var2[predur,epic.ind],rep(0, length(predur))), col = 'dark gray', border = NA)
      polygon(c(mardur,rev(mardur)), c(var2[mardur,epic.ind],rep(0, length(mardur))), col = 'red', border = NA)
   #   lines(tts, var1[tts,epic.ind], col = 'black')
      lines(tts, var2[tts,epic.ind], col = 'black', lty = 1)
#      if(ss==2) axis(1, at = 12*(seq(yrmin, yrmax, by = 5)-1900), label = NA)
      axis(2, at = seq(0,c(.2,.2)[ss], by = .1), las = 2)
      if(ss==1) { ## relationship arrows
        par(mar=c(1,5,1,.5))
        plot(0,0, type = 'n', bty = 'n', axes = F, xlab='',ylab='', xlim = c((yrmin-1900)*12, (yrmax-1900)*12), ylim = c(-6,14))
        par(xpd=NA)
        mloc <- 4
        lloc <- mloc+3.5
        dar <- 15
        segments(tmsd, lloc, tmar, lloc) #male ysa
        segments(tfsd, -lloc, tmar, -lloc)#f
        segments(tmar, -lloc, tmar, lloc)# cf
        ## segments(tmar,mloc, tint, mloc) # mdur
        ## segments(tmar,-mloc, tint, -mloc) #fdur
        polygon(c(tmar,tint,tint,tmar),c(mloc,mloc,-mloc,-mloc), col = 'dark gray', border = NA)
        arlen <- .03
        arwd <- 1.5
        arrows((tmar+tmsd)/2, dar, (tmar+tmsd)/2, lloc+2, col = 'dark gray', length = arlen, code = 2, lwd = arwd)
        arrows((tmar+tint)/2, dar, (tmar+tint)/2, mloc+2, col = 'red', length = arlen, code = 2, lwd = arwd)
        arrows((tmar+tfsd)/2, -dar, (tmar+tfsd)/2, -lloc-2, col = 'dark gray', length = arlen, code = 2, lwd = arwd)
        arrows((tmar+tint)/2, -dar, (tmar+tint)/2, -mloc-2, col = 'red', length = arlen, code = 2, lwd = arwd)
        arrows((tmar+tint)/2, -mloc, (tmar+tint)/2, mloc, col = 'blue', length = arlen, code = 3, lwd = arwd)
        text((tmar+tmsd)/2, (dar + lloc + 4)/2, expression(lambda[paste('M,B')]), pos = 2, col = 'dark gray')
        text((tmar+tfsd)/2, (-dar - lloc - 4)/2, expression(lambda[paste('F,B')]), pos = 2, col = 'dark gray')
        text((tmar+tint)/2-8, (mloc + 7), expression(lambda[paste('M,E')]), pos = 3, col = 'red')        
        text((tmar+tint)/2-8, (- mloc - 7), expression(lambda[paste('F,E')]), pos = 1, col = 'red')
        text((tmar+tint)/2+30, (mloc + 2), expression(lambda[paste('M,P')]), pos = 3, col = 'blue')
        text((tmar+tint)/2+30, (- mloc - 2), expression(lambda[paste('F,P')]), pos = 1, col = 'blue')    
      }
    }
    plot(0,0, type = 'n', bty = 'n', axes = F, xlab='',ylab='', xlim = c((yrmin-1900)*12, (yrmax-1900)*12), ylim = c(-6,14),
           pch = c(mort.pch), bty = 'n')
    par(xpd=F)    
    par(mar = c(4,5,0,.5), cex.axis = cex, cex.lab = cex, lend = 1)
    ys <- 1:ncpl        ## y-locations of each couple
    mys <- 1:ncpl + .23 ## male-specific symbol/line locations on each couple (pre-couple formation)
    fys <- 1:ncpl - .23 ## female
    mysc <- 1:ncpl + .05 ## same as above but while they're in couple tighten their lines
    fysc <- 1:ncpl - .05
    plot(0,0, xlim = c((yrmin-1900)*12, (yrmax-1900)*12), # initialize plots
         ylim = c(0, length(sel)+4), xlab = "", ylab = ylab1, xaxt = "n", bty = "n", yaxt = "n")
    axis(1, at = 12*(seq(yrmin, yrmax, by = 5)-1900), label = seq(yrmin, yrmax, by = 5)) ## x-axis
    segments(ev$tms, mys, ev$tmar, mys, lwd = .8, col = 'black')# col.mpre) ## pre-couple durations
    segments(ev$tfs, fys, ev$tmar, fys, lwd = .8, col = 'black') # = col.fpre)
    segments(ev$tmar, mys, ev$tmar, fys, lwd = .8, col = 'black') # = col.mpre) ## couple fromation
#    segments(ev$tmar, mys, ev$tmar, mysc) ## couple fromation
#    segments(ev$tmar, fys, ev$tmar, fysc)    
#    segments(ev$tmar, mysc, ev$taend, mysc, col = col.m) ## relationship
#    segments(ev$tmar, fysc, ev$taend, fysc, col = col.f)
    ## couple serostatus
    sel <- ev$ser==4   ## ccn
    segments(ev$tmar[sel], ys[sel], ev$taend[sel], ys[sel], col = col.ccn, lwd = lwd.cs)
    ## mSDC
    sel <- ev$msdc & ev$mdoi>=ev$tmar     ## msdc, inf after couple formation
    segments(ev$tmar[sel], ys[sel], ev$mdoi[sel], ys[sel], col = col.ccn, lwd = lwd.cs) ## marriage to male infection
    sel <- ev$msdc & !ev$ccc & ev$mdoi>=ev$tmar     ## msdc, inf after couple formation, couple never became ++
    segments(ev$mdoi[sel], ys[sel], ev$taend[sel], ys[sel], col = col.msdc, lwd = lwd.cs)
    sel <- ev$msdc & ev$ccc & ev$mdoi>=ev$tmar     ## msdc, inf after couple formation, couple  became ++    
    segments(ev$mdoi[sel], ys[sel], ev$fdoi[sel], ys[sel], col = col.msdc, lwd = lwd.cs) # time as mSDC
    segments(ev$fdoi[sel], ys[sel], ev$taend[sel], ys[sel], col = col.ccp, lwd = lwd.cs) # time as ++
    sel <- ev$msdc & !ev$ccc & ev$mdoi<ev$tmar ## msdc, inf before couple formation, never ++
    segments(ev$tmar[sel], ys[sel], ev$taend[sel], ys[sel], col = col.msdc, lwd = lwd.cs)
    sel <- ev$msdc & ev$ccc & ev$mdoi<ev$tmar & ev$fdoi>ev$tmar ## msdc w/ female infected after marriage 
    segments(ev$tmar[sel], ys[sel], ev$fdoi[sel], ys[sel], col = col.msdc, lwd = lwd.cs)
    segments(ev$fdoi[sel], ys[sel], ev$taend[sel], ys[sel], col = col.ccp, lwd = lwd.cs)
    ## fSDC
    sel <- ev$fsdc & ev$fdoi>=ev$tmar ## fsdc, inf after couple formation
    segments(ev$tmar[sel], ys[sel], ev$fdoi[sel], ys[sel], col = col.ccn, lwd = lwd.cs)
    sel <- ev$fsdc & !ev$ccc & ev$fdoi>=ev$tmar ## fsdc, inf after couple formation, couple never became ++
    segments(ev$fdoi[sel], ys[sel], ev$taend[sel], ys[sel], col = col.fsdc, lwd = lwd.cs)
    sel <- ev$fsdc & ev$ccc & ev$fdoi>=ev$tmar ## fsdc, inf after couple formation, couple  became ++ 
    segments(ev$fdoi[sel], ys[sel], ev$mdoi[sel], ys[sel], col = col.fsdc, lwd = lwd.cs)
    segments(ev$mdoi[sel], ys[sel], ev$taend[sel], ys[sel], col = col.ccp, lwd = lwd.cs)        
    sel <- ev$fsdc & !ev$ccc & ev$fdoi<ev$tmar     ## fsdc, inf before couple formation
    segments(ev$tmar[sel], ys[sel], ev$taend[sel], ys[sel], col = col.fsdc, lwd = lwd.cs)
    sel <- ev$fsdc & ev$ccc & ev$fdoi<ev$tmar & ev$mdoi>ev$tmar ## fsdc w/ male infected after marriage 
    segments(ev$tmar[sel], ys[sel], ev$mdoi[sel], ys[sel], col = col.fsdc, lwd = lwd.cs)
    segments(ev$mdoi[sel], ys[sel], ev$taend[sel], ys[sel], col = col.ccp, lwd = lwd.cs)        
    ## ++
    sel <- ev$ccc & !ev$msdc & !ev$fsdc ## they never were discordant
    segments(ev$tmar[sel], ys[sel], ev$taend[sel], ys[sel], col = col.ccp, lwd = lwd.cs)
    ## add symbols: infections before couple formation
    cex.inf.pt <- .9
    points(ev$mdoi[ev$mdoi<=ev$tmar], mys[ev$mdoi<=ev$tmar]+.1, pch = inf.pch.pre, col = inf.col.pre, cex = cex.inf.pt, crt = 180)
    points(ev$fdoi[ev$fdoi<=ev$tmar], fys[ev$fdoi<=ev$tmar]-.2, pch = inf.pch.pre, col = inf.col.pre, cex = cex.inf.pt)
    ## add symbols: infections after couple formation from outside (extra-couple)
    points(ev$mdoi[ev$mdoi>ev$tmar & ev$mcoi =='e'], mys[ev$mdoi>ev$tmar & ev$mcoi =='e']+.1, pch = inf.pch.extra,  col = inf.col.extra, cex = cex.inf.pt, crt = 180) 
    points(ev$fdoi[ev$fdoi>ev$tmar & ev$fcoi =='e'], fys[ev$fdoi>ev$tmar & ev$fcoi =='e']-.2, pch = inf.pch.extra,  col = inf.col.extra, cex = cex.inf.pt)
    ## add symbols: infections after couple formation from inside (extra-couple)                                                      
    points(ev$mdoi[ev$mdoi>ev$tmar & ev$mcoi =='p'], mys[ev$mdoi>ev$tmar & ev$mcoi =='p']+.1, pch = inf.pch.within, col = inf.col.within, cex = cex.inf.pt, crt = 180)
    points(ev$fdoi[ev$fdoi>ev$tmar & ev$fcoi =='p'], fys[ev$fdoi>ev$tmar & ev$fcoi =='p']-.2, pch = inf.pch.within, col = inf.col.within, cex = cex.inf.pt)
    ## add symbols: deaths
    if(mort.pch == '|') {
      points(ev$mdod[ev$mdod==ev$taend], .5*(mys+fys)[ev$mdod==ev$taend], pch = mort.pch, col = col.m, cex = .8)
      points(ev$fdod[ev$fdod==ev$taend], .5*(mys+fys)[ev$fdod==ev$taend], pch = mort.pch, col = col.f, cex = .8)
    }else{
      points(ev$mdod[ev$mdod==ev$taend], .5*(mys+fys)[ev$mdod==ev$taend], pch = mort.pch, col = 'black', cex = .8) #col.m)
      points(ev$fdod[ev$fdod==ev$taend], .5*(mys+fys)[ev$fdod==ev$taend], pch = mort.pch, col = 'black', cex = .8) #col.f)
    }
    ## age
    sel <- ev$mdod!=ev$taend & ev$fdod!=ev$taend & ev$taend < yrmax
    points(ev$taend[sel], ys[sel], pch = age.pch, col = 'black')#col.m)
    par(mar = c(4,5,0,.5))
    ## legend('topleft', c('--','+-','++'), col = c(col.ccn,col.fsdc,col.ccp), bty = 'n', pch = 15, title='couple serostatus')
    ## legend('topright', leg = c('AIDS death'), col = c('black'),
    plot(0,0)
    plot(0,0)    
  }
