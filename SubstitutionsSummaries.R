rm(list=ls())
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
source('PlotFunctions.r')
library(abind)
library(lattice)
library(multicore)
load('data files/epic.Rdata')
load('data files/ds.nm.all.Rdata')
load('data files/dframe.s.Rdata')
load('data files/pars.arr.Rdata')
load(file='data files/sdp ord.Rdata')
outdir <- file.path('results','SubstitutionAnalysis')
## source('SubstitutionsSummaries.R')

## create some data frame to store difference in SDP based on substitution
show <- c('group.ind',"s.epic", "s.demog", "s.bmb", "s.bfb", "s.bme", "s.bfe", "s.bmp", "s.bfp")
styps <- c('s.bmb','s.bme','s.bmp','s.epic','s.demog')
slabs <- c('pre-couple \ncontact coefficient','extra-couple \ncontact coefficient','transmission \nrate','prevalence','relationship \npatterns', 'summary')

                                        #s.labs <- LETTERS[1:6]                  # panel labeling by letters          
drs <- list.files(outdir, full.names = T)[grepl('Acute', list.files(outdir))] # directories with results

drs <- drs[1]
drs <- drs[c(1,8,3,7)]
dd <- drs[1]

fs <- list.files(dd, recursive = T, full.names = T, pattern='.Rdata') # find results
fs <- fs[!grepl('cfs',fs)]                   # make sure we're not grabbing a cfs (summary) file)
cfs <- coll(fs, give.ev = F, nc = 12)
attach(cfs)
already.done <- cframe$job
save(already.done, file=file.path('data files/SubAlreadyDone.Rdata'))

for(dd in drs) {
        s.country <- as.numeric(data.frame(apply(cframe[,show],1,unique))[2,]) # substitution country #'s
        sframe <- data.frame(cframe[,show], s.country = s.country,
                             ## time, SDP data subst country, SDP sim base country, SDP sim base country with sub input
                             tt.dat.s1 = NA, sdp.dat.s1 = NA, sdp.sim.b1 = NA, sdp.sim.s1 = NA, # 1st survey
                             tt.dat.s2 = NA, sdp.dat.s2 = NA, sdp.sim.b2 = NA, sdp.sim.s2 = NA, # 2nd survey
                             tt.dat.s3 = NA, sdp.dat.s3 = NA, sdp.sim.b3 = NA, sdp.sim.s3 = NA, # 2nd survey                             
                             x1 = NA, y1 = NA, x2 = NA, y2 = NA, x3 = NA, y3 = NA,
                             uvar = NA, evar = NA)
        for(ii in 1:nrow(sframe)) {      # fill it in; *b*ase countries (recipient country)
            tt.dat.s <- dframe.s[dframe.s$group == sframe$s.country[ii],'yr'] # time, data, substitution country (1-2 vals)
            sdp.dat.s <- dframe.s[dframe.s$group == sframe$s.country[ii],'psdc'] # SDP, data, sub country
            sframe$tt.dat.s1[ii] <- as.numeric(which.min(abs(t.arr[,'yr',ii] - tt.dat.s[1]))) # add 1st time
            tt <- sframe$tt.dat.s1[ii] ## not in cmc but in t.arr row #'s
            sframe$sdp.dat.s1[ii] <- sdp.dat.s[1]
            base.temp <- which(rowSums(cframe[,styps]!=sframe$group.ind[ii])==0) ## which styps are all base country
            sframe$sdp.sim.b1[ii] <- (t.arr[tt,'mm',base.temp] + t.arr[tt,'ff',base.temp]) / t.arr[tt,'inf.alive',base.temp] # sdp, sim, base
            sframe$sdp.sim.s1[ii] <- (t.arr[tt,'mm',ii] + t.arr[tt,'ff',ii]) / t.arr[tt,'inf.alive',ii]
            if(length(tt.dat.s)==2) { ## if there were two survey years in substitution (donor) country
                sframe$tt.dat.s2[ii] <- as.numeric(which.min(abs(t.arr[,'yr',ii] - tt.dat.s[2]))) # add 2st time
                tt <- sframe$tt.dat.s2[ii]
                sframe$sdp.dat.s2[ii] <- sdp.dat.s[2]
                sframe$sdp.sim.b2[ii] <- (t.arr[tt,'mm',base.temp] + t.arr[tt,'ff',base.temp]) / t.arr[tt,'inf.alive',base.temp] # sdp, sim, base
                sframe$sdp.sim.s2[ii] <- (t.arr[tt,'mm',ii] + t.arr[tt,'ff',ii]) / t.arr[tt,'inf.alive',ii]
            }
            if(length(tt.dat.s)==3) { ## if there were two survey years in substitution (donor) country
                sframe$tt.dat.s3[ii] <- as.numeric(which.min(abs(t.arr[,'yr',ii] - tt.dat.s[3]))) # add 3st time
                tt <- sframe$tt.dat.s3[ii]
                sframe$sdp.dat.s3[ii] <- sdp.dat.s[3]
                sframe$sdp.sim.b3[ii] <- (t.arr[tt,'mm',base.temp] + t.arr[tt,'ff',base.temp]) / t.arr[tt,'inf.alive',base.temp] # sdp, sim, base
                sframe$sdp.sim.s3[ii] <- (t.arr[tt,'mm',ii] + t.arr[tt,'ff',ii]) / t.arr[tt,'inf.alive',ii]
            }
            ## X1 (distance from base country simulation to substition country data)
            sframe$x1[ii] <- sframe$sdp.sim.b1[ii] - sframe$sdp.dat.s1[ii]
            ## Y1 (distance from base country simulation with susbtituted input to substition country data)            
            sframe$y1[ii] <- sframe$sdp.sim.s1[ii] - sframe$sdp.dat.s1[ii]            
            if(length(ts)==3) {           # if 3 surveys
                ## X2 (distance from base country simulation to substition country data)
                sframe$x2[ii] <- sframe$sdp.sim.b2[ii] - sframe$sdp.dat.s2[ii]
                ## Y2 (distance from base country simulation with susbtituted input to substition country data)            
                sframe$y2[ii] <- sframe$sdp.sim.s2[ii] - sframe$sdp.dat.s2[ii]
                ## X3 (distance from base country simulation to substition country data)
                sframe$x3[ii] <- sframe$sdp.sim.b3[ii] - sframe$sdp.dat.s3[ii]
                ## Y3 (distance from base country simulation with susbtituted input to substition country data)            
                sframe$y3[ii] <- sframe$sdp.sim.s3[ii] - sframe$sdp.dat.s3[ii]
                ## unexplained variance: gemoetric mean
                sframe$uvar[ii] <- (sframe$y1[ii]*sframe$y2[ii]*sframe$y3[ii])^(2/3) / (sframe$x1[ii]*sframe$x2[ii]*sframe$x2[ii])^(2/3)
            }else{
                if(length(ts)==2) {           # if 2 surveys
                    ## X2 (distance from base country simulation to substition country data)
                    sframe$x2[ii] <- sframe$sdp.sim.b2[ii] - sframe$sdp.dat.s2[ii]
                    ## Y2 (distance from base country simulation with susbtituted input to substition country data)            
                    sframe$y2[ii] <- sframe$sdp.sim.s2[ii] - sframe$sdp.dat.s2[ii]
                    ## unexplained variance: gemoetric mean of Y^2/X^2: (Y1^2/X1^2 * Y2^2/X2^2)^1/2 = (Y1*Y2) / (X1*X2)
                    sframe$uvar[ii] <- (sframe$y1[ii]*sframe$y2[ii]) / (sframe$x1[ii]*sframe$x2[ii])
                }else{                    # 1 survey
                    sframe$uvar[ii] <- (sframe$y1[ii] / sframe$x1[ii])^2 # (Y/X)^2
                }
            }
        }
        sframe$evar <- 1-sframe$uvar
        sfr <- sframe[,c(show,'s.country','uvar','evar')]
        sfr$var <- show[-1][apply(sfr[,show[-1]] == sfr[,1], 1, function(x) as.numeric(which(!x)[1]))]
        sfr$var <- factor(sfr$var)
        sfr <- sfr[sfr$s.country != sfr$group.ind,] # remove same country sims
        levels(sfr$var) <- c('pre-couple ','extra-couple ','transmission rate','relationship \npatterns','prevalence')
        sfr$var <- ordered(sfr$var, levels=c('pre-couple ','extra-couple ','transmission rate','prevalence','relationship \npatterns'))
        save(cfs, sfr, sframe, file = file.path(dd, 'cfs.Rdata'))
#        detach(cfs)
    }

## load(file = file.path(dd, 'cfs.Rdata'))
## attach(cfs)

####################################################################################################
## Plot Figures
col.pl <- 'white'
#col.pl <- 'black'
figdir <- file.path(outdir,'SubFigures')
if(!file.exists(figdir))        dir.create(figdir)
show.sdp <- T
for(one.file in c(F,T)) {
    ncount <- length(ds.nm)
    for(dd in drs) {# for each assumed acute phase
        ## ## collect results from directory
        if(one.file) pdf(file.path(figdir,paste('SubsAc',cframe$acute.sc[1],'Pts'[show.sdp],'.pdf', sep = '')),
                         w = 6.5, h = 4.5)
        for(gg in c(1:ncount)[-14]) {              # for each country   
            country <- ds.nm[gg]                        # country name
            ## initialize substitution data frame for gg-th country on dd-th simulation group
            if(!one.file) pdf(file.path(figdir,paste(ds.nm[gg],'SubsAcute',cframe$acute.sc[1],'Pts'[show.sdp],'.pdf', sep = '')),
                              w = 6.5, h = 4.5)
            layout(matrix(c(1,4,2,5,3,6,7,7),nr=2,nc=4), w = c(1,1,1,.6))
            par(mar = c(3,1.5,4,.3), oma = c(0,2.7,0,0), fg = col.pl, col.axis = col.pl, col.lab=col.pl, col.main = col.pl, col = col.pl)
            col <- col.pl
            leg <- ds.nm
            for(styp.i in 1:5) { ## for each substitution type, plot using plot.sdp()
                styp <- styps[styp.i]              
                js <- which(cframe[,styp] != cframe$group & cframe$group.ind==gg)
             js1 <- which(rowSums(cframe[,styps]!=gg)==0)
                js <- c(js1,js)
                main <- paste0(LETTERS[1:6], ') ',slabs)[styp.i]
                main <- slabs[styp.i]
                plot.sdp(js = js,  repl = styp, leg = leg,                      
                         main=main, ylab = '',
                         yaxt = (styp.i %in% c(1,4)), show.sdp = show.sdp, make.pdf = F, early.yr = 1985, browse=F,
                         cols =NA, col.pl = col, sep.leg = F, show.leg=F, cex.leg = .4, leg.col = 2)
    #            if(styp.i==3) {     ## show X & Y for explained variance calculations
                    par(xpd=NA)
                    sfr.styp <- sframe[js,]
                    arrows(2013.5, sfr.styp[sfr.styp[,styp]==3,'sdp.sim.b1'],                               
                           2013.5, sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'],
                           len = .01, angle = 90, code = 3, col = 'black')
                    text(x=2013.5, y=sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'], labels='X', col = 'black', pos = 3)
                    arrows(2015.5, sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'],                               
                           2015.5, sfr.styp[sfr.styp[,styp]==3,'sdp.sim.s1'],
                           len = .01, angle = 90, code = 3, col = rainbow(length(ds.nm))[3])
                    text(x=2015.5, y=sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'], labels='Y', col=rainbow(length(ds.nm))[3], pos = 3)
    #            }
            }
            ##  summary plot ready for all countries to put
            col <- col.pl          # plot color
            ymin <- -1
            ymax <- 1
            ylim <- c(ymin,ymax)
            par(col = col, fg = col, col.axis = col, col.main = col, col.lab=col, bty = 'l', mar = c(3,4,4,1.5))
            frame()
            plot.window(xlim = c(.58,6), ylim = ylim)
            grcol <- rgb(t(col2rgb(gray(.5))), alpha = 250,max=255)
            wher <- 5.5
            segments(.7,1,wher,1, col= grcol, lty = 1, lwd = 3)
            ## text(5.5,1, 'all variance \nexplained', pos=4)
            ## text(5.5,0, 'no variance \nexplained', pos=4)
            main <- paste0(LETTERS[1:6], ') ',slabs)[6]
            main <- slabs[6]
            segments(.7,0,wher,0, col = grcol, lty = 1, lwd = 3)
            boxplot(evar ~ var, data = sfr, bty = 'n', ylim = c(-.25,1), ylab = '', add = T, range = 1.5,
                    las = 2, axes = F, col = 'orange', lty = 1, main = main,
            outline = F) # boxplot of explained variance
            axis(1, at = 1:5, LETTERS[1:5], pos = ymin-.1, lwd = 1, las = 1)
            axis(2, at = seq(ymin,ymax, by = .5), labels = NA,mgp = c(3,1,-.2))
            axis(2, at = seq(ymin,ymax, by = 1), las = 2, mgp = c(3,1,-.2))
            mtext(expression(paste('deviance explained 1 - ', (Y^2/X^2))), side = 2, line = 1.5, adj = .5, cex = .7)
            par(xpd=NA)
          text(wher,1, 'all', las = 1, pos = 4) 
           text(wher,.5, 'deviance decreased', las = 2, pos = 4)
           text(wher,-.7, 'deviance increased', las = 2, pos = 4)
          text(wher,0, 'none', las = 1, pos = 4)            
                                        #mtext('serodiscordant proportion', side = 2, outer = T, adj = .5, cex = 1.5)
            par(mar=c(.2,0,0,0))
            plot.new()
            ord <- rev(order(cframe[which(rowSums(cframe[,styps]!=cframe$group.ind)==0),'sdp08']))
            cols <- rainbow(length(ds.nm))
            cols[gg] <- 'black'
            lwds <- rep(1,length(ds.nm))
            lwds[gg] <- 2
            legend('top', ds.nm[ord], col=cols[ord], lwd = lwds[ord], title = 'top to bottom', cex = .85,
                   inset = .1) #, bty = 'n')
            mtext('serodiscordant proportion (SDP)', side = 2.4, outer = T, adj = .5, line = 1)
            if(!one.file) dev.off()
        }                               ## loop through countries
    } ## end country loop
    if(one.file) dev.off()
    ## Summarize all donor-recipient pair relationships by substituted input
    ## only look at relevant columns
}

aggregate(sfr$evar, list(sfr$var), median)


####################################################################################################
## same figure but by donor country instead of recipient
####################################################################################################
figdir <- file.path(outdir,'SubFigures')
if(!file.exists(figdir))        dir.create(figdir)
show.sdp <- T
for(one.file in F) { #c(F,T)) {
    ncount <- length(ds.nm)
    for(dd in drs) {# for each assumed acute phase
        ## ## collect results from directory
        if(one.file) pdf(file.path(figdir,paste('Donor SubsAc',cframe$acute.sc[1],'Pts'[show.sdp],'.pdf', sep = '')),
                         w = 6.5, h = 4.5)
        for(gg in 15) { #c(1:ncount)[-14]) {              # for each country   
            country <- ds.nm[gg]                        # country name
            ## initialize substitution data frame for gg-th country on dd-th simulation group
            if(!one.file) pdf(file.path(figdir,paste('Donor', ds.nm[gg],'SubsAcute',cframe$acute.sc[1],'Pts'[show.sdp],'.pdf', sep = '')),
                              w = 6.5, h = 4.5)
            layout(matrix(c(1,4,2,5,3,6,7,7),nr=2,nc=4), w = c(1,1,1,.6))
            par(mar = c(3,1.5,4,.3), oma = c(0,2.7,0,0))
            col <- 'black'
            leg <- ds.nm
            for(styp.i in 1:5) { ## for each substitution type, plot using plot.sdp()
                styp <- styps[styp.i]              
                js <- which(cframe[,styp] == gg & cframe$group.ind!=gg)
                js1 <- which(rowSums(cframe[,styps]!=gg)==0)
                js <- c(js1,js)
                plot.sdp(js = js,  repl = styp, leg = leg,                      
                         main=paste0(LETTERS[1:6], ') ',slabs)[styp.i], ylab = '',
                         yaxt = (styp.i %in% c(1,4)), show.sdp = show.sdp, make.pdf = F, early.yr = 1985, browse=F,
                         cols =NA, col.pl = col, sep.leg = F, show.leg=F, cex.leg = .4, leg.col = 2)
    #            if(styp.i==3) {     ## show X & Y for explained variance calculations
                    par(xpd=NA)
                    sfr.styp <- sframe[js,]
                    ## arrows(2013.5, sfr.styp[sfr.styp[,styp]==3,'sdp.sim.b1'],                               
                    ##        2013.5, sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'],
                    ##        len = .01, angle = 90, code = 3, col = 'black')
                    ## text(x=2013.5, y=sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'], labels='X', col = 'black', pos = 3)
                    ## arrows(2015.5, sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'],                               
                    ##        2015.5, sfr.styp[sfr.styp[,styp]==3,'sdp.sim.s1'],
                    ##        len = .01, angle = 90, code = 3, col = rainbow(length(ds.nm))[3])
                    ## text(x=2015.5, y=sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'], labels='Y', col=rainbow(length(ds.nm))[3], pos = 3)
    #            }
            }
            ##  summary plot ready for all countries to put
            col <- 'black'          # plot color
            ymin <- -1
            ymax <- 1
            ylim <- c(ymin,ymax)
            par(col = col, fg = col, col.axis = col, col.main = col, col.lab=col, bty = 'l', mar = c(3,3,4,1.5))
            frame()
            plot.window(xlim = c(.58,6), ylim = ylim)
            grcol <- rgb(t(col2rgb(gray(.5))), alpha = 250,max=255)
            wher <- 5.5
            segments(.7,1,wher,1, col= grcol, lty = 1, lwd = 3)
            ## text(5.5,1, 'all variance \nexplained', pos=4)
            ## text(5.5,0, 'no variance \nexplained', pos=4)                  
            segments(.7,0,wher,0, col = grcol, lty = 1, lwd = 3)
            boxplot(evar ~ var, data = sfr, bty = 'n', ylim = c(-.25,1), ylab = '', add = T, range = 1.5,
                    las = 2, axes = F, col = 'orange', lty = 1, main = paste0(LETTERS[1:6], ') ',slabs)[6],
            outline = F) # boxplot of explained variance
            axis(1, at = 1:5, LETTERS[1:5], pos = ymin-.1, lwd = 1, las = 1)
            axis(2, at = seq(ymin,ymax, by = .5), labels = NA,mgp = c(3,1,-.2))
            axis(2, at = seq(ymin,ymax, by = 1), las = 2, mgp = c(3,1,-.2))
            mtext(expression(paste('deviance explained 1 - ', (Y^2/X^2))), side = 2, line = 1.5, adj = .5, cex = .7)
            par(xpd=NA)
          text(wher,1, 'all', las = 1, pos = 4) 
           text(wher,.5, 'deviance decreased', las = 2, pos = 4)
           text(wher,-.7, 'deviance increased', las = 2, pos = 4)
          text(wher,0, 'none', las = 1, pos = 4)            
                                        #mtext('serodiscordant proportion', side = 2, outer = T, adj = .5, cex = 1.5)
            par(mar=c(.2,0,0,0))
            plot.new()
            ord <- rev(order(cframe[which(rowSums(cframe[,styps]!=cframe$group.ind)==0),'sdp08']))
            cols <- rainbow(length(ds.nm))
            cols[gg] <- 'black'
            lwds <- rep(1,length(ds.nm))
            lwds[gg] <- 2
            legend('top', ds.nm[ord], col=cols[ord], lwd = lwds[ord], title = 'top to bottom', cex = .85,
                   inset = .1) #, bty = 'n')
            mtext('serodiscordant proportion (SDP)', side = 2.4, outer = T, adj = .5, line = 1)
            if(!one.file) dev.off()
        }                               ## loop through countries
    } ## end country loop
    if(one.file) dev.off()
    ## Summarize all donor-recipient pair relationships by substituted input
    ## only look at relevant columns
}
graphics.off()


####################################################################################################
## Sensitivity analysis of panel (F) to acute phase
drs <- list.files(outdir, full.names = T)[grepl('Acute', list.files(outdir))] # directories with results
drs <- drs[c(1,8,3,7)]                  # 1,7,25,50
col <- 'black'          # plot color
ymin <- -1
ymax <- 1
ylim <- c(ymin,ymax)
pdf(file.path(figdir,'Sub Sens Acute.pdf'), w = 6.5, h = 3)
par(mfrow = c(1,4), col = col, fg = col, col.axis = col, col.main = col, col.lab=col, bty = 'l', mar = c(3,1,2,0), oma = c(0,4,0,7.5))
for(dd in 1:length(drs)) {
  frame()
  plot.window(xlim = c(.58,6), ylim = ylim)
  grcol <- rgb(t(col2rgb(gray(.5))), alpha = 250,max=255)
  wher <- 5.5
  segments(.7,1,wher,1, col= grcol, lty = 1, lwd = 3)
  segments(.7,0,wher,0, col = grcol, lty = 1, lwd = 3)
  load(file = file.path(drs[dd], 'cfs.Rdata'))
  attach(cfs)
  boxplot(evar ~ var, data = sfr, bty = 'n', ylim = c(-.25,1), ylab = '', range = 1.5, add = T,
          las = 2, axes = F, col = 'orange', lty = 1, main = expression(alpha), outline = F) # boxplot of explained variance
  axis(1, at = 1:5, LETTERS[1:5], pos = ymin-.1, lwd = 1, las = 1)
  axis(2, at = seq(ymin,ymax, by = .5), labels = NA,mgp = c(3,1,.2))
  if(dd==1) axis(2, at = seq(ymin,ymax, by = 1), las = 2, mgp = c(3,1,.2))
  detach(cfs)
}
mtext('deviance explained', side = 2, outer = T, adj = .5, cex = 1, line = 2)
par(xpd=NA)
text(wher,1, 'all', las = 1, pos = 4) 
text(wher,.5, 'deviance decreased', las = 2, pos = 4)
text(wher,-.7, 'deviance increased', las = 2, pos = 4)
text(wher,0, 'none', las = 1, pos = 4)            
dev.off()


## ####################################################################################################
## styp.i <- 1
## gg <- 15
## show.sdp <- T
## col.pl <- 'white'
## col <- col.pl
## pdf(file.path(figdir,paste(ds.nm[gg],'SubsAcute-',gg,'-',cframe$acute.sc[1],'Pts'[show.sdp],'.pdf', sep = '')),
##     w = 5, h = 4)
## layout(matrix(c(1,2),nr=1,nc=2), w = c(1,.38))
## styp <- styps[styp.i]              
## js <- which(cframe[,styp] != cframe$group & cframe$group.ind==gg)
## js1 <- which(rowSums(cframe[,styps]!=gg)==0)
## js <- c(js1,js)
## main <- paste0(LETTERS[1:6], ') ',slabs)[styp.i]
## main <- slabs[styp.i]
## cols <- rainbow(length(ds.nm))
## par(fg = col.pl, col.axis = col.pl, col.lab=col.pl, col.main = col.pl, col = col.pl)
## plot.sdp(js = js,  repl = styp, leg = leg,                      
##          main=main, ylab = '', lwds = NA, show.all = T, show.all.pts = T, to.show = NA,
##          yaxt = (styp.i %in% c(1,4)), show.sdp = show.sdp, make.pdf = F, early.yr = 1985, browse=F,
##          cols =NA, col.pl = col, sep.leg = F, show.leg=F, cex.leg = .4, leg.col = 2)
##                                         #            if(styp.i==3) {     ## show X & Y for explained variance calculations
## par(xpd=NA)
## sfr.styp <- sframe[js,]
## ## arrows(2013.5, sfr.styp[sfr.styp[,styp]==3,'sdp.sim.b1'],                               
## ##        2013.5, sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'],
## ##        len = .01, angle = 90, code = 3, col = 'black')
## ## text(x=2013.5, y=sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'], labels='X', col = 'black', pos = 3)
## ## arrows(2015.5, sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'],                               
## ##        2015.5, sfr.styp[sfr.styp[,styp]==3,'sdp.sim.s1'],
## ##        len = .01, angle = 90, code = 3, col = rainbow(length(ds.nm))[3])
## ## text(x=2015.5, y=sfr.styp[sfr.styp[,styp]==3,'sdp.dat.s1'], labels='Y', col=rainbow(length(ds.nm))[3], pos = 3)
## par(mar=c(.2,0,0,0))
## plot.new()
## ord <- rev(order(cframe[which(rowSums(cframe[,styps]!=cframe$group.ind)==0),'sdp08']))
## cols <- rainbow(length(ds.nm))
## cols[gg] <- 'black'
## lwds <- rep(2,length(ds.nm))
## lwds[gg] <- 3
## legend('top', ds.nm[ord], col=cols[ord], lwd = lwds[ord], title = 'top to bottom', cex = .8,
##        inset = .1) #, bty = 'n')
## mtext('serodiscordant proportion (SDP)', side = 2.4, outer = T, adj = .5, line = 1)
## dev.off()

## ####################################################################################################
##            ##  summary plot ready for all countries to put
## pdf(file.path(figdir,paste('SubsAcute-Summary',cframe$acute.sc[1],'.pdf', sep = '')),
##     w = 6.5, h = 4.5)
## col.pl <- 'white'
## col <- col.pl          # plot color
## ymin <- -1
## ymax <- 1
## ylim <- c(ymin,ymax)
## par(col = col, fg = col, col.axis = col, col.main = col, col.lab=col, bty = 'l', mar = c(4,4,4,1.5))
## frame()
## plot.window(xlim = c(.58,6), ylim = ylim)
## grcol <- rgb(t(col2rgb(gray(.5))), alpha = 250,max=255)
## wher <- 5.5
## segments(.7,1,wher,1, col= grcol, lty = 1, lwd = 3)
## ## text(5.5,1, 'all variance \nexplained', pos=4)
## ## text(5.5,0, 'no variance \nexplained', pos=4)
## main <- paste0(LETTERS[1:6], ') ',slabs)[6]
## main <- slabs[6]
## segments(.7,0,wher,0, col = grcol, lty = 1, lwd = 3)
## boxplot(evar ~ var, data = sfr, bty = 'n', ylim = c(-.25,1), ylab = '', add = T, range = 1.5,
##         las = 2, axes = F, col = 'orange', lty = 1, main = '',
##         outline = F) # boxplot of explained variance
## axis(1, at = 1:5,NA, pos = ymin-.1, lwd = 1, las = 1)
## axis(2, at = seq(ymin,ymax, by = .5), labels = NA,mgp = c(3,1,-.2))
## axis(2, at = seq(ymin,ymax, by = 1), las = 2, mgp = c(3,1,-.2))
## mtext('explained deviance', side = 2, line = 2, adj = .5, cex = 1)
## par(xpd=NA)
## ## text(wher,1, 'all', las = 1, pos = 4) 
## ## text(wher,.5, 'deviance decreased', las = 2, pos = 4)
## ## text(wher,-.7, 'deviance increased', las = 2, pos = 4)
## ## text(wher,0, 'none', las = 1, pos = 4)            
## dev.off()
