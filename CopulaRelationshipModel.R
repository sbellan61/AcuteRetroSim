####################################################################################################
## Assess multivariate correlation between male & female ages at sexual debut, durations of sexual
## activities before marriage & their couple formation date for all data groups analyzed. Then use
## multivariate normal copulas to summarize multivariate correlations and simulate representative
## couples from each country, i.e. a pseudo-couple population.
####################################################################################################
## Steve Bellan, 2013
####################################################################################################
## Note that age at first intercourse is given in years in the DHS
## (though those that were at the time of marriage may be in months
## since marital date was given as a month).
####################################################################################################
rm(list=ls())                           # clear workspace
library(graphics); library(animation); library(abind); library(mvtnorm); library(copula); library(logspline); library(actuar)
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
source("SimulationFunctions.R")                   # load simulation/fitting functions
load("data files/epic.Rdata")                     ##  epidemic curves
load("data files/pars.arr.ac.Rdata") # fitted transmission coefficients (out.arr) across range of acute phase RH's (in.arr)
load("data files/csurv.Rdata")       # HIV survival times
load("data files/ds.nm.all.Rdata")     # country names
load("data files/allDHSAIS.Rdata")         # DHS data
outdir <- file.path('results','Copulas')
if(!file.exists(outdir))        dir.create(outdir)
## dat <- dat[,names(aisdat)]          # only use columns in AIS too
## dat <- rbind(dat, aisdat)           # add data sets
odat <- dat ## original data doesn't change ever
odat$bd <- apply(cbind(odat$tmar-odat$tms,odat$tmar-odat$tfs), 1, max) # before sex duration (max across partners)
odat$cd <- odat$tint - odat$tmar                                       # couple durations
## weight by survival probabilities to account for underrepresentation of couples who are likely to
## have died from HIV by the time of the DHS interview:
weight <- TRUE
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") # names of transmission coefficients (betas) = prevalence-standardized hazards
vars <- c("ams","afs","mdur", "fdur","tmar") # m/f sexual debuts; m/f duration of sex before marriage; marital date
nvar <- length(vars)                         # 5
lwd <- 1.2                                   # line width for plots
hist.col <- "dark gray"                      # histogram colors
fit.col <- "red"                             # fit distribution colors
labs <- c("M sexual debut", "F sexual debut", # label name
          "M duration of sex before marriage", "F duration of sex before marriage", "couple formation date")
start <- 0                              # start age (months)
end <- 60*12                            # end age (months); 60 is oldest age of anyone in DHS
## ds.nm <- levels(dat$group)
nds <- length(ds.nm)                         # country-groups
## source('CopulaRelationshipModel.R')

######################################################################
## The below chunk of code proceeds as follows: for each country,
## (1) for each variable fit a continuous kernel density estimator (using density or logspline) to
##      get marginal probability density functions (PDF) for each of 5 relationship variables
## (2) convert continuous PDFs to discrete probability distribution function (PDF) for the
##      probability that each variable takes a specific value (in monthly units)
## (3) plot histogram of variable, plot fitted PDF over it to assess fit.
######################################################################
## Set up data.frames for each discrete distribution (cols are each country, rows are months (mm),
## values are probability that that the variable (ams, etc...) is mm months for that country)
maxm <- c(rep(60,4),120)*12 ## first 4 variables go up to 60 yrs, marital date up to 120 (yr 2020)
for(vv in 1:nvar) {
    temp <- data.frame(matrix(NA, nr = maxm[vv] + 1, nc = 1 + nds))
    colnames(temp) <- c('mm', ds.nm)        # mm is months; rest of columns are country names
    temp$mm <- 0:maxm[vv]
    assign(paste0('p',vars[vv]), temp) # set pams, pafs,... to temp
}
tail(pams,5); tail(ptmar,5)
## set up pdf to save fits
file.nm <- ifelse(weight, "weighted marginal fits.pdf", "marginal fits.pdf") # pdf names
pdf(file.path(outdir,file.nm), w = 14, h = 10)
for(ii in 1:length(ds.nm)) { ## for each country
    print(paste("working on", ds.nm[ii]))
    group.ind <- ii
    dat <- odat[odat$group == ds.nm[group.ind],]
    dat <- add.vars(dat)
    par(mfrow=c(3,2))
    if(weight) { ## let's weight couples survival based on transmission parameters fit with
        ## assumption that acute RH = 7
        pars <- out.arr[hazs,2, which(in.arr[,ii,2]==7), group.ind] # median transmission coefficients for that country
        ## weight couples by infl fator = 1 / survival probabilities. Survival probabilities are
        ## calculated by fitted model's estimate that each observed couple observed to survive
        ## (Markov Chain state probability model)
        sim <- pcalc(pars, dat, sim = T, survive = T, browse = F, trace = F)
        surp <- rowSums(sim$pser.a)
        infl <- 1/surp
        infl1 <- infl/sum(infl)          # normalize for density()
    }
    for(jj in 1:5) { ## for each relationship variable
        temp <- dat[,vars][,jj]
        xmax <- max(temp)
        logd <- log(temp+1)             # fitting all variables on a log scale since they're all >=0
        if(jj<3) { # use normal kernel density estimator to smooth distribution of sexual debuts
            breaks <- seq(0,xmax+12, by = 12)
            bw <- 10
            if(weight) {
                den <- density(temp, from = start, to = xmax, bw = bw, weights = infl1)
            }else{
                den <- density(temp, from = start, to = xmax, bw = bw)
            }
            fden <- approxfun(den$x, den$y, yleft = 0, yright = 0) # create a function that interpolates the densities
            h.real <- hist(temp, breaks = breaks, plot = T, xlim = c(0,xmax), col = hist.col, # histogram of results
                           freq = F, main = labs[jj], xlab = "age in months", border = NA)
            if(vars[jj] == "ams") {     # smooth right edge, but truncate left edge
                ## no one was < 5 years old when they started having sex, use discret() to
                ## discretize a continuous density function and then assign to pams
                pams[,ii+1] <- discret(fden, leftz = 5*12, redge.smooth = T, vals = temp)$pp 
                lines(pams[,c(1,ii+1)], col = fit.col, lwd = lwd) # add line to histogram showing fit
            }
            if(vars[jj] == "afs") { ## same for afs
                pafs[,ii+1] <- discret(fden, leftz = 5*12, redge.smooth = T, vals = temp)$pp
                lines(pafs[,c(1,ii+1)], col = fit.col, lwd = lwd)
            }            
        }else{                        
            if(jj<5) { ## for durations of sex before marriage do logspline kernel density estimator
                ## (seems to deal better with numbers closer to 0)
                if(weight) { ## to weight we need to inverse-prob-of-survival inflate observed
                    ## couples since logspline() does not have a weight argument
                    ## numb: poisson inflation size, doesn't need to be infl1 since absolute # unimportant
                    numb <- rpois(length(logd), lambda = 50*infl) 
                    pslogd <- rep(logd, numb)                     # inflated log-ed variable
                    den <- logspline(pslogd, lbound = 0, ubound = max(pslogd +.1), # knots chosen empirically for best fit
                                     knots = c(log(c(1,2,5,9, 25, 80, 150))), maxknots = 7)
                }else{ ## no need to inflate since not weighting
                    den <- logspline(logd, lbound = 0, ubound = max(pslogd +.1),
                                     knots = c(log(c(1,2,5,9, 25, 80, 150))))
                }
                h.real <- hist(temp, breaks = c(0:1000), plot = F) # get y-max of real data
                h.fit <- hist(floor(exp(rlogspline(10000,den))), breaks = c(0:100000), plot = F) # get y-max of fitted spline
                ## plot histogram of real data
                h.real <- hist(temp, breaks = c(0:1000), plot = T, xlim = c(0,200), ylim = c(0,max(h.real$dens, h.fit$dens)),
                               main = labs[jj], xlab = "months", ylab = "Density", freq = F, col = hist.col, border = NA)
                if(vars[jj] == "mdur") { ## logspline gives CDF, transform to PDF
                    temp <- plogspline(log(1:(12*60+2)), den) # cdf
                    pmdur[,ii+1] <- temp[2:length(temp)] - temp[1:(length(temp)-1)] # pdf
                    lines(pmdur[,c(1,ii+1)], col = fit.col, lwd = lwd) # add fitted density function to plot
                }
                if(vars[jj] == "fdur") { ## same for female duration of sex before couple formation
                    temp <- plogspline(log(1:(12*60+2)), den) # cdf
                    pfdur[,ii+1] <- temp[2:length(temp)] - temp[1:(length(temp)-1)] #pdf
                    lines(pfdur[,c(1,ii+1)], col = fit.col, lwd = lwd)
                }            
            }else{ ## if jj==5, vars[vv]=tmar; get distribution of couple formation dates, just use density
                bw <- 25 ## bandwith for density (empirically determined)
                if(weight) {
                    den <- density(temp, from = start, to = 1.3*xmax, bw = bw, weights = infl1)
                }else{
                    den <- density(temp, from = start, to = 1.3*xmax, bw = bw)
                }
                fden <- approxfun(den$x, den$y, yleft = 0, yright = 0) # interpolate densities
                ptmar[,ii+1] <- discret(fden, range = c(0,120*12))$pp  # discretize into monthly intervals from 1900-2020
                h.real <- hist(temp, breaks = 30, plot = T, xlim = c(0,1500), col = hist.col,
                               freq = F, main = labs[jj], xlab = "months since 1900", border = NA)
                lines(den, col = fit.col, lwd = lwd) # add fitted distribution
            }
        }
    }
    plot(0,0, type = "n", axes = F, bty = "n", xlab = "", ylab = "") # add legend to last panel
    legend("topleft",c("data","fitted distribution"), col = c(hist.col,fit.col), pch = 19, bty = "n", cex = 2)
    mtext(ds.nm[ii], side = 3, outer = T, line = -3, cex = 2) # write country name at top
}
dev.off()
####################################################################################################

####################################################################################################
## transform PDF's to cumulative density functions's (CDF's) since we need these for the copula
## model
cnames <- paste0("c",vars)              # names of cfs
pnames <- paste0("p",vars)              # names of pdfs
for(jj in 1:5)  assign(cnames[jj], get(pnames[jj])) # copy data frame structure
for(ii in 1:length(ds.nm)) { ## for each country
    print(paste("working on", ds.nm[ii]))
    for(jj in 1:5) {
        tempp <- get(pnames[jj]) ## temp PDF
        for(rr in 1:nrow(tempp)) { ## for each row, sum up cumulative probability; inefficient loop but only doing once
            tempc <- get(cnames[jj])    # temp CDF
            tempc[rr, ii+1] <- sum(tempp[tempp$mm <= tempc$mm[rr], ii + 1]) # cumulative probability
            assign(cnames[jj], tempc)                                       # set CDF
        }
    }
}
####################################################################################################

####################################################################################################
## Find correlation between quantiles of the 5 variables from spline-fitted CDFs created above
for(ii in 1:length(ds.nm)) { ## for each country
    print(paste("working on", ds.nm[ii]))
    group.ind <- ii
    file.nm <- file.path(outdir, paste(ds.nm[group.ind],"couple parameter correlations")) # create file name
    dat <- odat[odat$group == ds.nm[group.ind],] ## select data set
    dat <- add.vars(dat) ## add transformed variables (ams,afs,mdur,fdur,tmar)
    cdat <- dat[,vars]   # cumulative probability behind value, just to get size of data frame to fill ni
    cdat[,] <- NA        ## remove all values
    for(jj in 1:5) {   ## transform data to empirical quantiles on [0,1] using CDF's created above
        temp <- dat[,vars][,jj]
        tcdf <- get(cnames[jj])
        cdat[,jj] <- tcdf[match(temp,tcdf$mm), ii + 1]
    }
    rdat <- apply(cdat, 2, qnorm) ## transform data to (-Inf, Inf) using normal quantile function
    colmns <- abs(apply(rdat,2,mean)) ## check col means ~= 0 (i.e. that mean quantile is ~50%)
    if(sum(colmns > .25) > 1) {
        print("error with normal quantile values: distr not centered around 0")
        print(colmns)
    }    
    if(weight) {## find covariance (need to weight based on survival later on)
        ## need to calculate the survival probability for the
        hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")
        ## weight couples by infl fator = 1 / survival probabilities (again using acute RH=7)
        pars <- out.arr[hazs,2, which(in.arr[,ii,2]==7), group.ind] # median transmission coefficients for that country
        sim <- pcalc(pars, dat, sim = T, survive = T, browse = F, trace = F)
        surp <- rowSums(sim$pser.a)
        infl <- 1/surp
        infl1 <- infl/sum(infl)         # normalize for density()
        tempcwt <- cov.wt(rdat, wt = infl1) # get weighted covariance-variance matrix
    }else{ ## get covariance-variance matrix without weights
        tempcwt <- cov.wt(rdat)
    }
    if(ii==1) { ## initialize sigmas
        sigmas <- tempcwt$cov
    }else{ ## append to sigmas into an array
        sigmas <- abind(sigmas, tempcwt$cov, along = 3)
    }
}
cop.sigmas <- sigmas ## rename so as not to conflict with 'sigmas' corresponding to MCMC multivariate block sampling
save(cop.sigmas, cnames, pnames, cams, cafs, cmdur, cfdur, ctmar, ## save everything into data files directory
     pams, pafs, pmdur, pfdur, ptmar, vars, labs, file = file.path('data files',"copula sigmas.Rdata"))
save.image(file = file.path(outdir,"copula workspace.Rdata"))
####################################################################################################

####################################################################################################
## Plot fitted copula models along with DHS data: do it as a function and send countries to
## different cores for speed
####################################################################################################
cophist <- function(country) {
    ## get axes ranges to make them the same for both plots
    breaks <- list(seq(0,60*12, by = 12), seq(0,60*12, by = 12), seq(0,60*12, by = 1), seq(0,60*12, by = 1),
                   seq(600, 1600, by = 12))
    lbreaks <- lapply(breaks, function(x) log(x+1)) # log breaks
    ylim <- data.frame(l = rep(0,5), u = c(.02, .04, .25, .3, .01)) # lower & upper y limits
    temp <- add.vars(odat)                                          ## original data
    rgs <- apply(temp[,vars], 2, range)                             ## range of all variables
    lrgs <- log(rgs+1)                                              ## on a log scale
    typs <- c(' real', ' copula')
    transfs <- c('untransformed', 'log-transformed')
    print(paste("working on", ds.nm[country]))
    group.ind <- country
    dat <- add.vars(odat[odat$group==ds.nm[group.ind],]) ## real data
    s.epic.nm <- NA                                      # to let rcop() know we're not substituting
    psdat <- add.vars(rcop(country, N = 10^5, tint = 2011)) ## call rcop to simulate pseudopopulation of 10^5 couples
    for(i.typ in 1:2) {                                     ## for real & fit
        for(i.transf in 1:2) { ## for untransformed & log-transformed
            typ <- typs[i.typ]
            transf <- transfs[i.transf]
            if(!file.exists(file.path(outdir,transf))) dir.create(file.path(outdir,transf)) # create directory
            jpeg(file.path(outdir,transf, paste0(ds.nm[country], typ, ' ',transf, ".jpg")), w = 1200, h = 1200)
            if(i.typ==1) tempd <- dat else tempd <- psdat
            sbpairs(tempd[,vars], browse = F, do.pdf = F, do.jpeg = F, rgs = rgs, hist.freq = F, yrg = ylim, breaks = breaks)
            mtext(paste(typ, transf), side = 3, line = 0, outer = T, cex = 2)
            dev.off()
        }
    }
    print("difference between variable means (real - copula) =") ## shouldn't be more than a few months
    print(apply(dat[,vars], 2, mean) - apply(psdat[,vars], 2, mean))
}
mclapply(1:length(ds.nm), cophist) ## send country-group jobs out to cores

######################################################################
## Explore truncation problem: want to see if the secular trends are
## messed up by the fact that we are fitting a continuous PDF to tmar
## even though the data are all necessarily truncated by the interview
## date. Let's compare secular trends & distributions between
## pseudopopulations of couples from sampling based on marital cohors
## (10 couples form per month for months 784:1336), based on sampling
## marital dates from the fitted marginal curve, and from the real
## data)
######################################################################
if(!file.exists(file.path(outdir, 'secular trends'))) dir.create(file.path(outdir, 'secular trends'))
sz <- 600
sec.fxn <- function(country) {    
    dtype <- c('synthetic cohort copula', 'full sample copula', 'real data')
    ## with & without sampling tmar
    test1 <- add.vars(rcop(country, sample.tmar = F, NN =10^4, tmar=784:1336, each = 10, tint = 1337, browse = F))
    test2 <- add.vars(rcop(country, sample.tmar = T, NN =10^4, tmar=784:1336, each = 10, tint = 1337, browse = F))
    real <- add.vars(odat[odat$group == ds.nm[country],])
    rgs <- apply(rbind(test1[,vars],test2[,vars],real[,vars]), 2 , range) ## get rgs for axis
    tempds <- list(test1, test2, real)
    for(dt.i in 1:3) {
        tempd <- tempds[[dt.i]]
        jpeg(file.path(outdir, 'secular trends', paste0(ds.nm[country], ' ', dtype[dt.i],".jpg")), width = sz, height = sz)
        par(mfrow=c(2,2))
        for(ii in 1:4) { ## for ams,afs,mdur,fdur
            scatter.smooth(tempd$tmar, tempd[,vars[ii]], xlab = "couple formation date", ylab = labs[ii],
                           bty = "n", pch = 19, cex = .5, ylim = rgs[,ii])
            abline(lm(tempd[,vars[ii]]~tempd$tmar), col = "red", lwd = 3)
        }
        dev.off()
    }
}
mclapply(1:length(ds.nm), sec.fxn) ## send country-group jobs out to
## cores Results look OK. But we'll still do a sensitivity analysis of
## all results to marital cohorts vs marital date sampling.
