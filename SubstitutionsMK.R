####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################
rm(list=ls())                                  # clear workspace
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
load('data files/ds.nm.all.Rdata')        # load country names
load('data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
load(file=file.path('data files/SubAlreadyDone.Rdata'))
## source('SubstitutionsMK.R')

####################################################################################################
####################################################################################################
##   SUBSTITUTION ANALYSIS
####################################################################################################
## To investigate possible drivers of country-to-country variation in serodiscordant proportions we
## substituted population-level HIV prevalence, couple-formation patterns or HIV transmission rates
## estimated from a “donor” country into a simulation model fit for a “recipient” country, and
## measured the extent to which a substitution shifted the serodiscordant proportion from that of
## the recipient country to that of the donor country.
####################################################################################################
countries <- c(1:length(ds.nm)) ## countries to do
ncount <- length(ds.nm)
# acute phase relative hazards to do (must match those fitted previously, since it determines the
# transmission coefficients used)
acutes <- c(1,5,7,10,25,30,40,50) 
each.val <- 200 #  number of couples per couple formation (marital) cohort
counterf.betas <- F                       # change betas in counterfactuals? if not change beta_within & c's (so beta_within affects all routes)
sub.betas <- F                           # substitute betas? if not beta_within & c's
nrb <- 0                     # initialize blocks data frame (rows = 0)
totn <- 0                    # total # jobs (starts at 0)
num.doing <- 0
nn <- 200 # number of simulations per country-acute combination (must be bigger than max(sel) later but is fine to leave big
substitute <- TRUE                      # doing a substitution analysis?
#################################################################################################### 
sink("Ac1SubstitutionsControlFile.txt")                      # create a control file to send to the cluster
outdir <- file.path('results','SubstitutionAnalysis') # set up directory for output
if(!file.exists(outdir))      dir.create(outdir) # create directory if necessary
for(aa in acutes)  {          # loop through acute phase relative hazard
    acdirnm <- file.path(outdir,paste0('Acute', aa)) # set up directory for each acute phase relative hazard
    if(!file.exists(acdirnm))      dir.create(acdirnm) # create directory if necessary
    for(cc in countries) {
        batchdirnm <- file.path(acdirnm, ds.nm[cc]) # setup directory for country
        if(!file.exists(batchdirnm))      dir.create(batchdirnm) # created if necessary
        if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) # setup directory to store Rout files
           ######################################################################
           ######################################################################           
           ## Set defaults for all parameters for each simulation, simulatin specific-values set later
           ######################################################################                      
           group <- rep(cc,nn)          # country group.
           ## set substitution country (donor country) indices to country by default & change later:
           ## s.epic, s.demog, s.bmb...
           for(svar in c('epic','demog', hazs)) assign(paste0('s.',svar), rep(cc,nn))
           ## set all phases to have same infectivity
           for(ph in c('acute','late', 'aids')) assign(paste0(ph,'.sc'), rep(1,nn))
           death <- rep(T,nn)       # include death
           ## all haz scalars (bmb.sc,etc...) set to 1
           for(hh in hazs) assign(paste0(hh,'.sc'), rep(1,nn))
           ## set heterogeneity defaults (i.e. het.b, het.b.sd, het.b.cor,etc...)
           for(ht in c('b','e','p','gen','beh')) {
             assign(paste0('het.',ht), rep(F,nn))
             assign(paste0('het.',ht,'.sd'), rep(0,nn))
             assign(paste0('het.',ht,'.cor'), rep(0,nn))
           }
           scale.by.sd <- rep(T,nn)     # scale by standard deviation?
           scale.adj <- rep(1,nn)       # arbitrary scalar if not doing that.
           infl.fac <- rep(200,nn)  # inflation factor for non-parametric couple pseudo-population builder
           maxN <- rep(10^5,nn)     # max pseudopopulation size
           sample.tmar <- rep(F,nn) # sample marital (couple formation) date from copulas?
           psNonPar <- rep(F,nn) #  use non-parametric couple pseudo-population builder?
           each <- rep(each.val, nn) # how many couples per marital (couple formation) cohort
           ######################################################################
           ######################################################################           
           ## ***set acute phase for within this aa-loop***
           acute.sc <- rep(aa,nn)
           ######################################################################
           ## Keep track of simulations we're doing in a data.frame called blocks
           if(cc==1 & aa==1) {          # initialize blocks for first aa-loop iteration
             blocks <- data.frame(start = 1, end = 1, acute = aa, cc = cc, country = ds.nm[cc], lab = 'as fitted', sub.country = NA)
           }else{                       # add to blocks 
             if(!(aa==1 & cc==1))       nrb <- nrow(blocks) # how long is blocks (so we know how to add appropriately
             blocks <- rbind(blocks, data.frame(start = 1+nrb, end = 1+nrb, acute = aa, cc = cc, country = ds.nm[cc], lab = 'as fitted', sub.country = NA))
           }
           ######################################################################
           ######################################################################
           ## Substitutions
           ######################################################################
           ## Epidemic curves (different for WA b/c that country-group has multiple epidemic curves)
           if(ds.nm[cc]!='WA') { # if not WA substitute all other countries epidemic curves in
             sel <- 2:(ncount-1)                 # selection index
             s.epic[sel] <- c(1:ncount)[-c(which(ds.nm=='WA'), cc)]
             blocks <- rbind(blocks, data.frame(start = min(sel)+nrb, end = max(sel)+nrb, acute = acute.sc[sel], cc = cc,
                                                country = ds.nm[cc], lab = 'epidemic curve',
                                                sub.country = c(1:ncount)[-c(which(ds.nm=='WA'), cc)]))
           }else{
             sel <- 2:ncount                    # selection index
             s.epic[sel] <- c(1:ncount)[-c(cc)] # if WA, do DRC at end intsead for WA
             blocks <- rbind(blocks, data.frame(start = min(sel)+nrb, end = max(sel)+nrb, acute = acute.sc[sel], cc = cc,
                                                country = ds.nm[cc], lab = 'epidemic curve', sub.country = c(1:ncount)[-cc]))
           }    
           ## Pre-couple transmission
           sel <- (sel[length(sel)] + 1):(sel[length(sel)] + ncount - 1) # update selection index, next 9 simulations
           s.bmb[sel] <- c(1:ncount)[-cc]
           s.bfb[sel] <- c(1:ncount)[-cc]
           blocks <- rbind(blocks, data.frame(start = min(sel)+nrb, end = max(sel)+nrb, acute = acute.sc[sel], cc = cc,
                                              country = ds.nm[cc], lab = 'pre-couple', sub.country = c(1:ncount)[-cc]))
           ## Extra-couple transmission
           sel <- (sel[length(sel)] + 1):(sel[length(sel)] + ncount - 1) # update selection index, next 9 simulations
           s.bme[sel] <- c(1:ncount)[-cc]
           s.bfe[sel] <- c(1:ncount)[-cc]
           blocks <- rbind(blocks, data.frame(start = min(sel)+nrb, end = max(sel)+nrb, acute = acute.sc[sel], cc = cc,
                                              country = ds.nm[cc], lab = 'extra-couple', sub.country = c(1:ncount)[-cc]))
           ## Within-couple transmission
           sel <- (sel[length(sel)] + 1):(sel[length(sel)] + ncount - 1) # update selection index, next 9 simulations
           s.bmp[sel] <- c(1:ncount)[-cc]
           s.bfp[sel] <- c(1:ncount)[-cc]
           blocks <- rbind(blocks, data.frame(start = min(sel)+nrb, end = max(sel)+nrb, acute = acute.sc[sel], cc = cc,
                                              country = ds.nm[cc], lab = 'within-couple', sub.country = c(1:ncount)[-cc]))
           ## Relationship patterns (from copulas)
           sel <- (sel[length(sel)] + 1):(sel[length(sel)] + ncount - 1) # update selection index, next 9 simulations
           s.demog[sel] <- c(1:ncount)[-cc]
           blocks <- rbind(blocks, data.frame(start = min(sel)+nrb, end = max(sel)+nrb, acute = acute.sc[sel], cc = cc,
                                              country = ds.nm[cc], lab = 'relationship patterns', sub.country = c(1:ncount)[-cc]))
           ######################################################################
           ## Create control file text to send to cluster.
           for(ii in 1:max(sel)) {
             jb <- ii                   # job num
             totn <- totn+1             # total jobs
             cmd <- paste("R CMD BATCH '--args jobnum=", totn, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
                          " group.ind=", group[ii], " substitute=", substitute, " s.epic=", s.epic[ii],  " s.demog=", s.demog[ii],
                          " sub.betas=", sub.betas, " counterf.betas=", counterf.betas,
                          " s.bmb=", s.bmb[ii], " s.bfb=", s.bfb[ii],
                          " s.bme=", s.bme[ii], " s.bfe=", s.bfe[ii],
                          " s.bmp=", s.bmp[ii], " s.bfp=", s.bfp[ii], 
                          " death=", death[ii],
                          " acute.sc=", acute.sc[ii], " late.sc=", late.sc[ii]," aids.sc=", aids.sc[ii],
                          " bmb.sc=", bmb.sc[ii], " bfb.sc=", bfb.sc[ii],
                          " bme.sc=", bme.sc[ii], " bfe.sc=", bfe.sc[ii],
                          " bmp.sc=", bmp.sc[ii], " bfp.sc=", bfp.sc[ii],
                          " het.b=", het.b[ii], " het.b.sd=", het.b.sd[ii], " het.b.cor=", het.b.cor[ii],
                          " het.e=", het.e[ii], " het.e.sd=", het.e.sd[ii], " het.e.cor=", het.e.cor[ii],
                          " het.p=", het.p[ii], " het.p.sd=", het.p.sd[ii], " het.p.cor=", het.p.cor[ii],                     
                          " het.gen=", het.gen[ii], " het.gen.sd=", het.gen.sd[ii], " het.gen.cor=", het.gen.cor[ii],
                          " het.beh=", het.beh[ii], " het.beh.sd=", het.beh.sd[ii], " het.beh.cor=", het.beh.cor[ii],
                          " scale.by.sd=", scale.by.sd[ii], " scale.adj=", scale.adj[ii],
                          " infl.fac=", infl.fac[ii], " maxN=", maxN[ii], " sample.tmar=", sample.tmar[ii],
                          " psNonPar=", psNonPar[ii], " seed=1 tmar=(65*12):(113*12) each=", each[ii],
                          " tint=113*12' SimulationStarter.R ", file.path(batchdirnm, "routs", paste0(ds.nm[group[ii]], jb, ".Rout")), sep='')
             if(ii > 0 & acute.sc[ii]==1 & !totn %in% already.done ) {
               num.doing <- num.doing + 1
                 cat(cmd)               # add command
                 cat('\n')              # add new line
               } } } }
sink()
blocks$sub.country.nm <- ds.nm[blocks$sub.country]
head(blocks,50)
tail(blocks,50)
print(totn)
print(num.doing)
save(blocks, file=file.path(outdir,'blocks.Rdata'))
####################################################################################################
