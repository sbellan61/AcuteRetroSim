## ####################################################################
## create a control file to send to the cluster
## ####################################################################
setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
batch <- 2

batchdirnm <- file.path('results',paste0('abcBatch',batch))
if(!file.exists(batchdirnm))      dir.create(batchdirnm) # create directory if necessary
if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) # setup directory to store Rout files


for(jj in 1:12) { ## separate batches into smaller ones so they start on queue faster
    sink(paste0("abcGo",jj,".txt"))
    ## ####################################################################
    to.do <- (jj-1)*25 + 1:25
    for(ii in to.do) {
        cmd <- paste("R CMD BATCH '--no-restore --no-save --args batch.dir=\"", batchdirnm, "\" seed=", ii, " batch=", batch,
                     "' ABCsimStarter.R ", file.path(batchdirnm, "routs", paste0('abc', ii, ".Rout")), sep='')
        cat(cmd)               # add command
        cat('\n')              # add new line
    }
    sink()
}
