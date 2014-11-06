## ####################################################################
## create a control file to send to the cluster
## ####################################################################
setwd('/home1/02413/sbellan/Rakai/AcuteRetroSim/')
batchdirnm <- file.path('results','abcBatch1')
if(!file.exists(batchdirnm))      dir.create(batchdirnm) # create directory if necessary
if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) # setup directory to store Rout files

sink("abcGo.txt") 
## ####################################################################
to.do <- 1:50
for(ii in to.do) {
  cmd <- paste("R CMD BATCH '--no-restore --no-save --args out.dir=\"", batchdirnm, "\" seed=", ii, "' ABCsimStarter.R ", file.path(batchdirnm, "routs", paste0('abc', ii, ".Rout")), sep='')
  cat(cmd)               # add command
  cat('\n')              # add new line
}
sink()
