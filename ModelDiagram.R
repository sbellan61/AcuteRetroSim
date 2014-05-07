####################################################################################################
## Make couple timeline model diagram figure
####################################################################################################
rm(list=ls())                           # clear workspace
graphics.off()
set.seed(1)
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/')     # setwd
source("SimulationFunctions.R")
source("ModDiagFxn.R")
load('data files/epic.Rdata')
#fs <- list.files('results/CounterFactual/Acute7/Zambia/', pattern='.Rdata', full.names = T)
#fs <- "results/CounterFactual/Acute7/Zambia/Zambia-110600-1.Rdata"

load(file = file.path('results','RakAcute','blocks.Rdata')) # these are country-acute phase specific blocks
names(blocks)
tp <- which(with(blocks, acute.sc==7 & dur.ac==2 & het.gen.sd==3 & dur.lt==10 & dur.aids==10 & late.sc==5))
tp


fs <- 'results/RakAcute/Uganda/'
fls <- list.files(fs, full.names=T)
head(fls[grepl('96200', fls)],5)
length(fls[grepl('96200', fls)])
fl <- fls[grepl(tp, fls)]
load(fl) ## use this for model diagram figure
                                        #if(output$simj != 1) stop('Not using an "as fitted" output file!')

outdir <- file.path('results','ModelDiagramFigure')
if(!file.exists(outdir)) dir.create(outdir)      # create it
## source('ModelDiagram.R')

## wrapper function for sending all the jobs to different cores for speed. Do this plot for 100
## subgroups & choose the best one for the paper in terms of clearly showing the concepts.
col.sdc <- 'orange'
col.ccp <- 'red'
wrp <- function(run, browse=F, examp) {
  print(paste0('working on run', run))
  pdf(file.path(outdir, paste0('Figure 1-',run,'.pdf')), w = 6.83, h = 10)
  cptime(output, ncpl = 40, inf.only = F, n.inf = 20, seed = run,  ylab2 = 'serodiscordant \nproportion',
         age.pch = NA, inf.pch.pre = NA, inf.pch.extra=NA, inf.pch.within=NA,
         yrmin = 1985, yrmax = 2000, browse = browse,
         col.m = col.sdc,  col.msdc = col.sdc,
         col.f = col.sdc,  col.fsdc = col.sdc,
         examp = examp,
         col.ccp = col.ccp)
  dev.off()
}
 
source("ModDiagFxn.R")
wrp(1, examp = 38)

wrp(56, browse=F)


## nruns <- 100
## mclapply(1:nruns, wrp) ## automatically sends it to all available cores



