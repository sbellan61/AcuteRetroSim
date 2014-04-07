library(mnormt);library(coda); library(stats4);library(plyr)#library(lme4)
rm(list=ls(all=TRUE))
setwd('/home1/02413/sbellan/Rakai/SDPSimulations/') # setwd
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
source('RakFunctions.R')                   # load simulation functions
load("data files/allDHSAIS.Rdata")         # DHS data
load("data files/pars.arr.ac.Rdata")
load("data files/ds.nm.all.Rdata")
load("data files/epic.Rdata")     # infectious HIV prevalence
to.plot <- T

head(dat)
