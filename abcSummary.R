sapply(c("SimulationFunctions.R","RakFunctions.R",'abcFunctions.R'), source) # load Rakai analysis simulation functions from script
fls <- list.files('results/abcBatch1/', pattern='Rdata', full.names=T)
rcohsListAll <- list()
for(ff in fls) {
    load(ff)
    rcohsListAll <- c(rcohsListAll, rcohsList)
}
length(rcohsListAll)
head(rcohsListAll[[1]]$rakll)

wtabSims <- mclapply(rcohsListAll, function(x) sbmod.to.wdat(x$rakll, excl.by.err = T, browse=F, giveLate=F, condRakai=T, giveProp=T),
                     mc.cores=12)

contTabsSim <- mclapply(wtabSims, function(x) {
    x <- within(x, {
        inct$ni <- with(inct, n-i)
        inct <- inct[,c('i','ni')]
        prevt$ni <- with(prevt, n-i)
        prevt <- prevt[,c('i','ni')]})}, mc.cores=12)

gS <- mclapply(contTabsSim, gSumStat, mc.cores = 12)
gVals <- unlist(lapply(gS, '[',1))

parmsMat <- matrix(unlist(lapply(rcohsList, '[[', 'pars')), nr = 59, nc = 9, byrow=T)
colnames(parmsMat) <- names(rcohsList[[1]]$pars)
parmsMat <- data.frame(parmsMat, gVals, job = 1:nrow(parmsMat))
pmat <- parmsMat[order(gVals),c('job','acute.sc','dur.ac','het.gen.sd','gVals')]

head(pmat)


gS[[46]]


tst <- contTabsSim[[1]]
gSumStat(tst)

x <- matrix(c(10,13,8,15),2,2)
chisq.test(x+.5)

chisq.test(matrix(c(10,13,8,15),2,2))$statistic
chisq.test(matrix(c(2,11,1,8),2,2))$statistic
