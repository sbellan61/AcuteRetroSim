#################################################################################################### 
## Showing that person-time calculations depend on hazard and partner serostatus
figdir <- 'FiguresAndTables'

n <- 10^5 ## number of couples
intervalLength <- 10 ## how many months long is the interval
hazs.inc <- .007 * c(1:50) ## the rate at which they infect their partner (example is chronic phase hazard in Rakai * RH[acute]=1,5,7 or 25)
pt.df.inc <- data.frame(all = rep(NA, length(hazs.inc)), serodiscordant = NA, concordantPos = NA) ## intialize person-time data frame
for(ii in 1:length(hazs.inc)) {
    haz <- hazs.inc[ii]
    timeIndexInf <- runif(n, 0,intervalLength) ## index partner's infection time is uniformly randomly distributed within an interval unconditional on any other information
    timeTilPartnerInf <- rexp(n, rate = haz) ## how long after index partner infection do they infect their partner?
    bothInfSameInterval <- (timeIndexInf + timeTilPartnerInf) < intervalLength ## was the second partner infected in the same interval?
    endOfPersonTimeAtRisk <- sapply(timeIndexInf + timeTilPartnerInf, function(x) min(intervalLength,x)) ## time of 2nd partner infection or end of interval
    personTimeAtRiskSecondPartner <-  endOfPersonTimeAtRisk - timeIndexInf 
    pt.df.inc$all[ii] <- mean(personTimeAtRiskSecondPartner)
    pt.df.inc$serodiscordant[ii] <- mean(personTimeAtRiskSecondPartner[!bothInfSameInterval])
    pt.df.inc$concordantPos[ii] <- mean(personTimeAtRiskSecondPartner[bothInfSameInterval])
}
pt.df.inc


#################################################################################################### 
## Non-incident intervals (prevalent couples)

n <- 10^5 ## number of couples
intervalLength <- 10 ## how many months long is the interval
hazs.prev <- .007 * c(.2,.5,1,2,5) ## the rate at which they infect their partner (example is chronic phase hazard in Rakai * 1/5 through 5
pt.df.prev <- data.frame(all = rep(NA, length(hazs.prev)), serodiscordant = NA, concordantPos = NA) ## intialize person-time data frame
for(ii in 1:length(hazs.prev)) {
    haz <- hazs.prev[ii]
    timeTilPartnerInf <- rexp(n, rate = haz) ## how long after index partner infection do they infect their partner?
    infected <- timeTilPartnerInf < intervalLength ## was the second partner infected in the same interval?
    personTimeAtRiskSecondPartner <- sapply(timeTilPartnerInf, function(x) min(intervalLength,x)) ## time of 2nd partner infection or end of interval
    pt.df.prev$all[ii] <- mean(personTimeAtRiskSecondPartner)
    pt.df.prev$serodiscordant[ii] <- mean(personTimeAtRiskSecondPartner[!infected])
    pt.df.prev$concordantPos[ii] <- mean(personTimeAtRiskSecondPartner[infected])
}
pt.df.prev


pdf(file.path(figdir, 'Person-Time Conditionality.pdf'), w = 7, h = 3)
par('ps'=8, lwd = 2, mar = c(4,5,1,.5), mfrow = c(1,2))
with(pt.df.prev, {
    plot(hazs.prev, all, xlab = 'hazard (per month)', ylab = 'mean person-time at risk (months)', type = 'l', xlim = c(0,max(hazs.prev)), ylim = c(0,10), bty = 'n',
         main = 'Prevalent Couples')
    lines(hazs.prev, serodiscordant, col = 'red')
    lines(hazs.prev, concordantPos, col = 'orange')
})
abline(h=5, lty = 2, col = 'orange')
abline(h=10, lty = 2, col = 'red')
par(xpd=NA)
arrows(.007, -2.5, .007, -.5, len = .05, angle = 30, code = 2)
text(.007, -2.5, 'Rakai \nestimate', pos = 1)
par(xpd=F)
with(pt.df.inc, {
    plot(hazs.inc, all, xlab = 'hazard (in multiples of chronic phase \nhazard 0.007 per month)', ylab = 'mean person-time at risk (months)', type = 'l', xlim = c(0,max(hazs.inc)), ylim = c(0,10), bty = 'n', xaxt = 'n', main = 'Incident Couples')
    axis(1, at = pretty(hazs.inc/.007,5)*.007, label = pretty(hazs.inc/.007,5))
    lines(hazs.inc, serodiscordant, col = 'red')
    lines(hazs.inc, concordantPos, col = 'orange')
})
abline(h=2.5, lty = 2, col = 'orange')
abline(h=5, lty = 2, col = 'red')
legend('topleft', leg = c('all', 'serodiscordant', 'concordant positive', 'Wawer assumption'), title = 'by serostatus at end of interval', col = c('black','red','orange','black'), lwd=2, lty = c(1,1,1,2), bty = 'n', cex = .8)
graphics.off()
