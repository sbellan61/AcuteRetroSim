library(plyr)
rm(list=ls(all=T))
####################################################################################################
## Jacquez et al. 1994

## One-group (homogenous) model
beta1 <- c(.08,.1,.12,.14)
beta24 <- .001           # p.9 (results of sim: homogenous populations)
beta5 <- .005
beta67 <- .01
dur.ac <- 2
rh.ac <- beta1/beta24
ehm.ac <- (rh.ac-1)*dur.ac
r1 <- c(.8,1,1.2,1.4) ## Table 4
r24 <- rep(.33,4)
r5 <- rep(.81,4)
r67 <- rep(1.25,4)
af.ac <- signif(r1 / (r1+r24+r5+r67) * 100,3)

onegr <- data.frame(rh.ac, dur.ac, ehm.ac, af.ac, data = 'epidemic curve', type = 'homogenous', estim = 'R0 %')


## Three-group model (each has different contact rates)
beta1 <- rep(c(.1,.15,.2,.25), each = 3)
beta24 <- .001           # p.9 (results of sim: homogenous populations)
beta5 <- .005
beta67 <- .01
dur.ac <- 2                             # months
rh.ac <- beta1/beta24
ehm.ac <- (rh.ac-1)*dur.ac

r1 <- c(.1,.8,3.2,.15,1.2,4.8,.2,1.6,6.4,.25,2,8) # Table 4
r24 <- rep(c(.33,.263,1.05), 4)
r5 <- rep(c(.81,.646,2.58), 4)
r67 <- rep(c(.126,1,4.02),4)
group <- rep(1:3, 4)
af.ac <- signif(r1 / (r1+r24+r5+r67) * 100,3)
temp <- data.frame(af.ac, ehm.ac, rh.ac, dur.ac)
## re-weight by population make up: group1 = 50% population, group2 = 40%, group3 = 10%
threegr = ddply(temp, .(ehm.ac), summarise, af.ac = af.ac[1]*.5 + af.ac[2]*.4 + af.ac[3]*.1, data = 'epidemic curve', type = 'three risk groups', estim = 'R0 %')
threegr <- data.frame(rh.ac=unique(rh.ac), dur.ac =unique(dur.ac), threegr)
jacq <- data.frame(ms = 'Jacquez et al. 1994', rbind(onegr, threegr)) 

####################################################################################################
## Pinkerton 1996, citing Jacquez
beta1 <- c(.05,.3)
beta2 <- .001
dur.ac <- 2
rh.ac <- beta1/beta2
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- c(20,55,20,87) ## Figure 7
type <- rep(c('100 sex acts per partner', '1 sex act per partner'), 2)
pink <- data.frame(ms = 'Pinkerton and Abramson 1996', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Jacquez', type = type, estim = 'R0 %')

####################################################################################################
## Koopman et al. 1997
## 
## "We chose absolute transmission probability values for early infection that generated epidemic
## curves compatible with those observed in the San Francisco Hepatitis B cohort study (Bacchetti and Moss 1989)."
## 
## Model 1: (not enough info in text to get details from model 2 which deals with people switching between
## risk groups)
beta.ac <- .2
beta.ch <- .001
dur.ac <- 1.5
rh.ac <- beta.ac/beta.ch
ehm.ac <- (rh.ac - 1 ) * dur.ac
af.ac <- c(36, 47.4)                    # Table 1
type <- c('homogenous contact', 'age-peaked contact rate')
koop <- data.frame(ms = 'Koopman et al. 1997', rh.ac, dur.ac, ehm.ac, af.ac, data = 'epidemic curve', type = type, estim = 'R0 %')

####################################################################################################
## Kretzschmar & Dietz 1998, citing Jacquez
beta.ac <- .3
beta.ch <- .001
dur.ac <- 1/7.3*12 ## transition rate is 7.3 / yr for some reason (Figure 8 caption)
rh.ac <- beta.ac/beta.ch
ehm.ac <- signif((rh.ac-1)*dur.ac,3)
## at endemic equilibrium in Figure 8 "The fraction of new cases produced by stage 1 individuals is
## now monotone de- creasing for model M(+/v) until in the endemic state around 65% of cases are
## produced by stage 1 individuals."
af.ac <- 65
kret <- data.frame(ms = 'Kretzschmar & Dietz 1998', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Jacquez', type = '', estim = 'AF_acute')

####################################################################################################
## Xiridou, citing studies on viral load
dur.ac <- runif(10^5, 1,5) ## uniformly sampled from this range
rh.ac <- runif(10^5, 5,35) ## uniformly sampled from this range
ehm.ac <- (rh.ac-1)*dur.ac
rh.ac <- quantile(rh.ac, c(.025, .975))
dur.ac <- quantile(dur.ac, c(.025, .975))
ehm.ac <- quantile(ehm.ac, c(.025, .975)) ## get 95% quantiles to match the latin hypercube thingy
ehm.ac <- signif(ehm.ac,3)
af.ac <- c(2.71, 24.97)
xir <- data.frame(ms = "Xiridou et al. 2004", rh.ac, dur.ac, ehm.ac, af.ac, data = 'viral load studies', type = '', estim = 'AF_acute')


####################################################################################################
## Hayes et al., citing Wawer
rh.ac <- 7.25
dur.ac <- 5
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- c(41,23)
hayes <- data.frame(ms = 'Hayes et al. 2006', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Rakai', type = '', estim = 'crude')


####################################################################################################
## Pinkerton 2007
dur.ac <-  c(42,49,56)/30     ## months (49 days, varied 42-56)
rh.ac <- c(4.2,8.1,12) ## varied 4.2-12
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- c(2.5, 5.5, 9) ## diagonal of Table 2
pink2 <- data.frame(ms = 'Pinkerton 2007', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Pilcher and Rapatski', type = '', estim = 'AF_acute')

####################################################################################################
## Hollingsworth, citing Wawer
dur.ac <- 1.9
rh.ac <- 26
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- c(9, 31)
type <- c('serial monogamy', 'random mixing')
holl <- data.frame(ms = 'Hollingsworth et al. 2008', rh.ac, dur.ac, ehm.ac, af.ac, data = NA, type = '', estim = 'AF_acute')


####################################################################################################
## Abu-Radad, citing Wawer
dur.ac <- 2.5
rh.ac <- .107/.008
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- c(17, 25)
type <- c('Kisumu', 'Yaounde')
abu <- data.frame(ms = 'Abu-Raddad et al. 2008', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Rakai', type = '', estim = 'AF_acute')


####################################################################################################
## Salomon, citing Wawer
dur.ac <- 4
rh.ac <- .0082/.001
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- c(41, 23)                      # Fig 3 no treatment 
type <- c('single lifelong partner', 'multiple partners')
sal <- data.frame(ms = 'Salomon & Hogan 2008', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Rakai', type = '', estim = 'AF_acute')

####################################################################################################
## Prabhu, citing Pinkerton
dur.ac <-  49/30
rh.ac <- 8.1
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- 11.4
prab <- data.frame(ms = 'Prabhu et al. 2009', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Pinkerton', type = '', estim = 'AF_acute')

####################################################################################################
## Powers, citing Hollingsworth
dur.ac <-  4.8
rh.ac <- 30.3
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- 38.4
powers <- data.frame(ms = 'Powers et al. 2011', rh.ac, dur.ac, ehm.ac, af.ac, data = 'Rakai', type = '', estim = 'AF_acute')


####################################################################################################
## Williams, viral load
dur.ac <-  1
rh.ac <- 3
ehm.ac <- (rh.ac-1)*dur.ac
af.ac <- 2
will <- data.frame(ms = 'Cohen et al. 2013 (Williams)', rh.ac, dur.ac, ehm.ac, af.ac, data = 'viral load', type = '', estim = 'AF_acute')

out <- rbind(jacq, pink, koop, kret, xir, pink2, hayes, holl, abu, sal, prab, powers, will)
out$ms <- factor(out$ms)
out$col <- makeTransparent(rainbow(nlevels(out$ms))[out$ms], 70)

out <- ddply(out, .(ms), transform, tag = paste0(as.numeric(ms)[1], letters[1:length(ms)]))

pdf('AF vs EHM review.pdf', w = 3.5, h=4)
par(mar=c(4,4,1,1))
plot(out$ehm.ac, out$af.ac, col = out$col, pch = 19, xlim = c(0, 170), ylim = c(0,50), cex = 2,
     xlab = expression(EHM[acute]), ylab = expression(AF[acute]), bty = 'n')
text(out$ehm.ac, out$af.ac, out$tag, cex = .6)
##legend('topleft', leg = unique(out$ms), col = unique(out$ms), pch = 19)
graphics.off()

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

write.csv(out, file = 'Cohen rev.csv')
