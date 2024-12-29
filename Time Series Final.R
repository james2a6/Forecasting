library(tsDyn)
library(vars)
library(quantmod)
library(tseries)
library(urca)
library(dynlm)


#loading the data
df <- read.csv("FinalTS.csv", header = TRUE)

#### CONVERT TO TIME SERIES
dfts = ts(df,  start = c(1992, 1), frequency = 12)

plot.ts()

un<-dfts[,2]
re<-dfts[,3]
cpi<-dfts[,4]
sen<-dfts[,5]
plot(un)
plot(cpi)
plot(sen)
lre = log(re)



source(file="intord.R")
intord(un)    #I1
intord(lre)  #I1
intord(cpi)  #I1
intord(sen)   #I1

#first differencing to make them all I(0)
dun<- diff(un)   #I0
dlre<- diff(lre)  #I0
dcpi<- diff(cpi)  #I0
dsen<- diff(sen)  #I0


n<- length(lre)


#variables I1
ly<- cbind(un, cpi, sen)

# stationary variables -- I(0)
dy <- cbind(dun, dcpi, dsen)

library(urca)
library(vars)

# lag selection criteria
VARselect(dy, lag.max=12, type="const")

## cointegration - variables must be in levels! #
jc <- ca.jo(ly, type="eigen", ecdet="const",K=2) 
summary(jc)
jct <- ca.jo(ly, type="trace", ecdet="const",K=2) 
summary(jct)

# eigenvalues, etc., can be taken from the output of ca.jo
eigenvals <- jc@lambda

# cointegrating relationships
cointv <- jc@V
cointj <- cointv[,1]
yym <- as.matrix(ly)
ecmj <- yym%*%cointj[1:3] + cointj[4] # eigenvectors with data plus intercept
intord(ecmj)

# Engle-Granger
c1 <- dynlm(sen~un+cpi) # 
summary(c1)
ecm <- c1$residuals
#dev.off()
#par(mar = c(1, 1, 1, 1))
intord(ecm)


plot(ecm,type='l')
plot(ecmj,type='l')

#create a lag for ECM ter,
ec <- embed(ecm,2) # first lag
ecm1 <- ec[,2]
ecj <- embed(ecmj,2) # first lag
ecmj1 <- ecj[,2]



# VECM since appear cointegrated, except for possibly gov. spending

library(vars)
var3 <- VAR(dy, p=2, type="cons", season = 12) # type ="cons","trend","both"
summary(var3)

plot(var3, names = "dun")
plot(var3, names = "dcpi")
plot(var3, names = "dsen")


# VECM with Johansen
var4 <- VAR(dy, p=2, type="cons",exogen=ecmj1,season=12)
summary(var4)

# VECM with EG
var5 <- VAR(dy, p=2, type="cons",exogen=ecm1,season=12)
summary(var5)


plot(var4, names = "dun")
plot(var4, names = "dcpi")
plot(var4, names = "dsen")


# test for serial correlation for residuals

# BoX-Ljung Q Statistic for dun

resi = var4$varresult$dun$residuals
b = Box.test(resi,lag = 20, type="Ljung-Box")
b

blt = rep(0,20)
for (i in 1:20){
  b = Box.test(resi,lag = i, type="Ljung-Box")
  blt[i]=b$p.value
}
blt

# BoX-Ljung Q Statistic for dcpi

resi = var4$varresult$dcpi$residuals
b = Box.test(resi,lag = 20, type="Ljung-Box")
b

blt = rep(0,20)
for (i in 1:20){
  b = Box.test(resi,lag = i, type="Ljung-Box")
  blt[i]=b$p.value
}
blt


# BoX-Ljung Q Statistic for dsen

resi = var4$varresult$dsen$residuals
b = Box.test(resi,lag = 20, type="Ljung-Box")
b

blt = rep(0,20)
for (i in 1:20){
  b = Box.test(resi,lag = i, type="Ljung-Box")
  blt[i]=b$p.value
}
blt



#ly<- cbind(un, cpi, sen)
vecm1 <- ca.jo(ly, ecdet = "trend", type="eigen", K=2, spec="longrun",
               season=12)
ve1 <- cajools(vecm1,reg.number=2)

# restricted
lyr <- cbind(un,cpi)    
vecmr <- ca.jo(lyr, ecdet = "const", type="eigen", K=2, spec="longrun",
               season=12)
ve1r <- cajools(vecmr,reg.number=2)  

# Joint F-test for Granger causality
anova(ve1, ve1r, test="F") 



irfs1a <- irf(var4, impulse = "dun", response = c("dun","dsen"), boot =
              TRUE)
# par(mar = c(1, 1, 1, 1))
plot(irfs1a)

irfs1b <- irf(var4, impulse = "dun", response = c("dcpi"), boot =
               TRUE)
# par(mar = c(1, 1, 1, 1))
plot(irfs1b)

irfs1c <- irf(var4, impulse = "dun", response = c("dun","dcpi", "dsen"), boot =
                TRUE)
# par(mar = c(1, 1, 1, 1))
plot(irfs1c)

# IRFs
irfs2a <- irf(var4, impulse = "dcpi", response = c("dun","dcpi","dsen"), boot =
              TRUE)
# par(mar = c(1, 1, 1, 1))
plot(irfs2a)

irfs2b <- irf(var4, impulse = "dcpi", response = c("dun"), boot =
               TRUE)
plot(irfs2b)

irfs2c <- irf(var4, impulse = "dcpi", response = c("dcpi"), boot =
                TRUE)
# par(mar = c(1, 1, 1, 1))
plot(irfs2c)

irfs2d <- irf(var4, impulse = "dcpi", response = c("dsen"), boot =
                TRUE)
plot(irfs2d)


irfs3a <- irf(var4, impulse = "dsen", response = c("dun", "dcpi","dsen"), boot =
                TRUE)
plot(irfs3a)

irfs3b <- irf(var4, impulse = "dsen", response = c("dun"), boot =
                TRUE)
plot(irfs3b)

irfs3c <- irf(var4, impulse = "dsen", response = c("dcpi"), boot =
                TRUE)
plot(irfs3c)

irfs3d <- irf(var4, impulse = "dsen", response = c("dsen"), boot =
                TRUE)
plot(irfs3d)

# Variance Decompositions (contribution of each variable to predicting a variable)
vard <- fevd(var4, n.ahead=12)
vard
vard$dsen
plot(vard, col=1:3)


# Forecasting with a VECM
varf <- vec2var(vecm1, r = 1)
fcast <- predict(varf, n.ahead = 8, ci = 0.95) 
plot(fcast)
fcast

# plotting forecasts for one variable (lm)
lhfaf <- fcast$fcst$sen[,1]
lhfaflow <- fcast$fcst$sen[,2]
lhfafupp <- fcast$fcst$sen[,3]

ff <- cbind(lhfaf,lhfaflow,lhfafupp)
matplot(ff,col=c(1,2,2),lty=1,lwd=2,type='l')
