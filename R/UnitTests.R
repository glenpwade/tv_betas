#  Initialise  ----

rm(list=ls())
library(MTVGARCH)
library(knitr)

estCtrl = list(verbose=TRUE,calcSE=TRUE)

setwd("C:/Source/Repos/tv_betas")

# # Generate Data with known TV, GARCH & Corr ----
nr.obs = 2000
nr.series = 3

# Generate Noise Data:
set.seed(42)
noiseData <- list()
for(n in 1:10){
    # Normal Error/Noise Distribution (Default is Standard-Normal)
    noiseData[[n]] <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)    
}

u = noiseData[[5]]
# Step1: Create CCC Correlation 
# - - - CCC - - -
COR1 = ccc(nr.series)
P.sqrt <-  sqrt_mat1(COR1$P)
e <- t(P.sqrt %*% t(u))
plot(e[,3],type='l')

# Step2: Inject GARCH into Data
garchObj = garch(garchtype$general)

# Generate Discard Data:
discardObs <- 1500
# Normal Error/Noise Distribution (Default is Standard-Normal)
discardData <- matrix(rnorm(discardObs*nr.series),nrow=discardObs, ncol=nr.series)

garchObj$pars["omega",1] <- ( 1 - garchObj$pars["alpha",1] - garchObj$pars["beta",1] )
e <- rbind(discardData,e)
endRow <- discardObs + nr.obs

for (b in 1:nr.series){
    w <- z <- e[,b]
    ht_1 <- 1
    w[1] <- z[1]
    for (t in 2:endRow) {
        ht <- garchObj$pars["omega",1] + garchObj$pars["alpha",1]*(w[t-1])^2 + garchObj$pars["beta",1]*ht_1
        if(garchObj$type == garchtype$gjr) { ht <- ht + garchObj$pars["gamma",1]*(min(w[t-1],0)^2) }
        ht_1 <- ht
        w[t] <- sqrt(ht)*z[t]
    }
    e[,b] <- as.numeric(w)
}

# Discard the first 2000
startRow <- discardObs + 1
e <- e[(startRow:endRow), ]
plot(e[,3],type='l')

# Step3: Inject TV into Data
st = seq(0,1,length.out=nr.obs)
TV1 = tv(st,tvshape$single)
TV1$pars["deltaN",1] = 3
TV1$pars["locN1",1] = 0.3
#
TV2 = tv(st,tvshape$single)
TV2$pars["deltaN",1] = 3
TV2$pars["locN1",1] = 0.5
#
TV3 = tv(st,tvshape$single)
TV3$pars["deltaN",1] = 3
TV3$pars["locN1",1] = 0.7

gt1 <- get_g(TV1)
plot(gt1,type='l')
e[,1] <- e[,1]*sqrt(gt1)
#
gt2 <- get_g(TV2)
plot(gt2,type='l')
e[,2] <- e[,2]*sqrt(gt2)
#
gt3 <- get_g(TV3)
plot(gt3,type='l')
e[,3] <- e[,3]*sqrt(gt3)
plot(e[,3],type='l')

# Data is now generated! ----

#  Unit Test the Tests! ----

TV1 <- estimateTV(e[,1],TV1,estCtrl)
TV2 <- estimateTV(e[,2],TV1,estCtrl)
TV3 <- estimateTV(e[,3],TV1,estCtrl)
plot(TV1)
summary(TV1)

TVG1 <- tvgarch(TV1,garchtype$general)
TVG2 <- tvgarch(TV2,garchtype$general)
TVG3 <- tvgarch(TV3,garchtype$general)

TVG1 <- estimateTVGARCH(e[,1],TVG1)
TVG1 <- estimateTVGARCH(e[,1],TVG1)
TVG2 <- estimateTVGARCH(e[,2],TVG2)
TVG2 <- estimateTVGARCH(e[,2],TVG2)
TVG3 <- estimateTVGARCH(e[,3],TVG3)
TVG3 <- estimateTVGARCH(e[,3],TVG3)


tvg_list = list()
tvg_list[[1]] = TVG1
tvg_list[[2]] = TVG2
tvg_list[[3]] = TVG3

tvg_names = c("TVG1","TVG2","TVG3")
tvgObj <- ntvgarch(tvg_list,tvg_names)

Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
NullHyp = ccc(tvgObj@N,tvgObj)
NullHyp = estimateCCC(e,NullHyp,estCtrl)
testOrd = 1

ccctest <- test.CCCParsim(e,NullHyp,st,testOrd)
print(round(ccctest,5))

z <- filterData(e,tvgObj)
H_1 <- stcc1(z,tvgObj)

stcctest = test.CCCvSTCC1(e,NullHyp,H_1,testOrd)
print(round(stcctest,5))

summary(TVG1)
