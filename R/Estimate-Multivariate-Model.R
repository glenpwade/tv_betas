#  Initialise  ----

rm(list=ls())
library(MTVGARCH)
library(knitr)

estCtrl = list(verbose=TRUE,calcSE=TRUE)

setwd("C:/Source/Repos/tv_betas")

# Get Data  ----
# Get all the final estimated univariate models & compile them into a multivariate ntvgarch Object
if(FALSE){
    tvg_list = list()
    tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
    tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
    tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
    tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")
    tvg_list[[5]] = readRDS("Results/STW_Final_TVG_model.RDS")
    tvg_list[[6]] = readRDS("Results/SPY_Final_TVG_model.RDS")
    tvg_list[[7]] = readRDS("Results/PR_Final_TVG_model.RDS")
    tvg_list[[8]] = readRDS("Results/CRAU_Final_TVG_model.RDS")
    tvg_list[[9]] = readRDS("Results/CRUS_Final_TVG_model.RDS")
    
    # Make a multivaraite model object ----
    
    tvg_names = c("ANZ","CBA","NAB","WBC","STW","SPY","PR","CRAU","CRUS")
    tvBetas <- ntvgarch(tvg_list,tvg_names)
    
    saveRDS(tvBetas,"Results/multivar-spec.RDS")  
}

tvBetas = readRDS("Results/multivar-spec.RDS")

# Estimate the multivariate model ----

# "1/11/2019" -> Index 2165

# TV Obj ----

# TV only estimate, h(t) = 1
ptitle = "ANZ"
prices <- read.csv("data/tv_betas_prices.csv")
e_anz <- diff(log(as.numeric(prices$ANZ)) ) * 100  # Percentage Returns
e <- e_anz
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(1,-1,6,-5)
TV$pars["speedN",] = c(2,2,5,4)
TV$pars["locN1",] = c(0.35,0.4,0.7,0.75)
TV$optimcontrol$reltol = 1e-9
TV$optimcontrol$ndeps = rep(1e-6,13)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

TV <- estimateTV(e,TV,estCtrl)

plot(TV)
summary(TV)

# TVGARCH Obj ----

ANZ_tvg <- tvBetas$ANZ
# estimated tvpars
ANZ_tvg$Estimated$tv$delta0
ANZ_tvg$Estimated$tv$pars

# Default plot of tvgarch object
plot(tvBetas$ANZ)
# summary of tvgarch object
summary(tvBetas$ANZ)

# Plot the returns, with g & h overlay ----
tvBetas$ANZ$e_desc <- "ANZ % Returns"
ptitle = tvBetas$ANZ$e_desc
e = tvBetas$ANZ@e
plot(e,type='l',col="grey40",main=ptitle)
abline(v=seq(1,3200,100),col="grey80")    
lines(e,type='l',col="grey40")

g = tvBetas$ANZ$Estimated$g
lines(sqrt(g),type='l',col="red",lwd=2)

h = tvBetas$ANZ$Estimated$h
lines(sqrt(h),type='l',col="blue")

gh = sqrt(g) + sqrt(h)
lines(gh,type='l', col="green")

summary(gh)
var(e)


# Test for CCC: ----

tvg_list = list()
tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")

## Banks Only ----

#tvgBanks <- tvBetas[1:4]

tvg_names = c("ANZ","CBA","NAB","WBC")
tvgBanks <- ntvgarch(tvg_list,tvg_names)

e <- matrix()
for (n in 1:tvgBanks@N){
    if(n==1) e <- tvgBanks[[1]]@e
    else e = cbind(e,tvgBanks[[n]]@e)
}

Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
NullHyp = ccc(tvgBanks@N,tvgBanks)
NullHyp$P = diag(nrow = tvgBanks@N)
NullHyp = estimateCCC(e,NullHyp,estCtrl)
testOrd = 1

ccctest <- test.CCCParsim(e,NullHyp,st,testOrd)
print(round(ccctest,5))

stcctest = test.CCCvSTCC1(e,NullHyp,st,testOrd)
print(round(stcctest,5))

## Banks Plus STW ----

tvg_list = list()
tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")
tvg_list[[5]] = readRDS("Results/STW_Final_TVG_model.RDS")

tvg_names = c("ANZ","CBA","NAB","WBC","STW")
tvgBanks1 <- ntvgarch(tvg_list,tvg_names)

e <- matrix()
for (n in 1:tvgBanks1@N){
    if(n==1) e <- tvgBanks1[[1]]@e
    else e = cbind(e,tvgBanks1[[n]]@e)
}

Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
NullHyp = ccc(tvgBanks1@N,tvgBanks1)
NullHyp$P = diag(nrow = tvgBanks1@N)
NullHyp = estimateCCC(e,NullHyp,estCtrl)
testOrd = 1

ccctest <- test.CCCParsim(e,NullHyp,st,testOrd)
print(round(ccctest,5))

stcctest = test.CCCvSTCC1(e,NullHyp,st,testOrd)
print(round(stcctest,5))

# Estimate the STCC model ----

tvg_list = list()
tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")

## Banks Only ----

tvg_names = c("ANZ","CBA","NAB","WBC")
tvgBanks <- ntvgarch(tvg_list,tvg_names)

e <- matrix()
for (n in 1:tvgBanks@N){
    if(n==1) e <- tvgBanks[[1]]@e
    else e = cbind(e,tvgBanks[[n]]@e)
}

Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
z <- filterData(e,tvgBanks)
NullHyp = stcc1(z,tvgBanks)
NullHyp$P = diag(nrow = tvgBanks@N)
NullHyp$pars[1] <- 5
NullHyp$pars[2] <- 0.6
# Switch P1 & P2
NullHyp$P1 -> P1tmp
NullHyp$P2 -> P2tmp
NullHyp$P1 <- P2tmp
NullHyp$P2 <- P1tmp

NullHyp$optimcontrol$reltol = 1e-04
NullHyp = estimateSTCC1(z,NullHyp,estCtrl)
testOrd = 1

stcctest1 = test.TVCC1vTVCC2(e,NullHyp,testOrd)

# TODO: Tidy up ntvgarch object to avoid this shuffling!!
for(n in seq_along(NullHyp$ntvgarch)){
    NullHyp$ntvgarch[[n]]$Estimated = list()
    NullHyp$ntvgarch[[n]]$Estimated$tv = list()
    NullHyp$ntvgarch[[n]]$Estimated$tv$g <- NullHyp$ntvgarch[[n]]$tv$g
    
    NullHyp$ntvgarch[[n]]$Estimated$garch <- garch(garchtype$gjr)
    NullHyp$ntvgarch[[n]]$Estimated$garch$h <- NullHyp$ntvgarch[[n]]$garch$h
    
    #NullHyp$ntvgarch[[n]]$tv <- tvg@tvObj
    #NullHyp$ntvgarch[[n]]$tv$g <- tvg@tvObj@g
    
    #plot(NullHyp$ntvgarch[[n]]$tv)
}

stcctest2 = test.TVCC1vTVCC2(e,NullHyp,testOrd)
print(round(stcctest,5))



# Check Test Consistency ----

tvg_list = list()

tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/NAB_Final_TVG_model.RDS")

## 2 Banks Only ----

tvg_names = c("ANZ","NAB")
tvgBanks <- ntvgarch(tvg_list,tvg_names)

e <- matrix()
for (n in 1:tvgBanks@N){
    if(n==1) e <- tvgBanks[[1]]@e
    else e = cbind(e,tvgBanks[[n]]@e)
}

plot(tvgBanks[[1]]@e,type='l')

Tobs = NROW(e)
st = seq(-0.5,0.5,length.out=Tobs)
# Null Hypothesis
NullHyp = ccc(tvgBanks@N,tvgBanks)
NullHyp = estimateCCC(e,NullHyp,estCtrl)

# Alternate Hypothesis
z <- filterData(e,tvgBanks)
H_1 <- stcc1(z,tvgBanks)
#
testOrd = 1

ccctest <- test.CCCParsim(e,NullHyp,st,testOrd)
print(round(ccctest,5))

stcctest = test.CCCvSTCC1(e,NullHyp,st,testOrd)
print(round(stcctest,5))

all.equal(ccctest,stcctest)
# Tests are consistent


# Test for CEC: ----

tvg_list = list()
tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")

## Banks Only ----

#tvgBanks <- tvBetas[1:4]

tvg_names = c("ANZ","CBA","NAB","WBC")
tvgBanks <- ntvgarch(tvg_list,tvg_names)

e <- matrix()
for (n in 1:tvgBanks@N){
    if(n==1) e <- tvgBanks[[1]]@e
    else e = cbind(e,tvgBanks[[n]]@e)
}

Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
testOrd = 1
#
NullHyp = cec(tvgBanks@N,tvgBanks)
NullHyp = estimateCEC(e,NullHyp,estCtrl)
NullHyp$Estimated$P = diag(x=NullHyp$Estimated$pars,tvgBanks@N,tvgBanks@N)
cectest <- test.CCCParsim(e,NullHyp,st,testOrd)
print(round(cectest,5))
#
NullHyp = ccc(tvgBanks@N,tvgBanks)
NullHyp = estimateCCC(e,NullHyp,estCtrl)
ccctest <- test.CCCParsim(e,NullHyp,st,testOrd)
print(round(ccctest,5))

stcctest = test.CCCvSTCC1(e,NullHyp,st,testOrd)
print(round(stcctest,5))

