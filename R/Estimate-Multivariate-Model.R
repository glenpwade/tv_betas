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
    #tvg_list[[7]] = readRDS("Results/PR_Final_TVG_model.RDS")
    tvg_list[[7]] = readRDS("Results/PR_Final_TVG_model_StdGarch.RDS")
    tvg_list[[8]] = readRDS("Results/CRAU_Final_TVG_model.RDS")
    tvg_list[[9]] = readRDS("Results/CRUS_Final_TVG_model.RDS")
    
    # Make a multivaraite model object ----
    
    tvg_names = c("ANZ","CBA","NAB","WBC","STW","SPY","PR","CRAU","CRUS")
    tvBetas <- ntvgarch(tvg_list,tvg_names)
    
    saveRDS(tvBetas,"Results/multivar-spec.RDS")  
}

tvBetas = readRDS("Results/multivar-spec.RDS")

# Plot the returns, with g & h overlay ----
tvBetas$WBC$e_desc <- "WBC % Returns"
ptitle = tvBetas$WBC$e_desc
e = tvBetas$WBC@e
plot(e,type='l',col="grey40",main=ptitle)
abline(v=seq(1,3200,100),col="grey80")    
lines(e,type='l',col="grey40")

g = tvBetas$WBC$Estimated$g
lines(sqrt(g),type='l',col="red",lwd=2)

h = tvBetas$WBC$Estimated$h
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

tvg_names = c("ANZ","CBA","NAB","WBC")
tvgBanks <- ntvgarch(tvg_list,tvg_names)

Tobs = tvgBanks$Tobs
st = seq(0,1,length.out=Tobs)
NullHyp = ccc(tvgBanks)
NullHyp = estimateCCC(NullHyp,estCtrl)
round(NullHyp$Estimated$P,3)
testOrd = 1

ccctest <- test.CCCParsim(NullHyp,st,testOrd)
print(round(ccctest,5))

# Pairwise tvBetas:





## Banks Only ----

tvg_list = list()
tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")

tvg_names = c("ANZ","CBA","NAB","WBC")
tvgBanks <- ntvgarch(tvg_list,tvg_names)

e <- matrix()
for (n in 1:tvgBanks@N){
    if(n==1) e <- tvgBanks[[1]]@e
    else e = cbind(e,tvgBanks[[n]]@e)
}

Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
NullHyp = ccc(tvgBanks)
NullHyp = estimateCCC(NullHyp,estCtrl)
testOrd = 1

ccctest <- test.CCCParsim(NullHyp,st,testOrd)
print(round(ccctest,5))

stcctest = test.CCCvSTCC1(NullHyp,st,testOrd)
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

Tobs = tvgBanks1@N
st = seq(0,1,length.out=Tobs)
NullHyp = ccc(tvgBanks1)
NullHyp = estimateCCC(NullHyp,estCtrl)
testOrd = 1

ccctest <- test.CCCParsim(NullHyp,st,testOrd)
print(round(ccctest,5))

stcctest = test.CCCvSTCC1(NullHyp,st,testOrd)
print(round(stcctest,5))

# Estimate the multivariate model ----
# Estimate the STCC model ----

tvg_list = list()
tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")

## Banks Only ----

tvg_names = c("ANZ","CBA","NAB","WBC")
tvgBanks <- ntvgarch(tvg_list,tvg_names)

# e <- matrix()
# for (n in 1:tvgBanks@N){
#     if(n==1) e <- tvgBanks[[1]]@e
#     else e = cbind(e,tvgBanks[[n]]@e)
# }

Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
NullHyp = stcc1(tvgBanks)
NullHyp$pars[1] <- 4
NullHyp$pars[2] <- 0.6
# Switch P1 & P2  # TODO:  Why are we doing this??
NullHyp$P1 -> P1tmp
NullHyp$P2 -> P2tmp
NullHyp$P1 <- P2tmp
NullHyp$P2 <- P1tmp

NullHyp$optimcontrol$reltol = 1e-04
NullHyp = estimateSTCC1(NullHyp,estCtrl)

saveRDS(NullHyp,"Results/multivar-4banks-stcc1.RDS")

NullHyp <- readRDS("Results/multivar-4banks-stcc1.RDS")

testOrd = 1
stcctest1 = test.TVCC1vTVCC2(NullHyp,testOrd)
testOrd = 2
stcctest2 = test.TVCC1vTVCC2(NullHyp,testOrd)
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

