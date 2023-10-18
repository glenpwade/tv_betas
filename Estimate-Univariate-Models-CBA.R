#  Initialise  ----

rm(list=ls())
library(MTVGARCH)
library(knitr)

#setwd("D:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project4_NAO/R")   #Anna's Laptop - GOOGLE DRIVE
setwd("C:/Source/Repos/tv_betas")
#
estCtrl = list(verbose=TRUE,calcSE=TRUE)
noisedist = list(name='Normal', mean=0, sd=1)
#
# Set simulation controls
simcontrol = list()
simcontrol$numLoops = 1200     # Number of series to generate - to simulate the error distribution, for testing
simcontrol$numCores = 7
maxTestOrd = 3
simcontrol$maxTestorder=maxTestOrd

## noiseDist is a named-list describing the error-distribution and parameters
## e.g. noiseDist$name = 'Normal'     noiseDist$mean = 0  noiseDist$sd = 1
## or   noiseDist$name = 'Student-t'  noiseDist$df = 6    noiseDist$ncp = 0
#
## mtvgarch pkg needs an update to make this optional & std Normal by default

# Get Data  ----

prices <- read.csv("data/tv_betas_prices.csv")
e_cba <- diff(log(as.numeric(prices$CBA)) ) * 100  # Percentage Returns
#
# Replace all NA's with the previous valid entry if required:
if(FALSE){
    na_pos = which(is.na(e_cba) )
    for(n in seq_along(na_pos)){
        t_pos <- na_pos[n]
        e_cba[t_pos] <- e_cba[t_pos-1]
    }
}

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_cba[1350:2250]
Tobs = NROW(e)
ptitle = "CBA stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_cba[1400:1800]
# Est se1 sig
# omega 0.072941 NaN    
# alpha 0.070489 NaN    
# beta  0.852360 NaN 

# e_cba[2500:3150]
# Est      se1 sig
# omega 0.132530 0.088212    
# alpha 0.061153 0.027582 ** 
# beta  0.839795 0.083945 ***

# Try the rolling window method:
if(FALSE){
    estCtrl$vartargetWindow = 300
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    # win=300, a=0.097559, b=0.847240, a+b=0.945
    # win=400, a=0.082415, b=0.894717, a+b=0.977
    # win=500, a=0.078079, b=0.905088, a+b=0.983
    
    # Method:  MLE, variance-targetting a rolling Window of 300 observations 
    # Est      se1 sig
    # omega 0.073229       NA    
    # alpha 0.097559 0.023529 ***
    # beta  0.847240 0.049561 ***
}

# Next, We need a standard TV object to generate the data:
e <- e_cba
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_CBA <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist, seed=1)
ptitle = "CBA"
saveRDS(refData_CBA,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_cba
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_CBA

# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# RESTART HERE with an updated TV Model specification:   ####

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
RefTests = list()
for (n in  1:maxTestOrd){
  RefTests[[n]] <- list()
  RefTests[[n]] <- getTestStats(e,TV,n)
}

## Calculate p_value by simulating the distrbution:

# Note: testStatDist() will generate a simulated test distribution and calculate p_values as follows:
# 1. For each column of refData, estimate the null TV model, then calculate test statistics (TR2 & Robust) for each requested TestOrder
# 2. Then compare actual RefTest statistics with the TestStat-distributions (i.e. get p-values)
# 3. This takes approx. 3 minutes for the TV$delta0 specification running in local-cpu parallel mode on an i7 3GHz 4-core/8-processor CPU,  Win10PC
SIMRESULT = testStatDist(refData,TV,RefTests,simcontrol)

# Print the test Results:
testResults = data.frame()
for(n in 1:maxTestOrd){
  testResults[n,1] = n
  testResults[n,2] = SIMRESULT[[n]]$pVal_TR2
  testResults[n,3] = SIMRESULT[[n]]$pVal_ROB
}
knitr::kable(testResults,'pipe',digits=3,col.names = c('Test Ord','TR2','Robust'))
  

# View the Distributions, if you like: 
if(FALSE){
  
  hist(SIMRESULT$TestStatDist[,"Stat_TR2.2"],breaks = 20)
  abline(v=RefTests[[2]]$TR2,col="red")
  
  hist(SIMRESULT$TestStatDist[,"Stat_Robust.2"],breaks = 20)
  abline(v=RefTests[[2]]$Robust,col="red")
  
  #View(SIMRESULT$TestStatDist)     # Have a look here for Troubleshooting
}

## SAVE Results
saveFile = paste0("Results/TestResults_",simcontrol$saveAs)
saveRDS(SIMRESULT,saveFile)

## RELOAD Results
#SIMRESULT = readRDS(saveFile)

## RESULTS Section ----

## P-Values from TEST Results, TV-delta0 only, CBA[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.002|  0.000|
# |        2| 0.013|  0.005|
# |        3| 0.000|  0.005|

## Conclusion: 
## Evidence of a bump first and/or third order.
## We can try to identify & estimate a single transition, or split the data and retest:
## Let's try estimating a single order transition...

e <- e_cba
Tobs = NROW(e)
#plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_CBA
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$single)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-1_Trans, CBA[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.141|  0.342|
# |        2| 0.007|  0.045|
# |        3| 0.000|  0.003|


## Conclusion: 
## Evidence of another transition exits
## Let's try estimating a model with 2 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-2_Trans, CBA[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.050|  0.375|
# |        2| 0.012|  0.013|
# |        3| 0.010|  0.000|

## Conclusion: 
## Still Evidence of another transition exits
## Let's try estimating a model with 3 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## Default starting values don't work well for this 3-Trans model, so...
## Looking at the plot, the missing transition seems to be high-to-low around Obs 2500
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(1,8,-7)
TV$pars["speedN",] = c(4,6,6)
TV$pars["locN1",] = c(0.2,0.6,0.8)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-3_Trans, CBA[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.829|  0.347|
# |        2| 0.836|  0.449|
# |        3| 0.405|  0.027|

#
# Log-likelihood value(TV):  -4910.034

## Conclusion: 
## Small Evidence of another transition exits - Robust TestOrd-3
## Let's try estimating a model with 4 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
abline(v=seq(0,3200,by=100),col="grey")
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## Default starting values don't work well for this 4-Trans model, so...
## Looking at the plot, the missing transition seems to be around Obs 1000 
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(1,-1,6,-5)
TV$pars["speedN",] = c(2,2,5,4)
TV$pars["locN1",] = c(0.35,0.4,0.7,0.75)
TV$optimcontrol$reltol = 1e-7
TV$optimcontrol$ndeps = rep(1e-6,13)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
#
## P-Values from TEST Results, TV-4_Trans, CBA[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.236|  0.012|
# |        2| 0.081|  0.004|
# |        3| 0.061|  0.011|

    
## Conclusion: 
## Some Evidence of another transition exits.
## Maybe a chance with Test Ord 2 being 0.004, but...
## it was difficult to find starting pars & optim controls to estimate this model
## so we will stop here.  Can't visually identify where another transition might be.

## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 1 1 1 
# 
# Estimation Results:
# 
# Delta0 = 1.735253    se0 = 0.256453*** 
# 
#             st1      se1 sig       st2      se2 sig.1      st3      se3 sig.2       st4      se4 sig.3
# deltaN 1.176916 0.173243 *** -2.534550 0.118584   *** 7.958872 0.855070   *** -7.159399 0.857084   ***
# speedN 4.464657 0.248553 ***  1.461050 0.107982   *** 6.999998 0.000000   ***  5.171512 0.359494   ***
# locN1  0.303670 0.006313 ***  0.313629 0.065846   *** 0.708183 0.000991   ***  0.732232 0.002721   ***
# locN2        NA      NaN            NA      NaN             NA      NaN              NA      NaN      
# 
# Log-likelihood value(TV):  -4891.686



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full CBA data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

prices <- read.csv("data/tv_betas_prices.csv")
e_cba <- diff(log(as.numeric(prices$CBA)) ) * 100  # Percentage Returns
e <- e_cba
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "CBA"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "CBA Std.% Returns"

TVG$tvpars[,1] = c(0.8,2.0,0.33,NA)
TVG$tvOptimcontrol$reltol = 5e-04
TVG$tvOptimcontrol$ndeps = rep(1e-04,length(TVG$tvOptimcontrol$ndeps))
TVG$garchpars[,1] = c(0.05,0.02,0.7,0.05)
TVG$garchOptimcontrol$reltol = 1e-04
#TVG$garchOptimcontrol$parscale = c(5,1,50,10)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
# Re-run the estimator until fully converged:
TVG <- estimateTVGARCH(e,TVG,estCtrl)
# Now Save!
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

summary(TVG1)
plot(TVG1)

