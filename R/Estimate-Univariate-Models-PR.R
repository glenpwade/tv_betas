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
ptitle = "PR"
allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_pr <- allData$PR

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_pr
e <- e_pr[2290:3152]
Tobs = NROW(e)
ptitle = "PR stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# Method:  MLE 
#            Est      se1 sig
# omega 0.044692 0.034313    
# alpha 0.063415 0.026330 ** 
# beta  0.902903 0.047262 ***
#     
#     Log-likelihood value(GARCH):  -1323.746

# Try the rolling window method:
if(FALSE){
    e <- e_pr
    Tobs = NROW(e)
    estCtrl$vartargetWindow = 400
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    
    # Method:  MLE, variance-targetting a rolling Window of 400 observations 
    #            Est      se1 sig
    # omega 0.055363       NA    
    # alpha 0.046571 0.009038 ***
    # beta  0.911268 0.020036 ***
    #     
    #     Log-likelihood value(GARCH):  -4737.781
}



# Next, We need a standard TV object to generate the data:
e <- e_pr
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_PR <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "PR"
saveRDS(refData_PR,paste0('RefData/',ptitle,'_RefData.RDS'))

#refData = readRDS(paste0('RefData/',ptitle,'_RefData.RDS'))
#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_pr
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_PR

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

## SAVE Results
saveFile = paste0("Results/TestResults_",simcontrol$saveAs)
saveRDS(SIMRESULT,saveFile)

## RELOAD Results
#SIMRESULT = readRDS(saveFile)

## RESULTS Section ----

## P-Values from TEST Results, TV-delta0 only, PR[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.799|  0.798|
# |        2| 0.923|  0.932|
# |        3| 0.790|  0.813|

## Conclusion: 
## No Evidence of a bump, though the plot suggests otherwise
## Let's try splitting the series to detect transitions...

e <- e_pr[800:2200]
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_PR
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-delta0 only,, PR[800:2200]

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.005|  0.001|
#     |        2| 0.008|  0.001|
#     |        3| 0.001|  0.002|

## Conclusion: 
## Lots of Evidence of a transitions exits 
## Based on the plot, Let's try estimating a model with 3 x single order transitions...

## Let's try to estimate a TV-3_Trans model on the full dataset:

e <- e_pr
Tobs = NROW(e)
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$double1loc, tvshape$single))
TV$delta0 = 2
TV$pars["deltaN",] = c(-0.5,-1,1)
TV$pars["speedN",] = c(4,4,5)
TV$pars["locN1",] = c(0.2,0.4,0.7)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
lines(sqrt(TV@g),col="red")
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")
# The above specification works, but the last transition may be a double...
# Log-likelihood value(TV):  -4768.274

# Shape guessed by observing the plot:
TV <- tv(st,c(tvshape$single,tvshape$double1loc,tvshape$double1loc))
TV$delta0 = 5
TV$pars["deltaN",] = c(-1,-1,-1)
TV$pars["speedN",] = c(4,4,4)
TV$pars["locN1",] = c(0.2,0.4,0.75)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-3Trans,, PR[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.991|  0.958|
# |        2| 0.495|  0.305|
# |        3| 0.541|  0.381|

## Conclusion: 
## No Evidence of another transition!

## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 3 3 
# 
# Estimation Results:
#     
#     Delta0 = 8.09636    se0 = 0.585468*** 
#     
#     st1      se1 sig       st2      se2 sig.1       st3      se3 sig.2
# deltaN -1.375687 0.365601 *** -3.755534 0.403204   *** -2.374920 0.237568   ***
#     speedN  3.224559 0.281768 ***  5.501066 0.147246   ***  4.276504 0.274627   ***
#     locN1   0.111321 0.033232 ***  0.377540 0.004398   ***  0.804204 0.009305   ***
#     locN2         NA      NaN            NA      NaN              NA      NaN      
# 
# Log-likelihood value(TV):  -4756.7



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full PR data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.


e <- e_pr
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "PR"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "PR % Returns"

#TVG$tvOptimcontrol$reltol = 1e-05
#TVG$tvOptimcontrol$ndeps = rep(1e-05,length(TVG$tvOptimcontrol$ndeps))
TVG$garchpars[,1] = c(0.1,0.05,0.85,0.01)
TVG$garchOptimcontrol$reltol = 1e-04
#TVG$garchOptimcontrol$ndeps = c(1e-05,1e-05,1e-04,1e-05)
TVG$garchOptimcontrol$parscale = c(4,1,40,2)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG@garchObj)
summary(TVG@tvObj)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
# Note: The Std Garch produces the better model specification for PR
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model_StdGarch.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model_StdGarch.RDS'))
plot(TV,main=ptitle)
summary(TVG1)
plot(TVG1)

# Try Std Garch Specification:
TVG <- tvgarch(TV,garchtype$general)
TVG$e_desc = "PR % Returns"
TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  

