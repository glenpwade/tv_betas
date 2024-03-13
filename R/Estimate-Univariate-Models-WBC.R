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
ptitle = "WBC"
allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_wbc <- allData$WBC
#

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_wbc
e <- e_wbc[1320:2225]
Tobs = NROW(e)
ptitle = "WBC stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_wbc
# Method:  MLE 
# Est      se1 sig
# omega 0.103720 0.021855 ***
#     alpha 0.140434 0.020762 ***
#     beta  0.806624 0.027887 ***
#     
#     Log-likelihood value(GARCH):  -5150.218


# No suitable stable subset, so...
# Try the rolling window method:
if(FALSE){
    e <- e_wbc
    estCtrl$vartargetWindow = 500
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    
    #     Method:  MLE, variance-targetting a rolling Window of 500 observations 
    # Est      se1 sig
    # omega 0.100310       NA    
    # alpha 0.151057 0.019162 ***
    # beta  0.787219 0.030467 ***
    #     
    #     Log-likelihood value(GARCH):  -5171.531
}


# Next, We need a standard TV object to generate the data:
e <- e_wbc
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_WBC <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "WBC"
saveRDS(refData_WBC,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_wbc
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_WBC

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

## P-Values from TEST Results, TV-delta0 only, WBC[1:3152]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.373|  0.284|
#     |        2| 0.653|  0.572|
#     |        3| 0.123|  0.206|

## Conclusion: 
## No Evidence of a bump - but many visible in the plot
## Let's try estimating a model with 3 x single order transitions...
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## Default starting values don't work well for this 3-Trans model, so...
## Looking at the plot, the missing transition seems to be high-to-low around Obs 2500
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-05
TV$optimcontrol$ndeps = rep(1e-05,10)
#TV$optimcontrol$parscale = c(5,3,1,1,30,10,1,30,10,1)
TV$pars["deltaN",] = c(-1,18,-17)
TV$pars["speedN",] = c(2,6,5)
TV$pars["locN1",] = c(0.25,0.7,0.75)
## The above estimates the model well, but Tests won't run...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
plot(TV)

## Log-likelihood value(TV):  -5281.867


## 4-TRANS ####

## Let's try estimating a model with 3 x single order transitions...
## Default starting values don't work well for this 4-Trans model, so...
## Looking at the plot, the missing transition seems to be high-to-low around Obs 2500
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-04
TV$optimcontrol$ndeps = rep(1e-07,length(TV$optimcontrol$ndeps))
TV$optimcontrol$parscale = c(5,3,1,1,30,10,1,5,3,1,30,10,1)
TV$pars["deltaN",] = c(1,-1,18,-17)
TV$pars["speedN",] = c(3,3,6,5)
TV$pars["locN1",] = c(0.1,0.25,0.7,0.75)
## The above estimates the model well, but Tests won't run...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
plot(TV)

# Log-likelihood value(TV):  -5265.986

## 5-TRANS ####

## Let's try estimating a model with 5 x single order transitions...
## Default starting values don't work well for this 5-Trans model, so...
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-05
TV$optimcontrol$ndeps = rep(1e-07,length(TV$optimcontrol$ndeps))
TV$optimcontrol$parscale = c(5,3,1,1,5,3,1,5,3,1,30,10,1,30,10,1)
TV$pars["deltaN",] = c(-1,2,-2,16,-15)
TV$pars["speedN",] = c(3,4,4,6,5)
TV$pars["locN1",] = c(0.2,0.3,0.4,0.7,0.75)
## The above estimates the model well, but Tests won't run...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
plot(TV)

saveRDS(TV,paste0('Results/',ptitle,'_Final_TV_model.RDS'))


## Final Model Specification:  ----
# 
# TV OBJECT
# 
# Transition Shapes: 1 1 1 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 2.579733    se0 = 0.293072*** 
#     
#     st1      se1 sig      st2      se2 sig.1       st3      se3 sig.2       st4      se4 sig.3        st5      se5 sig.4
#       deltaN -2.927694 0.466200 *** 3.501768 0.619137   *** -2.007502 0.568469   *** 17.305049 2.140172   *** -16.947159 2.077289   ***
#       speedN  2.743619 0.296067 *** 3.723701 0.264668   ***  6.999675 0.402992   ***  6.998945 0.628456   ***   4.629745 0.115255   ***
#       locN1   0.188679 0.014303 *** 0.320066 0.016233   ***  0.418380 0.002156   ***  0.712607 0.000532   ***   0.735961 0.002909   ***
# 
# Log-likelihood value(TV):  -5177.808

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full WBC data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

e <- e_wbc
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "WBC"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "WBC % Returns"

TVG$tvOptimcontrol$reltol = 1e-05
TVG$garchOptimcontrol$reltol = 1e-04

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG@garchObj)
summary(TVG@tvObj)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
plot(e,type='l')
lines(sqrt(TVG$Estimated$g),col="red")

# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

plot(TV,main=ptitle)
summary(TVG1)




