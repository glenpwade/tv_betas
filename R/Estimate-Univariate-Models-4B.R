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
ptitle = "FourB"
allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_fourB <- allData$FourB
#

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_fourB
e <- e_fourB[1:2200]
Tobs = NROW(e)
ptitle = "FourB stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_fourB[1:2200]
# Method:  MLE 
# Est      se1 sig
# omega 0.029398 0.009632 ***
#     alpha 0.081690 0.016093 ***
#     beta  0.893205 0.021847 ***
#     
#     Log-likelihood value(GARCH):  -3144.83

# No suitable stable subset, so...
# Try the rolling window method:
if(FALSE){
    e <- e_fourB[1:2200]
    estCtrl$vartargetWindow = 300
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    
    # Method:  MLE, variance-targetting a rolling Window of 300 observations 
    # Est      se1 sig
    # omega 0.118569       NA    
    # alpha 0.117745 0.024330 ***
    #     beta  0.777439 0.050894 ***
    #     
    #     Log-likelihood value(GARCH):  -3163.089
}


# Next, We need a standard TV object to generate the data:
e <- e_fourB
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_FourB <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "FourB"
saveRDS(refData_FourB,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_fourB
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_FourB

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

## P-Values from TEST Results, TV-delta0 only, FourB[1:3152]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.147|  0.065|
#     |        2| 0.283|  0.182|
#     |        3| 0.011|  0.052|

## Conclusion: 
## Evidence of a bump - at least 3 - consistent with plot (which suggests 4 or 5)
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
TV$optimcontrol$parscale = c(5,3,1,1,30,10,1,30,10,1)
TV$pars["deltaN",] = c(-1,18,-17)
TV$pars["speedN",] = c(2,6,5)
TV$pars["locN1",] = c(0.4,0.7,0.75)
## The above estimates the model well, but Tests won't run...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
plot(TV)

#  Log-likelihood value(TV):  -4863.825

## P-Values from TEST Results, TV-3Trans only, FourB[1:3152]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.033|  0.000|
#     |        2| 0.051|  0.001|
#     |        3| 0.063|  0.001|

## Conclusion: There is at least one more Transition in there


## 4-TRANS ####

## Let's try estimating a model with 4 x single order transitions...
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


## 5-TRANS ####

## Let's try estimating a model with 5 x single order transitions...
## Default starting values don't work well for this 5-Trans model, so...
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single,tvshape$single))

TV$optimcontrol$reltol = 1e-05
TV$optimcontrol$ndeps = rep(1e-07,length(TV$optimcontrol$ndeps))
TV$delta0 = 3
TV$pars["deltaN",] = c(-2,2,-2,21,-20)
TV$pars["speedN",] = c(4,4,4,6,5)
TV$pars["locN1",] = c(0.2,0.3,0.4,0.7,0.75)
## The above estimates the model well, but Tests won't run...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
plot(TV)

saveRDS(TV,paste0('Results/',ptitle,'_Final_TV_model.RDS'))


## Let's try estimating a model with 6 x single order transitions...
## Default starting values don't work well for this 6-Trans model, so...
TV <- tv(st,c(tvshape$single,tvshape$double1loc,tvshape$single,tvshape$single,tvshape$single,tvshape$single))

TV$optimcontrol$reltol = 1e-05
TV$optimcontrol$ndeps = rep(1e-07,length(TV$optimcontrol$ndeps))
TV$delta0 = 3
TV$pars["deltaN",] = c(-2,1,2,-2,21,-20)
TV$pars["speedN",] = c(4,4,4,4,6,5)
TV$pars["locN1",] = c(0.1,0.2,0.3,0.4,0.7,0.75)
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

# Delta0 = 2.527624    se0 = 0.342116*** 
#     
#     st1      se1 sig      st2      se2 sig.1       st3      se3 sig.2       st4      se4 sig.3        st5      se5 sig.4
# deltaN -1.777747 0.338554 *** 1.660221 0.611926   *** -1.492315 0.613621   **  23.585797 0.332149   *** -23.267901 0.325150   ***
#     speedN  6.585951 1.180486 *** 5.312794 0.356804   ***  6.811791 2.065288   ***  6.469602 0.318602   ***   4.908007 0.116052   ***
#     locN1   0.054008 0.002174 *** 0.326807 0.004908   ***  0.418146 0.003847   ***  0.711973 0.000828   ***   0.729068 0.002596   ***
#     locN2         NA      NaN           NA      NaN              NA      NaN              NA      NaN               NA      NaN      
# 
# Log-likelihood value(TV):  -4743.887

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full FourB data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

e <- e_fourB
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "FourB"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "FourB % Returns"

TVG$tvOptimcontrol$reltol = 1e-05
TVG$garchOptimcontrol$reltol = 1e-04

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG@garchObj)
summary(TVG@tvObj)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model_1.RDS'))
plot(e,type='l')
lines(sqrt(TVG$Estimated$g),col="red")

# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

plot(TV,main=ptitle)
summary(TVG1)




