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
ptitle = "CRAU"
allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_crau <- allData$CRAU
# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_crau
e <- e_crau[1430:2240]
Tobs = NROW(e)
ptitle = "CRAU stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# Method:  MLE 
#            Est      se1 sig
# omega 0.043531 0.061400    
# alpha 0.008945 0.011538    
# beta  0.935761 0.081705 ***
#     
#     Log-likelihood value(GARCH):  -1059.747


# Try the rolling window method:
if(FALSE){
    e <- e_crau[1:2220]
    estCtrl$vartargetWindow = 500
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    
    # Method:  MLE, variance-targetting a rolling Window of 500 observations 
    #                  Est      se1 sig
    #       omega 0.022534       NA    
    #       alpha 0.066584 0.014676 ***
    #       beta  0.907564 0.023719 ***
    #     
    #     Log-likelihood value(GARCH):  -3046.391
    
}

# Next, We need a standard TV object to generate the data:
e <- e_crau
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_CRAU <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "CRAU"
saveRDS(refData_CRAU,paste0('RefData/',ptitle,'_RefData.RDS'))

#refData = readRDS(paste0('RefData/',ptitle,'_RefData.RDS'))
#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_crau
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_CRAU

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

## P-Values from TEST Results, TV-delta0 only, CRAU[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.139|  0.080|
# |        2| 0.323|  0.197|
# |        3| 0.458|  0.258|

## Conclusion: 
## No Evidence of a bump, despite clear visual indicator from chart
## Let's try testing a subset

e <- e_crau[1:1400]
Tobs = NROW(e)
refData = refData_CRAU
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$delta0only))
TV <- estimateTV(e,TV,estCtrl)

## P-Values from TEST Results, TV-delta0 only, CRAU[1:1400]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.062|  0.098|
#     |        2| 0.041|  0.119|
#     |        3| 0.062|  0.007|

## Conclusion: 
## Evidence of 2 or 3 transitions in this subset, and at least one more in remaining sample
## Let's try to estimate a TV-3_Trans model on the full dataset:


## Now Let's try estimating a TV_3Trans model
e <- e_crau
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
#
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(1.5,-1.4,1.6)
TV$pars["speedN",] = c(5,4,5)
TV$pars["locN1",] = c(0.2,0.4,0.7)
TV$optimcontrol$reltol = 1e-5
TV$optimcontrol$ndeps = rep(1e-7,10)
TV$optimcontrol$parscale = c(2,20,60,1,20,60,2,20,60,2)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-3Trans CRAU[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.279|  0.992|
# |        2| 0.402|  0.974|
# |        3| 0.220|  0.171|

## Conclusion: 
## No Evidence of another transition!

## Final TV Model Specification:  ----
# TV OBJECT
# 
# Transition Shapes: 1 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 0.864655    se0 = 0.037235*** 
#     
#     st1      se1 sig       st2      se2 sig.1      st3      se3 sig.2
# deltaN 1.629320 0.240167 *** -1.698121 0.241103   *** 1.579615 0.117380   ***
#     speedN 6.343883 1.889447 ***  4.095089 0.251514   *** 6.999991 0.569397   ***
#     locN1  0.345040 0.002051 ***  0.441625 0.008531   *** 0.711075 0.001301   ***
#     locN2        NA      NaN            NA      NaN             NA      NaN      
# 
# Log-likelihood value(TV):  -4878.357



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full CR_AU data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

e <- e_crau
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "CRAU"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "CRAU % Returns"

TVG$tvpars["speedN",] = c(2,5,5,5)
TVG$tvpars["locN1",] = c(0.1,0.3,0.45,0.63)
TVG$tvOptimcontrol$reltol = 5e-04
TVG$garchpars[,1] = c(0.1,0.01,0.8,0.02)
TVG$garchOptimcontrol$reltol = 5e-05
TVG$garchOptimcontrol$ndeps = c(1e-04,1e-09,1e-04,1e-04)
TVG$garchOptimcontrol$parscale = c(50,1,140,7)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG@garchObj)
summary(TVG@tvObj)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
summary(TVG1@garchObj)
summary(TVG1@tvObj)
plot(TVG1)



