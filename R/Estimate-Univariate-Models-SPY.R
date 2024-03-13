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
ptitle = "SPY"
allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_spy <- allData$SPY

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_spy
e <- e_spy[180:1099]
Tobs = NROW(e)
ptitle = "SPY stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# Method:  MLE 
# Est      se1 sig
# omega 0.078414 0.023503 ***
# alpha 0.117832 0.030850 ***
# beta  0.743659 0.057963 ***
#     
#     Log-likelihood value(GARCH):  -1011.895

# Try the rolling window method:
if(FALSE){
    e <- e_spy[200:2200]
    Tobs = NROW(e)
    estCtrl$vartargetWindow = 500
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    
    # Method:  MLE, variance-targetting a rolling Window of 500 observations 
    #            Est      se1 sig
    # omega 0.089715       NA    
    # alpha 0.173298 0.021142 ***
    # beta  0.726568 0.035248 ***
    #     
    #     Log-likelihood value(GARCH):  -2215.578
}

# Next, We need a standard TV object to generate the data:
e <- e_spy
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_SPY <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "SPY"
saveRDS(refData_SPY,paste0('RefData/',ptitle,'_RefData.RDS'))

#refData = readRDS(paste0('RefData/',ptitle,'_RefData.RDS'))
#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_spy
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_SPY

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

## P-Values from TEST Results, TV-delta0 only, SPY[1:3152]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.150|  0.106|
#     |        2| 0.218|  0.001|
#     |        3| 0.020|  0.000|

## Conclusion: 
## Strong Evidence of a bump, most likely 3
## Let's try estimating a 3 transition model...

TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
## Default starting values don't work perfectly for this 3-Trans model, so...
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-03
TV$optimcontrol$ndeps = rep(1e-06,length(TV$optimcontrol$ndeps))
TV$pars["deltaN",] = c(-1,1,5,-4)
TV$pars["speedN",] = c(4,4,6,5)
TV$pars["locN1",] = c(0.15,0.4,0.7,0.75)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
# A TV-3_Trans model fitted well, even with default starting params

## P-Values from TEST Results, TV-3_Trans, SPY[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.163|  0.023|
# |        2| 0.166|  0.032|
# |        3| 0.152|  0.050|


## Conclusion: 
## Still Evidence of another transition exits
## Let's try estimating a model based on the plot

# By observing the plot, we can try the following shape:
e <- e_spy
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$double1loc,tvshape$single,tvshape$single,tvshape$double1loc))
TV$delta0 = 4
TV$pars["deltaN",] = c(1,9,-7,-2)
TV$pars["speedN",] = c(5,6,5,6)
TV$pars["locN1",] = c(0.5,0.7,0.75,0.9)
#
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
lines(sqrt(TV@g),col="red")
plot(TV)
abline(v=seq(0,3200,by=200),col="grey80")

## P-Values from TEST Results, TV-4_Trans, SPY[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.760|  0.404|
# |        2| 0.705|  0.564|
# |        3| 0.382|  0.249|

## Conclusion: 
## No Evidence of another transition!


## Final Model Specification:  ----

# Transition Shapes: 3 1 1 3 
# 
# Estimation Results:
#     
#     Delta0 = 4.110256    se0 = 0.678409*** 
#     
#     st1      se1 sig       st2      se2 sig.1        st3      se3 sig.2       st4      se4 sig.3
#     deltaN 1.569631 0.078639 *** 13.278953 4.219621   *** -13.634410 4.135510   *** -4.730295 0.677332   ***
#     speedN 5.819993 0.139181 ***  6.263004 0.224640   ***   4.649898 0.296112   ***  6.954825 0.185826   ***
#     locN1  0.482457 0.003236 ***  0.710209 0.000933   ***   0.726438 0.005047   ***  0.904094 0.002608   ***
#     locN2        NA      NaN            NA      NaN               NA      NaN              NA      NaN      
# 
# Log-likelihood value(TV):  -4258.319


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full SPY data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

e <- e_spy
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "SPY"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "SPY % Returns"

#TVG$garchpars[,1] = c(0.02,0.002,0.8,0.15)
#TVG$garchOptimcontrol$reltol = 1e-03
TVG$garchOptimcontrol$ndeps = c(1e-07,1e-03,1e-07,1e-07)
TVG$garchOptimcontrol$parscale = c(10,1,50,10)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG@garchObj)
summary(TVG@tvObj)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
plot(TV,main=ptitle)
summary(TVG1)
plot(TVG1)




