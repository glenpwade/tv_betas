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
ptitle = "STW"
allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_stw <- allData$STW
#

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_stw
e <- e_stw[1320:2225]
Tobs = NROW(e)
ptitle = "STW stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# Method:  MLE 
# Est      se1 sig
# omega 0.028801 0.012700 ** 
# alpha 0.093439 0.025424 ***
# beta  0.849340 0.042595 ***
#     
#     Log-likelihood value(GARCH):  -925.0433

# Try the rolling window method:
if(FALSE){
    e <- e_stw
    estCtrl$vartargetWindow = 300
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn

    # Method:  MLE, variance-targetting a rolling Window of 300 observations 
    # Est      se1 sig
    # omega 0.045592       NA    
    # alpha 0.102220 0.014145 ***
    # beta  0.827896 0.028061 ***
    #     
    #     Log-likelihood value(GARCH):  -3955.634
}

# Next, We need a standard TV object to generate the data:
e <- e_stw
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_STW <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "STW"
saveRDS(refData_STW,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_stw
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_STW

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

## P-Values from TEST Results, TV-delta0 only, STW[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.872|  0.856|
# |        2| 0.904|  0.660|
# |        3| 0.068|  0.069|

## Conclusion: 
## No Evidence of a bump - but the plot suggests otherwise: 3 or 4
## We can try to identify & estimate a single transition, or split the data and retest:
## Let's try estimating a TV-3Trans model...

TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
lines(sqrt(TV@g),col="red")
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

# Transition Shapes: 1 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 0.955045    se0 = 0.039536*** 
#     
#     st1      se1 sig      st2      se2 sig.1       st3      se3 sig.2
#       deltaN -0.493840 0.045876 *** 3.862118 0.463357   *** -3.619124 0.461356   ***
#       speedN  4.206573 0.316292 *** 6.999630 0.403541   ***  4.908192 0.203547   ***
#       locN1   0.427292 0.010621 *** 0.711101 0.000650   ***  0.736789 0.004028   ***
#      
# Log-likelihood value(TV):  -4047.496

# P-Values from TEST Results, TV-3Trans only, STW[1:3152]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.819|  0.862|
#     |        2| 0.827|  0.550|
#     |        3| 0.724|  0.156|

## Conclusion: 
## No Evidence of another bump - but the plot suggests otherwise
## We could try alternate specifications...


# By observing the plot, we can try the following shape:
TV <- tv(st,c(tvshape$single,tvshape$double1loc,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-05
TV$optimcontrol$ndeps = rep(1e-05,length(TV$optimcontrol$ndeps))
TV$optimcontrol$parscale = c(3,1,6,1,1,6,1,12,6,1,10,6,1)
TV$delta0 = 3
TV$pars["deltaN",] = c(-0.9,-1,12,-11)
TV$pars["speedN",] = c(5,5,6,5)
TV$pars["locN1",] = c(0.05,0.3,0.7,0.73)
#
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
plot(TV)

## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 3 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 3.514712    se0 = 0.312669*** 
#     
#     st1      se1 sig       st2      se2 sig.1       st3      se3 sig.2        st4      se4 sig.3
# deltaN -1.250001 0.162690 *** -1.752758 0.276607   *** 14.523741 3.714487   *** -14.300991 3.710504   ***
#     speedN  6.439409 0.717728 ***  6.470635 0.287741   ***  6.999987      NaN         5.566730 0.136536   ***
#     locN1   0.058050 0.001513 ***  0.360354 0.005211   ***  0.712740 0.000446   ***   0.723508 0.002155   ***
#     locN2         NA      NaN            NA      NaN              NA      NaN               NA      NaN      
# 
# Log-likelihood value(TV):  -3953.356


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full STW data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

e <- e_stw
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "STW"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "STW % Returns"

TVG$tvpars[,1] = c(-1.5,5,0.1,NA)
# TVG$tvpars[,2] = c(-0.7,2.5,0.25,NA)
# TVG$tvOptimcontrol$reltol = 1e-04
# TVG$tvOptimcontrol$ndeps = rep(1e-06,length(TVG$tvOptimcontrol$ndeps))
TVG$garchpars[,1] = c(0.02,0.002,0.85,0.1)
TVG$garchOptimcontrol$reltol = 1e-04
TVG$garchOptimcontrol$ndeps = c(1e-04,1e-07,1e-04,1e-04)
TVG$garchOptimcontrol$parscale = c(50,1,300,65)

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




