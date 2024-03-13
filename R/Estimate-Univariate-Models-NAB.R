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

ptitle = "NAB"
allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_nab <- allData$NAB
#

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_nab[1350:2044]
Tobs = NROW(e)
ptitle = "NAB stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_nab[1350:2044]
# Method:  MLE 
#                  Est      se1 sig
#       omega 0.123536 0.054356 ** 
#       alpha 0.121999 0.037354 ***
#       beta  0.732413 0.086044 ***

# Try the rolling window method:
if(FALSE){
    e <- e_nab[1:2200]
    estCtrl$vartargetWindow = 500
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn

    # Method:  MLE, variance-targetting a rolling Window of 500 observations 
    #            Est      se1 sig
    # omega 0.047447       NA    
    # alpha 0.086569 0.016788 ***
    # beta  0.862922 0.028792 ***
}

# Next, We need a standard TV object to generate the data:
e <- e_nab
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_NAB <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "NAB"
saveRDS(refData_NAB,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_nab
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_NAB

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

## P-Values from TEST Results, TV-delta0 only, NAB[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.339|  0.258|
# |        2| 0.639|  0.322|
# |        3| 0.053|  0.061|

## Conclusion: 
## Some Evidence of a bump - third order.
## We can try to split the data and retest, but transitions areclearly visible in the chart
## Let's try estimating a 3 transition model...


TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## Default starting values don't work great for this 3-Trans model, so...
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(1,8,-7)
TV$pars["speedN",] = c(4,6,6)
TV$pars["locN1",] = c(0.2,0.6,0.8)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(e,type='l')
lines(TV@g,col="red")
lines(sqrt(TV@g),col="red")
plot(TV)

## P-Values from TEST Results, TV-3_Trans, NAB[1:3152]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.967|  0.941|
# |        2| 0.758|  0.415|
# |        3| 0.550|  0.373|

#

## Conclusion: 
## No Evidence of another transition!


## Final TV Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 1.656822    se0 = 0.065028*** 
#     
#     st1      se1 sig       st2      se2 sig.1        st3      se3 sig.2
# deltaN -0.717305 0.080684 *** 21.188481 3.891130   *** -20.696662 3.875164   ***
#     speedN  5.769810 2.370178 **   6.622214 0.217393   ***   4.701197 0.116965   ***
#     locN1   0.421499 0.011489 ***  0.711615 0.000784   ***   0.734605 0.003236   ***
#     locN2         NA      NaN            NA      NaN               NA      NaN      
# 
# Log-likelihood value(TV):  -5084.701



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full NAB data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

e <- e_nab
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "NAB"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "NAB % Returns"

TVG$tvpars[,1] = c(-0.4,6,0.1,NA)
TVG$tvpars[,2] = c(8,6,0.6,NA)
TVG$tvpars[,3] = c(-7,6,0.75,NA)
TVG$tvOptimcontrol$reltol = 1e-07
TVG$tvOptimcontrol$ndeps = rep(1e-06,length(TVG$tvOptimcontrol$ndeps))
TVG$tvOptimcontrol$parscale = c(5,70,1,90,70,7,80,60,7)
TVG$garchpars[,1] = c(0.1,0.01,0.7,0.07)
#TVG$garchOptimcontrol$reltol = 1e-04

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

