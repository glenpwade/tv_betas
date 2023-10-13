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
e_crau <- tail(as.numeric(prices$CR_AU),-1)  # Percentage Returns - drop first (null) observation

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_crau
e <- e_crau[1385:2250]
Tobs = NROW(e)
ptitle = "CRAU stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# Next, We need a standard TV object to generate the data:
e <- e_crau
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_CRAU <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist, seed=1)
ptitle = "CRAU"
saveRDS(refData_CRAU,paste0('RefData/',ptitle,'_RefData.RDS'))

#refData = readRDS(paste0('RefData/',ptitle,'_RefData.RDS'))
##  End of Ref Data Generation ----


e <- e_crau
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_CRAU

# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

## RESTART HERE with an updated TV Model specification:   ####

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

## P-Values from TEST Results, TV-delta0 only, CRAU[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.000|      0|
# |        2| 0.001|      0|
# |        3| 0.003|      0|

## Conclusion: 
## Lots of Evidence of a bump
## Let's try estimating a TV_1Trans model

e <- e_crau
Tobs = NROW(e)
refData = refData_CRAU
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

#Log-likelihood value(TV):  -5026.424

## Now Let's try estimating a TV_2Trans model
TV <- tv(st,c(tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(0.5,1)
TV$pars["speedN",] = c(4,3)
TV$pars["locN1",] = c(0.3,0.7)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
#Log-likelihood value(TV):  -4921.398

## Now Let's try estimating a TV_3Trans model
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(0.5,-1,1)
TV$pars["speedN",] = c(4,3,5)
TV$pars["locN1",] = c(0.2,0.4,0.7)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
# Not easy to estimate - let's test for more than 2 transitions

## P-Values from TEST Results, TV-2Trans CRAU[1:3153]:

# | Test Ord| TR2| Robust|
# |--------:|---:|------:|
# |        1|   0|      0|
# |        2|   0|      0|
# |        3|   0|      0|

## Conclusion: 
## Strong Evidence of another transition!


## Let's try to estimate a TV-3_Trans model on the full dataset:

e <- e_crau
Tobs = NROW(e)
refData = refData_CRAU
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single, tvshape$double1loc))
TV$delta0 = 3
TV$pars["deltaN",] = c(1,-0.5,-1)
TV$pars["speedN",] = c(4,4,4)
TV$pars["locN1",] = c(0.25,0.4,0.75)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")
# The above specification works, but the last transition may be a double...


# Shape guessed by observing the plot:
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single, tvshape$single))
TV$delta0 = 1
TV$pars["deltaN",] = c(2,-2,2)
TV$pars["speedN",] = c(2,2,2)
TV$pars["locN1",] = c(0.25,0.35,0.65)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

# Shape guessed by observing the plot:
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single, tvshape$single))
TV$delta0 = 2
TV$pars["deltaN",] = c(-1,2,-1,1)
TV$pars["speedN",] = c(3,3,3,4)
TV$pars["locN1",] = c(0.05,0.35,0.45,0.65)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-4Trans, CRAU[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.058|  0.168|
# |        2| 0.046|  0.214|
# |        3| 0.040|  0.078|

## Conclusion: 
## Some Evidence of another transition, based on TR2 test...
## Given difficulties estimating this we will stop here

## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 1 1 1 
# 
# Estimation Results:
#     
# Delta0 = 1.210009    se0 = 0.244214*** 
#     
#              st1      se1 sig      st2      se2 sig.1       st3      se3 sig.2      st4      se4 sig.3
# deltaN -1.639075 0.627587 *** 3.334272 0.378988   *** -2.108204 0.906644   **  1.742438 0.134088   ***
# speedN  2.655870 0.286956 *** 4.545300 0.331021   ***  4.389454 0.627183   *** 6.999988 0.000044   ***
# locN1   0.201110 0.082037 **  0.290866 0.003387   ***  0.423108 0.028230   *** 0.710160 0.001318   ***
# locN2         NA      NaN           NA      NaN              NA      NaN             NA      NaN      
# 
# Log-likelihood value(TV):  -4786.416
