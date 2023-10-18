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
e_wbc <- diff(log(as.numeric(prices$WBC)) ) * 100  # Percentage Returns
#
# Replace all NA's with the previous valid entry if required:
if(FALSE){
    na_pos = which(is.na(e_wbc) )
    for(n in seq_along(na_pos)){
        t_pos <- na_pos[n]
        e_wbc[t_pos] <- e_wbc[t_pos-1]
    }
}

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_wbc[1320:2225]
Tobs = NROW(e)
ptitle = "WBC stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_wbc[1400:1800]
# Est se1 sig
# omega 0.072941 NaN    
# alpha 0.070489 NaN    
# beta  0.852360 NaN 

# e_wbc[2500:3150]
# Est      se1 sig
# omega 0.132530 0.088212    
# alpha 0.061153 0.027582 ** 
# beta  0.839795 0.083945 ***


# Next, We need a standard TV object to generate the data:
e <- e_wbc
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_WBC <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist, seed=1)
ptitle = "WBC"
saveRDS(refData_WBC,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_wbc
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
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

## P-Values from TEST Results, TV-delta0 only, WBC[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.085|  0.035|
# |        2| 0.211|  0.114|
# |        3| 0.000|  0.003|

## Conclusion: 
## Evidence of a bump - first or third order.
## We can try to identify & estimate a single transition, or split the data and retest:
## Let's try estimating a single order transition...

e <- e_wbc
Tobs = NROW(e)
#plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_WBC
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$single)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-1_Trans, WBC[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.004|  0.014|
# |        2| 0.003|  0.012|
# |        3| 0.001|  0.001|

## Conclusion: 
## Evidence of another transition exits
## Let's try estimating a model with 2 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-7
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-2_Trans, WBC[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.144|  0.062|
# |        2| 0.008|  0.013|
# |        3| 0.002|  0.005|


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
TV$optimcontrol$reltol = 1e-04
TV$optimcontrol$ndeps = rep(1e-06,10)
TV$pars["deltaN",] = c(-1,8,-7)
TV$pars["speedN",] = c(2,6,5)
TV$pars["locN1",] = c(0.45,0.7,0.75)
## The above estimates the model well, but Tests won't run...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

# So - let's tweak the Optim controls:
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-04
TV$optimcontrol$ndeps = rep(1e-06,10)
TV$optimcontrol$parscale = c(1,3,3,1,15,7,1,15,5,1)
TV$pars["deltaN",] = c(-1,8,-7)
TV$pars["speedN",] = c(2,6,5)
TV$pars["locN1",] = c(0.15,0.7,0.75)
#
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
plot(e_wbc,type='l')

# Note: Migh need to revisit WBC univar TV estimate...

## P-Values from TEST Results, TV-3_Trans, WBC[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.078|  0.057|
# |        2| 0.086|  0.116|
# |        3| 0.097|  0.206|

#

## Conclusion: 
## No Evidence of another transition!


## Final Model Specification:  ----
# 
# TV OBJECT
# 
# Transition Shapes: 1 1 1 
# 
# Estimation Results:
#     
# Delta0 = 1.540203    se0 = 0.047555*** 
#     
#             st1      se1 sig      st2      se2 sig.1       st3      se3 sig.2
# deltaN 4.725349 0.636648 *** 5.019944 1.925816   *** -9.778395 2.245770   ***
# speedN 4.293456 0.403468 *** 4.927609 0.591673   ***  4.913416 0.245546   ***
# locN1  0.706207 0.008651 *** 0.706462 0.004997   ***  0.749147 0.006449   ***
# locN2        NA      NaN           NA      NaN              NA      NaN      
# 
# Log-likelihood value(TV):  -5331.788


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full NAB data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

prices <- read.csv("data/tv_betas_prices.csv")
e_wbc <- diff(log(as.numeric(prices$WBC)) ) * 100  # Percentage Returns
e <- e_wbc
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "WBC"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "WBC Std.% Returns"

TVG$tvpars[,1] = c(-0.5,3,0.01,NA)
TVG$tvOptimcontrol$reltol = 1e-05
TVG$garchOptimcontrol$reltol = 1e-04

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

plot(TV,main=ptitle)

summary(TVG1)

