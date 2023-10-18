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
e_anz <- diff(log(as.numeric(prices$ANZ)) ) * 100  # Percentage Returns
#
# Replace all NA's with the previous valid entry if required:
if(FALSE){
    e_anz_na_pos = which(is.na(e_anz) )
    for(n in seq_along(e_anz_na_pos)){
        t_pos <- e_anz_na_pos[n]
        e_anz[t_pos] <- e_anz[t_pos-1]
    }
}

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_anz
e <- e_anz[2500:3153]
Tobs = NROW(e)
ptitle = "ANZ stable subest"
plot(e,type='l',main=ptitle)
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_anz[1400:1800]
# Est se1 sig
# omega 0.072941 NaN    
# alpha 0.070489 NaN    
# beta  0.852360 NaN 

# e_anz[2500:3150]
# Est      se1 sig
# omega 0.132530 0.088212    
# alpha 0.061153 0.027582 ** 
# beta  0.839795 0.083945 ***

# Next, We need a standard TV object to generate the data:
e <- e_anz
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_ANZ <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist, seed=1)
ptitle = "ANZ"
saveRDS(refData_ANZ,paste0('RefData/',ptitle,'_RefData.RDS'))

# Reload refdata:
refData = readRDS(paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation  #

# START - Test for any transition ----

e <- e_anz
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
refData = refData_ANZ

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
#saveRDS(SIMRESULT,saveFile)

## RELOAD Results
#SIMRESULT = readRDS(saveFile)

## RESULTS Section ----

## P-Values from TEST Results, TV-delta0 only, ANZ[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.048|  0.008|
# |        2| 0.078|  0.039|
# |        3| 0.002|  0.018|

## Conclusion: 
## Evidence of a bump first and/or third order.
## We can try to identify & estimate a single transition, or split the data and retest:
## Let's try estimating a single order transition...

e <- e_anz
Tobs = NROW(e)
#plot(e,type='l',main=ptitle)
refData = refData_ANZ
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$single)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-1_Trans, ANZ[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.019|  0.100|
# |        2| 0.005|  0.020|
# |        3| 0.000|  0.005|

## Conclusion: 
## Evidence of another transition exits
## Let's try estimating a model with 2 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-2_Trans, ANZ[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.105|  0.683|
# |        2| 0.029|  0.093|
# |        3| 0.001|  0.000|

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
TV$pars["deltaN",] = c(-0.5,2,-1)
TV$pars["speedN",] = c(4,6,6)
TV$pars["locN1",] = c(0.4,0.6,0.8)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-3_Trans, ANZ[1:3153]:

# | Test Ord| TR2| Robust|
# |--------:|---:|------:|
# |        1|   0|  0.144|
# |        2|   0|  0.176|
# |        3|   0|  0.001|
#
# Log-likelihood value(TV):  -5287.732

## Conclusion: 
## Still Evidence of another transition exits
## Let's try estimating a model with 4 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## Default starting values don't work well for this 4-Trans model, so...
## Looking at the plot, the missing transition seems to be around Obs 1000 
e <- e_anz
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(1,-1,6,-5)
TV$pars["speedN",] = c(2,2,5,4)
TV$pars["locN1",] = c(0.35,0.4,0.7,0.75)
TV$optimcontrol$reltol = 1e-9
TV$optimcontrol$ndeps = rep(1e-6,13)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
#
## P-Values from TEST Results, TV-4_Trans, ANZ[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.437|  0.216|
# |        2| 0.182|  0.023|
# |        3| 0.194|  0.050|
    
## Conclusion: 
## No Evidence of another transition exits!
## OK, so maybe a very slight chance with Test Ord 3 being 5%, but...
## it was difficult to find starting pars & optim controls to estimate this model
## so we will stop here.

## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 1 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 1.501141    se0 = 0.164958*** 
#     
#             st1      se1           st2      se2             st3      se3              st4      se4 
# deltaN 1.945325 0.434912 *** -2.432078 0.395698   *** 18.181107 2.594712   *** -17.892599 2.585121   ***
# speedN 4.788028 0.717099 ***  2.388774 0.362497   ***  6.999997 0.000000   ***   4.625933 0.104238   ***
# locN1  0.313708 0.010694 ***  0.342392 0.031867   ***  0.712862 0.000429   ***   0.738316 0.002317   ***
# locN2        NA      NaN            NA      NaN              NA      NaN               NA      NaN      
# 
# Log-likelihood value(TV):  -5130.843



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full ANZ data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.
ptitle = "ANZ"
prices <- read.csv("data/tv_betas_prices.csv")
e_anz <- diff(log(as.numeric(prices$ANZ)) ) * 100  # Percentage Returns
e <- e_anz
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
# TV = get final specification from code above
#
TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "ANZ Std.% Returns"

TVG$tvpars["speedN",1] = 4
TVG$tvOptimcontrol$reltol = 1e-07
TVG$tvOptimcontrol$ndeps = rep(1e-05,length(TVG$tvOptimcontrol$ndeps))
TVG$garchpars[,1] = c(0.1,0.02,0.7,0.05)
TVG$garchOptimcontrol$reltol = 1e-04
TVG$garchOptimcontrol$reltol = 1e-05
TVG$garchOptimcontrol$parscale = c(4,1,80,10)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))


