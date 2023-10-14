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
e_pr <- tail(as.numeric(prices$PR),-1)  # Percentage Returns - drop first (null) observation

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_pr
e <- e_pr[2300:3153]
Tobs = NROW(e)
ptitle = "PR stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# Next, We need a standard TV object to generate the data:
e <- e_pr
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_PR <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist, seed=1)
ptitle = "PR"
saveRDS(refData_PR,paste0('RefData/',ptitle,'_RefData.RDS'))

#refData = readRDS(paste0('RefData/',ptitle,'_RefData.RDS'))
#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_pr
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

## P-Values from TEST Results, TV-delta0 only, PR[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.846|  0.842|
# |        2| 0.955|  0.957|
# |        3| 0.887|  0.910|

## Conclusion: 
## No Evidence of a bump, though the plot suggests otherwise
## Let's try splitting the series to detect transitions...

e <- e_pr[300:2200]
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_PR
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-delta0 only,, PR[300:2200]

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.532|  0.440|
# |        2| 0.082|  0.047|
# |        3| 0.180|  0.055|

## Conclusion: 
## Some Evidence of a transition exits 
## Let's try estimating a model with 1 x single order transitions...

e <- e_pr[300:2200]
Tobs = NROW(e)
refData = refData_PR
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$single)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-1_Trans, PR[300:2200]

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.009|  0.000|
# |        2| 0.028|  0.008|
# |        3| 0.050|  0.013|

## Conclusion: 
## Still Evidence of another transition exits
## Let's try estimating a model with 2 x single order transitions...

e <- e_pr[300:2200]
Tobs = NROW(e)
refData = refData_PR
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(1,-1)
TV$pars["speedN",] = c(4,4)
TV$pars["locN1",] = c(0.3,0.7)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.548|  0.416|
# |        2| 0.622|  0.698|
# |        3| 0.598|  0.619|

## Conclusion: 
## No Evidence of another transition in this section
## Let's examine the rest of the series:

e <- e_pr[2200:3153]
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_PR
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-delta0 only,, PR[2200:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.198|  0.303|
# |        2| 0.400|  0.103|
# |        3| 0.335|  0.059|

## Conclusion: 
## No Evidence of another transition!

## Let's try to estimate the TV-2_Trans model on the full dataset:

e <- e_pr
Tobs = NROW(e)
refData = refData_PR
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(0.5,-0.5)
TV$pars["speedN",] = c(4,3)
TV$pars["locN1",] = c(0.3,0.75)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## Let's try to estimate a TV-3_Trans model on the full dataset:

e <- e_pr
Tobs = NROW(e)
refData = refData_PR
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$double1loc, tvshape$single))
TV$delta0 = 2
TV$pars["deltaN",] = c(-0.5,-1,1)
TV$pars["speedN",] = c(4,4,5)
TV$pars["locN1",] = c(0.2,0.4,0.7)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")
# The above specification works, but the last transition may be a double...


# Shape guessed by observing the plot:
TV <- tv(st,c(tvshape$single,tvshape$double1loc,tvshape$double1loc))
TV$delta0 = 5
TV$pars["deltaN",] = c(-1,-1,-1)
TV$pars["speedN",] = c(4,4,4)
TV$pars["locN1",] = c(0.2,0.4,0.75)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-delta0 only,, PR[1:3153]:

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
# Delta0 = 8.0297    se0 = 0.556298*** 
#     
#              st1      se1 sig       st2     se2 sig.1       st3      se3 sig.2
# deltaN -1.312477 0.317171 *** -3.733073 0.39928   *** -2.392260 0.240142   ***
# speedN  3.247389 0.272583 ***  5.501283 0.14505   ***  4.293124 0.266072   ***
# locN1   0.116495 0.030588 ***  0.377289 0.00439   ***  0.804533 0.009102   ***
# locN2         NA      NaN            NA     NaN              NA      NaN      
# 
# Log-likelihood value(TV):  -4758.008



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full PR data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

TVG <- tvgarch(TV,garchtype$gjr)

#TVG$tvOptimcontrol$reltol = 1e-05
#TVG$tvOptimcontrol$ndeps = rep(1e-05,length(TVG$tvOptimcontrol$ndeps))
TVG$garchpars[,1] = c(0.1,0.05,0.5,0.015)
#TVG$garchOptimcontrol$reltol = 1e-03
#TVG$garchOptimcontrol$ndeps = c(1e-05,1e-05,1e-04,1e-05)
TVG$garchOptimcontrol$parscale = c(4,1,40,2)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG)
plot(TVG,main=ptitle)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

plot(TV,main=ptitle)
summary(TVG1)

identical(TVG,TVG1)



