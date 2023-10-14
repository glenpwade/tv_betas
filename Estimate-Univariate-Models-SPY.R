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
e_spy <- tail(as.numeric(prices$r_SPY),-1)  # Percentage Returns - drop first (null) observation

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_spy[160:1099]
Tobs = NROW(e)
ptitle = "SPY stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# Next, We need a standard TV object to generate the data:
e <- e_spy
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_SPY <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist, seed=1)
ptitle = "SPY"
saveRDS(refData_SPY,paste0('RefData/',ptitle,'_RefData.RDS'))

#refData = readRDS(paste0('RefData/',ptitle,'_RefData.RDS'))
#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_spy
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
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

## P-Values from TEST Results, TV-delta0 only, SPY[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.062|  0.032|
# |        2| 0.085|  0.000|
# |        3| 0.002|  0.000|

## Conclusion: 
## Strong Evidence of a bump
## Let's try estimating a single order transition...

e <- e_spy
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_SPY
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$single)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-1_Trans, SPY[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.009|  0.037|
# |        2| 0.008|  0.137|
# |        3| 0.002|  0.015|

## Conclusion: 
## Evidence of another transition exits 
## Let's try estimating a model with 2 x single order transitions...

e <- e_spy
Tobs = NROW(e)
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-03
TV$optimcontrol$ndeps = rep(1e-03,length(TV$optimcontrol$ndeps))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

# A TV-2_Trans model cannot be fitted well.
# let's try a TV-3_Trans (Since strongest test evidence points to a 3rd order)

TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
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
## Let's try estimating a model with 2 x single order transitions...

e <- e_spy
Tobs = NROW(e)
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## Default starting values don't work perfectly for this 4-Trans model, so...
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-03
TV$optimcontrol$ndeps = rep(1e-06,length(TV$optimcontrol$ndeps))
TV$pars["deltaN",] = c(-1,1,5,-4)
TV$pars["speedN",] = c(4,4,6,5)
TV$pars["locN1",] = c(0.15,0.4,0.7,0.75)
## Can't seem to fit a 3-transition model...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

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
plot(TV)
abline(v=seq(0,3200,by=100),col="grey")

## P-Values from TEST Results, TV-4_Trans, SPY[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.710|  0.320|
# |        2| 0.627|  0.459|
# |        3| 0.327|  0.188|

## Conclusion: 
## No Evidence of another transition!


## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 3 1 1 3 
# 
# Estimation Results:
#     
# Delta0 = 4.311583    se0 = 0.724877*** 
#     
#            st1      se1 sig       st2      se2 sig.1        st3      se3 sig.2       st4      se4 sig.3
# deltaN 1.55987 0.077772 *** 12.128165 2.455598   *** -12.488008 2.402175   *** -4.926207 0.723412   ***
# speedN 5.81061 0.137957 ***  6.316793 0.185997   ***   4.618113 0.209851   ***  6.935429 0.170850   ***
# locN1  0.48237 0.003227 ***  0.709922 0.000774   ***   0.726778 0.003843   ***  0.904372 0.002577   ***
# locN2       NA      NaN            NA      NaN               NA      NaN              NA      NaN      
# 
# Log-likelihood value(TV):  -4259.378


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full SPY data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

TVG <- tvgarch(TV,garchtype$gjr)

#TVG$tvpars[,1] = c(-1.5,5,0.05,NA)
#TVG$tvpars[,2] = c(-0.7,2.5,0.25,NA)
TVG$tvOptimcontrol$reltol = 1e-05
TVG$tvOptimcontrol$ndeps = rep(1e-05,length(TVG$tvOptimcontrol$ndeps))
#TVG$garchpars[,1] = c(0.02,0.002,0.8,0.15)
#TVG$garchOptimcontrol$reltol = 1e-03
TVG$garchOptimcontrol$ndeps = c(1e-09,1e-03,1e-09,1e-09)
TVG$garchOptimcontrol$parscale = c(20,1,500,100)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG)
plot(TVG,main=ptitle)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

plot(TVG1,main=ptitle)
summary(TVG1)

identical(TVG,TVG1)



