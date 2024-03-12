#  Initialise  ----

rm(list=ls())
library(MTVGARCH)
library(knitr)

#setwd("D:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project4_NAO/R")   #Anna's Laptop - GOOGLE DRIVE
setwd("C:/Source/Repos/tv_betas")
#
estCtrl = list(verbose=TRUE,calcSE=TRUE)
estCtrl = list(verbose=TRUE,calcSE=FALSE)
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

# prices <- read.csv("data/tv_betas_prices.csv")
# e_cba <- diff(log(as.numeric(prices$CBA)) ) * 100  # Percentage Returns
# #
# # Replace all NA's with the previous valid entry if required:
# if(FALSE){
#     na_pos = which(is.na(e_cba) )
#     for(n in seq_along(na_pos)){
#         t_pos <- na_pos[n]
#         e_cba[t_pos] <- e_cba[t_pos-1]
#     }
# }

allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_cba <- allData$CBA

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_cba[2400:3150]
Tobs = NROW(e)
ptitle = "CBA stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1,estCtrl)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_cba[2500:3150]
#              Est      se1 sig
#   omega 0.428363 0.280113    
#   alpha 0.070826 0.033635 ** 
#   beta  0.653320 0.201732 ***

# Try the rolling window method:
if(FALSE){
    e <- e_cba[1:2200]
    estCtrl$vartargetWindow = 400
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    # win=300, a=0.117185, b=0.764039, a+b=0.945
    # win=400, a=0.075521, b=0.839063, a+b=0.915

    # Method:  MLE, variance-targetting a rolling Window of 400 observations 
    #              Est      se1 sig
    #   omega 0.096290       NA    
    #   alpha 0.075521 0.024454 ***
    #   beta  0.839063 0.058954 ***
    #     
    #     Log-likelihood value(GARCH):  -3183.426
    
}

# Next, We need a standard TV object to generate the data:
e <- e_cba
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_CBA <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
ptitle = "CBA"
saveRDS(refData_CBA,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_cba
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_CBA

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

## P-Values from TEST Results, TV-delta0 only, CBA[1:3152]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.013|  0.000|
#     |        2| 0.027|  0.010|
#     |        3| 0.000|  0.009|

## Conclusion: 
## Evidence of a bump first,second and/or third order.
## We can try to identify & estimate a single transition, or split the data and retest:
## Let's try estimating a single order transition...

e <- e_cba
Tobs = NROW(e)
#plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=200),col="grey")
refData = refData_CBA
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$single)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-1_Trans, CBA[1:3153]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.311|  0.163|
#     |        2| 0.012|  0.043|
#     |        3| 0.001|  0.007|

## Conclusion: 
## Evidence of another transition exits
## Let's try estimating a model with 2 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-2_Trans, CBA[1:3153]:

#     | Test Ord|   TR2| Robust|
#     |--------:|-----:|------:|
#     |        1| 0.050|  0.345|
#     |        2| 0.017|  0.003|
#     |        3| 0.019|  0.000|

## Conclusion: 
## Still Evidence of another transition exits - most likely 2nd or 3rd Order
## Let's try estimating a model with 4 x single order transitions...
#
## Didn't work - went back to 3-Trans


## Default starting values don't work well for this 3-Trans model, so...
## Looking at the plot, the missing transition seems to be high-to-low around Obs 2500
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$pars["deltaN",] = c(-1,8,-7)
TV$pars["speedN",] = c(4,6,6)
TV$pars["locN1",] = c(0.4,0.6,0.8)
TV$optimcontrol$reltol = 1e-5
TV$optimcontrol$ndeps = rep(1e-5,10)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-3_Trans, CBA[1:3153]:
# Not Done!


## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 1.242041    se0 = 0.048803*** 
#     
#                  st1      se1 sig       st2      se2 sig.1        st3      se3 sig.2
#     deltaN -0.278313 0.068814 *** 16.439793 4.149025   *** -15.831147 4.176572   ***
#     speedN  4.862726 0.972682 ***  6.693379 0.272620   ***   5.276195 0.591347   ***
#     locN1   0.424439 0.016376 ***  0.710923 0.001560   ***   0.725823 0.005119   ***
#     locN2         NA      NaN            NA      NaN               NA      NaN      
# 
# Log-likelihood value(TV):  -4914.328



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full CBA data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

e <- e_cba
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
ptitle = "CBA"
estCtrl = list(verbose=TRUE,calcSE=TRUE)
#TV <- # Get model spec from "Final Model Specification" above

TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "CBA % Returns"

TVG$tvpars[,1] = c(0.1,2.0,0.33,NA)

TVG$tvpars[,2] = c(8,5.0,0.7,NA)
TVG$tvpars[,3] = c(-7,5.0,0.71,NA)
TVG$tvOptimcontrol$reltol = 5e-04
TVG$tvOptimcontrol$ndeps = rep(1e-04,length(TVG$tvOptimcontrol$ndeps))
TVG$garchpars[,1] = c(0.05,0.02,0.7,0.05)
TVG$garchOptimcontrol$reltol = 1e-04
#TVG$garchOptimcontrol$parscale = c(5,1,50,10)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG@garchObj)
summary(TVG@tvObj)
plot(TVG)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
# Re-run the estimator until fully converged:
TVG1 <- estimateTVGARCH(e,TVG,estCtrl)
# Now Save!
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

summary(TVG1)
plot(TVG1)

# TEST the $Estimated$pars & @obj$Estimated$pars
identical(TVG@tvObj$Estimated$delta0, TVG$Estimated$tv$delta0)
identical(TVG@tvObj$Estimated$delta0_se, TVG$Estimated$tv$delta0_se)
identical(TVG@tvObj$Estimated$pars, TVG$Estimated$tv$pars)
identical(TVG@tvObj$Estimated$se, TVG$Estimated$tv$se)
identical(TVG@tvObj@g, TVG$Estimated$g)
identical(TVG$Estimated$tv$g, TVG$Estimated$g)
#
identical(TVG@garchObj$Estimated$pars, TVG$Estimated$garch$pars)
identical(TVG@garchObj$Estimated$se, TVG$Estimated$garch$se)
identical(TVG@garchObj@h, TVG$Estimated$h)
identical(TVG$Estimated$garch$h, TVG$Estimated$h)
#
identical(TVG,TVG1)

# TEST 1:  Initial Run Only
# Pass

# TEST 2:  Second Run (No Improvement)
# pass

# TEST 3:  Second Run (No Improvement)  Is the TVG identical to previous version?
# pass

# TEST 4:  Second Run (With Improvement)






