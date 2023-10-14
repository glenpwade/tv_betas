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
e_stw <- diff(log(as.numeric(prices$STW)) ) * 100  # Percentage Returns
#
# Replace all NA's with the previous valid entry if required:
if(FALSE){
    na_pos = which(is.na(e_stw) )
    for(n in seq_along(na_pos)){
        t_pos <- na_pos[n]
        e_stw[t_pos] <- e_stw[t_pos-1]
    }
}

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_stw[1320:2225]
Tobs = NROW(e)
ptitle = "STW stable subest"
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
#
Garch1 <- garch(garchtype$general)
Garch1 = estimateGARCH(e,Garch1)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_stw[1400:1800]
# Est se1 sig
# omega 0.072941 NaN    
# alpha 0.070489 NaN    
# beta  0.852360 NaN 

# e_stw[2500:3150]
# Est      se1 sig
# omega 0.132530 0.088212    
# alpha 0.061153 0.027582 ** 
# beta  0.839795 0.083945 ***

# Try the rolling window method:
if(FALSE){
    estCtrl$vartargetWindow = 300
    Garch1 <- garch(garchtype$general)
    Garch1_rollWin <- estimateGARCH_RollingWindow(e,Garch1,estCtrl)
    summary(Garch1_rollWin)
    Garch1_rollWin$Estimated$pars["alpha",1] + Garch1_rollWin$Estimated$pars["beta",1]
    Garch1$pars <- Garch1_rollWin$Estimated$pars  # The starting pars are used in the generateRefData fn
    # win=300, a=0.097559, b=0.847240, a+b=0.945
    # win=400, a=0.082415, b=0.894717, a+b=0.977
    # win=500, a=0.078079, b=0.905088, a+b=0.983
    
    # Method:  MLE, variance-targetting a rolling Window of 300 observations 
    # Est      se1 sig
    # omega 0.073229       NA    
    # alpha 0.097559 0.023529 ***
    # beta  0.847240 0.049561 ***
}

# Next, We need a standard TV object to generate the data:
e <- e_stw
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_STW <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist, seed=1)
ptitle = "STW"
saveRDS(refData_STW,paste0('RefData/',ptitle,'_RefData.RDS'))

#  End of Ref Data Generation 

# START - Test for any transition ----

e <- e_stw
Tobs = NROW(e)
plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
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

## P-Values from TEST Results, TV-delta0 only, STW[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.870|  0.850|
# |        2| 0.895|  0.647|
# |        3| 0.086|  0.082|

## Conclusion: 
## No Evidence of a bump - but the plot suggests otherwise
## We can try to identify & estimate a single transition, or split the data and retest:
## Let's try estimating a single order transition...

e <- e_stw
Tobs = NROW(e)
#plot(e,type='l',main=ptitle)
abline(v=seq(0,3200,by=100),col="grey")
refData = refData_STW
# create TV object and estimate
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$single)
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-1_Trans, STW[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.008|  0.002|
# |        2| 0.013|  0.008|
# |        3| 0.003|  0.017|

## Conclusion: 
## Evidence of another transition exits
## Let's try estimating a model with 2 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-7
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-2_Trans, STW[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.086|  0.224|
# |        2| 0.179|  0.474|
# |        3| 0.013|  0.049|


## Conclusion: 
## Still Evidence of another transition exits
## Let's try estimating a model with 3 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## Default starting values don't work perfectly for this 3-Trans model, so...
TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-03
TV$optimcontrol$ndeps = rep(1e-06,10)
TV$pars["deltaN",] = c(-1,5,-4)
TV$pars["speedN",] = c(4,6,5)
TV$pars["locN1",] = c(0.45,0.7,0.75)
## Can't seem to fit a 3-transition model...
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

# By observing the plot, we can try the following shape:
TV <- tv(st,c(tvshape$single,tvshape$double1loc,tvshape$single,tvshape$single))
TV$optimcontrol$reltol = 1e-05
TV$optimcontrol$ndeps = rep(1e-07,length(TV$optimcontrol$ndeps))
TV$optimcontrol$parscale = c(1,1,3,1,1,3,1,5,5,1,5,5,1)
TV$pars["deltaN",] = c(-0.5,-0.5,9,-7)
TV$pars["speedN",] = c(5,2,6,5)
TV$pars["locN1",] = c(0.05,0.3,0.7,0.75)
#
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-3_Trans, STW[1:3153]:

# | Test Ord|   TR2| Robust|
# |--------:|-----:|------:|
# |        1| 0.313|  0.409|
# |        2| 0.350|  0.664|
# |        3| 0.358|  0.768|


## Conclusion: 
## No Evidence of another transition!


## Final Model Specification:  ----


# TV OBJECT
# 
# Transition Shapes: 1 3 1 1 
# 
# Estimation Results:
#     
# Delta0 = 2.530764    se0 = 0.396564*** 
#     
#              st1      se1 sig       st2      se2 sig.1       st3      se3 sig.2        st4      se4 sig.3
# deltaN -1.228953 0.417568 *** -0.792440 0.122200   *** 13.285247 1.855377   *** -13.099437 1.865727   ***
# speedN  6.980272 1.929797 ***  5.178263 0.243529   ***  6.996496 1.162175   ***   5.445603 0.103166   ***
# locN1   0.058608 0.001501 ***  0.364195 0.008154   ***  0.712345 0.001236   ***   0.724423 0.001583   ***
# locN2         NA      NaN            NA      NaN              NA      NaN               NA      NaN      
# 
# Log-likelihood value(TV):  -3972.32




##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full NAB data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.

TVG <- tvgarch(TV,garchtype$gjr)

TVG$tvpars[,1] = c(-1.5,5,0.05,NA)
TVG$tvpars[,2] = c(-0.7,2.5,0.25,NA)
TVG$tvOptimcontrol$reltol = 5e-05
TVG$tvOptimcontrol$ndeps = rep(1e-04,length(TVG$tvOptimcontrol$ndeps))
TVG$garchpars[,1] = c(0.02,0.002,0.8,0.15)
TVG$garchOptimcontrol$reltol = 1e-03
TVG$garchOptimcontrol$ndeps = c(1e-09,1e-03,1e-09,1e-09)
TVG$garchOptimcontrol$parscale = c(50,1,300,65)

TVG <- estimateTVGARCH(e,TVG,estCtrl)
summary(TVG)
plot(TVG,main=ptitle)   # Note: produces 2 plots: sqrt(g)  &  sqrt(h)  
saveRDS(TVG,paste0('Results/',ptitle,'_Final_TVG_model.RDS'))
#
# Reload the saved TVG object:
TVG1 <- readRDS(paste0('Results/',ptitle,'_Final_TVG_model.RDS'))

plot(TV,main=ptitle)
summary(TVG1)

identical(TVG,TVG1)



