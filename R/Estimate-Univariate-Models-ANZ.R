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

# prices <- read.csv("data/tv_betas_prices.csv")
# e_anz <- diff(log(as.numeric(prices$ANZ)) ) * 100  # Percentage Returns
# #
# # Replace all NA's with the previous valid entry if required:
# if(FALSE){
#     e_anz_na_pos = which(is.na(e_anz) )
#     for(n in seq_along(e_anz_na_pos)){
#         t_pos <- e_anz_na_pos[n]
#         e_anz[t_pos] <- e_anz[t_pos-1]
#     }
# }

allData <- readRDS("Data/Returns_USlagged_4B.RDS")
e_anz <- allData$ANZ

# Generate suitable reference Data (once only) ----
# Visually find the longest stable-variance subset of data - to determine best Garch pars
e <- e_anz
plot(e,type='l')
abline(v=c(1350,2200),col="red")
abline(v=c(2400,3152),col="red")
#e <- e_anz[1350:2200]
e <- e_anz[2400:3152]
Tobs = NROW(e)
ptitle = "ANZ stable subset"
plot(e,type='l',main=ptitle)
#
Garch1 <- garch(garchtype$general)
Garch1 <- estimateGARCH(e,Garch1,estCtrl)
summary(Garch1)
Garch1$pars <- Garch1$Estimated$pars  # The starting pars are used in the generateRefData fn

# e_anz[2400:3152]
# Est      se1 sig
# omega 0.032604 0.008911 ***
# alpha 0.074835 0.009786 ***
# beta  0.906791 0.012940 ***

# Next, We need a standard TV object to generate the data:
e <- e_anz
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
TV <- tv(st,tvshape$delta0only)

refData_ANZ <- generateRefData(simcontrol$numLoops,Tobs,TV,Garch1,corrObj = NULL, noiseDist = noisedist)
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

## P-Values from TEST Results, TV-delta0 only, ANZ[1:3152]:

##   | Test Ord|   TR2| Robust|
#    |--------:|-----:|------:|
#    |        1| 0.087|  0.028|
#    |        2| 0.129|  0.089|
#    |        3| 0.003|  0.034|

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

## P-Values from TEST Results, TV-1_Trans, ANZ[1:3152]:

#   | Test Ord|   TR2| Robust|
#   |--------:|-----:|------:|
#   |        1| 0.026|  0.148|
#   |        2| 0.005|  0.051|
#   |        3| 0.000|  0.010|

## Conclusion: 
## Evidence of another transition exits
## Let's try estimating a model with 2 x single order transitions...

TV <- tv(st,c(tvshape$single,tvshape$single))
TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")

## P-Values from TEST Results, TV-2_Trans, ANZ[1:3152]:

#    | Test Ord|   TR2| Robust|
#    |--------:|-----:|------:|
#    |        1| 0.242|  0.948|
#    |        2| 0.069|  0.162|
#    |        3| 0.004|  0.004|

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
TV$pars["deltaN",] = c(-0.5,6,-6)
TV$pars["speedN",] = c(4,5,5)
TV$pars["locN1",] = c(0.4,0.6,0.7)
TV$optimcontrol$reltol = 1e-07

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)

## P-Values from TEST Results, TV-3_Trans, ANZ[1:3152]:

#    | Test Ord|   TR2| Robust|
#    |--------:|-----:|------:|
#    |        1| 0.093|  0.195|
#    |        2| 0.101|  0.364|
#    |        3| 0.077|  0.038|
#
# Log-likelihood value(TV):  -5211.561

## Conclusion: 
## Some small Evidence of another transition exits
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
TV$pars["deltaN",] = c(0.5,3,7,-9)
TV$pars["speedN",] = c(5,3,5,4)
TV$pars["locN1",] = c(0.35,0.6,0.7,0.75)
TV$optimcontrol$reltol = 1e-7
TV$optimcontrol$ndeps = rep(1e-9,13)

TV <- estimateTV(e,TV,estCtrl)
summary(TV)
plot(TV)
#
## P-Values from TEST Results, TV-4_Trans, ANZ[1:3152]:

# Log-likelihood value(TV):  -5249.475
    
## Conclusion: 
## The 4-Trans model is hard to estimate with SE.  The ll_value is worse than 3-Trans.
## Considering the Test only barely indicated a 4th Transition, we will stop at 3.
## so we will stop here.

## Final Model Specification:  ----

# TV OBJECT
# 
# Transition Shapes: 1 1 1 1 
# 
# Delta0 = 1.765208    se0 = 0.073322*** 
#
#               st1      se1 sig       st2      se2 sig.1       st3      se3 sig.2
#  deltaN -0.577503 0.094753 *** 10.653926 1.053234   *** -9.906067 0.988327   ***
#  speedN  5.973245 1.576291 ***  6.999987      NaN        4.392746 0.127020   ***
#  locN1   0.420742 0.004207 ***  0.711766 0.000618   ***  0.742545 0.003914   ***
#  locN2         NA      NaN            NA      NaN              NA      NaN      

# Log-likelihood value(TV):  -5211.561



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# START tvgarch Specification ----

# We have 'g' ignoring Garch, now we need to find 'g' & 'h'

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

# Make sure 'e' is set to the full ANZ data set!
# Run the estimation using default starting params & optim-controls
# Use the results to fine tune above if needed.
ptitle = "ANZ"
e <- e_anz
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
# TV = get final specification from code above
#
TVG <- tvgarch(TV,garchtype$gjr)
TVG$e_desc = "ANZ % Returns"

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

## Final Model Specification ####

# -- TVGARCH Model Specification --
#     
#     Multiplicative Model Log-Likelihood Value:  -5032.714
# 
# TVGARCH Model Parameters:
#     GARCH OBJECT
# 
# Type:  GJR Garch
# Order: ( 1 , 1 )
# Estimation Results:
#     
#     Method:  MLE 
#                Est      se1 sig
#     omega 0.015541 0.003691 ***
#     alpha 0.014334 0.006879 ** 
#     beta  0.930682 0.010235 ***
#     gamma 0.075482 0.010818 ***
#     
#     Log-likelihood value(GARCH):  -5032.714
# 
# TV OBJECT
# 
# Transition Shapes: 1 1 1 
# 
# Estimation Results:
#     
#     Delta0 = 1.765208    se0 = NaN 
# 
#              st1      se1 sig      st2      se2 sig.1       st3      se3 sig.2
# deltaN -1.144102 0.423570 *** 3.724952 1.743657   **  -2.375671 1.720534      
# speedN  5.162228 0.740163 *** 4.179354 0.448528   ***  6.999983      NaN      
# locN1   0.513852 0.008638 *** 0.599475 0.019978   ***  0.629153 0.001805   ***
#     locN2         NA      NaN           NA      NaN              NA      NaN      
# 
# Log-likelihood value(TV):  -5032.96
# 
# -- End of TVGARCH Model Specification --