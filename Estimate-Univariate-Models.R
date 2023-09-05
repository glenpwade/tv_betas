#  Initialise  ----

rm(list=ls())
library(MTVGARCH)

#setwd("D:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project4_NAO/R")   #Anna's Laptop - GOOGLE DRIVE
setwd("C:/Source/Repos/tv_betas")

# Get Data  ----

anz <- read.csv("data/anz.csv")
e_anz <- diff(log(as.numeric(anz$Close)) )
#
# Replace all NA's with the previous valid entry:
e_anz_na_pos = which(is.na(e_anz) )
for(n in seq_along(e_anz_na_pos)){
    t_pos <- e_anz_na_pos[n]
    e_anz[t_pos] <- e_anz[t_pos-1]
}

# Generate suitable reference Data ----
estCtrl$vartargetWindow = 400
G1 <- garch(garchtype$general)
G1_rw <- estimateGARCH_RollingWindow(e,G1,estCtrl)
summary(G1_rw)
# a=0.096, b=0.854

# create TV object and estimate
Tobs = NROW(e)
st =  (1:Tobs)/Tobs
TV <- tv(st,tvshape$delta0only)

G1$pars <- G1_rw$Estimated$pars

anz_refdata <- generateRefData(100,6023,TV,G1,corrObj = NULL, noiseDist = "Normal", seed=2)


# Univariate tv estimations ----
e <- e_anz
ptitle = "ANZ"

if(TRUE)
{
  # create TV object and estimate
  Tobs = NROW(e)
  st =  (1:Tobs)/Tobs
  TV <- tv(st,tvshape$delta0only)
  estCtrl = list(verbose=TRUE,calcSE=TRUE)
  TV <- estimateTV(e,TV,estCtrl)
  
  ## Calculate p_value by simulating the distrbution:
  
  # Set simulation controls
  simcontrol = list()
  simcontrol$numLoops = 1000
  simcontrol$numCores = 6
  
  
  simcontrol$saveAs = paste0("Simdist_",ptitle,"_TV",TV@nr.transitions,"Trans.RDS")
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
  
  ## RESTART HERE with an updated TV Model specification: ####
  
  maxTestOrd = 3
  simcontrol$maxTestorder=maxTestOrd
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
 
  RefTests <- getTestStats(e,TV,maxTestOrd)
  
  RefTests = list()
  for (n in  1:maxTestOrd){
      RefTests[[n]] <- list()
      RefTests[[n]] <- getTestStats(e,TV,n)
  }
 
  # Note: testStatDist() will generate a simulated test distribution and calculate p_values as follows:
  # 1. For each column of refData, estimate the null TV model, then calculate test statistics (TR2 & Robust) for each requested TestOrder
  # 2. Then compare actual RefTest statistics with the TestStat-distributions (i.e. get p-values)
  # 3. This takes approx. 3 minutes for the TV$delta0 specification running in local-cpu parallel mode on an i7 3GHz 4-core/8-processor CPU,  Win10PC
  SIMRESULT = testStatDist(refData,TV,RefTests,simcontrol)
  
  # Print the test Results:
  if(TRUE){
      numDP = 5
      if(maxTestOrd >= 1) cat("\nTestOrder1:\n  P_Value(TR2): ",format(SIMRESULT$Order1$pVal_TR2,nsmall= numDP)," P_Value(Robust):",format(SIMRESULT$Order1$pVal_ROB,nsmall= numDP))
      if(maxTestOrd >= 2) cat("\nTestOrder2:\n  P_Value(TR2): ",format(SIMRESULT$Order2$pVal_TR2,nsmall= numDP)," P_Value(Robust):",format(SIMRESULT$Order2$pVal_ROB,nsmall= numDP))
      if(maxTestOrd >= 3) cat("\nTestOrder3:\n  P_Value(TR2): ",format(SIMRESULT$Order3$pVal_TR2,nsmall= numDP)," P_Value(Robust):",format(SIMRESULT$Order3$pVal_ROB,nsmall= numDP))
      if(maxTestOrd >= 4) cat("\nTestOrder4:\n  P_Value(TR2): ",format(SIMRESULT$Order4$pVal_TR2,nsmall= numDP)," P_Value(Robust):",format(SIMRESULT$Order4$pVal_ROB,nsmall= numDP))  
  }
  
  # View the Distributions, if you like: 
  if(FALSE){
      hist(SIMRESULT$TestStatDist[,"Stat_TR2.2"],breaks = 20)
      abline(v=RefTests[[3]]$TR2,col="red")
      hist(SIMRESULT$TestStatDist[,"Stat_Robust.2"],breaks = 20)
      abline(v=RefTests[[3]]$Robust,col="red")
      
      #View(SIMRESULT$TestStatDist)     # Have a look here for Troubleshooting
  }
  
  ## SAVE Results
  saveFile = paste0("Results/TestResults_",ptitle,ptime,"_TV",TV@nr.transitions,"Trans.RDS")
  #saveRDS(SIMRESULT,saveFile)
  
  ## RELOAD Results
  #SIMRESULT = readRDS(saveFile)

}

rm(TV,MTV,e)

