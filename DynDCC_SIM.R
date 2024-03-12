# # Make some data dynamic DCC:
# N = 2 or 5
# Tobs=1000
# rho = 0.5 or 0.9
# Q1 = diag(N)
# Q2 = toeplitz(rho^seq.int(0, N-1))
# a= 0.05
# b= 0.9
# speed= 2.5, 5
# loc = 0.5
gc()

library(MTVGARCH)
library(doParallel)
library(knitr)

setwd("C:/Source/Repos/tv_betas") 
dir.create("SimResults")

SimRuns = 2000
Tobs=1000
estCtrl = list(verbose=FALSE,calcSE=FALSE)

## Sim START ####

N = 2
# Q1 = 0.0
Q1 <- diag(N)
q1name="00"
# Q2 = toeplitz(rho^seq.int(0, N-1))
rho = 0.9
rhoname = "090"
Q2 <- toeplitz(rho^seq.int(0, N-1))
q2name="09"
Qbar <- toeplitz(rho^seq.int(0, N-1))
# Speed (Eta[0..7])
speed = 2.5
spdname = "25"
# Loc (0..1)
loc = 0.5
locname = "05"
# Select a = 0.1 , 0.05 , 0.025 AND aname = 100 , 050 , 025
a=0.050
aname = "050"
# Select b = 0.8 , 0.9 , 0.95 AND bname = 800 , 900 , 950
b=0.800
bname="800"
#
st = seq(-0.5,0.5,length.out=Tobs)

numCores = 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

timestamp()
dccTest= foreach(n = 1:SimRuns, .combine = rbind, .inorder = FALSE, .packages = "MTVGARCH")%dopar%{
    #e <- generateDCCRefData(N,Tobs,Qbar,a,b)
    e <- generateDynDCCRefData(N,Tobs,Q1,Q2,speed,loc,a,b)

    TV = list()
    MTV = list()
    Names = vector()
    for(i in 1:N){
        TV[[i]] <- tv(st,tvshape$delta0only)
        TV[[i]] <- estimateTV(e[,i],TV[[i]],estCtrl)
        MTV[[i]] = tvgarch(TV[[i]],garchtype$noGarch)
        MTV[[i]] <- estimateTVGARCH(e[,i],MTV[[i]],estCtrl)
        Names <- c(Names,i)
    }
    NTV = ntvgarch(MTV, as.character(Names))
    cccObj = ccc(NTV)
    cccObj = estimateCCC(cccObj,estimationCtrl = estCtrl)
    
    #testOrd=1
    #ccctest <- test.CCCParsim(cccObj,st,testOrd)
    testOrd=1
    stcctest1 = test.CCCvSTCC1(cccObj,st,testOrd)
    testOrd=2
    stcctest2 = test.CCCvSTCC1(cccObj,st,testOrd)
    # returns statistic and pvalue
    c(as.integer(n),stcctest1["Statistic",1],stcctest1["P_value",1],stcctest2["Statistic",1],stcctest2["P_value",1])  # stat1 ~ pval ~ stat2 ~ pval
    
}
timestamp()

stopCluster(cl)

saveName = paste0("Pilot_N",N,"_T",Tobs,"_Q1_diag","_Q2_Tplz",rhoname,"_a",aname,"_b",bname,".RDS")
saveRDS(dccTest,file=saveName)

# 
gc()
## Sim 2 ####
N = 5
# Q1 = 0.0
Q1 <- diag(N)
q1name="00"
# Q2 = toeplitz(rho^seq.int(0, N-1))
rho = 0.9
rhoname = "090"
Q2 <- toeplitz(rho^seq.int(0, N-1))
q2name="09"
Qbar <- toeplitz(rho^seq.int(0, N-1))
# Speed (Eta[0..7])
speed = 2.5
spdname = "25"
# Loc (0..1)
loc = 0.5
locname = "05"
# Select a = 0.1 , 0.05 , 0.025 AND aname = 100 , 050 , 025
a=0.050
aname = "050"
# Select b = 0.8 , 0.9 , 0.95 AND bname = 800 , 900 , 950
b=0.800
bname="800"
#
st = seq(-0.5,0.5,length.out=Tobs)

numCores = 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

timestamp()
dccTest= foreach(n = 1:SimRuns, .combine = rbind, .inorder = FALSE, .packages = "MTVGARCH")%dopar%{
    #e <- generateDCCRefData(N,Tobs,Qbar,a,b)
    e <- generateDynDCCRefData(N,Tobs,Q1,Q2,speed,loc,a,b)
    
    TV = list()
    MTV = list()
    Names = vector()
    for(i in 1:N){
        TV[[i]] <- tv(st,tvshape$delta0only)
        TV[[i]] <- estimateTV(e[,i],TV[[i]],estCtrl)
        MTV[[i]] = tvgarch(TV[[i]],garchtype$noGarch)
        MTV[[i]] <- estimateTVGARCH(e[,i],MTV[[i]],estCtrl)
        Names <- c(Names,i)
    }
    NTV = ntvgarch(MTV, as.character(Names))
    cccObj = ccc(NTV)
    cccObj = estimateCCC(cccObj,estimationCtrl = estCtrl)
    
    #testOrd=1
    #ccctest <- test.CCCParsim(cccObj,st,testOrd)
    testOrd=1
    stcctest1 = test.CCCvSTCC1(cccObj,st,testOrd)
    testOrd=2
    stcctest2 = test.CCCvSTCC1(cccObj,st,testOrd)
    # returns statistic and pvalue
    c(as.integer(n),stcctest1["Statistic",1],stcctest1["P_value",1],stcctest2["Statistic",1],stcctest2["P_value",1])  # stat1 ~ pval ~ stat2 ~ pval
    
}
timestamp()

stopCluster(cl)

saveName = paste0("Pilot_N",N,"_T",Tobs,"_Q1_diag","_Q2_Tplz",rhoname,"_a",aname,"_b",bname,".RDS")
saveRDS(dccTest,file=saveName)

# 
gc()
## Sim 3 ####
N = 10
# Q1 = 0.0
Q1 <- diag(N)
q1name="00"
# Q2 = toeplitz(rho^seq.int(0, N-1))
rho = 0.9
rhoname = "090"
Q2 <- toeplitz(rho^seq.int(0, N-1))
q2name="09"
Qbar <- toeplitz(rho^seq.int(0, N-1))
# Speed (Eta[0..7])
speed = 2.5
spdname = "25"
# Loc (0..1)
loc = 0.5
locname = "05"
# Select a = 0.1 , 0.05 , 0.025 AND aname = 100 , 050 , 025
a=0.050
aname = "050"
# Select b = 0.8 , 0.9 , 0.95 AND bname = 800 , 900 , 950
b=0.800
bname="800"
#
st = seq(-0.5,0.5,length.out=Tobs)

numCores = 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

timestamp()
dccTest= foreach(n = 1:SimRuns, .combine = rbind, .inorder = FALSE, .packages = "MTVGARCH")%dopar%{
    #e <- generateDCCRefData(N,Tobs,Qbar,a,b)
    e <- generateDynDCCRefData(N,Tobs,Q1,Q2,speed,loc,a,b)
    
    TV = list()
    MTV = list()
    Names = vector()
    for(i in 1:N){
        TV[[i]] <- tv(st,tvshape$delta0only)
        TV[[i]] <- estimateTV(e[,i],TV[[i]],estCtrl)
        MTV[[i]] = tvgarch(TV[[i]],garchtype$noGarch)
        MTV[[i]] <- estimateTVGARCH(e[,i],MTV[[i]],estCtrl)
        Names <- c(Names,i)
    }
    NTV = ntvgarch(MTV, as.character(Names))
    cccObj = ccc(NTV)
    cccObj = estimateCCC(cccObj,estimationCtrl = estCtrl)
    
    #testOrd=1
    #ccctest <- test.CCCParsim(cccObj,st,testOrd)
    testOrd=1
    stcctest1 = test.CCCvSTCC1(cccObj,st,testOrd)
    testOrd=2
    stcctest2 = test.CCCvSTCC1(cccObj,st,testOrd)
    # returns statistic and pvalue
    c(as.integer(n),stcctest1["Statistic",1],stcctest1["P_value",1],stcctest2["Statistic",1],stcctest2["P_value",1])  # stat1 ~ pval ~ stat2 ~ pval
    
}
timestamp()

stopCluster(cl)

saveName = paste0("Pilot_N",N,"_T",Tobs,"_Q1_diag","_Q2_Tplz",rhoname,"_a",aname,"_b",bname,".RDS")
saveRDS(dccTest,file=saveName)

# 






## Review the output:   ####

dccTest = readRDS("Pilot_N10_T1000_Q1_diag_Q2_Tplz090_a100_b800.RDS")
dccTest = readRDS("Pilot_N5_T1000_Q1_diag_Q2_Tplz090_a100_b800.RDS")
dccTest = readRDS("Pilot_N2_T1000_Q1_diag_Q2_Tplz090_a100_b800.RDS")

hist(dccTest[,2],breaks = 25)
sum(as.integer( dccTest[,3]<0.01) )/nrow(dccTest)
sum(as.integer( dccTest[,3]<0.05) )/nrow(dccTest)
sum(as.integer( dccTest[,3]<0.10) )/nrow(dccTest)
summary(dccTest[,3])


